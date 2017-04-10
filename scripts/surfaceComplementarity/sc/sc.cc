//////////////////////////////////////////////////////////////////////
// GPU-SC: GPU-accelerated version of the original Lawrence &
// Coleman shape complementarity program from CCP4.
// Luki Goldschmidt <luki@mbi.ucla.edu>, March 2011 
//////////////////////////////////////////////////////////////////////
// Based on Fortran code:
// Sc (Version 2.0): A program for determining Shape Complementarity
// Copyright Michael Lawrence, Biomolecular Research	Institute
// 343 Royal Parade Parkville Victoria Australia
// This program is designed to compute Sc between two molecules.
//////////////////////////////////////////////////////////////////////

#include <iostream>
#include <ostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <stdlib.h>
#include "sc.h"

using namespace std;
std::vector<SCCALC_ATOM_RADIUS> ScCalc::radii;	// static

#define PI 3.14159265

////////////////////////////////////////////////////////////
// C functions

// Atom distance callback for sort 
PDBAtom *_atom_distance_ref = NULL;
int _atom_distance_cb(void *a1, void *a2)
{	
	float d1 = _atom_distance_ref->Distance(*((PDBAtom*)a1));
	float d2 = _atom_distance_ref->Distance(*((PDBAtom*)a2));
	return d1 < d2 ? -1: (d2 > d1 ? 1 : 0); 
}

////////////////////////////////////////////////////////////
// Public class functions

ScCalc::ScCalc()
{
	memset(&settings, 0, sizeof(settings));
	memset(&results, 0, sizeof(results));

	// Set defaults
	settings.sc_radii = "sc_radii.lib";

	// Defaults from sc source
	settings.rp = 1.7;
	settings.density = 15.0;
	settings.band = 1.5;
	settings.sep = 8.0;
	settings.weight = 0.5;
	settings.binwidth_dist = 0.02;
	settings.binwidth_norm = 0.02;

#ifdef GPU
	settings.gpu.device = -1;
#endif

	// Streams for messages
	verbose = NULL;
	error = NULL;
}

int ScCalc::Init() {
	if(radii.empty())
		ReadScRadii(settings.sc_radii);
#ifdef GPU
	gpuInit();
#endif
	return 1;
}

ScCalc::~ScCalc() {
}

// Reset data for another calculation
void ScCalc::Reset()
{
	memset(&results, 0, sizeof(results));
	atoms.clear();
	dots[0].clear();
	dots[1].clear();
	probes.clear();
}

// Read atom radius definitions from file
int ScCalc::ReadScRadii(const char *fn)
{
	SCCALC_ATOM_RADIUS radius;
	FILE *in = fopen(fn, "r");
	
	if(!in)
		throw ScCalcException("Cannot read radii definition file: %s", fn);
		
	radii.clear();
	while(!feof(in)) {
		if(fscanf(in, "%4c %4c %f\n",
			radius.residue, 
			radius.atom, 
			&radius.radius
		) == 3)
			radii.push_back(radius);
	}
	fclose(in);

	VERBOSE("Atom radii read: " << radii.size() << std::endl);

	return !radii.empty();
}

// Read atom coordinates from PDB file
std::vector<PDBAtom> ScCalc::ReadPdb(const char *fn)
{
	PDBAtom atom;
	std::vector<PDBAtom> atoms;
	char line[80];
	FILE *in; 
	
	in = fopen(fn, "r");
	if(!in)
		throw ScCalcException("Cannot read PDB file: %s", fn);

	atoms.clear();
	while(fgets(line, sizeof(line), in)) {
		if(sscanf(line, "ATOM %5d %3c %3c %c %4d %f %f %f %f %f",
			&atom.natom, 
			atom.atom, 
			atom.residue, 
			&atom.chain, 
			&atom.nresidue, 
			&atom.x, 
			&atom.y, 
			&atom.z, 
			&atom.occ, 
			&atom.b 
		) == 10) {
			atoms.push_back(atom);
		}
	}
	
	fclose(in);

	return atoms;
}

// Add an atom to a molecule for computation,
// also looks-up the atom radius and density
int ScCalc::AddAtom(
		int molecule,
		PDBAtom &atom)
{
	if(AssignAtomRadius(atom)) {
		molecule = (molecule == 1);
		atom.density = settings.density;
		atom.molecule = molecule;
		atom.natom = ++results.nAtoms;
		atom.access = 0;
		atoms.push_back(atom);
		results.surface[molecule].nAtoms++;
		/*
		printf("AddAtom[%d] %d: %s:%s (%10.4f, %10.4f, %10.4f) = %.4f\n", molecule,
		       results.surface[molecule].nAtoms,
		       atom.residue, atom.atom, atom.x, atom.y, atom.z, atom.radius);
		*/
		return 1;
		
	} else {
		if(error)
		*error << "Failed to assign atom radius for residue " << atom.chain << ":" << atom.nresidue
			<< " = " << atom.residue << ":" << atom.atom
			<< ". Skipping atom!" << std::endl;
	}
	return 0;
}

// Run the SC calculation for defined molecules
SCCALC_RESULTS *ScCalc::CalcSc()
{
	if(atoms.empty())
		throw ScCalcException("No atoms defined");
	if(!results.surface[0].nAtoms)
		throw ScCalcException("No atoms defined for molecule 1");
	if(!results.surface[1].nAtoms)
		throw ScCalcException("No atoms defined for molecule 2");

	// Determine assign the attention numbers for each atom		
	AssignAttentionNumbers(atoms);
	
	// Now compute the surface for the atoms in the interface and its neighbours
	VERBOSE("Generating molecular surface, " << settings.density << " dots/A^2" << std::endl);

	CalcDotsForAllAtoms(atoms);
	VERBOSE("		   Convex dots: " << results.dots.convex << std::endl);
	VERBOSE("		 Toroidal dots: " << results.dots.toroidal << std::endl);
	VERBOSE("		  Concave dots: " << results.dots.concave << std::endl);
	VERBOSE("Total surface dots (1): " << dots[0].size() << std::endl);
	VERBOSE("Total surface dots (2): " << dots[1].size() << std::endl);
	VERBOSE("	Total surface dots: " << (dots[0].size()+dots[1].size()) << std::endl);

	// Cut away the periphery of each surface
	VERBOSE("Trimming peripheral band, " << settings.band << "A range" << std::endl);
	vector<const SCCALC_DOT*> trimmed_dots[2];
	for(int i = 0; i < 2; i++) {
		results.surface[i].trimmedArea = TrimPeripheralBand(dots[i], trimmed_dots[i]);
		results.surface[i].nTrimmedDots = trimmed_dots[i].size();
		results.surface[i].nAllDots = dots[i].size();
	}

	// Compute distance arrays and histograms for each surface
	VERBOSE("Computing surface separation and vectors" << std::endl);

	CalcNeighborDistance(0, trimmed_dots[0], trimmed_dots[1]);
	CalcNeighborDistance(1, trimmed_dots[1], trimmed_dots[0]);

	results.surface[2].d_mean = (results.surface[0].d_mean + results.surface[1].d_mean) / 2;
	results.surface[2].d_median = (results.surface[0].d_median + results.surface[1].d_median) / 2;
	results.surface[2].s_mean = (results.surface[0].s_mean + results.surface[1].s_mean) / 2;
	results.surface[2].s_median = (results.surface[0].s_median + results.surface[1].s_median) / 2;

	results.surface[2].nAtoms = (results.surface[0].nAtoms + results.surface[1].nAtoms);
	results.surface[2].nBuriedAtoms = (results.surface[0].nBuriedAtoms + results.surface[1].nBlockedAtoms);
	results.surface[2].nBlockedAtoms = (results.surface[0].nBuriedAtoms + results.surface[1].nBuriedAtoms);
	results.surface[2].nAllDots = (results.surface[0].nAllDots + results.surface[1].nAllDots);
	results.surface[2].nTrimmedDots = (results.surface[0].nTrimmedDots + results.surface[1].nTrimmedDots);
	//results.surface[2].nBuriedDots = (results.surface[0].nBuriedDots + results.surface[1].nBuriedDots);
	//results.surface[2].nAccessibleDots = (results.surface[0].nAccessibleDots + results.surface[1].nAccessibleDots);
	results.surface[2].trimmedArea = results.surface[0].trimmedArea + results.surface[1].trimmedArea;

	results.sc = results.surface[2].s_median;
	results.separation = results.surface[2].d_median;
	results.area = results.surface[2].trimmedArea;

	return &results;
}

////////////////////////////////////////////////////////////
// Protected class functions

// Look up the atom radius for an atom
int ScCalc::AssignAtomRadius(PDBAtom &atom)
{
	std::vector<SCCALC_ATOM_RADIUS>::const_iterator radius;

	// Assign radius with wildcard matching
	for(radius = radii.begin(); radius != radii.end(); radius++) {
		if(WildcardMatch(atom.residue, radius->residue) &&
			WildcardMatch(atom.atom, radius->atom)) {
				atom.radius = radius->radius;
				return 1;
		}
	}

	return 0;
}

// Inline residue and atom name matching function
int ScCalc::WildcardMatch(
		const char *r,
		const char *pattern)
{
	while(*pattern && *r) {
		if(*pattern != '*' && *r != *pattern)
			return 0;
		r++;
		pattern++;
	}
	return 1;
}

// Determine assign the attention numbers for each atom
int ScCalc::AssignAttentionNumbers(std::vector<PDBAtom> &atoms)
{
	std::vector<PDBAtom>::iterator pAtom1, pAtom2;

	for(pAtom1 = atoms.begin(); pAtom1 < atoms.end(); pAtom1++) {
		// find nearest neighbour in other molecule
		float dist_min = 99999.0, r;
		for(pAtom2 = atoms.begin(); pAtom2 < atoms.end(); pAtom2++) {
			if(pAtom1->molecule == pAtom2->molecule)
				continue;
			r = pAtom1->Distance(*pAtom2);
			if(r < dist_min)
				dist_min = r;
		}
		
		// check if within separator distance
		if(dist_min >= settings.sep) {
			// too far away from other molecule, blocker atom only
			pAtom1->atten = ATTEN_BLOCKER;
			results.surface[pAtom1->molecule].nBlockedAtoms++;
		} else {
			// potential interface or neighbouring atom	
			pAtom1->atten = ATTEN_BURIED_FLAGGED;
			results.surface[pAtom1->molecule].nBuriedAtoms++;
		}
	}

	return 1;
}

////////////////////////////////////////////////////////////
// Molecular surface calculation

// Compute the surface for the atoms in the interface and its neighbours
// (routine mds)
int ScCalc::CalcDotsForAllAtoms(std::vector<PDBAtom> &atoms)
{
	// Calc maximum radius for atoms list
	radmax = 0.0;
	for(std::vector<PDBAtom>::const_iterator pAtom1 = atoms.begin(); pAtom1 < atoms.end(); pAtom1++) {
		if(pAtom1->radius > radmax)
			radmax = pAtom1->radius;
	}

	// Add dots for each atom in the list	
	for(std::vector<PDBAtom>::iterator pAtom1 = atoms.begin(); pAtom1 < atoms.end(); pAtom1++) {
		PDBAtom &atom1 = *pAtom1;
		if(atom1.atten <= 0)
			continue;

		// Find neighbor
		if(!FindNeighbordsAndBuriedAtoms(atom1))
			continue;

		if(!atom1.access)
			continue;
		if(atom1.atten <= ATTEN_BLOCKER)
			continue;
		if(atom1.atten == ATTEN_6 && atom1.buried.empty())
			continue;

		// Generate convex surface
		GenerateConvexSurface(atom1);
	}

	// Concave surface generation
	if(settings.rp > 0)
		GenerateConcaveSurface();
			
	return 1;
}

// Calculate surface dots around a single atom (main loop in original code)
int ScCalc::FindNeighbordsAndBuriedAtoms(PDBAtom &atom1)
{
	if(!FindNeighborsForAtom(atom1))
		return 0;
	
	// sort neighbors by distance from atom1
	_atom_distance_ref = &atom1; 
	std::sort(atom1.neighbors.begin(), atom1.neighbors.end(), _atom_distance_cb);
	_atom_distance_ref = NULL;

	SecondLoop(atom1);

	return atom1.neighbors.size();
}

// Make a list of neighboring atoms from atom1
// (loop to label 100)
int ScCalc::FindNeighborsForAtom(PDBAtom &atom1)
{
	std::vector<PDBAtom>::iterator iAtom2;
	std::vector<PDBAtom*> &neighbors = atom1.neighbors;
	float d2;
	float bridge;
	float bb2 = pow(4 * radmax + 4 * settings.rp, 2);
	int nbb = 0;

	for(iAtom2 = atoms.begin(); iAtom2 < atoms.end(); iAtom2++) {
		PDBAtom &atom2 = *iAtom2;
		if(atom1 == atom2 || atom2.atten <= 0)
			continue;

		if(atom1.molecule == atom2.molecule) {
			
			d2 = atom1.Distance2(atom2);

			if(d2 <= 0.0001)
				throw ScCalcException("Coincident atoms: %d == %d", atom1.natom, atom2.natom);
				
			bridge = atom1.radius + atom2.radius + 2 * settings.rp;
			if (d2 >= bridge * bridge)
				continue;

			neighbors.push_back(&atom2);
					
		} else {
			
			if(atom2.atten < ATTEN_BURIED_FLAGGED)
				continue;

			d2 = atom1.Distance2(atom2);
			if (d2 < bb2)
				nbb++;

			bridge = atom1.radius + atom2.radius + 2 * settings.rp;
			if (d2 >= bridge * bridge)
				continue;

			atom1.buried.push_back(&atom2);
		}
	}

	if(atom1.atten == ATTEN_6 && !nbb)
		return 0;

	if(neighbors.empty()) {
		// no neighbors
		atom1.access = 1;
		// no convex surface generation
		return 0;
	}

	return neighbors.size();
}

// second loop (per original code)
int ScCalc::SecondLoop(PDBAtom &atom1)
{
	CVec3 uij, tij;
	float erj, eri, rij, density, dij, asymm, far, contain;
	int between;
	std::vector<PDBAtom*> &neighbors = atom1.neighbors;

	eri = atom1.radius + settings.rp;
	
	for(std::vector<PDBAtom*>::iterator iAtom2 = neighbors.begin(); iAtom2 < neighbors.end(); iAtom2++) {
		PDBAtom &atom2 = **iAtom2;
		
		if(atom2 <= atom1)
			continue; 
			
		erj = atom2.radius + settings.rp;
		density = (atom1.density + atom2.density) / 2;
		dij = atom1.Distance(atom2);
		
		uij = (atom2 - atom1) / dij;
		asymm = (eri * eri - erj * erj) / dij;
		between = (ABS(asymm) < dij);

		tij = ((atom1 + atom2) * 0.5) + (uij * (asymm * 0.5));
		
		far = (eri + erj) * (eri + erj) - dij * dij;
		if (far <= 0.0)
			continue;	// return?
				
		far = sqrt(far);

		contain = dij * dij - ((atom1.radius - atom2.radius) * (atom1.radius - atom2.radius));
		if (contain <= 0.0)
			continue;
				
		contain = sqrt(contain);
		rij = 0.5 * far * contain / dij;

		if (neighbors.size() <= 1) {
			atom1.access = 1;
			atom2.access = 1;
			break;
		} 

		ThirdLoop(atom1, atom2, uij, tij, rij);

		if(atom1.atten > ATTEN_BLOCKER || (atom2.atten > ATTEN_BLOCKER && settings.rp > 0.0))
			GenerateToroidalSurface(atom1, atom2, uij, tij, rij, between);
	}

	return 1;
}

// third loop (per original code)
int ScCalc::ThirdLoop(
		PDBAtom& atom1,
		PDBAtom& atom2,
		const CVec3 &uij,
		const CVec3 &tij,
		float rij)
{
	std::vector<PDBAtom*> &neighbors = atom1.neighbors;
	float eri, erj, erk, djk, dik;
	float asymm, dt, dtijk2, hijk, isign, wijk, swijk, rkp2;
	int is0;
	CVec3 uik, dijk, pijk, utb, iujk, tv, tik, bijk, uijk;

	eri = atom1.radius + settings.rp;
	erj = atom2.radius + settings.rp;

	for(std::vector<PDBAtom*>::iterator iAtom3 = neighbors.begin();
			iAtom3 < neighbors.end(); iAtom3++) {
		PDBAtom &atom3 = **iAtom3;
		if(atom3 <= atom2)
			continue;

		erk = atom3.radius + settings.rp;
		djk = atom2.Distance(atom3);
		if(djk >= erj+erk)
			continue;

		dik = atom1.Distance(atom3);
		if(dik >= eri+erk)
			continue;

		if(atom1.atten <= ATTEN_BLOCKER && atom2.atten <= ATTEN_BLOCKER && atom3.atten <= ATTEN_BLOCKER)
			continue;

		uik = (atom3 - atom1) / dik;
		dt = uij.Dot(uik);
		wijk = acosf(dt);
		swijk = sinf(wijk);

		if(dt >= 1.0 || dt <= -1.0 || wijk <= 0.0 || swijk <= 0.0) {
			// collinear and other
			dtijk2 = tij.Distance(atom3);
			rkp2 = erk * erk - rij * rij;
			if(dtijk2 < rkp2)
				// 600
				return 0;
			continue;
		}

		uijk = uij.Cross(uik) / swijk;
		utb = uijk.Cross(uij);
		asymm = (eri*eri - erk*erk) / dik;
		tik = (atom1 + atom3)*0.5 + uik*asymm*0.5;
		tv = uik * (tik - tij);
		dt = tv.x + tv.y + tv.z;
		bijk = tij + utb * dt / swijk;
		hijk = eri*eri - bijk.Distance2(atom1);
		if(hijk <= 0.0)
			// no height, skip
			continue;

		hijk = sqrt(hijk);
		for(is0 = 1; is0 <= 2; is0++) {
			isign = 3 - 2 * is0;
			pijk = bijk + uijk * hijk * isign;

			// check for collision
			if(CheckAtomCollision2(pijk, atom2, atom3, neighbors))
				continue;

			// new probe position
			SCCALC_PROBE probe;
			if(isign > 0) {
				probe.pAtoms[0] = &atom1;
				probe.pAtoms[1] = &atom2;
				probe.pAtoms[2] = &atom3;
			} else {
				probe.pAtoms[0] = &atom2;
				probe.pAtoms[1] = &atom1;
				probe.pAtoms[2] = &atom3;
			}
			probe.height = hijk;
			probe.point = pijk;
			probe.alt = uijk * isign;
			probes.push_back(probe);

			atom1.access = 1;
			atom2.access = 1;
			atom3.access = 1;
		}
	}

	return 1;
}

// Check two atoms against a list of neighbors for collision
int ScCalc::CheckAtomCollision2(
		const CVec3 &pijk,
		const PDBAtom &atom1,
		const PDBAtom &atom2,
		const std::vector<PDBAtom*> &atoms)
{
	for(std::vector<PDBAtom*>::const_iterator ineighbor = atoms.begin();
			ineighbor < atoms.end(); ineighbor++) {
		const PDBAtom &neighbor = **ineighbor;
		if(&atom1 == &neighbor || &atom2 == &neighbor)
			continue;
		if(pijk.Distance2(neighbor) <= pow(neighbor.radius + settings.rp, 2))
			// collision detected
			return 1;
	}
	return 0;
}

// Generate convex surface for a specific atom
int ScCalc::GenerateConvexSurface(const PDBAtom &atom1)
{
	const std::vector<PDBAtom*> &neighbors = atom1.neighbors;
	const PDBAtom *neighbor;
	CVec3 north(0, 0, 1);
	CVec3 south(0, 0, -1);
	CVec3 eqvec(1, 0, 0);
	CVec3 vtemp, vql, uij, tij, pij, cen, pcen;
	float dt, erj, dij, eri, far, contain;
	float area, asymm, cs, ps, rad, ri, rij, rj;

	ri = atom1.radius;
	eri = (atom1.radius + settings.rp);

	if(!neighbors.empty()) {
		// use first neighbor
		neighbor = neighbors[0];

		north = atom1 - *neighbor;
		north.Normalize();

		vtemp.x = north.y*north.y + north.z*north.z;
		vtemp.y = north.x*north.x + north.z*north.z;
		vtemp.z = north.x*north.x + north.y*north.y;
		vtemp.Normalize();

		dt = vtemp.Dot(north);
		if(ABS(dt) > 0.99)
			vtemp = CVec3(1, 0, 0);

		eqvec = north.Cross(vtemp);
		eqvec.Normalize();
		vql = eqvec.Cross(north);

		rj = neighbor->radius;
		erj = neighbor->radius + settings.rp;
		dij = atom1.Distance(*neighbor);
		uij = (*neighbor - atom1) / dij;

		asymm = (eri*eri - erj*erj) / dij;
		tij = ((atom1 + *neighbor) * 0.5) + (uij * (asymm * 0.5));
		far = pow(eri + erj, 2) - dij*dij;
		if(far <= 0.0)
			throw ScCalcException("Imaginary far for atom %d, neighbor atom %d", atom1.natom, neighbor->natom);
		far = sqrt(far);

		contain = pow(dij, 2) - pow(ri - rj, 2);
		if(contain <= 0.0)
			throw ScCalcException("Imaginary contain for atom %d, neighbor atom %d", atom1.natom, neighbor->natom);
		contain = sqrt(contain);
		rij = 0.5 * far * contain / dij;
		pij = tij + (vql * rij);
		south = (pij - atom1) / eri;

		if(north.Cross(south).Dot(eqvec) <= 0.0)
			throw ScCalcException("Non-positive frame for atom %d, neighbor atom %d", atom1.natom, neighbor->natom);
	}

	// Generate subdivided arc
	std::vector<CVec3> lats;
	CVec3 o(0, 0, 0);
	cs = SubArc(o, ri, eqvec, atom1.density, north, south, lats);

	if(lats.empty())
		return 0;		
	
	// Project onto north vector
	std::vector<CVec3> points;
	for(vector<CVec3>::iterator ilat = lats.begin(); ilat < lats.end(); ilat++) {
		dt = ilat->Dot(north);
		cen = atom1 + (north*dt);
		rad = ri*ri - dt*dt;
		if(rad <= 0.0)
			continue;
		rad = sqrt(rad);

		points.clear();
		ps = SubCir(cen, rad, north, atom1.density, points);
		if(points.empty())
			continue;
		area = ps * cs;

		for(vector<CVec3>::iterator point = points.begin();
				point < points.end(); point++) {

			pcen = atom1 + ((*point - atom1) * (eri/ri));

			// Check for collision
			if(CheckPointCollision(pcen, neighbors))
				continue;

			// No collision, put point
			results.dots.convex++;
			AddDot(atom1.molecule, 1, *point, area, pcen, atom1);
		}
	}
	return 1;
}

// Check a point for collision against a list of atoms
int ScCalc::CheckPointCollision(
		const CVec3 &pcen,
		const std::vector<PDBAtom*> &atoms)
{
	for(std::vector<PDBAtom*>::const_iterator ineighbor = atoms.begin()+1;
		ineighbor < atoms.end(); ineighbor++) {
		if(pcen.Distance(**ineighbor) <= ((*ineighbor)->radius + settings.rp))
			// collision detected
			return 1;
	}
	return 0;
}

// Generate toroidal surface between two atoms
int ScCalc::GenerateToroidalSurface(
		PDBAtom &atom1,
		PDBAtom &atom2,
		const CVec3 uij,
		const CVec3 tij,
		float rij,
		int between)
{
	std::vector<PDBAtom*> &neighbors = atom1.neighbors;
	float density, ri, rj, rb, rci, rcj, rs, e, edens, eri, erj, erl, dtq, pcusp, anglei, anglej, dt, ts, ps, area;
	CVec3 pi, pj, axis, dij, pqi, pqj, qij, qjk, qj;
	
	vector<CVec3> subs;
	
	// following Fortran original
	// will be optimized by compiler
	density = (atom1.density + atom2.density) / 2;
	ri = atom1.radius;
	rj = atom2.radius;
	eri = (atom1.radius + settings.rp);
	erj = (atom2.radius + settings.rp);
	rci = rij * atom1.radius / eri;		
	rcj = rij * atom2.radius / erj;
	rb = rij - settings.rp;
	
	if(rb <= 0.0)
		rb = 0.0;
		
	rs = (rci + 2 * rb + rcj) / 4;
	e = rs / rij;
	edens = e * e * density;

	ts = SubCir(tij, rij, uij, edens, subs);
	if(subs.empty())
		return 0;

	for(vector<CVec3>::iterator iSub = subs.begin(); iSub < subs.end(); iSub++) {
		CVec3 &sub = *iSub;

		// check for collision
		int tooclose = 0;
		float d2 = 0;
		for(std::vector<PDBAtom*>::iterator ineighbor = neighbors.begin();
				!tooclose && ineighbor < neighbors.end() && !tooclose;
				ineighbor++) {
			const PDBAtom &neighbor = **ineighbor;	// for readability
			if(atom2 == neighbor)
				continue;
			erl = neighbor.radius + settings.rp;
			d2 = sub.Distance2(neighbor);
			tooclose = d2 < (erl * erl);
		}
		if(tooclose)
			continue;

		// no collision, toroidal arc generation
		CVec3 &pij = sub;
		atom1.access = 1;
		atom2.access = 1;
	
		if(atom1.atten == ATTEN_6 && atom2.atten == ATTEN_6&& atom1.buried.empty())
			continue;
			
		pi = (atom1 - pij) / eri;
		pj = (atom2 - pij) / erj;
		axis = pi.Cross(pj);
		axis.Normalize();
		
		dtq = pow(settings.rp, 2) - pow(rij, 2);
		pcusp = dtq > 0 && between;
		if(pcusp) {
			// point cusp -- two shortened arcs
			dtq = sqrt(dtq);
			qij = tij - uij * dtq;
			qjk = tij + uij * dtq;
			pqi = (qij - pij) / settings.rp;
			pqj = CVec3(0, 0, 0);

		} else {
			// no cusp
			pqi = pi + pj;		
			pqi.Normalize();
			pqj = pqi;
		}
		
		dt = pqi.Dot(pi);
		if(dt >= 1.0 || dt <= -1.0)
			return 0;
		anglei = acosf(dt);

		dt = pqj.Dot(pj);
		if(dt >= 1.0 || dt <= -1.0)
			return 0;
		anglej = acosf(dt);
		
		// convert two arcs to points
		if(atom1.atten >= ATTEN_2) {
			vector<CVec3> points;
			ps = SubArc(pij, settings.rp, axis, density, pi, pqi, points);
			for(std::vector<CVec3>::iterator point = points.begin(); point < points.end(); point++) {
				area = ps * ts * DistancePointToLine(tij, uij, *point) / rij;
				results.dots.toroidal++;
				AddDot(atom1.molecule, 2, *point, area, pij, atom1);
			}
		}

		if(atom2.atten >= ATTEN_2) {
			vector<CVec3> points;
			ps = SubArc(pij, settings.rp, axis, density, pqj, pj, points);
			for(std::vector<CVec3>::iterator point = points.begin(); point < points.end(); point++) {
				area = ps * ts * DistancePointToLine(tij, uij, *point) / rij;
				results.dots.toroidal++;
				AddDot(atom1.molecule, 2, *point, area, pij, atom2);
			}
		}
	}
	return 1;
}

// Generate concave surface for all probes
int ScCalc::GenerateConcaveSurface()
{
	vector<const SCCALC_PROBE*> lowprobs, nears;

	// collect low probes
	for(vector<SCCALC_PROBE>::iterator probe = probes.begin();
			probe < probes.end(); probe++) {
		if(probe->height < settings.rp)
			lowprobs.push_back(&(*probe));
	}

	for(vector<SCCALC_PROBE>::iterator probe = probes.begin();
			probe < probes.end(); probe++) {

		if(	probe->pAtoms[0]->atten == ATTEN_6 &&
			probe->pAtoms[1]->atten == ATTEN_6 &&
			probe->pAtoms[2]->atten == ATTEN_6) {
			continue;
		}

		CVec3 &pijk = probe->point, &uijk = probe->alt;
		float hijk = probe->height;
		float density = (
				probe->pAtoms[0]->density +
				probe->pAtoms[1]->density +
				probe->pAtoms[2]->density ) / 3;

		// gather nearby low probes
		nears.clear();
		for(vector<const SCCALC_PROBE*>::const_iterator lprobe = lowprobs.begin();
				lprobe < lowprobs.end(); lprobe++) {
			if(&(*probe) == *lprobe)
				continue;

			float d2 = pijk.Distance2((*lprobe)->point);
			if(d2 > 4 * pow(settings.rp, 2))
				continue;

			nears.push_back(*lprobe);
		}

		// set up vectors from probe center to three atoms
		CVec3 vp[3], vectors[3];
		for(int i = 0; i < 3; i++) {
			vp[i] = *(probe->pAtoms[i]) - pijk;
			vp[i].Normalize();
		}

		// set up vectors to three cutting planes
		vectors[0] = vp[0].Cross(vp[1]);
		vectors[1] = vp[1].Cross(vp[2]);
		vectors[2] = vp[2].Cross(vp[0]);
		vectors[0].Normalize();
		vectors[1].Normalize();
		vectors[2].Normalize();

		// find latitude of highest vertex of triangle
		float dm = -1.0;
		int mm = 0;
		for(int i = 0; i < 3; i++) {
			float dt = uijk.Dot(vp[i]);
			if(dt > dm) {
				dm = dt;
				mm = i;
			}
		}

		// create arc for selecting latitudes
		CVec3 south = -uijk;
		CVec3 axis = vp[mm].Cross(south);
		axis.Normalize();

		std::vector<CVec3> lats;
		CVec3 o(0, 0, 0);
		float cs;

		cs = SubArc(o, settings.rp, axis, density, vp[mm], south, lats);
		if(lats.empty())
			continue;

		std::vector<CVec3> points;
		for(vector<CVec3>::iterator ilat = lats.begin();
				ilat < lats.end(); ilat++) {
			float dt, area, rad, ps;
			CVec3 cen;

			dt = ilat->Dot(south);
			cen = south * dt;
			rad = pow(settings.rp, 2) - pow(dt, 2);
			if(rad <= 0.0)
				continue;
			rad = sqrtf(rad);

			points.clear();
			ps = SubCir(cen, rad, south, density, points);
			if(points.empty())
				continue;

			area = ps * cs;

			for(vector<CVec3>::iterator point = points.begin();
					point < points.end(); point++) {
				// check against 3 planes
				int bail = 0;
				for(int i = 0; i < 3; i++) {
					float dt = point->Dot(vectors[i]);
					if(dt >= 0.0) {
						bail = 1;
						break;
					}
				}
				if(bail)
					continue;

				*point += pijk;

				if((hijk < settings.rp && !nears.empty()) &&
					CheckProbeCollision(*point, nears, pow(settings.rp, 2)))
						continue;

				// determine which atom the surface point is closest to
				int mc = 0;
				float dmin = 2 * settings.rp;
				for(int i = 0; i < 3; i++) {
					float d = point->Distance(*(probe->pAtoms[i])) -
							probe->pAtoms[i]->radius;
					if(d < dmin) {
						dmin = d;
						mc = i;
					}
				}

				// No collision, put point
				results.dots.concave++;
				AddDot(probe->pAtoms[mc]->molecule, 3, *point, area, pijk, *probe->pAtoms[mc]);
			}
		}
	}
	return 1;
}

// Check a point against a set of probes for collision within radius^2
int ScCalc::CheckProbeCollision(
		const CVec3 &point,
		const std::vector<const SCCALC_PROBE*> nears,
		const float r2)
{
	for(vector<const SCCALC_PROBE*>::const_iterator near = nears.begin();
			near < nears.end(); near++) {
		if(point.Distance2((*near)->point) < r2)
			// Collision
			return 1;
	}
	return 0;
}

// Add a molecular dot
void ScCalc::AddDot(
		const int molecule,
		const int type,
		const CVec3 coor,
		float area,
		const CVec3 pcen,
		const PDBAtom &atom)
{
	SCCALC_DOT dot = { coor, 0, type, area, CVec3(), &atom };
	float pradius = settings.rp, erl;

	// calculate outward pointing unit normal vector
	if(pradius <= 0)
		dot.outnml = coor - atom;
	else
		dot.outnml = (pcen - coor) / pradius;

	// determine whether buried

#if not MULTI_THREADED
	// This shortcut was in the original Fortran source, but
	// it will not work in a multi-threaded application.
	static CVec3 pprev;
	static int pburied =0;

	// first check whether probe changed
	if(pcen.Distance2(pprev) <= 0.0) {
		dot.buried = pburied;
	} else {
#endif

		// check for collision with neighbors in other molecules
		dot.buried = 0;
		for(std::vector<PDBAtom*>::const_iterator iNeighbor = atom.buried.begin();
				iNeighbor < atom.buried.end();
				iNeighbor++) {
			erl = (*iNeighbor)->radius + pradius;
			float d = pcen.Distance2(**iNeighbor);
			if(d <= erl*erl) {
				dot.buried = 1;
				break;
			}

		}

#if not MULTI_THREADED
		pprev = pcen;
		pburied = dot.buried;
	}
#endif

/*
	printf("#%4d (%d): %8.4f, %8.4f, %8.4f=  %8.4f  =  %8.4f, %8.4f, %8.4f :: flags:%04X\n",
			atom.natom, type, coor.x, coor.y, coor.z, area,
			dot.outnml.x, dot.outnml.y, dot.outnml.z, dot.buried);
*/

	dots[molecule].push_back(dot);
}


////////////////////////////////////////////////////////////
// Elementary functions

// Calculate distance from point to line
float ScCalc::DistancePointToLine(
		const CVec3 &cen,
		const CVec3 &axis,
		const CVec3 &pnt)
{
	CVec3 vec = pnt - cen;
	float dt = vec.Dot(axis);
	float d2 = vec.Magnitude2() - pow(dt, 2);
	return d2 < 0.0 ? 0.0 : sqrt(d2);
}

// Generate sub arc of molecular dots centered around a defined point
float ScCalc::SubArc(
		const CVec3 &cen,
		const float rad,
		const CVec3 &axis,
		const float density,
		const CVec3 &x,
		const CVec3 &v,
		std::vector<CVec3> &points)
{
	CVec3 y;
	float angle;
	float dt1, dt2;

	y = axis.Cross(x);
	dt1 = v.Dot(x);
	dt2 = v.Dot(y);
	angle = atan2(dt2, dt1);

	if(angle < 0.0)
		angle = angle + 2*PI;

	return SubDiv(cen, rad, x, y, angle, density, points);
}

// Subdivide defined arc and generate molecular dots
float ScCalc::SubDiv(
		const CVec3 &cen,
		const float rad,
		const CVec3 &x,
		const CVec3 &y,
		const float angle,
		const float density,
		std::vector<CVec3> &points)
{
	float delta, a, c, s, ps;
	int i;

	delta = 1.0 / (sqrt(density) * rad);
	a = - delta / 2;

	for(i = 0; i < MAX_SUBDIV; i++) {
		a = a + delta;
		if(a > angle)
			break;

		c = rad * cosf(a);
		s = rad * sinf(a);
		points.push_back(CVec3(cen + x*c + y*s));
	}

	if (a + delta < angle)
		throw ScCalcException("Too many subdivisions");

	if (!points.empty())
		ps = rad * angle / points.size();
	else
		ps = 0.0;

	return ps;
}

// Generate an arbitrary unit vector perpendicular to axis
float ScCalc::SubCir(
		const CVec3 &cen,
		const float rad,
		const CVec3 &axis,
		const float density,
		std::vector<CVec3> &points)
{
	CVec3 v1, v2, x, y;
	float dt;

	v1.x = pow(axis.y, 2) + pow(axis.z, 2);
	v1.y = pow(axis.x, 2) + pow(axis.z, 2);
	v1.z = pow(axis.x, 2) + pow(axis.y, 2);
	v1.Normalize();
	dt = v1.Dot(axis);

	if(ABS(dt) > 0.99) {
		v1.x = 1.0;
		v1.y = 0.0;
		v1.z = 0.0;
	}

	v2 = axis.Cross(v1);
	v2.Normalize();
	x = axis.Cross(v2);
	x.Normalize();
	y = axis.Cross(x);

	return SubDiv(cen, rad, x, y, 2.0 * PI, density, points);
}

////////////////////////////////////////////////////////////
// SC molecular dot trimming, vector dot product calculation and statistics

// Trim dots and retain only the peripheral band
float ScCalc::TrimPeripheralBand(
		const std::vector<SCCALC_DOT> &sdots,
		std::vector<const SCCALC_DOT*> &trimmed_dots)
{
	float area = 0;

#ifdef GPU
	if(settings.gpu.ready) {
		area = gpuTrimPeripheralBand(sdots, trimmed_dots);
	} else {
#endif

	// Loop over one surface
	// If a point is buried then see if there is an accessible point within distance band

	for(std::vector<SCCALC_DOT>::const_iterator iDot = sdots.begin(); iDot < sdots.end(); iDot++) {
		const SCCALC_DOT &dot = *iDot;
		// Paralelleizable kernel function
		if(dot.buried && TrimPeripheralBandCheckDot(dot, sdots)) {
			area += dot.area;
			trimmed_dots.push_back(&dot);
		}
	}

#ifdef GPU
	}
#endif

	return area;
}

// Test a dot against a set of dots for collision
// NOTE: ~75% of time is spent in this function
int ScCalc::TrimPeripheralBandCheckDot(
		const SCCALC_DOT &dot,
		const std::vector<SCCALC_DOT> &sdots)
{
	// Caching of r2 only brings 0.5% speed boost
	float r2 = pow(settings.band, 2);

	for(std::vector<SCCALC_DOT>::const_iterator iDot2 = sdots.begin(); iDot2 < sdots.end(); iDot2++) {
		const SCCALC_DOT &dot2 = *iDot2;
		if(&dot == &dot2)
			continue;
		if(dot2.buried)
			continue;
		if(dot.coor.Distance2(dot2.coor) <= r2)
	 		return 0;
	}
	return 1;
}

////////////////////////////////////////////////////////////

// Calculate separation distance and and normal vector dot product (shape)
// distributions, mean and median

int ScCalc::CalcNeighborDistance(
		const int molecule,
		const std::vector<const SCCALC_DOT*> &my_dots,
		const std::vector<const SCCALC_DOT*> &their_dots)
{
	map<int,int> dbins;	// Distance bins
	map<int,int> sbins;	// Vector dot product bins (sc)
	float norm_sum = 0.0, distmin_sum = 0.0;
	int ibin;
	float total = 0.0;

        if(my_dots.empty() || their_dots.empty())
                return 0;

	for(std::vector<const SCCALC_DOT*>::const_iterator iDot = my_dots.begin();
			iDot < my_dots.end(); iDot++) {
		//if((*iDot)->buried)
		//	results.surface[molecule].nBuriedDots++;
		//else
		//	results.surface[molecule].nAccessibleDots++;
		total += (*iDot)->area;
	}

#ifdef GPU
	std::vector<const SCCALC_DOT*> neighbors;
	std::vector<const SCCALC_DOT*>::const_iterator iNeighbor;

	if(settings.gpu.ready) {
		gpuFindClosestNeighbors(my_dots, their_dots, neighbors);
		iNeighbor = neighbors.begin();
	}
#endif

	for(std::vector<const SCCALC_DOT*>::const_iterator iDot = my_dots.begin();
			iDot < my_dots.end(); iDot++) {
		const SCCALC_DOT &dot1 = **iDot;

		float distmin, r;
		const SCCALC_DOT *neighbor = NULL;

#ifdef GPU
		if(settings.gpu.ready)
			neighbor = *iNeighbor++;
		else
#endif
		neighbor = CalcNeighborDistanceFindClosestNeighbor(dot1, their_dots);

		if(!neighbor)
			continue;

		// having looked at all possible neighbours now accumulate stats
		distmin = neighbor->coor.Distance(dot1.coor);
		distmin_sum += distmin;
		// decide which bin to put it into and then add to distance histogram
		ibin = (int)(distmin / settings.binwidth_dist);
		dbins[ibin]++;

		//work out dot product
		r = dot1.outnml.Dot(neighbor->outnml);

		// weight dot product
		// cpjx I think the weighting factor is the denominator 2 below?
		// cpjx		  r = r * exp( - (distmin**2) / 2.)
		r = r * exp( - pow(distmin, 2) * settings.weight );
		// rounding errors a problem, so ensure abs(r) <1
		r = MIN(0.999, MAX(r, -0.999));
		norm_sum += r;

		// left_trunc float to int ibin
		// otherwise: (int)-0.9 = 0.
		r /= settings.binwidth_norm;
		if(r >= 0)
			ibin = (int)r;
		else
			ibin = (int)r -1;
		sbins[ibin]++;
	}

	// Determine the last distance bin that has anything in it
	// Accumulate percentages and area from all filled distance bins
	float abin, cumarea =0, cumperc = 0, perc, c;
	float rleft =0, rmedian =0;
	map<int,int>::const_iterator it;

	VERBOSE("\nDistance between surfaces D(" << (molecule+1) << "->" << (molecule+1)%2+1 << "):" << std::endl);
	VERBOSE("From - To\tArea\tCum. Area\t%\tCum. %\n");

	for(it = dbins.begin(); it != dbins.end(); it++) {
		abin = total * (it->second) / my_dots.size();
		cumarea += abin;
		perc = abin * 100 / total;
		c = cumperc + perc;
		if(cumperc <= 50 && c >= 50) {
			rleft = (it->first) * settings.binwidth_dist;
			rmedian = rleft + (50 - cumperc) * settings.binwidth_dist / ( c - cumperc );
		}
		cumperc = c;

		if(verbose) {
			char buf[128];
			snprintf(buf, sizeof(buf),
				"%.2f - %.2f\t%.1f\t%.1f\t%.1f\t%.1f",
				(float)it->first * settings.binwidth_dist,
				(float)it->first * settings.binwidth_dist + settings.binwidth_dist,
				abin, cumarea,
				perc, cumperc);
			*verbose << buf << std::endl;
		}
	}

	results.surface[molecule].d_mean = distmin_sum / my_dots.size();
	results.surface[molecule].d_median = rmedian;

	VERBOSE("\nSurface complementarity S(" << (molecule+1) << "->" << (molecule+1)%2+1 << "):" << std::endl);
	VERBOSE("From - To\tNumber\t%\tCumm. %\n");

	cumperc = 0;
	for(it = sbins.begin(); it != sbins.end(); it++) {
		perc = (float)(it->second) * 100 / my_dots.size();
		c = cumperc + perc;
		if(cumperc <= 50 && c >= 50) {
			rleft = (float)(it->first) * settings.binwidth_norm;
			rmedian = rleft + (50 - cumperc) * settings.binwidth_norm / ( c - cumperc );
		}
		cumperc = c;

		if(verbose) {
			char buf[128];
			snprintf(buf, sizeof(buf),
				"%.2f - %.2f\t%d\t%.1f\t%.1f",
				(float)-it->first * settings.binwidth_norm - settings.binwidth_norm,
				(float)-it->first * settings.binwidth_norm,
				it->second, perc, cumperc);
			*verbose << buf << std::endl;
		}
	}

	results.surface[molecule].s_mean= -norm_sum / my_dots.size();
	results.surface[molecule].s_median = -rmedian;

	return 1;
}

// Find closest neighbor dot for a given dot
// NOTE: ~20% of time is spent in this function
const SCCALC_DOT *ScCalc::CalcNeighborDistanceFindClosestNeighbor(
		const SCCALC_DOT &dot1,
		const std::vector<const SCCALC_DOT*> &their_dots)
{
	float distmin = 999999.0, d;
	const SCCALC_DOT *neighbor = NULL;

	// Loop over the entire surface: find and flag neighbour of each point
	// that we're interested in and store nearest neighbour pointer

	for(std::vector<const SCCALC_DOT*>::const_iterator iDot2 = their_dots.begin();
			iDot2 < their_dots.end(); iDot2++) {
		const SCCALC_DOT &dot2 = **iDot2;
		if(!dot2.buried)
			continue;
		d = dot2.coor.Distance2(dot1.coor);
		if(d <= distmin) {
			distmin = d;
			neighbor = &dot2;
		}
	}
	return neighbor;
}

#ifdef GPU
#define gpuAssert(err) gpuThrowException(err, __FILE__, __LINE__)
#define UPPER_MULTIPLE(n,d) (((n)%(d)) ? (((n)/(d)+1)*(d)) : (n))
#endif

float inline GetTimerMs(clock_t &start, int reset)
{
	clock_t now = clock();
	float d = (now - start)/(CLOCKS_PER_SEC/1000);
	if(reset)
		start = now;
	return d;
}

#ifdef CUDA_GPU

void ScCalc::cudaThrowException(cudaError_t err, const char *fn, int line)
{
	if (cudaSuccess != err)
		throw ScCalcException("CUDA Exception at %s:%d: %s", fn, line, cudaGetErrorString(err));
}

void ScCalc::cudaInit()
{
	if(settings.gpu.ready || !settings.gpu.device)
		return;

	cudaDeviceProp deviceProp;

	if(settings.gpu.device < 0) {
		// Detect GPU and available threads
		int dev, deviceCount;
		if (cudaGetDeviceCount(&deviceCount) == cudaSuccess) {
			for (dev = 0; dev < deviceCount; ++dev) {
				cudaGetDeviceProperties(&deviceProp, dev);
				if(deviceProp.maxThreadsPerBlock*deviceProp.multiProcessorCount > 
					settings.gpu.threads*settings.gpu.proc) {
					settings.gpu.proc = deviceProp.multiProcessorCount;
					settings.gpu.threads = settings.gpu.threads ? 
						MIN(settings.gpu.threads, deviceProp.maxThreadsPerBlock) : 
						deviceProp.maxThreadsPerBlock;
					settings.gpu.device = dev + 1;
				}
			}
		}
	} else if(settings.gpu.device > 0) {
	   	gpuAssert( cudaGetDeviceProperties(&deviceProp, settings.gpu.device -1) );
		settings.gpu.proc = deviceProp.multiProcessorCount;
		settings.gpu.threads = settings.gpu.threads ? 
			MIN(settings.gpu.threads, deviceProp.maxThreadsPerBlock) : 
			deviceProp.maxThreadsPerBlock;
	}

	if(!settings.gpu.threads) {
		if(settings.gpu.device < 0 && *error)
			*error << "No GPU available. GPU acceleration disabled." << std::endl;
	} else {
		VERBOSE("GPU support enabled: " << deviceProp.name << 
				" [" << (deviceProp.clockRate/1000) << " MHz, capability " << 
				deviceProp.major << "." << deviceProp.minor << "] with " <<
				deviceProp.multiProcessorCount << " processors, " <<
				settings.gpu.threads <<" threads." <<
				std::endl);
		settings.gpu.ready = 1;
	}
}

float ScCalc::cudaTrimPeripheralBand(
		const std::vector<SCCALC_DOT> &dots,
		std::vector<const SCCALC_DOT*> &trimmed_dots)
{
	int n, nBur, nAcc;
	int threads;
	float area = 0;
	clock_t timer = clock();

	threads = MIN(1024, settings.gpu.threads);
	n = dots.size();

	// Host and device (GPU) memory pointers for dot coordinates and results
	float3 *hAccDotCoords, *phAccDotCoords;
	float3 *hBurDotCoords, *phBurDotCoords;
	float3 *dAccDotCoords;
	float3 *dBurDotCoords;
	char *hDotColl;
	char *dDotColl;

	cudaEvent_t start, stop;

	// Allocate host memory
	// Use cudaHostAlloc for DMA zero-copy
	gpuAssert( cudaHostAlloc((void**)&hAccDotCoords, n * sizeof(*hAccDotCoords), cudaHostAllocDefault) );
	gpuAssert( cudaHostAlloc((void**)&hBurDotCoords, n * sizeof(*hBurDotCoords), cudaHostAllocDefault) );

	// Make CUDA copy of (x, y, z) buried and accessible coordinates
	phAccDotCoords = hAccDotCoords;
	phBurDotCoords = hBurDotCoords;
	for(std::vector<SCCALC_DOT>::const_iterator iDot = dots.begin();
			iDot < dots.end(); iDot++) {
		// Quick copy
		if(iDot->buried)
			*phBurDotCoords++ = *((float3*)&iDot->coor);
		else	
			*phAccDotCoords++ = *((float3*)&iDot->coor);
	}
	nBur = phBurDotCoords - hBurDotCoords;
	nAcc = phAccDotCoords - hAccDotCoords;

	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	// Allocate host memory for results (detected collisions)
	gpuAssert( cudaHostAlloc((void**)&hDotColl, nBur * sizeof(*hDotColl), cudaHostAllocDefault) );
	for(int i =0; i < nBur; i++) hDotColl[i] = 0;

	// Allocate GPU memory
	gpuAssert( cudaMalloc((void **)&dBurDotCoords, UPPER_MULTIPLE(nBur, threads) * sizeof(*dBurDotCoords)) );
	gpuAssert( cudaMalloc((void **)&dAccDotCoords, nAcc * sizeof(*dAccDotCoords)) );
	gpuAssert( cudaMalloc((void **)&dDotColl, UPPER_MULTIPLE(nBur, threads) * sizeof(*dDotColl)) );

	// Copy data from host to GPU
	gpuAssert( cudaMemcpy(dBurDotCoords, hBurDotCoords, nBur * sizeof(*dBurDotCoords), cudaMemcpyHostToDevice) );
	gpuAssert( cudaMemcpy(dAccDotCoords, hAccDotCoords, nAcc * sizeof(*dAccDotCoords), cudaMemcpyHostToDevice) );

	for(int i =0; i < nBur; i++) hDotColl[i] = -1;

	// Run kernel in multi-threaded blocks on GPU
	cudaEventRecord( start, 0 );
	_cuda_TrimPeripheralBand(UPPER_MULTIPLE(nBur, threads)/threads, threads,
				dAccDotCoords, nAcc,
				dBurDotCoords,
				dDotColl, pow(settings.band, 2));
	cudaEventRecord( stop, 0 );

	// Wait for threads to finish and copy back memory from GPU to host
	gpuAssert( cudaThreadSynchronize() );
	gpuAssert( cudaMemcpy(hDotColl, dDotColl,
			nBur * sizeof(*dDotColl), cudaMemcpyDeviceToHost) );

	// Make a new list of dots that have no collisions
	char *p = hDotColl;
	for(std::vector<SCCALC_DOT>::const_iterator iDot = dots.begin();
			iDot < dots.end(); iDot++) {
		const SCCALC_DOT &dot1 = *iDot;
		if(!iDot->buried)
			continue;	
		if(!*p++) {
			area += dot1.area;
			trimmed_dots.push_back(&dot1);
		}
	}	
	
	// Free GPU and host memory
	gpuAssert( cudaFree(dBurDotCoords) );
	gpuAssert( cudaFree(dAccDotCoords) );
	gpuAssert( cudaFree(dDotColl) );

	gpuAssert( cudaFreeHost(hBurDotCoords) );
	gpuAssert( cudaFreeHost(hAccDotCoords) );
	gpuAssert( cudaFreeHost(hDotColl) );

	float time_kernel;
	float time_total = GetTimerMs(timer, 0);
	cudaEventElapsedTime( &time_kernel, start, stop );
	VERBOSE("Peripheral trimming GPU processing time: " << time_kernel << " ms kernel, " << time_total << " ms total" << std::endl);

	gpuAssert( cudaEventDestroy(start) );
	gpuAssert( cudaEventDestroy(stop) );

	return area;
}

int ScCalc::cudaFindClosestNeighbors(
	const std::vector<const SCCALC_DOT*> &my_dots,
	const std::vector<const SCCALC_DOT*> &their_dots,
	std::vector<const SCCALC_DOT*> &neighbors)
{
	int nMyDots, nTheirDots, nNeighbors;
	int threads;
	clock_t timer = clock();

	threads = MIN(512, settings.gpu.threads);

	// Memory pointers for my and their dot coordinate arrays, CPU and GPU
	float3 *hMyDotCoords, *phMyDotCoords;
	float3 *hTheirDotCoords, *phTheirDotCoords;
	float3 *dMyDotCoords, *dTheirDotCoords;

	// Dot point pointer map
	const SCCALC_DOT **hTheirDots, **phTheirDots;

	// Neighbor ID memory pointers
	uint *hNeighbors;
	uint *dNeighbors;

	// Timers
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	nMyDots = my_dots.size();
	nTheirDots = their_dots.size();
	nNeighbors = nMyDots;

	// Allocate host memory
	gpuAssert( cudaHostAlloc((void**)&hMyDotCoords, nMyDots * sizeof(*hMyDotCoords), cudaHostAllocDefault) );
	gpuAssert( cudaHostAlloc((void**)&hTheirDotCoords, nTheirDots * sizeof(*hTheirDotCoords), cudaHostAllocDefault) );
	gpuAssert( cudaHostAlloc((void**)&hTheirDots, nTheirDots * sizeof(*hTheirDots), cudaHostAllocDefault) );
	gpuAssert( cudaHostAlloc((void**)&hNeighbors, nNeighbors * sizeof(*hNeighbors), cudaHostAllocDefault) );

	// Make CUDA copy of (x, y, z) dot coordinates for my dots
	phMyDotCoords = hMyDotCoords;
	for(std::vector<const SCCALC_DOT*>::const_iterator iDot = my_dots.begin(); iDot < my_dots.end(); iDot++) {
		// Quick copy
		*phMyDotCoords++ = *((float3*)&(*iDot)->coor);
	}
	nMyDots = phMyDotCoords - hMyDotCoords;

	// Make CUDA copy of (x, y, z) dot coordinates for their dots and keep a map
	phTheirDotCoords = hTheirDotCoords;
	phTheirDots = hTheirDots;
	for(std::vector<const SCCALC_DOT*>::const_iterator iDot = their_dots.begin(); iDot < their_dots.end(); iDot++) {
		// Quick copy
		if(!(*iDot)->buried)
			continue;
		*phTheirDotCoords++ = *((float3*)&(*iDot)->coor);
		*phTheirDots++ = *iDot;
	}
	nTheirDots = phTheirDotCoords - hTheirDotCoords;

	// Allocate GPU memory and copy data there
	gpuAssert( cudaMalloc((void **)&dMyDotCoords, UPPER_MULTIPLE(nMyDots, threads) * sizeof(*dMyDotCoords)) );
	gpuAssert( cudaMalloc((void **)&dTheirDotCoords, nTheirDots * sizeof(*dTheirDotCoords)) );
	gpuAssert( cudaMalloc((void **)&dNeighbors, UPPER_MULTIPLE(nNeighbors, threads) * sizeof(*dNeighbors)) );
	gpuAssert( cudaMemcpy(dMyDotCoords, hMyDotCoords, nMyDots * sizeof(*dMyDotCoords), cudaMemcpyHostToDevice) );
	gpuAssert( cudaMemcpy(dTheirDotCoords, hTheirDotCoords, nTheirDots * sizeof(*dTheirDotCoords), cudaMemcpyHostToDevice) );

	// Run kernel in multi-threaded blocks on GPU
	cudaEventRecord( start, 0 );
	_cuda_FindClosestNeighbor(UPPER_MULTIPLE(nMyDots, threads)/threads, threads,
				dMyDotCoords,
				dTheirDotCoords, nTheirDots,
				dNeighbors);
	cudaEventRecord( stop, 0 );

	// Wait for threads to finish and copy back memory from GPU to host
	gpuAssert( cudaThreadSynchronize() );
	gpuAssert( cudaMemcpy(hNeighbors, dNeighbors, nMyDots * sizeof(*hNeighbors), cudaMemcpyDeviceToHost) );

	for(int i = 0; i < nNeighbors; i++)
		neighbors.push_back( hTheirDots[hNeighbors[i]] );

	// Free memory
	gpuAssert( cudaFree(dMyDotCoords) );
	gpuAssert( cudaFree(dTheirDotCoords) );
	gpuAssert( cudaFree(dNeighbors) );

	gpuAssert( cudaFreeHost(hMyDotCoords) );
	gpuAssert( cudaFreeHost(hTheirDotCoords) );
	gpuAssert( cudaFreeHost(hTheirDots) );
	gpuAssert( cudaFreeHost(hNeighbors) );

	float time_kernel;
	float time_total = GetTimerMs(timer, 0);
	cudaEventElapsedTime( &time_kernel, start, stop );
	VERBOSE("Find Neighbors GPU processing time: " << time_kernel << " ms kernel, " << time_total << " ms total" << std::endl);

	gpuAssert( cudaEventDestroy(start) );
	gpuAssert( cudaEventDestroy(stop) );

	return 1;
}
#endif

#ifdef OPENCL_GPU
typedef cl_float4  float4;

const char *errstr(int errorCode) {
	switch(errorCode) {
		case CL_DEVICE_NOT_FOUND:
			return "CL_DEVICE_NOT_FOUND";
		case CL_DEVICE_NOT_AVAILABLE:
			return "CL_DEVICE_NOT_AVAILABLE";
		case CL_COMPILER_NOT_AVAILABLE:
			return "CL_COMPILER_NOT_AVAILABLE";
		case CL_MEM_OBJECT_ALLOCATION_FAILURE:
			return "CL_MEM_OBJECT_ALLOCATION_FAILURE";
		case CL_OUT_OF_RESOURCES:
			return "CL_OUT_OF_RESOURCES";
		case CL_OUT_OF_HOST_MEMORY:
			return "CL_OUT_OF_HOST_MEMORY";
		case CL_PROFILING_INFO_NOT_AVAILABLE:
			return "CL_PROFILING_INFO_NOT_AVAILABLE";
		case CL_MEM_COPY_OVERLAP:
			return "CL_MEM_COPY_OVERLAP";
		case CL_IMAGE_FORMAT_MISMATCH:
			return "CL_IMAGE_FORMAT_MISMATCH";
		case CL_IMAGE_FORMAT_NOT_SUPPORTED:
			return "CL_IMAGE_FORMAT_NOT_SUPPORTED";
		case CL_BUILD_PROGRAM_FAILURE:
			return "CL_BUILD_PROGRAM_FAILURE";
		case CL_MAP_FAILURE:
			return "CL_MAP_FAILURE";
		case CL_MISALIGNED_SUB_BUFFER_OFFSET:
			return "CL_MISALIGNED_SUB_BUFFER_OFFSET";
		case CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST:
			return "CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST";
		case CL_INVALID_VALUE:
			return "CL_INVALID_VALUE";
		case CL_INVALID_DEVICE_TYPE:
			return "CL_INVALID_DEVICE_TYPE";
		case CL_INVALID_PLATFORM:
			return "CL_INVALID_PLATFORM";
		case CL_INVALID_DEVICE:
			return "CL_INVALID_DEVICE";
		case CL_INVALID_CONTEXT:
			return "CL_INVALID_CONTEXT";
		case CL_INVALID_QUEUE_PROPERTIES:
			return "CL_INVALID_QUEUE_PROPERTIES";
		case CL_INVALID_COMMAND_QUEUE:
			return "CL_INVALID_COMMAND_QUEUE";
		case CL_INVALID_HOST_PTR:
			return "CL_INVALID_HOST_PTR";
		case CL_INVALID_MEM_OBJECT:
			return "CL_INVALID_MEM_OBJECT";
		case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR:
			return "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR";
		case CL_INVALID_IMAGE_SIZE:
			return "CL_INVALID_IMAGE_SIZE";
		case CL_INVALID_SAMPLER:
			return "CL_INVALID_SAMPLER";
		case CL_INVALID_BINARY:
			return "CL_INVALID_BINARY";
		case CL_INVALID_BUILD_OPTIONS:
			return "CL_INVALID_BUILD_OPTIONS";
		case CL_INVALID_PROGRAM:
			return "CL_INVALID_PROGRAM";
		case CL_INVALID_PROGRAM_EXECUTABLE:
			return "CL_INVALID_PROGRAM_EXECUTABLE";
		case CL_INVALID_KERNEL_NAME:
			return "CL_INVALID_KERNEL_NAME";
		case CL_INVALID_KERNEL_DEFINITION:
			return "CL_INVALID_KERNEL_DEFINITION";
		case CL_INVALID_KERNEL:
			return "CL_INVALID_KERNEL";
		case CL_INVALID_ARG_INDEX:
			return "CL_INVALID_ARG_INDEX";
		case CL_INVALID_ARG_VALUE:
			return "CL_INVALID_ARG_VALUE";
		case CL_INVALID_ARG_SIZE:
			return "CL_INVALID_ARG_SIZE";
		case CL_INVALID_KERNEL_ARGS:
			return "CL_INVALID_KERNEL_ARGS";
		case CL_INVALID_WORK_DIMENSION:
			return "CL_INVALID_WORK_DIMENSION";
		case CL_INVALID_WORK_GROUP_SIZE:
			return "CL_INVALID_WORK_GROUP_SIZE";
		case CL_INVALID_WORK_ITEM_SIZE:
			return "CL_INVALID_WORK_ITEM_SIZE";
		case CL_INVALID_GLOBAL_OFFSET:
			return "CL_INVALID_GLOBAL_OFFSET";
		case CL_INVALID_EVENT_WAIT_LIST:
			return "CL_INVALID_EVENT_WAIT_LIST";
		case CL_INVALID_EVENT:
			return "CL_INVALID_EVENT";
		case CL_INVALID_OPERATION:
			return "CL_INVALID_OPERATION";
		case CL_INVALID_GL_OBJECT:
			return "CL_INVALID_GL_OBJECT";
		case CL_INVALID_BUFFER_SIZE:
			return "CL_INVALID_BUFFER_SIZE";
		case CL_INVALID_MIP_LEVEL:
			return "CL_INVALID_MIP_LEVEL";
		case CL_INVALID_GLOBAL_WORK_SIZE:
			return "CL_INVALID_GLOBAL_WORK_SIZE";
		default:
			return "unknown error code";
    }
}

void ScCalc::clThrowException(int err, const char *fn, int line)
{
	if (err != CL_SUCCESS)
		throw ScCalcException("CUDA Exception at %s:%d (error %d): %s", fn, line, err, errstr(err));
}

void ScCalc::clInit()
{
	int err;
	
	if(settings.gpu.ready || !settings.gpu.device)
		return;

	// Find hardware
	cl_uint numPlatforms;
	cl_platform_id *platforms = NULL;

	err = clGetPlatformIDs(0, NULL, &numPlatforms);
	if(err == CL_SUCCESS)
		platforms = new cl_platform_id[numPlatforms];
	if(err == CL_SUCCESS)
		err = clGetPlatformIDs(numPlatforms, platforms, NULL);

	// We always use the first GPU; selection via settings.gpu.device not supported
	if(err == CL_SUCCESS)
		err = clGetDeviceIDs(platforms[0], CL_DEVICE_TYPE_GPU, 1, &gpu.device_id, NULL);
	if(platforms)
		delete platforms;
	
	if(err != CL_SUCCESS) {
		*error << "No GPU available. GPU acceleration disabled." << std::endl;
		settings.gpu.device = 0;
		return;
	}

	// Locate CL program
	char c, *p, kernel_source_fn[128];
	char *source;

	if(readlink("/proc/self/exe", kernel_source_fn, sizeof(kernel_source_fn)) > 0) {
		p = strrchr(kernel_source_fn, '/');
		if(p) p[1] = 0;
		strcat(kernel_source_fn, "sc.cl");

		FILE *in = fopen(kernel_source_fn, "r");
		if(!in) {
			*error << "No CL program found at " << kernel_source_fn << ". GPU support disabled." << std::endl;
			settings.gpu.device = 0;
			return;
		}
		
		fseek(in, 0, 2);
		int size = ftell(in);
		fseek(in, 0, 0);
		
		source = new char[size+2];
		size = fread(source, sizeof(*source), size, in);
		source[size] = 0;
		fclose(in);
	}
	
	// Setup context and command queue
	gpu.context = clCreateContext(0, 1, &gpu.device_id, NULL, NULL, &err);
	if(!gpu.context)
		gpuAssert(err);

	gpu.queue = clCreateCommandQueue(gpu.context, gpu.device_id, CL_QUEUE_PROFILING_ENABLE, &err);
	if(!gpu.queue)
		gpuAssert(err);

	// Build CL program
	gpu.program = clCreateProgramWithSource(gpu.context, 1, (const char **) &source, NULL, &err);
	if(!gpu.program)
		gpuAssert(err);
	
	err = clBuildProgram(gpu.program, 0, NULL, "-cl-single-precision-constant -cl-mad-enable -cl-no-signed-zeros -cl-fast-relaxed-math -w", NULL, NULL);
	if(err != CL_SUCCESS) {
		size_t len_status, len_options, len_log;
		char buffer_options[81920];
		char buffer_log[81920];

		cl_build_status status;
		clGetProgramBuildInfo(gpu.program, gpu.device_id, CL_PROGRAM_BUILD_STATUS, sizeof(cl_build_status), &status, &len_status);
		VERBOSE("CL_PROGRAM_BUILD_STATUS: '" << status << "'" << std::endl);

		clGetProgramBuildInfo(gpu.program, gpu.device_id, CL_PROGRAM_BUILD_OPTIONS, sizeof(buffer_options), buffer_options, &len_options);
		VERBOSE("CL_PROGRAM_BUILD_OPTIONS: '" << std::endl);
		for(uint i=0;i<len_options;++i)
			VERBOSE(buffer_options[i]);
		VERBOSE("'" << std::endl);

		clGetProgramBuildInfo(gpu.program, gpu.device_id, CL_PROGRAM_BUILD_LOG, sizeof(buffer_log), buffer_log, &len_log);
		VERBOSE("---- start build log ----" << std::endl);
		for(uint i=0;i<len_log;++i) 
			VERBOSE(buffer_log[i]);
		VERBOSE(std::endl << "----  end build log  ----" << std::endl);

		throw ScCalcException("Error: Failed to build program executable!");
	}
    
	clMakeKernel("TrimPeripheralBand");
	clMakeKernel("FindClosestNeighbor");
	delete source;

	// Get device info
	char devProp_name[128];
	size_t deviceProp_threads;
	cl_uint deviceProp_clockRate, deviceProp_multiProcessorCount;
	
	gpuAssert( clGetDeviceInfo(gpu.device_id, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(deviceProp_threads), &deviceProp_threads, NULL) );
	clGetDeviceInfo(gpu.device_id, CL_DEVICE_NAME, sizeof(devProp_name), devProp_name, NULL);
	clGetDeviceInfo(gpu.device_id, CL_DEVICE_MAX_CLOCK_FREQUENCY, sizeof(deviceProp_clockRate), &deviceProp_clockRate, NULL);
	clGetDeviceInfo(gpu.device_id, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(deviceProp_multiProcessorCount), &deviceProp_multiProcessorCount, NULL);
	settings.gpu.threads = deviceProp_threads;

	VERBOSE("GPU support enabled: " << devProp_name << 
		" [" << (deviceProp_clockRate) << " MHz] with " <<
		deviceProp_multiProcessorCount << " processors, " <<
		settings.gpu.threads <<" threads." <<
		std::endl);

	settings.gpu.ready = 1;
}

void ScCalc::clMakeKernel(const char *kname) {
	int err;
	gpu.kernels[kname] = clCreateKernel(gpu.program, kname, &err);
	if(!(gpu.kernels[kname]) || err != CL_SUCCESS) {
		VERBOSE("Error: Failed to create compute kernel " << kname << std::endl);
		gpuAssert(err);
	}
}

float ScCalc::clTrimPeripheralBand(
		const std::vector<SCCALC_DOT> &dots,
		std::vector<const SCCALC_DOT*> &trimmed_dots)
{
	int n, nBur, nAcc;
	int threads;
	float area = 0;
	clock_t timer;

	threads = MIN(512, settings.gpu.threads);
	n = dots.size();

	timer = clock();
	
	// Host and device (GPU) memory pointers for dot coordinates and results
	float4 *hAccDotCoords, *phAccDotCoords;
	float4 *hBurDotCoords, *phBurDotCoords;
	char *hDotColl;

	cl_mem dAccDotCoords;
	cl_mem dBurDotCoords;
	cl_mem dDotColl;

	// Timers
	cl_ulong start, end;
	cl_event kernelEvent;

	// Allocate host memory
	// Use cudaHostAlloc for DMA zero-copy
	hAccDotCoords = new float4[n];
	hBurDotCoords = new float4[n];
	// Allocate host memory for results (detected collisions)
	hDotColl = new char[n];

	if(!hAccDotCoords || !hBurDotCoords || !hDotColl)
		throw ScCalcException("Out of host memory!");

	// Make CUDA copy of (x, y, z) buried and accessible coordinates
	phAccDotCoords = hAccDotCoords;
	phBurDotCoords = hBurDotCoords;
	for(std::vector<SCCALC_DOT>::const_iterator iDot = dots.begin();
			iDot < dots.end(); iDot++) {
		// Quick copy
		if(iDot->buried)
			*phBurDotCoords++ = *((float4*)&iDot->coor);
		else	
			*phAccDotCoords++ = *((float4*)&iDot->coor);
	}
	nBur = phBurDotCoords - hBurDotCoords;
	nAcc = phAccDotCoords - hAccDotCoords;
	
	// Allocate GPU memory
	dBurDotCoords = clCreateBuffer(gpu.context, CL_MEM_READ_ONLY, UPPER_MULTIPLE(nBur, threads) * sizeof(*hBurDotCoords), NULL, NULL);
	dAccDotCoords = clCreateBuffer(gpu.context, CL_MEM_READ_ONLY, UPPER_MULTIPLE(nAcc, threads) * sizeof(*hAccDotCoords), NULL, NULL);
	dDotColl = clCreateBuffer(gpu.context, CL_MEM_READ_WRITE, UPPER_MULTIPLE(nBur, threads) * sizeof(*hDotColl), NULL, NULL);

	if(!dAccDotCoords || !dBurDotCoords || !dDotColl)
		throw ScCalcException("Out of GPU memory!");

	// Copy data from host to GPU
	gpuAssert( clEnqueueWriteBuffer(gpu.queue, dBurDotCoords, CL_TRUE, 0, nBur * sizeof(*hBurDotCoords), hBurDotCoords, 0, NULL, NULL) );
	gpuAssert( clEnqueueWriteBuffer(gpu.queue, dAccDotCoords, CL_TRUE, 0, nAcc * sizeof(*hAccDotCoords), hAccDotCoords, 0, NULL, NULL) );

	// Run kernel in multi-threaded blocks on GPU
	cl_kernel kernel = gpu.kernels["TrimPeripheralBand"];
	uint N = nAcc;
	size_t ldim = threads;
	size_t gdim = UPPER_MULTIPLE(nBur, threads);
	float r2 = pow(settings.band, 2);
	
	gpuAssert( clSetKernelArg(kernel, 0, sizeof(dAccDotCoords), &dAccDotCoords) );
	gpuAssert( clSetKernelArg(kernel, 1, sizeof(uint), (void*) &N) );
	gpuAssert( clSetKernelArg(kernel, 2, sizeof(dBurDotCoords), &dBurDotCoords) );
	gpuAssert( clSetKernelArg(kernel, 3, sizeof(hDotColl), &dDotColl) );
	gpuAssert( clSetKernelArg(kernel, 4, sizeof(float), (void*)&r2) );

	gpuAssert( clEnqueueNDRangeKernel(gpu.queue, kernel, 1, NULL, &gdim, &ldim, 0, NULL, &kernelEvent ) );

	// Read back data from GPU
	gpuAssert( clEnqueueReadBuffer( gpu.queue, dDotColl, CL_TRUE, 0, nBur, hDotColl, 0, NULL, NULL ) );
	
	// Wait for threads to finish and copy back memory from GPU to host
	clFinish(gpu.queue);

	clGetEventProfilingInfo(kernelEvent, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, 0);
	clGetEventProfilingInfo(kernelEvent, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, 0);

	// Make a new list of dots that have no collisions
	char *p = hDotColl;
	for(std::vector<SCCALC_DOT>::const_iterator iDot = dots.begin();
			iDot < dots.end(); iDot++) {
		const SCCALC_DOT &dot1 = *iDot;
		if(!iDot->buried)
			continue;	
		if(!*p++) {
			area += dot1.area;
			trimmed_dots.push_back(&dot1);
		}
	}
	
	float time_kernel = (end-start)*1.0e-6f;
	float time_total = GetTimerMs(timer, 0);

	VERBOSE("Peripheral trimming GPU processing time: " << time_kernel << " ms kernel, " << time_total << " ms total" << std::endl);

	// Free GPU and host memory
	gpuAssert( clReleaseMemObject(dBurDotCoords) );
	gpuAssert( clReleaseMemObject(dAccDotCoords) );
	gpuAssert( clReleaseMemObject(dDotColl) );

	delete hBurDotCoords;
	delete hAccDotCoords;
	delete hDotColl;

	return area;
}

int ScCalc::clFindClosestNeighbors(
	const std::vector<const SCCALC_DOT*> &my_dots,
	const std::vector<const SCCALC_DOT*> &their_dots,
	std::vector<const SCCALC_DOT*> &neighbors)
{
	int nMyDots, nTheirDots, nNeighbors;
	int threads;
	clock_t timer;

	timer = clock();
	threads = MIN(512, settings.gpu.threads);

	// Memory pointers for my and their dot coordinate arrays, CPU and GPU
	float4 *hMyDotCoords, *phMyDotCoords;
	float4 *hTheirDotCoords, *phTheirDotCoords;
	cl_mem dMyDotCoords, dTheirDotCoords;

	// Dot point pointer map
	SCCALC_DOT **hTheirDots, **phTheirDots;

	// Neighbor ID memory pointers
	uint *hNeighbors;
	cl_mem dNeighbors;

	// Timers
	cl_ulong start, end;
	cl_event kernelEvent;

	nMyDots = my_dots.size();
	nTheirDots = their_dots.size();
	nNeighbors = nMyDots;

	// Allocate host memory
	hMyDotCoords = new float4 [nMyDots];
	hTheirDotCoords = new float4 [nTheirDots];
	hTheirDots = (SCCALC_DOT**) new void* [nTheirDots];
	hNeighbors = new uint[nNeighbors];
	
	// Make CUDA copy of (x, y, z) dot coordinates for my dots
	phMyDotCoords = hMyDotCoords;
	for(std::vector<const SCCALC_DOT*>::const_iterator iDot = my_dots.begin(); iDot < my_dots.end(); iDot++) {
		// Quick copy
		phMyDotCoords->x = (*iDot)->coor.x;
		phMyDotCoords->y = (*iDot)->coor.y;
		phMyDotCoords->z = (*iDot)->coor.z;
		phMyDotCoords++;
	}
	nMyDots = phMyDotCoords - hMyDotCoords;

	// Make CUDA copy of (x, y, z) dot coordinates for their dots and keep a map
	phTheirDotCoords = hTheirDotCoords;
	phTheirDots = hTheirDots;
	for(std::vector<const SCCALC_DOT*>::const_iterator iDot = their_dots.begin(); iDot < their_dots.end(); iDot++) {
		// Quick copy
		if(!(*iDot)->buried)
			continue;

		phTheirDotCoords->x = (*iDot)->coor.x;
		phTheirDotCoords->y = (*iDot)->coor.y;
		phTheirDotCoords->z = (*iDot)->coor.z;
		phTheirDotCoords++;
		
		*phTheirDots++ = (SCCALC_DOT*)*iDot;
	}
	nTheirDots = phTheirDotCoords - hTheirDotCoords;
	
	// Allocate GPU memory and copy data there
	dMyDotCoords = clCreateBuffer(gpu.context, CL_MEM_READ_ONLY, UPPER_MULTIPLE(nMyDots, threads) * sizeof(*hMyDotCoords), NULL, NULL);
	dTheirDotCoords = clCreateBuffer(gpu.context, CL_MEM_READ_ONLY, UPPER_MULTIPLE(nTheirDots, threads) * sizeof(*hTheirDotCoords), NULL, NULL);
	dNeighbors = clCreateBuffer(gpu.context, CL_MEM_READ_WRITE, UPPER_MULTIPLE(nNeighbors, threads) * sizeof(*hNeighbors), NULL, NULL);
	
	if(!dMyDotCoords || !dTheirDotCoords || !dNeighbors)
		throw ScCalcException("Out of GPU mempory!");

	gpuAssert( clEnqueueWriteBuffer(gpu.queue, dMyDotCoords, CL_TRUE, 0, nMyDots * sizeof(*hMyDotCoords), hMyDotCoords, 0, NULL, NULL) );
	gpuAssert( clEnqueueWriteBuffer(gpu.queue, dTheirDotCoords, CL_TRUE, 0, nTheirDots * sizeof(*hTheirDotCoords), hTheirDotCoords, 0, NULL, NULL) );

	// Run kernel in multi-threaded blocks on GPU
	cl_kernel kernel = gpu.kernels["FindClosestNeighbor"];
	size_t ldim = threads;
	size_t gdim = UPPER_MULTIPLE(nMyDots, threads);
	uint N = nTheirDots;
	
	gpuAssert( clSetKernelArg(kernel, 0, sizeof(dMyDotCoords), &dMyDotCoords) );
	gpuAssert( clSetKernelArg(kernel, 1, sizeof(dTheirDotCoords), &dTheirDotCoords) );
	gpuAssert( clSetKernelArg(kernel, 2, sizeof(nTheirDots), &nTheirDots) );
	gpuAssert( clSetKernelArg(kernel, 3, sizeof(dNeighbors), &dNeighbors) );

	gpuAssert( clEnqueueNDRangeKernel(gpu.queue, kernel, 1, NULL, &gdim, &ldim, 0, NULL, &kernelEvent ) );

	// Read back data from GPU
	gpuAssert( clEnqueueReadBuffer( gpu.queue, dNeighbors, CL_TRUE, 0, nMyDots * sizeof(*hNeighbors), hNeighbors, 0, NULL, NULL ) );
	
	// Wait for threads to finish and copy back memory from GPU to host
	clFinish(gpu.queue);

	for(int i = 0; i < nNeighbors; i++)
		neighbors.push_back( hTheirDots[hNeighbors[i]] );

	clGetEventProfilingInfo(kernelEvent, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, 0);
	clGetEventProfilingInfo(kernelEvent, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, 0);

	float time_kernel = (end-start)*1.0e-6f;
	float time_total = GetTimerMs(timer, 0);

	VERBOSE("Peripheral trimming GPU processing time: " << time_kernel << " ms kernel, " << time_total << " ms total" << std::endl);

	// Free memory
	gpuAssert( clReleaseMemObject(dMyDotCoords) );
	gpuAssert( clReleaseMemObject(dTheirDotCoords) );
	gpuAssert( clReleaseMemObject(dNeighbors) );

	delete hMyDotCoords;
	delete hTheirDotCoords;
	delete hTheirDots;
	delete hNeighbors;

	return 1;
}

#endif

// END //
