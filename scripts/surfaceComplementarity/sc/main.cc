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
#include <getopt.h>
#include "sc.h"

////////////////////////////////////////////////////////////

using namespace std;
char sc_radii[256] = "sc_radii.lib";

void assign_molecules(vector<PDBAtom> &atoms, ScCalc &sc, char *chains[2]);
void init(int argc, char * const argv[], ScCalc &sc, int &fn_index, char *chains[2]);
void report_result(const char *fn, SCCALC_RESULTS &r, ScCalc &sc);

int main(int argc, char * const argv[])
{
	ScCalc sc;
	SCCALC_RESULTS *r = NULL;
	int fn_index = 0;
	const char *fn = NULL;
	char _chain[2][2] = { "\0", "\0" };
	char *chains[2] = { _chain[0], _chain[1] };
	
	init(argc, argv, sc, fn_index, chains);

	if(fn_index >= argc) {
		cerr << "This program computes the Lawrence & Coleman shape complementarity" << endl;
		cerr << "between two molecules (based on the sc code from CCP4)." << endl;
		cerr << "Luki Goldschmidt <luki@mbi.ucla.edu>, March 2011" << endl;
		cerr << endl;
		cerr << "Usage: [options] <PDB files>" << endl;
		cerr << "\nOptions:" << endl;
		cerr << "-1:\tMolecule 1 chain IDs (default: first chain)" << endl;
		cerr << "-2:\tMolecule 1 chain IDs (default: second chain)" << endl;
		cerr << "-H:\tInclude hydrogen atoms in the calculation (default: skip hydrogens)" << endl;
		cerr << "-q:\tQuick mode settings (lower accuracy)" << endl;
		cerr << "-v:\tShow verbose output" << endl;
		cerr << "-l:\tFilename of the atom radius library (default: sc_radii.lib)" << endl;
		cerr << "-r:\tProbe radius (default: 1.7A)" << endl;
		cerr << "-t:\tTrimming distance for the peripheral shell (default: 1.5A)" << endl;
		cerr << "-w:\tWeight factor using in sc calculation (default: 0.5A)" << endl;
		cerr << "-d:\tDot density for molecular surface (default: 15 dots/A^2)" << endl;
		cerr << "-s:\tInterface separation distance (default: 8A)" << endl;
#ifdef CUDA_GPU
		cerr << "-g:\tEnable CUDA GPU acceleration (default: use if available)" << endl;
		cerr << "-T:\tGPU thread limit (default: max supported by hardware, max 1024)" << endl;
#endif
		return 1;
	}

	while(fn_index < argc) {
		fn = argv[fn_index++];
		try {
			vector<PDBAtom> atoms;
			sc.Init();
			sc.Reset();

			atoms = sc.ReadPdb(fn);
			assign_molecules(atoms, sc, chains);
			r = sc.CalcSc();
			if(r)
				report_result(fn, *r, sc);

		} catch(ScCalcException e) {
			cout << "Failed to process " << fn << ": " << e.error << endl;
		        continue;
		}
	}
	
	return 0;
}

// Simple front-end interface to change settings for a this stand alone
// application.
void init(int argc, char * const argv[], ScCalc &sc, int &fn_index, char *chains[2])
{
	char c, *p;

	if(readlink("/proc/self/exe", sc_radii, sizeof(sc_radii)) > 0) {
		p = strrchr(sc_radii, '/');
		if(p) p[1] = 0;
		strcat(sc_radii, "sc_radii.lib");
		sc.settings.sc_radii = sc_radii;
	}

	sc.error = &std::cerr;
	while ( (c = getopt(argc, argv, "vl:r:t:w:d:s:q1:2:g:T:H")) != -1 ) {
		switch(c) {
		case 'v':
			sc.verbose = &std::cout;
			break;
		case 'l':
			sc.settings.sc_radii = argv[optind-1];
			break;
		case 'r':
			sc.settings.rp = atof(argv[optind-1]);
			break;
		case 't':
			sc.settings.band = atof(argv[optind-1]);
			break;
		case 'w':
			sc.settings.weight = atof(argv[optind-1]);
			break;
		case 'd':
			sc.settings.density = atof(argv[optind-1]);
			break;
		case 's':
			sc.settings.sep = atof(argv[optind-1]);
			break;
		case 'q':
			sc.settings.density = 5.0;
			break;
		case '1':
		case '2':
			chains[c - '1'] = argv[optind-1];
			break;
		case 'H':
			sc.settings.hydrogens = 1;
#ifdef GPU
		case 'g':
			sc.settings.gpu.device = atoi(argv[optind-1]);
			break;
		case 'T':
			sc.settings.gpu.threads = atoi(argv[optind-1]);
			break;
#else
		case 'g':
			cerr << "This build does not support GPU acceleration." << std::endl;
			break;
#endif
		}
	}

	fn_index = optind;
}

// Simple molecule assignment by chain only
void assign_molecules(vector<PDBAtom> &atoms, ScCalc &sc, char *chains[2])
{
	vector<PDBAtom>::iterator atom;
	for(atom = atoms.begin(); atom < atoms.end(); atom++) {
		// Skip hydrogens
		if(!sc.settings.hydrogens) {
			int i;
			for(i = 0; i < sizeof(atom->atom); i++) {
				if(atom->atom[i] < '0' || atom->atom[i] > '9')
					break;
			}
			if(atom->atom[i] == 'H')
				continue;
		}
		// Auto assign chain IDs
		if(!*chains[0])
			*chains[0] = atom->chain;
		if(!*chains[1] && atom->chain != *chains[0])
			*chains[1] = atom->chain;

		if(strchr(chains[0], atom->chain))
			sc.AddAtom(0, *atom);
		else if(strchr(chains[1], atom->chain))
			sc.AddAtom(1, *atom);
	}
}

void report_result(const char *fn, SCCALC_RESULTS &r, ScCalc &sc)
{
	if(!sc.verbose) {
		printf("%s\t%.4f\t%.4f\t%.f\n", fn, r.sc, r.separation, r.area);
		return;
	}

	// Verbose view
	cout << "==================================================" << endl;
	cout << endl;
	for(int i = 0; i <= 2; i++) {
		if(i < 2)
			cout << "Molecule " << (i+1) << ":" << endl;
		else
			cout << "Total/Average for both molecules:" << endl;

		cout << "          Total Atoms: " << r.surface[i].nAtoms << endl;
		cout << "         Buried Atoms: " << r.surface[i].nBuriedAtoms << endl;
		cout << "        Blocked Atoms: " << r.surface[i].nBlockedAtoms << endl;
		cout << "           Total Dots: " << r.surface[i].nAllDots << endl;
//		cout << "          Buried Dots: " << r.surface[i].nBuriedDots << endl;
//		cout << "      Accessible Dots: " << r.surface[i].nAccessibleDots << endl;
		cout << " Trimmed Surface Dots: " << r.surface[i].nTrimmedDots << endl;
		cout << "         Trimmed Area: " << r.surface[i].trimmedArea << endl;
		cout << endl;
	}
	cout << endl;

	for(int i = 0; i <= 2; i++) {
		if(i < 2)
			cout << "Molecule " << (i+1) << "->" << ((i+1)%2+1) << ": " << endl;
		else	
			cout << "Average for both molecules:" << endl;
		cout << "      Mean Separation: " << r.surface[i].d_mean << endl;
		cout << "    Median Separation: " << r.surface[i].d_median << endl;
		cout << "    Mean Shape Compl.: " << r.surface[i].s_mean << endl;
		cout << "  Median Shape Compl.: " << r.surface[i].s_median << endl;
		cout << endl;
	}

	cout << "==================================================" << endl;
	cout << "Shape Complementarity:          " << r.sc << endl;
	cout << "Interface separation (A):       " << r.separation << endl;
	cout << "Area buried in interface (A^2): " << r.area << endl;
	cout << "==================================================" << endl;
}
