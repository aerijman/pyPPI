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
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include "vec3.h"

//////////////////////////////////////////////////////////////////////
// GPU acceleration defs

// To enable GPU support, set define below and compile code with nvcc
// No additional libraries are needed, but nvcc will link against
// libcudart.so which is needed at runtime.

// #define CUDA_GPU
// #define OPENCL_GPU

#ifdef CUDA_GPU
#define GPU
#include <time.h>
#include <cuda_runtime_api.h>
float GetTimerMs(clock_t &start, int reset =0);
void _cuda_TrimPeripheralBand(int x, int y, float3 *dAccDotCoords, uint nAcc, float3 *dBurDotCoords, char *dDotColl, float r2);
void _cuda_FindClosestNeighbor(int x, int y, float3 *dMyDots, float3 *dTheirDotCoords, uint nTheirDotCoords, uint *dNeighbors);
#endif

#ifdef OPENCL_GPU
#define GPU
#ifdef MAC
#include <OpenCL/cl_platform.h>
#include <OpenCL/opencl.h>
#include <OpenCL/cl.h>
#else
#include <CL/cl_platform.h>
#include <CL/opencl.h>
#include <CL/cl.h>
#endif
#endif

////////////////////////////////////////////////////////////
// Defs from sc source

enum {
	ATTEN_BLOCKER = 1,
	ATTEN_2  = 2,
	ATTEN_BURIED_FLAGGED = 5,
	ATTEN_6 = 6
};

#define MAX_SUBDIV 100

////////////////////////////////////////////////////////////
// Types

#if not defined MIN
#define MIN(a,b) ((a) < (b) ? (a): (b))
#define MAX(a,b) ((a) > (b) ? (a): (b))
#endif

#define ABS(a) (((a) < 0) ? (-a) : (a))
#define VERBOSE(x) if(verbose) *verbose << x;
 
typedef struct _SCCALC_SETTINGS {
	const char *sc_radii;

	// From sc source	
	float rp;
	float density;
	float band;
	float sep;
	float weight;
	float binwidth_dist;
	float binwidth_norm;

	int hydrogens;
#ifdef GPU
	struct {
		int ready;
		int device;
		int threads;
		int proc;
	} gpu;
#endif
	
} SCCALC_SETTINGS;

typedef struct _SCCALC_RESULTS {
	float sc;
	float separation;
	float area;
	int nAtoms;
	struct {
		float d_mean;
		float d_median;
		float s_mean;
		float s_median;
		int nAtoms;
		int nBuriedAtoms;
		int nBlockedAtoms;
		int nAllDots;
		int nTrimmedDots;
		int nBuriedDots;
		int nAccessibleDots;
		float trimmedArea;
	} surface[3];
	struct {
		int convex;
		int concave;
		int toroidal;
	} dots;
} SCCALC_RESULTS;
	
// Atom radius definition
typedef struct _SCCALC_ATOM_RADIUS {
	char residue[5];
	char atom[5];
	float radius;
} SCCALC_ATOM_RADIUS;

// Molecular dot
class PDBAtom;
typedef struct _SCCALC_DOT {
	CVec3 coor;
	int buried;
	int type;
	float area;
	CVec3 outnml;
	const PDBAtom *atom;
} SCCALC_DOT;

// Molecular probe
typedef struct _SCCALC_PROBE {
	const PDBAtom *pAtoms[3];
	float height;
	CVec3 point;
	CVec3 alt;
} SCCALC_PROBE;

// PDB Atom
class PDBAtom : public CVec3 {
	public:
	float occ;
	float b;
	int natom;
	int nresidue;
	char atom[4];
	char residue[4];
	char chain;
	
	protected:
	friend class ScCalc;
	
	int molecule;
	float radius;
	float density;
	int atten;
	int access;
	std::vector<PDBAtom*> neighbors;
	std::vector<PDBAtom*> buried;
	
	public:
	PDBAtom() : CVec3() {
		occ = 0;
		b = 0;
		natom = 0;
		nresidue = 0;
		chain = ' ';

		molecule = 0;
		radius = 0;
		density = 0;
		atten = 0;
		access = 0;

		memset(atom, 0, sizeof(atom));
		memset(residue, 0, sizeof(residue));
	}
		
	operator std::string() {
		char buf[32];
		snprintf(buf, sizeof(buf), "%c:%d(%s):%s", chain, nresidue, residue,
				((std::string)atom).c_str());
		return std::string(buf);
	}

	int operator ==(const PDBAtom &atom2) {
		// For speed reasons, rather than comparing the coordinates,
		// we'll compare pointers as all atoms are part of the same atoms
		// vector, so the same atom will have the same address
		return this == &atom2;
	}
	int operator <=(PDBAtom &atom2) {
		return this <= &atom2;
	}
};

////////////////////////////////////////////////////////////
// Shape Complementarity Calculator class definition

class ScCalc {

	public:
	std::ostream *error, *verbose;
	SCCALC_SETTINGS settings;

	ScCalc();
	~ScCalc();
	int Init();
	void Reset();

	std::vector<PDBAtom> ReadPdb(const char *fn);
	int AddAtom(int molecule, PDBAtom &atom);
	SCCALC_RESULTS *CalcSc();

	const SCCALC_RESULTS& GetResults() { return results; }
	const std::vector<PDBAtom>& GetAtoms() { return atoms; }
	
	protected:
	static std::vector<SCCALC_ATOM_RADIUS> radii;
	std::vector<PDBAtom> atoms;
	std::vector <SCCALC_DOT> dots[2];
	std::vector<SCCALC_PROBE> probes;
	SCCALC_RESULTS results;

	int AssignAtomRadius(PDBAtom &atom);
	int WildcardMatch(const char *r, const char *pattern);
	int ReadScRadii(const char *fn);
	void AddDot(const int molecule, const int type, const CVec3 coor, const float area, const CVec3 pcen, const PDBAtom &atom);

	private:
	// variable from original Fortran code, keeping original names
	float radmax;

	int AssignAttentionNumbers(std::vector<PDBAtom>& atom);
	int CalcDotsForAllAtoms(std::vector<PDBAtom>& atoms);
	int CalcDotsForAtoms(std::vector<PDBAtom>& atoms);
	int FindNeighbordsAndBuriedAtoms(PDBAtom& atom);
	int FindNeighborsForAtom(PDBAtom& atom1);

	int GenerateToroidalSurface(PDBAtom& atom1, PDBAtom& atom2, const CVec3 uij, const CVec3 tij, float rij, int between);
	int GenerateConvexSurface(const PDBAtom& atom1);
	int GenerateConcaveSurface();

	// Function names similar to original source
	int SecondLoop(PDBAtom &pAtom1);
	int ThirdLoop(PDBAtom &pAtom1, PDBAtom &pAtom, const CVec3 &uij, const CVec3 &tij, const float rij);
	int CheckAtomCollision2(const CVec3 &pijk, const PDBAtom &atom1, const PDBAtom &atom2, const std::vector<PDBAtom*> &atoms);
	int CheckPointCollision(const CVec3 &pcen, const std::vector<PDBAtom*> &atoms);
	int CheckProbeCollision(const CVec3 &point, const std::vector<const SCCALC_PROBE*> nears, const float r2);

	// Dot trimming
	float TrimPeripheralBand(const std::vector<SCCALC_DOT> &sdots, std::vector<const SCCALC_DOT*> &trimmed_dots);
	int TrimPeripheralBandCheckDot(const SCCALC_DOT &dot, const std::vector<SCCALC_DOT> &sdots);

	// Elementary functions
	float DistancePointToLine(const CVec3 &cen, const CVec3 &axis, const CVec3 &pnt);
	float SubArc(const CVec3 &cen, const float rad, const CVec3 &axis, const float density,	const CVec3 &x, const CVec3 &v, std::vector<CVec3> &points);
	float SubDiv(const CVec3 &cen, const float rad, const CVec3 &x, const CVec3 &y, float angle, float density, std::vector<CVec3> &points);
	float SubCir(const CVec3 &cen, const float rad, const CVec3 &north, const float density, std::vector<CVec3> &points);

	// Kernel functions for parallalization
	int CalcNeighborDistance(const int molecule, const std::vector<const SCCALC_DOT*> &my_dots, const std::vector<const SCCALC_DOT*> &their_dots);
	const SCCALC_DOT *CalcNeighborDistanceFindClosestNeighbor(const SCCALC_DOT &dot1, const std::vector<const SCCALC_DOT*> &their_dots);

#ifdef CUDA_GPU
	#define GPU
	#define gpuFindClosestNeighbors cudaFindClosestNeighbors
	#define gpuTrimPeripheralBand cudaTrimPeripheralBand
	#define gpuInit	cudaInit
	#define gpuThrowException cudaThrowException

	protected:
	void cudaThrowException(cudaError_t err, const char *fn, int line);
	float cudaTrimPeripheralBand(const std::vector<SCCALC_DOT> &dots, std::vector<const SCCALC_DOT*> &trimmed_dots);
	int cudaFindClosestNeighbors(const std::vector<const SCCALC_DOT*> &my_dots, const std::vector<const SCCALC_DOT*> &their_dots, std::vector<const SCCALC_DOT*> &neighbors);

	public:
	void cudaInit();
#endif

#ifdef OPENCL_GPU
	#define GPU
	#define gpuFindClosestNeighbors clFindClosestNeighbors
	#define gpuTrimPeripheralBand clTrimPeripheralBand
	#define gpuInit	clInit
	#define gpuThrowException clThrowException

	private:
	struct {
		std::map<std::string, cl_kernel> kernels;
		cl_device_id device_id;
		cl_context context;
		cl_command_queue queue;
		cl_program program;
	} gpu;

	protected:
	void clThrowException(int err, const char *fn, int line);
	float clTrimPeripheralBand(const std::vector<SCCALC_DOT> &dots, std::vector<const SCCALC_DOT*> &trimmed_dots);
	int clFindClosestNeighbors(const std::vector<const SCCALC_DOT*> &my_dots, const std::vector<const SCCALC_DOT*> &their_dots, std::vector<const SCCALC_DOT*> &neighbors);
	void clMakeKernel(const char *kname);

	public:
	void clInit();
#endif

};

class ScCalcException {

	public:
	std::string error;

	ScCalcException(const char *err, ...) {
		va_list p;
   	char buf[256];
		va_start(p, err);
		vsnprintf(buf, sizeof(buf), err, p);
		va_end(p);
		error = buf;
	}
};
