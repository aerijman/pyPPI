//////////////////////////////////////////////////////////////////////
// GPU-SC: GPU-accelerated version of the original Lawrence &
// Coleman shape complementarity program from CCP4.
// Luki Goldschmidt <luki@mbi.ucla.edu>, March 2011
//////////////////////////////////////////////////////////////////////
// GPU kernel functions and stubs
//////////////////////////////////////////////////////////////////////

#include "sc.h"

//////////////////////////////////////////////////////////////////////
// Collision checking GPU kernel for TrimPeripheralBand  

__global__ void _cuda_TrimPeripheralBand_kernel(
	float3 *dAccDotCoords, 
	uint nAcc, 
	float3 *dBurDotCoords, 
	char *dDotColl, 
	float r2)
{
	register int i, j, l;
	register float3 dot1;
	__shared__ char sColl[1024];
	__shared__ float3 sCoords[1024];

	sColl[threadIdx.x] = 0;
 	dot1 = dBurDotCoords[blockIdx.x*blockDim.x + threadIdx.x];

	for(i = 0; i < nAcc; i += blockDim.x) {
		__syncthreads();
		sCoords[threadIdx.x] = dAccDotCoords[i + threadIdx.x];
		__syncthreads();

		l = MIN(nAcc - i, blockDim.x);
		for(j = 0; j < l; j++) {
			register float3 dot2 = sCoords[j];
			dot2.x -= dot1.x;
			dot2.y -= dot1.y;
			dot2.z -= dot1.z;
			sColl[threadIdx.x] |= (dot2.x*dot2.x + dot2.y*dot2.y + dot2.z*dot2.z) <= r2;
		}
	}
	dDotColl[blockIdx.x*blockDim.x + threadIdx.x] = sColl[threadIdx.x];
}

//////////////////////////////////////////////////////////////////////
// Finding closest dot neighbor GPU kernel  

__global__ void _cuda_FindClosestNeighbor_kernel(
	float3 *dMyDotCoords, 
	float3 *dTheirDotCoords, 
	uint nTheirDots, 
	uint *dNeighbors)
{
	register int i, j, l;
	register float3 dot1;
	__shared__ uint sNeighbors[512];
	__shared__ float3 sCoords[512];
	float distmin = 99999.0, d2;

 	dot1 = dMyDotCoords[blockIdx.x*blockDim.x + threadIdx.x];

	for(i = 0; i < nTheirDots; i += blockDim.x) {
		__syncthreads();
		sCoords[threadIdx.x] = dTheirDotCoords[i + threadIdx.x];
		__syncthreads();

		l = MIN(nTheirDots - i, blockDim.x);
		for(j = 0; j < l; j++) {
			register float3 dot2 = sCoords[j];
			dot2.x -= dot1.x;
			dot2.y -= dot1.y;
			dot2.z -= dot1.z;
			d2 = dot2.x*dot2.x + dot2.y*dot2.y + dot2.z*dot2.z;
			if(d2 <= distmin) {
				distmin = d2;
				sNeighbors[threadIdx.x] = i+j;
			}
		}
	}
	dNeighbors[blockIdx.x*blockDim.x + threadIdx.x] = sNeighbors[threadIdx.x];
}

//////////////////////////////////////////////////////////////////////
// Stubs called from CPU code

void _cuda_TrimPeripheralBand(int x, int y, float3 *dAccDotCoords, uint nAcc, float3 *dBurDotCoords, char *dDotColl, float r2)
{
	_cuda_TrimPeripheralBand_kernel<<<x, y>>>(dAccDotCoords, nAcc, dBurDotCoords, dDotColl, r2);
}

void _cuda_FindClosestNeighbor(int x, int y, float3 *dMyDotCoords, float3 *dTheirDotCoords, uint nTheirDotsCoords, uint *dNeighbors)
{
	_cuda_FindClosestNeighbor_kernel<<<x, y>>>(dMyDotCoords, dTheirDotCoords, nTheirDotsCoords, dNeighbors);
}
