all:	sc gpu-sc gpu-sc-cl
sc:	sc.cc main.cc sc.h vec3.h 
	g++ -O3 sc.cc main.cc -o sc
gpu-sc:	sc.cc main.cc sc.h vec3.h sc.cu
	nvcc sc.cc main.cc sc.cu -D CUDA_GPU -O3 -o gpu-sc
gpu-sc-cl:	sc.cc main.cc sc.h vec3.h sc.cu
	g++ sc.cc main.cc -D OPENCL_GPU -O3 -o gpu-sc-cl -l OpenCL -I /usr/local/cuda/include/ -L /usr/local/cuda/lib64/
clean:
	rm sc
package:
	tar czf sc.tar.gz main.cc sc.cc sc.h vec3.h sc.cu sc.cl Makefile sc_radii.lib sc gpu-sc gpu-sc-cl
