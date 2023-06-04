nbody_serial: timer.h timer.c nbody_cpu_serial.c
	gcc -o nbody_serial timer.c nbody_cpu_serial.c -lm -O3 -march=native -mtune=native

nbody_openmp: nbody_cpu_openmp.c
	gcc -fopenmp -DisParallel -O3 -flto -march=native -mtune=native nbody_cpu_openmp.c -o nbody_openmp -lm

nbody_cuda: timer.h timer.c nbody_gpu_cuda.cu
	cp timer.c timer.cu
	nvcc -o nbody_cuda timer.cu nbody_gpu_cuda.cu -arch sm_70