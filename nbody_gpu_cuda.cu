#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include "timer.h"

#define SOFTENING 1e-9f
#define MAX_BLOCKS_PER_DIM 65535
#define MIN(a, b) (((a) < (b)) ? (a) : (b))

typedef struct { float x, y, z, vx, vy, vz; } Body;

void randomizeBodies(Body *p, int n) {
    for (int i = 0; i < n; i++) {
        p[i].x = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
        p[i].y = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
        p[i].z = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
        p[i].vx = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
        p[i].vy = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
        p[i].vz = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
    }
}

__global__ void bodyForce(Body *p, float dt, int n) {
  int tid0 = blockIdx.x * blockDim.x + threadIdx.x;
  for (int i = tid0; i < n; i += blockDim.x * gridDim.x) {
    float Fx = 0.0f; float Fy = 0.0f; float Fz = 0.0f;

    for (int j = 0; j < n; j++) {
      float dx = p[j].x - p[i].x;
      float dy = p[j].y - p[i].y;
      float dz = p[j].z - p[i].z;
      float distSqr = dx*dx + dy*dy + dz*dz + SOFTENING;
      float invDist = 1.0f / sqrtf(distSqr);
      float invDist3 = invDist * invDist * invDist;

      Fx += dx * invDist3; Fy += dy * invDist3; Fz += dz * invDist3;
    }
    p[i].vx += dt*Fx; p[i].vy += dt*Fy; p[i].vz += dt*Fz;
  }
  for (int i = tid0 ; i < n; i += blockDim.x * gridDim.x) { // integrate position
    p[i].x = p[i].x + p[i].vx*dt;
    p[i].y = p[i].y + p[i].vy*dt;
    p[i].z = p[i].z + p[i].vz*dt;
  }
}

int particle_positions_to_csv(FILE *datafile, int iter, Body *p, int nBodies) {
    for (int i = 0 ; i < nBodies; i++) {
        fprintf(datafile, "%i, %f, %f, %f\n", iter, p[i].x, p[i].y, p[i].z);
    }
    return 0;
}

void calcPosition(Body *p, float dt, int nBodies, int nIters, int nthreads_per_block)
{
  FILE *datafile = fopen("nbody.csv","w");
  double totalTime = 0.0;
  Body *p_cu;
  cudaMalloc((void**)&p_cu, nBodies*sizeof(Body));
  cudaMemset(p_cu, 0.0, nBodies*sizeof(Body));
  cudaMemcpy(p_cu, p, nBodies*sizeof(Body), cudaMemcpyHostToDevice);
  StartTimer();
  for (int iter = 1; iter <= nIters; iter++) {
    int nblocks = MIN(nBodies / nthreads_per_block + 1, MAX_BLOCKS_PER_DIM);


    bodyForce<<<nblocks,nthreads_per_block>>>(p_cu, dt, nBodies);       // compute interbody forces
    
    const double tElapsed = GetTimer() / 1000.0;
    if (iter > 1) {                      // First iter is warm up
      totalTime += tElapsed; 
    }

    // cudaMemcpy(p, p_cu, nBodies*sizeof(Body), cudaMemcpyDeviceToHost);
    // if(iter%100 == 0){
    // particle_positions_to_csv(datafile, iter, p, nBodies);
    // }
  }
  // totalTime = totalTime/1000;
  double avgTime = totalTime / (double)(nIters-1);

  printf("avgTime: %.10f   totTime: %.10f \n", avgTime, totalTime);
  fclose(datafile);
}

int main(const int argc, const char** argv) {
  int nBodies;
  int nthreads_per_block;
  if (argc > 1) nBodies = atoi(argv[1]);
  if (argc > 2) nthreads_per_block = atoi(argv[2]);

  float dt = 0.01f; // time step
  int nIters = 20;  // simulation iterations

  int bytes = nBodies*sizeof(Body);
  float *buf = (float*)malloc(bytes);
  Body *p = (Body*)buf;

  randomizeBodies(p, nBodies); // Init pos / vel data
  calcPosition(p, dt, nBodies, nIters, nthreads_per_block);

  free(buf);
}
