#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define SOFTENING 1e-9f

typedef struct { float x, y, z, vx, vy, vz; } Body;

void randomizeBodies(float *data, int n) {
  for (int i = 0; i < n; i++) {
    data[i] = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
  }
}

void bodyForce(Body *p, float dt, int n, int nthreads) {
  int i, j;
  #ifdef isParallel
  #pragma omp parallel for default(none) shared(p, dt, n) private(i, j) num_threads(nthreads) schedule(static)
  #endif
  for (i = 0; i < n; i++) {
    float Fx = 0.0f; float Fy = 0.0f; float Fz = 0.0f;
    for (j = 0; j < n; j++) {
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
}

int main(const int argc, const char** argv) {
  FILE *datafile;  
  int nBodies;
  int nthreads;

  if (argc > 1) nBodies = atoi(argv[1]);
  if (argc > 2) nthreads = atoi(argv[2]);

  const float dt = 0.01f; // time step
  const int nIters = 20;  // simulation iterations

  int bytes = nBodies*sizeof(Body);
  float *buf = (float*)malloc(bytes);
  Body *p = (Body*)buf;

  randomizeBodies(buf, 6*nBodies); // Init pos / vel data

  double totalTime = 0.0;

  datafile = fopen("nbody.csv","w");
  printf("%d %d\n", nBodies, nthreads);

  /* ------------------------------*/
  /*     MAIN LOOP                 */
  /* ------------------------------*/
  double start = omp_get_wtime();
  for (int iter = 1; iter <= nIters; iter++) {
    printf("iteration:%d\n", iter);

    bodyForce(p, dt, nBodies, nthreads);           // compute interbody forces

    #ifdef isParallel
    #pragma omp parallel for num_threads(128) schedule(static)
    #endif
    for (int i = 0 ; i < nBodies; i++) { // integrate position
      p[i].x += p[i].vx*dt;
      p[i].y += p[i].vy*dt;
      p[i].z += p[i].vz*dt;
    }
  }
  
  const double tElapsed = omp_get_wtime(); 
  fclose(datafile);
  totalTime = (tElapsed - start);
  double avgTime =  totalTime / (double)(nIters-1); 

  printf("avgTime: %f   totTime: %f \n", avgTime, totalTime);
  free(buf);
}
