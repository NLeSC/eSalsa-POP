#include <stdio.h>

/*
 * This file contains the necessary functions for interfacing Fortran with the CUDA runtime API
 *
 *
 */



extern "C" {


  void cudamallochost(void **hostptr, int *p_size);

}



int cuda_initialized = 0;

void cuda_init() {
  cudaSetDeviceFlags(cudaDeviceMapHost);
  cudaSetDevice(0);
  cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);
  cudaDeviceSynchronize();


}


//Fortran entry for allocating pinned memory
void cudamallochost(void **hostptr, int *p_size) {
  if (!cuda_initialized) {
    cuda_initialized = 1;
    cuda_init();
  }
 
  cudaError_t err;

  err = cudaHostAlloc((void **)hostptr, (*p_size)*sizeof(double), cudaHostAllocMapped);
  if (err != cudaSuccess) fprintf(stderr, "Error in cudaHostAlloc: %s\n", cudaGetErrorString( err ));
}
