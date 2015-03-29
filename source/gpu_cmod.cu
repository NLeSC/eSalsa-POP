#include <stdio.h>

/*
 * This file contains the necessary functions for interfacing Fortran with the CUDA runtime API
 *
 *
 */



extern "C" {

  void cuda_init(int *pmy_task);

  void cudamallochost(void **hostptr, int *p_size);

}



int my_task;
int cuda_initialized = 0;

//Fortran entry for initializing CUDA
void cuda_init(int *pmy_task) {
  if (cuda_initialized == 0) {
    cuda_initialized = 1;
    my_task = *pmy_task;
    int deviceCount = 0;

    cudaError_t err = cudaGetDeviceCount(&deviceCount);
    if (err != cudaSuccess) fprintf(stderr, "Error in cuda initialization: %s\n", cudaGetErrorString( err ));

    if (deviceCount < 1) {
      fprintf(stderr,"Error: less than 1 cuda capable device detected proc=%d\n", my_task);
    }

    int dev = my_task % deviceCount;
    //fprintf(stdout,"Process %d: using CUDA device %d\n",my_task,dev);

    cudaSetDeviceFlags(cudaDeviceMapHost);
    cudaSetDevice(dev);

    cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);
    cudaDeviceSynchronize();
  }
}


//Fortran entry for allocating pinned memory
void cudamallochost(void **hostptr, int *p_size) {
  if (!cuda_initialized) {
    printf("Error: cudamallochost called before cuda_init\n");
  }
 
  cudaError_t err;

  err = cudaHostAlloc((void **)hostptr, (*p_size)*sizeof(double), cudaHostAllocMapped);
  if (err != cudaSuccess) fprintf(stderr, "Error in cudaHostAlloc: %s\n", cudaGetErrorString( err ));
}
