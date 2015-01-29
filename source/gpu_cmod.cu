#include <stdio.h>





extern "C" {


  void cudamallochost(void **hostptr, int* size);

}






//Fortran entry for allocating pinned memory
void cudamallochost(void **hostptr, int *p_size) {
  cudaError_t err;

  err = cudaHostAlloc((void **)hostptr, (*p_size)*sizeof(double), cudaHostAllocMapped);
  if (err != cudaSuccess) fprintf(stderr, "Error in cudaHostAlloc: %s\n", cudaGetErrorString( err ));
}
