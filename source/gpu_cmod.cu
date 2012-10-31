#include <stdio.h>
#include <math.h>

#define KM 42
#define NSTREAMS 42

#define CUDA_CHECK_ERROR(errorMessage) do {                                 \
    cudaError_t err = cudaGetLastError();                                    \
    if( cudaSuccess != err) {                                                \
        fprintf(stderr, "Cuda error: %s in file '%s' in line %i : %s.\n",    \
                errorMessage, __FILE__, __LINE__, cudaGetErrorString( err) );\
        exit(EXIT_FAILURE);                                                  \
    }                                                                        \
    err = cudaThreadSynchronize();                                           \
    if( cudaSuccess != err) {                                                \
        fprintf(stderr, "Cuda error: %s in file '%s' in line %i : %s.\n",    \
                errorMessage, __FILE__, __LINE__, cudaGetErrorString( err) );\
        exit(EXIT_FAILURE);                                                  \
    } } while (0)


//declared in extern C block otherwise C++ compilers mangle the function names
extern "C" {

void cuda_init();

void my_cudamallochost(void **hostptr, int* size);

void cuda_state_initialize(double *constants, double *pressz,
        double *tmin, double *tmax, double *smin, double *smax);



}




//functions


int cuda_initialized = 0;

//Fortran entry for initializing CUDA
void cuda_init() {
  if (cuda_initialized == 0) {
    cuda_initialized = 1;

    cudaSetDeviceFlags(cudaDeviceMapHost);
    cudaSetDevice(0);
  }
}


//Fortran entry for allocating pinned memory
void my_cudamallochost(void **hostptr, int *p_size) {
  cudaError_t err;

//  err = cudaMallocHost((void **)hostptr, (*p_size)*sizeof(double));
//  if (err != cudaSuccess) fprintf(stderr, "Error in cudaMallocHost: %s\n", cudaGetErrorString( err ));

  err = cudaHostAlloc((void **)hostptr, (*p_size)*sizeof(double), cudaHostAllocMapped);
  if (err != cudaSuccess) fprintf(stderr, "Error in cudaHostAlloc: %s\n", cudaGetErrorString( err ));

}

//Fortran entry for initializing constants used in state computations
//constants passed as parameters
double  mwjfnums0t1, mwjfnums0t3, mwjfnums1t1, mwjfnums2t0, mwjfdens0t2,
        mwjfdens0t4, mwjfdens1t0, mwjfdens1t1, mwjfdens1t3, mwjfdensqt0,
        mwjfdensqt2;

//all constants needed on the GPU
__constant__ double d_tmax[KM];
__constant__ double d_tmin[KM];
__constant__ double d_smax[KM];
__constant__ double d_smin[KM];

__constant__ double d_mwjfnums0t0[KM];
__constant__ double d_mwjfnums0t2[KM];
__constant__ double d_mwjfnums1t0[KM];

__constant__ double d_mwjfdens0t0[KM]; __constant__ double d_mwjfdens0t1[KM]; __constant__ double d_mwjfdens0t3[KM];

//declare streams
int cuda_state_initialized = 0;
cudaStream_t stream[NSTREAMS];

void cuda_state_initialize(double *constants, double *pressz,
        double *tmin, double *tmax, double *smin, double *smax) {

  if (cuda_state_initialized == 0) {
  cuda_state_initialized = 1;
  
  //Perhaps we need to find some better way to pass constants around
  //for now we do this because we want the GPU/C code to use the exact same
  //values as the Fortran code even if both parts are compiled with
  //different compilers
  mwjfnums0t1 = constants[21]; //= mwjfnp0s0t1;
  mwjfnums0t3 = constants[23]; //= mwjfnp0s0t3;
  mwjfnums1t1 = constants[25]; //= mwjfnp0s1t1;
  mwjfnums2t0 = constants[26]; //= mwjfnp0s2t0;
  mwjfdens0t2 = constants[34]; //= mwjfdp0s0t2;
  mwjfdens0t4 = constants[36]; //= mwjfdp0s0t4;
  mwjfdens1t0 = constants[37]; //= mwjfdp0s1t0;
  mwjfdens1t1 = constants[38]; //= mwjfdp0s1t1;
  mwjfdens1t3 = constants[39]; //= mwjfdp0s1t3;
  mwjfdensqt0 = constants[40]; //= mwjfdp0sqt0;
  mwjfdensqt2 = constants[41]; //= mwjfdp0sqt2;

  double mwjfnp0s0t0 = constants[20];
  double mwjfnp0s0t2 = constants[22];
  double mwjfnp0s1t0 = constants[24];
  double mwjfnp1s0t0 = constants[27];
  double mwjfnp1s0t2 = constants[28];
  double mwjfnp1s1t0 = constants[29];
  double mwjfnp2s0t0 = constants[30];
  double mwjfnp2s0t2 = constants[31];

  double mwjfdp0s0t0 = constants[32];
  double mwjfdp0s0t1 = constants[33];
  double mwjfdp0s0t3 = constants[35];
  double mwjfdp1s0t0 = constants[42];
  double mwjfdp2s0t3 = constants[43];
  double mwjfdp3s0t1 = constants[44];

  double p;

  //initialize all constant arrays to be stored in constant memory on the GPU
  double h_mwjfnums0t0[KM];
  double h_mwjfnums0t2[KM];
  double h_mwjfnums1t0[KM];
  double h_mwjfdens0t0[KM];
  double h_mwjfdens0t1[KM];
  double h_mwjfdens0t3[KM];

  int k;
  for (k=0; k<KM; k++) {
      p = 10.0*pressz[k];

      // first calculate numerator of MWJF density [P_1(S,T,p)]
      h_mwjfnums0t0[k] = mwjfnp0s0t0 + p*(mwjfnp1s0t0 + p*mwjfnp2s0t0);
      h_mwjfnums0t2[k] = mwjfnp0s0t2 + p*(mwjfnp1s0t2 + p*mwjfnp2s0t2);
      h_mwjfnums1t0[k] = mwjfnp0s1t0 + p*mwjfnp1s1t0;

      // now calculate denominator of MWJF density [P_2(S,T,p)]
      h_mwjfdens0t0[k] = mwjfdp0s0t0 + p*mwjfdp1s0t0;
      h_mwjfdens0t1[k] = mwjfdp0s0t1 + (p*p*p) * mwjfdp3s0t1;   //used to be p**3 in FORTRAN
      h_mwjfdens0t3[k] = mwjfdp0s0t3 + (p*p) * mwjfdp2s0t3;     //used to be p**2 in FORTRAN
  }

  //bunch of memcpy to symbols go here
  cudaMemcpyToSymbol("d_tmax", tmax, KM*sizeof(double), 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol("d_tmin", tmin, KM*sizeof(double), 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol("d_smax", smax, KM*sizeof(double), 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol("d_smin", smin, KM*sizeof(double), 0, cudaMemcpyHostToDevice);

  cudaMemcpyToSymbol("d_mwjfnums0t0", h_mwjfnums0t0, KM*sizeof(double), 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol("d_mwjfnums0t2", h_mwjfnums0t2, KM*sizeof(double), 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol("d_mwjfnums1t0", h_mwjfnums1t0, KM*sizeof(double), 0, cudaMemcpyHostToDevice);

  cudaMemcpyToSymbol("d_mwjfdens0t0", h_mwjfdens0t0, KM*sizeof(double), 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol("d_mwjfdens0t1", h_mwjfdens0t1, KM*sizeof(double), 0, cudaMemcpyHostToDevice);
  cudaMemcpyToSymbol("d_mwjfdens0t3", h_mwjfdens0t3, KM*sizeof(double), 0, cudaMemcpyHostToDevice);

  //error checking
  cudaDeviceSynchronize();
  CUDA_CHECK_ERROR("After cudaMemcpyToSymbols");

  //setup streams
  for (k=0; k<NSTREAMS; k++) {
    cudaStreamCreate(&stream[k]);
  }


  }
}






