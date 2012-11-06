#include <stdio.h>
#include <math.h>

//need to find out how to organize the compile time setting of the domain size
//such that the size is only set in one location
#define KM 42
#define NX_BLOCK 904
#define NY_BLOCK 604
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


//declared in extern C block to prevent C++ compilers from mangling function names
extern "C" {

void cuda_init();

void my_cudamallochost(void **hostptr, int* size);

void cuda_state_initialize(double *constants, double *pressz,
        double *tmin, double *tmax, double *smin, double *smax);

//specific functions
void mwjf_state_gpu(double *TEMPK, double *SALTK,
        		double *RHOOUT, double *DRHODT, double *DRHODS, 
        		int *pn_outputs, int *pstart_k, int *pend_k);

__global__ void mwjf_state_1D(double *TEMPK, double *SALTK,
		double *RHOFULL, double *DRHODT, double *DRHODS,
		int n_outputs, int start_k, int end_k);

void gpu_compare (double *a1, double *a2, int N);


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
//all constants needed on the GPU
__constant__ double d_mwjfnums0t1;
__constant__ double d_mwjfnums0t3;
__constant__ double d_mwjfnums1t1;
__constant__ double d_mwjfnums2t0;
__constant__ double d_mwjfdens0t2;
__constant__ double d_mwjfdens0t4;
__constant__ double d_mwjfdens1t0;
__constant__ double d_mwjfdens1t1;
__constant__ double d_mwjfdens1t3;
__constant__ double d_mwjfdensqt0;
__constant__ double d_mwjfdensqt2;

__constant__ double d_tmax[KM];
__constant__ double d_tmin[KM];
__constant__ double d_smax[KM];
__constant__ double d_smin[KM];

__constant__ double d_mwjfnums0t0[KM];
__constant__ double d_mwjfnums0t2[KM];
__constant__ double d_mwjfnums1t0[KM];

__constant__ double d_mwjfdens0t0[KM];
__constant__ double d_mwjfdens0t1[KM];
__constant__ double d_mwjfdens0t3[KM];

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
  cudaMemcpyToSymbol("d_mwjfnums0t1", &constants[21], sizeof(double), 0, cudaMemcpyHostToDevice); //= mwjfnp0s0t1
  cudaMemcpyToSymbol("d_mwjfnums0t3", &constants[23], sizeof(double), 0, cudaMemcpyHostToDevice); //= mwjfnp0s0t3
  cudaMemcpyToSymbol("d_mwjfnums1t1", &constants[25], sizeof(double), 0, cudaMemcpyHostToDevice); //= mwjfnp0s1t1
  cudaMemcpyToSymbol("d_mwjfnums2t0", &constants[26], sizeof(double), 0, cudaMemcpyHostToDevice); //= mwjfnp0s2t0
  cudaMemcpyToSymbol("d_mwjfdens0t2", &constants[34], sizeof(double), 0, cudaMemcpyHostToDevice); //= mwjfdp0s0t2
  cudaMemcpyToSymbol("d_mwjfdens0t4", &constants[36], sizeof(double), 0, cudaMemcpyHostToDevice); //= mwjfdp0s0t4
  cudaMemcpyToSymbol("d_mwjfdens1t0", &constants[37], sizeof(double), 0, cudaMemcpyHostToDevice); //= mwjfdp0s1t0
  cudaMemcpyToSymbol("d_mwjfdens1t1", &constants[38], sizeof(double), 0, cudaMemcpyHostToDevice); //= mwjfdp0s1t1
  cudaMemcpyToSymbol("d_mwjfdens1t3", &constants[39], sizeof(double), 0, cudaMemcpyHostToDevice); //= mwjfdp0s1t3
  cudaMemcpyToSymbol("d_mwjfdensqt0", &constants[40], sizeof(double), 0, cudaMemcpyHostToDevice); //= mwjfdp0sqt0
  cudaMemcpyToSymbol("d_mwjfdensqt2", &constants[41], sizeof(double), 0, cudaMemcpyHostToDevice); //= mwjfdp0sqt2

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
  double c10 = constants[7];

  //initialize all constant arrays to be stored in constant memory on the GPU
  double h_mwjfnums0t0[KM];
  double h_mwjfnums0t2[KM];
  double h_mwjfnums1t0[KM];
  double h_mwjfdens0t0[KM];
  double h_mwjfdens0t1[KM];
  double h_mwjfdens0t3[KM];

  int k;
  for (k=0; k<KM; k++) {
      p = c10*pressz[k];

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


void mwjf_state_gpu(double *TEMPK, double *SALTK, 
        		double *RHOOUT, double *DRHODT, double *DRHODS,
        		int *pn_outputs, int *pstart_k, int *pend_k) {
  int n_outputs = *pn_outputs;
  int start_k = *pstart_k-1;
  int end_k = *pend_k; //no -1 here as we're going from including to excluding
  //cudaError_t err;
  
  //execution parameters
  dim3 threads(256,1);
  dim3 grid(1,1);
  grid.x = (int)ceilf(((float)(NX_BLOCK*NY_BLOCK) / (float)threads.x));
  grid.y = (KM);
  
  
  // perhaps not needed on Fermi GPUs, who knew?
  //corresponding device pointers
  double *d_SALTK;
  double *d_TEMPK;

  double *d_RHOOUT;
  double *d_DRHODT = NULL;
  double *d_DRHODS = NULL;
  
  //obtain device pointers for host mapped memory
  err = cudaHostGetDevicePointer((void**) &d_TEMPK, TEMPK, 0);
  if (err != cudaSuccess) fprintf(stderr, "Error retrieving device pointer: %s\n", cudaGetErrorString( err ));
  err = cudaHostGetDevicePointer((void**) &d_SALTK, SALTK, 0);
  if (err != cudaSuccess) fprintf(stderr, "Error retrieving device pointer: %s\n", cudaGetErrorString( err ));
  err = cudaHostGetDevicePointer((void**) &d_RHOOUT, RHOOUT, 0);
  if (err != cudaSuccess) fprintf(stderr, "Error retrieving device pointer: %s\n", cudaGetErrorString( err ));
  
  if (n_outputs == 3) {
    err = cudaHostGetDevicePointer((void**) &d_DRHODT, DRHODT, 0);
    if (err != cudaSuccess) fprintf(stderr, "Error retrieving device pointer: %s\n", cudaGetErrorString( err ));
    err = cudaHostGetDevicePointer((void**) &d_DRHODS, DRHODS, 0);
    if (err != cudaSuccess) fprintf(stderr, "Error retrieving device pointer: %s\n", cudaGetErrorString( err ));
  }
  //
  
  //this synchronize is a bit over-protective but currently left in for debugging purposes
  cudaDeviceSynchronize();
  CUDA_CHECK_ERROR("Before mwjf_state_1D kernel execution");
  
  mwjf_state_1D<<<grid,threads,0,stream[1]>>>(TEMPK, SALTK, RHOOUT, DRHODT, DRHODS,
        n_outputs, start_k, end_k);
  
  
  //synchronize because we currently don't know when inputs or outputs will be used by CPU
  //the more this sync can be delayed the more overlap with CPU execution can be exploited
  cudaDeviceSynchronize();
  CUDA_CHECK_ERROR("After mwjf_state_1D kernel execution");
  
}



/*
 * The idea behind this kernel is to create a thread for each element in the array and not just for one level.
 * This eliminates the need to have a loop and is also cacheline boundary oblivious, making it a perfect for
 * using device mapped host memory.
 */
__global__ void mwjf_state_1D(double *TEMPK, double *SALTK, 
		double *RHOOUT, double *DRHODT, double *DRHODS,
		int n_outputs, int start_k, int end_k) {

  //obtain global ids
  //int j = threadIdx.y + blockIdx.y * BLOCK_Y;
  //int i = threadIdx.x + blockIdx.x * BLOCK_X;
  //obtain global id
  int i = blockIdx.y * gridDim.x * blockDim.x + blockIdx.x * blockDim.x + threadIdx.x;
  int k = start_k + (i / NX_BLOCK*NY_BLOCK);
  //obtain array index
  int index = i + start_k*NX_BLOCK*NY_BLOCK;

  double tq, sq, sqr, work1, work2, work3, work4, denomk;

//emulating
//for (j=0; j< NY_BLOCK; j++) {
//for (i=0; i< NX_BLOCK; i++) {
//emulating
//  if (j < NY_BLOCK && i < NX_BLOCK) {
  if (i < NX_BLOCK*NY_BLOCK*(end_k-start_k)) {

//unrolled for (k=start_k; k < end_k; k++)

        tq = min(TEMPK[index],d_tmax[k]);
        tq = max(tq,d_tmin[k]);

        sq = min(SALTK[index],d_smax[k]);
        sq = 1000.0 * max(sq,d_smin[k]);

        sqr = sqrt(sq); //double precision sqrt round towards zero

        work1 = d_mwjfnums0t0[k] + tq * (d_mwjfnums0t1 + tq * (d_mwjfnums0t2[k] + d_mwjfnums0t3 * tq)) +
                              sq * (d_mwjfnums1t0[k] + d_mwjfnums1t1 * tq + d_mwjfnums2t0 * sq);

        work2 = d_mwjfdens0t0[k] + tq * (d_mwjfdens0t1[k] + tq * (d_mwjfdens0t2 +
           tq * (d_mwjfdens0t3[k] + d_mwjfdens0t4 * tq))) +
           sq * (d_mwjfdens1t0 + tq * (d_mwjfdens1t1 + tq*tq*d_mwjfdens1t3)+
           sqr * (d_mwjfdensqt0 + tq*tq*d_mwjfdensqt2));

        denomk = 1.0/work2;
//      if (present(RHOFULL)) then
        RHOOUT[index] = work1*denomk;
//      endif

        if (n_outputs == 3) { 
        
	//      if (present(DRHODT)) then
	        work3 = // dP_1/dT
	                 d_mwjfnums0t1 + tq * (2.0*d_mwjfnums0t2[k] +
	                 3.0*d_mwjfnums0t3 * tq) + d_mwjfnums1t1 * sq;
	
	        work4 = // dP_2/dT
	                 d_mwjfdens0t1[k] + sq * d_mwjfdens1t1 +
	                 tq * (2.0*(d_mwjfdens0t2 + sq*sqr*d_mwjfdensqt2) +
	                 tq * (3.0*(d_mwjfdens0t3[k] + sq * d_mwjfdens1t3) +
	                 tq *  4.0*d_mwjfdens0t4));
	
	        DRHODT[index] = (work3 - work1*denomk*work4)*denomk;
	
	//      endif
	//      if (present(DRHODS)) then
	        work3 = // dP_1/dS
	                 d_mwjfnums1t0[k] + d_mwjfnums1t1 * tq + 2.0*d_mwjfnums2t0 * sq;
	
	        work4 = d_mwjfdens1t0 +   // dP_2/dS
	                 tq * (d_mwjfdens1t1 + tq*tq*d_mwjfdens1t3) +
	                 1.5*sqr*(d_mwjfdensqt0 + tq*tq*d_mwjfdensqt2);
	
	        DRHODS[index] = (work3 - work1*denomk*work4)*denomk * 1000.0;
	//      endif

        } //end if n_outputs == 3


  } // end of if-statement

}




void gpu_compare (double *a1, double *a2, int N) {
  int i,res = 0;
  int print = 0;
  int zeros = 0;
  double eps = 0.0000001;

  for (i=0; i<N; i++) {
//    if (i<1840080 && i>1840075) { printf("values at i=%d, a1= %20.17e, a2= %20.17e\n", i, a1[i], a2[i]); }

    if (a1[i] < eps && a1[i] > -eps) {
      zeros++;
    }

    double diff = a1[i]-a2[i];
    if (diff > eps || diff < -eps) {
        res++;
        if (print < 10) {
          print++;
          printf("Error detected at i=%d, a1= %20.17e a2= %20.17e\n",i,a1[i],a2[i]);
        }
    }

  }

  if (zeros > N/2) {
    fprintf(stderr, "Error: more than 50% of array 1 contains zeros\n");
  }

  fprintf(stdout,"Number of errors in GPU result: %d\n",res);

}



