#include <stdio.h>
#include <math.h>

#include "gpu_domain.h"

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
        double *tmin, double *tmax, double *smin, double *smax, int *pmy_task, int *kmt);

//specific functions
void mwjf_state_gpu(double *TEMPK, double *SALTK,
        		double *RHOOUT, double *DRHODT, double *DRHODS, 
        		int *pn_outputs, int *pstart_k, int *pend_k);

void mwjf_statepd_gpu(double *TEMPK, double *SALTK,
        		double *RHOOUT, int *pstart_k, int *pend_k);

__global__ void mwjf_state_1D(double *TEMPK, double *SALTK,
		double *RHOFULL, double *DRHODT, double *DRHODS,
		int n_outputs, int start_k, int end_k);

__global__ void mwjf_statepd_1D(double *TEMPK, double *SALTK, 
		double *RHOOUT, int start_k, int end_k);

void gpu_compare(double *a1, double *a2, int *pN, int *pName);

void buoydiff_gpu(double *DBLOC, double *DBSFC, double *TEMP, double *SALT);

__global__ void buoydiff_kernel1D(double *DBLOC, double *DBSFC, double *TEMP, double *SALT, int *KMT, int start_k, int end_k);

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
  if (err != cudaSuccess) fprintf(stdout, "Error in cudaHostAlloc: %s\n", cudaGetErrorString( err ));
 
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

__constant__ double d_grav;

//__constant__ double d_tmax[KM];
//__constant__ double d_tmin[KM];
//__constant__ double d_smax[KM];
//__constant__ double d_smin[KM];

__constant__ double d_mwjfnums0t0[KM];
__constant__ double d_mwjfnums0t2[KM];
__constant__ double d_mwjfnums1t0[KM];

__constant__ double d_mwjfdens0t0[KM];
__constant__ double d_mwjfdens0t1[KM];
__constant__ double d_mwjfdens0t3[KM];

//device pointer for kmt in GPU global memory
int *d_kmt;

//declare streams
int cuda_state_initialized = 0;
cudaStream_t stream[NSTREAMS];

//debugging
int my_task;

void cuda_state_initialize(double *constants, double *pressz,
        double *tmin, double *tmax, double *smin, double *smax, int *pmy_task, int *kmt) {
  cudaError_t err;
	
  if (cuda_state_initialized == 0) {
  cuda_state_initialized = 1;
  
  //Perhaps we need to find some better way to pass constants around
  //for now we do this because we want the GPU/C code to use the exact same
  //values as the Fortran code even if both parts are compiled with
  //different compilers
  /* CUDA 4 style
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

    cudaMemcpyToSymbol("d_grav", &constants[45], sizeof(double), 0, cudaMemcpyHostToDevice); //= grav
  */

  //CUDA 5.0 style
  err = cudaMemcpyToSymbol(d_mwjfnums0t1, &constants[21], sizeof(double), 0, cudaMemcpyHostToDevice); //= mwjfnp0s0t1
  if (err != cudaSuccess) fprintf(stderr, "Error doing cudaMemcpyToSymbol d_mwjfnums0t1\n");
  err = cudaMemcpyToSymbol(d_mwjfnums0t3, &constants[23], sizeof(double), 0, cudaMemcpyHostToDevice); //= mwjfnp0s0t3
  if (err != cudaSuccess) fprintf(stderr, "Error doing cudaMemcpyToSymbol d_mwjfnums0t3\n");
  err = cudaMemcpyToSymbol(d_mwjfnums1t1, &constants[25], sizeof(double), 0, cudaMemcpyHostToDevice); //= mwjfnp0s1t1
  if (err != cudaSuccess) fprintf(stderr, "Error doing cudaMemcpyToSymbol d_mwjfnums1t1\n");
  err = cudaMemcpyToSymbol(d_mwjfnums2t0, &constants[26], sizeof(double), 0, cudaMemcpyHostToDevice); //= mwjfnp0s2t0
  if (err != cudaSuccess) fprintf(stderr, "Error doing cudaMemcpyToSymbol d_mwjfnums2t0\n");
  err = cudaMemcpyToSymbol(d_mwjfdens0t2, &constants[34], sizeof(double), 0, cudaMemcpyHostToDevice); //= mwjfdp0s0t2
  if (err != cudaSuccess) fprintf(stderr, "Error doing cudaMemcpyToSymbol d_mwjfdens0t2\n");
  err = cudaMemcpyToSymbol(d_mwjfdens0t4, &constants[36], sizeof(double), 0, cudaMemcpyHostToDevice); //= mwjfdp0s0t4
  if (err != cudaSuccess) fprintf(stderr, "Error doing cudaMemcpyToSymbol d_mwjfdens0t4\n");
  err = cudaMemcpyToSymbol(d_mwjfdens1t0, &constants[37], sizeof(double), 0, cudaMemcpyHostToDevice); //= mwjfdp0s1t0
  if (err != cudaSuccess) fprintf(stderr, "Error doing cudaMemcpyToSymbol d_mwjfdens1t0\n");
  err = cudaMemcpyToSymbol(d_mwjfdens1t1, &constants[38], sizeof(double), 0, cudaMemcpyHostToDevice); //= mwjfdp0s1t1
  if (err != cudaSuccess) fprintf(stderr, "Error doing cudaMemcpyToSymbol d_mwjfdens1t1\n");
  err = cudaMemcpyToSymbol(d_mwjfdens1t3, &constants[39], sizeof(double), 0, cudaMemcpyHostToDevice); //= mwjfdp0s1t3
  if (err != cudaSuccess) fprintf(stderr, "Error doing cudaMemcpyToSymbol d_mwjfdens1t3\n");
  err = cudaMemcpyToSymbol(d_mwjfdensqt0, &constants[40], sizeof(double), 0, cudaMemcpyHostToDevice); //= mwjfdp0sqt0
  if (err != cudaSuccess) fprintf(stderr, "Error doing cudaMemcpyToSymbol d_mwjfdensqt0\n");
  err = cudaMemcpyToSymbol(d_mwjfdensqt2, &constants[41], sizeof(double), 0, cudaMemcpyHostToDevice); //= mwjfdp0sqt2
  if (err != cudaSuccess) fprintf(stderr, "Error doing cudaMemcpyToSymbol d_mwjfdensqt2\n");

  err = cudaMemcpyToSymbol(d_grav, &constants[45], sizeof(double), 0, cudaMemcpyHostToDevice); //= grav
  if (err != cudaSuccess) fprintf(stderr, "Error doing cudaMemcpyToSymbol d_grav\n");

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

//no longer used
//  cudaMemcpyToSymbol("d_tmax", tmax, KM*sizeof(double), 0, cudaMemcpyHostToDevice);
//  cudaMemcpyToSymbol("d_tmin", tmin, KM*sizeof(double), 0, cudaMemcpyHostToDevice);
//  cudaMemcpyToSymbol("d_smax", smax, KM*sizeof(double), 0, cudaMemcpyHostToDevice);
//  cudaMemcpyToSymbol("d_smin", smin, KM*sizeof(double), 0, cudaMemcpyHostToDevice);

  //bunch of memcpy to symbols go here
  err = cudaMemcpyToSymbol(d_mwjfnums0t0, h_mwjfnums0t0, KM*sizeof(double), 0, cudaMemcpyHostToDevice);
  if (err != cudaSuccess) fprintf(stderr, "Error doing cudaMemcpyToSymbol d_mwjfnums0t0\n");
  err = cudaMemcpyToSymbol(d_mwjfnums0t2, h_mwjfnums0t2, KM*sizeof(double), 0, cudaMemcpyHostToDevice);
  if (err != cudaSuccess) fprintf(stderr, "Error doing cudaMemcpyToSymbol d_mwjfnums0t2\n");
  err = cudaMemcpyToSymbol(d_mwjfnums1t0, h_mwjfnums1t0, KM*sizeof(double), 0, cudaMemcpyHostToDevice);
  if (err != cudaSuccess) fprintf(stderr, "Error doing cudaMemcpyToSymbol d_mwjfnums1t0\n");

  err = cudaMemcpyToSymbol(d_mwjfdens0t0, h_mwjfdens0t0, KM*sizeof(double), 0, cudaMemcpyHostToDevice);
  if (err != cudaSuccess) fprintf(stderr, "Error doing cudaMemcpyToSymbol d_mwjfdens0t0\n");
  err = cudaMemcpyToSymbol(d_mwjfdens0t1, h_mwjfdens0t1, KM*sizeof(double), 0, cudaMemcpyHostToDevice);
  if (err != cudaSuccess) fprintf(stderr, "Error doing cudaMemcpyToSymbol d_mwjfdens0t1\n");
  err = cudaMemcpyToSymbol(d_mwjfdens0t3, h_mwjfdens0t3, KM*sizeof(double), 0, cudaMemcpyHostToDevice);
  if (err != cudaSuccess) fprintf(stderr, "Error doing cudaMemcpyToSymbol d_mwjfdens0t3\n");

  err = cudaMalloc(&d_kmt, NX_BLOCK*NY_BLOCK*sizeof(int));
  if (err != cudaSuccess) fprintf(stderr, "Error doing cudaMalloc d_kmt\n");
  err = cudaMemcpy(d_kmt, kmt, NX_BLOCK*NY_BLOCK*sizeof(int), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) fprintf(stderr, "Error doing cudaMemcpyHostToDevice KMT\n");

  
  //error checking
  cudaDeviceSynchronize();
  CUDA_CHECK_ERROR("After cudaMemcpyToSymbols");

  //setup streams
  for (k=0; k<NSTREAMS; k++) {
    cudaStreamCreate(&stream[k]);
  }
  
  //error checking
  cudaDeviceSynchronize();
  CUDA_CHECK_ERROR("After cudaStreamCreate");
  
  my_task = *pmy_task;
//debugging
//  if (my_task == 0) {
//    printf("GPU_CMOD using constants:\n");
//    printf("c0=%20.17e, c1=%20.17e, c2=%20.17e, c3=%20.17e, c4=%20.17e, c5=%20.17e, c8=%20.17e, c10=%20.17e\n",
//    		constants[0], constants[1], constants[2], constants[3],
//    		constants[4], constants[5], constants[6], constants[7]);
//    printf("c16=%20.17e, c1000=%20.17e, c10000=%20.17e, c1p5=%20.17e, p33=%20.17e, p5=%20.17e, p25=%20.17e, p125=%20.17e, p001=%20.17e\n",
//    		constants[8], constants[9], constants[10], constants[11],
//    		constants[12], constants[13], constants[14], constants[15], constants[16]);
//    
//    
//  }

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
  grid.y = (end_k-start_k);
  
//  if (my_task == 0)
//  printf("n_outputs=%d, start_k=%d, end_k=%d, tx=%d, ty=%d, gx=%d, gy=%d\n", 
//  		  n_outputs, start_k, end_k, threads.x, threads.y, grid.x, grid.y);
  
  //zero output array, for debugging purposes only
  //memset(RHOOUT, 0, NX_BLOCK*NY_BLOCK*KM*sizeof(double));

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

void mwjf_statepd_gpu(double *TEMPK, double *SALTK, 
        		double *RHOOUT, int *pstart_k, int *pend_k) {
  int start_k = *pstart_k-1;
  int end_k = *pend_k; //no -1 here as we're going from including to excluding
  //cudaError_t err;
  
  //execution parameters
  dim3 threads(256,1);
  dim3 grid(1,1);
  grid.x = (int)ceilf(((float)(NX_BLOCK*NY_BLOCK) / (float)threads.x));
  grid.y = (end_k-start_k);
  
  mwjf_statepd_1D<<<grid,threads,0,stream[1]>>>(TEMPK, SALTK, RHOOUT, start_k, end_k);
  
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
  int k = start_k + (i /  (NX_BLOCK*NY_BLOCK));
  //obtain array index
  int index = i + start_k*(NX_BLOCK*NY_BLOCK);

  double tq, sq, sqr, work1, work2, work3, work4, denomk;

//emulating
//for (j=0; j< NY_BLOCK; j++) {
//for (i=0; i< NX_BLOCK; i++) {
//emulating
//  if (j < NY_BLOCK && i < NX_BLOCK) {
  if (i < (NX_BLOCK*NY_BLOCK)*(end_k-start_k)) {

//unrolled for (k=start_k; k < end_k; k++)

	    //tmax, tmin, smax, smin not really used in MWJF, replace with -2 and 999
        tq = min(TEMPK[index], 999.0);	//d_tmax[k]
        tq = max(tq, -2.0);				//d_tmin[k]

        sq = min(SALTK[index], 0.999);	//d_smax[k]
        sq = 1000.0 * max(sq, 0.0);		//d_smin[k]

        sqr = sqrt(sq);

        work1 = d_mwjfnums0t0[k] + tq * (d_mwjfnums0t1 + tq * (d_mwjfnums0t2[k] + d_mwjfnums0t3 * tq)) +
                              sq * (d_mwjfnums1t0[k] + d_mwjfnums1t1 * tq + d_mwjfnums2t0 * sq);

        work2 = d_mwjfdens0t0[k] + tq * (d_mwjfdens0t1[k] + tq * (d_mwjfdens0t2 +
           tq * (d_mwjfdens0t3[k] + d_mwjfdens0t4 * tq))) +
           sq * (d_mwjfdens1t0 + tq * (d_mwjfdens1t1 + tq*tq*d_mwjfdens1t3)+
           sqr * (d_mwjfdensqt0 + tq*tq*d_mwjfdensqt2));

        denomk = 1.0/work2;
//      if (present(RHOFULL)) then
        RHOOUT[index] = work1*denomk;
        //RHOOUT[index] = 1337.0; //for debugging obviously
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

/*
 * The idea behind this kernel is to create a thread for each element in the array and not just for one level.
 * This eliminates the need to have a loop and is also cacheline boundary oblivious, making it a perfect for
 * using device mapped host memory.
 */
__global__ void mwjf_statepd_1D(double *TEMPK, double *SALTK, 
		double *RHOOUT, int start_k, int end_k) {

  //obtain global id
  int i = blockIdx.y * gridDim.x * blockDim.x + blockIdx.x * blockDim.x + threadIdx.x;
  //obtain array index
  int index = i + start_k*(NX_BLOCK*NY_BLOCK);

  double tq, sq, sqr, work1, work2, denomk;

  if (i < (NX_BLOCK*NY_BLOCK)*(end_k-start_k)) {

	    //tmax, tmin, smax, smin not really used in MWJF, replace with -2 and 999
        tq = min(TEMPK[index], 999.0);	//d_tmax[k]
        tq = max(tq, -2.0);				//d_tmin[k]

        sq = min(SALTK[index], 0.999);	//d_smax[k]
        sq = 1000.0 * max(sq, 0.0);		//d_smin[k]

        sqr = sqrt(sq);

        work1 = d_mwjfnums0t0[0] + tq * (d_mwjfnums0t1 + tq * (d_mwjfnums0t2[0] + d_mwjfnums0t3 * tq)) +
                              sq * (d_mwjfnums1t0[0] + d_mwjfnums1t1 * tq + d_mwjfnums2t0 * sq);

        work2 = d_mwjfdens0t0[0] + tq * (d_mwjfdens0t1[0] + tq * (d_mwjfdens0t2 +
           tq * (d_mwjfdens0t3[0] + d_mwjfdens0t4 * tq))) +
           sq * (d_mwjfdens1t0 + tq * (d_mwjfdens1t1 + tq*tq*d_mwjfdens1t3)+
           sqr * (d_mwjfdensqt0 + tq*tq*d_mwjfdensqt2));

        denomk = 1.0/work2;
        RHOOUT[index] = work1*denomk;

  } // end of if-statement

}

const char *var_names[6] = { "ERROR", "RHOOUT", "DBLOC", "DBSFC", "STEPMOD_RHO", "ADVT_PD" };


void gpu_compare (double *a1, double *a2, int *pN, int *pName) {
  int N = *pN;
  int vName = *pName;
  int i=0, res=0;
  int print = 0;
  int zero_one = 0;
  int zero_two = 0;
  double eps = 0.00000000001;

  for (i=0; i<N; i++) {

    if (a1[i] < eps && a1[i] > -eps) { zero_one++; }
    if (a2[i] < eps && a2[i] > -eps) { zero_two++; }

    if (isnan(a1[i]) || isnan(a2[i])) {
        res++;
        if (print < 10) {
          print++;
          fprintf(stderr, "Node %d: Error detected at i=%d, a1= %20.17e a2= %20.17e\n",my_task,i,a1[i],a2[i]);
        }
    }

    double diff = a1[i]-a2[i];
    if (diff > eps || diff < -eps) {
        res++;
        if (print < 10) {
          print++;
          fprintf(stderr, "Node %d: Error detected at i=%d, \t a1= \t %20.17e \t a2= \t %20.17e\n",my_task,i,a1[i],a2[i]);
        }
    }

  }

  if (zero_one > 9*(N/10)) { fprintf(stderr, "Node %d: Error: array1 contains %d zeros\n",my_task, zero_one); }
  if (zero_two > 9*(N/10)) { fprintf(stderr, "Node %d: Error: array2 contains %d zeros\n",my_task, zero_two); }

  if (zero_one != zero_two) {
    fprintf(stderr, "Node %d: Error: number of zeros in arrays dont correspond zero1=%d, zero2=%d\n",my_task, zero_one, zero_two);
  }

  if (res > 0) {
    if (vName == 0) {
      fprintf(stdout,"Node %d: Number of errors in GPU result: %d\n",my_task,res);
    } else {
	  fprintf(stdout,"Node %d: Number of errors in %s GPU result: %d\n",my_task,var_names[vName],res);
    }
  }
}

double *d_TEMP;
double *d_SALT;

void buoydiff_gpu(double *DBLOC, double *DBSFC, double *TEMP, double *SALT) {
    cudaError_t err;
    
    //completely unnecessary but testing anyway
//    double *d_DBLOC;
//    double *d_DBSFC;    
//    err = cudaHostGetDevicePointer((void**)&d_DBLOC, DBLOC, 0);
//    if (err != cudaSuccess) fprintf(stderr, "Error retrieving device pointer: %s\n", cudaGetErrorString( err ));
//    err = cudaHostGetDevicePointer((void**)&d_DBSFC, DBSFC, 0);
//    if (err != cudaSuccess) fprintf(stderr, "Error retrieving device pointer: %s\n", cudaGetErrorString( err ));

    
    //allocate space and copy TRCR to GPU
    //this will later be reused by other GPU kernels in vmix_kpp
    err = cudaMalloc((void **)&d_TEMP, NX_BLOCK*NY_BLOCK*KM*sizeof(double));
    if (err != cudaSuccess) fprintf(stderr, "Error in cudaMalloc d_TRCR %s\n", cudaGetErrorString( err ));
    err = cudaMemcpy(d_TEMP, TEMP, NX_BLOCK*NY_BLOCK*KM*sizeof(double), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) fprintf(stderr, "Error in cudaMemcpy host to device TRCR: %s\n", cudaGetErrorString( err ));
    
    err = cudaMalloc((void **)&d_SALT, NX_BLOCK*NY_BLOCK*KM*sizeof(double));
    if (err != cudaSuccess) fprintf(stderr, "Error in cudaMalloc d_TRCR %s\n", cudaGetErrorString( err ));
    err = cudaMemcpy(d_SALT, SALT, NX_BLOCK*NY_BLOCK*KM*sizeof(double), cudaMemcpyHostToDevice);
    if (err != cudaSuccess) fprintf(stderr, "Error in cudaMemcpy host to device TRCR: %s\n", cudaGetErrorString( err ));

    memset(DBLOC, 0, NX_BLOCK*NY_BLOCK*KM*sizeof(double));
    memset(DBSFC, 0, NX_BLOCK*NY_BLOCK*KM*sizeof(double));
    
    //setup execution parameters
    dim3 threads(256,1);
    dim3 grid(1,1);
    grid.x = (int)ceilf(((float)(NX_BLOCK*NY_BLOCK) / (float)threads.x));
    grid.y = (KM);

    //this sync is a bit over protective but currently here for debugging purposes
    cudaDeviceSynchronize();
    CUDA_CHECK_ERROR("Before buoydiff_gpu kernel execution");

    buoydiff_kernel1D<<<grid,threads,0,stream[1]>>>(DBLOC, DBSFC, d_TEMP, d_SALT, d_kmt, 0, KM);
    //buoydiff_kernel1D<<<grid,threads,0,stream[1]>>>(DBLOC, DBSFC, d_TRCR, d_TRCR+(NX_BLOCK*NY_BLOCK*KM), d_kmt, 0, KM);
    //debugging
    //buoydiff_kernel1D<<<grid,threads,0,stream[1]>>>(DBLOC, DBSFC, TRCR, TRCR+(NX_BLOCK*NY_BLOCK*KM), d_kmt, 0, KM);
    
    cudaDeviceSynchronize();
    CUDA_CHECK_ERROR("After buoydiff_gpu kernel execution");

    //this free() should be removed if we want to reuse TRCR info in
    //other kernels from vmix_kpp
    cudaFree(d_TEMP);
    cudaFree(d_SALT);
    
}

//device version of state for rho only used in buoydiff GPU kernel
__device__ double state(double temp, double salt, int k) {
  double tq, sq, sqr, work1, work2;//, denomk;

  tq = min(temp, 999.0);		//d_tmax[k]
  tq = max(tq, -2.0);			//d_tmin[k]

  sq = min(salt, 0.999);		//d_smax[k]
  sq = 1000.0 * max(sq, 0.0);	//d_smin[k]

  sqr = sqrt(sq);

  work1 = d_mwjfnums0t0[k] + tq * (d_mwjfnums0t1 + tq * (d_mwjfnums0t2[k] + d_mwjfnums0t3 * tq)) +
                        sq * (d_mwjfnums1t0[k] + d_mwjfnums1t1 * tq + d_mwjfnums2t0 * sq);

  work2 = d_mwjfdens0t0[k] + tq * (d_mwjfdens0t1[k] + tq * (d_mwjfdens0t2 +
     tq * (d_mwjfdens0t3[k] + d_mwjfdens0t4 * tq))) +
     sq * (d_mwjfdens1t0 + tq * (d_mwjfdens1t1 + tq*tq*d_mwjfdens1t3)+
     sqr * (d_mwjfdensqt0 + tq*tq*d_mwjfdensqt2));

  //denomk = 1.0/work2;
  //return work1*denomk;
  
  return work1/work2;
}

__global__ void buoydiff_kernel1D(double *DBLOC, double *DBSFC, double *TEMP, double *SALT, int *KMT, int start_k, int end_k) {

  int i = blockIdx.y * gridDim.x * blockDim.x + blockIdx.x * blockDim.x + threadIdx.x;
  int k = start_k + (i / (NX_BLOCK*NY_BLOCK));
  int sfci = i - (i / (NX_BLOCK*NY_BLOCK))*(NX_BLOCK*NY_BLOCK);
  //obtain array index
  int index = i + start_k*NX_BLOCK*NY_BLOCK;
  int indexmk = index-(NX_BLOCK*NY_BLOCK);

  double rho1, rhokm, rhok, dbloc;

  if (i < NX_BLOCK*NY_BLOCK*(end_k-start_k)) {

    //if k==0 we write DBSFC=0 and exit, as DBLOC for k=0 is computed by k=1
	if (k == 0) {
		DBSFC[index] = 0.0;
		return;
	} 
	if (k == KM-1) {
		DBLOC[index] = 0.0;
	}
	
	//double tempsfc = max(TEMP[sfci],-2.0);
	//double tempmk  = max(TEMP[indexmk],-2.0);
	//double tempk   = max(TEMP[index],-2.0);
	
	rho1  = state(TEMP[sfci], SALT[sfci], k);
	rhokm = state(TEMP[indexmk], SALT[indexmk], k);
	rhok  = state(TEMP[index], SALT[index], k);
	
	if (rhok != 0.0) { //prevent div by zero
		DBSFC[index]   = d_grav*(1.0 - rho1/rhok);
		//debug DBLOC[indexmk] = 1337.0;
		dbloc = d_grav*(1.0 - rhokm/rhok);
		//DBLOC[indexmk] = d_grav*(1.0 - rhokm/rhok);
	} else {
		DBSFC[index]   = 0.0;
		dbloc = 0.0;
		//DBLOC[indexmk] = 0.0;
	}
	
	if (k >= KMT[sfci]){ //-1 removed because FORTRAN array index starts at 1
		dbloc = 0.0;
		//DBLOC[indexmk] = 0.0;
	}
	DBLOC[indexmk] = dbloc;

	
  }
}
