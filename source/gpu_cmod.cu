#include <stdio.h>
#include <math.h>

#include "gpu_domain.h"

#define REUSE_TRCR 1
//#define USE_READ_ONLY_CACHE 1

#define TMIN -2.0
#define TMAX 999.0
#define SMIN 0.0
#define SMAX 0.999

#define RRHO0 2.55
#define DSFMAX 1.0

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

void cuda_init(int *pmy_task);
void devsync();

void my_cudamallochost(void **hostptr, int* size);

void cuda_state_initialize(double *constants, double *pressz,
        double *tmin, double *tmax, double *smin, double *smax, int *pnblocks, int *kmt);

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

void buoydiff_gpu(double *DBLOC, double *DBSFC, double *TRCR, int *pbid);
void ddmix_gpu(double *VDC, double *TRCR);

__global__ void buoydiff_kernel_onek(double *DBLOC, double *DBSFC, double *TEMP, double *SALT, int *KMT, int start_k, int end_k);
__global__ void ddmix_kernel_onek(double *VDC1, double *VDC2, double *TEMP, double *SALT, int start_k);

}




//functions

int cuda_initialized = 0;

//Fortran entry for initializing CUDA
void cuda_init(int *pmy_task) {
  if (cuda_initialized == 0) {
    cuda_initialized = 1;
    
    my_task = *pmy_task;
    
    int deviceCount = 0;
    cudaError_t err = cudaGetDeviceCount(&deviceCount);
    if (err != cudaSuccess) fprintf(stderr, "Error in cuda initialization: %s\n", cudaGetErrorString( err ));

    int dev = my_task%deviceCount;
    
    cudaSetDeviceFlags(cudaDeviceMapHost);
    cudaSetDevice(dev);

    
  }
}

void devsync() {
	cudaDeviceSynchronize();
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
cudaStream_t stream[KM];
cudaEvent_t event_htod[KM];
cudaEvent_t event_comp[KM];
cudaEvent_t event_dtoh[KM];


//debugging
int my_task;

void cuda_state_initialize(double *constants, double *pressz,
        double *tmin, double *tmax, double *smin, double *smax, int *pnblocks, int *kmt) {
  cudaError_t err;
	
  if (cuda_state_initialized == 0) {
  cuda_state_initialized = 1;
  
  //Perhaps we need to find some better way to pass constants around
  //for now we do this because we want the GPU/C code to use the exact same
  //values as the Fortran code even if both parts are compiled with
  //different compilers

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

  //error checking
  cudaDeviceSynchronize();
  CUDA_CHECK_ERROR("After cudaMemcpyToSymbols");
  
  int nblocks = *pnblocks;
  
  //printf("Node %d: nblocks=%d\n", my_task, nblocks);
  
  err = cudaMalloc(&d_kmt, NX_BLOCK*NY_BLOCK*nblocks*sizeof(int));
  if (err != cudaSuccess) fprintf(stderr, "Error doing cudaMalloc d_kmt\n");
  err = cudaMemcpy(d_kmt, kmt, NX_BLOCK*NY_BLOCK*nblocks*sizeof(int), cudaMemcpyHostToDevice);
  if (err != cudaSuccess) fprintf(stderr, "Error doing cudaMemcpyHostToDevice KMT\n");

  //error checking
  cudaDeviceSynchronize();
  CUDA_CHECK_ERROR("After cudaMemcpy KMT");

  //setup streams
  for (k=0; k<KM; k++) {
    cudaStreamCreate(&stream[k]);
  }
  //create cuda events
  for (k=0; k<KM; k++) {
    err = cudaEventCreate(&event_htod[k]);
    if (err != cudaSuccess) fprintf(stderr, "Error in cudaEventCreate htod: %s\n", cudaGetErrorString( err ));
    err = cudaEventCreate(&event_comp[k]);
    if (err != cudaSuccess) fprintf(stderr, "Error in cudaEventCreate comp: %s\n", cudaGetErrorString( err ));
    err = cudaEventCreate(&event_dtoh[k]);
    if (err != cudaSuccess) fprintf(stderr, "Error in cudaEventCreate dtoh: %s\n", cudaGetErrorString( err ));
  }
  
  //error checking
  cudaDeviceSynchronize();
  CUDA_CHECK_ERROR("After cudaStream and event creates");
  


  
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
  //debugging
  /*
  if (my_task == 0) {
    printf("GPU_CMOD using pressz:\n");
    for (k=0; k<KM; k++) {
	    printf("%d=%20.17e\n",k,pressz[k]);
    }
    printf("\n");
  }
  */

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
  //cudaDeviceSynchronize();
  //CUDA_CHECK_ERROR("Before mwjf_state_1D kernel execution");
  
  mwjf_state_1D<<<grid,threads,0,stream[1]>>>(TEMPK, SALTK, RHOOUT, DRHODT, DRHODS,
        n_outputs, start_k, end_k);
  
  //synchronize because we currently don't know when inputs or outputs will be used by CPU
  //the more this sync can be delayed the more overlap with CPU execution can be exploited
  
  //cudaDeviceSynchronize();
  //CUDA_CHECK_ERROR("After mwjf_state_1D kernel execution");
  
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
  
  //synchronize can be delayed to increase overlap with CPU execution
  //cudaDeviceSynchronize();
  //CUDA_CHECK_ERROR("After mwjf_state_1D kernel execution");
  
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
        tq = min(TEMPK[index], TMAX);	//d_tmax[k]
        tq = max(tq, TMIN);				//d_tmin[k]

        sq = min(SALTK[index], SMAX);	//d_smax[k]
        sq = 1000.0 * max(sq, SMIN);		//d_smin[k]

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
        tq = min(TEMPK[index], TMAX);	//d_tmax[k]
        tq = max(tq, TMIN);				//d_tmin[k]

        sq = min(SALTK[index], SMAX);	//d_smax[k]
        sq = 1000.0 * max(sq, SMIN);		//d_smin[k]

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

const char *var_names[7] = { "ERROR", "RHOOUT", "DBLOC", "DBSFC", "STEPMOD_RHO", "ADVT_PD", "VDC" };


void gpu_compare (double *a1, double *a2, int *pN, int *pName) {
  int N = *pN;
  int vName = *pName;
  int i=0, res=0;
  int print = 0;
  int zero_one = 0;
  int zero_two = 0;
  double eps = 0.00000000001;
  
  if (vName < 0 || vName > 6) { vName = 0; }

  for (i=0; i<N; i++) {

    if (a1[i] < eps && a1[i] > -eps) { zero_one++; }
    if (a2[i] < eps && a2[i] > -eps) { zero_two++; }

    if (isnan(a1[i]) || isnan(a2[i])) {
        res++;
        if (print < 10) {
          print++;
          fprintf(stderr, "Node %d: %s Error detected at i=%d, a1= %20.17e a2= %20.17e\n",my_task,var_names[vName],i,a1[i],a2[i]);
        }
    }

    double diff = a1[i]-a2[i];
    if (diff > eps || diff < -eps) {
        res++;
        if (print < 10) {
          print++;
          fprintf(stderr, "Node %d: %s Error detected at i=%d, \t a1= \t %20.17e \t a2= \t %20.17e\n",my_task,var_names[vName],i,a1[i],a2[i]);
        }
    }

  }

  //if (zero_one > 95*(N/100)) { fprintf(stderr, "Node %d: Error: array1 contains %d zeros\n",my_task, zero_one); }
  //if (zero_two > 95*(N/100)) { fprintf(stderr, "Node %d: Error: array2 contains %d zeros\n",my_task, zero_two); }

  //if (zero_one != zero_two) {
  //  fprintf(stderr, "Node %d: %s Error: number of zeros in arrays dont correspond zero1=%d, zero2=%d\n", my_task, var_names[vName], zero_one, zero_two);
  //}

  if (res > 0) {
    if (vName == 0) {
      fprintf(stdout,"Node %d: Number of errors in GPU result: %d\n",my_task,res);
    } else {
	  fprintf(stdout,"Node %d: Number of errors in %s GPU result: %d\n",my_task,var_names[vName],res);
    }
  }
}

double *d_TRCR  = (double *) 0;
double *d_DBLOC = (double *) 0;
double *d_DBSFC = (double *) 0;

void buoydiff_gpu(double *DBLOC, double *DBSFC, double *TRCR, int *pbid) {
	  cudaError_t err;

	  //allocate device memory
	  double *d_SALT;
	  double *SALT = TRCR+NX_BLOCK*NY_BLOCK*KM;
	  
	  int bid = (*pbid) - 1; //-1 because fortran indices start at 0
	  //printf("Node %d: bid=%d\n", my_task, bid);
      if (d_DBLOC == (double *) 0) {
	    err = cudaMalloc((void **)&d_DBLOC, NX_BLOCK*NY_BLOCK*KM*sizeof(double));
	    if (err != cudaSuccess) fprintf(stderr, "Error in cudaMalloc d_DBLOC: %s\n", cudaGetErrorString( err ));
      }
      if (d_DBSFC == (double *) 0) {
	    err = cudaMalloc((void **)&d_DBSFC, NX_BLOCK*NY_BLOCK*KM*sizeof(double));
	    if (err != cudaSuccess) fprintf(stderr, "Error in cudaMalloc d_DBSFC: %s\n", cudaGetErrorString( err ));
      }
	  
      if (d_TRCR == (double *) 0) {
	    err = cudaMalloc((void **)&d_TRCR, NX_BLOCK*NY_BLOCK*KM*2*sizeof(double));
	    if (err != cudaSuccess) fprintf(stderr, "Error in cudaMalloc d_TRCR %s\n", cudaGetErrorString( err ));
      }
	  d_SALT =d_TRCR+NX_BLOCK*NY_BLOCK*KM;
	  
	  //only used in debugging
	  //cudaDeviceSynchronize();
	  //CUDA_CHECK_ERROR("After buoydiff memory setup");

	  //setup execution parameters
	  dim3 threads(32,8);
	  dim3 grid(1,1);
	  grid.y = (int)ceilf((float)NY_BLOCK  / (float)(threads.y));
	  grid.x = (int)ceilf((float)NX_BLOCK  / (float)(threads.x));

	  int array_size = NX_BLOCK*NY_BLOCK;
	  int lps = 1; //levels to compute per stream
	  int k = 0;

	  for (k=0; k<KM; k+=lps) {
	    err = cudaMemcpyAsync(d_TRCR+k*array_size, TRCR+k*array_size, lps*array_size*sizeof(double), cudaMemcpyHostToDevice, stream[k]);
	    if (err != cudaSuccess) fprintf(stderr, "Error in cudaMemcpy host to device TRCR: %s\n", cudaGetErrorString( err ));
	    err = cudaMemcpyAsync(d_SALT+k*array_size, SALT+k*array_size, lps*array_size*sizeof(double), cudaMemcpyHostToDevice, stream[k]);
	    if (err != cudaSuccess) fprintf(stderr, "Error in cudaMemcpy host to device TRCR: %s\n", cudaGetErrorString( err ));

	    //record cuda event
	    err = cudaEventRecord (event_htod[k], stream[k]);
	    if (err != cudaSuccess) fprintf(stderr, "Error in cudaEventRecord htod: %s\n", cudaGetErrorString( err ));
	  }

	  for (k=0; k<KM; k+=lps) {
	    if (k > 0) {
	      //wait for memcpy in stream k=0 to be complete
	      err = cudaStreamWaitEvent(stream[k], event_htod[0], 0);
	      if (err != cudaSuccess) fprintf(stderr, "Error in cudaStreamWaitEvent htod k=0: %s\n", cudaGetErrorString( err ));
	      //wait for memcpy in stream k-1 to be complete
	      err = cudaStreamWaitEvent(stream[k], event_htod[k-lps], 0);
	      if (err != cudaSuccess) fprintf(stderr, "Error in cudaStreamWaitEvent htod k-1: %s\n", cudaGetErrorString( err ));
	    }

	    buoydiff_kernel_onek<<<grid,threads,0,stream[k]>>>(d_DBLOC, d_DBSFC, d_TRCR, d_SALT, d_kmt+(bid*NX_BLOCK*NY_BLOCK), k, k+lps);
	    
	    err = cudaEventRecord (event_comp[k], stream[k]);
	    if (err != cudaSuccess) fprintf(stderr, "Error in cudaEventRecord htod: %s\n", cudaGetErrorString( err ));
	  }

	  for (k=0; k<KM; k+=lps) {
	    err = cudaMemcpyAsync(DBSFC+k*array_size, d_DBSFC+k*array_size, lps*array_size*sizeof(double), cudaMemcpyDeviceToHost, stream[k]);
	    if (err != cudaSuccess) fprintf(stderr, "Error in cudaMemcpy device to host DBSFC: %s\n", cudaGetErrorString( err ));

	    if (k < KM-1) {
	      //wait for computation in stream k+1 to be complete
	      err = cudaStreamWaitEvent(stream[k], event_comp[k+lps], 0);
	      if (err != cudaSuccess) fprintf(stderr, "Error in cudaStreamWaitEvent comp k+1: %s\n", cudaGetErrorString( err ));
	    }
	    err = cudaMemcpyAsync(DBLOC+k*array_size, d_DBLOC+k*array_size, lps*array_size*sizeof(double), cudaMemcpyDeviceToHost, stream[k]);
	    if (err != cudaSuccess) fprintf(stderr, "Error in cudaMemcpy device to host DBLOC: %s\n", cudaGetErrorString( err ));
	  }
	  
      //wait for device to finish
	  cudaDeviceSynchronize();
	  CUDA_CHECK_ERROR("After buoydiff_gpu kernel execution");
	  
      //we might move the malloc-free pairs to initialization and finalization routines in the end.
	  //cudaFree(d_DBLOC);
	  //cudaFree(d_DBSFC);
	  
	  //TRCR values may remain on the GPU for other vmix_kpp routines
	#ifndef REUSE_TRCR
	  //not reusing trcr so free it
	  //err = cudaFree(d_TRCR);
	  //if (err != cudaSuccess) fprintf(stderr, "Error in buoydiff cudaFree() d_TRCR: %s\n", cudaGetErrorString( err ));
	#endif

}

//device version of state for rho only used in buoydiff GPU kernel
__device__ double state(double temp, double salt, int k) {
  double tq, sq, sqr, work1, work2;//, denomk;

  tq = min(temp, TMAX);		//d_tmax[k]
  tq = max(tq, TMIN);			//d_tmin[k]

  sq = min(salt, SMAX);		//d_smax[k]
  sq = 1000.0 * max(sq, SMIN);	//d_smin[k]

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

__global__ void buoydiff_kernel_onek(double *DBLOC, double *DBSFC, double *TEMP, double *SALT, int *KMT, int start_k, int end_k) {
  //obtain indices
  int j = threadIdx.y + blockIdx.y * blockDim.y;
  int i = threadIdx.x + blockIdx.x * blockDim.x;
  int k = start_k;

  if (j < NY_BLOCK && i < NX_BLOCK) {

    if (k == 0) {
       DBSFC[i + j*NX_BLOCK] = 0.0;
    } else {
        double rho1, rhokm, rhok;
        double tempsfc = TEMP[i + j*NX_BLOCK];
        double saltsfc = SALT[i + j*NX_BLOCK];
        int kmt = KMT[i + j*NX_BLOCK];
        
        if (k == KM-1) {
           DBLOC[i + j*NX_BLOCK+k*NY_BLOCK*NX_BLOCK] = 0.0;
        }
        
        rho1 = state(tempsfc, saltsfc, k);
        rhokm = state(TEMP[i + j*NX_BLOCK+(k-1)*NY_BLOCK*NX_BLOCK], SALT[i + j*NX_BLOCK+(k-1)*NY_BLOCK*NX_BLOCK], k);
        rhok =  state(TEMP[i + j*NX_BLOCK+k*NY_BLOCK*NX_BLOCK], SALT[i + j*NX_BLOCK+k*NY_BLOCK*NX_BLOCK], k);

        if (rhok != 0.0) { //prevent div by zero
           DBSFC[i + j*NX_BLOCK+k*NY_BLOCK*NX_BLOCK]     = d_grav*(1.0 - rho1/rhok);
           DBLOC[i + j*NX_BLOCK+(k-1)*NY_BLOCK*NX_BLOCK] = d_grav*(1.0 - rhokm/rhok);
        } else {
           DBSFC[i + j*NX_BLOCK+k*NY_BLOCK*NX_BLOCK]     = 0.0;
           DBLOC[i + j*NX_BLOCK+(k-1)*NY_BLOCK*NX_BLOCK] = 0.0;
        }

        //zero if on land
        //why DBSFC isnt zeroed here is a mystery to me
        if (k >= kmt){ //-1 removed because FORTRAN array index starts at 1
           DBLOC[i + j*NX_BLOCK+(k-1)*NY_BLOCK*NX_BLOCK] = 0.0;
        }

      } // end of if k==0, else block
  }//end of bounds check
}

double *d_VDC = (double *) 0;

void ddmix_gpu(double *VDC, double *TRCR) {
	  cudaError_t err;
	
	  //allocate device memory
	  double *d_VDC1;
	  double *d_VDC2;
	  
	  if (d_VDC == (double *) 0) {
	    err = cudaMalloc((void **)&d_VDC, NX_BLOCK*NY_BLOCK*(KM+2)*2*sizeof(double));
	    if (err != cudaSuccess) fprintf(stderr, "Error in ddmix cudaMalloc d_VDC: %s\n", cudaGetErrorString( err ));
	  }
	  d_VDC1 = d_VDC+(NX_BLOCK*NY_BLOCK); //skip first level in VDC1 and VDC2
	  d_VDC2 = d_VDC1+(NX_BLOCK*NY_BLOCK*(KM+2));
	  
	#ifndef REUSE_TRCR
	  if (d_TRCR == (double *) 0) {
	    err = cudaMalloc((void **)&d_TRCR, NX_BLOCK*NY_BLOCK*KM*2*sizeof(double));
	    if (err != cudaSuccess) fprintf(stderr, "Error in ddmix cudaMalloc d_TRCR %s\n", cudaGetErrorString( err ));
	  }
	#endif
	  
	  //only for debugging
	  //cudaDeviceSynchronize();
	  //CUDA_CHECK_ERROR("After ddmix memory setup");

	  //setup execution parameters
	  dim3 threads(16,16);
	  dim3 grid(1,1);
	  grid.y = (int)ceilf((float)NY_BLOCK  / (float)(threads.y));
	  grid.x = (int)ceilf((float)NX_BLOCK  / (float)(threads.x));

	  //stream specific stuff
	  int array_size = NX_BLOCK*NY_BLOCK;
	  int lps = 1; //levels to compute per stream
	  int k = 0;

	  //separate pointers for VDC1 and VDC2, skip first level because of k=0 in FORTRAN
	  double *VDC1 = VDC+(NX_BLOCK*NY_BLOCK);
	  double *VDC2 = VDC+(NX_BLOCK*NY_BLOCK)+(NX_BLOCK*NY_BLOCK*(KM+2));

	  for (k=0; k<KM; k+=lps) {
	    //cpy tracers
	#ifndef REUSE_TRCR
	    err = cudaMemcpyAsync(d_TRCR+k*array_size, TRCR+k*array_size, lps*array_size*sizeof(double), cudaMemcpyHostToDevice, stream[k]);
	    if (err != cudaSuccess) fprintf(stderr, "Error in ddmix cudaMemcpy host to device TRCR: %s\n", cudaGetErrorString( err ));
	    err = cudaMemcpyAsync(d_TRCR+(NX_BLOCK*NY_BLOCK*KM)+k*array_size, TRCR+(NX_BLOCK*NY_BLOCK*KM)+k*array_size, lps*array_size*sizeof(double), cudaMemcpyHostToDevice, stream[k]);
	    if (err != cudaSuccess) fprintf(stderr, "Error in ddmix cudaMemcpy host to device TRCR: %s\n", cudaGetErrorString( err ));

	    //record cuda event
	    err = cudaEventRecord (event_htod[k], stream[k]);
	    if (err != cudaSuccess) fprintf(stderr, "Error in ddmix cudaEventRecord: %s\n", cudaGetErrorString( err ));
	#endif

	    if (k<KM-1) {
	      //cpy vdc
	      err = cudaMemcpyAsync(d_VDC1+k*array_size, VDC1+k*array_size, lps*array_size*sizeof(double), cudaMemcpyHostToDevice, stream[k]);
	      if (err != cudaSuccess) fprintf(stderr, "Error in cudaMemcpy host to device VDC1: %s\n", cudaGetErrorString( err ));
	      err = cudaMemcpyAsync(d_VDC2+k*array_size, VDC2+k*array_size, lps*array_size*sizeof(double), cudaMemcpyHostToDevice, stream[k]);
	      if (err != cudaSuccess) fprintf(stderr, "Error in cudaMemcpy host to device VDC2: %s\n", cudaGetErrorString( err ));
	    }

	  }

	  //call the kernel
	  for (k=0; k<KM-1; k+=lps) {
	    //wait for memcpys in stream K+1 to complete to guarantee correctness
	#ifndef REUSE_TRCR
	    if (k+lps < KM) {
	      err = cudaStreamWaitEvent(stream[k], event_htod[k+lps], 0);
	      if (err != cudaSuccess) fprintf(stderr, "Error in cudaStreamWaitEvent: %s\n", cudaGetErrorString( err ));
	    }
	#endif

	    ddmix_kernel_onek<<<grid,threads,0,stream[k]>>>(d_VDC, d_VDC+(NX_BLOCK*NY_BLOCK*(KM+2)), d_TRCR, d_TRCR+(NX_BLOCK*NY_BLOCK*KM), k);
	  }

	  //device to host copies
	  for (k=0; k<KM-1; k+=lps) {
	    err = cudaMemcpyAsync(VDC1+k*array_size, d_VDC1+k*array_size, lps*array_size*sizeof(double), cudaMemcpyDeviceToHost, stream[k]);
	    if (err != cudaSuccess) fprintf(stderr, "Error in cudaMemcpy host to device VDC1: %s\n", cudaGetErrorString( err ));
	    err = cudaMemcpyAsync(VDC2+k*array_size, d_VDC2+k*array_size, lps*array_size*sizeof(double), cudaMemcpyDeviceToHost, stream[k]);
	    if (err != cudaSuccess) fprintf(stderr, "Error in cudaMemcpy host to device VDC2: %s\n", cudaGetErrorString( err ));
	  }

	  //wait for completion
	  //this sync is delayed to stimulate overlap with CPU computation of bldepth()
	  //cudaDeviceSynchronize();
	  //CUDA_CHECK_ERROR("After ddmix_gpu kernel execution");
	  
	  
	  
	  //cudaFree(d_VDC);
	  //whether or not trcr was reused, we could free it now
	  //err = cudaFree(d_TRCR);
	  //if (err != cudaSuccess) fprintf(stderr, "Error in ddmix cudaFree() d_TRCR: %s\n", cudaGetErrorString( err ));

}

__global__ void ddmix_kernel_onek(double *VDC1, double *VDC2, double *TEMP, double *SALT, int start_k) {
	double talpha, sbeta;
	double salt, temp, temp_kup, salt_kup, talpha_kup, sbeta_kup;
	double alphadt, betads, diffdd;
	double vdc1, vdc2;

	//obtain indices
	int j = threadIdx.y + blockIdx.y * blockDim.y;
	int i = threadIdx.x + blockIdx.x * blockDim.x;
	int k = start_k;

	//check bounds
	if (j < NY_BLOCK && i < NX_BLOCK) {
		//kup  = 1
		//knxt = 2

#ifdef USE_READ_ONLY_CACHE
		temp_kup = __ldg( TEMP+(i+(j*NX_BLOCK)+(k*NX_BLOCK*NY_BLOCK)) );
		salt_kup = __ldg( SALT+(i+(j*NX_BLOCK)+(k*NX_BLOCK*NY_BLOCK)) );
		temp     = __ldg( TEMP+(i+(j*NX_BLOCK)+((k+1)*NX_BLOCK*NY_BLOCK)) );
		salt     = __ldg( SALT+(i+(j*NX_BLOCK)+((k+1)*NX_BLOCK*NY_BLOCK)) );
#else
		temp_kup = TEMP[i+(j*NX_BLOCK)+(k*NX_BLOCK*NY_BLOCK)];
		salt_kup = SALT[i+(j*NX_BLOCK)+(k*NX_BLOCK*NY_BLOCK)];
		temp     = TEMP[i+(j*NX_BLOCK)+((k+1)*NX_BLOCK*NY_BLOCK)];
		salt     = SALT[i+(j*NX_BLOCK)+((k+1)*NX_BLOCK*NY_BLOCK)];
#endif

		vdc1 = VDC1[i + j*NX_BLOCK+(k+1)*NX_BLOCK*NY_BLOCK];
		vdc2 = VDC2[i + j*NX_BLOCK+(k+1)*NX_BLOCK*NY_BLOCK];

		//computed rrho here is actually not used and overwritten by next call to state
		//   call state(1, 1, prandtl, salt_kup, &
		//                    RHOFULL=rrho, &
		//                    DRHODT=talpha_kup, DRHODS=sbeta_kup)
		//inlined function state
		double tq, sq, sqr, work1, work2, work3, work4, denomk;
		tq = min(temp_kup, TMAX);
		tq = max(tq, TMIN);
		sq = min(salt_kup, SMAX);
		sq = 1000.0 * max(sq, SMIN);

		sqr = sqrt(sq);

		work1 = d_mwjfnums0t0[k] + tq * (d_mwjfnums0t1 + tq * (d_mwjfnums0t2[k] + d_mwjfnums0t3 * tq)) +
				sq * (d_mwjfnums1t0[k] + d_mwjfnums1t1 * tq + d_mwjfnums2t0 * sq);

		work2 = d_mwjfdens0t0[k] + tq * (d_mwjfdens0t1[k] + tq * (d_mwjfdens0t2 +
				tq * (d_mwjfdens0t3[k] + d_mwjfdens0t4 * tq))) +
				sq * (d_mwjfdens1t0 + tq * (d_mwjfdens1t1 + tq*tq*d_mwjfdens1t3)+
						sqr * (d_mwjfdensqt0 + tq*tq*d_mwjfdensqt2));

		denomk = 1.0/work2;
		//unused        rrho = work1*denomk;

		work3 = // dP_1/dT
				d_mwjfnums0t1 + tq * (2.0*d_mwjfnums0t2[k] +
						3.0*d_mwjfnums0t3 * tq) + d_mwjfnums1t1 * sq;

		work4 = // dP_2/dT
				d_mwjfdens0t1[k] + sq * d_mwjfdens1t1 +
				tq * (2.0*(d_mwjfdens0t2 + sq*sqr*d_mwjfdensqt2) +
						tq * (3.0*(d_mwjfdens0t3[k] + sq * d_mwjfdens1t3) +
								tq *  4.0*d_mwjfdens0t4));

		talpha_kup = (work3 - work1*denomk*work4)*denomk;

		work3 = // dP_1/dS
				d_mwjfnums1t0[k] + d_mwjfnums1t1 * tq + 2.0*d_mwjfnums2t0 * sq;

		work4 = d_mwjfdens1t0 +   // dP_2/dS
				tq * (d_mwjfdens1t1 + tq*tq*d_mwjfdens1t3) +
				1.5*sqr*(d_mwjfdensqt0 + tq*tq*d_mwjfdensqt2);

		sbeta_kup = (work3 - work1*denomk*work4)*denomk * 1000.0;
		//end of inlined function state

		//real work starts here
		//   for (k = start_k; k < end_k; k++) {
		{ 
			//   do k=1,KM
			double prandtl = 0.0, rrho = 0.0;

			if ( k < KM-1 ) {//changed to KM-1 because we start at 0

				prandtl = max(temp, TMIN);
				//PRANDTL = merge(-c2,TRCR(:,:,k+1,1),TRCR(:,:,k+1,1) < -c2)


				//         call state(k+1, k+1, prandtl, salt,              &
				//                              RHOFULL=rrho, DRHODT=talpha, &
				//                                            DRHODS=sbeta)
				//inlined function state (k+1)
				double tq, sq, sqr, work1, work2, work3, work4, denomk;
				tq = min(temp, TMAX);
				tq = max(tq, TMIN);
				sq = min(salt, SMAX);
				sq = 1000.0 * max(sq, SMIN);
				sqr = sqrt(sq);

				work1 = d_mwjfnums0t0[k+1] + tq * (d_mwjfnums0t1 + tq * (d_mwjfnums0t2[k+1] + d_mwjfnums0t3 * tq)) +
						sq * (d_mwjfnums1t0[k+1] + d_mwjfnums1t1 * tq + d_mwjfnums2t0 * sq);
				work2 = d_mwjfdens0t0[k+1] + tq * (d_mwjfdens0t1[k+1] + tq * (d_mwjfdens0t2 +
						tq * (d_mwjfdens0t3[k+1] + d_mwjfdens0t4 * tq))) +
						sq * (d_mwjfdens1t0 + tq * (d_mwjfdens1t1 + tq*tq*d_mwjfdens1t3)+
								sqr * (d_mwjfdensqt0 + tq*tq*d_mwjfdensqt2));
				denomk = 1.0/work2;
				rrho = work1*denomk;

				work3 = // dP_1/dT
						d_mwjfnums0t1 + tq * (2.0*d_mwjfnums0t2[k+1] +
								3.0*d_mwjfnums0t3 * tq) + d_mwjfnums1t1 * sq;
				work4 = // dP_2/dT
						d_mwjfdens0t1[k+1] + sq * d_mwjfdens1t1 +
						tq * (2.0*(d_mwjfdens0t2 + sq*sqr*d_mwjfdensqt2) +
								tq * (3.0*(d_mwjfdens0t3[k+1] + sq * d_mwjfdens1t3) +
										tq *  4.0*d_mwjfdens0t4));
				talpha = (work3 - work1*denomk*work4)*denomk;

				work3 = // dP_1/dS
						d_mwjfnums1t0[k+1] + d_mwjfnums1t1 * tq + 2.0*d_mwjfnums2t0 * sq;
				work4 = d_mwjfdens1t0 +   // dP_2/dS
						tq * (d_mwjfdens1t1 + tq*tq*d_mwjfdens1t3) +
						1.5*sqr*(d_mwjfdensqt0 + tq*tq*d_mwjfdensqt2);
				sbeta = (work3 - work1*denomk*work4)*denomk * 1000.0;
				//end of inlined function state

				alphadt = -0.5*(talpha_kup + talpha) * (temp_kup - temp);
				//         ALPHADT = -p5*(TALPHA(:,:,kup) + TALPHA(:,:,knxt)) &
				//                      *(TRCR(:,:,k,1) - TRCR(:,:,k+1,1))

				betads  = 0.5*(sbeta_kup + sbeta) * (salt_kup - salt);
				//         BETADS  = p5*( SBETA(:,:,kup) +  SBETA(:,:,knxt)) &
				//                     *(TRCR(:,:,k,2) - TRCR(:,:,k+1,2))

				//         kup  = knxt
				//         knxt = 3 - kup

			} else {
				alphadt = 0.0;
				betads  = 0.0;
			}

			//!-----------------------------------------------------------------------
			//!
			//!     salt fingering case
			//!
			//!-----------------------------------------------------------------------


			//where ( ALPHADT > BETADS .and. BETADS > c0 )
			if ((alphadt > betads) && (betads > 0.0)) {

				rrho       = min(alphadt/betads, RRHO0);
				//RRHO       = MIN(ALPHADT/BETADS, Rrho0)
				diffdd     = (1.0-((rrho-1.0)/(RRHO0-1.0)));
				diffdd     = diffdd * diffdd * diffdd;
				diffdd     = DSFMAX * diffdd;
				//DIFFDD     = dsfmax*(c1-(RRHO-c1)/(Rrho0-c1))**3

				vdc1 += 0.7*diffdd; //k+1 because k=0 is only zeros
				vdc2 += diffdd;
				//VDC(:,:,k,1) = VDC(:,:,k,1) + 0.7*DIFFDD
				//VDC(:,:,k,2) = VDC(:,:,k,2) + DIFFDD

			} //endwhere
			//!-----------------------------------------------------------------------
			//!
			//!     diffusive convection
			//!
			//!-----------------------------------------------------------------------

			//where ( ALPHADT < c0 .and. BETADS < c0 .and. ALPHADT > BETADS )
			if ((alphadt < 0.0) && (betads < 0.0) && (alphadt > betads)) {
				rrho    = alphadt/ betads;
				//RRHO    = ALPHADT / BETADS
				diffdd  = 1.5e-2 * 0.909 * exp(4.6*exp(-0.54*(1.0/rrho-1.0)));
				//DIFFDD  = 1.5e-2_c_double*0.909_c_double* &
				//          exp(4.6_c_double*exp(-0.54_c_double*(c1/RRHO-c1)))
				prandtl = 0.15 *rrho;
				//PRANDTL = 0.15_c_double*RRHO
			} else { //elsewhere
				rrho    = 0.0;
				diffdd  = 0.0;
				prandtl = 0.0;
			} //endwhere


			//where (RRHO > p5) PRANDTL = (1.85_c_double - 0.85_c_double/RRHO)*RRHO
			if (rrho > 0.5) {
				//prandtl = (1.85 - (0.85/rrho))*rrho;
				//simplyfied
				prandtl = 1.85*rrho - 0.85;
			}

			VDC1[i + j*NX_BLOCK+(k+1)*NX_BLOCK*NY_BLOCK] = vdc1 + diffdd;
			VDC2[i + j*NX_BLOCK+(k+1)*NX_BLOCK*NY_BLOCK] = vdc2 + prandtl*diffdd;

		} //end do-loop over k

	} //end bounds check

}


