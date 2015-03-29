
#include <stdio.h>
#include <time.h>
#include <stdlib.h>

using namespace std;

#define DBLOC(i,j,k) DBLOC[(i)+(j)*nx_block+(k)*nx_block*ny_block]
#define VISC(i,j,k) VISC[(i)+(j)*nx_block+(k)*nx_block*ny_block]
#define VDC(i,j,k,l) VDC[(i)+(j)*nx_block+(k)*nx_block*ny_block+(l)*nx_block*ny_block*(km+1)]
#define VVC(i,j,k) VVC[(i)+(j)*nx_block+(k)*nx_block*ny_block]
#define KMU(i,j,k) KMU[(i)+(j)*nx_block+(k)*nx_block*ny_block]
#define KPP_SRC(i,j,k,l,m) KPP_SRC[(i)+(j)*nx_block+(k)*nx_block*ny_block+(l)*nx_block*ny_block*km+(m)*nx_block*ny_block*km*nt]
#define STF(i,j,k) STF[(i)+(j)*nx_block+(k)*nx_block*ny_block]
#define dz(i) dz[(i)]
#define GHAT(i,j,k) GHAT[(i)+(j)*nx_block+(k)*nx_block*ny_block]
#define KMT(i,j,k) KMT[(i)+(j)*nx_block+(k)*nx_block*ny_block]
#define HMXL(i,j,k) HMXL[(i)+(j)*nx_block+(k)*nx_block*ny_block]
#define zt(i) zt[(i)]
#define DBSFC(i,j,k) DBSFC[(i)+(j)*nx_block+(k)*nx_block*ny_block]
#define UUU(i,j,k) UUU[(i)+(j)*nx_block+(k)*nx_block*ny_block]
#define VVV(i,j,k) VVV[(i)+(j)*nx_block+(k)*nx_block*ny_block]
#define SMF(i,j,k) SMF[(i)+(j)*nx_block+(k)*nx_block*ny_block]
#define TRCR(i,j,k,l) TRCR[(i)+(j)*nx_block+(k)*nx_block*ny_block+(l)*nx_block*ny_block*km]
#define SHF_QSW(i,j) SHF_QSW[(i)+(j)*nx_block]
#define RI_BULK(i) RI_BULK[(i)]
#define HBLT(i,j) HBLT[(i)+(j)*nx_block]
#define BFSFC(i,j) BFSFC[(i)+(j)*nx_block]
#define Ricr(i) Ricr[(i)]
#define KBL(i,j) KBL[(i)+(j)*nx_block]
#define GAT1(i) GAT1[(i)]
#define DAT1(i) DAT1[(i)]
#define STABLE(i,j) STABLE[(i)+(j)*nx_block]
#define DKM1(i) DKM1[(i)]
#define WM(i,j) WM[(i)+(j)*nx_block]
#define USTAR(i,j) USTAR[(i)+(j)*nx_block]
#define WS(i,j) WS[(i)+(j)*nx_block]
#define TALPHA(i) TALPHA[(i)]
#define SBETA(i) SBETA[(i)]
#define TEMPK(i) TEMPK[(i)]


//manually replaced
#define zgrid(i) zgrid[(i+1)] //offset is advanced to allow indexing from -1 
#define hwide(i) hwide[(i+1)] //offset is advanced to allow indexing from -1

#define DZT(i,j,k,l) DZT[(i)+(j)*nx_block+(k+1)*nx_block*ny_block+(l)*nx_block*ny_block*(km+2)]
#define DZU(i,j,k,l) DZU[(i)+(j)*nx_block+(k+1)*nx_block*ny_block+(l)*nx_block*ny_block*(km+2)]

//#define bckgrnd_vdc(i,j,k,l) bckgrnd_vdc[(i)+(j)*nx_block+(k)*nx_block*ny_block+(l)*nx_block*ny_block*km]
//#define bckgrnd_vvc(i,j,k,l) bckgrnd_vvc[(i)+(j)*nx_block+(k)*nx_block*ny_block+(l)*nx_block*ny_block*km]
#define bckgrnd_vdc(i,j,k,l) bckgrnd_vdc[(k)]
#define bckgrnd_vvc(i,j,k,l) bckgrnd_vvc[(k)]


//manual inserts

#define ARRAY_UGRID(i,j,k) ARRAY_UGRID[(i)+(j)*nx_block+(k)*nx_block*ny_block]
#define ARRAY_TGRID(i,j,k) ARRAY_TGRID[(i)+(j)*nx_block+(k)*nx_block*ny_block]
#define AU(i,j,b) AU[(i)+(j)*nx_block+(b)*nx_block*ny_block]
#define UAREA_R(i,j,b) UAREA_R[(i)+(j)*nx_block+(b)*nx_block*ny_block]

#define VSHEARU(i,j,k) VSHEARU[(i)+(j)*nx_block+(k)*nx_block*ny_block]

#define DEBUG(i,j,k) DEBUG[(i)+(j)*nx_block+(k)*nx_block*ny_block]
#define WORK3(i,j) WORK3[(i)+(j)*nx_block]

#define VISCT(i,j,k) VISCT[(i)+(j)*nx_block+(k)*nx_block*ny_block]


#define h_KMU(i,j,k) h_KMU[(i)+(j)*nx_block+(k)*nx_block*ny_block]
#define h_KMT(i,j,k) h_KMT[(i)+(j)*nx_block+(k)*nx_block*ny_block]


#define km 42
#define nx_block 64
#define ny_block 64
#define nt 2
#define max_blocks_clinic 7
#define nblocks_clinic 7

#define BLOCK_X 32
#define BLOCK_Y 16



#define sign(x, y) (y >= 0 ? abs(x) : -abs(x) )
#define nint(x) ((int) nearbyint(x))




/* constants for equation of state */

double      mwjfnp0s0t0 =   9.9984369900000e+2 * 0.001;
double      mwjfnp0s0t1 =   7.3521284000000e+0 * 0.001;
double      mwjfnp0s0t2 =  -5.4592821100000e-2 * 0.001;
double      mwjfnp0s0t3 =   3.9847670400000e-4 * 0.001;
double      mwjfnp0s1t0 =   2.9693823900000e+0 * 0.001;
double      mwjfnp0s1t1 =  -7.2326881300000e-3 * 0.001;
double      mwjfnp0s2t0 =   2.1238234100000e-3 * 0.001;
double      mwjfnp1s0t0 =   1.0400459100000e-2 * 0.001;
double      mwjfnp1s0t2 =   1.0397052900000e-7 * 0.001;
double      mwjfnp1s1t0 =   5.1876188000000e-6 * 0.001;
double      mwjfnp2s0t0 =  -3.2404182500000e-8 * 0.001;
double      mwjfnp2s0t2 =  -1.2386936000000e-11* 0.001;

double      mwjfdp0s0t0 =   1.0e+0;
double      mwjfdp0s0t1 =   7.2860673900000e-3;
double      mwjfdp0s0t2 =  -4.6083554200000e-5;
double      mwjfdp0s0t3 =   3.6839057300000e-7;
double      mwjfdp0s0t4 =   1.8080918600000e-10;
double      mwjfdp0s1t0 =   2.1469170800000e-3;
double      mwjfdp0s1t1 =  -9.2706248400000e-6;
double      mwjfdp0s1t3 =  -1.7834364300000e-10;
double      mwjfdp0sqt0 =   4.7653412200000e-6;
double      mwjfdp0sqt2 =   1.6341073600000e-9;
double      mwjfdp1s0t0 =   5.3084887500000e-6;
double      mwjfdp2s0t3 =  -3.0317512800000e-16;
double      mwjfdp3s0t1 =  -1.2793413700000e-17;

__constant__ double mwjfnums0t1;
__constant__ double mwjfnums0t3;
__constant__ double mwjfnums1t1;
__constant__ double mwjfnums2t0;
__constant__ double mwjfdens0t2;
__constant__ double mwjfdens0t4;
__constant__ double mwjfdens1t0;
__constant__ double mwjfdens1t1;
__constant__ double mwjfdens1t3;
__constant__ double mwjfdensqt0;
__constant__ double mwjfdensqt2;

__constant__ double mwjfnums0t0[km];
__constant__ double mwjfnums0t2[km];
__constant__ double mwjfnums1t0[km];

__constant__ double mwjfdens0t0[km];
__constant__ double mwjfdens0t1[km];
__constant__ double mwjfdens0t3[km];

#define TMIN -2.0
#define TMAX 999.0
#define SMIN 0.0
#define SMAX 0.999


/* end of constants for equation of state */

//constants for blmix
    __constant__ double cg;
    __constant__ double Vtc;



cudaStream_t stream[1] = { 0 }; 

extern "C" {

    void start_timer();
    void stop_timer(float *time);
    void fill_random(double *p, int *n);

    int compare_int(int *a1, int *a2, int N);
    int compare (double *a1, double *a2, int N, const char *);
    int compare_to (double *a1, double *a2, int N, double *a3);

    __global__ void interior_convection_kernel (double *DBLOC, double *DBSFC, int *KBL, double *STF, double *GHAT, double *VISC, double convect_diff, double convect_visc, int bid, double *VDC, double *VVC, double *KPP_SRC, double *HMXL, double *VISCT, double *DZT, int *KMT);
    __global__ void interior_convection_part2 (double *DBLOC, double *DBSFC, int *KBL, double *STF, double *GHAT, double *VISC, double convect_diff, double convect_visc, int bid, double *VDC, double *VVC, double *KPP_SRC, double *HMXL, double *AU, double *UAREA_R, double *VISCT, double *DZT, int *KMT, int *KMU);
    __global__ void compute_vshear_ugrid(double *VSHEARU, double *UUU, double *VVV, double *DZU, int bid);
    __global__ void ri_iwmix_kernel (double *DBLOC, double *VISC, double *VDC, double *UUU, double *VVV, double *RHOMIX, double convect_diff, double convect_visc, int bid, double *VSHEARU, double *DZT, int *KMT);
    __global__ void compute_ustaru(double *USTAR, double *SMF);
    __global__ void ugrid_to_tgrid_kernel (double *ARRAY_TGRID, double *ARRAY_UGRID, int k);
    __global__ void compute_vshearu_bldepth(double *VSHEARU, double *UUU, double *VVV, double *DZU, int bid);
    __global__ void bldepth_kernel (double *DBLOC, double *DBSFC, double *TRCR, double *UUU, double *VVV, double *STF, double *SHF_QSW, double *HBLT, double *USTAR, double *BFSFC, double *STABLE, int *KBL, int bid, double *SMF, double *VSHEARU, double *d_WORK1, double *d_WORK2, double *d_WORK3, double *DZT, double *DZU, int *KMT);
    __global__ void correct_stability_and_buoyancy(double *BFSFC, double *STABLE, double *BO, double *BOSOL, double *HBLT);
    __global__ void blmix_kernel (double *VISC, double *VDC, double *HBLT, double *USTAR, double *BFSFC, double *STABLE, int *KBL, double *GHAT, int bid, double *DZT);
    __global__ void ddmix_kernel (double *VDC, double *TRCR, int bid);
    __global__ void buoydiff_kernel (double *DBLOC, double *DBSFC, double *TRCR, int bid, int *KMT);
    __global__ void smooth_hblt_kernel(int bid, double *HBLT, int *KBL, double *d_WORK3, int *KMT, double *DZT);

    __device__ double wmscale (double sigma, double hbl, double ustar, double bfsfc);
    __device__ double wsscale (double sigma, double hbl, double ustar, double bfsfc);
    __device__ double ugrid_to_tgrid (double *ARRAY_UGRID, int i, int j, int k);
    __device__ double tgrid_to_ugrid(double *ARRAY_TGRID, double *AU, double *UAREA_R, int iblock, int i, int j, int k);
    __device__ double sw_absorb_frac(double depth);
    __device__ double compute_rho(double temp, double salt, int k);
    __device__ double mwjf_numerator(double tq, double sq, int k);
    __device__ double mwjf_denominator(double tq, double sq, double sqr, int k);
    __device__ double compute_drhodt(double tq, double sq, double sqr, int k, double nomk, double denomk);
    __device__ double compute_drhods(double tq, double sq, double sqr, int k, double nomk, double denomk);


    void init_global_variables( double *h_DZT, int *h_KMU, double *h_dz, double *h_zt, double *h_DZU, int *h_KMT, double *h_bckgrnd_vdc, double *h_bckgrnd_vvc, double *h_zgrid, double *h_Ricr, double *h_hwide, double *pressz, double *h_AU, double *h_UAREA_R);
    void vmix_coeffs_kpp_gpu_entry (double *VDC, double *VVC, double *TRCR, double *UUU, double *VVV, double *STF, double *SHF_QSW, int *pbid, double *pconvect_diff, double *pconvect_visc, double *SMF, double *HMXL, double *KPP_HBLT, double *KPP_SRC);

    void vmix_coeffs_kpp_gpu_entry_test (double *VDC, double *VVC, double *TRCR, double *UUU, double *VVV, double *STF, double *SHF_QSW, int *pbid, double *pconvect_diff, double *pconvect_visc, double *SMF, double *HMXL, double *KPP_HBLT, double *KPP_SRC);

    void interior_convection (double *DBLOC, double *DBSFC, int *KBL, double *STF, double *GHAT, double *VISC, double *convect_diff, double *convect_visc, int *bid, double *VDC, double *VVC, double *KPP_SRC, double *HMXL);
    void ri_iwmix (double *DBLOC, double *VISC, double *VDC, double *UUU, double *VVV, double *RHOMIX, double *convect_diff, double *convect_visc, int *bid);
    void blmix (double *VISC, double *VDC, double *HBLT, double *USTAR, double *BFSFC, double *STABLE, int *KBL, double *GHAT, int *bid);
    void ddmix (double *VDC, double *TRCR, int *bid);

    void buoydiff(double *DBLOC, double *DBSFC, double *TRCR, int *bid);

//    void init_vmix_kpp( double *h_DZT, int *h_KMU, double *h_dz, double *h_HMXL, double *h_zt, double *h_DZU, int *h_KMT, double *h_bckgrnd_vdc, double *h_bckgrnd_vvc, double *h_zgrid, double *h_Ricr, double *h_hwide, double *pressz, double *h_AU);

    void bldepth_test (double *DBLOC, double *DBSFC, double *TRCR, double *UUU, double *VVV, double *STF, double *SHF_QSW, double *HBLT, double *USTAR, double *BFSFC, double *STABLE, int *KBL, int *bid, double *SMF);
//    void bldepth (double *DBLOC, double *DBSFC, double *TRCR, double *UUU, double *VVV, double *STF, double *SHF_QSW, double *HBLT, double *USTAR, double *BFSFC, double *STABLE, int *KBL, int *bid, double *SMF);

}





/*
double      c0     =    0.0   ;
double      c1     =    1.0   ;
double      c2     =    2.0   ;
double      c3     =    3.0   ;
double      c4     =    4.0   ;
double      c5     =    5.0   ;
double      c8     =    8.0   ;
double      c10    =   10.0   ;
double      c16    =   16.0   ;
double      c1000  = 1000.0   ;
double      c10000 =10000.0   ;
double      c1p5   =    1.5   ;
double      p33    = c1/c3    ;
double      p5     = 0.500    ;
double      p25    = 0.250    ;
double      p125   = 0.125    ;
double      p001   = 0.001    ;
double      eps    = 1.0e-10  ;
double      eps2   = 1.0e-20  ;
double      bignum = 1.0e+30  ;
double      grav   = 980.6    ;
*/
    __constant__ double c0 = 0.0;
    __constant__ double c1 = 1.0;
    __constant__ double c2 = 2.0;
    __constant__ double c3 = 3.0;
    __constant__ double c4 = 4.0;
    __constant__ double c5 = 5.0;
    __constant__ double c8 = 8.0;
    __constant__ double c10 = 10.0;
    __constant__ double c16 = 16.0;
    __constant__ double c1000 = 1000.0;
    __constant__ double c10000 = 10000.0;
    __constant__ double c1p5 = 1.5;
    __constant__ double p33 = 1.0/3.0;
    __constant__ double p5 = 0.500;
    __constant__ double p25 = 0.250;
    __constant__ double p125 = 0.125;
    __constant__ double p001 = 0.001;
    __constant__ double eps = 1.0e-10;
    __constant__ double eps2 = 1.0e-20;
    __constant__ double bignum = 1.0e+30;

    __constant__ double grav = 980.6;
    __constant__ double vonkar = 0.4;

int *KMT = (int *) NULL;
double *DZU = (double *) NULL;
double *HMXL = (double *) NULL;
int *KMU = (int *) NULL;
double *DZT = (double *) NULL;
double *AU = (double *) NULL;
double *UAREA_R = (double *) NULL;

__constant__ double hwide[km+2];
__constant__ double zgrid[km+2];
__constant__ double Ricr[km];
__constant__ double zt[km];
__constant__ double dz[km];

__constant__ double bckgrnd_vvc[km];
__constant__ double bckgrnd_vdc[km];

extern int my_task;

void init_global_variables( double *h_DZT, int *h_KMU, double *h_dz, double *h_zt, double *h_DZU, int *h_KMT, double *h_bckgrnd_vdc, double *h_bckgrnd_vvc, double *h_zgrid, double *h_Ricr, double *h_hwide, double *pressz, double *h_AU, double *h_UAREA_R) {

  cudaError_t err;

  err = cudaStreamCreate(stream);
  if (err != cudaSuccess) { fprintf (stderr, "Error in cuda stream creation: %s\n", cudaGetErrorString (err)); }


  err = cudaMalloc ((void **) &DZT, nx_block * ny_block * (km+2) * max_blocks_clinic * sizeof(double));
  if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMalloc DZT: %s\n", cudaGetErrorString (err)); }

  err = cudaMemcpyAsync(DZT, h_DZT, nx_block * ny_block * (km+2) * max_blocks_clinic * sizeof(double), cudaMemcpyHostToDevice, 0);
  if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMemcpy host to device DZT: %s\n", cudaGetErrorString(err)); }

  err = cudaMalloc ((void **) &KMU, nx_block * ny_block * max_blocks_clinic * sizeof(int));
  if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMalloc KMU: %s\n", cudaGetErrorString (err)); }

  err = cudaMemcpyAsync(KMU, h_KMU, nx_block * ny_block * max_blocks_clinic * sizeof(int), cudaMemcpyHostToDevice, 0);
  if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMemcpy host to device KMU: %s\n", cudaGetErrorString(err)); }

  err = cudaMemcpyToSymbol(dz, h_dz, km * sizeof(double), 0, cudaMemcpyHostToDevice);
  if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMemcpy host to device dz: %s\n", cudaGetErrorString(err)); }

  err = cudaMemcpyToSymbol(zt, h_zt, km * sizeof(double), 0, cudaMemcpyHostToDevice);
  if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMemcpy host to device zt: %s\n", cudaGetErrorString(err)); }

  err = cudaMalloc ((void **) &DZU, nx_block * ny_block * (km+2) * max_blocks_clinic * sizeof(double));
  if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMalloc DZU: %s\n", cudaGetErrorString (err)); }

  err = cudaMemcpyAsync(DZU, h_DZU, nx_block * ny_block * (km+2) * max_blocks_clinic * sizeof(double), cudaMemcpyHostToDevice, 0);
  if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMemcpy host to device DZU: %s\n", cudaGetErrorString(err)); }

  err = cudaMalloc ((void **) &KMT, nx_block * ny_block * max_blocks_clinic * sizeof(int));
  if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMalloc KMT: %s\n", cudaGetErrorString (err)); }

  err = cudaMemcpyAsync(KMT, h_KMT, nx_block * ny_block * max_blocks_clinic * sizeof(int), cudaMemcpyHostToDevice, 0);
  if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMemcpy host to device KMT: %s\n", cudaGetErrorString(err)); }

//  err = cudaMalloc ((void **) &bckgrnd_vdc, nx_block * ny_block * km * nblocks_clinic * sizeof(double));
  err = cudaMalloc ((void **) &bckgrnd_vdc, km * sizeof(double));
  if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMalloc bckgrnd_vdc: %s\n", cudaGetErrorString (err)); }

//  err = cudaMemcpyAsync(bckgrnd_vdc, h_bckgrnd_vdc, nx_block * ny_block * km * nblocks_clinic * sizeof(double), cudaMemcpyHostToDevice, 0);
  double background[km];
  for (int k = 0; k<km; k++) {
    background[k] = h_bckgrnd_vdc[(k)*nx_block*ny_block];
  }

//  err = cudaMemcpyAsync(bckgrnd_vdc, h_bckgrnd_vdc, nx_block * ny_block * km * nblocks_clinic * sizeof(double), cudaMemcpyHostToDevice, 0);
  err = cudaMemcpyToSymbol(bckgrnd_vdc, background, km * sizeof(double), 0, cudaMemcpyHostToDevice);
  if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMemcpy host to device bckgrnd_vdc: %s\n", cudaGetErrorString(err)); }

//  err = cudaMalloc ((void **) &bckgrnd_vvc, nx_block * ny_block * km * nblocks_clinic * sizeof(double));
  err = cudaMalloc ((void **) &bckgrnd_vvc, km * sizeof(double));
  if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMalloc bckgrnd_vvc: %s\n", cudaGetErrorString (err)); }

  for (int k = 0; k<km; k++) {
    background[k] = h_bckgrnd_vvc[(k)*nx_block*ny_block];
  }

  err = cudaMemcpyToSymbol(bckgrnd_vvc, background, km * sizeof(double), 0, cudaMemcpyHostToDevice);
  if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMemcpy host to device bckgrnd_vvc: %s\n", cudaGetErrorString(err)); }

  err = cudaMemcpyToSymbol(zgrid, h_zgrid, (km+2) * sizeof(double), 0, cudaMemcpyHostToDevice);
  if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMemcpy host to device zgrid: %s\n", cudaGetErrorString(err)); }

  err = cudaMemcpyToSymbol(Ricr, h_Ricr, km * sizeof(double), 0, cudaMemcpyHostToDevice);
  if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMemcpy host to device Ricr: %s\n", cudaGetErrorString(err)); }

  err = cudaMemcpyToSymbol(hwide, h_hwide, (km+2) * sizeof(double), 0, cudaMemcpyHostToDevice);
  if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMemcpy host to device hwide: %s\n", cudaGetErrorString(err)); }

  err = cudaMalloc ((void **) &AU, nx_block * ny_block * max_blocks_clinic * sizeof(double));
  if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMalloc DZU: %s\n", cudaGetErrorString (err)); }

  err = cudaMemcpyAsync(AU, h_AU, nx_block * ny_block * max_blocks_clinic * sizeof(double), cudaMemcpyHostToDevice, 0);
  if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMemcpy host to device DZU: %s\n", cudaGetErrorString(err)); }

  err = cudaMalloc ((void **) &UAREA_R, nx_block * ny_block * max_blocks_clinic * sizeof(double));
  if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMalloc DZU: %s\n", cudaGetErrorString (err)); }

  err = cudaMemcpyAsync(UAREA_R, h_UAREA_R, nx_block * ny_block * max_blocks_clinic * sizeof(double), cudaMemcpyHostToDevice, 0);
  if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMemcpy host to device DZU: %s\n", cudaGetErrorString(err)); }

  //constants for blmix
  double cstar = 10.0;
  double vonkar = 0.4;
  double c_s = 98.96;
  double epssfc = 0.1;
  double h_cg = cstar*vonkar* pow(c_s*vonkar*epssfc, 1.0/3.0);
  double h_Vtc = sqrt (0.2 / c_s / epssfc) / pow(vonkar,2);

  err = cudaMemcpyToSymbol(cg, &h_cg, sizeof(double), 0, cudaMemcpyHostToDevice);
  if (err != cudaSuccess) { fprintf(stderr, "Error doing cudaMemcpyToSymbol cg %s\n", cudaGetErrorString( err )); }
  err = cudaMemcpyToSymbol(Vtc, &h_Vtc, sizeof(double), 0, cudaMemcpyHostToDevice);
  if (err != cudaSuccess) { fprintf(stderr, "Error doing cudaMemcpyToSymbol Vtc %s\n", cudaGetErrorString( err )); }

//printf("cg=%.16f\n", h_cg);
//printf("Vtc=%.16f\n", h_Vtc);

  //init phase for equation of state
  //initialize all constant arrays to be stored in constant memory on the GPU
  double p;

  err = cudaMemcpyToSymbol(mwjfnums0t1, &mwjfnp0s0t1, sizeof(double), 0, cudaMemcpyHostToDevice); //= mwjfnp0s0t1
  if (err != cudaSuccess) { fprintf(stderr, "Error doing cudaMemcpyToSymbol mwjfnums0t1 %s\n", cudaGetErrorString( err )); }
  err = cudaMemcpyToSymbol(mwjfnums0t3, &mwjfnp0s0t3, sizeof(double), 0, cudaMemcpyHostToDevice); //= mwjfnp0s0t3
  if (err != cudaSuccess) { fprintf(stderr, "Error doing cudaMemcpyToSymbol mwjfnums0t3 %s\n", cudaGetErrorString( err )); }
  err = cudaMemcpyToSymbol(mwjfnums1t1, &mwjfnp0s1t1, sizeof(double), 0, cudaMemcpyHostToDevice); //= mwjfnp0s1t1
  if (err != cudaSuccess) { fprintf(stderr, "Error doing cudaMemcpyToSymbol mwjfnums1t1 %s\n", cudaGetErrorString( err )); }
  err = cudaMemcpyToSymbol(mwjfnums2t0, &mwjfnp0s2t0, sizeof(double), 0, cudaMemcpyHostToDevice); //= mwjfnp0s2t0
  if (err != cudaSuccess) { fprintf(stderr, "Error doing cudaMemcpyToSymbol mwjfnums2t0 %s\n", cudaGetErrorString( err )); }
  err = cudaMemcpyToSymbol(mwjfdens0t2, &mwjfdp0s0t2, sizeof(double), 0, cudaMemcpyHostToDevice); //= mwjfdp0s0t2
  if (err != cudaSuccess) { fprintf(stderr, "Error doing cudaMemcpyToSymbol mwjfdens0t2 %s\n", cudaGetErrorString( err )); }
  err = cudaMemcpyToSymbol(mwjfdens0t4, &mwjfdp0s0t4, sizeof(double), 0, cudaMemcpyHostToDevice); //= mwjfdp0s0t4
  if (err != cudaSuccess) { fprintf(stderr, "Error doing cudaMemcpyToSymbol mwjfdens0t4 %s\n", cudaGetErrorString( err )); }
  err = cudaMemcpyToSymbol(mwjfdens1t0, &mwjfdp0s1t0, sizeof(double), 0, cudaMemcpyHostToDevice); //= mwjfdp0s1t0
  if (err != cudaSuccess) { fprintf(stderr, "Error doing cudaMemcpyToSymbol mwjfdens1t0 %s\n", cudaGetErrorString( err )); }
  err = cudaMemcpyToSymbol(mwjfdens1t1, &mwjfdp0s1t1, sizeof(double), 0, cudaMemcpyHostToDevice); //= mwjfdp0s1t1
  if (err != cudaSuccess) { fprintf(stderr, "Error doing cudaMemcpyToSymbol mwjfdens1t1 %s\n", cudaGetErrorString( err )); }
  err = cudaMemcpyToSymbol(mwjfdens1t3, &mwjfdp0s1t3, sizeof(double), 0, cudaMemcpyHostToDevice); //= mwjfdp0s1t3
  if (err != cudaSuccess) { fprintf(stderr, "Error doing cudaMemcpyToSymbol mwjfdens1t3 %s\n", cudaGetErrorString( err )); }
  err = cudaMemcpyToSymbol(mwjfdensqt0, &mwjfdp0sqt0, sizeof(double), 0, cudaMemcpyHostToDevice); //= mwjfdp0sqt0
  if (err != cudaSuccess) { fprintf(stderr, "Error doing cudaMemcpyToSymbol mwjfdensqt0 %s\n", cudaGetErrorString( err )); }
  err = cudaMemcpyToSymbol(mwjfdensqt2, &mwjfdp0sqt2, sizeof(double), 0, cudaMemcpyHostToDevice); //= mwjfdp0sqt2
  if (err != cudaSuccess) { fprintf(stderr, "Error doing cudaMemcpyToSymbol mwjfdensqt2 %s\n", cudaGetErrorString( err )); }

  double h_mwjfnums0t0[km];
  double h_mwjfnums0t2[km];
  double h_mwjfnums1t0[km];
  double h_mwjfdens0t0[km];
  double h_mwjfdens0t1[km];
  double h_mwjfdens0t3[km];

  int k;
  for (k=0; k<km; k++) {
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

  err = cudaMemcpyToSymbol(mwjfnums0t0, h_mwjfnums0t0, km*sizeof(double), 0, cudaMemcpyHostToDevice);
  if (err != cudaSuccess) { fprintf(stderr, "Error doing cudaMemcpyToSymbol mwjfnums0t0 %s\n", cudaGetErrorString( err )); }
  err = cudaMemcpyToSymbol(mwjfnums0t2, h_mwjfnums0t2, km*sizeof(double), 0, cudaMemcpyHostToDevice);
  if (err != cudaSuccess) { fprintf(stderr, "Error doing cudaMemcpyToSymbol mwjfnums0t2 %s\n", cudaGetErrorString( err )); }
  err = cudaMemcpyToSymbol(mwjfnums1t0, h_mwjfnums1t0, km*sizeof(double), 0, cudaMemcpyHostToDevice);
  if (err != cudaSuccess) { fprintf(stderr, "Error doing cudaMemcpyToSymbol mwjfnums1t0 %s\n", cudaGetErrorString( err )); }
  cudaMemcpyToSymbol(mwjfdens0t0, h_mwjfdens0t0, km*sizeof(double), 0, cudaMemcpyHostToDevice);
  if (err != cudaSuccess) { fprintf(stderr, "Error doing cudaMemcpyToSymbol mwjfdens0t0 %s\n", cudaGetErrorString( err )); }
  cudaMemcpyToSymbol(mwjfdens0t1, h_mwjfdens0t1, km*sizeof(double), 0, cudaMemcpyHostToDevice);
  if (err != cudaSuccess) { fprintf(stderr, "Error doing cudaMemcpyToSymbol mwjfdens0t1 %s\n", cudaGetErrorString( err )); }
  cudaMemcpyToSymbol(mwjfdens0t3, h_mwjfdens0t3, km*sizeof(double), 0, cudaMemcpyHostToDevice);
  if (err != cudaSuccess) { fprintf(stderr, "Error doing cudaMemcpyToSymbol mwjfdens0t3 %s\n", cudaGetErrorString( err )); }



  cudaDeviceSynchronize();
  err = cudaGetLastError();
  if (err != cudaSuccess) { fprintf (stderr, "Error in init global variables vmix KPP GPU: %s\n", cudaGetErrorString (err)); }


  if (my_task == 0) {
    printf("Initialized KPP GPU version using thread block size: %d x %d\n", BLOCK_X, BLOCK_Y);
  }


}







  //module vmix_kpp
    //! BOP
    //! MODULE: vmix_kpp
    //!
    //! DESCRIPTION:
    //! This module contains routines for initializing and computing
    //! vertical mixing coefficients for the KPP parameterization
    //! (see Large, McWilliams and Doney, Reviews of Geophysics, 32, 363
    //! November 1994) .
    //!
    //! REVISION HISTORY:
    //! SVN:$Id: vmix_kpp.F90 12674 2008-10-31 22:21:32Z njn01 $
    //! USES:
    //use POP_IOUnitsMod
    //use kinds_mod
    //use blocks
    //use distribution
    //use domain_size
    //use domain
    //use constants
    //use grid
    //use broadcast
    //use io
    //use state_mod
    //use exit_mod
    //use sw_absorption
    //use tavg, only: define_tavg_field, tavg_requested, accumulate_tavg_field
    //use io_types, only: stdout
    //use communicate, only: my_task, master_task
    //use tidal_mixing, only: TIDAL_COEF, tidal_mix_max, ltidal_mixing
    //use iso_c_binding
    //implicit none
    //private
    //save
    //! PUBLIC MEMBER FUNCTIONS:
    //public :: init_vmix_kpp, vmix_coeffs_kpp, vmix_coeffs_kpp_gpu, add_kpp_sources, smooth_hblt, linertial
    //! PUBLIC DATA MEMBERS:
    double *KPP_HBLT;
    double *BOLUS_SP;
    //! mixing parameterization
    __constant__ double bckgrnd_vdc2;
    //! EOP
    //! BOC
    //! -----------------------------------------------------------------------
    //!
    //! mixing constants
    //!
    //! -----------------------------------------------------------------------
    __constant__ double rich_mix = 50.0;
    __constant__ int num_v_smooth_Ri = 1;
    __constant__ double epssfc = 0.1;
    __constant__ double Prandtl = 10.0;
    double *FSTOKES;
    //! parameterization
    //! -----------------------------------------------------------------------
    //!
    //! non-local mixing (counter-gradient mixing) , treated as source term
    //!
    //! -----------------------------------------------------------------------
    double *KPP_SRC;
    //! -----------------------------------------------------------------------
    //!
    //! parameters for subroutine bldepth: computes bndy layer depth
    //!
    //! concv = ratio of interior buoyancy frequency to
    //! buoyancy frequency at entrainment depth
    //! parameter statement sets the minimum value.
    //!
    //! -----------------------------------------------------------------------
    //! as a function of vertical resolution
    __constant__ double cekman = 0.7;
    __constant__ double cmonob = 1.0;
    __constant__ double concv = 1.7;
    __constant__ double hbf = 1.0;
    //! -----------------------------------------------------------------------
    //!
    //! parameters for subroutine ri_iwmix which computes
    //! vertical mixing coefficients below boundary layer due to shear
    //! instabilities, internal waves and convection
    //!
    //! -----------------------------------------------------------------------
    __constant__ double Riinfty = 0.8;
    __constant__ double BVSQcon = 0.0;
    //! -----------------------------------------------------------------------
    //!
    //! parameters for subroutine ddmix (double-diffusive mixing)
    //!
    //! -----------------------------------------------------------------------
    __constant__ double Rrho0 = 2.55;
    __constant__ double dsfmax = 1.0;
    //! turbulent velocity shear coefficient (for bulk Ri no)
    //! -----------------------------------------------------------------------
    //!
    //! parameters for velocity scale function (from Large et al.)
    //!
    //! -----------------------------------------------------------------------
    __constant__ double zeta_m = -0.2;
    __constant__ double zeta_s = -1.0;
    __constant__ double c_m = 8.38;
    __constant__ double c_s = 98.96;
    __constant__ double a_m = 1.26;
    __constant__ double a_s = -28.86;
    //! -----------------------------------------------------------------------
    //!
    //! parameters for subroutine blmix: mixing within boundary layer
    //!
    //! cstar = proportionality coefficient for nonlocal transport
    //! cg = non-dimensional coefficient for counter-gradient term
    //!
    //! -----------------------------------------------------------------------
    __constant__ double cstar = 10.0;
    //! -----------------------------------------------------------------------
    //!
    //! common vertical grid arrays used by KPP mixing
    //!
    //! -----------------------------------------------------------------------
    //! -----------------------------------------------------------------------
    //!
    //! ids for tavg diagnostics computed in this module
    //!
    //! -----------------------------------------------------------------------
  
    //! EOC
    //! ***********************************************************************
    //contains
    //! ***********************************************************************
    //! BOP
    //! IROUTINE: init_vmix_kpp
    //! INTERFACE:
    //! ***********************************************************************
    //! BOP
    //! IROUTINE: vmix_coeffs_kpp
    //! INTERFACE:
    //! ***********************************************************************
    //! BOP
    //! IROUTINE: vmix_coeffs_kpp_gpu
    //! INTERFACE:
    //! ***********************************************************************
    //! BOP
    //! IROUTINE: vmix_coeffs_kpp_gpu_entry
    //! INTERFACE:





    void vmix_coeffs_kpp_gpu_entry (double *VDC, double *VVC, double *TRCR, double *UUU, double *VVV, double *STF, double *SHF_QSW, int *pbid, double *pconvect_diff, double *pconvect_visc, double *SMF, double *HMXL, double *KPP_HBLT, double *KPP_SRC) {

      cudaError_t err;

      
      //! DESCRIPTION:
      //! This is the main GPU driver routine which calculates the vertical
      //! mixing coefficients for the KPP mixing scheme as outlined in
      //! Large, McWilliams and Doney, Reviews of Geophysics, 32, 363
      //! (November 1994) . The non-local mixing is also computed here, but
      //! is treated as a source term in baroclinic. This routine will be
      //! translated to C.
      //!
      //! REVISION HISTORY:
      //! same as module
      //! INPUT PARAMETERS:
      static double *d_TRCR = (double *) NULL;
      if (d_TRCR == NULL) {
        err = cudaMalloc ((void **) &d_TRCR, nx_block * ny_block * km * nt * sizeof(double));
        if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMalloc d_TRCR: %s\n", cudaGetErrorString (err)); }
      }

      static double *d_UUU = (double *) NULL;
      if (d_UUU == NULL) {
        err = cudaMalloc ((void **) &d_UUU, nx_block * ny_block * km * sizeof(double));
        if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMalloc d_UUU: %s\n", cudaGetErrorString (err)); }
      }

      static double *d_VVV = (double *) NULL;
      if (d_VVV == NULL) {
        err = cudaMalloc ((void **) &d_VVV, nx_block * ny_block * km * sizeof(double));
        if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMalloc d_VVV: %s\n", cudaGetErrorString (err)); }
      }

      static double *d_STF = (double *) NULL;
      if (d_STF == NULL) {
        err = cudaMalloc ((void **) &d_STF, nx_block * ny_block * nt * sizeof(double));
        if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMalloc d_STF: %s\n", cudaGetErrorString (err)); }
      }

      static double *d_SHF_QSW = (double *) NULL;
      if (d_SHF_QSW == NULL) {
        err = cudaMalloc ((void **) &d_SHF_QSW, nx_block * ny_block * sizeof(double));
        if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMalloc d_SHF_QSW: %s\n", cudaGetErrorString (err)); }
      }

      static double *d_SMF = (double *) NULL;
      if (d_SMF == NULL) {
        err = cudaMalloc ((void **) &d_SMF, nx_block * ny_block * 2 * sizeof(double));
        if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMalloc d_SMF: %s\n", cudaGetErrorString (err)); }
      }

      double convect_diff = *pconvect_diff;
      double convect_visc = *pconvect_visc;
//      int fbid = *pbid;
      int bid = *pbid-1;//adjusted to c
      //! INPUT/OUTPUT PARAMETERS:
      static double *d_VDC = (double *) NULL;
      if (d_VDC == NULL) {
        err = cudaMalloc ((void **) &d_VDC, nx_block * ny_block * (km+1) * 2 * sizeof(double));
        if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMalloc d_VDC: %s\n", cudaGetErrorString (err)); }
        err = cudaMemset (d_VDC, 0, nx_block * ny_block * (km+1) * 2 * sizeof(double));
        if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMemset d_VDC: %s\n", cudaGetErrorString (err)); }
      }

      //! OUTPUT PARAMETERS:
      static double *d_VVC = (double *) NULL;
      if (d_VVC == NULL) {
        err = cudaMalloc ((void **) &d_VVC, nx_block * ny_block * km * sizeof(double));
        if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMalloc d_VVC: %s\n", cudaGetErrorString (err)); }
      }

      static double *d_HMXL = (double *) NULL;
      if (d_HMXL == NULL) {
        err = cudaMalloc ((void **) &d_HMXL, nx_block * ny_block * sizeof(double));
        if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMalloc d_HMXL: %s\n", cudaGetErrorString (err)); }
        err = cudaMemset (d_HMXL, 0, nx_block * ny_block * sizeof(double));
        if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMemset d_HMXL: %s\n", cudaGetErrorString (err)); }
      }

      static double *d_KPP_HBLT = (double *) NULL;
      if (d_KPP_HBLT == NULL) {
        err = cudaMalloc ((void **) &d_KPP_HBLT, nx_block * ny_block * sizeof(double));
        if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMalloc d_KPP_HBLT: %s\n", cudaGetErrorString (err)); }
      }

/*      static double *d_KPP_SRC = (double *) NULL;
      if (d_KPP_SRC == NULL) {
        err = cudaMalloc ((void **) &d_KPP_SRC, nx_block * ny_block * km * nt * max_blocks_clinic * sizeof(double));
        if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMalloc d_KPP_SRC: %s\n", cudaGetErrorString (err)); }
        err = cudaMemset (d_KPP_SRC, 0, nx_block * ny_block * km * nt * max_blocks_clinic * sizeof(double));
        if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMemset d_KPP_SRC: %s\n", cudaGetErrorString (err)); }
      }
*/
      //! EOP
      //! BOC
      //! -----------------------------------------------------------------------
      //!
      //! local variables
      //!
      //! -----------------------------------------------------------------------
      static int *d_KBL = (int *) NULL;
      if (d_KBL == NULL) {
        err = cudaMalloc ((void **) &d_KBL, nx_block * ny_block * sizeof(int));
        if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMalloc d_KBL: %s\n", cudaGetErrorString (err)); }
      }

      static double *d_USTAR = (double *) NULL;
      if (d_USTAR == NULL) {
        err = cudaMalloc ((void **) &d_USTAR, nx_block * ny_block * sizeof(double));
        if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMalloc d_USTAR: %s\n", cudaGetErrorString (err)); }
      }

      static double *d_BFSFC = (double *) NULL;
      if (d_BFSFC == NULL) {
        err = cudaMalloc ((void **) &d_BFSFC, nx_block * ny_block * sizeof(double));
        if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMalloc d_BFSFC: %s\n", cudaGetErrorString (err)); }
      }

      static double *d_STABLE = (double *) NULL;
      if (d_STABLE == NULL) {
        err = cudaMalloc ((void **) &d_STABLE, nx_block * ny_block * sizeof(double));
        if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMalloc d_STABLE: %s\n", cudaGetErrorString (err)); }
      }

      static double *d_DBLOC = (double *) NULL;
      if (d_DBLOC == NULL) {
        err = cudaMalloc ((void **) &d_DBLOC, nx_block * ny_block * km * sizeof(double));
        if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMalloc d_DBLOC: %s\n", cudaGetErrorString (err)); }
      }

      static double *d_DBSFC = (double *) NULL;
      if (d_DBSFC == NULL) {
        err = cudaMalloc ((void **) &d_DBSFC, nx_block * ny_block * km * sizeof(double));
        if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMalloc d_DBSFC: %s\n", cudaGetErrorString (err)); }
      }

      static double *d_GHAT = (double *) NULL;
      if (d_GHAT == NULL) {
        err = cudaMalloc ((void **) &d_GHAT, nx_block * ny_block * km * sizeof(double));
        if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMalloc d_GHAT: %s\n", cudaGetErrorString (err)); }
      }

      static double *d_VISC = (double *) NULL;
      if (d_VISC == NULL) {
        err = cudaMalloc ((void **) &d_VISC, nx_block * ny_block * (km+1) * sizeof(double));
        if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMalloc d_VISC: %s\n", cudaGetErrorString (err)); }
      }

      //! dummy variable that is not used by ri_iwmix but added here to have the
      //! same interface during testing
      static double *d_RHOMIX = (double *) NULL;
// rhomix is not used
//      if (d_RHOMIX == NULL) {
//        err = cudaMalloc ((void **) &d_RHOMIX, nx_block * ny_block * km * sizeof(double));
//        if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMalloc d_RHOMIX: %s\n", cudaGetErrorString (err)); }
//      }
//
      //! temporary variable for storing velocity shear squared on u grid, converted to t grid by ri_iwmix
      static double *d_VSHEARU = (double *) NULL;
      if (d_VSHEARU == NULL) {
        err = cudaMalloc ((void **) &d_VSHEARU, nx_block * ny_block * km * sizeof(double));
        if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMalloc d_VSHEARU: %s\n", cudaGetErrorString (err)); }
      }
      static double *d_WORK1 = (double *) NULL;
      if (d_WORK1 == NULL) {
        err = cudaMalloc ((void **) &d_WORK1, nx_block * ny_block * sizeof(double));
        if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMalloc d_WORK1: %s\n", cudaGetErrorString (err)); }
      }
      static double *d_WORK2 = (double *) NULL;
      if (d_WORK2 == NULL) {
        err = cudaMalloc ((void **) &d_WORK2, nx_block * ny_block * sizeof(double));
        if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMalloc d_WORK2: %s\n", cudaGetErrorString (err)); }
      }
      static double *d_WORK3 = (double *) NULL;
      if (d_WORK3 == NULL) {
        err = cudaMalloc ((void **) &d_WORK3, nx_block * ny_block * sizeof(double));
        if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMalloc d_WORK3: %s\n", cudaGetErrorString (err)); }
      }

      //copy inputs from Host to Device
      err = cudaMemcpyAsync(d_TRCR, TRCR, nx_block * ny_block * km * nt * sizeof(double), cudaMemcpyHostToDevice, stream[0]);
      if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMemcpy host to device d_TRCR: %s\n", cudaGetErrorString(err)); }

      err = cudaMemcpyAsync(d_UUU, UUU, nx_block * ny_block * km * sizeof(double), cudaMemcpyHostToDevice, stream[0]);
      if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMemcpy host to device d_UUU: %s\n", cudaGetErrorString(err)); }
      err = cudaMemcpyAsync(d_VVV, VVV, nx_block * ny_block * km * sizeof(double), cudaMemcpyHostToDevice, stream[0]);
      if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMemcpy host to device d_VVV: %s\n", cudaGetErrorString(err)); }

      err = cudaMemcpyAsync(d_STF, STF, nx_block * ny_block * nt * sizeof(double), cudaMemcpyHostToDevice, stream[0]);
      if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMemcpy host to device d_STF: %s\n", cudaGetErrorString(err)); }
      err = cudaMemcpyAsync(d_SHF_QSW, SHF_QSW, nx_block * ny_block * sizeof(double), cudaMemcpyHostToDevice, stream[0]);
      if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMemcpy host to device d_SHF_QSW: %s\n", cudaGetErrorString(err)); }
      err = cudaMemcpyAsync(d_SMF, SMF, nx_block * ny_block * 2 * sizeof(double), cudaMemcpyHostToDevice, stream[0]);
      if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMemcpy host to device d_SMF: %s\n", cudaGetErrorString(err)); }

//      err = cudaMemcpyAsync(d_VDC, VDC, nx_block * ny_block * (km+1) * 2 * sizeof(double), cudaMemcpyHostToDevice, stream[0]);
//      if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMemcpy host to device d_VDC: %s\n", cudaGetErrorString(err)); }

//      double *DEBUG = (double *) NULL;
//      float time;

      //setup threads and grid
      dim3 threads(BLOCK_X, BLOCK_Y);
      dim3 grid( (int)ceilf((float) nx_block / (float)(BLOCK_X)) , (int)ceilf((float) ny_block / (float)(BLOCK_Y)));

      //! -----------------------------------------------------------------------
      //!
      //! compute buoyancy differences at each vertical level.
      //!
      //! -----------------------------------------------------------------------
      buoydiff_kernel<<<grid, threads, 0, stream[0]>>>(d_DBLOC, d_DBSFC, d_TRCR, bid, KMT);

//      err = cudaPeekAtLastError();
//      if (err != cudaSuccess) { fprintf (stderr, "Error in buoydiff kernel: %s\n", cudaGetErrorString (err)); }

      //! -----------------------------------------------------------------------
      //!
      //! compute mixing due to shear instability, internal waves and
      //! convection
      //!
      //! -----------------------------------------------------------------------
      compute_vshear_ugrid<<<grid, threads, 0, stream[0]>>>(d_VSHEARU, d_UUU, d_VVV, DZU, bid);

      ri_iwmix_kernel<<<grid, threads, 0, stream[0]>>>(d_DBLOC, d_VISC, d_VDC, d_UUU, d_VVV, d_RHOMIX, convect_diff, convect_visc, bid, d_VSHEARU, DZT, KMT);

//      err = cudaPeekAtLastError();
//      if (err != cudaSuccess) { fprintf (stderr, "Error in ri_iwmix kernel: %s\n", cudaGetErrorString (err)); }

      //! -----------------------------------------------------------------------
      //!
      //! compute double diffusion if desired
      //!
      //! -----------------------------------------------------------------------
      ddmix_kernel<<<grid, threads, 0, stream[0]>>>(d_VDC, d_TRCR, bid);
   
//      err = cudaPeekAtLastError();
//      if (err != cudaSuccess) { fprintf (stderr, "Error in ddmix kernel: %s\n", cudaGetErrorString (err)); }

      //! -----------------------------------------------------------------------
      //!
      //! compute boundary layer depth
      //!
      //! -----------------------------------------------------------------------

      compute_ustaru<<<grid, threads, 0, stream[0]>>>(d_STABLE, d_SMF); //storing temporary result in stable
      ugrid_to_tgrid_kernel<<<grid, threads, 0, stream[0]>>>(d_USTAR, d_STABLE, 0);
      compute_vshearu_bldepth<<<grid, threads, 0, stream[0]>>>(d_VSHEARU, d_UUU, d_VVV, DZU, bid);
      bldepth_kernel<<<grid, threads, 0, stream[0]>>>(d_DBLOC, d_DBSFC, d_TRCR, d_UUU, d_VVV, d_STF, d_SHF_QSW, d_KPP_HBLT, d_USTAR, d_BFSFC, d_STABLE, d_KBL, bid, d_SMF, d_VSHEARU, d_WORK1, d_WORK2, d_WORK3, DZT, DZU, KMT);
      smooth_hblt_kernel<<<grid, threads, 0, stream[0]>>>(bid, d_KPP_HBLT, d_KBL, d_WORK3, KMT, DZT);
      correct_stability_and_buoyancy<<<grid, threads, 0, stream[0]>>>(d_BFSFC, d_STABLE, d_WORK1, d_WORK2, d_KPP_HBLT); 
    
//      err = cudaPeekAtLastError();
//      if (err != cudaSuccess) { fprintf (stderr, "Error in bldepth kernel: %s\n", cudaGetErrorString (err)); }
 
      err = cudaMemcpyAsync(KPP_HBLT, d_KPP_HBLT, nx_block * ny_block * sizeof(double), cudaMemcpyDeviceToHost, stream[0]);
      if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMemcpy device to host KPP_HBLT: %s\n", cudaGetErrorString (err)); }

      //! -----------------------------------------------------------------------
      //!
      //! compute boundary layer diffusivities
      //!
      //! -----------------------------------------------------------------------
      blmix_kernel<<<grid, threads, 0, stream[0]>>>(d_VISC, d_VDC, d_KPP_HBLT, d_USTAR, d_BFSFC, d_STABLE, d_KBL, d_GHAT, bid, DZT);

//      err = cudaPeekAtLastError();
//      if (err != cudaSuccess) { fprintf (stderr, "Error in blmix kernel: %s\n", cudaGetErrorString (err)); }
 
      //! -----------------------------------------------------------------------
      //!
      //! consider interior convection:
      //!
      //! -----------------------------------------------------------------------

      interior_convection_kernel<<<grid, threads, 0, stream[0]>>>(d_DBLOC, d_DBSFC, d_KBL, d_STF, d_GHAT, d_VISC, convect_diff, convect_visc, bid, d_VDC, d_VVC, KPP_SRC, d_HMXL, d_VSHEARU, DZT, KMT); //reuse VSHEARU to store VISC on T grid
      interior_convection_part2<<<grid, threads, 0, stream[0]>>>(d_DBLOC, d_DBSFC, d_KBL, d_STF, d_GHAT, d_VISC, convect_diff, convect_visc, bid, d_VDC, d_VVC, KPP_SRC, d_HMXL, AU, UAREA_R, d_VSHEARU, DZT, KMT, KMU);
    
//      err = cudaPeekAtLastError();
//      if (err != cudaSuccess) { fprintf (stderr, "Error in interior convection kernel: %s\n", cudaGetErrorString (err)); }
 
      //! -----------------------------------------------------------------------
      //! EOC
          //copy outputs from Device to Host
      err = cudaMemcpyAsync(VDC, d_VDC, nx_block * ny_block * (km+1) * 2 * sizeof(double), cudaMemcpyDeviceToHost, stream[0]);
      if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMemcpy device to host VDC: %s\n", cudaGetErrorString (err)); }
      err = cudaMemcpyAsync(VVC, d_VVC, nx_block * ny_block * km * sizeof(double), cudaMemcpyDeviceToHost, stream[0]);
      if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMemcpy device to host VVC: %s\n", cudaGetErrorString (err)); }
      err = cudaMemcpyAsync(HMXL, d_HMXL, nx_block * ny_block * sizeof(double), cudaMemcpyDeviceToHost, stream[0]);
      if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMemcpy device to host HMXL: %s\n", cudaGetErrorString (err)); }
//      err = cudaMemcpyAsync(KPP_SRC, d_KPP_SRC, nx_block * ny_block * km * nt * sizeof(double), cudaMemcpyDeviceToHost, stream[0]);
//      if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMemcpy device to host KPP_SRC: %s\n", cudaGetErrorString (err)); }

      cudaDeviceSynchronize();
      err = cudaGetLastError();
      if (err != cudaSuccess) { fprintf (stderr, "Error in GPU version of entrypoint: %s\n", cudaGetErrorString (err)); }

//      printf("finished execution of entrypoint\n");

    
}



/*
 * Test version of entry point goes here
 *
*/ 
    void vmix_coeffs_kpp_gpu_entry_test (double *VDC, double *VVC, double *TRCR, double *UUU, double *VVV, double *STF, double *SHF_QSW, int *pbid, double *pconvect_diff, double *pconvect_visc, double *SMF, double *HMXL, double *KPP_HBLT, double *KPP_SRC) {

      cudaError_t err;

      
      //! DESCRIPTION:
      //! This is the main GPU driver routine which calculates the vertical
      //! mixing coefficients for the KPP mixing scheme as outlined in
      //! Large, McWilliams and Doney, Reviews of Geophysics, 32, 363
      //! (November 1994) . The non-local mixing is also computed here, but
      //! is treated as a source term in baroclinic. This routine will be
      //! translated to C.
      //!
      //! REVISION HISTORY:
      //! same as module
      //! INPUT PARAMETERS:
      static double *d_TRCR = (double *) NULL;
      if (d_TRCR == NULL) {
        err = cudaMalloc ((void **) &d_TRCR, nx_block * ny_block * km * nt * sizeof(double));
        if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMalloc d_TRCR: %s\n", cudaGetErrorString (err)); }
      }

      static double *d_UUU = (double *) NULL;
      if (d_UUU == NULL) {
        err = cudaMalloc ((void **) &d_UUU, nx_block * ny_block * km * sizeof(double));
        if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMalloc d_UUU: %s\n", cudaGetErrorString (err)); }
      }

      static double *d_VVV = (double *) NULL;
      if (d_VVV == NULL) {
        err = cudaMalloc ((void **) &d_VVV, nx_block * ny_block * km * sizeof(double));
        if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMalloc d_VVV: %s\n", cudaGetErrorString (err)); }
      }

      static double *d_STF = (double *) NULL;
      if (d_STF == NULL) {
        err = cudaMalloc ((void **) &d_STF, nx_block * ny_block * nt * sizeof(double));
        if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMalloc d_STF: %s\n", cudaGetErrorString (err)); }
      }

      static double *d_SHF_QSW = (double *) NULL;
      if (d_SHF_QSW == NULL) {
        err = cudaMalloc ((void **) &d_SHF_QSW, nx_block * ny_block * sizeof(double));
        if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMalloc d_SHF_QSW: %s\n", cudaGetErrorString (err)); }
      }

      static double *d_SMF = (double *) NULL;
      if (d_SMF == NULL) {
        err = cudaMalloc ((void **) &d_SMF, nx_block * ny_block * 2 * sizeof(double));
        if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMalloc d_SMF: %s\n", cudaGetErrorString (err)); }
      }

      double convect_diff = *pconvect_diff;
      double convect_visc = *pconvect_visc;
      int fbid = *pbid;
      int bid = *pbid-1;//adjusted to c
      //! INPUT/OUTPUT PARAMETERS:
      static double *d_VDC = (double *) NULL;
      if (d_VDC == NULL) {
        err = cudaMalloc ((void **) &d_VDC, nx_block * ny_block * (km+1) * 2 * sizeof(double));
        if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMalloc d_VDC: %s\n", cudaGetErrorString (err)); }
      }

      //! OUTPUT PARAMETERS:
      static double *d_VVC = (double *) NULL;
      if (d_VVC == NULL) {
        err = cudaMalloc ((void **) &d_VVC, nx_block * ny_block * km * sizeof(double));
        if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMalloc d_VVC: %s\n", cudaGetErrorString (err)); }
      }

      static double *d_HMXL = (double *) NULL;
      if (d_HMXL == NULL) {
        err = cudaMalloc ((void **) &d_HMXL, nx_block * ny_block * sizeof(double));
        if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMalloc d_HMXL: %s\n", cudaGetErrorString (err)); }
        err = cudaMemset (d_HMXL, 0, nx_block * ny_block * sizeof(double));
        if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMemset d_HMXL: %s\n", cudaGetErrorString (err)); }
      }

      static double *d_KPP_HBLT = (double *) NULL;
      if (d_KPP_HBLT == NULL) {
        err = cudaMalloc ((void **) &d_KPP_HBLT, nx_block * ny_block * sizeof(double));
        if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMalloc d_KPP_HBLT: %s\n", cudaGetErrorString (err)); }
      }

      static double *d_KPP_SRC = (double *) NULL;
      if (d_KPP_SRC == NULL) {
        err = cudaMalloc ((void **) &d_KPP_SRC, nx_block * ny_block * km * nt * max_blocks_clinic * sizeof(double));
        if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMalloc d_KPP_SRC: %s\n", cudaGetErrorString (err)); }
        err = cudaMemset (d_KPP_SRC, 0, nx_block * ny_block * km * nt * max_blocks_clinic * sizeof(double));
        if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMemset d_KPP_SRC: %s\n", cudaGetErrorString (err)); }
        
      }

      //! EOP
      //! BOC
      //! -----------------------------------------------------------------------
      //!
      //! local variables
      //!
      //! -----------------------------------------------------------------------
      static int *d_KBL = (int *) NULL;
      if (d_KBL == NULL) {
        err = cudaMalloc ((void **) &d_KBL, nx_block * ny_block * sizeof(int));
        if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMalloc d_KBL: %s\n", cudaGetErrorString (err)); }
      }

      static double *d_USTAR = (double *) NULL;
      if (d_USTAR == NULL) {
        err = cudaMalloc ((void **) &d_USTAR, nx_block * ny_block * sizeof(double));
        if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMalloc d_USTAR: %s\n", cudaGetErrorString (err)); }
      }

      static double *d_BFSFC = (double *) NULL;
      if (d_BFSFC == NULL) {
        err = cudaMalloc ((void **) &d_BFSFC, nx_block * ny_block * sizeof(double));
        if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMalloc d_BFSFC: %s\n", cudaGetErrorString (err)); }
      }

      static double *d_STABLE = (double *) NULL;
      if (d_STABLE == NULL) {
        err = cudaMalloc ((void **) &d_STABLE, nx_block * ny_block * sizeof(double));
        if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMalloc d_STABLE: %s\n", cudaGetErrorString (err)); }
      }

      static double *d_DBLOC = (double *) NULL;
      if (d_DBLOC == NULL) {
        err = cudaMalloc ((void **) &d_DBLOC, nx_block * ny_block * km * sizeof(double));
        if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMalloc d_DBLOC: %s\n", cudaGetErrorString (err)); }
      }

      static double *d_DBSFC = (double *) NULL;
      if (d_DBSFC == NULL) {
        err = cudaMalloc ((void **) &d_DBSFC, nx_block * ny_block * km * sizeof(double));
        if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMalloc d_DBSFC: %s\n", cudaGetErrorString (err)); }
      }

      static double *d_GHAT = (double *) NULL;
      if (d_GHAT == NULL) {
        err = cudaMalloc ((void **) &d_GHAT, nx_block * ny_block * km * sizeof(double));
        if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMalloc d_GHAT: %s\n", cudaGetErrorString (err)); }
      }

      static double *d_VISC = (double *) NULL;
      if (d_VISC == NULL) {
        err = cudaMalloc ((void **) &d_VISC, nx_block * ny_block * (km+1) * sizeof(double));
        if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMalloc d_VISC: %s\n", cudaGetErrorString (err)); }
      }

      //! dummy variable that is not used by ri_iwmix but added here to have the
      //! same interface during testing
      static double *d_RHOMIX = (double *) NULL;
//      if (d_RHOMIX == NULL) {
//        err = cudaMalloc ((void **) &d_RHOMIX, nx_block * ny_block * km * sizeof(double));
//        if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMalloc d_RHOMIX: %s\n", cudaGetErrorString (err)); }
//      }

      //! temporary variable for storing velocity shear squared on u grid, converted to t grid by ri_iwmix
      static double *d_VSHEARU = (double *) NULL;
      if (d_VSHEARU == NULL) {
        err = cudaMalloc ((void **) &d_VSHEARU, nx_block * ny_block * km * sizeof(double));
        if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMalloc d_VSHEARU: %s\n", cudaGetErrorString (err)); }
      }
      static double *d_WORK1 = (double *) NULL;
      if (d_WORK1 == NULL) {
        err = cudaMalloc ((void **) &d_WORK1, nx_block * ny_block * sizeof(double));
        if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMalloc d_WORK1: %s\n", cudaGetErrorString (err)); }
      }
      static double *d_WORK2 = (double *) NULL;
      if (d_WORK2 == NULL) {
        err = cudaMalloc ((void **) &d_WORK2, nx_block * ny_block * sizeof(double));
        if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMalloc d_WORK2: %s\n", cudaGetErrorString (err)); }
      }
      static double *d_WORK3 = (double *) NULL;
      if (d_WORK3 == NULL) {
        err = cudaMalloc ((void **) &d_WORK3, nx_block * ny_block * sizeof(double));
        if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMalloc d_WORK3: %s\n", cudaGetErrorString (err)); }
      }

      //copy inputs from Host to Device
      err = cudaMemcpyAsync(d_TRCR, TRCR, nx_block * ny_block * km * nt * sizeof(double), cudaMemcpyHostToDevice, stream[0]);
      if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMemcpy host to device d_TRCR: %s\n", cudaGetErrorString(err)); }
      err = cudaMemcpyAsync(d_UUU, UUU, nx_block * ny_block * km * sizeof(double), cudaMemcpyHostToDevice, stream[0]);
      if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMemcpy host to device d_UUU: %s\n", cudaGetErrorString(err)); }
      err = cudaMemcpyAsync(d_VVV, VVV, nx_block * ny_block * km * sizeof(double), cudaMemcpyHostToDevice, stream[0]);
      if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMemcpy host to device d_VVV: %s\n", cudaGetErrorString(err)); }
      err = cudaMemcpyAsync(d_STF, STF, nx_block * ny_block * nt * sizeof(double), cudaMemcpyHostToDevice, stream[0]);
      if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMemcpy host to device d_STF: %s\n", cudaGetErrorString(err)); }
      err = cudaMemcpyAsync(d_SHF_QSW, SHF_QSW, nx_block * ny_block * sizeof(double), cudaMemcpyHostToDevice, stream[0]);
      if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMemcpy host to device d_SHF_QSW: %s\n", cudaGetErrorString(err)); }
      err = cudaMemcpyAsync(d_SMF, SMF, nx_block * ny_block * 2 * sizeof(double), cudaMemcpyHostToDevice, stream[0]);
      if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMemcpy host to device d_SMF: %s\n", cudaGetErrorString(err)); }
      err = cudaMemcpyAsync(d_VDC, VDC, nx_block * ny_block * (km+1) * 2 * sizeof(double), cudaMemcpyHostToDevice, stream[0]);
      if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMemcpy host to device d_VDC: %s\n", cudaGetErrorString(err)); }

      float time;

  //allocate host memory for reference 
  static double *VDC_REF = (double *) NULL;
  static double *VVC_REF = (double *) NULL;
  static double *TRCR_REF = (double *) NULL;
  static double *UUU_REF = (double *) NULL;
  static double *VVV_REF = (double *) NULL;
  static double *STF_REF = (double *) NULL;
  static double *SHF_QSW_REF = (double *) NULL;
  static double *SMF_REF = (double *) NULL;
  static double *HMXL_REF = (double *) NULL;
  static double *KPP_HBLT_REF = (double *) NULL;
  static double *KPP_SRC_REF = (double *) NULL;

if (TRCR_REF == NULL) {
  err = cudaHostAlloc ((void **) &TRCR_REF, nx_block * ny_block * km * nt * sizeof(double), cudaHostAllocMapped);
  if (err != cudaSuccess) { fprintf (stderr, "Error in cudaHostAlloc TRCR: %s\n", cudaGetErrorString (err)); }
  err = cudaHostAlloc ((void **) &UUU_REF, nx_block * ny_block * km * sizeof(double), cudaHostAllocMapped);
  if (err != cudaSuccess) { fprintf (stderr, "Error in cudaHostAlloc UUU: %s\n", cudaGetErrorString (err)); }
  err = cudaHostAlloc ((void **) &VVV_REF, nx_block * ny_block * km * sizeof(double), cudaHostAllocMapped);
  if (err != cudaSuccess) { fprintf (stderr, "Error in cudaHostAlloc VVV: %s\n", cudaGetErrorString (err)); }
  err = cudaHostAlloc ((void **) &STF_REF, nx_block * ny_block * nt * sizeof(double), cudaHostAllocMapped);
  if (err != cudaSuccess) { fprintf (stderr, "Error in cudaHostAlloc STF: %s\n", cudaGetErrorString (err)); }
  err = cudaHostAlloc ((void **) &SHF_QSW_REF, nx_block * ny_block * sizeof(double), cudaHostAllocMapped);
  if (err != cudaSuccess) { fprintf (stderr, "Error in cudaHostAlloc SHF_QSW: %s\n", cudaGetErrorString (err)); }
  err = cudaHostAlloc ((void **) &SMF_REF, nx_block * ny_block * 2 * sizeof(double), cudaHostAllocMapped);
  if (err != cudaSuccess) { fprintf (stderr, "Error in cudaHostAlloc SMF: %s\n", cudaGetErrorString (err)); }
  err = cudaHostAlloc ((void **) &VDC_REF, nx_block * ny_block * (km+1) * 2 * sizeof(double), cudaHostAllocMapped);
  if (err != cudaSuccess) { fprintf (stderr, "Error in cudaHostAlloc VDC: %s\n", cudaGetErrorString (err)); }
  err = cudaHostAlloc ((void **) &VVC_REF, nx_block * ny_block * km * sizeof(double), cudaHostAllocMapped);
  if (err != cudaSuccess) { fprintf (stderr, "Error in cudaHostAlloc VVC: %s\n", cudaGetErrorString (err)); }
  err = cudaHostAlloc ((void **) &HMXL_REF, nx_block * ny_block * sizeof(double), cudaHostAllocMapped);
  if (err != cudaSuccess) { fprintf (stderr, "Error in cudaHostAlloc HMXL: %s\n", cudaGetErrorString (err)); }
  err = cudaHostAlloc ((void **) &KPP_HBLT_REF, nx_block * ny_block * sizeof(double), cudaHostAllocMapped);
  if (err != cudaSuccess) { fprintf (stderr, "Error in cudaHostAlloc KPP_HBLT: %s\n", cudaGetErrorString (err)); }
  err = cudaHostAlloc ((void **) &KPP_SRC_REF, nx_block * ny_block * km * nt * max_blocks_clinic * sizeof(double), cudaHostAllocMapped);
  if (err != cudaSuccess) { fprintf (stderr, "Error in cudaHostAlloc KPP_SRC: %s\n", cudaGetErrorString (err)); }
}
  //make sure CPU and GPU version compute on same data
  err = cudaMemcpy(TRCR_REF, TRCR, nx_block * ny_block * km * nt * sizeof(double), cudaMemcpyHostToHost);
  if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMemcpyHostToHost %s\n", cudaGetErrorString (err)); }
  err = cudaMemcpy(UUU_REF, UUU, nx_block * ny_block * km * sizeof(double), cudaMemcpyHostToHost);
  if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMemcpyHostToHost %s\n", cudaGetErrorString (err)); }
  err = cudaMemcpy(VVV_REF, VVV, nx_block * ny_block * km * sizeof(double), cudaMemcpyHostToHost);
  if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMemcpyHostToHost %s\n", cudaGetErrorString (err)); }
  err = cudaMemcpy(STF_REF, STF, nx_block * ny_block * nt * sizeof(double), cudaMemcpyHostToHost);
  if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMemcpyHostToHost %s\n", cudaGetErrorString (err)); }
  err = cudaMemcpy(SHF_QSW_REF, SHF_QSW, nx_block * ny_block * sizeof(double), cudaMemcpyHostToHost);
  if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMemcpyHostToHost %s\n", cudaGetErrorString (err)); }
  err = cudaMemcpy(SMF_REF, SMF, nx_block * ny_block * 2 * sizeof(double), cudaMemcpyHostToHost);
  if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMemcpyHostToHost %s\n", cudaGetErrorString (err)); }
  err = cudaMemcpy(VDC_REF, VDC, nx_block * ny_block * (km+1) * 2 * sizeof(double), cudaMemcpyHostToHost);
  if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMemcpyHostToHost %s\n", cudaGetErrorString (err)); }
  err = cudaMemcpy(VVC_REF, VVC, nx_block * ny_block * km * sizeof(double), cudaMemcpyHostToHost);
  if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMemcpyHostToHost %s\n", cudaGetErrorString (err)); }
  err = cudaMemcpy(HMXL_REF, HMXL, nx_block * ny_block * sizeof(double), cudaMemcpyHostToHost);
  if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMemcpyHostToHost %s\n", cudaGetErrorString (err)); }
  err = cudaMemcpy(KPP_HBLT_REF, KPP_HBLT, nx_block * ny_block * sizeof(double), cudaMemcpyHostToHost);
  if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMemcpyHostToHost %s\n", cudaGetErrorString (err)); }
  err = cudaMemcpy(KPP_SRC_REF, KPP_SRC, nx_block * ny_block * km * nt * max_blocks_clinic * sizeof(double), cudaMemcpyHostToHost);
  if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMemcpyHostToHost %s\n", cudaGetErrorString (err)); }

  static int *KBL_REF = (int *) NULL;
  static double *USTAR_REF = (double *) NULL;
  static double *BFSFC_REF = (double *) NULL;
  static double *STABLE_REF = (double *) NULL;
  static double *DBLOC_REF = (double *) NULL;
  static double *DBSFC_REF = (double *) NULL;
  static double *GHAT_REF = (double *) NULL;
  static double *VISC_REF = (double *) NULL;
  static double *RHOMIX_REF = (double *) NULL;
  static double *VSHEARU_REF = (double *) NULL;

if (KBL_REF == NULL) {
  cudaHostAlloc ((void **) &KBL_REF, nx_block * ny_block * sizeof(int), cudaHostAllocMapped);
  cudaHostAlloc ((void **) &USTAR_REF, nx_block * ny_block * sizeof(double), cudaHostAllocMapped);
  cudaHostAlloc ((void **) &BFSFC_REF, nx_block * ny_block * sizeof(double), cudaHostAllocMapped);
  cudaHostAlloc ((void **) &STABLE_REF, nx_block * ny_block * sizeof(double), cudaHostAllocMapped);
  cudaHostAlloc ((void **) &DBLOC_REF, nx_block * ny_block * km * sizeof(double), cudaHostAllocMapped);
  cudaHostAlloc ((void **) &DBSFC_REF, nx_block * ny_block * km * sizeof(double), cudaHostAllocMapped);
  cudaHostAlloc ((void **) &GHAT_REF, nx_block * ny_block * km * sizeof(double), cudaHostAllocMapped);
  cudaHostAlloc ((void **) &VISC_REF, nx_block * ny_block * (km+1) * sizeof(double), cudaHostAllocMapped);
  cudaHostAlloc ((void **) &RHOMIX_REF, nx_block * ny_block * km * sizeof(double), cudaHostAllocMapped);
  cudaHostAlloc ((void **) &VSHEARU_REF, nx_block * ny_block * km * sizeof(double), cudaHostAllocMapped);
}

  static int *KBL = (int *) NULL;
  static double *USTAR = (double *) NULL;
  static double *BFSFC = (double *) NULL;
  static double *STABLE = (double *) NULL;
  static double *DBLOC = (double *) NULL;
  static double *DBSFC = (double *) NULL;
  static double *GHAT = (double *) NULL;
  static double *VISC = (double *) NULL;
  static double *RHOMIX = (double *) NULL;
  static double *VSHEARU = (double *) NULL;

if (KBL == NULL) {
  cudaHostAlloc ((void **) &KBL, nx_block * ny_block * sizeof(int), cudaHostAllocMapped);
  cudaHostAlloc ((void **) &USTAR, nx_block * ny_block * sizeof(double), cudaHostAllocMapped);
  cudaHostAlloc ((void **) &BFSFC, nx_block * ny_block * sizeof(double), cudaHostAllocMapped);
  cudaHostAlloc ((void **) &STABLE, nx_block * ny_block * sizeof(double), cudaHostAllocMapped);
  cudaHostAlloc ((void **) &DBLOC, nx_block * ny_block * km * sizeof(double), cudaHostAllocMapped);
  cudaHostAlloc ((void **) &DBSFC, nx_block * ny_block * km * sizeof(double), cudaHostAllocMapped);
  cudaHostAlloc ((void **) &GHAT, nx_block * ny_block * km * sizeof(double), cudaHostAllocMapped);
  cudaHostAlloc ((void **) &VISC, nx_block * ny_block * (km+1) * sizeof(double), cudaHostAllocMapped);
  cudaHostAlloc ((void **) &RHOMIX, nx_block * ny_block * km * sizeof(double), cudaHostAllocMapped);
  cudaHostAlloc ((void **) &VSHEARU, nx_block * ny_block * km * sizeof(double), cudaHostAllocMapped);
}

  //output parameter for debug
//  double *DEBUG;
//  double *DEBUG_REF;
//  cudaHostAlloc ((void **) &DEBUG, nx_block * ny_block * km * sizeof(double), cudaHostAllocMapped);
//  cudaHostAlloc ((void **) &DEBUG_REF, nx_block * ny_block * km * sizeof(double), cudaHostAllocMapped);

  cudaDeviceSynchronize();
  err = cudaGetLastError();
  if (err != cudaSuccess) { fprintf (stderr, "Error in cudaHostAllocs: %s\n", cudaGetErrorString (err)); }


      //setup threads and grid
      dim3 threads(BLOCK_X, BLOCK_Y);
      dim3 grid( (int)ceilf((float) nx_block / (float)(BLOCK_X)) , (int)ceilf((float) ny_block / (float)(BLOCK_Y)));

      //! -----------------------------------------------------------------------
      //!
      //! compute buoyancy differences at each vertical level.
      //!
      //! -----------------------------------------------------------------------
      cudaDeviceSynchronize();
      start_timer();
      buoydiff_kernel<<<grid, threads, 0, stream[0]>>>(d_DBLOC, d_DBSFC, d_TRCR, bid, KMT);
      cudaDeviceSynchronize();
      stop_timer(&time);
      printf("buoydiff: %.6f ms\n", time);

      err = cudaPeekAtLastError();
      if (err != cudaSuccess) { fprintf (stderr, "Error in buoydiff kernel: %s\n", cudaGetErrorString (err)); }

      cudaDeviceSynchronize();
      start_timer();
      buoydiff(DBLOC_REF, DBSFC_REF, TRCR_REF, &fbid);
      cudaDeviceSynchronize();
      stop_timer(&time);
      printf("buoydiff cpu: %.6f ms\n", time);
 
      err = cudaMemcpy(DBLOC, d_DBLOC, nx_block * ny_block * km * sizeof(double), cudaMemcpyDeviceToHost);
      if (err != cudaSuccess) { fprintf (stderr, "Error in buoydiff kernel memcpy: %s\n", cudaGetErrorString (err)); }
      err = cudaMemcpy(DBSFC, d_DBSFC, nx_block * ny_block * km * sizeof(double), cudaMemcpyDeviceToHost);
      if (err != cudaSuccess) { fprintf (stderr, "Error in buoydiff kernel memcpy: %s\n", cudaGetErrorString (err)); }
      compare(DBLOC, DBLOC_REF, nx_block * ny_block * km, "buoydiff DBLOC");
      compare(DBSFC, DBSFC_REF, nx_block * ny_block * km, "buoydiff DBSFC");


      //! -----------------------------------------------------------------------
      //!
      //! compute mixing due to shear instability, internal waves and
      //! convection
      //!
      //! -----------------------------------------------------------------------
      start_timer();
      compute_vshear_ugrid<<<grid, threads, 0, stream[0]>>>(d_VSHEARU, d_UUU, d_VVV, DZU, bid);
      cudaDeviceSynchronize();
      stop_timer(&time);
      printf("vshear ugrid: %.6f ms\n", time);

      err = cudaPeekAtLastError();
      if (err != cudaSuccess) { fprintf (stderr, "Error in vshear ugrid kernel: %s\n", cudaGetErrorString (err)); }
 
      start_timer();
      ri_iwmix_kernel<<<grid, threads, 0, stream[0]>>>(d_DBLOC, d_VISC, d_VDC, d_UUU, d_VVV, d_RHOMIX, convect_diff, convect_visc, bid, d_VSHEARU, DZT, KMT);
      cudaDeviceSynchronize();
      stop_timer(&time);
      printf("ri_iwmix: %.6f ms\n", time);

      err = cudaPeekAtLastError();
      if (err != cudaSuccess) { fprintf (stderr, "Error in ri_iwmix kernel: %s\n", cudaGetErrorString (err)); }

      start_timer();
      ri_iwmix(DBLOC_REF, VISC_REF, VDC_REF, UUU_REF, VVV_REF, RHOMIX_REF, &convect_diff, &convect_visc, &fbid);
      cudaDeviceSynchronize();
      stop_timer(&time);
      printf("ri_iwmix cpu: %.6f ms\n", time);

      err = cudaMemcpy(VDC, d_VDC, nx_block * ny_block * (km+1) * 2 * sizeof(double), cudaMemcpyDeviceToHost);
      if (err != cudaSuccess) { fprintf (stderr, "Error in ri_iwmix kernel memcpy: %s\n", cudaGetErrorString (err)); }
      err = cudaMemcpy(VISC, d_VISC, nx_block * ny_block * (km+1) * sizeof(double), cudaMemcpyDeviceToHost);
      if (err != cudaSuccess) { fprintf (stderr, "Error in ri_iwmix kernel memcpy: %s\n", cudaGetErrorString (err)); }
      compare(VDC, VDC_REF, nx_block * ny_block * (km+1) * 2, "ri_iwmix VDC");
      compare(VISC, VISC_REF, nx_block * ny_block * (km+1), "ri_iwmix VISC");
 
      //! -----------------------------------------------------------------------
      //!
      //! compute double diffusion if desired
      //!
      //! -----------------------------------------------------------------------
      start_timer();
      ddmix_kernel<<<grid, threads, 0, stream[0]>>>(d_VDC, d_TRCR, bid);
      cudaDeviceSynchronize();
      stop_timer(&time);
      printf("ddmix: %.6f ms\n", time);
   
      err = cudaPeekAtLastError();
      if (err != cudaSuccess) { fprintf (stderr, "Error in ddmix kernel: %s\n", cudaGetErrorString (err)); }

      start_timer();
      ddmix(VDC_REF, TRCR_REF, &fbid);
      cudaDeviceSynchronize();
      stop_timer(&time);
      printf("ddmix cpu: %.6f ms\n", time);
 
      err = cudaMemcpy(VDC, d_VDC, nx_block * ny_block * (km+1) * 2 * sizeof(double), cudaMemcpyDeviceToHost);
      if (err != cudaSuccess) { fprintf (stderr, "Error in ddmix kernel memcpy: %s\n", cudaGetErrorString (err)); }
      compare(VDC, VDC_REF, nx_block * ny_block * (km+1) * 2, "ddmix VDC");

      //! -----------------------------------------------------------------------
      //!
      //! compute boundary layer depth
      //!
      //! -----------------------------------------------------------------------

      start_timer();
      compute_ustaru<<<grid, threads, 0, stream[0]>>>(d_STABLE, d_SMF); //storing temporary result in stable
      cudaDeviceSynchronize();
      stop_timer(&time);
      printf("compute ustaru: %.6f ms\n", time);

      start_timer();
      ugrid_to_tgrid_kernel<<<grid, threads, 0, stream[0]>>>(d_USTAR, d_STABLE, 0);
      cudaDeviceSynchronize();
      stop_timer(&time);
      printf("ugrid_to_grid kernel: %.6f ms\n", time);

      start_timer();
      compute_vshearu_bldepth<<<grid, threads, 0, stream[0]>>>(d_VSHEARU, d_UUU, d_VVV, DZU, bid);
      cudaDeviceSynchronize();
      stop_timer(&time);
      printf("vshearu bldepth: %.6f ms\n", time);

      err = cudaPeekAtLastError();
      if (err != cudaSuccess) { fprintf (stderr, "Error in vshearu bldepth kernel: %s\n", cudaGetErrorString (err)); }
 
      start_timer();
      bldepth_kernel<<<grid, threads, 0, stream[0]>>>(d_DBLOC, d_DBSFC, d_TRCR, d_UUU, d_VVV, d_STF, d_SHF_QSW, d_KPP_HBLT, d_USTAR, d_BFSFC, d_STABLE, d_KBL, bid, d_SMF, d_VSHEARU, d_WORK1, d_WORK2, d_WORK3, DZT, DZU, KMT);
      cudaDeviceSynchronize();
      stop_timer(&time);
      printf("bldepth: %.6f ms\n", time);

      err = cudaPeekAtLastError();
      if (err != cudaSuccess) { fprintf (stderr, "Error in bldepth kernel: %s\n", cudaGetErrorString (err)); }
 
      start_timer();
      smooth_hblt_kernel<<<grid, threads, 0, stream[0]>>>(bid, d_KPP_HBLT, d_KBL, d_WORK3, KMT, DZT);
      cudaDeviceSynchronize();
      stop_timer(&time);
      printf("smooth hblt: %.6f ms\n", time);

      err = cudaPeekAtLastError();
      if (err != cudaSuccess) { fprintf (stderr, "Error in smooth hblt kernel: %s\n", cudaGetErrorString (err)); }
 
      start_timer();
      correct_stability_and_buoyancy<<<grid, threads, 0, stream[0]>>>(d_BFSFC, d_STABLE, d_WORK1, d_WORK2, d_KPP_HBLT); 
      cudaDeviceSynchronize();
      stop_timer(&time);
      printf("correct stable bfsfc: %.6f ms\n", time);
    
      err = cudaPeekAtLastError();
      if (err != cudaSuccess) { fprintf (stderr, "Error in correct stable bfsfc kernel: %s\n", cudaGetErrorString (err)); }

      start_timer();
      bldepth_test(DBLOC_REF, DBSFC_REF, TRCR_REF, UUU_REF, VVV_REF, STF_REF, SHF_QSW_REF, KPP_HBLT_REF, USTAR_REF, BFSFC_REF, STABLE_REF, KBL_REF, &fbid, SMF_REF);
      cudaDeviceSynchronize();
      stop_timer(&time);
      printf("bldepth cpu: %.6f ms\n", time);

      err = cudaMemcpy(BFSFC, d_BFSFC, nx_block * ny_block * sizeof(double), cudaMemcpyDeviceToHost);
      if (err != cudaSuccess) { fprintf (stderr, "Error in bldepth kernel memcpy: %s\n", cudaGetErrorString (err)); }
      err = cudaMemcpy(KPP_HBLT, d_KPP_HBLT, nx_block * ny_block * sizeof(double), cudaMemcpyDeviceToHost);
      if (err != cudaSuccess) { fprintf (stderr, "Error in bldepth kernel memcpy: %s\n", cudaGetErrorString (err)); }
      err = cudaMemcpy(STABLE, d_STABLE, nx_block * ny_block * sizeof(double), cudaMemcpyDeviceToHost);
      if (err != cudaSuccess) { fprintf (stderr, "Error in bldepth kernel memcpy: %s\n", cudaGetErrorString (err)); }
      err = cudaMemcpy(USTAR, d_USTAR, nx_block * ny_block * sizeof(double), cudaMemcpyDeviceToHost);
      if (err != cudaSuccess) { fprintf (stderr, "Error in bldepth kernel memcpy: %s\n", cudaGetErrorString (err)); }
      err = cudaMemcpy(KBL, d_KBL, nx_block * ny_block * sizeof(int), cudaMemcpyDeviceToHost);
      if (err != cudaSuccess) { fprintf (stderr, "Error in bldepth kernel memcpy: %s\n", cudaGetErrorString (err)); }

      compare(BFSFC, BFSFC_REF, nx_block * ny_block, "bldepth BFSFC");
      compare(KPP_HBLT, KPP_HBLT_REF, nx_block * ny_block, "bldepth HBLT");
      compare(STABLE, STABLE_REF, nx_block * ny_block, "bldepth STABLE");
      compare(USTAR, USTAR_REF, nx_block * ny_block, "bldepth USTAR");
      compare_int(KBL, KBL_REF, nx_block * ny_block);
 
      //! -----------------------------------------------------------------------
      //!
      //! compute boundary layer diffusivities
      //!
      //! -----------------------------------------------------------------------
      start_timer();
      blmix_kernel<<<grid, threads, 0, stream[0]>>>(d_VISC, d_VDC, d_KPP_HBLT, d_USTAR, d_BFSFC, d_STABLE, d_KBL, d_GHAT, bid, DZT);
      cudaDeviceSynchronize();
      stop_timer(&time);
      printf("blmix: %.6f ms\n", time);

      err = cudaPeekAtLastError();
      if (err != cudaSuccess) { fprintf (stderr, "Error in blmix kernel: %s\n", cudaGetErrorString (err)); }
 
      start_timer();
      blmix(VISC_REF, VDC_REF, KPP_HBLT_REF, USTAR_REF, BFSFC_REF, STABLE_REF, KBL_REF, GHAT_REF, &fbid);
      cudaDeviceSynchronize();
      stop_timer(&time);
      printf("blmix cpu: %.6f ms\n", time);

      err = cudaMemcpy(VDC, d_VDC, nx_block * ny_block * (km+1) * 2 * sizeof(double), cudaMemcpyDeviceToHost);
      if (err != cudaSuccess) { fprintf (stderr, "Error in blmix kernel memcpy: %s\n", cudaGetErrorString (err)); }
      err = cudaMemcpy(VISC, d_VISC, nx_block * ny_block * (km+1) * sizeof(double), cudaMemcpyDeviceToHost);
      if (err != cudaSuccess) { fprintf (stderr, "Error in blmix kernel memcpy: %s\n", cudaGetErrorString (err)); }
      err = cudaMemcpy(GHAT, d_GHAT, nx_block * ny_block * km * sizeof(double), cudaMemcpyDeviceToHost);
      if (err != cudaSuccess) { fprintf (stderr, "Error in blmix kernel memcpy: %s\n", cudaGetErrorString (err)); }

      compare(VDC, VDC_REF, nx_block * ny_block * (km+1) * 2, "blmix VDC");
      compare(VISC, VISC_REF, nx_block * ny_block * (km+1), "blmix VISC");
      compare(GHAT, GHAT_REF, nx_block * ny_block * km, "blmix GHAT");

      //! -----------------------------------------------------------------------
      //!
      //! consider interior convection:
      //!
      //! -----------------------------------------------------------------------

      cudaDeviceSynchronize();
      start_timer();
      interior_convection_kernel<<<grid, threads, 0, stream[0]>>>(d_DBLOC, d_DBSFC, d_KBL, d_STF, d_GHAT, d_VISC, convect_diff, convect_visc, bid, d_VDC, d_VVC, d_KPP_SRC, d_HMXL, d_VSHEARU, DZT, KMT); //reuse VSHEARU to store VISC on T grid
      cudaDeviceSynchronize();
      stop_timer(&time);
      printf("interior part1: %.6f ms\n", time);

      err = cudaPeekAtLastError();
      if (err != cudaSuccess) { fprintf (stderr, "Error in interior part1 kernel: %s\n", cudaGetErrorString (err)); }
 
      start_timer();
      interior_convection_part2<<<grid, threads, 0, stream[0]>>>(d_DBLOC, d_DBSFC, d_KBL, d_STF, d_GHAT, d_VISC, convect_diff, convect_visc, bid, d_VDC, d_VVC, d_KPP_SRC, d_HMXL, AU, UAREA_R, d_VSHEARU, DZT, KMT, KMU);
      cudaDeviceSynchronize();
      stop_timer(&time);
      printf("interior part2: %.6f ms\n", time);
    
      err = cudaPeekAtLastError();
      if (err != cudaSuccess) { fprintf (stderr, "Error in interior part2 kernel: %s\n", cudaGetErrorString (err)); }
 
      start_timer();
      interior_convection(DBLOC_REF, DBSFC_REF, KBL_REF, STF_REF, GHAT_REF, VISC_REF, &convect_diff, &convect_visc, &fbid, VDC_REF, VVC_REF, KPP_SRC_REF, HMXL_REF);
      cudaDeviceSynchronize();
      stop_timer(&time);
      printf("interior convection cpu: %.6f ms\n", time);

      err = cudaMemcpy(VVC, d_VVC, nx_block * ny_block * km * sizeof(double), cudaMemcpyDeviceToHost);
      if (err != cudaSuccess) { fprintf (stderr, "Error in blmix kernel memcpy: %s\n", cudaGetErrorString (err)); }
      err = cudaMemcpy(VDC, d_VDC, nx_block * ny_block * (km+1) * 2 * sizeof(double), cudaMemcpyDeviceToHost);
      if (err != cudaSuccess) { fprintf (stderr, "Error in blmix kernel memcpy: %s\n", cudaGetErrorString (err)); }
      err = cudaMemcpy(HMXL, d_HMXL, nx_block * ny_block * sizeof(double), cudaMemcpyDeviceToHost);
      if (err != cudaSuccess) { fprintf (stderr, "Error in blmix kernel memcpy: %s\n", cudaGetErrorString (err)); }
      err = cudaMemcpy(KPP_SRC, d_KPP_SRC, nx_block * ny_block * km * nt * max_blocks_clinic * sizeof(double), cudaMemcpyDeviceToHost);
      if (err != cudaSuccess) { fprintf (stderr, "Error in blmix kernel memcpy: %s\n", cudaGetErrorString (err)); }

      compare(VVC, VVC_REF, nx_block * ny_block * km, "interior VVC");
      compare(VDC, VDC_REF, nx_block * ny_block * (km+1) * 2, "interior VDC");
      compare(HMXL, HMXL_REF, nx_block * ny_block, "interior HMXL");
      compare(KPP_SRC, KPP_SRC_REF, nx_block * ny_block * km * nt * max_blocks_clinic, "interior KPP_SRC");


      //! -----------------------------------------------------------------------
      //! EOC
          //copy outputs from Device to Host
      err = cudaMemcpyAsync(VDC, d_VDC, nx_block * ny_block * (km+1) * 2 * sizeof(double), cudaMemcpyDeviceToHost, stream[0]);
      if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMemcpy device to host VDC: %s\n", cudaGetErrorString (err)); }
      err = cudaMemcpyAsync(VVC, d_VVC, nx_block * ny_block * km * sizeof(double), cudaMemcpyDeviceToHost, stream[0]);
      if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMemcpy device to host VVC: %s\n", cudaGetErrorString (err)); }
      err = cudaMemcpyAsync(HMXL, d_HMXL, nx_block * ny_block * sizeof(double), cudaMemcpyDeviceToHost, stream[0]);
      if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMemcpy device to host HMXL: %s\n", cudaGetErrorString (err)); }
      err = cudaMemcpyAsync(KPP_HBLT, d_KPP_HBLT, nx_block * ny_block * sizeof(double), cudaMemcpyDeviceToHost, stream[0]);
      if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMemcpy device to host KPP_HBLT: %s\n", cudaGetErrorString (err)); }
      err = cudaMemcpyAsync(KPP_SRC, d_KPP_SRC, nx_block * ny_block * km * nt * sizeof(double), cudaMemcpyDeviceToHost, stream[0]);
      if (err != cudaSuccess) { fprintf (stderr, "Error in cudaMemcpy device to host KPP_SRC: %s\n", cudaGetErrorString (err)); }

      cudaDeviceSynchronize();
      err = cudaGetLastError();
      if (err != cudaSuccess) { fprintf (stderr, "Error in test version of entrypoint: %s\n", cudaGetErrorString (err)); }

      printf("finished execution of test version of entrypoint\n");

    }
  








  
    __global__ void interior_convection_kernel (double *DBLOC, double *DBSFC, int *KBL, double *STF, double *GHAT, double *VISC, double convect_diff, double convect_visc, int bid, double *VDC, double *VVC, double *KPP_SRC, double *HMXL, double *VISCT, double *DZT, int *KMT) {

      int j = threadIdx.y + blockIdx.y * blockDim.y;
      int i = threadIdx.x + blockIdx.x * blockDim.x;

      if (j < ny_block && i < nx_block) {
  
      //! DESCRIPTION:
      //!
      //! This routine implements the final part of vmix_coeffs_kpp.
      //! It is rewritten into a separate routine for the purpose of creating a
      //! GPU implementation of the entire vertical mixing scheme.
      //! Previously USTAR and BFSFC were reused as temps by this final part,
      //! these were renamed to WORK3 and WORK4 respectively.
      //!
      //! PARAMETERS
      //! INPUTS:
      //! fake INOUTS:
      //! these are actually inputs but are listed as inouts because they are also
      //! written by part of the computation, we are however not interested in their
      //! final value after the execution of this subroutine
      //! true INOUTS:
      //! OUTPUTS:
      //! local variables:
      double WORK1;
      double FCON;
      int k;
    
      //! -----------------------------------------------------------------------
      //! compute function of Brunt-Vaisala squared for convection.
      //!
      //! use either a smooth
      //!
      //! WORK1 = N**2, FCON is function of N**2
      //! FCON = 0 for N**2 > 0
      //! FCON = [1- (1-WORK1/BVSQcon) **2]**3 for BVSQcon < N**2 < 0
      //! FCON = 1 for N**2 < BVSQcon
      //!
      //! or a step function. The smooth function has been used with
      //! BVSQcon = -0.2e-4_dbl_kind.
      //!
      //! after convection, average viscous coeffs to U-grid and reset sea
      //! floor values
      //!
      //! -----------------------------------------------------------------------
      for (k = 0; k<km-1; k++) {
        
        WORK1 = DBLOC(i,j,k) / (p5* (DZT(i,j,k,bid) + DZT(i,j,k+1,bid) ) );
      
        if (WORK1 > c0) { //converted from where statement
          FCON = c0;
        } else {
          FCON = c1;
        }
      
        //! *** add convection and reset sea floor values to zero
        //removed for (j=1;j<ny_block;j+=1)
        //removed for (i=1;i<nx_block;i+=1)

        if ( k >= KBL (i,j)-1 ) { //adjusted to c
          VISC(i,j,k) = VISC(i,j,k) + convect_visc * FCON;
          VDC(i,j,k,0) = VDC(i,j,k,0) + convect_diff * FCON; //adjusted to c
          VDC(i,j,k,1) = VDC(i,j,k,1) + convect_diff * FCON;//adjusted to c
        }
    
        if (k >= KMT (i,j,bid)-1 ) { //adjusted to c
          VISC(i,j,k) = c0;
          VDC(i,j,k,0) = c0;//adjusted to c
          VDC(i,j,k,1) = c0;//adjusted to c
        }
      
	VISCT(i,j,k) = VISC(i,j,k);
      }
    
      VISCT(i,j,km-1) = VISC(i,j,km-1);
      VISCT(i,j,km) = VISC(i,j,km);

      VDC(i,j,km-1,0) = c0; //adjusted to c
      VDC(i,j,km-1,1) = c0; //adjusted to c

      }
//split kernel here
    }

    __global__ void interior_convection_part2 (double *DBLOC, double *DBSFC, int *KBL, double *STF, double *GHAT, double *VISC, double convect_diff, double convect_visc, int bid, double *VDC, double *VVC, double *KPP_SRC, double *HMXL, double *AU, double *UAREA_R, double *VISCT, double *DZT, int *KMT, int *KMU) {

      int j = threadIdx.y + blockIdx.y * blockDim.y;
      int i = threadIdx.x + blockIdx.x * blockDim.x;
      int k;
      int n;
    
      if (j < ny_block && i < nx_block) {

      double WORK2;
      double WORK3;
      double WORK4;
      double WORK5;

      //! *** now average visc to U grid
      for (k = 0; k<km-1; k++) {
        WORK2 = tgrid_to_ugrid(VISCT, AU, UAREA_R, bid, i, j, k);
     
        VVC(i,j,k) = ((k < KMU(i,j,bid)-1 ) ? WORK2 : c0); //adjusted k to k+1 in c

      }
      VVC(i,j,km-1) = c0;//adjusted to c

      //! -----------------------------------------------------------------------
      //!
      //! add ghatp term from previous computation to right-hand-side
      //! source term on current row
      //!
      //! -----------------------------------------------------------------------
      for (n = 0; n<nt; n++) {
        
        KPP_SRC(i,j,0,n,bid) = STF(i,j,n) /dz(0) * (-VDC(i,j,0,n) *GHAT(i,j,0) ); //adjusted to c
      
        for (k = 1; k<km; k++) {
//          KPP_SRC(:,:,k,n,bid) = STF(:,:,n) /DZT(:,:,k,bid) * ( VDC(:,:,k-1,mt2) *GHAT(:,:,k-1) -VDC(:,:,k ,mt2) *GHAT (:,:,k) )
          KPP_SRC(i,j,k,n,bid) = STF(i,j,n) /DZT(i,j,k,bid) * ( VDC(i,j,k-1,n)   *GHAT(i,j,k-1) -VDC(i,j,k,n) *GHAT(i,j,k) );
        }
      }

      //! -----------------------------------------------------------------------
      //!
      //! compute diagnostic mixed layer depth (cm) using a max buoyancy
      //! gradient criterion. This part previously used USTAR and BFSFC as temps.
      //! These have been renamed to WORK3 and WORK4
      //!
      //! -----------------------------------------------------------------------
      WORK3 = c0;
      double hmxl = 0.0;
    
      if (KMT(i,j,bid) == 1) { //converted from where statement
        hmxl = zt(0); //adjusted to c
      } else {
        hmxl = c0;
      }
    
      for (k = 1; k<km; k++) {
        
        if (k <= KMT(i,j,bid)-1 ) { //converted from where statement  //adjusted to c
          WORK5 = zt(k-1) + p5* (DZT(i,j,k-1,bid) + DZT(i,j,k,bid) );
          WORK3 = max (DBSFC(i,j,k) /WORK5,WORK3);
          hmxl = WORK5;
        }
      }
    
      VISC(i,j,0) = c0;
    
      for (k = 1; k<km; k++) {
        
        if (WORK3 > c0 ) { //converted from where statement
          VISC(i,j,k) = (DBSFC(i,j,k) -DBSFC(i,j,k-1) ) / (p5* (DZT(i,j,k,bid) + DZT(i,j,k-1,bid) ) );
        }
      
        if ( (VISC(i,j,k) >= WORK3) && ((VISC(i,j,k)-VISC(i,j,k-1)) != c0) && (WORK3 > c0) ) { //converted from where statement
          WORK4 = (VISC(i,j,k) - WORK3) / (VISC(i,j,k)-VISC(i,j,k-1) );
          //! tqian
          //! HMXL (:,:,bid) = (zt (k-1) + p5*DZT (:,:,k-1,bid) ) * (c1-WORK4) ! + (zt (k-1) - p5*DZT (:,:,k-1,bid) ) *WORK4
          hmxl = (zt(k-1) + p25* (DZT(i,j,k-1,bid) +DZT(i,j,k,bid) ) ) * (c1-WORK4) + (zt(k-1) - p25* (DZT(i,j,k-2,bid) +DZT(i,j,k-1,bid) ) ) *WORK4;
          WORK3 = c0;
        }
      }

      HMXL(i,j,bid) = hmxl;


      }
    }

    /* this function computes the velocity shear squared on the u grid
     * it will be averaged to the t grid by ri_iwmix
     */
    __global__ void compute_vshear_ugrid(double *VSHEARU, double *UUU, double *VVV, double *DZU, int bid) {
      int j = threadIdx.y + blockIdx.y * blockDim.y;
      int i = threadIdx.x + blockIdx.x * blockDim.x;
      int k;
 
      if (j < ny_block && i < nx_block) {

      double tmp;
      double u;
      double v;
      double dz;

      for (k = 0; k < km-1; k++) {
        //! -----------------------------------------------------------------------
        //!
        //! compute velocity shear squared on U grid
        //! VSHEAR = (UUU (k) -UUU (k+1) ) **2+ (VVV (k) -VVV (k+1) ) **2
        //!
        //! -----------------------------------------------------------------------
        u = UUU(i,j,k)-UUU(i,j,k+1);
        v = VVV(i,j,k)-VVV(i,j,k+1);
        tmp = u*u + v*v;
        dz = p5* (DZU(i,j,k,bid) + DZU(i,j,k+1,bid));
        
        VSHEARU(i,j,k) = tmp / (dz*dz);
      }
      VSHEARU(i,j,km-1) = c0;
      
      }
    }

  
    //! ***********************************************************************
    //! BOP
    //! IROUTINE: ri_iwmix
    //! INTERFACE:
    __global__ void ri_iwmix_kernel (double *DBLOC, double *VISC, double *VDC, double *UUU, double *VVV, double *RHOMIX, double convect_diff, double convect_visc, int bid, double *VSHEARU, double *DZT, int *KMT) {

      int j = threadIdx.y + blockIdx.y * blockDim.y;
      int i = threadIdx.x + blockIdx.x * blockDim.x;

      if (j < ny_block && i < nx_block) {

      //! DESCRIPTION:
      //! Computes viscosity and diffusivity coefficients for the interior
      //! ocean due to shear instability (richardson number dependent) ,
      //! internal wave activity, and to static instability (Ri < 0) .
      //!
      //! REVISION HISTORY:
      //! Ben changed VDC from in out to output, because it was not read
      //! but only written
      //! INPUT PARAMETERS:
      //! OUTPUT PARAMETERS:
      //! EOP
      //! BOC
      //! -----------------------------------------------------------------------
      //!
      //! local variables
      //!
      //! -----------------------------------------------------------------------
      int k;
      int n;
      double VSHEAR;
      double RI_LOC;
      double FRI;
    
      //! -----------------------------------------------------------------------
      //!
      //! compute mixing at each level
      //!
      //! -----------------------------------------------------------------------
    
      //! VISC (:,:,0) = c0
      for (k  = 0; k < km; k ++) {

        //average velocity shear squared from u grid to t grid
        if (k < km-1) {//adjusted to c
          VSHEAR = ugrid_to_tgrid(VSHEARU,i,j,k);
        } else {
          VSHEAR = 0.0;
        }


        //! -----------------------------------------------------------------------
        //!
        //! compute local richardson number
        //! use visc array as temporary Ri storage to be smoothed
        //!
        //! -----------------------------------------------------------------------
        if (k < km-1) {//adjusted to c
          double pterm = p5* (DZT(i,j,k,bid) + DZT(i,j,k+1,bid));
          RI_LOC = DBLOC(i,j,k) / ( VSHEAR + eps / (pterm*pterm) )
                                / pterm ;
        }
        else {
          double pterm = p5*DZT(i,j,k,bid);
          RI_LOC = DBLOC(i,j,k) / (VSHEAR + eps / (pterm*pterm)) 
                                / pterm ;
        }

        if (k == 0) {//adjusted to c
          VISC(i,j,k) = (k+1 <= KMT(i,j,bid) ? RI_LOC : c0); //adjusted to c
        }
        else {
          VISC(i,j,k) = (k+1 <= KMT(i,j,bid) ? RI_LOC : VISC(i,j,k-1)); //adjusted to c
        }

      }
    


      //! -----------------------------------------------------------------------
      //!
      //! vertically smooth Ri num_v_smooth_Ri times with 1-2-1 weighting
      //! result again stored temporarily in VISC and use RI_LOC and FRI
      //! as temps
      //!
      //! -----------------------------------------------------------------------
      for (n = 0; n < num_v_smooth_Ri; n++) {
        
        FRI = p25 * VISC(i,j,0);  //adjusted to c
        VISC(i,j,km) = VISC(i,j,km-1);//adjusted to c
      
        for (k = 0; k<km; k++) {
          
          //! DIR$ NODEP
          //removed for (j=1;j<ny_block;j+=1)
          //! DIR$ NODEP
          //removed for (i=1;i<nx_block;i+=1)
          
          RI_LOC = VISC(i,j,k);
        
          if (KMT (i,j,bid) >= 3) {
            VISC(i,j,k) = FRI + p5*RI_LOC + p25*VISC(i,j,k+1);
          }
        
          FRI = p25*RI_LOC;
        }
      }



      //! -----------------------------------------------------------------------
      //!
      //! now that we have a smoothed Ri field, finish computing coeffs
      //! at each level
      //!
      //! -----------------------------------------------------------------------
      for (k = 0; k < km; k++) {
        
        //! -----------------------------------------------------------------------
        //!
        //! if Ri-number mixing requested,
        //! evaluate function of Ri for shear instability:
        //! for 0 < Ri < Riinfty, function = (1 - (Ri/Riinfty) **2) **3
        //! for Ri > Riinfty, function = 0
        //! for Ri < 0 , function = 1
        //! compute contribution due to shear instability
        //! VISC holds smoothed Ri at k, but replaced by real VISC
        //!
        //! otherwise only use iw
        //! convection is added later
        //!
        //! -----------------------------------------------------------------------
        //if ( k < km-1 ) {//adjusted to c
        //  KVMIX = bckgrnd_vdc(i,j,k,bid);
        //}
      
        FRI = min ( (max (VISC(i,j,k), c0)) / Riinfty, c1);

        double oneminusfri = (1.0 - (FRI*FRI));
        double oneminusfricubed = oneminusfri*oneminusfri*oneminusfri;
        VISC(i,j,k) = bckgrnd_vvc(i,j,k,bid) + rich_mix * oneminusfricubed;
      
        if ( k < km-1 ) {//adjusted to c
          VDC(i,j,k,1) = bckgrnd_vdc(i,j,k,bid) + rich_mix * oneminusfricubed; //adjusted to c
          VDC(i,j,k,0) = VDC(i,j,k,1); //adjusted to c
        }
      
        //! -----------------------------------------------------------------------
        //!
        //! set seafloor values to zero
        //!
        //! -----------------------------------------------------------------------
        //! DIR$ NODEP
        //! DIR$ COLLAPSE
        //removed for (j=1;j<ny_block;j+=1)
        //! DIR$ NODEP
        //removed for (i=1;i<nx_block;i+=1)
        
        if ( k+1 >= KMT (i,j,bid) ) {//adjusted to c
          VISC(i,j,k) = c0;
          VDC(i,j,k,0) = c0;//adjusted to c
          VDC(i,j,k,1) = c0;//adjusted to c
        }
      
        //! if (tavg_requested (tavg_KVMIX) ) then
        //! call accumulate_tavg_field (KVMIX,tavg_KVMIX,bid,k)
        //! end if
        //! if (tavg_requested (tavg_TPOWER) ) then
        //! WORK1 (:,:) = KVMIX (:,:) *RHOMIX (:,:,k) *DBLOC (:,:,k) / ! (zgrid (k) - zgrid (k+1) )
        //! call accumulate_tavg_field (WORK1,tavg_TPOWER,bid,k)
        //! end if
        //! -----------------------------------------------------------------------
        //!
        //! move to next level.
        //!
        //! -----------------------------------------------------------------------
      }


    
      //! -----------------------------------------------------------------------
      //!
      //! fill extra coefficients for blmix
      //!
      //! -----------------------------------------------------------------------
      //! VISC (:,:,0 ) = c0 no longer necessary Ben
      //! VDC (:,:,0,:) = c0 no longer necessary Ben
      VISC(i,j,km) = c0;//adjusted to c
      VDC(i,j,km-1,0) = c0;//adjusted to c
      VDC(i,j,km-1,1) = c0;//adjusted to c
      VDC(i,j,km,0) = c0;//adjusted to c
      VDC(i,j,km,1) = c0;//adjusted to c
  
  
      //! -----------------------------------------------------------------------
      //! EOC


      }
    }
  


   __global__ void compute_ustaru(double *USTAR, double *SMF) {

      int j = threadIdx.y + blockIdx.y * blockDim.y;
      int i = threadIdx.x + blockIdx.x * blockDim.x;

      if (j < ny_block && i < nx_block) {

      //! -----------------------------------------------------------------------
      //!
      //! compute friction velocity USTAR on U-grid
      //!
      //! -----------------------------------------------------------------------
      USTAR(i,j) = sqrt (sqrt (SMF(i,j,0)*SMF(i,j,0) + SMF(i,j,1)*SMF(i,j,1)) ); //adjusted to c
    
      }
    }



    /**
     * This method computes the velocity shear squared on the U grid, as needed by bldepth.
     *
     * This method computes the square of the velocity differences with the first level on the tracer grid and are taken
     * in bldepth as the maximum of the neighboring velocity points.
     *
     */
    __global__ void compute_vshearu_bldepth(double *VSHEARU, double *UUU, double *VVV, double *DZU, int bid) {

      int j = threadIdx.y + blockIdx.y * blockDim.y;
      int i = threadIdx.x + blockIdx.x * blockDim.x;
      int kl;

      if (j < ny_block && i < nx_block) {

      double dzusfc = DZU(i,j,0,bid);

      for (kl = 1; kl < km; kl++) {    
        double u = UUU(i,j,0)-UUU(i,j,kl);
        double v = VVV(i,j,0)-VVV(i,j,kl);

        double work = u*u + v*v;
        double denom = (-zgrid(kl-1) + p5* (DZU(i,j,kl,bid) + DZU(i,j,kl-1,bid) - dzusfc ) );
        VSHEARU(i,j,kl) = work / (denom*denom);

      }

      }
    }


    //! ***********************************************************************
    //! BOP
    //! IROUTINE: bldepth
    //! INTERFACE:
    __global__ void bldepth_kernel (double *DBLOC, double *DBSFC, double *TRCR, double *UUU, double *VVV, double *STF, double *SHF_QSW, double *HBLT, double *USTAR, double *BFSFC, double *STABLE, int *KBL, int bid, double *SMF, double *VSHEARU, double *d_WORK1, double *d_WORK2, double *d_WORK3, double *DZT, double *DZU, int *KMT) {

      int j = threadIdx.y + blockIdx.y * blockDim.y;
      int i = threadIdx.x + blockIdx.x * blockDim.x;
      
      if (j < ny_block && i < nx_block) {
      //! DESCRIPTION:
      //! This routine computes the ocean boundary layer depth defined as
      //! the shallowest depth where the bulk Richardson number is equal to
      //! a critical value, Ricr.
      //!
      //! NOTE: bulk richardson numbers are evaluated by computing
      //! differences between values at zgrid (kl) $< 0$ and surface
      //! reference values. currently, the reference values are equal
      //! to the values in the surface layer. when using higher
      //! vertical grid resolution, these reference values should be
      //! computed as the vertical averages from the surface down to
      //! epssfc*zgrid (kl) .
      //!
      //! This routine also computes where surface forcing is stable
      //! or unstable (STABLE)
      //!
      //! REVISION HISTORY:
      //! same as module
      //! INPUT PARAMETERS:
      //! *** either one or the other (not
      //! *** both) should be passed
      //! OUTPUT PARAMETERS:
      //! EOP
      //! BOC
      //! -----------------------------------------------------------------------
      //!
      //! local variables
      //!
      //! -----------------------------------------------------------------------
      int kupper;
      int kup;
      int kdn;
      int ktmp;
      int kl;
      double VSHEAR;
      double SIGMA;
      double WM;
      double WS;
      double BO;
      double BOSOL;
      double WORK;
      double ZKL;
      double B_FRQNCY;
      //! (= min (HEKMAN,HMONOB) )
      double RI_BULK[3];
      double absorb_frac;
      double sqrt_arg;
      double z_upper;
      double z_up;
      double a_co;
      double b_co;
      double c_co;
      //! $ (a_{co}z^2+b_{co}|z|+c_{co}=Ri_b) used to
      //! find the boundary layer depth. when
      //! finding the roots, c_co = c_co - Ricr
      double slope_up;
    
      //! at zup. this is used as a boundary
      //! condition to determine the coefficients.
 
      //! -----------------------------------------------------------------------
      //!
      //! average friction velocity (USTAR) computed on U-grid to T-grid.
      //!
      //! -----------------------------------------------------------------------
       //USTAR has been precomputed by other kernels


      //! -----------------------------------------------------------------------
      //!
      //! compute density and expansion coefficients at surface
      //!
      //! -----------------------------------------------------------------------
      WORK = (TRCR(i,j,0,0) < -c2 ? -c2 : TRCR(i,j,0,0)); //adjusted to c

      double tq = max(min(WORK,TMAX),TMIN);
      double sq = 1000.0 * max(min(TRCR(i,j,0,1),SMAX),SMIN);
      double sqr = sqrt(sq);

      double nomk = mwjf_numerator(tq, sq, 0);
      double denomk = mwjf_denominator(tq, sq, sqr, 0);
      double RHO1 = nomk * denomk;
      double TALPHA = compute_drhodt(tq, sq, sqr, 0, nomk, denomk);
      double SBETA = compute_drhods(tq, sq, sqr, 0, nomk, denomk);

      //! -----------------------------------------------------------------------
      //!
      //! compute turbulent and radiative sfc buoyancy forcing
      //!
      //! -----------------------------------------------------------------------
      //removed for (j=1;j<ny_block;j+=1)
      //removed for (i=1;i<nx_block;i+=1)

      if (RHO1 != c0) {
        BO = grav* (-TALPHA *STF(i,j,0) - SBETA *STF(i,j,1) ) /RHO1;   //adjusted to c
        BOSOL = -grav*TALPHA *SHF_QSW(i,j) /RHO1;
      }
      else {
        BO = c0;
        BOSOL = c0;
      }

      //! -----------------------------------------------------------------------
      //!
      //! Find bulk Richardson number at every grid level until > Ricr
      //! max values when Ricr never satisfied are KBL = KMT and
      //! HBLT = -zgrid (KMT)
      //!
      //! NOTE: the reference depth is -epssfc/2.*zgrid (i,k) , but the
      //! reference u,v,t,s values are simply the surface layer
      //! values and not the averaged values from 0 to 2*ref.depth,
      //! which is necessary for very fine grids (top layer < 2m
      //! thickness)
      //!
      //!
      //! Initialize hbl and kbl to bottomed out values
      //! Initialize HEKMAN and HLIMIT (= HMONOB until reset) to model bottom
      //! Initialize Monin Obukhov depth to value at z_up
      //! Set HMONOB=-zgrid (km) if unstable
      //!
      //! -----------------------------------------------------------------------
      kupper = 0; //adjusted to c
      kup = 1; //adjusted to c
      kdn = 2; //adjusted to c
      z_upper = c0;
      z_up = zgrid(0); //adjusted to c
      RI_BULK(kupper) = c0;
      RI_BULK(kup) = c0;
      KBL(i,j) = ((KMT(i,j,bid) > 1) ? KMT(i,j,bid) : 1);
    
      for (kl = 0; kl<km; kl++) {
        
        if (kl > 0) {
          ZKL = -zgrid(kl-1) + p5* (DZT(i,j,kl,bid) + DZT(i,j,kl-1,bid) );
        }
        else {
          ZKL = -zgrid(0);
        }
      
        //! DIR$ COLLAPSE
        //removed for (j=1;j<ny_block;j+=1)
        //removed for (i=1;i<nx_block;i+=1)
        
        if (kl == KBL (i,j)-1 ) {  //adjusted to c
          HBLT(i,j) = ZKL;
        }
      }


      //! -----------------------------------------------------------------------
      //!
      //! compute velocity shear squared on U-grid and use the maximum
      //! of the four surrounding U-grid values for the T-grid.
      //!
      //! -----------------------------------------------------------------------
      for (kl = 1; kl < km; kl++) {        
        VSHEAR = c0;
      
        if (j >= 1) { //if statement replaced for loop
          
          VSHEAR = c0;
        
          if (i >= 1) { //if statement replaced for loop

	    //requires communication between different threads, split into separate kernel            
            VSHEAR =max( max(VSHEARU(i,j  ,kl), VSHEARU(i-1,j  ,kl)), 
                         max(VSHEARU(i,j-1,kl), VSHEARU(i-1,j-1,kl)) );
          }
        }
      
        ZKL = -zgrid(kl-1) + p5* (DZT(i,j,kl,bid) + DZT(i,j,kl-1,bid) );

//bldepth correct sofar  

        //! -----------------------------------------------------------------------
        //!
        //! compute bfsfc= Bo + radiative contribution down to hbf * hbl
        //! add epsilon to BFSFC to ensure never = 0
        //!
        //! -----------------------------------------------------------------------
        //removed for (j=1;j<ny_block;j+=1)
        //removed for (i=1;i<nx_block;i+=1)
        
        absorb_frac = sw_absorb_frac(ZKL);

      
        BFSFC(i,j) = BO + BOSOL * (c1 - absorb_frac);

        STABLE(i,j) = (BFSFC(i,j) >= c0 ? c1 : c0);
        BFSFC(i,j) = BFSFC(i,j) + STABLE(i,j)*eps;
      

        //! -----------------------------------------------------------------------
        //!
        //! compute the Ekman and Monin Obukhov depths using above stability
        //!
        //! -----------------------------------------------------------------------
        //! -----------------------------------------------------------------------
        //!
        //! compute velocity scales at sigma, for hbl = -zgrid (kl)
        //!
        //! -----------------------------------------------------------------------
        SIGMA = epssfc;
      
        WS = wsscale (SIGMA, ZKL, USTAR(i,j), BFSFC(i,j));
      

        //! -----------------------------------------------------------------------
        //!
        //! compute the turbulent shear contribution to RI_BULK and store
        //! in WM.
        //!
        //! -----------------------------------------------------------------------
        if (kl < km-1) {//adjusted to c
          B_FRQNCY = sqrt ( p5* (DBLOC(i,j,kl) + abs (DBLOC(i,j,kl) ) + eps2) / (p5* (DZT(i,j,kl,bid) + DZT(i,j,kl+1,bid) ) ) );
        }
        else {
          B_FRQNCY = sqrt ( p5* (DBLOC(i,j,kl) + abs (DBLOC(i,j,kl) ) + eps2) / DZT(i,j,kl,bid) );
        }
      
        WM = ZKL*WS*B_FRQNCY* ( (Vtc/Ricr(kl) ) *max (2.1 - 200.0*B_FRQNCY,concv) );
      
        //! -----------------------------------------------------------------------
        //!
        //! compute bulk Richardson number at new level
        //!
        //! -----------------------------------------------------------------------

	double denom = (-zgrid(kl-1) + p5* (DZT(i,j,kl-1,bid) + DZT(i,j,kl,bid) - DZT(i,j,0,bid) ) ); //adjusted to c
        WORK = (KMT(i,j,bid) >= kl+1 ? DBSFC(i,j,kl) / denom : c0); //adjusted to c

	WM = WM / (denom*denom);

	denom = (-zgrid(kl-1) + p5* (DZU(i,j,kl,bid) + DZU(i,j,kl-1,bid) - DZU(i,j,0,bid) ) ); //adjusted to c
        RI_BULK(kdn) = WORK/ (VSHEAR+WM+eps / (denom*denom));
      
        //! -----------------------------------------------------------------------
        //!
        //! find hbl where Rib = Ricr. if possible, use a quadratic
        //! interpolation. if not, linearly interpolate. the quadratic
        //! equation coefficients are determined using the slope and
        //! Ri_bulk at z_up and Ri_bulk at zgrid (kl) . the slope at
        //! z_up is computed linearly between z_upper and z_up.
        //!
        //! compute Langmuir depth always
        //! -----------------------------------------------------------------------
        //removed for (j=1;j<ny_block;j+=1)
        //removed for (i=1;i<nx_block;i+=1)
        
        if ( KBL(i,j) == KMT (i,j,bid) && RI_BULK (kdn) > Ricr (kl) ) {
          slope_up = (RI_BULK(kupper) - RI_BULK(kup) ) / (z_up - z_upper);
          a_co = (RI_BULK(kdn) - RI_BULK(kup) - slope_up* (ZKL + z_up) ) / ((z_up + ZKL)*(z_up + ZKL));
          b_co = slope_up + c2 * a_co * z_up;
          c_co = RI_BULK(kup) + z_up* (a_co*z_up + slope_up) - Ricr(kl);
          sqrt_arg = b_co*b_co - c4*a_co*c_co;
        
          if ( ( abs (b_co) > eps && abs (a_co) /abs (b_co) <= eps ) || sqrt_arg <= c0 ) {
            HBLT(i,j) = -z_up + (z_up + ZKL ) * (Ricr(kl) - RI_BULK(kup) ) / (RI_BULK(kdn) - RI_BULK(kup) );
          }
          else {
            HBLT(i,j) = (-b_co + sqrt (sqrt_arg) ) / (c2*a_co);
          }
        
          KBL(i,j) = kl +1; //adjusted to c
          //RSH_HBLT = (VSHEAR *Ricr(kl) / (DBSFC(i,j,kl) +eps) ) /HBLT(i,j);
        }
      
        //! -----------------------------------------------------------------------
        //!
        //! swap klevel indices and move to next level
        //!
        //! -----------------------------------------------------------------------
        ktmp = kupper;
        kupper = kup;
        kup = kdn;
        kdn = ktmp;
        z_upper = z_up;
        z_up = zgrid(kl);
      }

    
      //here we have split up bldepth because the smoothing of HBLT has to happen in a separate kernel
      //that what comes after the smoothing of HBLT also goes in another kernel,
      //although this last kernel could in theory be merged with the smooth
      d_WORK1[i+j*nx_block] = BO; //store BO and BOSOL as needed for the correction of BFSFC and STABLE after smoothing
      d_WORK2[i+j*nx_block] = BOSOL;
      d_WORK3[i+j*nx_block] = HBLT(i,j);

      //! -----------------------------------------------------------------------
      //! EOC
      }
    }

    /* split from bldepth, STABLE and BFSFC have to be corrected after HBLT was smoothed */
    __global__ void correct_stability_and_buoyancy(double *BFSFC, double *STABLE, double *BO, double *BOSOL, double *HBLT) {

      int j = threadIdx.y + blockIdx.y * blockDim.y;
      int i = threadIdx.x + blockIdx.x * blockDim.x;

      if (j < ny_block && i < nx_block) {

      double bfsfc;
      double stable;

      //! -----------------------------------------------------------------------
      //!
      //! correct stability and buoyancy forcing for SW up to boundary layer
      //!
      //! -----------------------------------------------------------------------
      //removed for (j = 1;j <ny_block;j +=1)
      //removed for (i = 1;i <nx_block;i +=1)
      
      double absorb_frac = sw_absorb_frac(HBLT(i,j));
    
      bfsfc = BO[i+j*nx_block] + BOSOL[i+j*nx_block] * (c1 - absorb_frac);
    
      //removed tavg_requested for tavg_QSW_HBL in GPU version
      //! if (tavg_requested (tavg_QSW_HBL) ) then
      //! WORK = SHF_QSW* (c1-absorb_frac) /hflux_factor
      //! call accumulate_tavg_field (WORK,tavg_QSW_HBL,bid,1)
      //! end if
      stable = (bfsfc >= c0 ? c1 : c0);
      bfsfc = bfsfc + stable * eps;
    
      BFSFC(i,j) = bfsfc;
      STABLE(i,j) = stable;

      }
    }

  
    //! ***********************************************************************
    //! BOP
    //! IROUTINE: blmix
    //! INTERFACE:
    __global__ void blmix_kernel (double *VISC, double *VDC, double *HBLT, double *USTAR, double *BFSFC, double *STABLE, int *KBL, double *GHAT, int bid, double *DZT) {

      int j = threadIdx.y + blockIdx.y * blockDim.y;
      int i = threadIdx.x + blockIdx.x * blockDim.x;

      if (j < ny_block && i < nx_block) {

      //! DESCRIPTION:
      //! This routine computes mixing coefficients within boundary layer
      //! which depend on surface forcing and the magnitude and gradient
      //! of interior mixing below the boundary layer (matching) . These
      //! quantities have been computed in other routines.
      //!
      //! Caution: if mixing bottoms out at hbl = -zgrid (km) then
      //! fictitious layer at km+1 is needed with small but finite width
      //! hwide (km+1) .
      //!
      //! REVISION HISTORY:
      //! same as module
      //! INPUT/OUTPUT PARAMETERS:
      //! combined interior/bndy layer coeff output
      //! INPUT PARAMETERS:
      //! OUTPUT PARAMETERS:
      //! EOP
      //! BOC
      //! -----------------------------------------------------------------------
      //!
      //! local variables
      //!
      //! -----------------------------------------------------------------------
      int k;
      int KN;
      double BLMC[3];
      double GAT1[3];
      double DAT1[3];
      double DKM1[3];
      double WM;
      double WS;
      double CASEA;
      double SIGMA;
      double VISCH;
      double DIFTH;
      double DIFSH;
      double DELHAT;
      double R;
      double DVDZUP;
      double DVDZDN;
      double VISCP;
      double DIFTP;
      double DIFSP;
      double F1;
      double WORK1;
      double WORK2;

      //! -----------------------------------------------------------------------
      //!
      //! compute velocity scales at hbl
      //!
      //! -----------------------------------------------------------------------
      SIGMA = epssfc;
    
      WM = wmscale(SIGMA, HBLT(i,j), USTAR(i,j), BFSFC(i,j));
      WS = wsscale(SIGMA, HBLT(i,j), USTAR(i,j), BFSFC(i,j));

      //! -----------------------------------------------------------------------
      //!
      //! determine caseA = 0 if closer to KBL than KBL-1
      //! KN is then the closest klevel to HBLT
      //!
      //! -----------------------------------------------------------------------
      //! DIR$ COLLAPSE
      //removed for (j=1;j<ny_block;j+=1)
      //removed for (i=1;i<nx_block;i+=1)
      
      k = KBL(i,j)-1;//adjusted to c
      k = max(k,0); //sanity check
    
      if (k == 0) { //adjusted to c
        CASEA = p5 + sign (p5, -zgrid(-1) -HBLT(i,j) ); //adjusted to c
      }
      else {
        CASEA = p5 + sign (p5, -zgrid(k-1) +p5*DZT(i,j,k-1,bid) -HBLT(i,j) );
      }
    
      KN = nint (CASEA) * (KBL(i,j)-1) + (1-nint (CASEA) ) *KBL(i,j);
    
      //! -----------------------------------------------------------------------
      //!
      //! find the interior viscosities and derivatives at hbl by
      //! interpolating derivative values at vertical interfaces. compute
      //! matching conditions for shape function.
      //!
      //! -----------------------------------------------------------------------
      F1 = STABLE(i,j)*c5*BFSFC(i,j)/ (pow(USTAR(i,j),4)+eps);
    
      k = KN-1; //adjusted to c
      k = min(max(k,0),km-1); //sanity check
        
      if (k == 0) {//adjusted to c
        WORK1 = c0;
      }
      else {
        WORK1 = DZT(i,j,k-1,bid);
      }
      
      if (k == km-1) {//adjusted to c
        WORK2 = eps;
      }
      else {
        WORK2 = DZT(i,j,k+1,bid);
      }
      
      //! DIR$ COLLAPSE
      //removed for (j=1;j<ny_block;j+=1)
      //removed for (i=1;i<nx_block;i+=1)

      DELHAT = - zgrid(k-1) + DZT(i,j,k,bid) + p5*WORK1 - HBLT(i,j);
      R = c1 - DELHAT /DZT(i,j,k,bid);
        
      if (k == 0) { //adjusted to c
        DVDZUP = (0.0 - VISC(i,j,k) ) /DZT(i,j,k,bid);
      }
      else {
        DVDZUP = (VISC(i,j,k-1) - VISC(i,j,k) ) /DZT(i,j,k,bid);
      }
        
      DVDZDN = (VISC(i,j,k) - VISC(i,j,k+1) ) /WORK2;
      VISCP = p5* ( (c1-R ) * (DVDZUP + abs (DVDZUP ) ) + R * (DVDZDN + abs (DVDZDN ) ) );
        
      if (k == 0) {//adjusted to c
        DVDZUP = (0.0 - VDC(i,j,k,1) ) /DZT(i,j,k,bid); //adjusted to c
      }
      else {
        DVDZUP = (VDC(i,j,k-1,1) - VDC(i,j,k,1) ) /DZT(i,j,k,bid); //adjusted to c
      }
        
      DVDZDN = (VDC(i,j,k,1) - VDC(i,j,k+1,1) ) /WORK2; //adjusted to c
      DIFSP = p5* ( (c1-R ) * (DVDZUP + abs (DVDZUP ) ) + R * (DVDZDN + abs (DVDZDN ) ) );
        
      if (k == 0) { //adjusted to c
        DVDZUP = (0.0 - VDC(i,j,k,0) ) /DZT(i,j,k,bid); // adjusted to c
      }
      else {
        DVDZUP = (VDC(i,j,k-1,0) - VDC(i,j,k,0) ) /DZT(i,j,k,bid); //adjusted to c
      }
        
      DVDZDN = (VDC(i,j,k,0) - VDC(i,j,k+1,0) ) /WORK2; //adjusted to c
      DIFTP = p5* ( (c1-R ) * (DVDZUP + abs (DVDZUP ) ) + R * (DVDZDN + abs (DVDZDN ) ) );
      VISCH = VISC(i,j,k) + VISCP *DELHAT;
      DIFSH = VDC(i,j,k,1) + DIFSP *DELHAT; //adjusted to c
      DIFTH = VDC(i,j,k,0) + DIFTP *DELHAT; //adjusted to c
      GAT1(0) = VISCH / HBLT(i,j) / (WM +eps); //adjusted to c
      DAT1(0) = -VISCP / (WM +eps) + F1*VISCH; //adjusted to c
      GAT1(1) = DIFSH / HBLT(i,j) / (WS +eps); //adjusted to c
      DAT1(1) = -DIFSP / (WS +eps) + F1*DIFSH; //adjusted to c
      GAT1(2) = DIFTH / HBLT(i,j) / (WS +eps); //adjusted to c
      DAT1(2) = -DIFTP / (WS +eps) + F1*DIFTH; //adjusted to c
    
      DAT1[0] = min (DAT1[0], 0.0); //adjusted to c
      DAT1[1] = min (DAT1[1], 0.0); //adjusted to c
      DAT1[2] = min (DAT1[2], 0.0); //adjusted to c    
//      if (DAT1[0] >= 0) { DAT1[0] = 0.0; }
//      if (DAT1[1] >= 0) { DAT1[1] = 0.0; }
//      if (DAT1[2] >= 0) { DAT1[2] = 0.0; }





//correct so far

      //! -----------------------------------------------------------------------
      //!
      //! find diffusivities at kbl-1 grid level
      //!
      //! -----------------------------------------------------------------------
      //! DIR$ COLLAPSE
      //removed for (j=1;j<ny_block;j+=1)
      //removed for (i=1;i<nx_block;i+=1)
      
      k = KBL(i,j) - 2;//adjusted to c
      k = max(k,-1); //sanity check

      SIGMA = -zgrid(k) /HBLT(i,j);
    
      F1 = min (SIGMA, epssfc);
    
      WM = wmscale (F1, HBLT(i,j), USTAR(i,j), BFSFC(i,j));
      WS = wsscale (F1, HBLT(i,j), USTAR(i,j), BFSFC(i,j));

    
      //! DIR$ COLLAPSE
      //removed for (j=1;j<ny_block;j+=1)
      //removed for (i=1;i<nx_block;i+=1)
      
      DKM1(0) = HBLT(i,j) *WM *SIGMA * (1.0+SIGMA * ( (SIGMA -2.0) + (3.0-2.0*SIGMA ) *GAT1(0) + (SIGMA -1.0) *DAT1(0) ) ); //adjusted to c
      DKM1(1) = HBLT(i,j) *WS *SIGMA * (1.0+SIGMA * ( (SIGMA -2.0) + (3.0-2.0*SIGMA ) *GAT1(1) + (SIGMA -1.0) *DAT1(1) ) ); //adjusted to c
      DKM1(2) = HBLT(i,j) *WS *SIGMA * (1.0+SIGMA * ( (SIGMA -2.0) + (3.0-2.0*SIGMA ) *GAT1(2) + (SIGMA -1.0) *DAT1(2) ) ); //adjusted to c





      //! -----------------------------------------------------------------------
      //!
      //! compute the dimensionless shape functions and diffusivities
      //! at the grid interfaces. also compute function for non-local
      //! transport term (GHAT) .
      //!
      //! -----------------------------------------------------------------------
      for (k = 0; k < km; k++) {
        
        if (k > 0) { //adjusted to c
          SIGMA = (-zgrid(k-1) + p5*DZT(i,j,k-1,bid) + DZT(i,j,k,bid) ) / HBLT(i,j);
        }
        else {
          SIGMA = (-zgrid(k) + p5*hwide(k) ) / HBLT(i,j);
        }
      
        F1 = min (SIGMA,epssfc);
      
        WM = wmscale (F1, HBLT(i,j), USTAR(i,j), BFSFC(i,j));
        WS = wsscale (F1, HBLT(i,j), USTAR(i,j), BFSFC(i,j));

        //! DIR$ COLLAPSE
        //removed for (j=1;j<ny_block;j+=1)
        //removed for (i=1;i<nx_block;i+=1)
        
/* This computation has some difficulty with numerical stability, still experimenting
        double sig3 = SIGMA*SIGMA*SIGMA;
        double sig2 = SIGMA*SIGMA;

        double a = (sig3 - 2.0*sig2);
        double b = (3.0*sig2-2.0*sig3);
        double c = (sig3 - sig2);

        double term1 = SIGMA + a + b * GAT1(0) + c * DAT1(0);
        BLMC[0] = HBLT(i,j) * WM * term1;
        double term2 = SIGMA + a + b * GAT1(1) + c * DAT1(1);
        BLMC[1] = HBLT(i,j) * WS * term2;
        double term3 = SIGMA + a + b * GAT1(2) + c * DAT1(2);
        BLMC[2] = HBLT(i,j) * WS * term3;
*/
          BLMC[0] = HBLT(i,j) * WM * SIGMA * (1.0+SIGMA * ( (SIGMA -2.0) + (3.0-2.0*SIGMA ) *GAT1(0) + (SIGMA -1.0) *DAT1(0) ) ); //adjusted to c
          BLMC[1] = HBLT(i,j) * WS * SIGMA * (1.0+SIGMA * ( (SIGMA -2.0) + (3.0-2.0*SIGMA ) *GAT1(1) + (SIGMA -1.0) *DAT1(1) ) ); //adjusted to c
          BLMC[2] = HBLT(i,j) * WS * SIGMA * (1.0+SIGMA * ( (SIGMA -2.0) + (3.0-2.0*SIGMA ) *GAT1(2) + (SIGMA -1.0) *DAT1(2) ) ); //adjusted to c


        GHAT(i,j,k) = (c1-STABLE(i,j) ) * cg/ (WS * HBLT(i,j) + eps);
    
        //! -----------------------------------------------------------------------
        //!
        //! compute the enhanced mixing
        //!
        //! -----------------------------------------------------------------------
        //! DIR$ NOVECTOR
        if (k < km-1) { //fused with previous for-loop from 0 to k<km
        
          if (k == 0) { //adjusted to c
            WORK1 = -p5*DZT(i,j,k,bid);
          }
          else {
            WORK1 = zgrid(k-1) - p5* (DZT(i,j,k-1,bid) + DZT(i,j,k,bid) );
          }
      
          if (k == (KBL(i,j) - 2) ) { //converted from where statement  //adjusted to c
            DELHAT = (HBLT(i,j) + WORK1) / (p5* (DZT(i,j,k,bid) + DZT(i,j,k+1,bid) ) );
          }
      
          //! DIR$ COLLAPSE
          //removed for (j=1;j<ny_block;j+=1)
          //removed for (i=1;i<nx_block;i+=1)
        
          if (k == (KBL (i,j) - 2) ) { //adjusted to c
  	    double oneminusdelhatsq = (c1-DELHAT)*(c1-DELHAT);
	    double delhatsq = DELHAT*DELHAT;
            BLMC[0] = (c1-DELHAT ) *VISC(i,j,k)  + DELHAT * ( oneminusdelhatsq *DKM1(0) + delhatsq * (CASEA *VISC(i,j,k)  + (c1-CASEA ) *BLMC[0] ) ); //adjusted to c
            BLMC[1] = (c1-DELHAT ) *VDC(i,j,k,1) + DELHAT * ( oneminusdelhatsq *DKM1(1) + delhatsq * (CASEA *VDC(i,j,k,1) + (c1-CASEA ) *BLMC[1] ) ); //adjusted to c
            BLMC[2] = (c1-DELHAT ) *VDC(i,j,k,0) + DELHAT * ( oneminusdelhatsq *DKM1(2) + delhatsq * (CASEA *VDC(i,j,k,0) + (c1-CASEA ) *BLMC[2] ) ); //adjusted to c
            GHAT(i,j,k) = (c1-CASEA ) * GHAT(i,j,k);
          }
        }

        //! -----------------------------------------------------------------------
        //!
        //! combine interior and boundary layer coefficients and nonlocal term
        //!
        //! -----------------------------------------------------------------------
        
        //fused for-loop over k with previous
        //removed for (j=1;j<ny_block;j+=1)
        //removed for (i=1;i<nx_block;i+=1)
        
        if (k < KBL (i,j)-1 ) { //adjusted to c
          VISC(i,j,k) = BLMC[0]; //adjusted to c
          VDC(i,j,k,1) = BLMC[1]; //adjusted to c
          VDC(i,j,k,0) = BLMC[2]; //adjusted to c
        }
        else {
          GHAT(i,j,k) = c0;
        }
 

      } // end of for-loop over k (fused three for-loops)



    
      //! -----------------------------------------------------------------------
      //! EOC
      }
    }
  

/*
 * Split version of wscale 
 * this function computes momentum, which equals wscale with m_or_s = 1
 */

    __device__ __forceinline__ double wmscale (double sigma, double hbl, double ustar, double bfsfc) {

      double zeta;
      double zetah;
      double wm;

      //! -----------------------------------------------------------------------
      //!
      //! compute zetah and zeta - surface layer is special case
      //!
      //! -----------------------------------------------------------------------
      zetah = sigma*hbl*vonkar*bfsfc;
      double ustar_cube = ustar*ustar*ustar;
      zeta = zetah / (ustar_cube + eps);
    
      //! -----------------------------------------------------------------------
      //!
      //! compute velocity scales for momentum
      //!
      //! -----------------------------------------------------------------------

        if (zeta >= c0) {
          wm = vonkar*ustar / (c1 + c5*zeta);
        } else if (zeta >= zeta_m) {
          wm = vonkar*ustar * pow(c1-c16*zeta,p25);
        } else {
          wm = vonkar* pow(a_m* ustar_cube -c_m*zetah,p33);
        }

      return wm;
    }

    __device__ __forceinline__ double wsscale (double sigma, double hbl, double ustar, double bfsfc) {
      double zeta;
      double zetah;
      double ws;

      //! -----------------------------------------------------------------------
      //!
      //! compute zetah and zeta - surface layer is special case
      //!
      //! -----------------------------------------------------------------------
      zetah = sigma*hbl*vonkar*bfsfc;
      double ustar_cube = ustar*ustar*ustar;
      zeta = zetah / (ustar_cube + eps);
    
      //! -----------------------------------------------------------------------
      //!
      //! compute velocity scales for momentum
      //!
      //! -----------------------------------------------------------------------

        if (zeta >= c0) {
          ws = vonkar*ustar / (c1 + c5*zeta );        
        } else if (zeta >= zeta_s) {
          ws = vonkar*ustar * sqrt(c1 - c16*zeta);
        } else {
          ws = vonkar * pow(a_s* ustar_cube -c_s*zetah, p33);
        }

      return ws;
    }



  
    //! ***********************************************************************
    //! BOP
    //! IROUTINE: ddmix
    //! INTERFACE:
    __global__ void ddmix_kernel (double *VDC, double *TRCR, int bid) {

      int j = threadIdx.y + blockIdx.y * blockDim.y;
      int i = threadIdx.x + blockIdx.x * blockDim.x;

      if (j < ny_block && i < nx_block) {  
      //! DESCRIPTION:
      //! $R_\rho$ dependent interior flux parameterization.
      //! Add double-diffusion diffusivities to Ri-mix values at blending
      //! interface and below.
      //!
      //! REVISION HISTORY:
      //! same as module
      //! INPUT PARAMETERS:
      //! INPUT/OUTPUT PARAMETERS:
      //! EOP
      //! BOC
      //! -----------------------------------------------------------------------
      //!
      //! local variables
      //!
      //! -----------------------------------------------------------------------
      int k;
      double ALPHADT;
      double BETADS;
      double RRHO;
      double DIFFDD;
      double PRANDTL;
      double talpha_kup;
      double talpha_knxt;
      double sbeta_kup;
      double sbeta_knxt;
    
      //! -----------------------------------------------------------------------
      //!
      //! compute alpha*DT and beta*DS at interfaces. use RRHO and
      //! PRANDTL for temporary storage for call to state
      //!
      //! -----------------------------------------------------------------------
      double temp_kup = TRCR(i,j,0,0);
      double salt_kup = TRCR(i,j,0,1);
      double temp;
      double salt;

      //PRANDTL = (temp_kup < -c2 ? -c2 : temp_kup); //adjusted to c

      double tq = max(min(temp_kup,TMAX),TMIN);
      double sq = 1000.0 * max(min(salt_kup,SMAX),SMIN);
      double sqr = sqrt(sq);

      double nomk = mwjf_numerator(tq, sq, 0);
      double denomk = mwjf_denominator(tq, sq, sqr, 0);
      talpha_kup = compute_drhodt(tq, sq, sqr, 0, nomk, denomk);
      sbeta_kup = compute_drhods(tq, sq, sqr, 0, nomk, denomk);
    
      for (k = 0; k<km; k++) {
        double vdc1 = VDC(i,j,k,0);
        double vdc2 = VDC(i,j,k,1);

        if ( k < km-1 ) { //adjusted to c
          temp = TRCR(i,j,k+1,0);
          salt = TRCR(i,j,k+1,1);

          //PRANDTL = (temp < -c2 ? -c2 : temp); //adjusted to c
          
          tq = max(min(temp,TMAX),TMIN);
          sq = 1000.0 * max(min(salt,SMAX),SMIN);
          sqr = sqrt(sq);          
          nomk = mwjf_numerator(tq, sq, k+1);
          denomk = mwjf_denominator(tq, sq, sqr, k+1);
          talpha_knxt = compute_drhodt(tq, sq, sqr, k+1, nomk, denomk);
          sbeta_knxt = compute_drhods(tq, sq, sqr, k+1, nomk, denomk);
          RRHO = nomk * denomk;

          ALPHADT = -p5* (talpha_kup + talpha_knxt ) * (temp_kup - temp );
          BETADS = p5* ( sbeta_kup + sbeta_knxt ) * (salt_kup - salt ); 

        } else {        
          ALPHADT = c0;
          BETADS = c0;
        }
      
        //! -----------------------------------------------------------------------
        //!
        //! salt fingering case
        //!
        //! -----------------------------------------------------------------------
        if ( ALPHADT > BETADS && BETADS > c0 ) { //converted from where statement
          RRHO = min (ALPHADT/BETADS, Rrho0);
	  double tmp = (c1-(RRHO-c1)/(Rrho0-c1));
          DIFFDD = dsfmax* (tmp*tmp*tmp);
          vdc1 = vdc1 + 0.7*DIFFDD;  //adjusted to c
          vdc2 = vdc2 + DIFFDD; //adjusted to c
        }
      
        //! -----------------------------------------------------------------------
        //!
        //! diffusive convection
        //!
        //! -----------------------------------------------------------------------
        if ( ALPHADT < c0 && BETADS < c0 && ALPHADT > BETADS ) { //converted from where statement
          RRHO = ALPHADT / BETADS;
          DIFFDD = 1.5e-2*0.909* exp (4.6*exp (-0.54* (c1/RRHO-c1) ) );
          PRANDTL = 0.15*RRHO;
        } else {
          RRHO = c0;
          DIFFDD = c0;
          PRANDTL = c0;
        }
      
        if (RRHO > p5) { //converted from where statement
          //PRANDTL = (1.85 - 0.85/RRHO) *RRHO;
          PRANDTL = 1.85*RRHO - 0.85;
        }
      
        VDC(i,j,k,0) = vdc1 + DIFFDD;  //adjusted to c
        VDC(i,j,k,1) = vdc2 + PRANDTL*DIFFDD;  //adjusted to c

        //reuse these for next round
        temp_kup = temp;
        salt_kup = salt;
        talpha_kup = talpha_knxt;
        sbeta_kup = sbeta_knxt;

      }
    
      //! -----------------------------------------------------------------------
      //! EOC
      }
    }
  
    //! ***********************************************************************
    //! BOP
    //! IROUTINE: buoydiff
    //! INTERFACE:
    __global__ void buoydiff_kernel (double *DBLOC, double *DBSFC, double *TRCR, int bid, int *KMT) {

      int j = threadIdx.y + blockIdx.y * blockDim.y;
      int i = threadIdx.x + blockIdx.x * blockDim.x;
      
      if (j < ny_block && i < nx_block) {
      //! DESCRIPTION:
      //! This routine calculates the buoyancy differences at model levels.
      //!
      //! REVISION HISTORY:
      //! same as module
      //! INPUT PARAMETERS:
      //! OUTPUT PARAMETERS:
      //! EOP
      //! BOC
      //! -----------------------------------------------------------------------
      //!
      //! local variables
      //!
      //! -----------------------------------------------------------------------
      int k;
      double RHO1;
      double RHOKM;
      double RHOK;
      double TEMPSFC = TRCR(i,j,0,0);
      double SALTSFC = TRCR(i,j,0,1);
      int kmt = KMT(i,j,bid);
    
      //! -----------------------------------------------------------------------
      //!
      //! calculate density and buoyancy differences at surface
      //!
      //! -----------------------------------------------------------------------
      DBSFC(i,j,0) = 0.0; //adjusted to c
    
      //! -----------------------------------------------------------------------
      //!
      //! calculate DBLOC and DBSFC for all other levels
      //!
      //! -----------------------------------------------------------------------
      for (k = 1; k < km; k++) {
        
        RHO1 = compute_rho(TEMPSFC, SALTSFC, k);
      
        RHOKM = compute_rho(TRCR(i,j,k-1,0), TRCR(i,j,k-1,1), k);
      
        RHOK = compute_rho(TRCR(i,j,k,0), TRCR(i,j,k,1), k);
      
        //removed for (j=1;j<ny_block;j+=1)
        //removed for (i=1;i<nx_block;i+=1)
        
        if (RHOK != 0.0) {
          DBSFC(i,j,k) = grav* (1.0 - RHO1 / RHOK );
          DBLOC(i,j,k-1) = grav* (1.0 - RHOKM / RHOK );
        }
        else {
          DBSFC(i,j,k) = 0.0;
          DBLOC(i,j,k-1) = 0.0;
        }
      
        if (k >= kmt) { //adjusted to c
          DBLOC(i,j,k-1) = 0.0;
        }

      }
    
      DBLOC(i,j,km-1) = 0.0;//adjusted to c
    
      //! -----------------------------------------------------------------------
      //! EOC
      }
    }
  
    //! ***********************************************************************
    //! BOP
    //! IROUTINE: add_kpp_sources
    //! INTERFACE:
    //! ***********************************************************************
    //! BOP
    //! IROUTINE: smooth_hblt
    //! INTERFACE:
    __global__ void smooth_hblt_kernel(int bid, double *HBLT, int *KBL, double *WORK3, int *KMT, double *DZT) {

      int j = threadIdx.y + blockIdx.y * blockDim.y;
      int i = threadIdx.x + blockIdx.x * blockDim.x;
      

      if (j < ny_block && i < nx_block) {


      //! DESCRIPTION:
      //! This subroutine uses a 1-1-4-1-1 Laplacian filter one time
      //! on HBLT or HMXL to reduce any horizontal two-grid-point noise.
      //! If HBLT is overwritten, KBL is adjusted after smoothing.
      //!
      //! REVISION HISTORY:
      //! same as module
    
      //! INPUT/OUTPUT PARAMETERS:
      //! OUTPUT PARAMETERS:
      //! EOP
      //! BOC
      //! -----------------------------------------------------------------------
      //!
      //! local variables
      //!
      //! -----------------------------------------------------------------------
    
      int k;
      double WORK2;
      double cc;
      double cw;
      double ce;
      double cn;
      double cs;
      double ztmp;
    
      //! -----------------------------------------------------------------------
      //!
      //! perform one smoothing pass since we cannot do the necessary
      //! boundary updates for multiple passes.
      //!
      //! -----------------------------------------------------------------------
      WORK2 = WORK3(i,j);
    
      if (j > 0 && j < ny_block-1 && i > 0 && i < nx_block-1) {
          if ( KMT(i,j,bid) != 0 ) {
            cw = p125;
            ce = p125;
            cn = p125;
            cs = p125;
            cc = p5;
          
            if ( KMT(i-1,j,bid) == 0 ) {
              cc = cc + cw;
              cw = c0;
            }
          
            if ( KMT(i+1,j,bid) == 0 ) {
              cc = cc + ce;
              ce = c0;
            }
          
            if ( KMT(i,j-1,bid) == 0 ) {
              cc = cc + cs;
              cs = c0;
            }
          
            if ( KMT(i,j+1,bid) == 0 ) {
              cc = cc + cn;
              cn = c0;
            }

            WORK2 = cw * WORK3(i-1,j) + ce * WORK3(i+1,j) + cs * WORK3(i,j-1) + cn * WORK3(i,j+1) + cc * WORK3(i,j);
          }
    
        k = KMT (i,j,bid)-1; //adjusted to c
        if (KMT (i,j,bid) > 0) {
          ztmp = -zgrid(k-1) + p5* (DZT(i,j,k-1,bid) + DZT(i,j,k,bid) );
          if (WORK2 > ztmp) {
            WORK2 = ztmp;
          }
        }

     
        int kbl = KBL[i+j*nx_block];
        for (k = 0; k<km; k++) {
          ztmp = -zgrid(k-1) + p5* (DZT(i,j,k-1,bid) + DZT(i,j,k,bid) );

          if ( KMT (i,j,bid) != 0 && ( WORK2 > -zgrid (k-1) ) && ( WORK2 <= ztmp ) ) {
            kbl = k + 1;
          }
        }
        KBL[i+j*nx_block] = kbl; //adjusted to c

      }

      HBLT(i,j) = WORK2;

      //! -----------------------------------------------------------------------
      }
    }  
    //! ***********************************************************************



  //! |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||;




    //! ***********************************************************************

    __global__ void ugrid_to_tgrid_kernel (double *ARRAY_TGRID, double *ARRAY_UGRID, int k) {

      int j = threadIdx.y + blockIdx.y * blockDim.y;
      int i = threadIdx.x + blockIdx.x * blockDim.x;

      if (j < ny_block && i < nx_block) {

      //! DESCRIPTION:
      //! Interpolates values at U points on a B grid to T points.
      //! Note that ghost cells are not updated.
      //! Also note that the input array is assumed to be in the baroclinic
      //! distribution (where the stencil weights are defined) .
      //!
      //! REVISION HISTORY:
      //! same as module
      //! INPUT PARAMETERS:
      //! OUTPUT PARAMETERS:
      //! EOP
      //! BOC
      //! -----------------------------------------------------------------------
      //!
      //! local variables
      //!
      //! -----------------------------------------------------------------------
      //! -----------------------------------------------------------------------
      //!
      //! southwest 4pt average
      //!
      //! -----------------------------------------------------------------------
      double at0 = 0.25;
      double ats = 0.25;
      double atw = 0.25;
      double atsw = 0.25;
      double tgrid = 0.0;

      if (j > 0) { //if statement replaced for loop
        if (i > 0) { //if statement replaced for loop
          tgrid = at0 *ARRAY_UGRID(i,j,k) +
		  ats *ARRAY_UGRID(i,j-1,k) +
		  atw *ARRAY_UGRID(i-1,j,k) +
		  atsw*ARRAY_UGRID(i-1,j-1,k);
        }
      }
  

      ARRAY_TGRID(i,j,k) = tgrid;  
      }
    }
  


    __device__ double ugrid_to_tgrid (double *ARRAY_UGRID, int i, int j, int k) {

      //! DESCRIPTION:
      //! Interpolates values at U points on a B grid to T points.
      //! Note that ghost cells are not updated.
      //! Also note that the input array is assumed to be in the baroclinic
      //! distribution (where the stencil weights are defined) .
      //!
      //! REVISION HISTORY:
      //! same as module
      //! INPUT PARAMETERS:
      //! OUTPUT PARAMETERS:
      //! EOP
      //! BOC
      //! -----------------------------------------------------------------------
      //!
      //! local variables
      //!
      //! -----------------------------------------------------------------------
      //! -----------------------------------------------------------------------
      //!
      //! southwest 4pt average
      //!
      //! -----------------------------------------------------------------------
      double at0 = 0.25;
      double ats = 0.25;
      double atw = 0.25;
      double atsw = 0.25;
      double tgrid = 0.0;

      if (j > 0) { //if statement replaced for loop
        if (i > 0) { //if statement replaced for loop
          tgrid = at0 *ARRAY_UGRID(i,j,k) +
		  ats *ARRAY_UGRID(i,j-1,k) +
		  atw *ARRAY_UGRID(i-1,j,k) +
		  atsw*ARRAY_UGRID(i-1,j-1,k);
        }
      }
  
      return tgrid;  
    }
  
    //! ***********************************************************************
    //! BOP
    //! IROUTINE: tgrid_to_ugrid
    //! INTERFACE:
/*
 * GPU version uses modified coefficients:
 *
 * original version uses:
   AU0  = TAREA
   AUN  = eoshift(TAREA,dim=2,shift=+1)
   AUE  = eoshift(TAREA,dim=1,shift=+1)
   AUNE = eoshift(AUE  ,dim=2,shift=+1)

   AU0  = AU0 *p25*UAREA_R
   AUN  = AUN *p25*UAREA_R
   AUE  = AUE *p25*UAREA_R
   AUNE = AUNE*p25*UAREA_R
 *
 * GPU version uses
   AU = TAREA(i,j) *p25 
   ugrid *= UAREA_R
 *
 */
   __device__ double tgrid_to_ugrid(double *ARRAY_TGRID, double *AU, double *UAREA_R, int iblock, int i, int j, int k) {

      double ugrid = 0.0;
      
      //! DESCRIPTION:
      //! Interpolates values at T points on a B grid to U points.
      //! Note that ghost cells are not updated.
      //! Also note that the input array is assumed to be in the baroclinic
      //! distribution (where the stencil weights are defined).
      //!
      //! -----------------------------------------------------------------------
      //!
      //! northeast 4pt average
      //!
      //! -----------------------------------------------------------------------
      if (i < nx_block-1) {
        if (j < ny_block-1) {
          double uarea_r = UAREA_R(i,j,iblock);
          ugrid =	AU(i,j,iblock)    * ARRAY_TGRID(i,j,k) +
			AU(i,j+1,iblock)  * ARRAY_TGRID(i,j+1,k) +
			AU(i+1,j,iblock)  * ARRAY_TGRID(i+1,j,k) +
			AU(i+1,j+1,iblock)* ARRAY_TGRID(i+1,j+1,k);
          ugrid *= uarea_r;
        }
      }

      return ugrid;
  
      //! -----------------------------------------------------------------------
    }
  
    //! ***********************************************************************
    __device__ double sw_absorb_frac(double depth) {
      
      //! DESCRIPTION:
      //! Computes fraction of solar short-wave flux penetrating to
      //! specified depth due to exponential decay in Jerlov water type.
      //! Reference : two band solar absorption model of Simpson and
      //! Paulson (1977)
      //! Note: below 200m the solar penetration gets set to zero,
      //! otherwise the limit for the exponent ($+/- 5678$) needs to be
      //! taken care of.
      //!
      //! REVISION HISTORY:
      //! same as module
      //! INPUT PARAMETERS:
      //! OUTPUT PARAMETERS:
      double sw_absorb_fraction;
      //! EOP
      //! BOC
      //! -----------------------------------------------------------------------
      //!
      //! local variables
      //!
      //! -----------------------------------------------------------------------
      double depth_cutoff = -200.0;
      double depth_neg_meters;
      //! -----------------------------------------------------------------------
      //!
      //! define Jerlov water properties with rfac, depth1, depth2
      //! Jerlov water type : I IA IB II III
      //! jerlov_water_type : 1 2 3 4 5
      //!
      //! -----------------------------------------------------------------------
      double rfac = 0.67;
      double depth1 = 1.00;
      double depth2 = 17.0;
    
      //! -----------------------------------------------------------------------
      //!
      //! compute absorption fraction
      //!
      //! -----------------------------------------------------------------------
      depth_neg_meters = -depth*0.01;//mpercm
    
      //! change sign
      if (depth_neg_meters < depth_cutoff) {
        sw_absorb_fraction = 0.0;
      } else {
        sw_absorb_fraction = rfac * exp(depth_neg_meters/depth1) + (c1 - rfac) * exp(depth_neg_meters/depth2);
      }
    
      //! -----------------------------------------------------------------------
      //! EOC
      return sw_absorb_fraction;
    }



  //! |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||;


  __device__ __forceinline__ double compute_rho(double temp, double salt, int k) {
     double tq = max(min(temp,TMAX),TMIN);
     double sq = 1000.0 * max(min(salt,SMAX),SMIN);
     return mwjf_numerator(tq, sq, k) * mwjf_denominator(tq, sq, sqrt(sq), k);
  }

  __device__ __forceinline__ double mwjf_numerator(double tq, double sq, int k) {
        return mwjfnums0t0[k] + tq * (mwjfnums0t1 + tq * (mwjfnums0t2[k] + mwjfnums0t3 * tq)) +
                              sq * (mwjfnums1t0[k] + mwjfnums1t1 * tq + mwjfnums2t0 * sq);
  }

  __device__ __forceinline__ double mwjf_denominator(double tq, double sq, double sqr, int k) {
        double work2 = mwjfdens0t0[k] + tq * (mwjfdens0t1[k] + tq * (mwjfdens0t2 +
           tq * (mwjfdens0t3[k] + mwjfdens0t4 * tq))) +
           sq * (mwjfdens1t0 + tq * (mwjfdens1t1 + tq*tq*mwjfdens1t3)+
           sqr * (mwjfdensqt0 + tq*tq*mwjfdensqt2));
        return 1.0/work2;
  }

  __device__ __forceinline__ double compute_drhodt(double tq, double sq, double sqr, int k, double nomk, double denomk) {
          double work3 = // dP_1/dT
                 mwjfnums0t1 + tq * (2.0*mwjfnums0t2[k] +
                 3.0*mwjfnums0t3 * tq) + mwjfnums1t1 * sq;

          double work4 = // dP_2/dT
                 mwjfdens0t1[k] + sq * mwjfdens1t1 +
                 tq * (2.0*(mwjfdens0t2 + sq*sqr*mwjfdensqt2) +
                 tq * (3.0*(mwjfdens0t3[k] + sq * mwjfdens1t3) +
                 tq *  4.0*mwjfdens0t4));

          return (work3 - nomk*denomk*work4)*denomk;
  }

  __device__ __forceinline__ double compute_drhods(double tq, double sq, double sqr, int k, double nomk, double denomk) {
          double work3 = // dP_1/dS
                 mwjfnums1t0[k] + mwjfnums1t1 * tq + 2.0*mwjfnums2t0 * sq;

          double work4 = mwjfdens1t0 +   // dP_2/dS
                 tq * (mwjfdens1t1 + tq*tq*mwjfdens1t3) +
                 1.5*sqr*(mwjfdensqt0 + tq*tq*mwjfdensqt2);

          return (work3 - nomk*denomk*work4)*denomk * 1000.0;
  }




//fortran entry
void fill_random(double *p, int *n) {
  int end = *n;
  int i;
  for (i=0;i<end;i++) {
    p[i] = ( rand() % 100000 ) / 100000.0; 
  }
}


void fillrandom(double *p, int n) {
  int i;
  for (i=0;i<n;i++) {
    p[i] = (1+(rand() % 100000)) / 100000.0; 
  }
}


void fillrandom_int(int *p, int n) {
  int i;
  for (i=0;i<n;i++) {
    p[i] = ( rand() % km );
  }
}


/*
int main() {

  cudaError_t err;

  cudaSetDeviceFlags(cudaDeviceMapHost|cudaDeviceScheduleSpin);
  cudaSetDevice(0);
  cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);
  cudaDeviceSynchronize();

//  srand(0);//time(NULL));
srand(time(NULL));

  printf("cuda initialized\n");

  double *h_DZT		= (double *) malloc(nx_block * ny_block * (km+2) * sizeof(double));
  double *h_HMXL	= (double *) malloc(nx_block * ny_block * max_blocks_clinic * sizeof(double));
  double *h_DZU		= (double *) malloc(nx_block * ny_block * (km+2) * sizeof(double));
  double *h_bckgrnd_vdc = (double *) malloc(nx_block * ny_block * km * max_blocks_clinic * sizeof(double));
  double *h_bckgrnd_vvc = (double *) malloc(nx_block * ny_block * km * max_blocks_clinic * sizeof(double));
  double *h_Ricr	= (double *) malloc(nx_block * ny_block * km * sizeof(double));
  double *h_AU		= (double *) malloc(nx_block * ny_block * max_blocks_clinic * sizeof(double));
  int *h_KMU		= (int *) malloc(nx_block * ny_block * max_blocks_clinic * sizeof(int));
  int *h_KMT		= (int *) malloc(nx_block * ny_block * max_blocks_clinic * sizeof(int));

  double *pressz	= (double *) malloc(km * sizeof(double));
  double *h_dz		= (double *) malloc(km * sizeof(double));
  double *h_zt		= (double *) malloc(km * sizeof(double));
  double *h_hwide	= (double *) malloc((km+2) * sizeof(double));
  double *h_zgrid	= (double *) malloc((km+2) * sizeof(double));

  printf("allocated host data\n");

  fillrandom(h_DZT, nx_block * ny_block * (km+2));
  fillrandom(h_HMXL, nx_block * ny_block * max_blocks_clinic);
  fillrandom(h_DZU, nx_block * ny_block * (km+2));
  fillrandom(h_bckgrnd_vdc, nx_block * ny_block * km * max_blocks_clinic);
  fillrandom(h_bckgrnd_vdc, nx_block * ny_block * km * max_blocks_clinic);
  fillrandom(h_Ricr, nx_block * ny_block * km);
  fillrandom(h_AU, nx_block * ny_block * max_blocks_clinic);
//  fillrandom_int(h_KMT, nx_block * ny_block * max_blocks_clinic);

  //construct a random yet somewhat realistic KMT
  int i=0,j=0,b=0;
  int bl=8;
  for (b=0; b < max_blocks_clinic; b++) {
    for (j=0; j < ny_block; j+=bl) {
      for (i=0; i < nx_block; i+=bl) {
        int r = rand() % km;
        int ib=i,jb=j;
        for (jb=j; jb < j+bl && jb < ny_block; jb++) {
          for (ib=i; ib < i+bl && ib < nx_block; ib++) {
            h_KMT(ib,jb,b) = r;
          }
        }
      }
    }
  }

  //smooth KMT
  for (int f=0; f < 10; f++) {
  for (b=0; b < max_blocks_clinic; b++) {
    for (j=1; j < ny_block-1; j++) {
      for (i=1; i < nx_block-1; i++) {
        h_KMT(i,j,b) = min( ( h_KMT(i,j,b) + h_KMT(i+1,j,b) + h_KMT(i,j+1,b) + h_KMT(i-1,j,b) + h_KMT(i,j-1,b) ) / 5 , km );
      }
    }
  }
  }

  printf("KMT\n");
  for (j=0; j < 32; j++) {
    for (i=0; i < 32; i++) {
      printf("%.2d ", h_KMT(i,j,0));
    }
    printf("\n");
  }
  printf("\n");

  //fill KMU
  memcpy(h_KMU, h_KMT, nx_block * ny_block * max_blocks_clinic * sizeof(int));
  for (b=0; b < max_blocks_clinic; b++) {
    for (j=0; j < ny_block-1; j++) {
      for (i=0; i < nx_block-1; i++) {
        h_KMU(i,j,b) = min ( min ( h_KMT(i, j  , b), h_KMT(i+1, j  , b) ) ,
                             min ( h_KMT(i, j+1, b), h_KMT(i+1, j+1, b) ) );
      }
    }
  }

  fillrandom(pressz, km);
  fillrandom(h_dz, km);
  fillrandom(h_zt, km);
  fillrandom(h_hwide, (km+2));
  fillrandom(h_zgrid, (km+2));

  printf("filled globals in host memory with test data\n");

  printf ("before h_bckgrnd_vdc[15]=%.20f\n",h_bckgrnd_vdc[15]);

  init_vmix_kpp(h_DZT, h_KMU, h_dz, h_HMXL, h_zt, h_DZU, h_KMT, h_bckgrnd_vdc, h_bckgrnd_vvc, h_zgrid, h_Ricr, h_hwide, pressz, h_AU);

  printf ("after h_bckgrnd_vdc[15]=%.20f\n",h_bckgrnd_vdc[15]);
  printf ("after h_bckgrnd_vvc[15]=%.20f\n",h_bckgrnd_vvc[15]);

  printf("allocated and copied global data to and from Fortran memory\n");

  init_global_variables(h_DZT, h_KMU, h_dz, h_HMXL, h_zt, h_DZU, h_KMT, h_bckgrnd_vdc, h_bckgrnd_vvc, h_zgrid, h_Ricr, h_hwide, pressz, h_AU);

  printf("allocated and copied global GPU data\n");


  double *VDC;
  double *VVC;
  double *TRCR;
  double *UUU;
  double *VVV;
  double *STF;
  double *SHF_QSW;
  double *SMF;
  double *HMXL;
  double *KPP_HBLT;
  double *KPP_SRC;

  err = cudaHostAlloc ((void **) &TRCR, nx_block * ny_block * km * nt * sizeof(double), cudaHostAllocMapped);
  if (err != cudaSuccess) { fprintf (stderr, "Error in cudaHostAlloc TRCR: %s\n", cudaGetErrorString (err)); }
  err = cudaHostAlloc ((void **) &UUU, nx_block * ny_block * km * sizeof(double), cudaHostAllocMapped);
  if (err != cudaSuccess) { fprintf (stderr, "Error in cudaHostAlloc UUU: %s\n", cudaGetErrorString (err)); }
  err = cudaHostAlloc ((void **) &VVV, nx_block * ny_block * km * sizeof(double), cudaHostAllocMapped);
  if (err != cudaSuccess) { fprintf (stderr, "Error in cudaHostAlloc VVV: %s\n", cudaGetErrorString (err)); }
  err = cudaHostAlloc ((void **) &STF, nx_block * ny_block * nt * sizeof(double), cudaHostAllocMapped);
  if (err != cudaSuccess) { fprintf (stderr, "Error in cudaHostAlloc STF: %s\n", cudaGetErrorString (err)); }
  err = cudaHostAlloc ((void **) &SHF_QSW, nx_block * ny_block * sizeof(double), cudaHostAllocMapped);
  if (err != cudaSuccess) { fprintf (stderr, "Error in cudaHostAlloc SHF_QSW: %s\n", cudaGetErrorString (err)); }
  err = cudaHostAlloc ((void **) &SMF, nx_block * ny_block * 2 * sizeof(double), cudaHostAllocMapped);
  if (err != cudaSuccess) { fprintf (stderr, "Error in cudaHostAlloc SMF: %s\n", cudaGetErrorString (err)); }
  err = cudaHostAlloc ((void **) &VDC, nx_block * ny_block * (km+1) * 2 * sizeof(double), cudaHostAllocMapped);
  if (err != cudaSuccess) { fprintf (stderr, "Error in cudaHostAlloc VDC: %s\n", cudaGetErrorString (err)); }
  err = cudaHostAlloc ((void **) &VVC, nx_block * ny_block * km * sizeof(double), cudaHostAllocMapped);
  if (err != cudaSuccess) { fprintf (stderr, "Error in cudaHostAlloc VVC: %s\n", cudaGetErrorString (err)); }
  err = cudaHostAlloc ((void **) &HMXL, nx_block * ny_block * sizeof(double), cudaHostAllocMapped);
  if (err != cudaSuccess) { fprintf (stderr, "Error in cudaHostAlloc HMXL: %s\n", cudaGetErrorString (err)); }
  err = cudaHostAlloc ((void **) &KPP_HBLT, nx_block * ny_block * sizeof(double), cudaHostAllocMapped);
  if (err != cudaSuccess) { fprintf (stderr, "Error in cudaHostAlloc KPP_HBLT: %s\n", cudaGetErrorString (err)); }
  err = cudaHostAlloc ((void **) &KPP_SRC, nx_block * ny_block * km * nt * sizeof(double), cudaHostAllocMapped);
  if (err != cudaSuccess) { fprintf (stderr, "Error in cudaHostAlloc KPP_SRC: %s\n", cudaGetErrorString (err)); }

  int bid = 1;
  double convect_diff = 1000.0;
  double convect_visc = 1000.0;

  printf("finished allocation of host data\n");

  fillrandom(TRCR, nx_block * ny_block * km * nt);
  fillrandom(UUU, nx_block * ny_block * km);
  fillrandom(VVV, nx_block * ny_block * km);
  fillrandom(STF, nx_block * ny_block * nt);
  fillrandom(SHF_QSW, nx_block * ny_block);
  fillrandom(SMF, nx_block * ny_block * 2);
  fillrandom(VDC, nx_block * ny_block * (km+1) * 2);
  fillrandom(VVC, nx_block * ny_block * km);
  fillrandom(HMXL, nx_block * ny_block);
  fillrandom(KPP_HBLT, nx_block * ny_block);
  fillrandom(KPP_SRC, nx_block * ny_block * km * nt);

  printf("finished filling host data with random values\n");



  vmix_coeffs_kpp_gpu_entry_test(VDC, VVC, TRCR, UUU, VVV, STF, SHF_QSW, &bid, &convect_diff, 
				&convect_visc, SMF, HMXL, KPP_HBLT, KPP_SRC);



  return 0;
}

*/

int compare_int (int *a1, int *a2, int N) {
  int i=0,res=0,print=0;

  for (i=0; i<N; i++) {
    if (a1[i] != a2[i]) {
      res++;
      if (print < 10) {
        print++;
        fprintf(stderr, "Error detected at i=%d, \t a1= \t %d \t a2= \t %d \n",i,a1[i],a2[i]);
      }

    }
  }

  if (res > 0) { fprintf(stdout,"Number of errors in GPU result: %d\n",res); }

  return res;
}




int compare (double *a1, double *a2, int N, const char * str) {
  int i=0, res=0;
  int print = 0;
  int zero_one = 0;
  int zero_two = 0;

  double eps = 0.000001;

  for (i=0; i<N; i++) {

    if (a1[i] < eps && a1[i] > -eps) { zero_one++; }
    if (a2[i] < eps && a2[i] > -eps) { zero_two++; }

    if (isnan(a1[i]) || isnan(a2[i])) {
        res++;
        if (print < 10) {
          print++;
          fprintf(stderr, "Error in %s isnan at i=%d, a1= %30.27e a2= %30.27e\n",str,i,a1[i],a2[i]);
        }
    }

    double diff1 = abs( a1[i]-a2[i] ) / a2[i];
    double diff2 = abs( a1[i]-a2[i] ) / a1[i];
    double diff = max(diff1, diff2);
    if (diff > eps) {
        res++;
        if (print < 10) {
          print++;
         /*
          unsigned long long int int_a1 = *(unsigned long long int *)(a1+i);
          unsigned long long int int_a2 = *(unsigned long long int *)(a2+i);
          unsigned long long int dist = (unsigned long long int)0;
          if (int_a1 > int_a2) {
            dist = int_a1 - int_a2;
          } else {
            dist = int_a2 - int_a1;
          }
          fprintf(stderr, "Error detected at i=%d, \t a1= \t %30.27e \t a2= \t %30.27e \t ulp_dist=\t %llu\n",i,a1[i],a2[i],dist);
         */

          fprintf(stderr, "Error in %s at i=%d, \t a1= \t %30.27e \t a2= \t %30.27e \n",str,i,a1[i],a2[i]);
        }
    }

  }

  if (zero_one > 3*(N/4)) { fprintf(stderr, "Error: array1 contains %d zeros\n", zero_one); }
  if (zero_two > 3*(N/4)) { fprintf(stderr, "Error: array2 contains %d zeros\n", zero_two); }

  if (zero_one != zero_two) {
    fprintf(stderr, "Error in %s number of zeros in arrays dont correspond zero1=%d, zero2=%d\n", str, zero_one, zero_two);
  }

  if (res > 0) { fprintf(stdout,"Number of errors in GPU result %s: %d\n",str,res); }

  return res;
}




int compare_to (double *a1, double *a2, int N, double *a3) {
  int i=0, res=0;
  int print = 0;
  int zero_one = 0;
  int zero_two = 0;

  double eps = 0.000001;

  for (i=0; i<N; i++) {

    if (a1[i] < eps && a1[i] > -eps) { zero_one++; }
    if (a2[i] < eps && a2[i] > -eps) { zero_two++; }

    if (isnan(a1[i]) || isnan(a2[i])) {
        res++;
        if (print < 10) {
          print++;
          fprintf(stderr, "Error detected isnan at i=%d, a1= %30.27e a2= %30.27e, a3=%30.27e\n",i,a1[i],a2[i],a3[i]);
        }
    }

    double diff1 = abs( a1[i]-a2[i] ) / a2[i];
    double diff2 = abs( a1[i]-a2[i] ) / a1[i];
    double diff = max(diff1, diff2);
    if (diff > eps) {
        res++;
        if (print < 10) {
          print++;
         /*
          unsigned long long int int_a1 = *(unsigned long long int *)(a1+i);
          unsigned long long int int_a2 = *(unsigned long long int *)(a2+i);
          unsigned long long int dist = (unsigned long long int)0;
          if (int_a1 > int_a2) {
            dist = int_a1 - int_a2;
          } else {
            dist = int_a2 - int_a1;
          }
          fprintf(stderr, "Error detected at i=%d, \t a1= \t %30.27e \t a2= \t %30.27e \t ulp_dist=\t %llu\n",i,a1[i],a2[i],dist);
         */

          fprintf(stderr, "Error detected at i=%d, a1= %30.27e a2= %30.27e a3= %30.27e\n",i,a1[i],a2[i],a3[i]);
        }
    }

  }

  if (zero_one > 3*(N/4)) { fprintf(stderr, "Error: array1 contains %d zeros\n", zero_one); }
  if (zero_two > 3*(N/4)) { fprintf(stderr, "Error: array2 contains %d zeros\n", zero_two); }

  if (zero_one != zero_two) {
    fprintf(stderr, "Error: number of zeros in arrays dont correspond zero1=%d, zero2=%d\n", zero_one, zero_two);
  }

  if (res > 0) { fprintf(stdout,"Number of errors in GPU result: %d\n",res); }

  return res;
}


