!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  module gpu_mod

!BOP
! MODULE: gpu_mod
! DESCRIPTION:
!  This module contains the routines for execution on GPU accelerators
!
! REVISION HISTORY:
!  gpu_mod.F90 2012-10-11 19:30:26Z B. van Werkhoven

! USES:

   use kinds_mod
   use io
   use exit_mod
   use constants
   use grid
   use state_mod     ! access to preszz, tmin, tmax, smin, smax, etc
   use domain_size   ! included for use of nx_block,ny_block,km,
   use prognostic    ! include for reference to TRACER
   use vertical_mix  ! include for ference to VDC and VVC
   use iso_c_binding

   implicit none
   private
   save

   ! include compile time domain info
   #include "gpu_domain.h"
   ! include the interface to the C code
   #include "gpu_mod.fh"

! !PUBLIC MEMBER FUNCTIONS:

   public :: init_gpu_mod, &
             gpumod_compare, &
             gpumod_devsync, &
             gpumod_mwjf_state, &
             gpumod_mwjf_statePD, &
             gpumod_buoydiff, &
             gpumod_ddmix

!             vmix_coeffs,                          &
!             vdifft, vdiffu,                       &
!             impvmixt, impvmixt_correct, impvmixu, &
!             convad

! !PUBLIC DATA MEMBERS:

   logical (log_kind), public :: &
      use_gpu_state   ! flag for computing density functions on gpu

!EOP
!BOC

!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: init_vertical_mix
! !INTERFACE:

 subroutine init_gpu_mod

! !DESCRIPTION:
!  Initializes various mixing quantities and calls initialization
!  routines for specific parameterizations.
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   type (block) :: &
      this_block          ! block information for current block

   integer (int_kind) ::  &
      k,                  &! vertical level index
      nu,                 &! i/o unit
      nml_error,          &! namelist i/o error flag
      bid                  ! block id

   character (char_len) :: &
      vmix_choice          ! input choice for desired parameterization

   character (char_len) :: &
      convection_type      ! input choice for method for convection

   type(c_ptr) :: &
      cptr                  ! ptr used for alllocating arrays

   namelist /gpu_mod_nml/ use_gpu_state, use_verify_results

   ! array used for sending constants to the GPU
   ! this prevents issues with compilers rounding double precision constants
   ! differently and thus reduces rounding errors
   real (r8), dimension(46) :: &
        constants = (/ c0, c1, c2, c3, c4, c5, c8, c10, c16, c1000, &
                       c10000, c1p5, p33, p5, p25, p125, p001, eps, &
                       eps2, bignum, mwjfnp0s0t0, mwjfnp0s0t1, mwjfnp0s0t2, &
                       mwjfnp0s0t3, mwjfnp0s1t0, mwjfnp0s1t1, mwjfnp0s2t0, &
                       mwjfnp1s0t0, mwjfnp1s0t2, mwjfnp1s1t0, mwjfnp2s0t0, &
                       mwjfnp2s0t2, mwjfdp0s0t0, mwjfdp0s0t1, mwjfdp0s0t2, &
                       mwjfdp0s0t3, mwjfdp0s0t4, mwjfdp0s1t0, mwjfdp0s1t1, &
                       mwjfdp0s1t3, mwjfdp0sqt0, mwjfdp0sqt2, mwjfdp1s0t0, &
                       mwjfdp2s0t3, mwjfdp3s0t1, 0.0_r8 /)

   constants(46) = grav

   this_block = get_block(blocks_clinic(1),1)

   bid = this_block%local_id

!-----------------------------------------------------------------------
!
!  read input namelist and set mixing type
!
!-----------------------------------------------------------------------

   use_gpu_state = .true.
   use_verify_results = .false.

   if (my_task == master_task) then
      open (nml_in, file=nml_filename, status='old',iostat=nml_error)
      if (nml_error /= 0) then
         nml_error = -1
      else
         nml_error =  1
      endif
      do while (nml_error > 0)
         read(nml_in, nml=gpu_mod_nml, iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
   endif

   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
      call exit_POP(sigAbort, &
                    'ERROR reading gpu_mod_nml')
   endif

   if (my_task == master_task) then
      write(stdout,blank_fmt)
      write(stdout,ndelim_fmt)
      write(stdout,'(a23)') 'GPU Accelerator options'
      write(stdout,blank_fmt)
      write(stdout,delim_fmt)

      if (use_gpu_state) then
         write(stdout,'(a43)') ' GPU usage for density computations enabled'
      else
         write(stdout,'(a44)') ' GPU usage for density computations disabled'
      endif

      if (use_verify_results) then
         write(stdout,'(a39)') ' GPU results verfication by CPU enabled'
      else
         write(stdout,'(a40)') ' GPU results verfication by CPU disabled'
      endif

   endif

   call broadcast_scalar(use_gpu_state, master_task)


!-----------------------------------------------------------------------
!
! Proceed to initialize gpu, if necessary
!-----------------------------------------------------------------------
  if (use_gpu_state) then

  !-----------------------------------------------------------------------
  !
  ! Initialize CUDA
  !
  !-----------------------------------------------------------------------

    call cuda_init()


  !-----------------------------------------------------------------------
  !
  !  allocate arrays
  !
  !-----------------------------------------------------------------------

    call my_cudaMallocHost(cptr, (nx_block*ny_block*km*nt*3*max_blocks_clinic))
    call c_f_pointer(cptr, TRACER, (/ nx_block,ny_block,km,nt,3,max_blocks_clinic /))

    call my_cudaMallocHost(cptr, (nx_block*ny_block*km*3*max_blocks_clinic))
    call c_f_pointer(cptr, RHO, (/ nx_block,ny_block,km,3,max_blocks_clinic /))

    call my_cudaMallocHost(cptr, (nx_block*ny_block*km))
    call c_f_pointer(cptr, RHOP, (/ nx_block,ny_block,km /))

!       real (r8), dimension(nx_block,ny_block,km) :: &
!      DBLOC,      &! buoyancy difference between adjacent levels
!      DBSFC,      &! buoyancy difference between level and surface

    call my_cudaMallocHost(cptr, (nx_block*ny_block*km))
    call c_f_pointer(cptr, DBLOC, (/ nx_block,ny_block,km /))

    call my_cudaMallocHost(cptr, (nx_block*ny_block*km))
    call c_f_pointer(cptr, DBSFC, (/ nx_block,ny_block,km /))

      !allocate (VDC(nx_block,ny_block,0:km+1,2,nblocks_clinic), &
      !          VVC(nx_block,ny_block,km,      nblocks_clinic))
    call my_cudaMallocHost(cptr, (nx_block*ny_block*(km+2)*2*max_blocks_clinic))
    call c_f_pointer(cptr, VDC, (nx_block,ny_block,0:km+1,2,max_blocks_clinic))

    call my_cudaMallocHost(cptr, (nx_block*ny_block*km*max_blocks_clinic))
    call c_f_pointer(cptr, VVC, (nx_block,ny_block,km,max_blocks_clinic))


    ! arrays used for correctness checks
    if (use_verify_results) then
      allocate(RHOREF(nx_block,ny_block,km), &
               DBLOCREF(nx_block,ny_block,km), &
               DBSFCREF(nx_block,ny_block,km), &
               VDCREF(nx_block,ny_block,0:km+1,2,max_blocks_clinic), &
               VVCREF(nx_block,ny_block,km,max_blocks_clinic))

    endif

  !-----------------------------------------------------------------------
  !
  !  set up coefficients (such as state() constants)
  !  subroutine init_state_coeffs in state_mod must be called first
  !
  !-----------------------------------------------------------------------

    ! it is important that state_mod has already been initialized
    call cuda_state_initialize(constants, pressz, tmin, tmax, smin, smax, my_task, KMT(:,:,bid))

    !write(stdout, *) ' grav= ', constants(46)

  else
  !-----------------------------------------------------------------------
  !
  ! if (use_gpu_mod == .false.)
  ! allocate TRACER in host memory the normal way
  !
  !-----------------------------------------------------------------------
    allocate(TRACER(nx_block,ny_block,km,nt,3,max_blocks_clinic))
    allocate(RHO(nx_block,ny_block,km,3,max_blocks_clinic))
    allocate(DBLOC(nx_block,ny_block,km))
    allocate(DBSFC(nx_block,ny_block,km))

  endif ! use_gpu_state

!-----------------------------------------------------------------------
!EOC

 end subroutine init_gpu_mod


!cudaDeviceSynchronize()
 subroutine gpumod_devsync

   call devsync

 end subroutine gpumod_devsync


!***********************************************************************
!BOP
! !IROUTINE: state
! !INTERFACE:

 subroutine gpumod_mwjf_state(TEMP, SALT, start_k, end_k, &
                         RHOOUT, DRHODT, DRHODS)

! !DESCRIPTION:
!  GPU Accelerated version of MWJF state assumes
!  state_itype == state_type_mwjf .and.
!  state_range_iopt == state_range_enforce
!
!  Returns the density of water at level k from equation of state
!  $\rho = \rho(d,\theta,S)$ where $d$ is depth, $\theta$ is
!  potential temperature, and $S$ is salinity.
!
!  This routine also computes derivatives of density with respect
!  to temperature and salinity at level k from equation of state
!  if requested (ie the optional arguments are present).
!
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      start_k,                    &! loop start (including) start index 1
      end_k                        ! loop end (including)

   real (r8), dimension(nx_block,ny_block,km), intent(in) :: &
      TEMP,             &! temperature at level k
      SALT               ! salinity    at level k

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,km), optional, intent(out) :: &
      RHOOUT,  &! perturbation density of water
      DRHODT,  &! derivative of density with respect to temperature
      DRHODS    ! derivative of density with respect to salinity

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables:
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      n_outputs = 1

!-----------------------------------------------------------------------
!
!  first check for requested functionality
!
!-----------------------------------------------------------------------
   n_outputs = 1

   if (.not. present(RHOOUT)) then
      ! throw an this is currently not support error
   endif


   if (present(DRHODT) .and. present(DRHODS)) then
      n_outputs = 3
   endif

!-----------------------------------------------------------------------
!
!  call C wrapper function for CUDA kernel
!
!-----------------------------------------------------------------------

   call mwjf_state_gpu(TEMP, SALT, RHOOUT, DRHODT, DRHODS, n_outputs, start_k, end_k)


 end subroutine gpumod_mwjf_state


 subroutine gpumod_mwjf_statePD(TEMP, SALT, start_k, end_k, RHOOUT)

   integer (int_kind), intent(in) :: &
      start_k,                    &! loop start (including) start index 1
      end_k                        ! loop end (including)

   real (r8), dimension(nx_block,ny_block,km), intent(in) :: &
      TEMP,             &! temperature at level k
      SALT               ! salinity    at level k

   real (r8), dimension(nx_block,ny_block,km), intent(out) :: &
      RHOOUT

   call mwjf_statepd_gpu(TEMP, SALT, RHOOUT, start_k, end_k)

 end subroutine gpumod_mwjf_statePD


 subroutine gpumod_compare(A, B, n, var_name)
    real (r8), dimension(nx_block,ny_block,KM), intent(in) :: &
      A,             &! array 1
      B               ! array 2

    integer (int_kind), intent(in) :: &
      n             ! number of elements

    integer (int_kind), intent(in), optional :: &
      var_name

   if (present(var_name)) then
     call gpu_compare(A, B, n, var_name)
   else
     call gpu_compare(A, B, n, 0)
   endif

 end subroutine


 subroutine gpumod_buoydiff(DBLOC, DBSFC, TRCR, this_block)

! !DESCRIPTION:
!  This routine calculates the buoyancy differences at model levels.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,km,2), intent(in) :: &
      TRCR              ! tracers at current time

   type (block), intent(in) :: &
      this_block          ! block information for current block

! !OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,km), intent(out) :: &
      DBLOC,         &! buoyancy difference between adjacent levels
      DBSFC           ! buoyancy difference between level and surface

!    if (my_task == master_task) then
!      write(stdout, *) ' tmin= ', tmin
!      write(stdout, *) ' tmax= ', tmax
!      write(stdout, *) ' smin= ', smin
!      write(stdout, *) ' smax= ', smax
!      write(stdout, *) ' pressz= ', pressz
!    endif

   call buoydiff_gpu(DBLOC, DBSFC, TRCR)

 end subroutine gpumod_buoydiff

 subroutine ddmix(VDC, TRCR, this_block)

! !DESCRIPTION:
!  $R_\rho$ dependent interior flux parameterization.
!  Add double-diffusion diffusivities to Ri-mix values at blending
!  interface and below.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,km,nt), intent(in) :: &
      TRCR                ! tracers at current time

   type (block), intent(in) :: &
      this_block          ! block information for current block

! !INPUT/OUTPUT PARAMETERS:

   real (r8), dimension(nx_block,ny_block,0:km+1,2),intent(inout) :: &
      VDC        ! diffusivity for tracer diffusion

!EOP
!BOC

   call ddmix_gpu(VDC, TRCR)

 end subroutine gpumod_ddmix

 end module gpu_mod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
