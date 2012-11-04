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
   use state_mod     ! access to preszz, tmin, tmax, smin, smax, etc
   use domain_size   ! included for use of nx_block,ny_block,km,
   use prognostic    ! include for reference to TRACER
   use iso_c_binding

   implicit none
   private
   save

   ! include the interface to the C code
   #include "gpu_mod.fh"

! !PUBLIC MEMBER FUNCTIONS:

   public :: init_gpu_mod, &
             mwjf_state, &
             gpumod_compare
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

   integer (int_kind) ::  &
      k,                  &! vertical level index
      nu,                 &! i/o unit
      nml_error            ! namelist i/o error flag

   character (char_len) :: &
      vmix_choice          ! input choice for desired parameterization

   character (char_len) :: &
      convection_type      ! input choice for method for convection

   type(c_ptr) :: &
      cptr                  ! ptr used for alllocating arrays

   ! array used for sending constants to the GPU
   ! this prevents issues with compilers rounding double precision constants
   ! differently and thus reduces rounding errors
   real (r8), parameter, dimension(45) :: &
        constants = (/ c0, c1, c2, c3, c4, c5, c8, c10, c16, c1000, &
                       c10000, c1p5, p33, p5, p25, p125, p001, eps, &
                       eps2, bignum, mwjfnp0s0t0, mwjfnp0s0t1, mwjfnp0s0t2, &
                       mwjfnp0s0t3, mwjfnp0s1t0, mwjfnp0s1t1, mwjfnp0s2t0, &
                       mwjfnp1s0t0, mwjfnp1s0t2, mwjfnp1s1t0, mwjfnp2s0t0, &
                       mwjfnp2s0t2, mwjfdp0s0t0, mwjfdp0s0t1, mwjfdp0s0t2, &
                       mwjfdp0s0t3, mwjfdp0s0t4, mwjfdp0s1t0, mwjfdp0s1t1, &
                       mwjfdp0s1t3, mwjfdp0sqt0, mwjfdp0sqt2, mwjfdp1s0t0, &
                       mwjfdp2s0t3, mwjfdp3s0t1 /)

   namelist /gpu_mod_nml/ use_gpu_state

!-----------------------------------------------------------------------
!
!  read input namelist and set mixing type
!
!-----------------------------------------------------------------------

   use_gpu_state = .true.

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

    allocate(RHOREF(nx_block,ny_block,km)) ! used for correctness checks

  !-----------------------------------------------------------------------
  !
  !  set up coefficients (such as state() constants)
  !  subroutine init_state_coeffs in state_mod must be called first
  !
  !-----------------------------------------------------------------------

    ! it is important that state_mod has already been initialized
    call cuda_state_initialize(constants, pressz, tmin, tmax, smin, smax)


  else
  !-----------------------------------------------------------------------
  !
  ! if (use_gpu_mod == .false.)
  ! allocate TRACER in host memory the normal way
  !
  !-----------------------------------------------------------------------
    allocate(TRACER(nx_block,ny_block,km,nt,3,max_blocks_clinic))
    allocate(RHO(nx_block,ny_block,km,3,max_blocks_clinic))


  endif ! use_gpu_state

!-----------------------------------------------------------------------
!EOC

 end subroutine init_gpu_mod


!***********************************************************************
!BOP
! !IROUTINE: state
! !INTERFACE:

 subroutine mwjf_state(TEMP, SALT, start_k, end_k, &
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
      n_outputs = 1;

!-----------------------------------------------------------------------
!
!  first check for requested functionality
!
!-----------------------------------------------------------------------

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


 end subroutine mwjf_state


 subroutine gpumod_compare(A, B, n)
    real (r8), dimension(nx_block,ny_block,km), intent(in) :: &
      A,             &! array 1
      B               ! array 2

    integer (int_kind), intent(in) :: &
      n             ! number of elements

   call gpu_compare(A, B, n)

 end subroutine


 end module gpu_mod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
