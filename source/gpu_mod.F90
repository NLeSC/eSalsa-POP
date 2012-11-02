!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  module gpu_mod

!BOP
! !MODULE: gpu_mod
! !DESCRIPTION:
!  This module contains the routines for computing several functions
!  on GPU accelerators
!
! !REVISION HISTORY:
!  gpu_mod.F90 2012-10-11 19:30:26Z B. van Werkhoven

! !USES:

   use kinds_mod
   use io
   use exit_mod
   use domain_size ! included for use of nx_block,ny_block,km,
   use prognostic  ! include for reference to TRACER

   implicit none
   private
   save

   ! include the interface to the C code
   #include "gpu_mod.fh"

! !PUBLIC MEMBER FUNCTIONS:

   public :: init_gpu_mod
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

    call my_cudaMallocHost(cptr, (nx_block*ny_block*km*nt*3*max_blocks_clinic));
    call c_f_pointer(cptr, TRACER, (nx_block,ny_block,km,nt,3,max_blocks_clinic))


  !-----------------------------------------------------------------------
  !
  !  set up coefficients (such as state() constants)
  !  subroutine init_state_coeffs in state_mod must be called first
  !
  !-----------------------------------------------------------------------



  endif ! use_gpu_state

!-----------------------------------------------------------------------
!EOC

 end subroutine init_gpu_mod


 end module gpu_mod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
