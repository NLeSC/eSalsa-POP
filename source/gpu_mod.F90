module gpu_mod
! MODULE: gpu_mod
! DESCRIPTION:
! This is the generic GPU module that contains the
! initialization routines for execution on GPUs
!
! REVISION HISTORY:
! gpu_mod.F90 2015-01-16 19:30:26Z B. van Werkhoven
! USES:

 use kinds_mod
 use io

 use iso_c_binding

 implicit none


 !PUBLIC MEMBER FUNCTIONS:
 public :: init_gpu_mod

 !PUBLIC DATA MEMBERS:
 logical (log_kind), public :: &
   use_gpu ! flag for computing density functions on gpu


 type(c_ptr), public :: &
   cptr ! temporary pointer used for allocating arrays


 contains


 subroutine init_gpu_mod 

   integer (int_kind) :: &
      nu, &! i/o unit
      nml_error ! namelist i/o error flag

   namelist /gpu_mod_nml/ use_gpu

   use_gpu = .true.  ! default setting

   if (my_task == master_task) then
     open (nml_in, file=nml_filename, status='old',iostat=nml_error)
     if (nml_error /= 0) then
       nml_error = -1
     else
       nml_error = 1
     endif
     do while (nml_error > 0)
       read(nml_in, nml=gpu_mod_nml, iostat=nml_error)
     end do
     if (nml_error == 0) then
       close(nml_in)
     endif
   endif

   call broadcast_scalar(nml_error, master_task)

   if (nml_error /= 0) then
     call exit_POP(sigAbort, 'ERROR reading gpu_mod_nml')
   endif

  if (my_task == master_task) then
    write(stdout,blank_fmt)
    write(stdout,ndelim_fmt)
    write(stdout,'(a23)') 'GPU Accelerator options'
    write(stdout,blank_fmt)
    write(stdout,delim_fmt)
    if (use_gpu) then
      write(stdout,'(a35)') ' GPU usage for computations enabled'
    else
      write(stdout,'(a36)') ' GPU usage for computations disabled'
    endif
  endif

  call broadcast_scalar(use_gpu, master_task)


  ! now we know if the gpu should be used or not





end subroutine init_gpu_mod
















end module gpu_mod
