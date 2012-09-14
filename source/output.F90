!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module output

!BOP
! !MODULE: output
! !DESCRIPTION:
!  Contains necessary routines, variables for large model output
!  files - restart, history, movies, drifter, time average files.
!  This module is primarily a driver for the individual output
!  modules.
!
! !REVISION HISTORY:
!  SVN:$Id: output.F90 12674 2008-10-31 22:21:32Z njn01 $
!
! !USES:

   use POP_KindsMod
   use POP_ErrorMod
   use POP_DomainSizeMod, only : POP_maxBlocksClinic

   use domain, only: distrb_clinic
   use restart, only: write_restart, init_restart
   use history, only: write_history, init_history
   use movie, only: write_movie, init_movie
   use tavg, only: write_tavg, init_tavg
   use floats, only: write_float_ascii, init_floats, write_float_netcdf,  &
                       float_fmt_out_public
   use timers, only: get_timer, timer_start, timer_stop

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: output_driver, &
             init_output


!EOP
!BOC
!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: output_driver
! !INTERFACE:

   subroutine output_driver(errorCode)

! !DESCRIPTION:
!  This is the main driver routine for all large model output routines.
!
! !REVISION HISTORY:
!  same as module

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode         ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   character (POP_charLength) :: &
      restart_type          ! type of restart being written - used
                            ! to pass restart info to tavg routines

   integer (POP_i4), save :: &
      timer_tavg,              &! timer for tavg
      timer_out,               &! timer for tavg
      timer_movie,             &! timer for movies
      timer_rest                ! timer for restart

   logical (POP_Logical), save :: &
      first_call = .true.       ! flag for initializing timers


!-----------------------------------------------------------------------
!
!  if this is the first call to output_driver, start some timers
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   if (first_call) then
      call get_timer(timer_out, 'OUTPUT'     ,POP_maxBlocksClinic,distrb_clinic%nprocs)
      call get_timer(timer_tavg,'OUTPUT TAVG',POP_maxBlocksClinic,distrb_clinic%nprocs)
      call get_timer(timer_rest,'OUTPUT REST',POP_maxBlocksClinic,distrb_clinic%nprocs)
      call get_timer(timer_movie,'OUTPUT MOVIE',POP_maxBlocksClinic,distrb_clinic%nprocs)
      first_call = .false.
   endif

   call timer_start(timer_out)

!-----------------------------------------------------------------------
!
!  write history, movie files - the decision when to write
!  is internal to each routine  
!  write these first so that if I/O fails, no restart is written
!
!-----------------------------------------------------------------------

#if (defined _NOIO) 
! Insufficient memory to write history files on Blue Gene
#else
   call write_history
   call timer_start(timer_movie)
   call write_movie
   call timer_stop(timer_movie)
#endif

!-----------------------------------------------------------------------
!
!  check for restart and write restart if required
!
!-----------------------------------------------------------------------

   call timer_start(timer_rest)

   call write_restart(restart_type, errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'init_output: error in write_restart')
      return
   endif

   call timer_stop (timer_rest)

!-----------------------------------------------------------------------
!
!  write tavg - the decision when to write
!  is internal to routine except for notifying tavg that a 
!  restart must be written
!
!-----------------------------------------------------------------------

   call timer_start(timer_tavg)
#if (defined _NOIO)
! Insufficient memory to restart tavg files on Blue Gene
#else
   call write_tavg(restart_type)
#endif
   call timer_stop (timer_tavg)

!-----------------------------------------------------------------------
!  write floats
!-----------------------------------------------------------------------

   if (float_fmt_out_public == 'ascii') then
      call write_float_ascii(restart_type)
   else
      call write_float_netcdf(restart_type)
   endif

   call timer_stop(timer_out)
!-----------------------------------------------------------------------
!EOC

 end subroutine output_driver

!***********************************************************************
!BOP
! !IROUTINE: init_output
! !INTERFACE:

 subroutine init_output(errorCode)

! !DESCRIPTION:
!  Initializes frequency of output and filenames for
!  various files by calling individual initialization routines
!
! !REVISION HISTORY:
!  same as module

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode         ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  call individual init routines
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   call init_restart(errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'init_output: error in init_restart')
      return
   endif

   call init_history
   call init_movie
   call init_tavg
   call init_floats

!-----------------------------------------------------------------------
!EOC

 end subroutine init_output


!***********************************************************************

 end module output

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
