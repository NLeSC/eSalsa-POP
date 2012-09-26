!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module floats

!BOP
! !MODULE: floats
! !DESCRIPTION:
!
! !REVISION HISTORY:
! CVS:$Id: float.F90,v 1.31 2003/12/23 22:32:16 pwjones Exp $
! CVS:$Name: $

! !USES:

   use kinds_mod
   use blocks
   use distribution
   use domain
   use constants
   use prognostic
   use grid
   use time_management
   use global_reductions
   use gather_floats
   use broadcast
   use io
   use exit_mod
   use timers, only : get_timer, timer_start, timer_stop

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: init_floats, &
             define_float_field, &
             float_move_driver, &
             calc_float_field_values,&
             update_float_buffer, &
             float_global_to_local, &
             float_requested, &
             write_float_ascii, &
             write_float_netcdf, &
             read_float_ascii, &
             read_float_netcdf

! !PUBLIC DATA MEMBERS:

   logical (log_kind), public :: &
      lfloat_on = .false., &! float file output wanted
      lfloat_restart = .false. ! float restart flag

   logical (log_kind), allocatable, dimension(:,:), public :: &
      float_on_this_block

   integer (int_kind), public :: &
      num_floats, &
      float_flag, & ! time flag for writing float files
      angle_bufloc

   real (r4), parameter, public :: &
      float_special_value = -999.0_r4

   character(char_len), public :: float_fmt_out_public ! public version of fmt

!EOP
!BOC
!-----------------------------------------------------------------------
!
! float field descriptor data type and array of such types
!
!-----------------------------------------------------------------------

   type :: float_field_desc
      character(char_len) :: short_name ! short name for field
      character(char_len) :: long_name ! long descriptive name
      character(char_len) :: units ! units
      real (r4) :: missing_value ! value on land pts
      real (r4), dimension(2) :: valid_range ! min/max
      integer (i4) :: buf_loc ! location in buffer
      integer (i4) :: field_loc ! grid location and field
      integer (i4) :: field_type ! type for io, ghost cells
   end type

   integer (int_kind), parameter :: &
      max_avail_float_fields = 200 ! limit on available fields - can
                                     ! be pushed as high as necessary

   integer (int_kind) :: &
      num_avail_float_fields = 0, &! current number of defined fields
      num_requested_float_fields ! number of fields requested

   type (float_field_desc), dimension(max_avail_float_fields) :: &
      avail_float_fields

   real (r4), allocatable, dimension(:,:), public :: &
      float_field_values_gathered, &
      float_ijk_global, &
      float_xyza_gathered

   real (r4), allocatable, dimension(:,:,:), public :: &
      float_field_values, &
      float_ijk, &
      float_xyza, &
      float_weights_Upoint, &
      float_weights_Tpoint

   real (r4), allocatable, dimension(:) :: &
      float_deploy_time

   character(char_len), allocatable, dimension(:) :: &
      float_type

   integer (int_kind), allocatable, dimension(:) :: &
      float_itype

!-----------------------------------------------------------------------
!
! buffers for holding running float variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      float_bufsize ! size of buffer for float fields

   real (r4), dimension(:,:,:,:,:), allocatable :: &
      FLOAT_BUF ! buffer for holding updated float fields

   real (r8), allocatable, dimension(:,:,:,:) :: &
      WVEL ! vertical velocity at U-points

  !-----------------------------------------------------------------------------
  ! define float field id for ANGLE
  !-----------------------------------------------------------------------------

   integer (int_kind) :: &
      float_ANGLE

!-----------------------------------------------------------------------
!
! variables for writing data
!
!-----------------------------------------------------------------------

   integer (i4), public :: &
      float_freq_iopt, &! frequency option for writing float
      float_freq ! frequency of float output

   character (char_len) :: &
      float_infile, & ! filename for restart input
      float_outfile, & ! root filename for float output
      float_fmt_in, & ! format (nc or bin) for reading restart
      float_fmt_out ! format (nc or bin) for writing

!-----------------------------------------------------------------------
!
! scalars
!
!-----------------------------------------------------------------------

   integer (int_kind), public :: &
      calc_float_timer, float_move_timer

   integer (int_kind) :: &
      float_advect_2d = 1, &
      float_advect_3d = 2

   real (r4), parameter :: &
      shallowest_depth_fraction = 0.1_r4

   character (10) :: &
      beg_date ! date on which the current accumulated sum
                     ! was started (not the tavg_start date)
!EOC
!***********************************************************************

 contains

!***********************************************************************
!EOP
! !IROUTINE: init_floats
! !INTERFACE:

 subroutine init_floats

! !DESCRIPTION:
! This routine initializes float options and reads in contents file to
! determine which fields for which the user wants float data.
!
! !REVISION HISTORY:
! same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!
! local variables
!
!-----------------------------------------------------------------------

   integer (i4) :: &
      n,isrch, &! dummy index
      iblock, &! local block index
      loc, &! location of field in buffer
      nu, &! unit for contents input file
      cindex, &! character index for manipulating strings
      nml_error, &! namelist i/o error flag
      contents_error ! error flag for contents file read

   character (char_len) :: &
      float_freq_opt, &! choice for frequency of float output
      float_init_opt, &! choice for initialization option
      float_init_file, &! name of init file
      float_init_file_fmt, &! format of init file
      float_contents, &! filename for choosing fields for output
      char_temp ! temporary for manipulating fields

   character (36), parameter :: &
      freq_fmt = "('float diagnostics every ',i6,a8)"

   character (47), parameter :: &
      start_fmt = "('float sums updated starting at ',a5,i8)"

   namelist /float_nml/ &
      float_freq_opt, float_freq, float_infile, &
      float_outfile, float_contents, &
      float_fmt_in, float_fmt_out, &
      float_init_opt, float_init_file, float_init_file_fmt

!-----------------------------------------------------------------------
!
! read float file output frequency and filenames from namelist
!
!-----------------------------------------------------------------------

   if (my_task == master_task) then
      write(stdout,delim_fmt)
      write(stdout,blank_fmt)
      write(stdout,'(a15)') 'Drifter options'
      write(stdout,blank_fmt)
      write(stdout,delim_fmt)
   endif

   float_freq_opt = 'never'
   float_freq_iopt = freq_opt_never
   float_freq = 100000
   float_infile = 'unknown_float_infile'
   float_outfile = 'float.'
   float_contents = 'unknown_float_contents'
   float_fmt_in = 'ascii'
   float_fmt_out = 'ascii'
   float_init_opt = 'file'
   float_init_file = 'unknown_float_init_file'
   float_init_file_fmt = 'unknown_float_init_file_fmt'

   if (my_task == master_task) then
      open (nml_in, file=nml_filename, status='old',iostat=nml_error)
      if (nml_error /= 0) then
         nml_error = -1
      else
         nml_error = 1
      endif
      do while (nml_error > 0)
         read(nml_in, nml=float_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
   endif

   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
      call exit_POP(sigAbort,'ERROR reading float_nml')
   endif

   if (my_task == master_task) then
      select case (float_freq_opt)
      case ('never')
         float_freq_iopt = freq_opt_never
         write(stdout,'(a20)') 'float diagnostics off'
      case ('nyear')
         float_freq_iopt = freq_opt_nyear
         write(stdout,freq_fmt) float_freq,' years  '
      case ('nmonth')
         float_freq_iopt = freq_opt_nmonth
         write(stdout,freq_fmt) float_freq,' months '
      case ('nday')
         float_freq_iopt = freq_opt_nday
         write(stdout,freq_fmt) float_freq,' days   '
      case ('nhour')
         float_freq_iopt = freq_opt_nhour
         write(stdout,freq_fmt) float_freq,' hours  '
      case ('nsecond')
         float_freq_iopt = freq_opt_nsecond
         write(stdout,freq_fmt) float_freq,' seconds'
      case ('nstep')
         float_freq_iopt = freq_opt_nstep
         write(stdout,freq_fmt) float_freq,' steps  '
      case default
         float_freq_iopt = -1000
      end select

   endif

   call broadcast_scalar(float_freq_iopt, master_task)

   if (float_freq_iopt == -1000) then
      call exit_POP(sigAbort,'unknown option for float file frequency')
   else if (float_freq_iopt /= freq_opt_never) then
      call broadcast_scalar(float_freq, master_task)
      call broadcast_scalar(float_infile, master_task)
      call broadcast_scalar(float_outfile, master_task)
      call broadcast_scalar(float_contents, master_task)
      call broadcast_scalar(float_fmt_in , master_task)
      call broadcast_scalar(float_fmt_out, master_task)
      call broadcast_scalar(float_init_opt, master_task)
      call broadcast_scalar(float_init_file, master_task)
      call broadcast_scalar(float_init_file_fmt, master_task)
      lfloat_on = .true.
   endif

   float_fmt_out_public = float_fmt_out

!-----------------------------------------------------------------------
!
! initialize time flag for writing float files
!
!-----------------------------------------------------------------------

   call init_time_flag('float', float_flag, default=.false., &
                                owner = 'init_floats', &
                                freq_opt = float_freq_iopt, &
                                freq = float_freq)

!-----------------------------------------------------------------------
!
! read contents file to determine which fields to dump
!
!-----------------------------------------------------------------------

   if (float_freq_iopt /= freq_opt_never) then

      call get_timer(calc_float_timer, 'UPDATE_FLOAT', &
                     nblocks_clinic, distrb_clinic%nprocs)
      call get_timer(float_move_timer, 'MOVE_FLOAT', &
                     nblocks_clinic, distrb_clinic%nprocs)
      call get_timer(gather_float_timer, 'GATHER_FLOAT', &
                     1, distrb_clinic%nprocs)

      float_bufsize = 0

!-----------------------------------------------------------------------
! we will always want ANGLE so do it now
!-----------------------------------------------------------------------

      call define_float_field(float_ANGLE,'ANGLE', &
                          long_name='grid angle interpolated to float position', &
                          units='radians', grid_loc='2221')
      char_temp = 'ANGLE'
      call request_float_field(trim(char_temp))
      angle_bufloc = float_bufsize

      call get_unit(nu)

      if (my_task == master_task) then
         open(nu, file=float_contents, status='old')
         read(nu,*) num_requested_float_fields, num_floats
         write(stdout,'(a41)') 'float diagnostics requested for fields:'
      endif

      call broadcast_scalar(num_requested_float_fields, master_task)
      call broadcast_scalar(num_floats, master_task)

      contents_error = 0

      do n=1,num_requested_float_fields
         if (my_task == master_task) then
            read(nu,'(a80)',iostat=contents_error) char_temp
            char_temp = adjustl(char_temp)
            cindex = index(char_temp,' ')
            char_temp(cindex:) = ' '
            write(stdout,*) '  ',trim(char_temp)
         endif

         call broadcast_scalar(contents_error, master_task)
         if (contents_error /= 0) then
            call exit_POP(sigAbort,'error reading float contents')
         endif

         call broadcast_scalar(char_temp, master_task)
         call request_float_field(trim(char_temp))
      end do

!-----------------------------------------------------------------------
!
! allocate float location arrays, then read in global logical
! locations for floats.
!
!-----------------------------------------------------------------------

      allocate( float_field_values_gathered(num_floats,float_bufsize), &
                float_field_values(num_floats,nblocks_clinic,float_bufsize) )
      allocate( float_ijk_global (num_floats,3) &
             , float_ijk (num_floats,nblocks_clinic,3) &
             , float_xyza_gathered(num_floats,4) &
             , float_xyza (num_floats,nblocks_clinic,4) &
             , float_weights_Upoint(num_floats,nblocks_clinic,3) &
             , float_weights_Tpoint(num_floats,nblocks_clinic,3) &
             , float_on_this_block(num_floats,nblocks_clinic) )
      allocate( float_deploy_time(num_floats))
      allocate( float_type(num_floats))
      allocate( float_itype(num_floats))
      float_ijk = c0
      float_ijk_global = c0
      float_xyza = float_special_value
      float_xyza_gathered = c0
      float_weights_Upoint = float_special_value
      float_weights_Tpoint = float_special_value
      float_on_this_block = .false.
      float_deploy_time = c0
      float_type = 'constant-depth'
      float_itype = float_advect_2d

      if (float_init_opt == 'file') then
         if(my_task == master_task) then ! continue reading contents file
            do n = 1, num_floats
               read(nu,*) float_ijk_global(n,1), float_ijk_global(n,2), &
                          float_ijk_global(n,3), float_deploy_time(n), &
                          float_type(n)
            enddo
         endif
      elseif (float_init_opt == 'restart') then
         if (float_fmt_in == 'ascii') then
            call read_float_ascii
         elseif (float_fmt_in == 'nc') then
            call read_float_netcdf
         else
            call exit_POP(sigAbort,'ERROR: unknown float_fmt_in value')
         endif
      else
         call exit_POP(sigAbort,'ERROR: unknown float_init_opt')
      endif

      if (my_task == master_task) close(nu)

      call release_unit(nu)

      !*** allocate and initialize float buffer

      allocate( &
         FLOAT_BUF(nx_block, ny_block, km, nblocks_clinic, float_bufsize) )
      FLOAT_BUF = c0

      allocate( WVEL(nx_block, ny_block, km, nblocks_clinic) )
      WVEL = c0

      call broadcast_array(float_ijk_global, master_task)
      call broadcast_array(float_deploy_time, master_task)
      do n = 1, num_floats
         call broadcast_scalar(float_type(n), master_task) ! no char array bcast
      enddo

      call float_global_to_local

      do n = 1, num_floats
         if (trim(float_type(n)) == 'constant-depth') then
            float_itype(n) = float_advect_2d
         elseif (trim(float_type(n)) == '3d') then
            float_itype(n) = float_advect_3d
         else
            call exit_POP(sigAbort,'ERROR: unknown float_type')
         endif
      enddo

   endif ! if iopt /= never

!-----------------------------------------------------------------------
!EOC

 end subroutine init_floats

!***********************************************************************
!BOP
! !IROUTINE: write_float_ascii
! !INTERFACE:

 subroutine write_float_ascii(restart_type)

! !DESCRIPTION:
! This routine writes requested float fields to a file. The fields are
! normalized by the time interval before writing.
!
! !REVISION HISTORY:
! same as module

! !INPUT PARAMETERS:

   character (POP_charLength), intent(in) :: &
      restart_type ! tells float whether to write restart

!EOP
!BOC
!-----------------------------------------------------------------------
!
! local variables
!
!-----------------------------------------------------------------------

   integer (i4) :: &
      nu, &! i/o unit for output file
      iblock, &! dummy block index
      nfield, &! dummy field index
      m, &! dummy float index
      k, &! dummy depth index
      dump_type, &! dummy restart type index
      loc ! buffer location for field

   character (char_len), dimension(2) :: &
      file_suffix ! suffix to append to float file name

   character (char_len) :: &
      hist_string, &! string containing file history
      float_filename, &! filename for float data
      float_pointer_file ! filename for pointer file containing
                             ! location/name of last restart file

   character (8) :: &
      date_created ! string with (real) date this file created

   character (10) :: &
      time_created ! string with (real) date this file created

   logical (log_kind), dimension(2) :: &
      lfloat_write ! time to write a file

   logical (log_kind) :: &
      lpositions, &! true if gathering up float positions
      already_gathered !

   real (r4), allocatable, dimension(:,:,:) :: &
      float_ijk_temp

!-----------------------------------------------------------------------
!
! is it time to write a file - if yes, create a file suffix
!
!-----------------------------------------------------------------------

   lfloat_write(:) = .false.
   file_suffix(:) = 'none'

   if (lfloat_on) then
      lfloat_write(1) = check_time_flag(float_flag)

      !*** regular float dump
      if (lfloat_write(1)) then
         call create_suffix_float(file_suffix(1))
      endif

      !*** float restart
      if (trim(restart_type) /= 'none') then
         if (.not. lfloat_write(2)) then
            lfloat_write(2) = .true.

            select case (trim(restart_type))
            case('even')
               file_suffix(2) = trim(runid)/&
                                         &/'.even'
            case('odd')
               file_suffix(2) = trim(runid)/&
                                         &/'.odd'
            case('end')
               file_suffix(2) = trim(runid)/&
                                         &/'.end'
            case default
               call create_suffix_float(file_suffix(2))
               file_suffix(2) = trim(file_suffix(2))/&
                                               &/'.restart'
            end select
         endif
      endif
   endif

!-----------------------------------------------------------------------
!
! do the rest only if it is time to do a float dump
!
!-----------------------------------------------------------------------

   already_gathered = .false.

   do dump_type = 1,2

!-----------------------------------------------------------------------
!
! create data file descriptor
!
!-----------------------------------------------------------------------

      call date_and_time(date=date_created, time=time_created)
      hist_string = char_blank
      write(hist_string,'(a26,a8,1x,a10)') &
         'POP FLOAT file created: ',date_created,time_created

      float_filename = trim(float_outfile)/&
                                               &/trim(file_suffix(dump_type))/&
                                                                    &/'.ascii'

!-----------------------------------------------------------------------
!
! Gather up output data from all processors.
!
!-----------------------------------------------------------------------

     if (lfloat_write(dump_type) .and. .not. already_gathered) then

!-----------------------------------------------------------------------
! float_ijk gets overwritten in gather_float_fields, so copy into a temp
!-----------------------------------------------------------------------

         lpositions = .true.
         allocate( float_ijk_temp(num_floats,nblocks_clinic,3) )
         float_ijk_temp = float_ijk
         call gather_float_fields(float_ijk_global, float_ijk_temp, &
                                    float_on_this_block, master_task, &
                                    distrb_clinic, float_special_value, &
                                    lpositions)
         deallocate(float_ijk_temp)

         lpositions = .false.
         call gather_float_fields(float_xyza_gathered, float_xyza, &
                                    float_on_this_block, master_task, &
                                    distrb_clinic, float_special_value, &
                                    lpositions)

         call gather_float_fields(float_field_values_gathered, &
                                    float_field_values, &
                                    float_on_this_block, master_task, &
                                    distrb_clinic, float_special_value, &
                                    lpositions)

         already_gathered = .true.

      endif ! if time to write and not already_gathered

!-----------------------------------------------------------------------
!
! open output file
!
!-----------------------------------------------------------------------

     if (lfloat_write(dump_type)) then

      call get_unit(nu)
      if (my_task == master_task) then

        write(stdout,*) 'Writing float file: ', trim(float_filename)
        open(nu,file=float_filename,status='unknown')
        write(nu,'(i5,i8,i10,f15.5,i8,i4,i5)') &
           num_requested_float_fields, num_floats, &
           nsteps_total,tday,iyear,imonth,iday

        do nfield = 1,num_avail_float_fields ! check all available fields
           loc = avail_float_fields(nfield)%buf_loc ! locate field in buffer
           if (loc > 0) then ! field is actually requested and in buffer
              write(nu,'(a,2x,a,2x,a)')trim(avail_float_fields(nfield)%short_name) &
                               ,trim(avail_float_fields(nfield)%units) &
                               ,trim(avail_float_fields(nfield)%long_name)
           endif
        enddo

        do m = 1, num_floats
          write(nu,'(4f15.6,2x,a)')float_ijk_global(m,1), &
                                   float_ijk_global(m,2), &
                                   float_ijk_global(m,3), &
                                   float_deploy_time(m), &
                                   trim(float_type(m))
        end do

        do m = 1, num_floats
          write(nu,'(4f15.6)')float_xyza_gathered(m,1), &
                              float_xyza_gathered(m,2), &
                              float_xyza_gathered(m,3), &
                              float_xyza_gathered(m,4)
        end do

!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------

        do nfield = 1,num_avail_float_fields ! check all available fields

           loc = avail_float_fields(nfield)%buf_loc ! locate field in buffer

           if (loc > 0) then ! field is actually requested and in buffer


              write(nu,'(a,2x,a,2x,a)') trim(avail_float_fields(nfield)%short_name), &
                               trim(avail_float_fields(nfield)%units), &
                               trim(avail_float_fields(nfield)%long_name)

              do m = 1, num_floats
                 write(nu,'(i10, e15.7)')m, float_field_values_gathered(m,loc)
              enddo

           endif
        end do

        write(nu,'(a)') trim(hist_string)

        close(nu)

      endif ! master_task

      call release_unit(nu)

!-----------------------------------------------------------------------
!
! if pointer files are used, write float filenames to pointer file
! do this only for float restarts - not float dumps
!
!-----------------------------------------------------------------------

      if (luse_pointer_files .and. lfloat_write(2)) then
         call get_unit(nu)
         if (my_task == master_task) then
            float_pointer_file = trim(pointer_filename)/&
                                                       &/'.float'

            open(nu,file=float_pointer_file,form='formatted', &
                    status='unknown')
            write(nu,'(a)') trim(float_filename)
            close(nu)
         endif
         call release_unit(nu)
      endif

   endif ! lfloat_write

   enddo ! dump_type

!-----------------------------------------------------------------------
!EOC

 end subroutine write_float_ascii

!***********************************************************************
!BOP
! !IROUTINE: read_float_ascii
! !INTERFACE:

 subroutine read_float_ascii

! !DESCRIPTION:
! This routine reads an ascii float restart dump.
!
! !REVISION HISTORY:
! same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!
! local variables
!
!-----------------------------------------------------------------------

   integer (i4) :: &
     nu, & ! i/o unit
     iblock, & ! dummy block index
     n,i,m,k, & ! dummy for indexing character string
     in_fields, & ! num of fields in restart file
     nfield, & ! dummy field counter
     hdr_error, & ! error file for reading restart hdr
     in_nsteps_total, & ! nsteps_total according to float file
     in_iyear, & ! iyear according to float file
     in_imonth, & ! imonth according to float file
     in_iday, & ! iday according to float file
     in_num_float_fields, & ! num_requested_float_fields according to float file
     in_num_floats, & ! num_floats according to float file
     loc ! buffer location

   real (r4) :: &
     in_tday, & ! tday according to float file
     in1,in2,in3,in4 ! placeholders for non-saved values read from float file

   character (char_len) :: &
     in_field_name , & ! temporary
     char_temp, & ! for string manipulation
     float_pointer_file ! filename for pointer file containing
                           ! location/name of last restart file

   logical (log_kind) :: restart_error

!-----------------------------------------------------------------------
!
! if pointer files are used, pointer file and must be read to get
! actual filenames
!
!-----------------------------------------------------------------------

   call get_unit(nu)

   if (luse_pointer_files) then

      if (my_task == master_task) then
         float_pointer_file = char_blank
         float_pointer_file = trim(pointer_filename)/&
                                                   &/'.float'
         write(stdout,*) 'Reading pointer file: ', &
                         trim(float_pointer_file)
         open(nu, file=trim(float_pointer_file), form='formatted', &
                  status='old')
         read(nu,'(a)') float_infile
         close(nu)
      endif
      call broadcast_scalar(float_infile, master_task)

   endif

   call release_unit(nu)

!-----------------------------------------------------------------------
!
! open input file
!
! check for consistency of restart file with float_contents file
!
!-----------------------------------------------------------------------

   call get_unit(nu)

   restart_error = .false.

   if (my_task == master_task) then

      open(nu, file = trim(float_infile), status = 'old')
      read(nu,*) in_num_float_fields, in_num_floats, &
                 in_nsteps_total

      !*** check nsteps total for validity
      if (in_nsteps_total /= nsteps_total) then
         write(stdout,'(i6,a32,i6,a35)') &
            in_nsteps_total,' nsteps_total in float restart; ', &
            nsteps_total, ' nsteps_total in current simulation'
         restart_error = .true.
         char_temp = 'FLOAT:restart file has wrong time step?'
      endif

      !*** check number of requested float fields for validity
      if (in_num_float_fields /= num_requested_float_fields) then
         write(stdout,'(i6,a44,i6,a47)') &
            in_num_float_fields, &
            'requested float fields in float restart; ', &
            num_requested_float_fields, &
            ' requested float fields in current simulation'
         restart_error = .true.
         char_temp = 'FLOAT:restart file has wrong number of fields'
      endif

      !*** check number of floats for validity
      if (in_num_floats /= num_floats) then
         write(stdout,'(i6,a35,i6,a40)') &
            in_num_floats,' floats in float restart; ', &
            num_floats, ' floats in current simulation'
         restart_error = .true.
         char_temp = 'FLOAT:restart file has wrong number of floats'
      endif

   endif ! master_task

   call broadcast_scalar(restart_error, master_task)
   call broadcast_scalar(char_temp , master_task)

   if (restart_error) &
      call exit_POP(sigAbort,trim(char_temp))

!-----------------------------------------------------------------------
! Now read in field names, float locations, etc.
!-----------------------------------------------------------------------

   if (my_task == master_task) then

      write(stdout,*)' the following fields are in the float restart file:'
      do nfield = 1, num_requested_float_fields
         read(nu,'(a)') in_field_name
         write(stdout,*)'   ',trim(in_field_name)
      enddo
!maltrud read in ANGLE
      read(nu,'(a)') in_field_name
      write(stdout,*)'   ',trim(in_field_name)
      write(stdout,*) &
        ' WARNING: float restart fields ARE NOT compared with float_contents file'

      do n = 1, num_floats ! float_ijk
         read(nu,*) float_ijk_global(n,1), float_ijk_global(n,2), &
                    float_ijk_global(n,3), float_deploy_time(n), &
                    float_type(n)
      enddo

      close(nu)

   endif ! master_task

   call release_unit(nu)

!-----------------------------------------------------------------------
!EOC

 end subroutine read_float_ascii

!***********************************************************************
!BOP
! !IROUTINE: calc_float_field_values
! !INTERFACE:

 subroutine calc_float_field_values

! !DESCRIPTION:
! This routine updates a float field. If the time average of the
! field is requested, it updates a time sum of a field by
! multiplying by the time step and updateng the sum into the
! float buffer array. If the min or max of a field is requested, it
! checks the current value and replaces the min, max if the current
! value is less than or greater than the stored value.
!
! !REVISION HISTORY:
! same as module

! !INPUT PARAMETERS:

!EOP
!BOC
!-----------------------------------------------------------------------
!
! local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      bufloc, &! location of field in float buffer
      field_grid_loc, &! location of requested field (T or U point)
      m, iblock, nfield ! loop counters

   real (r4) :: &
      interp_value ! value of FLOAT_BUF interpolated to float location

   real (r4), dimension(2) :: &
      interp_value_vec ! value of vector FLOAT_BUF interpolated to float

   logical (log_kind) :: &
      vector_field

!-----------------------------------------------------------------------
!
! get buffer location and field info from avail_float_field array
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
! get real space values for locations of floats
! Note that calc_float_xyz_position also calls calc_float_weights
!
!-----------------------------------------------------------------------

   do iblock = 1, nblocks_clinic

      call timer_start(calc_float_timer, block_id=iblock)

      call calc_float_xyz_position(iblock)

      nfield = 1
      do while (nfield <= num_avail_float_fields) ! check all available fields

         bufloc = avail_float_fields(nfield)%buf_loc
! if (bufloc <= 0) &
! call exit_POP(sigAbort, &
! 'float: attempt to update bad float field')
         if (bufloc > 0) then

         field_grid_loc = avail_float_fields(nfield)%field_loc

!-----------------------------------------------------------------------
!
! loop over floats and only interpolate if the float is
! on this block and already deployed.
!
!-----------------------------------------------------------------------
         vector_field = .false.
         if (trim(avail_float_fields(nfield )%short_name) =='UVEL'.and.&
             trim(avail_float_fields(nfield+1)%short_name) =='VVEL') &
             vector_field = .true.

         do m = 1, num_floats

            if (float_on_this_block(m,iblock) &
                .and. tday >= float_deploy_time(m)) then

               if (field_grid_loc == field_loc_center) then
                  call float_field_value_T(interp_value, iblock, m, bufloc)
               else if (field_grid_loc == field_loc_NEcorner) then
                  if (vector_field) then
                     call float_field_value_U_vec(interp_value_vec, iblock, m, bufloc)
                  else
                     call float_field_value_U (interp_value, iblock, m, bufloc)
                  endif
               else
! undefined field_grid_loc
               endif

!-----------------------------------------------------------------------
!
! update the field into the local block values array
!
!-----------------------------------------------------------------------

               if(vector_field) then
                  float_field_values(m,iblock,bufloc ) = interp_value_vec(1)
                  float_field_values(m,iblock,bufloc+1) = interp_value_vec(2)
               else
                  float_field_values(m,iblock,bufloc) = interp_value
               endif

            endif ! float on this block?

         enddo ! loop over floats
         if (vector_field) nfield = nfield + 1
         endif ! bufloc>0
         nfield = nfield + 1
      enddo ! loop over available fields

      call timer_stop(calc_float_timer, block_id=iblock)

   enddo ! loop over blocks

!-----------------------------------------------------------------------
!EOC

 end subroutine calc_float_field_values

!***********************************************************************
!***********************************************************************
!BOP
! !IROUTINE: calc_float_xyz_position
! !INTERFACE:

 subroutine calc_float_xyz_position(iblock)

! !DESCRIPTION:
! This routine updates a float field. If the time average of the
! field is requested, it updates a time sum of a field by
! multiplying by the time step and updateng the sum into the
! float buffer array. If the min or max of a field is requested, it
! checks the current value and replaces the min, max if the current
! value is less than or greater than the stored value.
!
! !REVISION HISTORY:
! same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      iblock ! local block address (in baroclinic distribution)

!EOP
!BOC
!-----------------------------------------------------------------------
!
! local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      bufloc, &! location of field in float buffer
      field_grid_loc, &! location of requested field (T or U point)
      il, jl, kl, ipl, jpl, kpl, &! indices
      m ! loop counter

   real (r4) :: &
      radian_r4, &! r4 version of radian
      interp_value ! value of ARRAY interpolated to float location

   real (r8) :: &
      dx,dy,dk,dth1,dth2,dth3,dph1,dph2,dph3 ! geometrical factors

   type (block) :: &
     this_block ! block info for current block

!-----------------------------------------------------------------------
!
! get buffer location and field info from avail_float_field array
!
!-----------------------------------------------------------------------

   call timer_start(calc_float_timer, block_id=iblock)

   radian_r4 = radian

!-----------------------------------------------------------------------
!
! need to calculate weights since last call was corrector step.
!
!-----------------------------------------------------------------------

   this_block = get_block(blocks_clinic(iblock),iblock)

   call calc_float_weights(float_ijk, iblock, this_block)

!-----------------------------------------------------------------------
!
! loop over floats and only interpolate if the float is on this block.
!
!-----------------------------------------------------------------------

   do m = 1, num_floats

      if (float_on_this_block(m,iblock)) then ! even if not deployed yet

!-----------------------------------------------------------------------
! Find lat-lon location of each float.
!-----------------------------------------------------------------------

         il = int(float_ijk(m,iblock,1))
         jl = int(float_ijk(m,iblock,2))
         kl = int(float_ijk(m,iblock,3))

         ipl = il + 1
         jpl = jl + 1
         kpl = kl + 1

         kpl = min(kpl,km)

         dk = float_ijk(m,iblock,3) - kl

         dx = float_weights_Upoint(m,iblock,1)
         dy = float_weights_Upoint(m,iblock,2)

         dth1 = ULAT(ipl,jl ,iblock) - ULAT(il ,jl ,iblock)
         dth2 = ULAT(il ,jpl,iblock) - ULAT(il ,jl ,iblock)
         dth3 = ULAT(ipl,jpl,iblock) - ULAT(ipl,jl ,iblock) - dth2

         dph1 = ULON(ipl,jl ,iblock) - ULON(il ,jl ,iblock)
         dph2 = ULON(il ,jpl,iblock) - ULON(il ,jl ,iblock)
         dph3 = ULON(ipl,jpl,iblock) - ULON(ipl,jl ,iblock)

         if (dph1 > c3*pih) dph1 = dph1 - pi2
         if (dph2 > c3*pih) dph2 = dph2 - pi2
         if (dph3 > c3*pih) dph3 = dph3 - pi2
         if (dph1 < -c3*pih) dph1 = dph1 + pi2
         if (dph2 < -c3*pih) dph2 = dph2 + pi2
         if (dph3 < -c3*pih) dph3 = dph3 + pi2

         dph3 = dph3 - dph2

         float_xyza(m,iblock,1) = ULON(il,jl,iblock) &
                                   + dph1*dx + dph2*dy + dph3*dx*dy
         float_xyza(m,iblock,2) = ULAT (il,jl,iblock) &
                                   + dth1*dx + dth2*dy + dth3*dx*dy

!-----------------------------------------------------------------------
! Convert to degrees
!-----------------------------------------------------------------------

         float_xyza(m,iblock,1) = float_xyza(m,iblock,1)*radian_r4
         if(float_xyza(m,iblock,1) > 360.0_r4) &
            float_xyza(m,iblock,1) = float_xyza(m,iblock,1) - 360.0_r4
         if(float_xyza(m,iblock,1) < 0.0_r4) &
            float_xyza(m,iblock,1) = float_xyza(m,iblock,1) + 360.0_r4
         float_xyza(m,iblock,2) = float_xyza(m,iblock,2)*radian_r4

!-----------------------------------------------------------------------
! find vertical location in meters
!-----------------------------------------------------------------------
         if (kl == 0) then
            float_xyza(m,iblock,3) = dk*dz(kpl)
         else
            float_xyza(m,iblock,3) = dk*dz(kpl) + zw(kl)
         endif
         float_xyza(m,iblock,3) = float_xyza(m,iblock,3)*0.01_r4 ! cm to m

!-----------------------------------------------------------------------
! Interpolate ANGLE to float location.
!-----------------------------------------------------------------------

         call float_field_value_U(interp_value, iblock, m, angle_bufloc)
         float_xyza(m,iblock,4) = interp_value

      endif ! float on this block?

   enddo ! loop over floats

   call timer_stop(calc_float_timer, block_id=iblock)

!-----------------------------------------------------------------------
!EOC

 end subroutine calc_float_xyz_position

!***********************************************************************
!BOP
! !IROUTINE: define_float_field
! !INTERFACE:

 subroutine define_float_field(id, short_name, &
                                  long_name, units, &
                                  grid_loc, missing_value, valid_range, &
                                  field_loc, field_type)

! !DESCRIPTION:
! Initializes description of an available field and returns location
! in the available fields array for use in later float calls.
!
! !REVISION HISTORY:
! same as module

! !OUTPUT PARAMETERS:

   integer (int_kind), intent(out) :: &
      id ! location in avail_fields array for use in
                        ! later float routines

! !INPUT PARAMETERS:

   character(*), intent(in) :: &
      short_name ! short name for field

   integer (i4), intent(in), optional :: &
      field_loc, &! location in grid
      field_type ! type of field (scalar, vector, angle)

   character(*), intent(in), optional :: &
      long_name, &! long descriptive name for field
      units ! physical units for field

   character(4), intent(in), optional :: &
      grid_loc ! location in grid (in 4-digit code)

   real (r4), intent(in), optional :: &
      missing_value ! value on land pts

   real (r4), dimension(2), intent(in), optional :: &
      valid_range ! min/max

!EOP
!BOC
!-----------------------------------------------------------------------
!
! increment the number of defined fields and make sure it does not
! exceed the maximum
! return the id as the current number
!
!-----------------------------------------------------------------------

   num_avail_float_fields = num_avail_float_fields + 1
   if (num_avail_float_fields > max_avail_float_fields) then
      call exit_POP(sigAbort, &
                    'float: defined float fields > max allowed')
   endif

   id = num_avail_float_fields

!-----------------------------------------------------------------------
!
! now fill the field descriptor
!
!-----------------------------------------------------------------------

   avail_float_fields(id)%short_name = short_name
   avail_float_fields(id)%buf_loc = 0 ! will be reset later

   if (present(long_name)) then
      avail_float_fields(id)%long_name = long_name
   else
      avail_float_fields(id)%long_name = char_blank
   endif

   if (present(units)) then
      avail_float_fields(id)%units = units
   else
      avail_float_fields(id)%units = char_blank
   endif

   if (present(missing_value)) then
      avail_float_fields(id)%missing_value = missing_value
   else
      avail_float_fields(id)%missing_value = undefined
   endif

   if (present(valid_range)) then
      avail_float_fields(id)%valid_range = valid_range
   else
      avail_float_fields(id)%valid_range = undefined
   endif

   !*** set field location, field type used by i/o, ghost cell update
   !*** and other communication routines. because ghost cells for float
   !*** fields are not typically used, the default is field_xxx_noupdate

   if (present(field_loc)) then
      avail_float_fields(id)%field_loc = field_loc
   else
      !*** try to decode field location from grid_loc
      if (grid_loc(2:2) == '1' .and. grid_loc(3:3) == '1') then
         avail_float_fields(id)%field_loc = field_loc_center
      else if (grid_loc(2:2) == '2' .and. grid_loc(3:3) == '2') then
         avail_float_fields(id)%field_loc = field_loc_NEcorner
      else if (grid_loc(2:2) == '1' .and. grid_loc(3:3) == '2') then
         avail_float_fields(id)%field_loc = field_loc_Nface
      else if (grid_loc(2:2) == '2' .and. grid_loc(3:3) == '1') then
         avail_float_fields(id)%field_loc = field_loc_Eface
      else
         avail_float_fields(id)%field_loc = field_loc_noupdate
      endif
   endif

   if (present(field_type)) then
      avail_float_fields(id)%field_type = field_type
   else
      avail_float_fields(id)%field_type = field_type_noupdate
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine define_float_field

!***********************************************************************
!BOP
! !IROUTINE: request_float_field
! !INTERFACE:

 subroutine request_float_field(short_name)

! !DESCRIPTION:
! This field marks an available field as requested and computes
! the location in the float buffer array.
!
! !REVISION HISTORY:
! same as module

! !INPUT PARAMETERS:

   character (*), intent(in) :: &
      short_name ! the short name of the field

!EOP
!BOC
!-----------------------------------------------------------------------
!
! local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      n, &! loop index
      id ! location of field in avail_fields array

!-----------------------------------------------------------------------
!
! search for field with same name
!
!-----------------------------------------------------------------------

   id = 0
   srch_loop: do n=1,num_avail_float_fields
      if (trim(avail_float_fields(n)%short_name) == short_name) then
         id = n
         exit srch_loop
      endif
   end do srch_loop

   if (id == 0) then
      if (my_task == master_task) write(stdout,*) 'Requested ', &
                                                  trim(short_name)
      call exit_POP(sigAbort,'float: requested field unknown')
   endif

!-----------------------------------------------------------------------
!
! set the position in the buffer and advance the buffer position
! for the next field
!
!-----------------------------------------------------------------------

   float_bufsize = float_bufsize + 1
   avail_float_fields(id)%buf_loc = float_bufsize

!-----------------------------------------------------------------------
!
! if float is on, but not started yet, set the buf_loc to -buf_loc
! to show that it is requested, but will not return true for
! requested_float_field
!
!-----------------------------------------------------------------------

   if (.not. lfloat_on) then
      avail_float_fields(id)%buf_loc = -avail_float_fields(id)%buf_loc
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine request_float_field

!***********************************************************************
!BOP
! !IROUTINE: float_requested
! !INTERFACE:

 function float_requested(id)

! !DESCRIPTION:
! This function determines whether an available (defined) float field
! has been requested by a user (through the input contents file) and
! returns true if it has. Note that if floats have been turned off,
! the function will always return false.
!
! !REVISION HISTORY:
! same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      id ! id returned by the define function which
                           ! gives the location of the field

! !OUTPUT PARAMETERS:

   logical (log_kind) :: &
      float_requested ! result of checking whether the field has
                         ! been requested

!EOP
!BOC
!-----------------------------------------------------------------------
!
! check the buffer location - if zero, the field has not been
! requested
!
!-----------------------------------------------------------------------

   if (id < 1 .or. id > num_avail_float_fields) then
      call exit_POP(sigAbort,'float_requested: invalid float id')
   endif

   if (avail_float_fields(id)%buf_loc > 0) then
      float_requested = .true.
   else
      float_requested = .false.
   endif

!-----------------------------------------------------------------------
!EOC

 end function float_requested

!***********************************************************************
!BOP
! !IROUTINE: create_suffix_float
! !INTERFACE:

 subroutine create_suffix_float(file_suffix)

! !DESCRIPTION:
! Creates a suffix to append to output filename based on frequency
! option and averaging interval.
!
! !REVISION HISTORY:
! same as module

! !OUTPUT PARAMETERS:

   character (char_len), intent(out) :: &
      file_suffix ! suffix to append to root filename

!EOP
!BOC
!-----------------------------------------------------------------------
!
! local variable
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      cindx1, cindx2, &! indices into character strings
      len_date ! length of date string

   character (char_len) :: &
      char_temp ! temp character space (for removing spaces)

   character (10) :: &
      cstep_beg, &! beginning step of this particular average
      cstep_end, &! ending step of this particular average
      cdate ! character string with yyyymmdd and optional
                      ! separator (eg yyyy-mm-dd)

   character (4) :: &
      cyear_beg, &! beginning year of this particular average
      cyear_end ! end year of this particular average

   character (2) :: &
      cmonth_beg, &! beginning month of this particular average
      cmonth_end, &! end month of this particular average
      cday_beg, &! beginning day of this particular average
      cday_end ! end day of this particular average

!-----------------------------------------------------------------------
!
! start suffix with runid
!
!-----------------------------------------------------------------------

   file_suffix = char_blank
   cindx2 = len_trim(runid) + 1
   file_suffix(1:cindx2) = trim(runid)/&
                                       &/'.'
   cindx1 = cindx2 + 1

!-----------------------------------------------------------------------
!
! extract beginning year, month, day or time step from beg_date
! and determine end date
!
!-----------------------------------------------------------------------

   !***
   !*** use step numbers if float freq option is nstep
   !***

   write(cstep_end,'(i10)') nsteps_total + 1000000000
   cdate = adjustl(cstep_end)
   cstep_end = trim(cdate)

   call time_stamp('last', 'ymd', date_string=cdate) ! last date

   if (date_separator == ' ') then ! no date separator
      cyear_end = cdate(1:4)
      cmonth_end = cdate(5:6)
      cday_end = cdate(7:8)
   else
      cyear_end = cdate(1:4)
      cmonth_end = cdate(6:7)
      cday_end = cdate(9:10)
   endif

!-----------------------------------------------------------------------
!
! create time portion of suffix based on frequency option
! note that concatenation operator split across lines to avoid
! problems with some cpp preprocessors
!
!-----------------------------------------------------------------------

   select case (float_freq_iopt)
   case (freq_opt_nyear, freq_opt_nmonth, freq_opt_nday)
      cindx2 = cindx1 + 7
      file_suffix(cindx1:cindx2) = cyear_end/&
                                 &/cmonth_end/&
                                 &/cday_end

   case (freq_opt_nstep)
! cindx2 = cindx1 + len_trim(cstep_end) - 1
      cindx2 = cindx1 + len_trim(cstep_end) - 2
      file_suffix(cindx1:cindx2) = cstep_end(2:10)

   case default ! use nstep for other options
      cindx2 = cindx1 + len_trim(cstep_end) - 1
      file_suffix(cindx1:cindx2) = trim(cstep_end)

   end select

!-----------------------------------------------------------------------
!EOC

 end subroutine create_suffix_float

!***********************************************************************
!BOP
! !IROUTINE: calc_float_weights
! !INTERFACE:

 subroutine calc_float_weights(float_local_ijk, iblock, this_block)

! !DESCRIPTION:
! Calculates U and T point bilinear interpolation weights from nearest neighbors to
! float locations.
!
! !REVISION HISTORY:
! same as module

! !OUTPUT PARAMETERS:

!EOP
!BOC

!-----------------------------------------------------------------------
!
! arguments
!
!-----------------------------------------------------------------------

   real (r4) :: float_local_ijk(num_floats,nblocks_clinic,3)

   type (block) :: &
     this_block ! block info for current block

!-----------------------------------------------------------------------
!
! local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: m,i,j,ip,jp,il,ipl,jl,jpl,iblock,kl,kp

   real (r8) :: dx,dy,dth1,dth2,dth3,dph1,dph2,dph3,di,dj,dk
   real (r4) imin,imax,jmin,jmax,interp_value

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

   this_block = get_block(blocks_clinic(iblock),iblock)

!-----------------------------------------------------------------------
! Cycle over all floats.
! Note that we calculate weights even if the float is not deployed
! in case we are just about to output values.
!-----------------------------------------------------------------------

   do m = 1,num_floats

      if (float_on_this_block(m,iblock)) then

         il = int(float_local_ijk(m,iblock,1))
         jl = int(float_local_ijk(m,iblock,2))

         di = float_local_ijk(m,iblock,1) - il
         dj = float_local_ijk(m,iblock,2) - jl

         kl = int(float_local_ijk(m,iblock,3))
         kl = min(max(kl,0),km) ! kl can be zero
         dk = float_local_ijk(m,iblock,3) - kl

!-----------------------------------------------------------------------
! store horizontal interpolation weights for U points and vertical
! weights for (cell bottom) W.
!-----------------------------------------------------------------------

         float_weights_Upoint(m,iblock,1) = di
         float_weights_Upoint(m,iblock,2) = dj
         float_weights_Upoint(m,iblock,3) = dk

!-----------------------------------------------------------------------
! horizontal weights for T points and vertical weights for T and U points
! need to be modified due to the staggered grid.
!-----------------------------------------------------------------------

         il = int(float_local_ijk(m,iblock,1) + p5)
         jl = int(float_local_ijk(m,iblock,2) + p5)

         di = float_local_ijk(m,iblock,1) - il + p5
         dj = float_local_ijk(m,iblock,2) - jl + p5

         float_weights_Tpoint(m,iblock,1) = di
         float_weights_Tpoint(m,iblock,2) = dj

         kl = int(float_local_ijk(m,iblock,3) + p5)

         kp = kl + 1
         kl = max(kl,1)
         kp = min(kp,km)
         if (dk < p5) then
            float_weights_Tpoint(m,iblock,3) = &
               (p5*dz(kl) + dk*dz(kp))*dzwr(kl)
         else
            float_weights_Tpoint(m,iblock,3) = &
               (dk - p5)*dz(kl)*dzwr(kl)
         endif

      endif ! float on this block

   enddo ! float loop

!-----------------------------------------------------------------------
!EOC

   end subroutine calc_float_weights

!***********************************************************************
!BOP
! !IROUTINE: float_move_driver
! !INTERFACE:

   subroutine float_move_driver(left_block)

! !DESCRIPTION:
! Driver for advecting float particles
!-----------------------------------------------------------------------
! Move each float using a predictor-corrector scheme. Trilinear
! interpolation is used to find velocity at float locations.
! Drifters are moved in logical space=> array float_ijk holds
! the logical positions of the floats. Only when float properties
! are dumped (in float_field_value) is the conversion from logical space
! to physical (lat-long-depth) performed. Algorithm provided
! by Doug Kothe.
!-----------------------------------------------------------------------
!
! !REVISION HISTORY:
! same as module

! !USES

   use operators, only : wcalc

! !OUTPUT PARAMETERS:

!EOP
!BOC

!-----------------------------------------------------------------------
!
! arguments
!
!-----------------------------------------------------------------------

   integer(int_kind) :: left_block

!-----------------------------------------------------------------------
!
! local variables
!
!-----------------------------------------------------------------------

   integer(int_kind) :: iblock, float, k

   real (r4) :: dt_move

   real (r4), allocatable, dimension(:,:,:) :: &
      float_predictor_ijk

   real (r8), allocatable, dimension(:,:,:) :: &
      WVELT ! vertical velocity at T-points

   type (block) :: &
     this_block ! block info for current block

   character (char_len) :: &
      predictor_or_corrector

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

   allocate (float_predictor_ijk(num_floats,nblocks_clinic,3))
   allocate (WVELT(nx_block,ny_block,km))
!maltrud debug
! float_predictor_ijk = float_special_value
   float_predictor_ijk = float_ijk
   WVELT = c0

   left_block = 0 ! add 1 for each float that leaves block

   do iblock = 1, nblocks_clinic

      this_block = get_block(blocks_clinic(iblock),iblock)

      call timer_start(float_move_timer, block_id=iblock)

!-----------------------------------------------------------------------
! Calculate weights based on current ijk
!-----------------------------------------------------------------------

      call calc_float_weights(float_ijk, iblock, this_block)

      call wcalc(WVELT, UVEL(:,:,:,curtime,iblock), &
                        VVEL(:,:,:,curtime,iblock), this_block)

      do k = 1, km
         call tgrid_to_ugrid(WVEL(:,:,k,iblock),WVELT(:,:,k), iblock)
      enddo

!-----------------------------------------------------------------------
! Predictor step
!-----------------------------------------------------------------------

      dt_move = p5*dtu
      if (avg_ts) then ! correct for averaging timesteps
         dt_move = p5*dt_move
      endif

      predictor_or_corrector = 'predictor'
      call float_move(float_predictor_ijk, float_ijk, iblock, &
                        this_block, dt_move, predictor_or_corrector)

!-----------------------------------------------------------------------
! Corrector step
! use updated positions to calculate weights, but corrected positions
! are with respect to the original positions.
!-----------------------------------------------------------------------

      dt_move = c2*dt_move

      call calc_float_weights(float_predictor_ijk, iblock, this_block)

      predictor_or_corrector = 'corrector'
      call float_move(float_predictor_ijk, float_ijk, iblock, &
                        this_block, dt_move, predictor_or_corrector)

!-----------------------------------------------------------------------
! Move fully updated positions into the float_ijk array
!-----------------------------------------------------------------------

      float_ijk = float_predictor_ijk

!maltrud need to take care of longitude wrap
!also, ny_global no longer need be land
!probably do this when converting from local to global indices
!or maybe here--i think if blocks dont fit exactly, then maybe
!a particle can have xp>nx_global without leaving the block

      do float = 1, num_floats
         if (float_on_this_block(float,iblock) &
             .and. tday >= float_deploy_time(float)) then
            float_ijk(float,iblock,3) = &
               min( float_ijk(float,iblock,3), (real(km) - 1.e-5_r4) )
            float_ijk(float,iblock,3) = &
               max( float_ijk(float,iblock,3), shallowest_depth_fraction )

            if ( float_ijk(float,iblock,1) < this_block%ib - 1 &
               .or. float_ijk(float,iblock,1) > this_block%ie &
               .or. float_ijk(float,iblock,2) < this_block%jb - 1 &
               .or. float_ijk(float,iblock,2) > this_block%je ) then
                   write(stdout,*)' Drifter #',float,' has left processor #' &
                 , my_task, 'block #', iblock, ' at cycle ',nsteps_total
                  left_block = left_block + 1
            endif
         endif
      enddo

      call timer_stop(float_move_timer, block_id=iblock)

   enddo ! loop over blocks

   deallocate (float_predictor_ijk)
   deallocate (WVELT)

!-----------------------------------------------------------------------

   end subroutine float_move_driver

!***********************************************************************
!BOP
! !IROUTINE: float_move
! !INTERFACE:

   subroutine float_move(float_ijk_new, float_ijk_old, iblock, &
                           this_block, dt_move, predictor_or_corrector)

! !DESCRIPTION:
!
! !REVISION HISTORY:
! same as module

! !OUTPUT PARAMETERS:

!EOP
!BOC

!-----------------------------------------------------------------------
!
! arguments
!
!-----------------------------------------------------------------------

   real(r4), dimension(num_floats,nblocks_clinic,3) :: &
      float_ijk_new, float_ijk_old

   real (r4) :: dt_move

   integer(int_kind) :: iblock

   type (block) :: &
     this_block ! block info for current block

   character (char_len) :: &
      predictor_or_corrector

!-----------------------------------------------------------------------
!
! local variables
!
!-----------------------------------------------------------------------

   logical (log_kind) :: use_lower_level_only

   real(r4) :: xp, yp, zp, di, dj, dk, wijk, wijpk, wijkp, wijpkp, &
               wipjk, wipjpk, wipjkp, wipjpkp, scale, up, vp, wp, &
               di_dx, dj_dy, dk_dz

   integer(int_kind) :: float, i, j, k, ip, jp, kp, ip2, jp2, kk, &
               kmt_c, kmt_n, kmt_s, kmt_e, kmt_w, kmt_ne, kmt_sw, &
               kmt_se, kmt_nw

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

   do float = 1, num_floats

      if (float_on_this_block(float,iblock) &
          .and. tday >= float_deploy_time(float)) then

         if (trim(predictor_or_corrector) == 'predictor') then
            xp = float_ijk_old(float, iblock, 1)
            yp = float_ijk_old(float, iblock, 2)
            zp = float_ijk_old(float, iblock, 3)
         else ! corrector
            xp = float_ijk_new(float, iblock, 1)
            yp = float_ijk_new(float, iblock, 2)
            zp = float_ijk_new(float, iblock, 3)
         endif

!-----------------------------------------------------------------------
! Get the cell indices containing the particle and surrounding cells
! Remember, these are LOCAL INDICES
!-----------------------------------------------------------------------

         i = int(xp)
         j = int(yp)
         k = int(zp)
         k = min(max(k,0),km) ! k can be zero

         ip = i + 1
         jp = j + 1
         kp = k + 1

         kp = min(max(kp,1),km)

         di = float_weights_Upoint(float,iblock,1)
         dj = float_weights_Upoint(float,iblock,2)
         dk = float_weights_Tpoint(float,iblock,3)

         wijk = (1. - di)*(1. - dj)*(1. - dk)
         wijpk = (1. - di)*dj*(1. - dk)
         wijkp = (1. - di)*(1. - dj)*dk
         wijpkp = (1. - di)*dj*dk
         wipjk = di*(1. - dj)*(1. - dk)
         wipjpk = di*dj*(1. - dk)
         wipjkp = di*(1. - dj)*dk
         wipjpkp = di*dj*dk

!-----------------------------------------------------------------------
! check for land in all driections in lower level.
! if so, use only lower level quantities in the interpolation.
!-----------------------------------------------------------------------

         ip2 = ip + 1
         ip2 = min(ip2,nx_block) ! do not reach off block
         jp2 = jp + 1
         jp2 = min(jp2,ny_block) ! do not reach off block

         kmt_c = KMT(ip ,jp , iblock)
         kmt_n = KMT(ip ,jp2, iblock)
         kmt_e = KMT(ip2,jp , iblock)
         kmt_w = KMT(i ,jp , iblock)
         kmt_s = KMT(ip ,j , iblock)
         kmt_ne = KMT(ip2,jp2, iblock)
         kmt_nw = KMT(i ,jp2, iblock)
         kmt_se = KMT(ip2,j , iblock)
         kmt_sw = KMT(i ,j , iblock)

         kk = int(zp + p5)
         kk = max(kk,1)

         use_lower_level_only = .false.

         if (dk > p5) then
            if (kmt_n == kk .and. kmt_c > kk .and. dj > p5) &
              use_lower_level_only = .true.
            if (kmt_e == kk .and. kmt_c > kk .and. di > p5) &
              use_lower_level_only = .true.
            if (kmt_s == kk .and. kmt_c > kk .and. dj < p5) &
              use_lower_level_only = .true.
            if (kmt_w == kk .and. kmt_c > kk .and. di < p5) &
              use_lower_level_only = .true.

            if (kmt_ne == kk .and. kmt_c > kk .and. di > p5 .and. dj > p5) &
              use_lower_level_only = .true.
            if (kmt_nw == kk .and. kmt_c > kk .and. di < p5 .and. dj > p5) &
              use_lower_level_only = .true.
            if (kmt_se == kk .and. kmt_c > kk .and. di > p5 .and. dj < p5) &
              use_lower_level_only = .true.
            if (kmt_sw == kk .and. kmt_c > kk .and. di < p5 .and. dj < p5) &
              use_lower_level_only = .true.
         endif

         if( use_lower_level_only) then
            scale = wijkp + wijpkp + wipjkp + wipjpkp
            if(scale /= c0) then
               scale = c1/scale
            else
               scale = c0
            endif

            up = scale * ( wijkp*UVEL(i ,j ,kp,curtime,iblock) + &
                           wijpkp*UVEL(i ,jp,kp,curtime,iblock) + &
                           wipjkp*UVEL(ip,j ,kp,curtime,iblock) + &
                          wipjpkp*UVEL(ip,jp,kp,curtime,iblock) )
            vp = scale * ( wijkp*VVEL(i ,j ,kp,curtime,iblock) + &
                           wijpkp*VVEL(i ,jp,kp,curtime,iblock) + &
                           wipjkp*VVEL(ip,j ,kp,curtime,iblock) + &
                          wipjpkp*VVEL(ip,jp,kp,curtime,iblock) )
         else ! use both levels
            up = wijk*UVEL(i ,j ,k ,curtime,iblock) + &
                   wijpk*UVEL(i ,jp,k ,curtime,iblock) + &
                   wijkp*UVEL(i ,j ,kp,curtime,iblock) + &
                  wijpkp*UVEL(i ,jp,kp,curtime,iblock) + &
                   wipjk*UVEL(ip,j ,k ,curtime,iblock) + &
                  wipjpk*UVEL(ip,jp,k ,curtime,iblock) + &
                  wipjkp*UVEL(ip,j ,kp,curtime,iblock) + &
                 wipjpkp*UVEL(ip,jp,kp,curtime,iblock)

            vp = wijk*VVEL(i ,j ,k ,curtime,iblock) + &
                   wijpk*VVEL(i ,jp,k ,curtime,iblock) + &
                   wijkp*VVEL(i ,j ,kp,curtime,iblock) + &
                  wijpkp*VVEL(i ,jp,kp,curtime,iblock) + &
                   wipjk*VVEL(ip,j ,k ,curtime,iblock) + &
                  wipjpk*VVEL(ip,jp,k ,curtime,iblock) + &
                  wipjkp*VVEL(ip,j ,kp,curtime,iblock) + &
                 wipjpkp*VVEL(ip,jp,kp,curtime,iblock)
         endif ! lower level only

         if (float_itype(float) == float_advect_3d) then

!maltrud use u-point W so only vertical weights are different

            dk = float_weights_Tpoint(float,iblock,3)

            wijk = (1. - di)*(1. - dj)*(1. - dk)
            wijpk = (1. - di)*dj*(1. - dk)
            wijkp = (1. - di)*(1. - dj)*dk
            wijpkp = (1. - di)*dj*dk
            wipjk = di*(1. - dj)*(1. - dk)
            wipjpk = di*dj*(1. - dk)
            wipjkp = di*(1. - dj)*dk
            wipjpkp = di*dj*dk

            if (k == 0) then
               wp = wijk*WVEL(i ,j ,k ,iblock) + &
                      wijpk*WVEL(i ,jp,k ,iblock) + &
                      wijkp*WVEL(i ,j ,kp,iblock) + &
                     wijpkp*WVEL(i ,jp,kp,iblock) + &
                      wipjk*WVEL(ip,j ,k ,iblock) + &
                     wipjpk*WVEL(ip,jp,k ,iblock) + &
                     wipjkp*WVEL(ip,j ,kp,iblock) + &
                    wipjpkp*WVEL(ip,jp,kp,iblock)

            elseif(k == km) then ! should not happen
               wp = c0
            else
               wp = wijk*WVEL(i ,j ,k ,iblock) + &
                      wijpk*WVEL(i ,jp,k ,iblock) + &
                      wijkp*WVEL(i ,j ,kp,iblock) + &
                     wijpkp*WVEL(i ,jp,kp,iblock) + &
                      wipjk*WVEL(ip,j ,k ,iblock) + &
                     wipjpk*WVEL(ip,jp,k ,iblock) + &
                     wipjkp*WVEL(ip,j ,kp,iblock) + &
                    wipjpkp*WVEL(ip,jp,kp,iblock)
            endif
         else ! constant depth
            wp = c0
         endif

!-----------------------------------------------------------------------
! Compute the derivative of the metric
! Remember that z is positive up and k increases down, so
! dk_dz is actually negative.
!-----------------------------------------------------------------------

         di_dx = (1. - di)/HUS(i,jp,iblock) + di/HUS(ip,jp,iblock)
         dj_dy = (1. - dj)/HUW(ip,j,iblock) + dj/HUW(ip,jp,iblock)
         dk_dz = (1. - dk)*dzwr(k) + dk*dzwr(kp)

         dk_dz = - dk_dz
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

         float_ijk_new(float, iblock, 1) = &
            float_ijk_old(float, iblock, 1) + dt_move*up*di_dx
         float_ijk_new(float, iblock, 2) = &
            float_ijk_old(float, iblock, 2) + dt_move*vp*dj_dy
         float_ijk_new(float, iblock, 3) = &
            float_ijk_old(float, iblock, 3) + dt_move*wp*dk_dz

      endif ! float on this block

   enddo ! loop over floats

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------

   end subroutine float_move

!***********************************************************************

!***********************************************************************
!BOP
! !IROUTINE: float_field_value_T
! !INTERFACE:

 subroutine float_field_value_T(interp_value, iblock, m, bufloc)

! !DESCRIPTION:
! Calculates interpolated value of T-point variables at float location.
!
! !REVISION HISTORY:
! same as module

! !OUTPUT PARAMETERS:

!EOP
!BOC

!-----------------------------------------------------------------------
!
! local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: i,j,ip,jp,il,ipl,jl,jpl,m,k,iblock, bufloc, kl, kpl

   real (r4) :: & ! interpolation weights
      wijk, wipjk, wijpk, wipjpk, wijkp, wipjkp, wijpkp, wipjpkp, weight_sum

   real (r4) :: interp_value

   type (block) :: &
     this_block ! block info for current block

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

   il = int(float_ijk(m,iblock,1) + p5)
   jl = int(float_ijk(m,iblock,2) + p5)
   kl = int(float_ijk(m,iblock,3) + p5)

   ipl = il + 1
   jpl = jl + 1
   kpl = kl + 1

   kl = max(kl,1)
   kpl = min(kpl,km)

   this_block = get_block(iblock,iblock)

!-----------------------------------------------------------------------
! create the trilinear weights from the float_weights_Tpoint array
!-----------------------------------------------------------------------

   wijk = (1.0_r4 - float_weights_Tpoint(m,iblock,1)) &
           *(1.0_r4 - float_weights_Tpoint(m,iblock,2)) &
           *(1.0_r4 - float_weights_Tpoint(m,iblock,3))

   wipjk = float_weights_Tpoint(m,iblock,1) &
           *(1.0_r4 - float_weights_Tpoint(m,iblock,2)) &
           *(1.0_r4 - float_weights_Tpoint(m,iblock,3))

   wijpk = (1.0_r4 - float_weights_Tpoint(m,iblock,1)) &
           *float_weights_Tpoint(m,iblock,2) &
           *(1.0_r4 - float_weights_Tpoint(m,iblock,3))

   wipjpk = float_weights_Tpoint(m,iblock,1) &
           *float_weights_Tpoint(m,iblock,2) &
           *(1.0_r4 - float_weights_Tpoint(m,iblock,3))

   wijkp = (1.0_r4 - float_weights_Tpoint(m,iblock,1)) &
           *(1.0_r4 - float_weights_Tpoint(m,iblock,2)) &
           *float_weights_Tpoint(m,iblock,3)

   wipjkp = float_weights_Tpoint(m,iblock,1) &
           *(1.0_r4 - float_weights_Tpoint(m,iblock,2)) &
           *float_weights_Tpoint(m,iblock,3)

   wijpkp = (1.0_r4 - float_weights_Tpoint(m,iblock,1)) &
           *float_weights_Tpoint(m,iblock,2) &
           *float_weights_Tpoint(m,iblock,3)

   wipjpkp= float_weights_Tpoint(m,iblock,1) &
           *float_weights_Tpoint(m,iblock,2) &
           *float_weights_Tpoint(m,iblock,3)

!-----------------------------------------------------------------------
! For tracers, we do not want to include values on land, so we need
! to renormalize the sum of the weights after setting any to zero.
!-----------------------------------------------------------------------

   if(KMT(il,jl,iblock) < kl) wijk = c0
   if(KMT(ipl,jl,iblock) < kl) wipjk = c0
   if(KMT(il,jpl,iblock) < kl) wijpk = c0
   if(KMT(ipl,jpl,iblock)< kl) wipjpk = c0

   if(KMT(il,jl,iblock) < kpl) wijkp = c0
   if(KMT(ipl,jl,iblock) < kpl) wipjkp = c0
   if(KMT(il,jpl,iblock) < kpl) wijpkp = c0
   if(KMT(ipl,jpl,iblock)< kpl) wipjpkp = c0

   weight_sum = wijk + wipjk + wijpk + wipjpk &
              + wijkp + wipjkp + wijpkp + wipjpkp
   if(weight_sum == c0) weight_sum = eps2 ! small value

   if(weight_sum < c1) then
     wijk = wijk /weight_sum
     wipjk = wipjk /weight_sum
     wijpk = wijpk /weight_sum
     wipjpk = wipjpk /weight_sum
     wijkp = wijkp /weight_sum
     wipjkp = wipjkp /weight_sum
     wijpkp = wijpkp /weight_sum
     wipjpkp = wipjpkp/weight_sum
   endif

   interp_value = wijk *FLOAT_BUF(il, jl, kl, iblock, bufloc) &
                + wipjk *FLOAT_BUF(ipl, jl, kl, iblock, bufloc) &
                + wijpk *FLOAT_BUF(il, jpl, kl, iblock, bufloc) &
                + wipjpk *FLOAT_BUF(ipl, jpl, kl, iblock, bufloc) &
                + wijkp *FLOAT_BUF(il, jl, kpl, iblock, bufloc) &
                + wipjkp *FLOAT_BUF(ipl, jl, kpl, iblock, bufloc) &
                + wijpkp *FLOAT_BUF(il, jpl, kpl, iblock, bufloc) &
                + wipjpkp*FLOAT_BUF(ipl, jpl, kpl, iblock, bufloc)

!-----------------------------------------------------------------------

  end subroutine float_field_value_T

!***********************************************************************
!BOP
! !IROUTINE: float_field_value_U
! !INTERFACE:

 subroutine float_field_value_U(interp_value, iblock, m, bufloc)

! !DESCRIPTION:
! Calculates interpolated value of U-point variables at float location.
!
! !REVISION HISTORY:
! same as module

! !OUTPUT PARAMETERS:

!EOP
!BOC

!-----------------------------------------------------------------------
!
! local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: i,j,ip,jp,il,ipl,jl,jpl,m,k,iblock,bufloc, kl, kpl

   real (r4) :: & ! interpolation weights
      wijk, wipjk, wijpk, wipjpk, wijkp, wipjkp, wijpkp, wipjpkp, weight_sum

   real (r4) :: interp_value

   type (block) :: &
     this_block ! block info for current block

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

   il = int(float_ijk(m,iblock,1))
   jl = int(float_ijk(m,iblock,2))

   ipl = il + 1
   jpl = jl + 1

   kl = int(float_ijk(m,iblock,3) + p5)
   kpl = kl + 1
   kl = max(kl,1)
   kpl = min(kpl,km)

!-----------------------------------------------------------------------
! Calculate trilinear weights and interpolate.
! Assume that we want to include zero values for land points.
!-----------------------------------------------------------------------

   wijk = (1.0_r4 - float_weights_Upoint(m,iblock,1)) &
           *(1.0_r4 - float_weights_Upoint(m,iblock,2)) &
           *(1.0_r4 - float_weights_Upoint(m,iblock,3))

   wipjk = float_weights_Upoint(m,iblock,1) &
           *(1.0_r4 - float_weights_Upoint(m,iblock,2)) &
           *(1.0_r4 - float_weights_Upoint(m,iblock,3))

   wijpk = (1.0_r4 - float_weights_Upoint(m,iblock,1)) &
           *float_weights_Upoint(m,iblock,2) &
           *(1.0_r4 - float_weights_Upoint(m,iblock,3))

   wipjpk = float_weights_Upoint(m,iblock,1) &
           *float_weights_Upoint(m,iblock,2) &
           *(1.0_r4 - float_weights_Upoint(m,iblock,3))

   wijkp = (1.0_r4 - float_weights_Upoint(m,iblock,1)) &
           *(1.0_r4 - float_weights_Upoint(m,iblock,2)) &
           *float_weights_Upoint(m,iblock,3)

   wipjkp = float_weights_Upoint(m,iblock,1) &
           *(1.0_r4 - float_weights_Upoint(m,iblock,2)) &
           *float_weights_Upoint(m,iblock,3)

   wijpkp = (1.0_r4 - float_weights_Upoint(m,iblock,1)) &
           *float_weights_Upoint(m,iblock,2) &
           *float_weights_Upoint(m,iblock,3)

   wipjpkp= float_weights_Upoint(m,iblock,1) &
           *float_weights_Upoint(m,iblock,2) &
           *float_weights_Upoint(m,iblock,3)

   interp_value = wijk *FLOAT_BUF(il, jl, kl, iblock, bufloc) &
                + wipjk *FLOAT_BUF(ipl, jl, kl, iblock, bufloc) &
                + wijpk *FLOAT_BUF(il, jpl, kl, iblock, bufloc) &
                + wipjpk *FLOAT_BUF(ipl, jpl, kl, iblock, bufloc) &
                + wijkp *FLOAT_BUF(il, jl, kpl, iblock, bufloc) &
                + wipjkp *FLOAT_BUF(ipl, jl, kpl, iblock, bufloc) &
                + wijpkp *FLOAT_BUF(il, jpl, kpl, iblock, bufloc) &
                + wipjpkp*FLOAT_BUF(ipl, jpl, kpl, iblock, bufloc)

!-----------------------------------------------------------------------

  end subroutine float_field_value_U

!-----------------------------------------------------------------------
!EOC

!***********************************************************************
!BOP
! !IROUTINE: float_field_value_U_vec
! !INTERFACE:

 subroutine float_field_value_U_vec(interp_value, iblock, m, bufloc)

! !DESCRIPTION:
! Calculates interpolated value of U-point variables at float location.
!
! !REVISION HISTORY:
! same as module

! !OUTPUT PARAMETERS:

!EOP
!BOC

!-----------------------------------------------------------------------
!
! local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: i,j,ip,jp,il,ipl,jl,jpl,m,k,iblock,bufloc, kl, kpl

   real (r4) :: & ! interpolation weights
      wijk, wipjk, wijpk, wipjpk, wijkp, wipjkp, wijpkp, wipjpkp, weight_sum

   real (r4) :: & ! vector field
      vijk, vipjk, vijpk, vipjpk, vijkp, vipjkp, vijpkp, vipjpkp

   real (r4), dimension(2) :: interp_value

   type (block) :: &
     this_block ! block info for current block

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

   il = int(float_ijk(m,iblock,1))
   jl = int(float_ijk(m,iblock,2))

   ipl = il + 1
   jpl = jl + 1

   kl = int(float_ijk(m,iblock,3) + p5)
   kpl = kl + 1
   kl = max(kl,1)
   kpl = min(kpl,km)

!-----------------------------------------------------------------------
! Calculate trilinear weights and interpolate.
! Assume that we want to include zero values for land points.
!-----------------------------------------------------------------------

   wijk = (1.0_r4 - float_weights_Upoint(m,iblock,1)) &
           *(1.0_r4 - float_weights_Upoint(m,iblock,2)) &
           *(1.0_r4 - float_weights_Upoint(m,iblock,3))

   wipjk = float_weights_Upoint(m,iblock,1) &
           *(1.0_r4 - float_weights_Upoint(m,iblock,2)) &
           *(1.0_r4 - float_weights_Upoint(m,iblock,3))

   wijpk = (1.0_r4 - float_weights_Upoint(m,iblock,1)) &
           *float_weights_Upoint(m,iblock,2) &
           *(1.0_r4 - float_weights_Upoint(m,iblock,3))

   wipjpk = float_weights_Upoint(m,iblock,1) &
           *float_weights_Upoint(m,iblock,2) &
           *(1.0_r4 - float_weights_Upoint(m,iblock,3))

   wijkp = (1.0_r4 - float_weights_Upoint(m,iblock,1)) &
           *(1.0_r4 - float_weights_Upoint(m,iblock,2)) &
           *float_weights_Upoint(m,iblock,3)

   wipjkp = float_weights_Upoint(m,iblock,1) &
           *(1.0_r4 - float_weights_Upoint(m,iblock,2)) &
           *float_weights_Upoint(m,iblock,3)

   wijpkp = (1.0_r4 - float_weights_Upoint(m,iblock,1)) &
           *float_weights_Upoint(m,iblock,2) &
           *float_weights_Upoint(m,iblock,3)

   wipjpkp= float_weights_Upoint(m,iblock,1) &
           *float_weights_Upoint(m,iblock,2) &
           *float_weights_Upoint(m,iblock,3)

!-----------------------------------------------------------------------
! now calculate the vector field at all nearest neighbors
! bufloc should hold the i-component, bufloc+1 the j-component
! uzonal = ui*cos(angle) + vi*sin(-angle)
!-----------------------------------------------------------------------

   vijk = FLOAT_BUF(il, jl, kl, iblock, bufloc) * &
             cos( FLOAT_BUF(il, jl, kl, iblock, angle_bufloc)) + &
             FLOAT_BUF(il, jl, kl, iblock, bufloc+1) * &
             sin(-FLOAT_BUF(il, jl, kl, iblock, angle_bufloc))

   vipjk = FLOAT_BUF(ipl, jl, kl, iblock, bufloc) * &
             cos( FLOAT_BUF(ipl, jl, kl, iblock, angle_bufloc)) + &
             FLOAT_BUF(ipl, jl, kl, iblock, bufloc+1) * &
             sin(-FLOAT_BUF(ipl, jl, kl, iblock, angle_bufloc))

   vijpk = FLOAT_BUF(il, jpl, kl, iblock, bufloc) * &
             cos( FLOAT_BUF(il, jpl, kl, iblock, angle_bufloc)) + &
             FLOAT_BUF(il, jpl, kl, iblock, bufloc+1) * &
             sin(-FLOAT_BUF(il, jpl, kl, iblock, angle_bufloc))

   vipjpk = FLOAT_BUF(ipl, jpl, kl, iblock, bufloc) * &
             cos( FLOAT_BUF(ipl, jpl, kl, iblock, angle_bufloc)) + &
             FLOAT_BUF(ipl, jpl, kl, iblock, bufloc+1) * &
             sin(-FLOAT_BUF(ipl, jpl, kl, iblock, angle_bufloc))

   vijkp = FLOAT_BUF(il, jl, kpl, iblock, bufloc) * &
             cos( FLOAT_BUF(il, jl, kpl, iblock, angle_bufloc)) + &
             FLOAT_BUF(il, jl, kpl, iblock, bufloc+1) * &
             sin(-FLOAT_BUF(il, jl, kpl, iblock, angle_bufloc))

   vipjkp = FLOAT_BUF(ipl, jl, kpl, iblock, bufloc) * &
             cos( FLOAT_BUF(ipl, jl, kpl, iblock, angle_bufloc)) + &
             FLOAT_BUF(ipl, jl, kpl, iblock, bufloc+1) * &
             sin(-FLOAT_BUF(ipl, jl, kpl, iblock, angle_bufloc))

   vijpkp = FLOAT_BUF(il, jpl, kpl, iblock, bufloc) * &
             cos( FLOAT_BUF(il, jpl, kpl, iblock, angle_bufloc)) + &
             FLOAT_BUF(il, jpl, kpl, iblock, bufloc+1) * &
             sin(-FLOAT_BUF(il, jpl, kpl, iblock, angle_bufloc))

   vipjpkp = FLOAT_BUF(ipl, jpl, kpl, iblock, bufloc) * &
             cos( FLOAT_BUF(ipl, jpl, kpl, iblock, angle_bufloc)) + &
             FLOAT_BUF(ipl, jpl, kpl, iblock, bufloc+1) * &
             sin(-FLOAT_BUF(ipl, jpl, kpl, iblock, angle_bufloc))

   interp_value(1) = wijk * vijk &
                   + wipjk * vipjk &
                   + wijpk * vijpk &
                   + wipjpk * vipjpk &
                   + wijkp * vijkp &
                   + wipjkp * vipjkp &
                   + wijpkp * vijpkp &
                   + wipjpkp* vipjpkp

!-----------------------------------------------------------------------
! now calculate the meridional component at all nearest neighbors
! vmerid = vi*cos(angle) - ui*sin(-angle)
!-----------------------------------------------------------------------

   vijk = FLOAT_BUF(il, jl, kl, iblock, bufloc+1) * &
             cos( FLOAT_BUF(il, jl, kl, iblock, angle_bufloc)) - &
             FLOAT_BUF(il, jl, kl, iblock, bufloc) * &
             sin(-FLOAT_BUF(il, jl, kl, iblock, angle_bufloc))

   vipjk = FLOAT_BUF(ipl, jl, kl, iblock, bufloc+1) * &
             cos( FLOAT_BUF(ipl, jl, kl, iblock, angle_bufloc)) - &
             FLOAT_BUF(ipl, jl, kl, iblock, bufloc) * &
             sin(-FLOAT_BUF(ipl, jl, kl, iblock, angle_bufloc))

   vijpk = FLOAT_BUF(il, jpl, kl, iblock, bufloc+1) * &
             cos( FLOAT_BUF(il, jpl, kl, iblock, angle_bufloc)) - &
             FLOAT_BUF(il, jpl, kl, iblock, bufloc) * &
             sin(-FLOAT_BUF(il, jpl, kl, iblock, angle_bufloc))

   vipjpk = FLOAT_BUF(ipl, jpl, kl, iblock, bufloc+1) * &
             cos( FLOAT_BUF(ipl, jpl, kl, iblock, angle_bufloc)) - &
             FLOAT_BUF(ipl, jpl, kl, iblock, bufloc) * &
             sin(-FLOAT_BUF(ipl, jpl, kl, iblock, angle_bufloc))

   vijkp = FLOAT_BUF(il, jl, kpl, iblock, bufloc+1) * &
             cos( FLOAT_BUF(il, jl, kpl, iblock, angle_bufloc)) - &
             FLOAT_BUF(il, jl, kpl, iblock, bufloc) * &
             sin(-FLOAT_BUF(il, jl, kpl, iblock, angle_bufloc))

   vipjkp = FLOAT_BUF(ipl, jl, kpl, iblock, bufloc+1) * &
             cos( FLOAT_BUF(ipl, jl, kpl, iblock, angle_bufloc)) - &
             FLOAT_BUF(ipl, jl, kpl, iblock, bufloc) * &
             sin(-FLOAT_BUF(ipl, jl, kpl, iblock, angle_bufloc))

   vijpkp = FLOAT_BUF(il, jpl, kpl, iblock, bufloc+1) * &
             cos( FLOAT_BUF(il, jpl, kpl, iblock, angle_bufloc)) - &
             FLOAT_BUF(il, jpl, kpl, iblock, bufloc) * &
             sin(-FLOAT_BUF(il, jpl, kpl, iblock, angle_bufloc))

   vipjpkp = FLOAT_BUF(ipl, jpl, kpl, iblock, bufloc+1) * &
             cos( FLOAT_BUF(ipl, jpl, kpl, iblock, angle_bufloc)) - &
             FLOAT_BUF(ipl, jpl, kpl, iblock, bufloc) * &
             sin(-FLOAT_BUF(ipl, jpl, kpl, iblock, angle_bufloc))

   interp_value(2) = wijk * vijk &
                   + wipjk * vipjk &
                   + wijpk * vijpk &
                   + wipjpk * vipjpk &
                   + wijkp * vijkp &
                   + wipjkp * vipjkp &
                   + wijpkp * vijpkp &
                   + wipjpkp* vipjpkp

!-----------------------------------------------------------------------

  end subroutine float_field_value_U_vec

!-----------------------------------------------------------------------
!EOC

!***********************************************************************
!BOP
! !IROUTINE: update_float_buffer
! !INTERFACE:

 subroutine update_float_buffer(ARRAY,field_id,block)

! !DESCRIPTION:
! This routine updates a float field to the current value.
!
! !REVISION HISTORY:
! same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
      block, &! local block address (in baroclinic distribution)
      field_id ! index into available fields for float field info

   real (r8), dimension(nx_block,ny_block,km), intent(in) :: &
      ARRAY ! array of data for this block update float buffer

!EOP
!BOC
!-----------------------------------------------------------------------
!
! local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
      bufloc ! location of field in float buffer

!-----------------------------------------------------------------------
!
! get buffer location and field info from avail_float_field array
!
!-----------------------------------------------------------------------

   bufloc = avail_float_fields(field_id)%buf_loc
   if (bufloc <= 0) &
     call exit_POP(sigAbort, &
                    'float: attempt to update bad float field')

!-----------------------------------------------------------------------
!
! update the field into the float buffer
!
!-----------------------------------------------------------------------

   FLOAT_BUF(:,:,:,block,bufloc) = ARRAY

!-----------------------------------------------------------------------
!EOC

 end subroutine update_float_buffer

!***********************************************************************
!BOP
! !IROUTINE: write_float_netcdf
! !INTERFACE:

 subroutine write_float_netcdf(restart_type)

! !DESCRIPTION:
! This routine writes requested float fields to a file. The fields are
! normalized by the time interval before writing.
!
! !REVISION HISTORY:
! same as module

   use netcdf

! !INPUT PARAMETERS:

   character (POP_charLength), intent(in) :: &
      restart_type ! tells float whether to write restart

!EOP
!BOC
!-----------------------------------------------------------------------
!
! local variables
!
!-----------------------------------------------------------------------

   integer (i4) :: &
      nu, &! i/o unit for output file
      iblock, &! dummy block index
      nfield, &! dummy field index
      m, &! dummy float index
      k, &! dummy depth index
      loc ! buffer location for field

   character (char_len), dimension(2) :: &
      file_suffix ! suffix to append to float file name

   character (char_len) :: &
      hist_string, &! string containing file history
      float_filename, &! filename for float data
      float_pointer_file ! filename for pointer file containing
                             ! location/name of last restart file

   character (8) :: &
      date_created ! string with (real) date this file created

   character (10) :: &
      time_created ! string with (real) date this file created

   type (io_field_desc), dimension(:), allocatable :: &
      float_fields

   logical (log_kind), dimension(2) :: &
      lfloat_write ! time to write a file

   logical (log_kind) :: &
      lpositions, &! true if gathering up float positions
      already_gathered !

   real (r4), allocatable, dimension(:) :: &
      temp_1d

   character (char_len), allocatable, dimension(:) :: &
      tempc_1d

   character (255) :: &
      path, &! temp to use for filename
      work_line ! temp to use for character manipulation

   integer (i4) :: &
      ncid, &! netCDF id for file
      iostat, &! status flag for netCDF function calls
      itmp ! integer temp for equivalent logical attribute

   integer (i4) :: & ! dimension and variable ids
      num_id, numv_id, time_id, timev_id, xloc_id, yloc_id, angle_id, &
      iloc_id, jloc_id, zloc_id, kloc_id, deploy_id, type_id, dump_type

   integer (i4), allocatable, dimension(:) :: &
      field_id

   integer (i4), dimension(2) :: &
      start, count, dimids

   logical (log_kind) :: &
      netcdf_error ! error flag for reading attributes

   real (r4), parameter :: float_missing_value = -1.e34

   real (r4), allocatable, dimension(:,:,:) :: &
      float_ijk_temp

!-----------------------------------------------------------------------
!
! is it time to write a file - if yes, create a file suffix
!
!-----------------------------------------------------------------------

   lfloat_write(:) = .false.
   file_suffix(:) = 'none'

   if (lfloat_on) then
      lfloat_write(1) = check_time_flag(float_flag)

      !*** regular float dump
      if (lfloat_write(1)) then
         call create_suffix_float(file_suffix(1))
      endif

      !*** float restart
      if (trim(restart_type) /= 'none') then
         if (.not. lfloat_write(2)) then
            lfloat_write(2) = .true.

            select case (trim(restart_type))
            case('even')
               file_suffix(2) = trim(runid)/&
                                         &/'.even'
            case('odd')
               file_suffix(2) = trim(runid)/&
                                         &/'.odd'
            case('end')
               file_suffix(2) = trim(runid)/&
                                         &/'.end'
            case default
               call create_suffix_float(file_suffix(2))
               file_suffix(2) = trim(file_suffix(2))/&
                                               &/'.restart'
            end select
         endif
      endif
   endif

!-----------------------------------------------------------------------
!
! do the rest only if it is time to do a float dump
!
!-----------------------------------------------------------------------

   already_gathered = .false.

   do dump_type = 1,2

!-----------------------------------------------------------------------
!
! create data file descriptor
!
!-----------------------------------------------------------------------

      float_filename = trim(float_outfile)/&
                                               &/trim(file_suffix(dump_type))/&
                                                                              &/'.nc'

!-----------------------------------------------------------------------
!
! Gather up output data from all processors.
!
!-----------------------------------------------------------------------

     if (lfloat_write(dump_type) .and. .not. already_gathered) then

!-----------------------------------------------------------------------
! float_ijk gets overwritten in gather_float_fields, so copy into a temp
!-----------------------------------------------------------------------

         lpositions = .true.
         allocate( float_ijk_temp(num_floats,nblocks_clinic,3) )
         float_ijk_temp = float_ijk
         call gather_float_fields(float_ijk_global, float_ijk_temp, &
                                    float_on_this_block, master_task, &
                                    distrb_clinic, float_special_value, &
                                    lpositions)
         deallocate(float_ijk_temp)

         lpositions = .false.
         call gather_float_fields(float_xyza_gathered, float_xyza, &
                                    float_on_this_block, master_task, &
                                    distrb_clinic, float_special_value, &
                                    lpositions)

         call gather_float_fields(float_field_values_gathered, &
                                    float_field_values, &
                                    float_on_this_block, master_task, &
                                    distrb_clinic, float_special_value, &
                                    lpositions)

         already_gathered = .true.

      endif ! if time to write and not already_gathered

!-----------------------------------------------------------------------
!
! open output netcdf file
!
!-----------------------------------------------------------------------

     if (lfloat_write(dump_type)) then

      iostat = nf90_noerr

      if (my_task == master_task) then
         write(stdout,*) 'Writing float file: ', trim(float_filename)
         path = trim(float_filename)
         iostat = nf90_create(path=trim(path), cmode=nf90_clobber, ncid=ncid)
         call check(iostat)
      endif

      call broadcast_scalar(iostat, master_task)
      if (iostat /= nf90_noerr) call exit_POP(sigAbort, &
                                              'Error opening file')

      netcdf_error = .false.

      if (my_task == master_task) then

         !*** standard attributes

         path = 'title'
         iostat = nf90_put_att(ncid, NF90_GLOBAL, 'title', &
                               trim(path))
         call check(iostat)
         if (iostat /= nf90_noerr) then
            write(stdout,*) 'Error writing TITLE to float netCDF file'
            netcdf_error = .true.
         endif

         call date_and_time(date=date_created, time=time_created)
         hist_string = char_blank
         write(hist_string,'(a26,a8,1x,a10)') &
            'POP FLOAT file created: ',date_created,time_created

         iostat = nf90_put_att(ncid, NF90_GLOBAL, 'history', &
                               trim(hist_string))
         call check(iostat)
         if (iostat /= nf90_noerr) then
            write(stdout,*) 'Error writing HISTORY to float netCDF file'
            netcdf_error = .true.
         endif

         iostat = nf90_put_att(ncid, NF90_GLOBAL, 'num_float_fields', &
                               num_requested_float_fields)
         call check(iostat)
         if (iostat /= nf90_noerr) then
            write(stdout,*) &
               'Error writing num_requested_float_fields to float netCDF file'
            netcdf_error = .true.
         endif

         iostat = nf90_put_att(ncid, NF90_GLOBAL, 'num_floats', &
                               num_floats)
         call check(iostat)
         if (iostat /= nf90_noerr) then
            write(stdout,*) 'Error writing num_floats to float netCDF file'
            netcdf_error = .true.
         endif

         iostat = nf90_put_att(ncid, NF90_GLOBAL, 'nsteps_total', &
                               nsteps_total)
         call check(iostat)
         if (iostat /= nf90_noerr) then
            write(stdout,*) 'Error writing nsteps_total to float netCDF file'
            netcdf_error = .true.
         endif

         iostat = nf90_put_att(ncid, NF90_GLOBAL, 'time', &
                               tday)
         call check(iostat)
         if (iostat /= nf90_noerr) then
            write(stdout,*) 'Error writing tday to float netCDF file'
            netcdf_error = .true.
         endif

! write(nu,'(i5,i8,f20.8,i10,f15.5,i8,i4,i5)') &
! num_requested_float_fields, num_floats, &
! nsteps_total,tday,iyear,imonth,iday
      endif ! master task

      call broadcast_scalar(netcdf_error, master_task)
      if (netcdf_error) call exit_POP(sigAbort, &
                                      'Error writing file attributes')

      allocate (field_id(num_avail_float_fields))

      if (my_task == master_task) then

! define dimensions
      iostat = NF90_DEF_DIM (ncid=ncid, name='num_floats', &
                             len=num_floats, dimid=num_id)

! define dimension variables
!-----------------------------------------------------------------------
! num_floats
!-----------------------------------------------------------------------
      dimids(1) = num_id
      iostat = NF90_DEF_VAR (ncid=ncid, name='num_floats', &
                             xtype=NF90_INT, dimids=dimids(1), &
                             varid=numv_id)

      !*** long_name
      iostat = NF90_PUT_ATT(ncid=ncid, varid=numv_id, &
                            name='long_name', values='number of floats')
      call check(iostat)
      if (iostat /= NF90_NOERR) netcdf_error = .true.

      !*** units
      iostat = NF90_PUT_ATT(ncid=ncid, varid=numv_id, &
                            name='units', values='float_index')
      call check(iostat)
      if (iostat /= NF90_NOERR) netcdf_error = .true.

! define variables
!-----------------------------------------------------------------------
! time
!-----------------------------------------------------------------------
      dimids(1) = num_id
      iostat = NF90_DEF_VAR (ncid=ncid, name='model_time', &
                             xtype=NF90_FLOAT, dimids=dimids(1), &
                             varid=time_id)

      !*** long_name
      iostat = NF90_PUT_ATT(ncid=ncid, varid=time_id, &
                            name='long_name', values='absolute model time')
      call check(iostat)
      if (iostat /= NF90_NOERR) netcdf_error = .true.

      !*** units
      iostat = NF90_PUT_ATT(ncid=ncid, varid=time_id, &
                            name='units', values='days')
      call check(iostat)
      if (iostat /= NF90_NOERR) netcdf_error = .true.

!-----------------------------------------------------------------------
! xloc
!-----------------------------------------------------------------------
      dimids(1) = num_id
      iostat = NF90_DEF_VAR (ncid=ncid, name='xloc', &
                             xtype=NF90_FLOAT, dimids=dimids(1), &
                             varid=xloc_id)

      !*** long_name
      iostat = NF90_PUT_ATT(ncid=ncid, varid=xloc_id, &
                            name='long_name', values='longitude of floats')
      call check(iostat)
      if (iostat /= NF90_NOERR) netcdf_error = .true.

      !*** units
      iostat = NF90_PUT_ATT(ncid=ncid, varid=xloc_id, &
                            name='units', values='degrees')
      call check(iostat)
      if (iostat /= NF90_NOERR) netcdf_error = .true.

!-----------------------------------------------------------------------
! yloc
!-----------------------------------------------------------------------
      iostat = NF90_DEF_VAR (ncid=ncid, name='yloc', &
                             xtype=NF90_FLOAT, dimids=dimids(1), &
                             varid=yloc_id)

      !*** long_name
      iostat = NF90_PUT_ATT(ncid=ncid, varid=yloc_id, &
                            name='long_name', values='latitude of floats')
      call check(iostat)
      if (iostat /= NF90_NOERR) netcdf_error = .true.

      !*** units
      iostat = NF90_PUT_ATT(ncid=ncid, varid=yloc_id, &
                            name='units', values='degrees')
      call check(iostat)
      if (iostat /= NF90_NOERR) netcdf_error = .true.

!-----------------------------------------------------------------------
! zloc
!-----------------------------------------------------------------------
      iostat = NF90_DEF_VAR (ncid=ncid, name='zloc', &
                             xtype=NF90_FLOAT, dimids=dimids(1), &
                             varid=zloc_id)

      !*** long_name
      iostat = NF90_PUT_ATT(ncid=ncid, varid=zloc_id, &
                            name='long_name', values='depth of floats')
      call check(iostat)
      if (iostat /= NF90_NOERR) netcdf_error = .true.

      !*** units
      iostat = NF90_PUT_ATT(ncid=ncid, varid=zloc_id, &
                            name='units', values='m')
      call check(iostat)
      if (iostat /= NF90_NOERR) netcdf_error = .true.

!-----------------------------------------------------------------------
! angle
!-----------------------------------------------------------------------
      iostat = NF90_DEF_VAR (ncid=ncid, name='angle', &
                             xtype=NF90_FLOAT, dimids=dimids(1), &
                             varid=angle_id)

      !*** long_name
      iostat = NF90_PUT_ATT(ncid=ncid, varid=angle_id, &
               name='long_name', values='grid angle at location  of floats')
      call check(iostat)
      if (iostat /= NF90_NOERR) netcdf_error = .true.

      !*** units
      iostat = NF90_PUT_ATT(ncid=ncid, varid=angle_id, &
                            name='units', values='radians')
      call check(iostat)
      if (iostat /= NF90_NOERR) netcdf_error = .true.

!-----------------------------------------------------------------------
! iloc
!-----------------------------------------------------------------------
      iostat = NF90_DEF_VAR (ncid=ncid, name='iloc', &
                             xtype=NF90_FLOAT, dimids=dimids(1), &
                             varid=iloc_id)

      !*** long_name
      iostat = NF90_PUT_ATT(ncid=ncid, varid=iloc_id, &
                            name='long_name', values='i location of floats')
      call check(iostat)
      if (iostat /= NF90_NOERR) netcdf_error = .true.

      !*** units
      iostat = NF90_PUT_ATT(ncid=ncid, varid=iloc_id, &
                            name='units', values='grid index')
      call check(iostat)
      if (iostat /= NF90_NOERR) netcdf_error = .true.

!-----------------------------------------------------------------------
! jloc
!-----------------------------------------------------------------------
      iostat = NF90_DEF_VAR (ncid=ncid, name='jloc', &
                             xtype=NF90_FLOAT, dimids=dimids(1), &
                             varid=jloc_id)

      !*** long_name
      iostat = NF90_PUT_ATT(ncid=ncid, varid=jloc_id, &
                            name='long_name', values='j location of floats')
      call check(iostat)
      if (iostat /= NF90_NOERR) netcdf_error = .true.

      !*** units
      iostat = NF90_PUT_ATT(ncid=ncid, varid=jloc_id, &
                            name='units', values='grid index')
      call check(iostat)
      if (iostat /= NF90_NOERR) netcdf_error = .true.

!-----------------------------------------------------------------------
! kloc
!-----------------------------------------------------------------------
      iostat = NF90_DEF_VAR (ncid=ncid, name='kloc', &
                             xtype=NF90_FLOAT, dimids=dimids(1), &
                             varid=kloc_id)

      !*** long_name
      iostat = NF90_PUT_ATT(ncid=ncid, varid=kloc_id, &
                            name='long_name', values='k location of floats')
      call check(iostat)
      if (iostat /= NF90_NOERR) netcdf_error = .true.

      !*** units
      iostat = NF90_PUT_ATT(ncid=ncid, varid=kloc_id, &
                            name='units', values='grid index')
      call check(iostat)
      if (iostat /= NF90_NOERR) netcdf_error = .true.

!-----------------------------------------------------------------------
! deployment time
!-----------------------------------------------------------------------
      iostat = NF90_DEF_VAR (ncid=ncid, name='deploy_time', &
                             xtype=NF90_FLOAT, dimids=dimids(1), &
                             varid=deploy_id)

      !*** long_name
      iostat = NF90_PUT_ATT(ncid=ncid, varid=deploy_id, &
                            name='long_name', values='deployment time of floats')
      call check(iostat)
      if (iostat /= NF90_NOERR) netcdf_error = .true.

      !*** units
      iostat = NF90_PUT_ATT(ncid=ncid, varid=deploy_id, &
                            name='units', values='model time in days')
      call check(iostat)
      if (iostat /= NF90_NOERR) netcdf_error = .true.

!-----------------------------------------------------------------------
! float type
!-----------------------------------------------------------------------
      iostat = NF90_DEF_VAR (ncid=ncid, name='float_itype', &
                             xtype=NF90_INT, dimids=dimids(1), &
                             varid=type_id)

      !*** long_name
      iostat = NF90_PUT_ATT(ncid=ncid, varid=type_id, &
                            name='long_name', values='type of floats')
      call check(iostat)
      if (iostat /= NF90_NOERR) netcdf_error = .true.

      !*** units
      iostat = NF90_PUT_ATT(ncid=ncid, varid=type_id, &
                            name='units', values='none')
      call check(iostat)
      if (iostat /= NF90_NOERR) netcdf_error = .true.

!-----------------------------------------------------------------------
! all of the field variables
!-----------------------------------------------------------------------

      do nfield = 1,num_avail_float_fields ! check all available fields
         loc = avail_float_fields(nfield)%buf_loc ! locate field in buffer
         if (loc > 0) then ! field is actually requested and in buffer

            iostat = NF90_DEF_VAR (ncid=ncid, &
                                   name=avail_float_fields(nfield)%short_name, &
                                   xtype=NF90_FLOAT, dimids=dimids(1), &
                                   varid=field_id(nfield))

            !*** long_name
            iostat = NF90_PUT_ATT(ncid=ncid, varid=field_id(nfield), &
                                  name='long_name', &
                                  values=avail_float_fields(nfield)%long_name)
            call check(iostat)
            if (iostat /= NF90_NOERR) netcdf_error = .true.

            !*** units
            iostat = NF90_PUT_ATT(ncid=ncid, varid=field_id(nfield), &
                                  name='units', &
                                  values=avail_float_fields(nfield)%units)
            call check(iostat)
            if (iostat /= NF90_NOERR) netcdf_error = .true.

            !*** missing value
            iostat = NF90_PUT_ATT(ncid=ncid, varid=field_id(nfield), &
                                  name='missing_value', &
                                  values=float_missing_value)
            call check(iostat)
            if (iostat /= NF90_NOERR) netcdf_error = .true.

         endif ! field is requested
      enddo ! loop over available fields

      iostat = nf90_enddef(ncid)
      call check(iostat)

      endif ! master task

      call broadcast_scalar(netcdf_error, master_task)
      if (netcdf_error) call exit_POP(sigAbort, &
                                      'Error writing file attributes')

!-----------------------------------------------------------------------
! write variables
!-----------------------------------------------------------------------

      allocate (temp_1d(num_floats))
      allocate (tempc_1d(num_floats))

      if (my_task == master_task) then

         start(1) = 1
         count(1) = num_floats
         iostat = NF90_PUT_VAR (ncid=ncid, &
                                varid=numv_id, &
                                values=(/ (m, m=1,num_floats) /), &
                                start=start(:), count=count(:))

         temp_1d = float_xyza_gathered(:,1)
         iostat = NF90_PUT_VAR (ncid=ncid, &
                                varid=xloc_id, &
                                values=temp_1d, &
                                start=start(:), count=count(:))

         temp_1d = float_xyza_gathered(:,2)
         iostat = NF90_PUT_VAR (ncid=ncid, &
                                varid=yloc_id, &
                                values=temp_1d, &
                                start=start(:), count=count(:))

         temp_1d = float_xyza_gathered(:,3)
         iostat = NF90_PUT_VAR (ncid=ncid, &
                                varid=zloc_id, &
                                values=temp_1d, &
                                start=start(:), count=count(:))

         temp_1d = float_xyza_gathered(:,4)
         iostat = NF90_PUT_VAR (ncid=ncid, &
                                varid=angle_id, &
                                values=temp_1d, &
                                start=start(:), count=count(:))

         temp_1d = float_ijk_global(:,1)
         iostat = NF90_PUT_VAR (ncid=ncid, &
                                varid=iloc_id, &
                                values=temp_1d, &
                                start=start(:), count=count(:))

         temp_1d = float_ijk_global(:,2)
         iostat = NF90_PUT_VAR (ncid=ncid, &
                                varid=jloc_id, &
                                values=temp_1d, &
                                start=start(:), count=count(:))

         temp_1d = float_ijk_global(:,3)
         iostat = NF90_PUT_VAR (ncid=ncid, &
                                varid=kloc_id, &
                                values=temp_1d, &
                                start=start(:), count=count(:))

         temp_1d = float_deploy_time(:)
         iostat = NF90_PUT_VAR (ncid=ncid, &
                                varid=deploy_id, &
                                values=temp_1d, &
                                start=start(:), count=count(:))

         iostat = NF90_PUT_VAR (ncid=ncid, &
                                varid=type_id, &
                                values=float_itype, &
                                start=start(:), count=count(:))

         temp_1d = tday + c0*float_deploy_time(:) ! easy way to generate tday array
         iostat = NF90_PUT_VAR (ncid=ncid, &
                                varid=time_id, &
                                values=temp_1d, &
                                start=start(:), count=count(:))

         do nfield = 1,num_avail_float_fields ! check all available fields

            loc = avail_float_fields(nfield)%buf_loc ! locate field in buffer

            if (loc > 0) then ! field is actually requested and in buffer

                  temp_1d = float_field_values_gathered(:,loc)
                  start(1) = 1
                  count(1) = num_floats
                  iostat = NF90_PUT_VAR (ncid=ncid, &
                                         varid=field_id(nfield), &
                                         values=temp_1d, &
                                         start=start(:), count=count(:))

            endif
         end do

      iostat = nf90_close(ncid)

      endif ! master_task

      deallocate (temp_1d, tempc_1d, field_id)

!-----------------------------------------------------------------------
!
! if pointer files are used, write float filenames to pointer file
! do this only for float restarts - not float dumps
!
!-----------------------------------------------------------------------

      if (luse_pointer_files .and. lfloat_write(2)) then
         call get_unit(nu)
         if (my_task == master_task) then
            float_pointer_file = trim(pointer_filename)/&
                                                       &/'.float'

            open(nu,file=float_pointer_file,form='formatted', &
                    status='unknown')
            write(nu,'(a)') trim(float_filename)
            close(nu)
         endif
         call release_unit(nu)
      endif

   endif ! lfloat_write

   enddo ! dump_type

!-----------------------------------------------------------------------
!EOC

 end subroutine write_float_netcdf

!***********************************************************************
!BOP
! !IROUTINE: read_float_netcdf
! !INTERFACE:

 subroutine read_float_netcdf

! !DESCRIPTION:
! This routine reads a time average restart dump to continue
! running time averages of requested fields.
!
! !REVISION HISTORY:
! same as module

   use netcdf

!EOP
!BOC
!-----------------------------------------------------------------------
!
! local variables
!
!-----------------------------------------------------------------------

   integer (i4) :: &
     nu, & ! i/o unit
     iblock, & ! dummy block index
     n,i, & ! dummy for indexing character string
     in_fields, & ! num of fields in restart file
     nfield, & ! dummy field counter
     hdr_error, & ! error file for reading restart hdr
     in_nsteps_total, & ! nsteps_total according to float file
     in_iyear, & ! iyear according to float file
     in_imonth, & ! imonth according to float file
     in_iday, & ! iday according to float file
     in_num_float_fields, & ! num_requested_float_fields according to float file
     in_num_floats, & ! num_floats according to float file
     loc ! buffer location

   real (r4) :: &
     in_tday, & ! tday according to float file
     in1,in2,in3 ! placeholders for non-saved values read from float file

   character (char_len) :: &
     in_field_name , & ! temporary
     char_temp, & ! for string manipulation
     float_pointer_file ! filename for pointer file containing
                           ! location/name of last restart file

   logical (log_kind) :: restart_error

   integer (i4), dimension(2) :: start, count
   integer (i4), dimension(5) :: field_id

   integer (i4) :: iostat, ncid, m

   real (r4), allocatable, dimension(:) :: &
      temp_1d

!-----------------------------------------------------------------------
!
! if pointer files are used, pointer file and must be read to get
! actual filenames
!
!-----------------------------------------------------------------------

   call get_unit(nu)

   if (luse_pointer_files) then

      if (my_task == master_task) then
         float_pointer_file = char_blank
         float_pointer_file = trim(pointer_filename)/&
                                                   &/'.float'
         write(stdout,*) 'Reading pointer file: ', &
                         trim(float_pointer_file)
         open(nu, file=trim(float_pointer_file), form='formatted', &
                  status='old')
         read(nu,'(a)') float_infile
         close(nu)
      endif
      call broadcast_scalar(float_infile, master_task)

   endif

   call release_unit(nu)

!-----------------------------------------------------------------------
!
! open input file
!
! check for consistency of restart file with float_contents file
!
!-----------------------------------------------------------------------

   iostat = nf90_noerr

   if (my_task == master_task) then
      char_temp = trim(float_infile)
      iostat = nf90_open(path=trim(char_temp), mode=nf90_nowrite, ncid=ncid)
      call check(iostat)
   endif

   call broadcast_scalar(iostat, master_task)
   if (iostat /= nf90_noerr) &
      call exit_POP(sigAbort,'error opening netCDF file for reading')

   call broadcast_scalar(ncid, master_task)

!-----------------------------------------------------------------------
!
! read global file attributes
!
!-----------------------------------------------------------------------

   restart_error = .false.

   if (my_task == master_task) then

      iostat = nf90_get_att(ncid=ncid, varid=NF90_GLOBAL, &
               name='num_float_fields',values=in_num_float_fields)
      call check(iostat)
      if (iostat /= nf90_noerr) then
         restart_error = .true.
         write(stdout,*) 'FLOAT:error reading num_float_fields'
      endif

      iostat = nf90_get_att(ncid=ncid, varid=NF90_GLOBAL, &
               name='num_floats',values=in_num_floats)
      call check(iostat)
      if (iostat /= nf90_noerr) then
         restart_error = .true.
         write(stdout,*) 'FLOAT:error reading num_floats'
      endif

      iostat = nf90_get_att(ncid=ncid, varid=NF90_GLOBAL, &
               name='nsteps_total',values=in_nsteps_total)
      call check(iostat)
      if (iostat /= nf90_noerr) then
         restart_error = .true.
         write(stdout,*) 'FLOAT:error reading nsteps_total'
      endif

   endif ! master_task

   call broadcast_scalar(restart_error, master_task)

   if (restart_error) &
      call exit_POP(sigAbort,'FLOAT: error reading global attributes')

   if (my_task == master_task) then

      !*** check nsteps total for validity
      if (in_nsteps_total /= nsteps_total) then
         write(stdout,'(i6,a32,i6,a35)') &
            in_nsteps_total,' nsteps_total in float restart; ', &
            nsteps_total, ' nsteps_total in current simulation'
         restart_error = .true.
         write(stdout,*)'FLOAT:restart file has wrong time step?'
      endif

      !*** check number of requested float fields for validity
      if (in_num_float_fields /= num_requested_float_fields) then
         write(stdout,'(i6,a44,i6,a47)') &
            in_num_float_fields, &
            'requested float fields in float restart; ', &
            num_requested_float_fields, &
            ' requested float fields in current simulation'
         restart_error = .true.
         write(stdout,*)'FLOAT:restart file has wrong number of fields'
      endif

      !*** check number of floats for validity
      if (in_num_floats /= num_floats) then
         write(stdout,'(i6,a35,i6,a40)') &
            in_num_floats,' float locations in float restart; ', &
            num_floats, ' float locations in current simulation'
         restart_error = .true.
         write(stdout,*)'FLOAT:restart file has wrong number of floats'
      endif

   endif ! master_task

   call broadcast_scalar(restart_error, master_task)

   if (restart_error) &
      call exit_POP(sigAbort,'FLOAT: restart file failed consistency checks')

!-----------------------------------------------------------------------
! Now read in values
!-----------------------------------------------------------------------

   if (my_task == master_task) then

      in_field_name = 'iloc'
      iostat = NF90_INQ_VARID(ncid=ncid, &
                              name=trim(in_field_name), &
                              varid=field_id(1))
      call check(iostat)
      if (iostat /= nf90_noerr) then
         write(stdout,*) 'FLOAT:could not find field: ', trim(in_field_name)
         restart_error = .true.
         char_temp = 'FLOAT:error finding float field'
      endif

      in_field_name = 'jloc'
      iostat = NF90_INQ_VARID(ncid=ncid, &
                              name=trim(in_field_name), &
                              varid=field_id(2))
      call check(iostat)
      if (iostat /= nf90_noerr) then
         write(stdout,*) 'FLOAT:could not find field: ', trim(in_field_name)
         restart_error = .true.
         char_temp = 'FLOAT:error finding float field'
      endif

      in_field_name = 'kloc'
      iostat = NF90_INQ_VARID(ncid=ncid, &
                              name=trim(in_field_name), &
                              varid=field_id(3))
      call check(iostat)
      if (iostat /= nf90_noerr) then
         write(stdout,*) 'FLOAT:could not find field: ', trim(in_field_name)
         restart_error = .true.
         char_temp = 'FLOAT:error finding float field'
      endif

      in_field_name = 'deploy_time'
      iostat = NF90_INQ_VARID(ncid=ncid, &
                              name=trim(in_field_name), &
                              varid=field_id(4))
      call check(iostat)
      if (iostat /= nf90_noerr) then
         write(stdout,*) 'FLOAT:could not find field: ', trim(in_field_name)
         restart_error = .true.
         char_temp = 'FLOAT:error finding float field'
      endif

      in_field_name = 'float_itype'
      iostat = NF90_INQ_VARID(ncid=ncid, &
                              name=trim(in_field_name), &
                              varid=field_id(5))
      call check(iostat)
      if (iostat /= nf90_noerr) then
         write(stdout,*) 'FLOAT:could not find field: ', trim(in_field_name)
         restart_error = .true.
         char_temp = 'FLOAT:error finding float field'
      endif

      allocate (temp_1d(num_floats))

      start(1) = 1
      count(1) = num_floats
      start(2) = 1
      count(2) = 1
      iostat = NF90_GET_VAR(ncid=ncid, &
                            varid=field_id(1), &
                            values=temp_1d)
      call check(iostat)
      if (iostat /= nf90_noerr) then
         write(stdout,*) 'FLOAT:could not read field: iloc'
         restart_error = .true.
         char_temp = 'FLOAT:error reading float field'
      endif
      float_ijk_global(:,1) = temp_1d

      iostat = NF90_GET_VAR(ncid=ncid, &
                            varid=field_id(2), &
                            values=temp_1d)
      call check(iostat)
      if (iostat /= nf90_noerr) then
         write(stdout,*) 'FLOAT:could not read field: jloc'
         restart_error = .true.
         char_temp = 'FLOAT:error reading float field'
      endif
      float_ijk_global(:,2) = temp_1d

      iostat = NF90_GET_VAR(ncid=ncid, &
                            varid=field_id(3), &
                            values=temp_1d)
      call check(iostat)
      if (iostat /= nf90_noerr) then
         write(stdout,*) 'FLOAT:could not read field: kloc'
         restart_error = .true.
         char_temp = 'FLOAT:error reading float field'
      endif
      float_ijk_global(:,3) = temp_1d

      iostat = NF90_GET_VAR(ncid=ncid, &
                            varid=field_id(4), &
                            values=float_deploy_time)
      call check(iostat)
      if (iostat /= nf90_noerr) then
         write(stdout,*) 'FLOAT:could not read field: deploy_time'
         restart_error = .true.
         char_temp = 'FLOAT:error reading float field'
      endif

      start(1) = 1
      count(1) = num_floats
      iostat = NF90_GET_VAR(ncid=ncid, &
                            varid=field_id(5), &
                            values=float_itype)
      call check(iostat)
      if (iostat /= nf90_noerr) then
         write(stdout,*) 'FLOAT:could not read field: float_type'
         restart_error = .true.
         char_temp = 'FLOAT:error reading float field'
      endif

      deallocate (temp_1d)

   endif ! master_task

   call broadcast_scalar(restart_error, master_task)
   call broadcast_scalar(char_temp , master_task)

   if (restart_error) &
      call exit_POP(sigAbort,trim(char_temp))

!-----------------------------------------------------------------------
! convert float_itype to float_type
!-----------------------------------------------------------------------

   call broadcast_array(float_itype, master_task)
   do n = 1, num_floats
      if (float_itype(n) == float_advect_2d) then
         float_type(n) = 'constant-depth'
      elseif (float_itype(n) == float_advect_3d) then
         float_type(n) = '3d'
      else
         call exit_POP(sigAbort,'ERROR: unknown float_itype')
      endif
   enddo

!-----------------------------------------------------------------------
!EOC

 end subroutine read_float_netcdf

!***********************************************************************
!EOP
! !IROUTINE: float_global_to_local
! !INTERFACE:

 subroutine float_global_to_local

! !DESCRIPTION:
!
! !REVISION HISTORY:
! same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!
! local variables
!
!-----------------------------------------------------------------------

   integer (i4) :: &
      iblock, imin, imax, jmin, jmax, float

   type (block) :: &
     this_block ! block info for current block

!-----------------------------------------------------------------------
!
! convert from global i,j,k to local block i,j,k
!
!-----------------------------------------------------------------------

   do iblock = 1, nblocks_clinic

      this_block = get_block(blocks_clinic(iblock),iblock)

      do float = 1, num_floats
         imin = this_block%i_glob(this_block%ib) - 1
         imax = this_block%i_glob(this_block%ie)
         jmin = this_block%j_glob(this_block%jb) - 1
         jmax = this_block%j_glob(this_block%je)
         if ( float_ijk_global(float,1) > imin &
           .and. float_ijk_global(float,1) <= imax &
           .and. float_ijk_global(float,2) > jmin &
           .and. float_ijk_global(float,2) <= jmax ) then

            float_on_this_block(float,iblock) = .true.
            float_ijk(float,iblock,1) = &
               float_ijk_global(float,1) - imin + nghost
            float_ijk(float,iblock,2) = &
               float_ijk_global(float,2) - jmin + nghost
            float_ijk(float,iblock,3) = float_ijk_global(float,3)

!maltrud debug
! write(stdout,*)' found float: #, proc, block, global_ijk, local_ijk = ', &
! float, my_task, iblock, float_ijk_global(float,:), &
! float_ijk(float,iblock,:)
         else
            float_on_this_block(float,iblock) = .false.
            float_ijk(float,iblock,:) = float_special_value
         endif

!maltrud i do not know if these are necessary
!maltrud longitude wrap
! if (float_local_ijk(m,iblock,1) < c0) &
! float_local_ijk(m,iblock,1) = &
! float_local_ijk(m,iblock,1) + nx_global
!maltrud fix bottom row where this_block%j_glob(1) == 0
! if (this_block%j_glob(1) == 0) &
! float_local_ijk(m,iblock,2) = &
! float_local_ijk(m,iblock,2) + nghost


      enddo
   enddo

!-----------------------------------------------------------------------
!EOC

 end subroutine float_global_to_local

!***********************************************************************
!BOP
! !IROUTINE: check
! !INTERFACE:

 subroutine check(status)

! !DESCRIPTION:
! This exception handler subroutine can be used to check error status
! after a netcdf call. It prints out a text message assigned to
! an error code but does not exit because this routine is typically
! only called from a single process.
!
! !REVISION HISTORY:
! same as module

   use netcdf

! !INPUT PARAMETERS:

   integer (i4), intent (in) :: &
      status ! status returned by netCDF call

!EOP
!BOC
!-----------------------------------------------------------------------
!
! call netCDF routine to return error message
!
!-----------------------------------------------------------------------

   if (status /= nf90_noerr) then
      write(stdout,*) trim(nf90_strerror(status))
   end if

!-----------------------------------------------------------------------
!EOC

 end subroutine check

!***********************************************************************

!***********************************************************************

 end module floats

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
