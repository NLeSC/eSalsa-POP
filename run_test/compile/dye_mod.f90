!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

module dye_mod

!BOP
! !MODULE: dye_mod
!
! !DESCRIPTION:
!
! !REVISION HISTORY:

! !USES:

   use POP_KindsMod
   use POP_IOUnitsMod
   use POP_DomainSizeMod
   use POP_ErrorMod

   use blocks, only: nx_block, ny_block
   use domain_size, only: max_blocks_clinic, km
   use domain, only: nblocks_clinic
   use exit_mod, only: sigAbort, exit_POP
   use communicate, only: my_task, master_task
   use prognostic, only: tracer_field
   use kinds_mod
   use constants, only: c0, c1, char_blank, delim_fmt, field_loc_center, &
       field_type_scalar, blank_fmt
   use io, only: data_set
   use io_types, only: stdout, nml_in, nml_filename, io_field_desc, io_dim, &
       construct_file, construct_io_dim, construct_io_field, rec_type_dbl, &
       datafile, destroy_io_field, destroy_file
   use io_tools, only: document
   use passive_tracer_tools, only: ind_name_pair, tracer_read, &
       rest_read_tracer_block, file_read_tracer_block
   use movie, only: define_movie_field, movie_requested, update_movie_field
   implicit none
   private

! !PUBLIC MEMBER FUNCTIONS:

   public :: dye_tracer_cnt, &
             dye_init, &
             dye_set_interior, &
             dye_set_sflux, &
             dye_movie, &
             dye_reset

!EOP
!BOC

!-----------------------------------------------------------------------
! module variables required by passive_tracer
!-----------------------------------------------------------------------

   integer(int_kind), parameter :: &
      dye_tracer_cnt = 1

!-----------------------------------------------------------------------
! relative tracer indices
!-----------------------------------------------------------------------

   integer(int_kind), parameter :: &
      dye_ind = 1 ! dye index

!-----------------------------------------------------------------------
! derived type & parameter for tracer index lookup
!-----------------------------------------------------------------------

   type(ind_name_pair), dimension(dye_tracer_cnt) :: &
      ind_name_table = (/ ind_name_pair(dye_ind, 'DYE') /)

!-----------------------------------------------------------------------
! surface flux masks
!-----------------------------------------------------------------------

   character(char_len) :: &
      dye_region_file, & ! filename for surface region
      dye_region_file_fmt ! file format for surface region

   logical(POP_Logical), dimension(:,:,:), allocatable, save :: &
       LAND_MASK

   real (POP_r8), dimension(:,:,:,:), allocatable, save :: &
       DYE_SURFACE_REGION_MASK

!-----------------------------------------------------------------------------
! define movie id for DYE
!-----------------------------------------------------------------------------

   integer (POP_i4), dimension(0:km,dye_tracer_cnt) :: &
      movie_DYE ! movie id for all DYEs
                   ! k=0 holds depth-integrated field

!EOC
!*****************************************************************************

contains

!*****************************************************************************
!BOP
! !IROUTINE: dye_init
! !INTERFACE:

 subroutine dye_init(init_ts_file_fmt, read_restart_filename, &
                      tracer_d_module, TRACER_MODULE, tadvect_ctype, errorCode)

! !DESCRIPTION:
! Initialize dye tracer module. This involves setting metadata, reading
! the module namelist and setting initial conditions.
!
! !REVISION HISTORY:
! same as module

! !USES:

   use broadcast, only: broadcast_scalar
   use prognostic, only: curtime, oldtime
   use grid, only: KMT

! !INPUT PARAMETERS:

   character (*), intent(in) :: &
      init_ts_file_fmt, & ! format (bin or nc) for input file
      read_restart_filename ! file name for restart file

! !INPUT/OUTPUT PARAMETERS:

   type (tracer_field), dimension(dye_tracer_cnt), intent(inout) :: &
      tracer_d_module ! descriptors for each tracer

   real(r8), dimension(nx_block,ny_block,km,dye_tracer_cnt,3,max_blocks_clinic), &
      intent(inout) :: TRACER_MODULE

! !OUTPUT PARAMETERS:

   character (char_len), dimension(dye_tracer_cnt), intent(out) :: &
      tadvect_ctype ! advection method for dye tracers

   integer (POP_i4), intent(out) :: &
      errorCode ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
! local variables
!-----------------------------------------------------------------------

   character(*), parameter :: subname = 'dye_mod:dye_init'

   character(char_len) :: &
      init_dye_option, & ! option for initialization of dye
      init_dye_init_file, & ! filename for option 'file'
      init_dye_init_file_fmt, & ! file format for option 'file'
      dye_tadvect_ctype ! advection method for ecosys tracers

   logical(log_kind) :: &
      lnml_found ! Was dye_nml found ?

   integer(int_kind) :: &
      n, & ! index for looping over tracers
      k, & ! index for looping over depth levels
      iblock, & ! index for looping over blocks
      nml_error ! namelist i/o error flag

! l, & ! index for looping over time levels

   type(tracer_read), dimension(dye_tracer_cnt) :: &
      tracer_init_ext ! namelist variable for initializing tracers

   namelist /dye_nml/ &
      init_dye_option, init_dye_init_file, tracer_init_ext, &
      init_dye_init_file_fmt, dye_region_file, dye_region_file_fmt, &
      dye_tadvect_ctype

   character (char_len) :: &
      dye_restart_filename ! modified file name for restart file

!-----------------------------------------------------------------------
! initialize tracer_d values
!-----------------------------------------------------------------------

   errorCode = POP_Success

   tracer_d_module(dye_ind)%short_name = 'DYE'
   tracer_d_module(dye_ind)%long_name = 'Dye'
   tracer_d_module(dye_ind)%units = 'conc'
   tracer_d_module(dye_ind)%tend_units = 'conc/s'
   tracer_d_module(dye_ind)%flux_units = 'cm conc/s'

!-----------------------------------------------------------------------
! default namelist settings
!-----------------------------------------------------------------------

   init_dye_option = 'unknown'
   init_dye_init_file = 'unknown'
   init_dye_init_file_fmt = 'bin'
   dye_region_file = 'unknown'
   dye_region_file_fmt = 'bin'
   dye_tadvect_ctype = 'base_model'

   do n = 1,dye_tracer_cnt
      tracer_init_ext(n)%mod_varname = 'unknown'
      tracer_init_ext(n)%filename = 'unknown'
      tracer_init_ext(n)%file_varname = 'unknown'
      tracer_init_ext(n)%scale_factor = c1
      tracer_init_ext(n)%default_val = c0
      tracer_init_ext(n)%file_fmt = 'bin'
   end do

   if (my_task == master_task) then
      open (nml_in, file=nml_filename, status='old',iostat=nml_error)
      if (nml_error /= 0) then
         nml_error = -1
      else
         nml_error = 1
      endif
      do while (nml_error > 0)
         read(nml_in, nml=dye_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
   endif

   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
      call document(subname, 'dye_nml not found')
      call exit_POP(sigAbort, 'stopping in ' /&
                           &/ subname)
   endif

!-----------------------------------------------------------------------
! broadcast all namelist variables
!-----------------------------------------------------------------------

   call broadcast_scalar(init_dye_option , master_task)
   call broadcast_scalar(init_dye_init_file, master_task)
   call broadcast_scalar(init_dye_init_file_fmt, master_task)
   call broadcast_scalar(dye_region_file, master_task)
   call broadcast_scalar(dye_region_file_fmt, master_task)
   call broadcast_scalar(dye_tadvect_ctype, master_task)
   tadvect_ctype(:) = dye_tadvect_ctype

   do n = 1,dye_tracer_cnt
      call broadcast_scalar(tracer_init_ext(n)%mod_varname, master_task)
      call broadcast_scalar(tracer_init_ext(n)%filename, master_task)
      call broadcast_scalar(tracer_init_ext(n)%file_varname, master_task)
      call broadcast_scalar(tracer_init_ext(n)%scale_factor, master_task)
      call broadcast_scalar(tracer_init_ext(n)%default_val, master_task)
      call broadcast_scalar(tracer_init_ext(n)%file_fmt, master_task)
   end do

!-----------------------------------------------------------------------
! initialize tracers
!-----------------------------------------------------------------------

   select case (init_dye_option)

   case ('startup', 'zero', 'startup_spunup')
      TRACER_MODULE = c0
      if (my_task == master_task) then
          write(stdout,delim_fmt)
          write(stdout,*) ' Initial 3-d Dye concentration set to all zeros'
          write(stdout,delim_fmt)
          call POP_IOUnitsFlush(POP_stdout)
      endif

   case ('restart', 'continue', 'branch', 'hybrid' )

      dye_restart_filename = char_blank

      if (init_dye_init_file == 'same_as_TS') then
         if (read_restart_filename == 'undefined') then
            call document(subname, 'no restart file to read dye from')
            call exit_POP(sigAbort, 'stopping in ' /&
                                 &/ subname)
         endif
         dye_restart_filename = read_restart_filename
         init_dye_init_file_fmt = init_ts_file_fmt

      else ! do not read from TS restart file

         dye_restart_filename = trim(init_dye_init_file)

      endif

      call rest_read_tracer_block(init_dye_init_file_fmt, &
                                  dye_restart_filename, &
                                  tracer_d_module, &
                                  TRACER_MODULE)

   case ('file')
      call document(subname, 'dye being read from separate file')

      call file_read_tracer_block(init_dye_init_file_fmt, &
                                  init_dye_init_file, &
                                  tracer_d_module, &
                                  ind_name_table, &
                                  tracer_init_ext, &
                                  TRACER_MODULE)

   case default
      call document(subname, 'unknown init_dye_option = ', init_dye_option)
      call exit_POP(sigAbort, 'stopping in ' /&
                           &/ subname)

   end select

!-----------------------------------------------------------------------
! apply land mask to tracers
!-----------------------------------------------------------------------

   do iblock=1,nblocks_clinic
      do n = 1,dye_tracer_cnt
         do k = 1,km
            where (k > KMT(:,:,iblock))
               TRACER_MODULE(:,:,k,n,curtime,iblock) = c0
               TRACER_MODULE(:,:,k,n,oldtime,iblock) = c0
            end where
         end do
      end do
   enddo

!-----------------------------------------------------------------------
! call other initialization subroutines
!-----------------------------------------------------------------------

   call dye_init_sflux(tracer_d_module)
   call dye_init_movie

    !call init_dye_tavg
    !call init_dye_mooring


!-----------------------------------------------------------------------
!EOC

 end subroutine dye_init

!***********************************************************************
!BOP
! !IROUTINE: dye_set_interior
! !INTERFACE:

 subroutine dye_set_interior(k, DTRACER_MODULE)

! !DESCRIPTION:
! set interior source/sink term for ideal age tracer
!
! !REVISION HISTORY:
! same as module

! !USES:

   use time_management, only: seconds_in_year

! !INPUT PARAMETERS:

   integer(int_kind), intent(in) :: &
      k ! vertical level index

! !OUTPUT PARAMETERS:

   real(r8), dimension(nx_block,ny_block,dye_tracer_cnt), intent(out) :: &
      DTRACER_MODULE ! computed source/sink term

!EOP
!BOC

!-----------------------------------------------------------------------

    DTRACER_MODULE = c0

!-----------------------------------------------------------------------
!EOC

 end subroutine dye_set_interior

!***********************************************************************
!BOP
! !IROUTINE: dye_reset
! !INTERFACE:

 subroutine dye_reset(TRACER_MODULE)

! !DESCRIPTION:
! reset surface value for ideal age tracer
!
! !REVISION HISTORY:
! same as module

! !INPUT/OUTPUT PARAMETERS:

   real(r8), dimension(nx_block,ny_block,km,dye_tracer_cnt), intent(inout) :: &
      TRACER_MODULE ! ideal age tracer

!EOP
!BOC

!-----------------------------------------------------------------------

! no resetting

!-----------------------------------------------------------------------
!EOC

 end subroutine dye_reset

!***********************************************************************
!***********************************************************************
!BOP
! !IROUTINE: dye_init_sflux
! !INTERFACE:

 subroutine dye_init_sflux(tracer_d_module)

! !DESCRIPTION:
! initialize surfcae flux fields for dye tracer
!
! !REVISION HISTORY:
! same as module

   use grid, only: KMT

! !INPUT/OUTPUT PARAMETERS:

   type (tracer_field), dimension(dye_tracer_cnt), intent(inout) :: &
      tracer_d_module ! descriptors for each tracer

!---------------------------------------------------------------------------
! local variables
!---------------------------------------------------------------------------

   type (io_field_desc) :: &
      io_field ! io field descriptors for input Tracer

   type (datafile) :: &
      io_file ! io file descriptor

   type (io_dim) :: &
      i_dim, j_dim, k_dim ! dimension descriptors

   character (char_len) :: &
      region_filename ! modified file name for surface region file

   integer (int_kind) :: &
      n

!EOP
!BOC

!---------------------------------------------------------------------------
! read in surface forcing region masks
!---------------------------------------------------------------------------

  if (my_task == master_task) then
      write(stdout,delim_fmt)
      write(stdout,*) ' Reading DYE surface regions from ' /&
                               &/trim(region_filename)
      write(stdout,delim_fmt)
      call POP_IOUnitsFlush(POP_stdout)
  endif

    io_file = construct_file(dye_region_file_fmt, &
                      full_name=trim(dye_region_file), &
                      record_length = rec_type_dbl, &
                      recl_words=POP_nxGlobal*POP_nyGlobal)
    call data_set(io_file,'open_read')

    i_dim = construct_io_dim('i',POP_nxGlobal)
    j_dim = construct_io_dim('j',POP_nyGlobal)
    k_dim = construct_io_dim('k',dye_tracer_cnt)

    allocate (DYE_SURFACE_REGION_MASK(nx_block,ny_block, &
                                      dye_tracer_cnt,max_blocks_clinic) )
    DYE_SURFACE_REGION_MASK = c0

!maltrud not general for more than 1 dye tracer
    do n = 1, dye_tracer_cnt
    io_field = &
        construct_io_field(trim(tracer_d_module(n)%short_name), &
        dim1=i_dim, dim2=j_dim, dim3=k_dim, &
        field_loc = field_loc_center, &
        field_type = field_type_scalar, &
        d3d_array=DYE_SURFACE_REGION_MASK)

    call data_set(io_file,'define',io_field)

    call data_set(io_file,'read' ,io_field)

    call destroy_io_field(io_field)
    enddo

    call data_set(io_file,'close')
    call destroy_file(io_file)

    if (my_task == master_task) then
       write(stdout,blank_fmt)
       write(stdout,'(a12,a)') ' file read: ', &
          trim(dye_region_file)
    endif

    !---------------------------------------------------------------------------
    ! allocate and initialize LAND_MASK (true for ocean points)
    !---------------------------------------------------------------------------

    allocate( LAND_MASK(nx_block,ny_block,max_blocks_clinic) )
    LAND_MASK = merge(.true., .false., KMT > 0)

!-----------------------------------------------------------------------
!EOC

 end subroutine dye_init_sflux

!***********************************************************************
!***********************************************************************
!BOP
! !IROUTINE: dye_set_sflux
! !INTERFACE:

 subroutine dye_set_sflux(STF_MODULE)

! !DESCRIPTION:
! set surface flux for dye tracer
!
! !REVISION HISTORY:
! same as module

! !INPUT/OUTPUT PARAMETERS:

 real(POP_r8), dimension(nx_block,ny_block,dye_tracer_cnt,max_blocks_clinic), &
    intent(inout) :: STF_MODULE

!---------------------------------------------------------------------------
! local variables
!---------------------------------------------------------------------------

 integer(POP_i4) :: &
    n, iblock ! loop indices

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

 STF_MODULE = c0

 !$OMP PARALLEL DO PRIVATE(iblock)
 do iblock = 1, nblocks_clinic
    do n = 1, dye_tracer_cnt
       where (LAND_MASK(:,:,iblock)) &
          STF_MODULE(:,:,n,iblock) = DYE_SURFACE_REGION_MASK(:,:,n,iblock)
    enddo
 enddo
 !$OMP END PARALLEL DO

!EOP
!BOC

!EOC

 end subroutine dye_set_sflux

!***********************************************************************

!BOP
! !IROUTINE: dye_init_movie
! !INTERFACE:

 subroutine dye_init_movie

! !DESCRIPTION:
!---------------------------------------------------------------------------
! set up variables for movie access
! 1) allocate single precision data buffers
! 2) initialize buffers to zero
! 3) register buffers so that movie can access them
!---------------------------------------------------------------------------
!
! !REVISION HISTORY:
! same as module

!EOP
!BOC
!-----------------------------------------------------------------------
! local variables
!-----------------------------------------------------------------------

  integer (POP_i4) :: n,k ! loop indicies
  character (len = 3) :: char_region ! character version of DYE region number
  character (POP_CharLength) :: &
     appended_long_name, & ! long name with region number appended
     appended_short_name ! short name with region number appended

!---------------------------------------------------------------------------
! prognostic variables
!---------------------------------------------------------------------------

  do n = 1, dye_tracer_cnt
     write(char_region,'(i3)') 100 + n
     appended_short_name = 'DYE'/&
                        &/char_region(2:3)

     appended_long_name = ' Dye Number '/&
                        &/char_region(2:3)
     do k = 1, km
        call define_movie_field(movie_DYE(k,n),trim(appended_short_name),k, &
                             long_name=trim(appended_long_name), &
                             units='years', grid_loc='3111')
     enddo

! now do depth integrated DYEs located in k = 0 id slot

     k = 0
     appended_short_name = 'DYE'/&
                        &/char_region(2:3)/&
                        &/'_zint'
     appended_long_name = 'Depth-Integrated Dye Number '/&
                        &/char_region(2:3)
     call define_movie_field(movie_DYE(k,n),trim(appended_short_name),k, &
                          long_name=trim(appended_long_name), &
                          units='years', grid_loc='2111')
  enddo

!-----------------------------------------------------------------------
!EOC

  end subroutine dye_init_movie

!***********************************************************************

!BOP
! !IROUTINE: dye_movie
! !INTERFACE:

  subroutine dye_movie(bid, TRACER_MODULE)

! !DESCRIPTION:
! set the index bounds of a single passive tracer module
!
! !REVISION HISTORY:
! same as module

! !USES:

  use grid, only : dzr

! !INPUT PARAMETERS:

  integer(POP_i4), intent(in) :: bid ! block id

  real(POP_r8), dimension(nx_block,ny_block,km,dye_tracer_cnt), intent(in) :: &
         TRACER_MODULE

!EOP
!BOC
!-----------------------------------------------------------------------
!
! local variables
!
!-----------------------------------------------------------------------

  integer(POP_i4) :: k, n ! loop indicies

  real(POP_r8), dimension(nx_block,ny_block) :: WORK

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

  do n = 1, dye_tracer_cnt
     do k = 1, km
        if (movie_requested(movie_DYE(k,n))) then
           call update_movie_field(TRACER_MODULE(:,:,k,n), &
                                    movie_DYE(k,n),bid,k)
        endif
     enddo
     if (movie_requested(movie_DYE(0,n))) then ! depth integrated
        WORK = c0
        do k = 2, km ! do not included top level
           WORK = WORK + TRACER_MODULE(:,:,k,n)*dzr(k)
        enddo
        call update_movie_field(WORK,movie_DYE(0,n),bid,0)
     endif
  enddo

!-----------------------------------------------------------------------
!EOC

 end subroutine dye_movie

!***********************************************************************

 end module dye_mod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
