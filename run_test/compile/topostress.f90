!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module topostress

!BOP
! !MODULE: topostress

! !DESCRIPTION:
! This module contains routines necessary for computing stress
! due to bottom topography.
!
! !REVISION HISTORY:
! CVS:$Id: topostress.F90,v 1.10 2002/12/02 13:45:11 pwjones Exp $
! CVS:$Name: $

! !USES:

   use POP_KindsMod
   use POP_ErrorMod
   use POP_HaloMod
   use POP_GridHorzMod
   use POP_FieldMod

   use domain
   use blocks
   use distribution
   use POP_DomainSizeMod, only : POP_nxGlobal, POP_nyGlobal
   use constants
   use io
   use grid
   use broadcast

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: init_topostress

! !PUBLIC DATA MEMBERS:

   logical (POP_logical), public :: &
     ltopostress ! true if topographic stress desired

   real (POP_r8), dimension(:,:,:), allocatable, public :: &
     TSU, TSV ! topographic stress velocities

!EOP
!BOC
!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: init_topostress
! !INTERFACE:

 subroutine init_topostress(errorCode)

! !DESCRIPTION:
! This routine allocates stress arrays if topographic stress is
! chosen and initializes topo stress parameters.
!
! !REVISION HISTORY:
! same as module
!
! !OUTPUT PARAMETERS:

   integer(POP_i4), intent(out) :: &
      errorCode ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
! input namelist to choose topostress
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      nsmooth_topo ! number of passes of topography smoother

   integer (POP_i4) :: nml_error ! namelist i/o error flag

   character (POP_charLength) :: &
      ustar_opt, &! option for specifying ustar
      ustar_file, &! input file for reading ustar
      ustar_file_fmt ! input file format for reading ustar

   namelist /topostress_nml/ltopostress, nsmooth_topo, &
                            ustar_opt, ustar_file, ustar_file_fmt

!-----------------------------------------------------------------------
!
! read namelist to see if topostress desired
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   ltopostress = .false.
   nsmooth_topo= 0
   ustar_opt = 'internal'
   ustar_file = 'unknown_ustar_file'
   ustar_file_fmt = 'bin'

   if (my_task == master_task) then
      open (nml_in, file=nml_filename, status='old',iostat=nml_error)
      if (nml_error /= 0) then
         nml_error = -1
      else
         nml_error = 1
      endif
      do while (nml_error > 0)
         read(nml_in, nml=topostress_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
   endif

   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
      call POP_ErrorSet(errorCode, &
         'POP_TopostressInit: error reading configuration info')
      return
   endif

   if (my_task == master_task) then
      write(stdout,blank_fmt)
      write(stdout,delim_fmt)
      write(stdout,blank_fmt)
      write(stdout,'(a27)') ' Topographic stress options'
      write(stdout,blank_fmt)
      write(stdout,delim_fmt)
      if (ltopostress) then
         write(stdout,'(a27)') ' Topographic stress enabled'
         if (trim(ustar_opt) == 'internal') then
            write(stdout,'(a27)') ' USTAR generated internally'
            if (nsmooth_topo > 0) then
               write(stdout,'(a26,i2,a8)') ' Topography smoothed with ', &
                                           nsmooth_topo,' passes.'
            else
               write(stdout,'(a25)') ' Topography not smoothed.'
            endif
         else ! ustar read in from file
            write(stdout,*)' USTAR to be read in from file: ', &
               trim(ustar_file)
         endif
      else
         write(stdout,'(a28)') ' Topographic stress disabled'
      endif
   endif

   call broadcast_scalar(ltopostress, master_task)
   call broadcast_scalar(nsmooth_topo, master_task)
   call broadcast_scalar(ustar_opt, master_task)
   call broadcast_scalar(ustar_file, master_task)
   call broadcast_scalar(ustar_file_fmt, master_task)

!-----------------------------------------------------------------------
!
! allocate the topographic stress velocity arrays if required
!
!-----------------------------------------------------------------------

   if (ltopostress) then
      allocate (TSU(nx_block,ny_block,nblocks_clinic), &
                TSV(nx_block,ny_block,nblocks_clinic))

      select case (ustar_opt)
      case ('internal')

         call topo_stress_internal(nsmooth_topo, errorCode)

         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_TopostressInit: error computing internal topostress')
            return
         endif

      case ('file')

         call read_ustar(ustar_file, ustar_file_fmt, errorCode)

         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'POP_TopostressInit: error reading ustar')
            return
         endif

      case default

         call POP_ErrorSet(errorCode, &
            'POP_TopostressInit: unknown topostress ustar option')
         return

      end select

   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine init_topostress

!***********************************************************************
!BOP
! !IROUTINE: topo_stress_internal
! !INTERFACE:

 subroutine topo_stress_internal(nsmooth_topo, errorCode)

! !DESCRIPTION:
! Calculate topographic stress (maximum entropy) velocities. These
! are time-independent 2d fields given by:
! \begin{eqnarray}
! \Psi^* &=& -f L^2 H \! H U^* &=& -\nabla_y(\Psi^*) \! H V^* &=& +\nabla_x(\Psi^*)


! \end{eqnarray}
!
! !REVISION HISTORY:
! same as module

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
      nsmooth_topo ! number of passes of topography smoother

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
! local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      i,j,iter,iblock, &! local iteration counters
      blkError ! block level error code

   real (POP_r8), dimension(nx_block,ny_block) :: &
      TSP, &! topo stress streamfunction
      SCALE, &! scale length
      HTOLD ! old topography

   real (POP_r8), dimension(nx_block,ny_block,nblocks_clinic) :: &
      HTNEW ! smoothed topography

   real (POP_r8), parameter :: &
      tslse = 12.0e5_POP_r8, &!
      tslsp = 3.0e5_POP_r8 !

   type (block) :: &
      this_block ! block information for current block

!-----------------------------------------------------------------------
!
! smooth topography if requested
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   do iblock = 1,nblocks_clinic
      HTNEW(:,:,iblock) = HT(:,:,iblock) ! initialize
   end do

   do iter = 1, nsmooth_topo

      do iblock = 1,nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)

         HTOLD = HTNEW(:,:,iblock)

         call smooth_topo2(HTOLD,HTNEW(:,:,iblock),this_block, blkError)

         if (blkError /= POP_Success) errorCode = POP_Fail

      end do

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_TopostressInternal: error in smoothing')
         return
      endif

      call POP_HaloUpdate(HTNEW, POP_haloClinic, &
                          POP_gridHorzLocCenter, POP_fieldKindScalar, &
                          errorCode, fillValue = 0.0_POP_r8)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_TopostressInternal: error updating htNew')
         return
      endif
   enddo

!-----------------------------------------------------------------------
!
! calculate the topographic stress equilibrium stream function.
!
!-----------------------------------------------------------------------

   do iblock = 1,nblocks_clinic
      this_block = get_block(blocks_clinic(iblock),iblock)

      SCALE = tslsp + (tslse - tslsp)* &
                      (p5 + p5*cos(c2*TLAT(:,:,iblock)))

      where (KMT(:,:,iblock) /= 0)
         TSP = -FCORT(:,:,iblock)*SCALE*SCALE*HTNEW(:,:,iblock)
      elsewhere
         TSP = c0
      end where

!-----------------------------------------------------------------------
!
! calculate the topographic stress velocities from stream function.
! compute gradient with 4 point stencil
!
!-----------------------------------------------------------------------

      TSU(:,:,iblock) = c0
      TSV(:,:,iblock) = c0

      do j=this_block%jb,this_block%je
      do i=this_block%ib,this_block%ie
         TSV(i,j,iblock) = DXUR(i,j,iblock)*p5*HUR(i,j,iblock)* &
                            (TSP(i+1,j+1) - TSP(i ,j) - &
                             TSP(i ,j+1) + TSP(i+1,j))
         TSU(i,j,iblock) = -DYUR(i,j,iblock)*p5*HUR(i,j,iblock)* &
                            (TSP(i+1,j+1) - TSP(i ,j) + &
                             TSP(i ,j+1) - TSP(i+1,j))
      end do
      end do

      where (KMU(:,:,iblock) == 0)
         TSV(:,:,iblock) = c0 ! zero at land points
         TSU(:,:,iblock) = c0
      endwhere

      ! apply only in 'deep' water
      ! where (KMU(:,:,iblock) <= 3)
      ! TSU(:,:,iblock) = c0
      ! TSV(:,:,iblock) = c0
      ! endwhere

   end do ! block loop

   call POP_HaloUpdate(TSU, POP_haloClinic, &
                       POP_gridHorzLocNECorner, POP_fieldKindVector, &
                       errorCode, fillValue = 0.0_POP_r8)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_TopostressInternal: error updating TSU')
      return
   endif

   call POP_HaloUpdate(TSV, POP_haloClinic, &
                       POP_gridHorzLocNECorner, POP_fieldKindVector, &
                       errorCode, fillValue = 0.0_POP_r8)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'POP_TopostressInternal: error updating TSV')
      return
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine topo_stress_internal

!***********************************************************************
!BOP
! !IROUTINE: smooth_topo2
! !INTERFACE:

 subroutine smooth_topo2(HTOLD,HTNEW,this_block, errorCode)

! !DESCRIPTION:
! This routine smooths topography using a 9-point averaging stencil
! given by
! \begin{equation}
! \begin{array}{ccccc}
! 1 & -- & 2 & -- & 1 \! | & & | & & | \! 2 & -- & 4 & -- & 2 \! | & & | & & | \! 1 & -- & 2 & -- & 1




! \end{array} \nonumber
! \end{equation}
! Land points are not included in the smoothing, and the
! stencil is modified to include only ocean points in the
! averaging. This routine is nearly identical to the smooth
! topography routine in the grid module.
!
! !REVISION HISTORY:
! same as module

! !INPUT PARAMETERS:

   real (POP_r8), dimension(nx_block,ny_block), intent(in) :: &
      HTOLD ! old HT field to be smoothed

   type (block), intent(in) :: &
      this_block ! block info for this sub block

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode ! returned error code

   real (POP_r8), dimension(nx_block,ny_block), intent(out) :: &
      HTNEW ! smoothed HT field

!EOP
!BOC
!-----------------------------------------------------------------------
!
! local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      i,j, &! loop counters
      bid ! local block index

   integer (POP_i4), dimension(nx_block,ny_block) :: &
      NB, &! array to compute number of ocean neighbors
      IWORK ! local work space

   real (POP_r8), dimension(nx_block,ny_block) :: &
      WORK ! local work space

!-----------------------------------------------------------------------
!
! smooth topography
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   bid = this_block%local_id

   where (KMT(:,:,bid) > 0)
      IWORK = 1
      HTNEW = HTOLD
   elsewhere
      IWORK = 0
      HTNEW = c0
   endwhere

   do j=this_block%jb,this_block%je
   do i=this_block%ib,this_block%ie

      WORK(i,j) = c4*HTNEW(i,j) + &
                  c2*HTNEW(i+1,j ) + c2*HTNEW(i-1,j ) + &
                  c2*HTNEW(i ,j+1) + c2*HTNEW(i ,j-1) + &
                     HTNEW(i+1,j+1) + HTNEW(i+1,j-1) + &
                     HTNEW(i-1,j+1) + HTNEW(i-1,j-1)

      NB(i,j) = c4*IWORK(i,j) + &
                c2*IWORK(i+1,j ) + c2*IWORK(i-1,j ) + &
                c2*IWORK(i ,j+1) + c2*IWORK(i ,j-1) + &
                   IWORK(i+1,j+1) + IWORK(i+1,j-1) + &
                   IWORK(i-1,j+1) + IWORK(i-1,j-1)

   end do
   end do

!-----------------------------------------------------------------------
!
! new depth field
!
!-----------------------------------------------------------------------

   where ((KMT(:,:,bid) /= 0) .and. (NB /= 0))
      HTNEW = WORK/real(NB)
   elsewhere
      HTNEW = c0
   endwhere

!-----------------------------------------------------------------------
!EOC

 end subroutine smooth_topo2

!***********************************************************************

!***********************************************************************
!BOP
! !IROUTINE: read_ustar
! !INTERFACE:

 subroutine read_ustar(ustar_file, ustar_file_fmt, errorCode)

! !DESCRIPTION:
! Reads ustar from input file
!
! !REVISION HISTORY:
! same as module

! !INPUT PARAMETERS:

   character (*), intent(in) :: &
      ustar_file, &! filename of file containing ustar
      ustar_file_fmt ! file format of file containing ustar

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
! local variables
!
!-----------------------------------------------------------------------

   type (datafile) :: &
      ustar_data_file ! io file type for ustar file

   type (io_field_desc) :: &
      io_ustar, io_vstar ! descriptors for USTAR

   type (io_dim) :: &
      i_dim, j_dim ! dimension descriptors for horiz dims

!-----------------------------------------------------------------------
!
! construct io file type and open for reading
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   ustar_data_file = construct_file(ustar_file_fmt, full_name=ustar_file, &
                                     record_length=rec_type_dbl, &
                                     recl_words=POP_nxGlobal*POP_nyGlobal)

   call data_set(ustar_data_file, 'open_read')

!-----------------------------------------------------------------------
!
! define variables to be read
!
!-----------------------------------------------------------------------

   !*** define dimensions

   i_dim = construct_io_dim('i', POP_nxGlobal)
   j_dim = construct_io_dim('j', POP_nyGlobal)

   io_ustar = construct_io_field('USTAR', dim1=i_dim, dim2=j_dim, &
                   long_name='topostress ustar', &
                   units ='cm/s', &
                   grid_loc ='2221', &
                   field_loc = field_loc_NEcorner, &
                   field_type = field_type_vector, &
                   d2d_array = TSU)

   io_vstar = construct_io_field('VSTAR', dim1=i_dim, dim2=j_dim, &
                   long_name='topostress vstar', &
                   units ='cm/s', &
                   grid_loc ='2221', &
                   field_loc = field_loc_NEcorner, &
                   field_type = field_type_vector, &
                   d2d_array = TSV)

   call data_set (ustar_data_file, 'define', io_ustar)
   call data_set (ustar_data_file, 'define', io_vstar)

!-----------------------------------------------------------------------
!
! read arrays then clean up
!
!-----------------------------------------------------------------------

   call data_set (ustar_data_file, 'read', io_ustar)
   call data_set (ustar_data_file, 'read', io_vstar)

   call destroy_io_field (io_ustar)
   call destroy_io_field (io_vstar)

   call data_set (ustar_data_file, 'close')
   call destroy_file(ustar_data_file)

!-----------------------------------------------------------------------
!EOC

 end subroutine read_ustar

!***********************************************************************

 end module topostress

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
