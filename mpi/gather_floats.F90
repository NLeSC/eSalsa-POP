!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!BOP
! !MODULE: gather_floats

 module gather_floats

! !DESCRIPTION:
!  This module contains routines for gathering data to a single
!  processor from a distributed array, specifically for floats.
!
! !REVISION HISTORY:
!  CVS: $Id:  $
!  CVS: $Name:  $

! !USES:

   use kinds_mod
   use communicate
   use constants
   use blocks
   use distribution
   use domain
   use domain_size
   use exit_mod
   use timers, only : timer_start, timer_stop

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: gather_float_fields, float_test_sum

   integer (int_kind), public ::  &
      gather_float_timer

!EOP
!BOC

!-----------------------------------------------------------------------
!
!  module variables
!
!-----------------------------------------------------------------------

!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
!***********************************************************************

 subroutine gather_float_fields(float_fields_global, float_fields, &
    float_on_this_block, dst_task, src_dist, float_special_value, lpositions)

!-----------------------------------------------------------------------
!
!  This subroutine gathers float information from all blocks to
!  the processor dst_task.
!
!-----------------------------------------------------------------------

   use blocks

   include 'mpif.h'

!-----------------------------------------------------------------------
!
!  input variables
!
!-----------------------------------------------------------------------

   integer (int_kind), intent(in) :: &
     dst_task       ! task to which array should be gathered

   type (distrb) :: &
     src_dist       ! distribution of blocks in the source array

   real (r4), dimension(:,:,:) :: &
     float_fields          ! array containing distributed field

   real (r4), intent(in) :: &
     float_special_value  !  special value indicates float not on this block

   logical (log_kind), intent(in) :: &
     lpositions

   logical (log_kind), dimension(:,:), intent(in) :: &
     float_on_this_block

!-----------------------------------------------------------------------
!
!  output variables
!
!-----------------------------------------------------------------------

   real (r4), dimension(:,:), intent(inout) :: &
     float_fields_global        ! array containing global field on dst_task

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: &
     num_floats   ,&! first dimension of the input array
     nblocks        ,&! second dimension of the input array
     num_fields     ,&! third dimension of the input array
     m,iblock,field  ,&! dummy loop counters
     i_local,j_local,&! local indices
     nsends         ,&! number of actual sends
     src_block      ,&! block locator for send
     ierr             ! MPI error flag

   integer (int_kind), dimension(MPI_STATUS_SIZE) :: &
     status

   integer (int_kind), dimension(:), allocatable :: &
     snd_request

   integer (int_kind), dimension(:,:), allocatable :: &
     snd_status

   real (r4) :: &
      del_i, del_j

   real (r4), dimension(:,:), allocatable :: &
     msg_buffer

   type (block) :: &
     this_block ! block info for current block

!-----------------------------------------------------------------------
!
!  if this task is the dst_task, copy local blocks into the global 
!  array and post receives for non-local blocks.
!
!-----------------------------------------------------------------------

   call timer_start(gather_float_timer)

   float_fields_global = float_special_value

   num_floats = SIZE(float_fields,DIM=1)
   nblocks      = SIZE(float_fields,DIM=2)
   num_fields   = SIZE(float_fields,DIM=3)

   if(lpositions) then
      do iblock = 1, nblocks
         this_block = get_block(blocks_clinic(iblock),iblock)
         do m=1,num_floats
            if (float_on_this_block(m,iblock) ) then
!maltrud debug
!write(*,*)' converting local indices: ',float_fields(m,iblock,1:2)
               i_local = int(float_fields(m,iblock,1)) 
               j_local = int(float_fields(m,iblock,2)) 
               del_i = float_fields(m,iblock,1) - i_local
               del_j = float_fields(m,iblock,2) - j_local
               float_fields(m, iblock, 1) =   &
                  this_block%i_glob(i_local) + del_i
               float_fields(m, iblock, 2) =   &
                  this_block%j_glob(j_local) + del_j
               if (float_fields(m, iblock, 1) > float(nx_global) )  &
                   float_fields(m, iblock, 1) =   &
                   float_fields(m, iblock, 1) - float(nx_global)
               if (float_fields(m, iblock, 1) < 0.0_r4 )  &
                   float_fields(m, iblock, 1) =   &
                   float_fields(m, iblock, 1) + float(nx_global)
!maltrud debug
!write(*,*)' to global indices: ',float_fields(m,iblock,1:2)
            endif
         enddo
      enddo
   endif

   do iblock = 1, nblocks
      do m=1,num_floats
         if (.not. float_on_this_block(m,iblock) )   &
            float_fields(m, iblock, 1:num_fields) = float_special_value
      enddo
   enddo

   if (my_task == dst_task) then

     do iblock=1,nblocks_tot

       !*** copy local blocks

       if (src_dist%proc(iblock) == my_task+1) then

         src_block = src_dist%local_block(iblock)

         do field=1,num_fields
         do m=1,num_floats
           if(float_on_this_block(m,src_block) )   &
              float_fields_global(m,field) = float_fields(m,src_block,field) 
!maltrud debug
!write(*,*)'*',my_task,block,src_block,m,field,float_fields(m,src_block,field),float_fields_global(m,field)
         end do
         end do

       endif

     end do

     !*** receive blocks to fill up the rest

     allocate (msg_buffer(num_floats,num_fields))

     do iblock=1,nblocks_tot
       if (src_dist%proc(iblock) > 0 .and. &
           src_dist%proc(iblock) /= my_task+1) then

         call MPI_RECV(msg_buffer, size(msg_buffer), &
                       mpi_real, src_dist%proc(iblock)-1, 7*mpitag_gs+iblock, &
                       src_dist%communicator, status, ierr)

         do field=1,num_fields
         do m=1,num_floats
           if(msg_buffer(m,field) /= float_special_value)   &
              float_fields_global(m,field) = msg_buffer(m,field) 
!maltrud debug
!write(*,*)my_task,block,src_block,m,field,msg_buffer(m,field),float_fields_global(m,field)
         end do
         end do
       endif
     end do

     deallocate(msg_buffer)

!-----------------------------------------------------------------------
!
!  otherwise send data to dst_task
!
!-----------------------------------------------------------------------

   else

     allocate(snd_request(nblocks_tot), &
              snd_status (MPI_STATUS_SIZE, nblocks_tot))
     allocate (msg_buffer(num_floats,num_fields))
     msg_buffer = float_special_value

     nsends = 0
     do iblock=1,nblocks_tot
       if (src_dist%proc(iblock) == my_task+1) then

         src_block = src_dist%local_block(iblock)

         do field=1,num_fields
         do m=1,num_floats
!maltrud debug
         if(float_on_this_block(m,src_block))  &
            msg_buffer(m,field) = float_fields(m,src_block,field)
!if(my_task == 3)  &
!write(*,*)my_task,block,src_block,m,field,float_fields(m,src_block,field)
         end do
         end do

         nsends = nsends + 1
         call MPI_ISEND(msg_buffer, num_floats*num_fields, &
                     mpi_real, dst_task, 7*mpitag_gs+iblock, &
                     src_dist%communicator, snd_request(nsends), ierr)
       endif
     end do

     if (nsends > 0) &
       call MPI_WAITALL(nsends, snd_request, snd_status, ierr)
     deallocate(snd_request, snd_status)
     deallocate(msg_buffer)

   endif

   call timer_stop(gather_float_timer)

!-----------------------------------------------------------------------

 end subroutine gather_float_fields

!EOC
!***********************************************************************
!BOP
! !IROUTINE: float_test_sum
! !INTERFACE:

!maltrud should be able to replace this with a function in global_reductions soon

   function float_test_sum(local_int_sum,src_dist)

   include 'mpif.h'   ! MPI Fortran include file

   type (distrb) :: &
     src_dist       ! distribution of blocks in the source array

   integer (int_kind) ::    &
      float_test_sum, local_int_sum, ierr

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

   call MPI_ALLREDUCE(local_int_sum, float_test_sum, 1,  &
                      MPI_INTEGER, MPI_SUM, src_dist%communicator, ierr)

!-----------------------------------------------------------------------
   end function float_test_sum

!EOC
!***********************************************************************

!***********************************************************************

 end module gather_floats

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
