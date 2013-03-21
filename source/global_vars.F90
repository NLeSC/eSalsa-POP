!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module global_vars

!BOP
! !MODULE: global_vars
! !DESCRIPTION:
!  This module contains all the global variables used by GPU extensions
!  in POP.
!
! !REVISION HISTORY:
!  CVS:$Id: prognostic.F90,v 1.5 2013/03/21 22:11:40 ben Exp $
!  CVS:$Name: POP_2_0_1 $

! !USES:

   use kinds_mod
   use blocks
   use domain_size
   use domain
   use constants

   implicit none
   public
   save

! !PUBLIC DATA TYPES:


! !PUBLIC DATA MEMBERS:

   real (r8), dimension(:,:,:,:,:), pointer :: &
      VDC, &             ! tracer diffusivity - public to allow possible modification by Gent-McWilliams horizontal mixing parameterization
      VDCREF             ! copy for GPU result verification

   real (r8), dimension(:,:,:,:), pointer :: &
      VVC,         &! momentum viscosity
      VVCREF        ! copy for GPU result verification


!EOP
!BOC
!EOC
!***********************************************************************


 end module global_vars

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
