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
      VDC,         &! pointer used to point to pinned memory allocation
      VDC_ALLOC     ! pointer used for allocation

   real (r8), dimension(:,:,:,:), pointer :: &
      VVC,         &! momentum viscosity
      VVCREF,      &! copy for GPU result verification
      VDCREF        ! copy for GPU result verification

    real (r8), dimension(:,:,:), pointer :: &
      DBLOC,      &! buoyancy difference between adjacent levels
      DBSFC,      &! buoyancy difference between level and surface
      DBLOCREF,   &! for verification
      DBSFCREF     ! for verification

   real (r8), dimension(:,:,:), pointer :: &
      VISC,       &! local temp for viscosity
      VISCREF      ! for verification

   real (r8), dimension(:,:,:), pointer :: &
      SHF_QSW      ! incoming short wave

   real (r8), dimension(:,:,:,:), pointer :: &
      SMF,    &!  surface momentum fluxes (wind stress)
      SMFT     !  surface momentum fluxes on T points if avail

   real (r8), dimension(:,:,:,:), pointer :: &
      STF          !  surface tracer fluxes

   real (r8), dimension(:,:,:), pointer :: &
      KPP_HBLT     ! boundary layer depth

   real (r8), dimension(:,:), pointer :: &
      HBLTREF,    &! for verification
      BFSFC,       &! surface buoyancy forcing
      BFSFCREF,    &! surface buoyancy forcing
      USTAR,       &! surface friction velocity
      USTARREF,    &! surface friction velocity
      STABLE,      &! = 1 for stable forcing; = 0 for unstable forcing
      STABLEREF     ! = 1 for stable forcing; = 0 for unstable forcing

   integer (int_kind), dimension(:,:), pointer :: &
      KBL,        &! index of first lvl below hbl
      KBLREF       ! for verification

   real (r8), dimension(:,:,:), pointer :: &
      GHAT,        &! non-local mixing coefficient
      GHATREF       ! for verification


!EOP
!BOC
!EOC
!***********************************************************************


 end module global_vars

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
