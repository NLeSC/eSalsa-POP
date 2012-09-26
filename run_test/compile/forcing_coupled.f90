!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module forcing_coupled

!BOP
!MODULE: forcing_coupled

! !DESCRIPTION:
! This module contains all the routines necessary for coupling POP to
! atmosphere and sea ice models using the NCAR CCSM flux coupler. To
! enable the routines in this module, the coupled ifdef option must
! be specified during the make process.
!
! !REVISION HISTORY:
! SVN:$Id: forcing_coupled.F90 17212 2009-07-20 23:01:42Z njn01 $
!
! !USES:

   use POP_KindsMod
   use POP_ErrorMod
   use POP_CommMod
   use POP_FieldMod
   use POP_GridHorzMod
   use POP_HaloMod

   use kinds_mod
   use blocks, only: nx_block, ny_block, block, get_block
   use domain_size
   use domain
   use io_types, only: stdout, nml_in

   use communicate
   use global_reductions
   use constants
   use io
   use time_management
   use grid
   use prognostic
   use exit_mod
   use ice, only: tfreez, tmelt, liceform,QFLUX, QICE, AQICE, tlast_ice
   use forcing_shf
   use forcing_sfwf
   use forcing_ws, only: ws_data_type
   use forcing_fields
   use timers
   implicit none
   save
!EOP
!BOC
!-----------------------------------------------------------------------
!
! module variables
!
!-----------------------------------------------------------------------
!EOC
!***********************************************************************
 contains
!***********************************************************************
!BOP
! !IROUTINE: pop_init_coupled
! !INTERFACE:
 subroutine pop_init_coupled
! !DESCRIPTION:
! This routine sets up everything necessary for coupling with
! the CCSM3 flux coupler, version 6 (cpl6). Initial cpl information
! is received from the coupler, but pop information is not sent
! from this routine
!
! !REVISION HISTORY:
! same as module
!EOP
!BOC
 end subroutine pop_init_coupled
!***********************************************************************
!BOP
! !IROUTINE: pop_init_partially_coupled
! !INTERFACE:
 subroutine pop_init_partially_coupled
! !DESCRIPTION:
! This routine initializes and allocates arrays for the partially-coupled
! option
!
! !REVISION HISTORY:
! same as module
!EOP
!BOC
 end subroutine pop_init_partially_coupled
!***********************************************************************
!BOP
! !IROUTINE: pop_set_coupled_forcing
! !INTERFACE:
 subroutine pop_set_coupled_forcing
! !DESCRIPTION:
! This routine is called immediately following the receipt of fluxes
! from the coupler. It combines fluxes received from the coupler into
! the STF array and converts from W/m**2 into model units. It also
! balances salt/freshwater in marginal seas and sets SHF_QSW_RAW
! and SHF_COMP. Compute QSW_COSZ_WGHT_NORM if needed.
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
!-----------------------------------------------------------------------
!EOC
 end subroutine pop_set_coupled_forcing
!***********************************************************************
!BOP
! !IROUTINE: set_combined_forcing
! !INTERFACE:
 subroutine set_combined_forcing (STF,FW,TFW)
! !DESCRIPTION:
!
! This routine combines heat flux components into the STF array and
! converts from W/m**2, then combines terms when the "partially-coupled"
! has been selected
!
! !REVISION HISTORY:
! same as module
! !INPUT/OUTPUT PARAMETERS:
   real (r8), dimension(nx_block,ny_block,nt,max_blocks_clinic), &
      intent(inout) :: &
      STF, &! surface tracer fluxes at current timestep
      TFW ! tracer concentration in water flux
   real (r8), dimension(nx_block,ny_block,max_blocks_clinic), &
      intent(inout) :: &
      FW ! fresh water flux
!EOP
!BOC
!-----------------------------------------------------------------------
!
! local variables
!
!-----------------------------------------------------------------------
  integer (int_kind) :: &
     iblock, &! local address of current block
     n ! index
!-----------------------------------------------------------------------
!EOC
 end subroutine set_combined_forcing
!***********************************************************************
!BOP
! !IROUTINE: tavg_coupled_forcing
! !INTERFACE:
 subroutine tavg_coupled_forcing
! !DESCRIPTION:
! This routine accumulates tavg diagnostics related to forcing_coupled
! forcing.
!
! !REVISION HISTORY:
! same as module
!EOP
!BOC
!-----------------------------------------------------------------------
!EOC
 end subroutine tavg_coupled_forcing
!***********************************************************************
!BOP
! !IROUTINE: update_ghost_cells_coupler_fluxes
! !INTERFACE:
   subroutine update_ghost_cells_coupler_fluxes(errorCode)
! !DESCRIPTION:
! This routine accumulates tavg diagnostics related to forcing_coupled
! forcing.
!
! !REVISION HISTORY:
! same as module
! !OUTPUT PARAMETERS:
   integer (POP_i4), intent(out) :: errorCode
!EOP
!BOC
!-----------------------------------------------------------------------
!
! update halos for all coupler fields
!
!-----------------------------------------------------------------------
   errorCode = POP_Success
!-----------------------------------------------------------------------
!EOC
 end subroutine update_ghost_cells_coupler_fluxes
!***********************************************************************
!BOP
! !IROUTINE: rotate_wind_stress
! !INTERFACE:
   subroutine rotate_wind_stress (WORK1,WORK2)
! !DESCRIPTION:
! This subroutine rotates true zonal/meridional wind stress into local
! coordinates, converts to dyne/cm**2, and shifts SMFT to the U grid
!
! !REVISION HISTORY:
! same as module
! !INPUT PARAMETERS:
   real (r8), dimension(nx_block,ny_block,max_blocks_clinic), intent(in) :: &
      WORK1, WORK2 ! contains taux and tauy from coupler
!EOP
!BOC
!-----------------------------------------------------------------------
!EOC
 end subroutine rotate_wind_stress
!***********************************************************************
!***********************************************************************
 end module forcing_coupled
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
