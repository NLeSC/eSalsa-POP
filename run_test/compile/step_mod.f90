!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module step_mod

!BOP
! !MODULE: step_mod

! !DESCRIPTION:
! Contains the routine for stepping the model forward one timestep
!
! !REVISION HISTORY:
! SVN:$Id: step_mod.F90 16641 2009-06-12 17:00:31Z njn01 $
!
! !USES:

   use POP_KindsMod
   use POP_ErrorMod
   use POP_CommMod
   use POP_FieldMod
   use POP_GridHorzMod
   use POP_HaloMod
   use POP_BlocksMod

   use POP_DomainSizeMod, only : POP_km, POP_nt, POP_maxBlocksClinic
   use blocks
   use domain
   use constants
   use prognostic
   use timers
   use grid
   use diagnostics
   use state_mod, only: state
   use time_management
   use baroclinic
   use barotropic
   use surface_hgt
   use tavg
   use forcing_fields
   use forcing
   use ice
   use passive_tracers
   use registry
   use communicate
   use io_types
   use floats, only: float_move_driver, float_global_to_local, &
                       calc_float_field_values, float_ijk_global, &
                       float_ijk, float_on_this_block, lfloat_on, &
                       float_special_value, float_freq_iopt, float_freq, &
                       float_flag
   use gather_floats, only: gather_float_fields, float_test_sum

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: step, init_step

!----------------------------------------------------------------------
!
! module variables
!
!----------------------------------------------------------------------
   integer (POP_i4), private :: &
      tavg_flag ! flag to access tavg frequencies

   integer (POP_i4), private :: &
      timer_step, &! timer number for step
      timer_baroclinic, &! timer for baroclinic parts of step
      timer_barotropic, &! timer for barotropic part of step
      timer_3dupdate ! timer for the 3D update after baroclinic component

   integer (POP_i4) :: ierr
!EOP
!BOC
!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: step
! !INTERFACE:

 subroutine step(errorCode)

! !DESCRIPTION:
! This routine advances the simulation on timestep.
! It controls logic for leapfrog and/or Matsuno timesteps and performs
! time-averaging if necessary. Prognostic variables are updated for
! the next timestep near the end of the routine.
! On Matsuno steps, the time (n) velocity and tracer arrays
! UBTROP,VBTROP,UVEL,VVEL,TRACER contain the predicted new
! velocities from the 1st pass for use in the 2nd pass.
!
! !REVISION HISTORY:
! same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!
! local or common variables:
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      errorCode, istat

   integer (POP_i4) :: &
      i,j,k,n, &! loop indices
      tmptime, &! temp space for time index swapping
      iblock, &! block counter
      ipass, &! pass counter
      num_floats_leaving_block, &!
      num_passes ! number of passes through time step
                          ! (Matsuno requires two)

   real (POP_r8), dimension(nx_block,ny_block,POP_maxBlocksClinic) :: &
      ZX,ZY, &! vertically integrated forcing terms
      DH,DHU ! time change of surface height minus
                          ! freshwater flux at T, U points
   real (POP_r8), allocatable, dimension(:,:) :: &
      FILTER_TEMP ! temporary used with modified Robert filter

   logical (POP_logical), save :: &
      first_call = .true. ! flag for initializing timers

   type (block) :: &
      this_block ! block information for current block






!-----------------------------------------------------------------------
!
! start step timer
!
!-----------------------------------------------------------------------







   call timer_start(timer_step)

   errorCode = POP_Success

!-----------------------------------------------------------------------
!
! read fields for surface forcing
!
!-----------------------------------------------------------------------

   call set_surface_forcing

!-----------------------------------------------------------------------
!
! update timestep counter, set corresponding model time, set
! time-dependent logical switches to determine program flow.
!
!-----------------------------------------------------------------------

   call time_manager(registry_match('lcoupled'), liceform)

!-----------------------------------------------------------------------
!
! compute and initialize some time-average diagnostics
!
!-----------------------------------------------------------------------

   call tavg_set_flag
   call tavg_forcing
   if (nt > 2) call passive_tracers_tavg_sflux(STF)
   call movie_forcing


!-----------------------------------------------------------------------
!
! set timesteps and time-centering parameters for leapfrog or
! matsuno steps.
!
!-----------------------------------------------------------------------

   mix_pass = 0
   if (matsuno_ts) then
      num_passes = 2
   else
      num_passes = 1
   endif


   do ipass = 1,num_passes







      if (matsuno_ts) mix_pass = mix_pass + 1

      if (leapfrogts) then ! leapfrog (and averaging) timestep
         mixtime = oldtime
         beta = alpha
         do k = 1,POP_km
            c2dtt(k) = c2*dt(k)
         enddo
         c2dtu = c2*dtu
         c2dtp = c2*dtp ! barotropic timestep = baroclinic timestep
         c2dtq = c2*dtu ! turbulence timestep = mean flow timestep
      else
         mixtime = curtime
         beta = theta
         do k = 1,POP_km
            c2dtt(k) = dt(k)
         enddo
         c2dtu = dtu
         c2dtp = dtp ! barotropic timestep = baroclinic timestep
         c2dtq = dtu ! turbulence timestep = mean flow timestep
      endif

!-----------------------------------------------------------------------
!
! on 1st pass of matsuno, set time (n-1) variables equal to
! time (n) variables.
!
!-----------------------------------------------------------------------


      if (mix_pass == 1) then

         !$OMP PARALLEL DO PRIVATE(iblock)
         do iblock = 1,nblocks_clinic
            UBTROP(:,:,oldtime,iblock) = UBTROP(:,:,curtime,iblock)
            VBTROP(:,:,oldtime,iblock) = VBTROP(:,:,curtime,iblock)
            UVEL(:,:,:,oldtime,iblock) = UVEL(:,:,:,curtime,iblock)
            VVEL(:,:,:,oldtime,iblock) = VVEL(:,:,:,curtime,iblock)
            RHO (:,:,:,oldtime,iblock) = RHO (:,:,:,curtime,iblock)
            TRACER(:,:,:,:,oldtime,iblock) = &
            TRACER(:,:,:,:,curtime,iblock)
         end do
         !$OMP END PARALLEL DO

      endif


!-----------------------------------------------------------------------
!
! initialize diagnostic flags and sums
!
!-----------------------------------------------------------------------

      call diag_init_sums

!-----------------------------------------------------------------------
!
! calculate change in surface height dh/dt from surface pressure
!
!-----------------------------------------------------------------------

      call dhdt(DH,DHU, errorCode)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'step: error in dhdt')
         return
      endif

!-----------------------------------------------------------------------
!
! Integrate baroclinic equations explicitly to find tracers and
! baroclinic velocities at new time. Update ghost cells for
! forcing terms leading into the barotropic solver.
!
!-----------------------------------------------------------------------







      if(profile_barrier) call POP_Barrier
      call timer_start(timer_baroclinic)
      call baroclinic_driver(ZX,ZY,DH,DHU, errorCode)
      if(profile_barrier) call POP_Barrier
      call timer_stop(timer_baroclinic)
      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'step: error in baroclinic driver')
         return
      endif







      call POP_HaloUpdate(ZX, POP_haloClinic, POP_gridHorzLocNECorner, &
                              POP_fieldKindVector, errorCode, &
                              fillValue = 0.0_POP_r8)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'step: error updating halo for ZX')
         return
      endif

      call POP_HaloUpdate(ZY, POP_haloClinic, POP_gridHorzLocNECorner, &
                              POP_fieldKindVector, errorCode, &
                              fillValue = 0.0_POP_r8)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'step: error updating halo for ZY')
         return
      endif







!-----------------------------------------------------------------------
!
! Solve barotropic equations implicitly to find surface pressure
! and barotropic velocities.
!
!-----------------------------------------------------------------------

      if(profile_barrier) call POP_Barrier
      call timer_start(timer_barotropic)
      call barotropic_driver(ZX,ZY,errorCode)
      if(profile_barrier) call POP_Barrier
      call timer_stop(timer_barotropic)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'POP_Step: error in barotropic solver')
         return
      endif







!-----------------------------------------------------------------------
!
! update tracers using surface height at new time
! also peform adjustment-like physics (convection, ice formation)
!
!-----------------------------------------------------------------------

      call timer_start(timer_baroclinic)
      call baroclinic_correct_adjust(errorCode)
      call timer_stop(timer_baroclinic)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'step: error in baroclinic correct-adjust')
         return
      endif







      if(profile_barrier) call POP_Barrier
      call timer_start(timer_3dupdate)

      call POP_HaloUpdate(UBTROP(:,:,newtime,:), &
                                  POP_haloClinic, &
                                  POP_gridHorzLocNECorner, &
                                  POP_fieldKindVector, errorCode, &
                                  fillValue = 0.0_POP_r8)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'step: error updating halo for UBTROP')
         return
      endif

      call POP_HaloUpdate(VBTROP(:,:,newtime,:), &
                                  POP_haloClinic, &
                                  POP_gridHorzLocNECorner, &
                                  POP_fieldKindVector, errorCode, &
                                  fillValue = 0.0_POP_r8)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'step: error updating halo for VBTROP')
         return
      endif

      call POP_HaloUpdate(UVEL(:,:,:,newtime,:), &
                                POP_haloClinic, &
                                POP_gridHorzLocNECorner, &
                                POP_fieldKindVector, errorCode, &
                                fillValue = 0.0_POP_r8)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'step: error updating halo for UVEL')
         return
      endif

      call POP_HaloUpdate(VVEL(:,:,:,newtime,:), &
                                POP_haloClinic, &
                                POP_gridHorzLocNECorner, &
                                POP_fieldKindVector, errorCode, &
                                fillValue = 0.0_POP_r8)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'step: error updating halo for VVEL')
         return
      endif

      call POP_HaloUpdate(RHO(:,:,:,newtime,:), &
                               POP_haloClinic, &
                               POP_gridHorzLocCenter, &
                               POP_fieldKindScalar, errorCode, &
                               fillValue = 0.0_POP_r8)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'step: error updating halo for RHO')
         return
      endif

      call POP_HaloUpdate(TRACER(:,:,:,:,newtime,:), POP_haloClinic, &
                                  POP_gridHorzLocCenter, &
                                  POP_fieldKindScalar, errorCode, &
                                  fillValue = 0.0_POP_r8)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'step: error updating halo for TRACER')
         return
      endif

      call POP_HaloUpdate(QICE(:,:,:), &
                               POP_haloClinic, &
                               POP_gridHorzLocCenter, &
                               POP_fieldKindScalar, errorCode, &
                               fillValue = 0.0_POP_r8)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'step: error updating halo for QICE')
         return
      endif

      call POP_HaloUpdate(AQICE(:,:,:), &
                               POP_haloClinic, &
                               POP_gridHorzLocCenter, &
                               POP_fieldKindScalar, errorCode, &
                               fillValue = 0.0_POP_r8)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'step: error updating halo for AQICE')
         return
      endif


      if(profile_barrier) call POP_Barrier
      call timer_stop(timer_3dupdate)







!-----------------------------------------------------------------------
!
! add barotropic to baroclinic velocities at new time
!
!-----------------------------------------------------------------------

      !$OMP PARALLEL DO PRIVATE(iblock,k,i,j)
      do iblock = 1,nblocks_clinic

!CDIR NOVECTOR
         do k=1,POP_km
            do j=1,ny_block
            do i=1,nx_block
               if (k <= KMU(i,j,iblock)) then
                  UVEL(i,j,k,newtime,iblock) = &
                  UVEL(i,j,k,newtime,iblock) + UBTROP(i,j,newtime,iblock)
                  VVEL(i,j,k,newtime,iblock) = &
                  VVEL(i,j,k,newtime,iblock) + VBTROP(i,j,newtime,iblock)
               endif
            enddo
            enddo
         enddo

!-----------------------------------------------------------------------
!
! on matsuno mixing steps update variables and cycle for 2nd pass
! note: first step is forward only.
!
!-----------------------------------------------------------------------

         if (mix_pass == 1) then

            UBTROP(:,:,curtime,iblock) = UBTROP(:,:,newtime,iblock)
            VBTROP(:,:,curtime,iblock) = VBTROP(:,:,newtime,iblock)
            UVEL(:,:,:,curtime,iblock) = UVEL(:,:,:,newtime,iblock)
            VVEL(:,:,:,curtime,iblock) = VVEL(:,:,:,newtime,iblock)
            RHO (:,:,:,curtime,iblock) = RHO (:,:,:,newtime,iblock)
            TRACER(:,:,:,:,curtime,iblock) = &
            TRACER(:,:,:,:,newtime,iblock)

         endif
      enddo ! block loop
      !$OMP END PARALLEL DO







   end do ! ipass: cycle for 2nd pass in matsuno step







!-----------------------------------------------------------------------
!
! extrapolate next guess for pressure from three known time levels
!
!-----------------------------------------------------------------------

   !$OMP PARALLEL DO PRIVATE(iblock)
   do iblock = 1,nblocks_clinic
      PGUESS(:,:,iblock) = c3*(PSURF(:,:,newtime,iblock) - &
                               PSURF(:,:,curtime,iblock)) + &
                               PSURF(:,:,oldtime,iblock)
   end do
   !$OMP END PARALLEL DO







!-----------------------------------------------------------------------
!
! compute some global diagnostics
! before updating prognostic variables
!
!-----------------------------------------------------------------------

   call diag_global_preupdate(DH,DHU)

!-----------------------------------------------------------------------
!
! if floats are active, advect them forward, then check to see if
! any have left their home block. if yes, then gather-scatter the
! positions.
! if it is time for output, interpolate properties to the current
! float positions. actual output is done in routine output.
!
!-----------------------------------------------------------------------

   if (lfloat_on) then
      call float_move_driver(num_floats_leaving_block)
      i = float_test_sum(num_floats_leaving_block,distrb_clinic)
      if (i > 0) then
         call gather_float_fields(float_ijk_global, float_ijk, &
                                    float_on_this_block, master_task, &
                                    distrb_clinic, float_special_value,&
                                    .true.)
         call broadcast_array(float_ijk_global, master_task)
         call float_global_to_local
      endif

      if (check_time_flag(float_flag)) &
         call calc_float_field_values

   endif

!-----------------------------------------------------------------------
!
! update prognostic variables for next timestep:
! on normal timesteps
! (n) -> (n-1)
! (n+1) -> (n)
! on averaging timesteps
! [(n) + (n-1)]/2 -> (n-1)
! [(n+1) + (n)]/2 -> (n)
!
!-----------------------------------------------------------------------







   if (avg_ts .or. back_to_back) then ! averaging step

      !$OMP PARALLEL DO PRIVATE(iblock,this_block,k,n)

      do iblock = 1,nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)

         !*** avg 2-d fields

         UBTROP(:,:,oldtime,iblock) = p5*(UBTROP(:,:,oldtime,iblock) + &
                                          UBTROP(:,:,curtime,iblock))
         VBTROP(:,:,oldtime,iblock) = p5*(VBTROP(:,:,oldtime,iblock) + &
                                          VBTROP(:,:,curtime,iblock))
         UBTROP(:,:,curtime,iblock) = p5*(UBTROP(:,:,curtime,iblock) + &
                                          UBTROP(:,:,newtime,iblock))
         VBTROP(:,:,curtime,iblock) = p5*(VBTROP(:,:,curtime,iblock) + &
                                          VBTROP(:,:,newtime,iblock))
         GRADPX(:,:,oldtime,iblock) = p5*(GRADPX(:,:,oldtime,iblock) + &
                                          GRADPX(:,:,curtime,iblock))
         GRADPY(:,:,oldtime,iblock) = p5*(GRADPY(:,:,oldtime,iblock) + &
                                          GRADPY(:,:,curtime,iblock))
         GRADPX(:,:,curtime,iblock) = p5*(GRADPX(:,:,curtime,iblock) + &
                                          GRADPX(:,:,newtime,iblock))
         GRADPY(:,:,curtime,iblock) = p5*(GRADPY(:,:,curtime,iblock) + &
                                          GRADPY(:,:,newtime,iblock))
         FW_OLD(:,:,iblock) = p5*(FW(:,:,iblock) + FW_OLD(:,:,iblock))

         !*** avg 3-d fields

         UVEL(:,:,:,oldtime,iblock) = p5*(UVEL(:,:,:,oldtime,iblock) + &
                                          UVEL(:,:,:,curtime,iblock))
         VVEL(:,:,:,oldtime,iblock) = p5*(VVEL(:,:,:,oldtime,iblock) + &
                                          VVEL(:,:,:,curtime,iblock))
         UVEL(:,:,:,curtime,iblock) = p5*(UVEL(:,:,:,curtime,iblock) + &
                                          UVEL(:,:,:,newtime,iblock))
         VVEL(:,:,:,curtime,iblock) = p5*(VVEL(:,:,:,curtime,iblock) + &
                                          VVEL(:,:,:,newtime,iblock))

         do n=1,POP_nt

            do k=2,POP_km
               TRACER(:,:,k,n,oldtime,iblock) = &
                          p5*(TRACER(:,:,k,n,oldtime,iblock) + &
                              TRACER(:,:,k,n,curtime,iblock))
               TRACER(:,:,k,n,curtime,iblock) = &
                          p5*(TRACER(:,:,k,n,curtime,iblock) + &
                              TRACER(:,:,k,n,newtime,iblock))
            end do
         end do

         if (sfc_layer_type == sfc_layer_varthick) then

            do n = 1,POP_nt

               TRACER(:,:,1,n,oldtime,iblock) = &
                   p5*((dz(1) + PSURF(:,:,oldtime,iblock)/grav)* &
                       TRACER(:,:,1,n,oldtime,iblock) + &
                       (dz(1) + PSURF(:,:,curtime,iblock)/grav)* &
                       TRACER(:,:,1,n,curtime,iblock) )
               TRACER(:,:,1,n,curtime,iblock) = &
                   p5*((dz(1) + PSURF(:,:,curtime,iblock)/grav)* &
                       TRACER(:,:,1,n,curtime,iblock) + &
                       (dz(1) + PSURF(:,:,newtime,iblock)/grav)* &
                       TRACER(:,:,1,n,newtime,iblock) )
            end do ! nt

            PSURF(:,:,oldtime,iblock) = p5*(PSURF(:,:,oldtime,iblock) + &
                                            PSURF(:,:,curtime,iblock))
            PSURF(:,:,curtime,iblock) = p5*(PSURF(:,:,curtime,iblock) + &
                                            PSURF(:,:,newtime,iblock))
            do n = 1,POP_nt

               TRACER(:,:,1,n,oldtime,iblock) = &
               TRACER(:,:,1,n,oldtime,iblock)/(dz(1) + &
                                      PSURF(:,:,oldtime,iblock)/grav)
               TRACER(:,:,1,n,curtime,iblock) = &
               TRACER(:,:,1,n,curtime,iblock)/(dz(1) + &
                                      PSURF(:,:,curtime,iblock)/grav)
            enddo

         else

            do n=1,POP_nt

               TRACER(:,:,1,n,oldtime,iblock) = &
                          p5*(TRACER(:,:,1,n,oldtime,iblock) + &
                              TRACER(:,:,1,n,curtime,iblock))
               TRACER(:,:,1,n,curtime,iblock) = &
                          p5*(TRACER(:,:,1,n,curtime,iblock) + &
                              TRACER(:,:,1,n,newtime,iblock))
            end do

            PSURF (:,:,oldtime,iblock) = &
                                  p5*(PSURF (:,:,oldtime,iblock) + &
                                      PSURF (:,:,curtime,iblock))
            PSURF (:,:,curtime,iblock) = &
                                  p5*(PSURF (:,:,curtime,iblock) + &
                                      PSURF (:,:,newtime,iblock))

         endif

         do k = 1,POP_km ! recalculate densities from averaged tracers
            call state(k,k,TRACER(:,:,k,1,oldtime,iblock), &
                           TRACER(:,:,k,2,oldtime,iblock), &
                           this_block, &
                         RHOOUT=RHO(:,:,k,oldtime,iblock))
            call state(k,k,TRACER(:,:,k,1,curtime,iblock), &
                           TRACER(:,:,k,2,curtime,iblock), &
                           this_block, &
                         RHOOUT=RHO(:,:,k,curtime,iblock))
         enddo

         !*** correct after avg
         PGUESS(:,:,iblock) = p5*(PGUESS(:,:,iblock) + &
                                   PSURF(:,:,newtime,iblock))
      end do ! block loop
      !$OMP END PARALLEL DO

   elseif (robert_filter_standard) then ! Standard Robert-Asselin filter

      !$OMP PARALLEL DO PRIVATE(iblock,this_block,k,n)

      do iblock = 1,nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)

         !*** filter 2-d fields

         UBTROP(:,:,curtime,iblock) = UBTROP(:,:,curtime,iblock) + &
                robert_curtime*(UBTROP(:,:,oldtime,iblock) + &
                                UBTROP(:,:,newtime,iblock) - &
                             c2*UBTROP(:,:,curtime,iblock))
         VBTROP(:,:,curtime,iblock) = VBTROP(:,:,curtime,iblock) + &
                robert_curtime*(VBTROP(:,:,oldtime,iblock) + &
                                VBTROP(:,:,newtime,iblock) - &
                             c2*VBTROP(:,:,curtime,iblock))
         GRADPX(:,:,curtime,iblock) = GRADPX(:,:,curtime,iblock) + &
                robert_curtime*(GRADPX(:,:,oldtime,iblock) + &
                                GRADPX(:,:,newtime,iblock) - &
                             c2*GRADPX(:,:,curtime,iblock))
         GRADPY(:,:,curtime,iblock) = GRADPY(:,:,curtime,iblock) + &
                robert_curtime*(GRADPY(:,:,oldtime,iblock) + &
                                GRADPY(:,:,newtime,iblock) - &
                             c2*GRADPY(:,:,curtime,iblock))

         !*** avg 3-d fields

         UVEL(:,:,:,curtime,iblock) = UVEL(:,:,:,curtime,iblock) + &
                robert_curtime*(UVEL(:,:,:,oldtime,iblock) + &
                                UVEL(:,:,:,newtime,iblock) - &
                             c2*UVEL(:,:,:,curtime,iblock))
         VVEL(:,:,:,curtime,iblock) = VVEL(:,:,:,curtime,iblock) + &
                robert_curtime*(VVEL(:,:,:,oldtime,iblock) + &
                                VVEL(:,:,:,newtime,iblock) - &
                             c2*VVEL(:,:,:,curtime,iblock))

         do n=1,POP_nt
            do k=2,POP_km
               TRACER(:,:,k,n,curtime,iblock) = TRACER(:,:,k,n,curtime,iblock) + &
                    robert_curtime*(TRACER(:,:,k,n,oldtime,iblock) + &
                                    TRACER(:,:,k,n,newtime,iblock) - &
                                 c2*TRACER(:,:,k,n,curtime,iblock))
            end do
         end do

         if (sfc_layer_type == sfc_layer_varthick) then

            do n = 1,POP_nt

               TRACER(:,:,1,n,curtime,iblock) = &
                   TRACER(:,:,1,n,curtime,iblock) * &
                      (dz(1) + PSURF(:,:,curtime,iblock)/grav) + &
                   robert_curtime*((dz(1) + PSURF(:,:,oldtime,iblock)/grav)* &
                                   TRACER(:,:,1,n,oldtime,iblock) + &
                                   (dz(1) + PSURF(:,:,newtime,iblock)/grav)* &
                                   TRACER(:,:,1,n,newtime,iblock) - &
                                c2*(dz(1) + PSURF(:,:,curtime,iblock)/grav)* &
                                   TRACER(:,:,1,n,curtime,iblock) )

            end do ! nt

            PSURF(:,:,curtime,iblock) = PSURF(:,:,curtime,iblock) + &
                robert_curtime*(PSURF(:,:,oldtime,iblock) + &
                                PSURF(:,:,newtime,iblock) - &
                             c2*PSURF(:,:,curtime,iblock))
            do n = 1,POP_nt

               TRACER(:,:,1,n,curtime,iblock) = &
               TRACER(:,:,1,n,curtime,iblock)/(dz(1) + &
                                      PSURF(:,:,curtime,iblock)/grav)
            enddo

         else

            do n=1,POP_nt

               TRACER(:,:,1,n,curtime,iblock) = TRACER(:,:,1,n,curtime,iblock) + &
                    robert_curtime*(TRACER(:,:,1,n,oldtime,iblock) + &
                                    TRACER(:,:,1,n,newtime,iblock) - &
                                 c2*TRACER(:,:,1,n,curtime,iblock))
            end do

            PSURF (:,:,curtime,iblock) = PSURF (:,:,curtime,iblock) + &
                    robert_curtime*(PSURF (:,:,oldtime,iblock) + &
                                    PSURF (:,:,newtime,iblock) - &
                                 c2*PSURF (:,:,curtime,iblock))

         endif

         do k = 1,POP_km ! recalculate densities from averaged tracers
            call state(k,k,TRACER(:,:,k,1,curtime,iblock), &
                           TRACER(:,:,k,2,curtime,iblock), &
                           this_block, &
                           RHOOUT=RHO(:,:,k,curtime,iblock))
         enddo

! just like standard leapfrog
         FW_OLD(:,:,iblock) = FW(:,:,iblock)

         !*** correct after avg

! extrapolate again based on new PSURF(curtime)
         PGUESS(:,:,iblock) = c3*(PSURF(:,:,newtime,iblock) - &
                                  PSURF(:,:,curtime,iblock)) + &
                                  PSURF(:,:,oldtime,iblock)
      end do ! block loop

! now reorder just like standard leapfrog
      tmptime = oldtime
      oldtime = curtime
      curtime = newtime
      newtime = tmptime

   elseif (robert_filter_modified) then ! Modified Robert-Asselin filter

      allocate( FILTER_TEMP(POP_nxBlock,POP_nyBlock), stat=istat)

      if (istat > 0) then
         call POP_ErrorSet(errorCode, &
            'step: error allocating temporary array FILTER_TEMP')
         return
      endif

      !$OMP PARALLEL DO PRIVATE(iblock,this_block,k,n)

      do iblock = 1,nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)

         !*** filter 2-d fields

         FILTER_TEMP = UBTROP(:,:,oldtime,iblock) + &
                       UBTROP(:,:,newtime,iblock) - &
                    c2*UBTROP(:,:,curtime,iblock)
         UBTROP(:,:,curtime,iblock) = UBTROP(:,:,curtime,iblock) + &
                                      robert_curtime*FILTER_TEMP
         UBTROP(:,:,newtime,iblock) = UBTROP(:,:,newtime,iblock) + &
                                      robert_newtime*FILTER_TEMP

         FILTER_TEMP = VBTROP(:,:,oldtime,iblock) + &
                       VBTROP(:,:,newtime,iblock) - &
                    c2*VBTROP(:,:,curtime,iblock)
         VBTROP(:,:,curtime,iblock) = VBTROP(:,:,curtime,iblock) + &
                                      robert_curtime*FILTER_TEMP
         VBTROP(:,:,newtime,iblock) = VBTROP(:,:,newtime,iblock) + &
                                      robert_newtime*FILTER_TEMP

         FILTER_TEMP = GRADPX(:,:,oldtime,iblock) + &
                       GRADPX(:,:,newtime,iblock) - &
                    c2*GRADPX(:,:,curtime,iblock)
         GRADPX(:,:,curtime,iblock) = GRADPX(:,:,curtime,iblock) + &
                                      robert_curtime*FILTER_TEMP
         GRADPX(:,:,newtime,iblock) = GRADPX(:,:,newtime,iblock) + &
                                      robert_newtime*FILTER_TEMP

         FILTER_TEMP = GRADPY(:,:,oldtime,iblock) + &
                       GRADPY(:,:,newtime,iblock) - &
                    c2*GRADPY(:,:,curtime,iblock)
         GRADPY(:,:,curtime,iblock) = GRADPY(:,:,curtime,iblock) + &
                                      robert_curtime*FILTER_TEMP
         GRADPY(:,:,newtime,iblock) = GRADPY(:,:,newtime,iblock) + &
                                      robert_newtime*FILTER_TEMP

         !*** avg 3-d fields

         do k=1,POP_km
            FILTER_TEMP = UVEL(:,:,k,oldtime,iblock) + &
                          UVEL(:,:,k,newtime,iblock) - &
                       c2*UVEL(:,:,k,curtime,iblock)
            UVEL(:,:,k,curtime,iblock) = UVEL(:,:,k,curtime,iblock) + &
                                         robert_curtime*FILTER_TEMP
            UVEL(:,:,k,newtime,iblock) = UVEL(:,:,k,newtime,iblock) + &
                                         robert_newtime*FILTER_TEMP

            FILTER_TEMP = VVEL(:,:,k,oldtime,iblock) + &
                          VVEL(:,:,k,newtime,iblock) - &
                       c2*VVEL(:,:,k,curtime,iblock)
            VVEL(:,:,k,curtime,iblock) = VVEL(:,:,k,curtime,iblock) + &
                                         robert_curtime*FILTER_TEMP
            VVEL(:,:,k,newtime,iblock) = VVEL(:,:,k,newtime,iblock) + &
                                         robert_newtime*FILTER_TEMP
         end do

         do n=1,POP_nt

            do k=2,POP_km
               FILTER_TEMP = TRACER(:,:,k,n,oldtime,iblock) + &
                             TRACER(:,:,k,n,newtime,iblock) - &
                          c2*TRACER(:,:,k,n,curtime,iblock)
               TRACER(:,:,k,n,curtime,iblock) = TRACER(:,:,k,n,curtime,iblock) &
                                              + robert_curtime*FILTER_TEMP
               TRACER(:,:,k,n,newtime,iblock) = TRACER(:,:,k,n,newtime,iblock) &
                                              + robert_newtime*FILTER_TEMP

            end do
         end do

         if (sfc_layer_type == sfc_layer_varthick) then

            do n = 1,POP_nt

               FILTER_TEMP = (dz(1) + PSURF(:,:,oldtime,iblock)/grav)* &
                             TRACER(:,:,1,n,oldtime,iblock) + &
                             (dz(1) + PSURF(:,:,newtime,iblock)/grav)* &
                             TRACER(:,:,1,n,newtime,iblock) - &
                          c2*(dz(1) + PSURF(:,:,curtime,iblock)/grav)* &
                             TRACER(:,:,1,n,curtime,iblock)
               TRACER(:,:,1,n,curtime,iblock) = &
                   TRACER(:,:,1,n,curtime,iblock) * &
                      (dz(1) + PSURF(:,:,curtime,iblock)/grav) + &
                       robert_curtime*FILTER_TEMP
               TRACER(:,:,1,n,newtime,iblock) = &
                   TRACER(:,:,1,n,newtime,iblock) * &
                      (dz(1) + PSURF(:,:,newtime,iblock)/grav) + &
                       robert_newtime*FILTER_TEMP
            end do ! nt

            FILTER_TEMP = PSURF(:,:,oldtime,iblock) + &
                          PSURF(:,:,newtime,iblock) - &
                       c2*PSURF(:,:,curtime,iblock)
            PSURF(:,:,curtime,iblock) = PSURF(:,:,curtime,iblock) + &
                                        robert_curtime*FILTER_TEMP
            PSURF(:,:,newtime,iblock) = PSURF(:,:,newtime,iblock) + &
                                        robert_newtime*FILTER_TEMP

            do n = 1,POP_nt

               TRACER(:,:,1,n,curtime,iblock) = &
               TRACER(:,:,1,n,curtime,iblock)/(dz(1) + &
                                      PSURF(:,:,curtime,iblock)/grav)
               TRACER(:,:,1,n,newtime,iblock) = &
               TRACER(:,:,1,n,newtime,iblock)/(dz(1) + &
                                      PSURF(:,:,newtime,iblock)/grav)
            enddo

         else

            do n=1,POP_nt

               FILTER_TEMP = TRACER(:,:,1,n,oldtime,iblock) + &
                              TRACER(:,:,1,n,newtime,iblock) - &
                           c2*TRACER(:,:,1,n,curtime,iblock)
               TRACER(:,:,1,n,curtime,iblock) = TRACER(:,:,1,n,curtime,iblock) &
                                              + robert_curtime*FILTER_TEMP
               TRACER(:,:,1,n,newtime,iblock) = TRACER(:,:,1,n,newtime,iblock) &
                                              + robert_newtime*FILTER_TEMP
            end do

            FILTER_TEMP = PSURF(:,:,oldtime,iblock) + &
                          PSURF(:,:,newtime,iblock) - &
                       c2*PSURF(:,:,curtime,iblock)
            PSURF(:,:,curtime,iblock) = PSURF(:,:,curtime,iblock) + &
                                        robert_curtime*FILTER_TEMP
            PSURF(:,:,newtime,iblock) = PSURF(:,:,newtime,iblock) + &
                                        robert_newtime*FILTER_TEMP

         endif

         do k = 1,POP_km ! recalculate densities from averaged tracers
            call state(k,k,TRACER(:,:,k,1,curtime,iblock), &
                           TRACER(:,:,k,2,curtime,iblock), &
                           this_block, &
                           RHOOUT=RHO(:,:,k,curtime,iblock))
            call state(k,k,TRACER(:,:,k,1,newtime,iblock), &
                           TRACER(:,:,k,2,newtime,iblock), &
                           this_block, &
                           RHOOUT=RHO(:,:,k,newtime,iblock))
         enddo

! just like standard leapfrog
         FW_OLD(:,:,iblock) = FW(:,:,iblock)

         !*** correct after avg

! extrapolate again based on new PSURF(curtime)
         PGUESS(:,:,iblock) = c3*(PSURF(:,:,newtime,iblock) - &
                                  PSURF(:,:,curtime,iblock) + &
                                  PSURF(:,:,oldtime,iblock))
      end do ! block loop

! now reorder just like standard leapfrog
      tmptime = oldtime
      oldtime = curtime
      curtime = newtime
      newtime = tmptime

      deallocate( FILTER_TEMP )

   else ! non-averaging step

      !$OMP PARALLEL DO PRIVATE(iblock)
      do iblock = 1,nblocks_clinic

         if (mix_pass == 2) then ! reset time n variables on 2nd pass matsuno

            UBTROP(:,:,curtime,iblock) = UBTROP(:,:,oldtime,iblock)
            VBTROP(:,:,curtime,iblock) = VBTROP(:,:,oldtime,iblock)
            UVEL(:,:,:,curtime,iblock) = UVEL(:,:,:,oldtime,iblock)
            VVEL(:,:,:,curtime,iblock) = VVEL(:,:,:,oldtime,iblock)
            TRACER(:,:,:,:,curtime,iblock) = &
                                     TRACER(:,:,:,:,oldtime,iblock)
            RHO(:,:,:,curtime,iblock) = RHO(:,:,:,oldtime,iblock)

         endif

         FW_OLD(:,:,iblock) = FW(:,:,iblock)

      end do ! block loop
      !$OMP END PARALLEL DO


      tmptime = oldtime
      oldtime = curtime
      curtime = newtime
      newtime = tmptime

   endif







!-----------------------------------------------------------------------
!
! end of timestep, all variables updated
! compute and print some more diagnostics
!
!-----------------------------------------------------------------------

!maltrud deal with qflux later
! if (registry_match('lcoupled')) then
! if ( liceform .and. check_time_flag(ice_cpl_flag) ) then
! call tavg_increment_sum_qflux(const=tlast_ice)
! !$OMP PARALLEL DO PRIVATE(iblock)
! do iblock = 1,nblocks_clinic
! call ice_flx_to_coupler(TRACER(:,:,:,:,curtime,iblock),iblock)
! if (tavg_requested(tavg_id('QFLUX')) ) &
! call accumulate_tavg_field(QFLUX(:,:,iblock), tavg_id('QFLUX'), &
! iblock,1,const=tlast_ice)
!
! end do ! block loop
! !$OMP END PARALLEL DO
!-----------------------------------------------------------------------
! time-averaging for ice formation related quantities
!-----------------------------------------------------------------------
! if (nt > 2) call passive_tracers_tavg_FvICE(cp_over_lhfusion, QICE)
! endif
! endif

   call diag_global_afterupdate
   call diag_print
   call diag_transport

!maltrud merge ncar does moorings at the end of the step
! if ( eod .and. ldiag_velocity) then
! call diag_velocity
! endif

!-----------------------------------------------------------------------
!
! stop step timer
!
!-----------------------------------------------------------------------

  call timer_stop(timer_step)

!-----------------------------------------------------------------------
!EOC






   end subroutine step

!***********************************************************************

!BOP
! !IROUTINE: init_step
! !INTERFACE:

 subroutine init_step

! !DESCRIPTION:
! This routine initializes timers and flags used in subroutine step.
!
! !REVISION HISTORY:
! added 17 August 2007 njn01

!EOP
!BOC

!-----------------------------------------------------------------------
!
! get tavg time-flag handle and initialize timers
!
!-----------------------------------------------------------------------

   tavg_flag = get_time_flag_id('tavg')

   call get_timer(timer_step,'STEP',1,distrb_clinic%nprocs)
   call get_timer(timer_baroclinic,'BAROCLINIC',1,distrb_clinic%nprocs)
   call get_timer(timer_barotropic,'BAROTROPIC',1,distrb_clinic%nprocs)
   call get_timer(timer_3dupdate,'3D-UPDATE',1,distrb_clinic%nprocs)


!-----------------------------------------------------------------------
!EOC

   end subroutine init_step
 end module step_mod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
