!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 module baroclinic

!BOP
! !MODULE: baroclinic
!
! !DESCRIPTION:
!  Contains main driver routines and variables for computing the 
!  baroclinic velocities and tracer fields.
!
! !REVISION HISTORY:
!  SVN:$Id: baroclinic.F90 14725 2009-03-04 22:50:06Z njn01 $

! !USES:

   use POP_KindsMod
   use POP_ErrorMod
   use POP_CommMod
   use POP_FieldMod 
   use POP_GridHorzMod
   use POP_HaloMod 

   use POP_DomainSizeMod, only : POP_km, POP_nt, POP_maxBlocksClinic

   use blocks, only: nx_block, ny_block, block, get_block
!   use distribution, only: 
   use domain, only: nblocks_clinic, blocks_clinic, POP_haloClinic
   use constants, only: delim_fmt, blank_fmt, p5, field_loc_center,          &
       field_type_scalar, c0, c1, c2, grav, ndelim_fmt,                      &
       hflux_factor, salinity_factor, salt_to_ppt
   use prognostic, only: TRACER, UVEL, VVEL, mixtime, &
       RHO, newtime, oldtime, curtime, PSURF
   use broadcast, only: broadcast_scalar
   use communicate, only: my_task, master_task
   use grid, only: FCOR, DZU, HUR, KMU, KMT, sfc_layer_type,                 &
       sfc_layer_varthick, partial_bottom_cells, dz, DZT, CALCT, dzw, dzr,   &
       FCORT, TAREA_R
   use advection, only: advu, advt, comp_flux_vel_ghost
   use pressure_grad, only: lpressure_avg, gradp
   use horizontal_mix, only: hdiffu, hdifft
   use vertical_mix, only: vmix_coeffs, implicit_vertical_mix, vdiffu,       &
       vdifft, impvmixt, impvmixu, impvmixt_correct, convad
   use vmix_kpp, only: add_kpp_sources
   use diagnostics, only: ldiag_cfl, cfl_check, ldiag_global,                &
       DIAG_KE_ADV_2D, DIAG_KE_PRESS_2D, DIAG_KE_HMIX_2D, DIAG_KE_VMIX_2D,   &
       DIAG_TRACER_HDIFF_2D, DIAG_PE_2D, DIAG_TRACER_ADV_2D,                 &
       DIAG_TRACER_SFC_FLX, DIAG_TRACER_VDIFF_2D, DIAG_TRACER_SOURCE_2D
   use movie, only: define_movie_field, movie_requested, update_movie_field
   use state_mod
   use ice, only: liceform, ice_formation, increment_tlast_ice
   use time_management, only: mix_pass, leapfrogts, impcor, c2dtu, beta,     &
       gamma, c2dtt
   use io_types, only: nml_in, nml_filename, stdout
   use tavg, only: define_tavg_field, tavg_requested, accumulate_tavg_field, &
       tavg_method_max, tavg_method_min
   use forcing_fields, only: STF, SMF, lsmft_avail, SMFT, TFW
   use forcing_shf, only: SHF_QSW
   use forcing_sfwf, only: lfw_as_salt_flx
   use sw_absorption, only:  add_sw_absorb
   use forcing_pt_interior, only: set_pt_interior
   use forcing_s_interior, only: set_s_interior
   use passive_tracers, only: set_interior_passive_tracers,  &
       reset_passive_tracers, tavg_passive_tracers,          &
       movie_passive_tracers
   use floats, only: define_float_field, float_requested, update_float_buffer
   use operators, only: zcurl
   use exit_mod, only: sigAbort, exit_pop

   use gpu_mod

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: init_baroclinic,          &
             baroclinic_driver,        &
             baroclinic_correct_adjust

! !PUBLIC DATA MEMBERS:

   logical (POP_logical) :: &
      reset_to_freezing   ! flag to prevent very cold water

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  ids for tavg diagnostics computed from baroclinic
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      tavg_UDP,          &! tavg id for pressure grad work
      tavg_TEMP,         &! tavg id for temperature
      tavg_TEMP_MAX,     &! tavg id for maximum temperature
      tavg_TEMP_MIN,     &! tavg id for maximum temperature
      tavg_dTEMP_POS_3D, &! tavg id for positive temperature timestep difference
      tavg_dTEMP_POS_2D, &! tavg id for positive temperature timestep difference
      tavg_dTEMP_NEG_3D, &! tavg id for negative temperature timestep difference
      tavg_dTEMP_NEG_2D, &! tavg id for negative temperature timestep difference
      tavg_SALT,         &! tavg id for salinity
      tavg_TEMP2,        &! tavg id for temperature squared
      tavg_SALT2,        &! tavg id for salinity    squared
      tavg_UVEL,         &! tavg id for U velocity
      tavg_VVEL,         &! tavg id for V velocity
      tavg_KE,           &! tavg id for kinetic energy
      tavg_ST,           &! tavg id for salt*temperature
      tavg_RHO,          &! tavg id for in-situ density
      tavg_UV,           &! tavg id for u times v
      tavg_RESID_T,      &! free-surface residual flux (T)
      tavg_RESID_S        ! free-surface residual flux (S)

!-----------------------------------------------------------------------
!
!  ids for movie diagnostics computed from baroclinic
!
!-----------------------------------------------------------------------

   integer (POP_i4), dimension(0:POP_km) :: &
      movie_TEMP,         &! movie id for temperature
      movie_SALT,         &! movie id for salinity
      movie_UVEL,         &! movie id for U velocity
      movie_VVEL,         &! movie id for V velocity
      movie_RHO            ! movie id for in-situ density

!-----------------------------------------------------------------------
!  
!  ids for float diagnostics computed from baroclinic
!  
!-----------------------------------------------------------------------
   
   integer (POP_i4) :: & 
      float_UVEL,       &! tavg id for u
      float_VVEL,       &! tavg id for v
      float_TEMP,       &! tavg id for temperature
      float_SALT,       &! float id for salinity
      float_PV,         &! float id for potential vorticity
      float_PD           ! float id for potential density ref to surface

!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: init_baroclinic
! !INTERFACE:

 subroutine init_baroclinic(errorCode)

! !DESCRIPTION:
!  Initializes some baroclinic options.
!
! !REVISION HISTORY:
!  same as module

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode             ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      nml_error          ! namelist i/o error flag
   integer (POP_i4) :: k

   namelist /baroclinic_nml/ reset_to_freezing

!-----------------------------------------------------------------------
!
!  read options from namelist and broadcast
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   reset_to_freezing = .true.

   if (my_task == master_task) then
      open (nml_in, file=nml_filename, status='old',iostat=nml_error)
      if (nml_error /= 0) then
         nml_error = -1
      else
         nml_error =  1
      endif
      do while (nml_error > 0)
         read(nml_in, nml=baroclinic_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nml_in)
   endif

   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
      call exit_POP(sigAbort,'ERROR reading baroclinic_nml')
   endif

   if (my_task == master_task) then
      write(stdout,blank_fmt)
      write(stdout,ndelim_fmt)
      write(stdout,blank_fmt)
      write(stdout,*) ' Baroclinic:'
      write(stdout,blank_fmt)
      write(stdout,*) ' baroclinic_nml namelist settings:'
      write(stdout,blank_fmt)
      write(stdout,baroclinic_nml)
      write(stdout,blank_fmt)
      write(stdout,delim_fmt)

      write(stdout,'(a18)') 'Baroclinic options'
      write(stdout,blank_fmt)
      write(stdout,delim_fmt)
      if (reset_to_freezing) then
         write(stdout,'(a40)') &
                          'Surface temperature reset to freezing on'
      else
         write(stdout,'(a41)') &
                          'Surface temperature reset to freezing off'
      endif
   endif

   call broadcast_scalar(reset_to_freezing, master_task)

!-----------------------------------------------------------------------
!
!  define tavg fields computed from baroclinic driver routines
!
!-----------------------------------------------------------------------

   call define_tavg_field(tavg_UDP,'UDP',3,                            &
                          long_name='Pressure work',                   &
                          units='erg', grid_loc='3221')

   call define_tavg_field(tavg_UVEL,'UVEL',3,                          &
                          long_name='Velocity in grid-x direction',    &
                          units='centimeter/s', grid_loc='3221')

   call define_tavg_field(tavg_VVEL,'VVEL',3,                          &
                          long_name='Velocity in grid-y direction',    &
                          units='centimeter/s', grid_loc='3221')

   call define_tavg_field(tavg_KE,'KE',3,                              &
                          long_name='Horizontal Kinetic Energy',       &
                          units='centimeter^2/s^2', grid_loc='3221')

   call define_tavg_field(tavg_TEMP,'TEMP',3,                          &
                          long_name='Potential Temperature',           &
                          units='degC', grid_loc='3111')

   call define_tavg_field(tavg_TEMP_MAX,'TEMP_MAX',3,                  &
                          tavg_method=tavg_method_max,                 &
                          long_name='Maximum Potential Temperature',   &
                          units='degC', grid_loc='3111')

   call define_tavg_field(tavg_TEMP_MIN,'TEMP_MIN',3,                  &
                          tavg_method=tavg_method_min,                 &
                          long_name='Minimum Potential Temperature',   &
                          units='degC', grid_loc='3111')

   call define_tavg_field(tavg_dTEMP_POS_3D,'dTEMP_POS_3D',3,          &
                          tavg_method=tavg_method_max,                 &
                          long_name='max pos temperature timestep diff', &
                          units='degC', grid_loc='3111')

   call define_tavg_field(tavg_dTEMP_POS_2D,'dTEMP_POS_2D',2,          &
                          tavg_method=tavg_method_max,                 &
                          long_name='max pos column temperature timestep diff', &
                          units='degC', grid_loc='2110')

   call define_tavg_field(tavg_dTEMP_NEG_3D,'dTEMP_NEG_3D',3,          &
                          tavg_method=tavg_method_min,                 &
                          long_name='min neg temperature timestep diff', &
                          units='degC', grid_loc='3111')

   call define_tavg_field(tavg_dTEMP_NEG_2D,'dTEMP_NEG_2D',2,          &
                          tavg_method=tavg_method_min,                 &
                          long_name='min neg column temperature timestep diff', &
                          units='degC', grid_loc='2110')

   call define_tavg_field(tavg_SALT,'SALT',3,                          &
                          long_name='Salinity',                        &
                          units='gram/kilogram', grid_loc='3111')

   call define_tavg_field(tavg_TEMP2,'TEMP2',3,                        &
                          long_name='Temperature**2',                  &
                          units='degC^2', grid_loc='3111')

   call define_tavg_field(tavg_SALT2,'SALT2',3,                        &
                          long_name='Salinity**2 ',                    &
                          units='(gram/gram)^2', grid_loc='3111')

   call define_tavg_field(tavg_ST,'ST',3,                              &
                          long_name='Temperature*Salinity',            &
                          units='degC*gram/gram', grid_loc='3111')

   call define_tavg_field(tavg_RHO,'RHO',3,                            &
                          long_name='In-Situ Density',                 &
                          units='gram/centimeter^3', grid_loc='3111')

   call define_tavg_field(tavg_UV,'UV',3,                              &
                          long_name='UV velocity product',             &
                          units='centimeter^2/s^2', grid_loc='3221')

   call define_tavg_field(tavg_RESID_T,'RESID_T',2,                           &
                    long_name='Free-Surface Residual Flux (T)',               &
                          units='watt/m^2', grid_loc='2110')

   call define_tavg_field(tavg_RESID_S,'RESID_S',2,                           &
                    long_name='Free-Surface Residual Flux (S)',               &
                          units='kg/m^2/s', grid_loc='2110')

!-----------------------------------------------------------------------
!
!  define movie fields computed from baroclinic driver routines
!
!-----------------------------------------------------------------------

   do k = 1, POP_km

      call define_movie_field(movie_UVEL(k),'UVEL',k,                  &
                          long_name='Zonal Velocity',                  &
                          units='cm/s', grid_loc='3221')

      call define_movie_field(movie_VVEL(k),'VVEL',k,                  &
                          long_name='Meridional Velocity',             &
                          units='cm/s', grid_loc='3221')

      call define_movie_field(movie_TEMP(k),'TEMP',k,                  &
                          long_name='Potential Temperature',           &
                          units='degC', grid_loc='3111')

      call define_movie_field(movie_SALT(k),'SALT',k,                  &
                          long_name='Salinity',                        &
                          units='psu', grid_loc='3111')

      call define_movie_field(movie_RHO(k),'RHO',k,                    &
                          long_name='In-situ density',                 &
                          units='sigma units', grid_loc='3111')

   enddo

!  now define depth integrated fields located in k = 0 id slot

  k = 0
  call define_movie_field(movie_UVEL(k),'UVEL_zint',k,                    &
                           long_name='Depth-Integrated Zonal Velocity',   &
                           units='cm^2/s', grid_loc='2221')

  call define_movie_field(movie_VVEL(k),'VVEL_zint',k,                    &
                          long_name='Depth-Integrated Meridional Velocity',  &
                          units='cm^2/s', grid_loc='2221')

!-----------------------------------------------------------------------
!
!  define float fields computed from baroclinic driver routines
!
!-----------------------------------------------------------------------

      call define_float_field(float_TEMP,'TEMP',                    &
                          long_name='Potential Temperature',           &
                          units='degC', grid_loc='3111')

      call define_float_field(float_SALT,'SALT',                    &
                          long_name='Salinity',                        &
                          units='psu', grid_loc='3111')

      call define_float_field(float_UVEL,'UVEL',                    &
                          long_name='Zonal Velocity',                  &
                          units='cm/s', grid_loc='3221')

      call define_float_field(float_VVEL,'VVEL',                    &
                          long_name='Meridional Velocity',             &
                          units='cm/s', grid_loc='3221')

      call define_float_field(float_PV,'PV',                              &
                          long_name='Potential Vorticity',             &
                          units='1/s', grid_loc='3111')

      call define_float_field(float_PD,'PD',                              &
                          long_name='Potential Density Ref to Surface',&
                          units='gram/centimeter^3', grid_loc='3111')

!-----------------------------------------------------------------------
!EOC

 end subroutine init_baroclinic

!***********************************************************************
!BOP
! !IROUTINE: baroclinic_driver
! !INTERFACE:

 subroutine baroclinic_driver(ZX,ZY,DH,DHU, errorCode)

! !DESCRIPTION:
!  This routine is the main driver for the explicit time integration of
!  baroclinic velocities $(u',v')$ and tracer fields $T$. 
!
!  Tracer equations:
!  \begin{equation}
!     (T^{n+1}-T^{n-1})/(2 \Delta t) = -L(T^n) + D_H(T^{n-1}) + 
!                                      D_V(T^{n-1}) + S
!  \end{equation}
!  where $S$ are source terms, $L$ is the advection operator and
!  $D$ is the diffusion operator in the horizontal and vertical
!  directions.
!
!  Momentum equations:
!  \begin{eqnarray}
!     (u'-u^{n-1})/(2 \Delta t) - f*\alpha*(v'-v^{n-1}) &=& F_x \\
!     (v'-v^{n-1})/(2 \Delta t) + f*\alpha*(u'-u^{n-1}) &=& F_y
!  \end{eqnarray}
!  \begin{eqnarray}
!     \tilde{u}' &=& u' - {1\over{H_U}}\sum_{k=1}^{km}dz_k u_k' \\
!     \tilde{v}' &=& v' - {1\over{H_U}}\sum_{k=1}^{km}dz_k v_k'
!  \end{eqnarray}
!
!  This routine calculates baroclinic velocities and tracers at
!  the new time level and stores them in arrays UVEL,VVEL,and TRACER
!  with time index newtime.  The arrays UVEL,VVEL, and TRACER
!  with time indices curtime and oldtime are not updated for the next
!  timestep until near the end of routine step.
!
!  The above equations are written for the case of (leapfrog)
!  implicit treatment of the coriolis terms.  if these terms are
!  treated explicitly then the coriolis terms appear only in the 
!  forcing terms $(F_x,F_y)$, which are calculated in clinic. 
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   real (POP_r8), dimension(nx_block,ny_block,POP_maxBlocksClinic), &
      intent(in) :: &
      DH, DHU              ! change in surface height at T,U points

! !OUTPUT PARAMETERS:

   real (POP_r8), dimension(nx_block,ny_block,POP_maxBlocksClinic), &
      intent(out) :: &
      ZX, ZY               ! vertical integrals of forcing

   integer (POP_i4), intent(out) :: &
      errorCode            ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) ::  &
      i,j,                &! dummy indices for horizontal directions
      n,k,                &! dummy indices for vertical level, tracer
      iblock,             &! counter for block loops
      kp1,km1              ! level index for k+1, k-1 levels

   real (POP_r8), dimension(nx_block,ny_block) :: & 
      FX,FY,              &! sum of r.h.s. forcing terms
      WORK1,WORK2,        &! local work space
      WUK,                &! vertical velocity at top of U box
      WTK                  ! vertical velocity at top of T box

   real (POP_r8)    ::        &
      factor
   type (block) ::        &
      this_block           ! block information for current block

   real (POP_r8), dimension(:,:,:), allocatable :: &
      WORK3D1,WORK3D2              ! temporaries for float PV and PD

!-----------------------------------------------------------------------
!
!  compute flux velocities in ghost cells
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   call comp_flux_vel_ghost(DH, errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'baroclinic_driver: error in comp_flux_vel_ghost')
      return
   endif

!-----------------------------------------------------------------------
!
!  first block loop to update tracers
!
!-----------------------------------------------------------------------


!jw    !$OMP PARALLEL DO PRIVATE(iblock,k,kp1,km1,this_block,WTK,WORK1,factor)

   do iblock = 1,nblocks_clinic
      this_block = get_block(blocks_clinic(iblock),iblock)  

      do k = 1,POP_km 

         kp1 = k+1
         km1 = k-1
         if (k == 1) km1 = 1
         if (k == POP_km) kp1 = POP_km

!-----------------------------------------------------------------------
!
!        compute vertical viscosity and diffusion coeffs
!
!-----------------------------------------------------------------------

!maltrud merge lanl also passes curtime

         if (lsmft_avail) then
            call vmix_coeffs(k,TRACER (:,:,:,:,mixtime,iblock), &
                               UVEL   (:,:,:  ,mixtime,iblock), &
                               VVEL   (:,:,:  ,mixtime,iblock), &
                               RHO    (:,:,:  ,mixtime,iblock), &
                               STF    (:,:,:          ,iblock), &
                               SHF_QSW(:,:            ,iblock), &
                               this_block, SMFT=SMFT(:,:,:,iblock))
         else
            call vmix_coeffs(k,TRACER (:,:,:,:,mixtime,iblock), &
                               UVEL   (:,:,:  ,mixtime,iblock), &
                               VVEL   (:,:,:  ,mixtime,iblock), &
                               RHO    (:,:,:  ,mixtime,iblock), &
                               STF    (:,:,:          ,iblock), &
                               SHF_QSW(:,:            ,iblock), &
                               this_block, SMF=SMF(:,:,:,iblock))
         endif

!-----------------------------------------------------------------------
!
!        calculate level k tracers at new time
!
!-----------------------------------------------------------------------


         call tracer_update(k, WTK,                             &
                               TRACER (:,:,:,:,newtime,iblock), &
                               TRACER (:,:,:,:,oldtime,iblock), &
                               TRACER (:,:,:,:,mixtime,iblock), &
                               TRACER (:,:,:,:,curtime,iblock), &
                               UVEL   (:,:,:  ,curtime,iblock), &
                               VVEL   (:,:,:  ,curtime,iblock), &
                               UVEL   (:,:,:  ,mixtime,iblock), &
                               VVEL   (:,:,:  ,mixtime,iblock), &
                               RHO    (:,:,:  ,curtime,iblock), &
                               STF    (:,:,:          ,iblock), &
                               TFW    (:,:,:          ,iblock), &
                               SHF_QSW(:,:            ,iblock), &
                               DH     (:,:            ,iblock), &
                               PSURF  (:,:    ,oldtime,iblock), &
                               PSURF  (:,:    ,curtime,iblock), &
                               this_block, errorCode)

         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'baroclinic_driver: error in tracer_update')
            return
         endif

!-----------------------------------------------------------------------
!
!        accumulate some tavg diagnostics if requested
!
!-----------------------------------------------------------------------

         if (mix_pass /= 1) then

         if (tavg_requested(tavg_UVEL)) then
            call accumulate_tavg_field(UVEL(:,:,k,curtime,iblock), &
                                       tavg_UVEL,iblock,k)
         endif

         if (tavg_requested(tavg_VVEL)) then
            call accumulate_tavg_field(VVEL(:,:,k,curtime,iblock), &
                                       tavg_VVEL,iblock,k)
         endif

         if (tavg_requested(tavg_KE)) then
            call accumulate_tavg_field(p5*(UVEL(:,:,k,curtime,iblock)**2 + &
                                           VVEL(:,:,k,curtime,iblock)**2), &
                                       tavg_KE,iblock,k)
         endif

         if (tavg_requested(tavg_UV)) then
            call accumulate_tavg_field(UVEL(:,:,k,curtime,iblock)* &
                                       VVEL(:,:,k,curtime,iblock), &
                                       tavg_UV,iblock,k)
         endif

         if (tavg_requested(tavg_TEMP)) then
            call accumulate_tavg_field(TRACER(:,:,k,1,curtime,iblock), &
                                       tavg_TEMP,iblock,k)
         endif

         if (tavg_requested(tavg_TEMP_MAX)) then
            call accumulate_tavg_field(TRACER(:,:,k,1,curtime,iblock), &
                                       tavg_TEMP_MAX,iblock,k)
         endif

         if (tavg_requested(tavg_TEMP_MIN)) then
            call accumulate_tavg_field(TRACER(:,:,k,1,curtime,iblock), &
                                       tavg_TEMP_MIN,iblock,k)
         endif

         if (tavg_requested(tavg_dTEMP_POS_3D)) then
            call accumulate_tavg_field(max(TRACER(:,:,k,1,curtime,iblock) - &
                                           TRACER(:,:,k,1,oldtime,iblock), c0), &
                                       tavg_dTEMP_POS_3D,iblock,k)
         endif

         if (tavg_requested(tavg_dTEMP_POS_2D)) then
            call accumulate_tavg_field(max(TRACER(:,:,k,1,curtime,iblock) - &
                                           TRACER(:,:,k,1,oldtime,iblock), c0), &
                                       tavg_dTEMP_POS_2D,iblock,1)
         endif

         if (tavg_requested(tavg_dTEMP_NEG_3D)) then
            call accumulate_tavg_field(min(TRACER(:,:,k,1,curtime,iblock) - &
                                           TRACER(:,:,k,1,oldtime,iblock), c0), &
                                       tavg_dTEMP_NEG_3D,iblock,k)
         endif

         if (tavg_requested(tavg_dTEMP_NEG_2D)) then
            call accumulate_tavg_field(min(TRACER(:,:,k,1,curtime,iblock) - &
                                           TRACER(:,:,k,1,oldtime,iblock), c0), &
                                       tavg_dTEMP_NEG_2D,iblock,1)
         endif

         if (tavg_requested(tavg_TEMP2)) then
            call accumulate_tavg_field(TRACER(:,:,k,1,curtime,iblock)**2, &
                                       tavg_TEMP2,iblock,k)
         endif

         if (tavg_requested(tavg_SALT)) then
            call accumulate_tavg_field(TRACER(:,:,k,2,curtime,iblock), &
                                       tavg_SALT,iblock,k)
         endif

         if (tavg_requested(tavg_SALT2)) then
            call accumulate_tavg_field(TRACER(:,:,k,2,curtime,iblock)**2, &
                                       tavg_SALT2,iblock,k)
         endif

         if (tavg_requested(tavg_ST)) then
            call accumulate_tavg_field(TRACER(:,:,k,1,curtime,iblock)* &
                                       TRACER(:,:,k,2,curtime,iblock), &
                                       tavg_ST,iblock,k)
         endif

         if (tavg_requested(tavg_RHO)) then
            call accumulate_tavg_field(RHO(:,:,k,curtime,iblock), &
                                       tavg_RHO,iblock,k)
         endif

        if ( sfc_layer_type /= sfc_layer_varthick .and. k == 1) then
          if (tavg_requested(tavg_RESID_T)) then
              WORK1 = c0
              factor = c1/hflux_factor  ! converts to W/m^2
              where (CALCT(:,:,iblock))  &
                WORK1=DH(:,:,iblock)*TRACER(:,:,1,1,curtime,iblock)*factor
              call accumulate_tavg_field(WORK1,tavg_RESID_T,iblock,k)
          endif

          if (tavg_requested(tavg_RESID_S)) then
              WORK1 = c0
              factor = c1/salinity_factor  ! converts to kg(freshwater)/m^2/s
              where (CALCT(:,:,iblock)) &
                WORK1 = DH(:,:,iblock)*TRACER(:,:,k,2,curtime,iblock)*factor
              call accumulate_tavg_field(WORK1,tavg_RESID_S,iblock,k)
          endif
        endif  ! sfc_layer_type

         if (POP_nt > 2) call tavg_passive_tracers(iblock,k)

!-----------------------------------------------------------------------
!
!        update movie fields if requested
!
!-----------------------------------------------------------------------

         if (movie_requested(movie_UVEL(k))) then
            call update_movie_field(UVEL(:,:,k,curtime,iblock), &
                                       movie_UVEL(k),iblock,k)
         endif

         if (movie_requested(movie_VVEL(k))) then
            call update_movie_field(VVEL(:,:,k,curtime,iblock), &
                                       movie_VVEL(k),iblock,k)
         endif

         if (movie_requested(movie_TEMP(k))) then
            call update_movie_field(TRACER(:,:,k,1,curtime,iblock), &
                                       movie_TEMP(k),iblock,k)
         endif

         if (movie_requested(movie_SALT(k))) then
            call update_movie_field(  &
                 TRACER(:,:,k,2,curtime,iblock)*salt_to_ppt, &
                                       movie_SALT(k),iblock,k)  !  convert to psu
         endif

         if (movie_requested(movie_RHO(k))) then
            call update_movie_field(  &
               salt_to_ppt*(RHO(:,:,k,curtime,iblock)-c1), movie_RHO(k),iblock,k)

         endif

         endif ! mix_pass

!-----------------------------------------------------------------------

      enddo  ! k loop

!-----------------------------------------------------------------------
!  put vertical integral movie calls after k loop
!  only UVEL, VVEL defined since no obvious need for others
!-----------------------------------------------------------------------

      if (movie_requested(movie_UVEL(0)) .and. mix_pass /= 1) then
         WORK1 = c0
         do k = 1, POP_km
            WORK1 = WORK1 + UVEL(:,:,k,curtime,iblock)*dzr(k)
         enddo
         call update_movie_field(WORK1,movie_UVEL(0),iblock,0)
      endif
      if (movie_requested(movie_VVEL(0)) .and. mix_pass /= 1) then
         WORK1 = c0
         do k = 1, POP_km
            WORK1 = WORK1 + VVEL(:,:,k,curtime,iblock)*dzr(k)
         enddo
         call update_movie_field(WORK1,movie_VVEL(0),iblock,0)
      endif

!-----------------------------------------------------------------------
!  do all levels of passive tracer movies as well as vertical integrals
!-----------------------------------------------------------------------

      if ( POP_nt > 2 .and. mix_pass /= 1) call movie_passive_tracers(iblock)

!-----------------------------------------------------------------------
!
!        update float fields if requested
!
!-----------------------------------------------------------------------

      if (float_requested(float_TEMP)) then
         call update_float_buffer(TRACER(:,:,:,1,curtime,iblock), &
                                    float_TEMP,iblock)
      endif

      if (float_requested(float_SALT)) then
         call update_float_buffer(  &
              TRACER(:,:,:,2,curtime,iblock)*salt_to_ppt, &
                                    float_SALT,iblock)  !  convert to psu
      endif

      if (float_requested(float_UVEL)) then
         call update_float_buffer(UVEL(:,:,:,curtime,iblock), &
                                    float_UVEL,iblock)
      endif

      if (float_requested(float_VVEL)) then
         call update_float_buffer(VVEL(:,:,:,curtime,iblock), &
                                    float_VVEL,iblock)
      endif

!-----------------------------------------------------------------------
!
!     if no implicit vertical mixing, we now have updated tracers
!     if using implicit vertical mixing and rigid lid or old free
!        surface form, update all the tracers now using impvmix
!     if using implicit vertical mixing and a variable thickness
!        surface layer, only update T,S as predicted values to
!        use for pressure averaging - the full update will
!        occur after the barotropic solver
!
!-----------------------------------------------------------------------

      if (implicit_vertical_mix) then
         if (sfc_layer_type /= sfc_layer_varthick) then
            call impvmixt(TRACER(:,:,:,:,newtime,iblock), &
                          TRACER(:,:,:,:,oldtime,iblock), &
                          PSURF (:,:,    curtime,iblock), &
                          1, POP_nt, this_block)

         else if (lpressure_avg .and. leapfrogts) then

            !*** predictor update of T,S
            !*** with PSURF(curtime) on the LHS at k=1
 
            call impvmixt(TRACER(:,:,:,:,newtime,iblock), &
                          TRACER(:,:,:,:,oldtime,iblock), &
                          PSURF (:,:,curtime,iblock),     &
                          1, 2, this_block)

         endif
      endif

!-----------------------------------------------------------------------
!
!     end of first block loop
!
!-----------------------------------------------------------------------

   enddo ! first block loop

!jw   !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!  update non-prognostic float fields if necessary
!-----------------------------------------------------------------------

   if (float_requested(float_PD)) then

      allocate (WORK3D1(nx_block,ny_block,POP_km))
      do iblock = 1,nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)
         do k = 1, POP_km
            call state(k,1,TRACER(:,:,k,1,curtime,iblock),                         &
                        TRACER(:,:,k,2,curtime,iblock), this_block, &
                        RHOFULL=WORK1)
            WORK3D1(:,:,k) = salt_to_ppt*(WORK1 - c1)
         enddo
         call update_float_buffer(WORK3D1,float_PD,iblock)
      enddo

      deallocate (WORK3D1)
   endif

   if (float_requested(float_PV)) then

      allocate (WORK3D2(nx_block,ny_block,POP_km))
      do iblock = 1,nblocks_clinic
         this_block = get_block(blocks_clinic(iblock),iblock)
         do k = 1, POP_km

! the following is copied from calculation of tavg_PV in advection.F90

            if (k == 1 ) then
               FX = RHO(:,:,k,curtime,iblock)
            else
               call state(k-1,k,TRACER(:,:,k-1,1,curtime,iblock)  &
                               ,TRACER(:,:,k-1,2,curtime,iblock), this_block, &
                                RHOOUT=FY)
               FX = p5*(FY + RHO(:,:,k,curtime,iblock))
            endif

            if (k == POP_km) then
               FY = RHO(:,:,k,curtime,iblock)
            else
               call state(k+1,k,TRACER(:,:,k+1,1,curtime,iblock)  &
                               ,TRACER(:,:,k+1,2,curtime,iblock), this_block, &
                                RHOOUT=WORK2) 

               do j=1, ny_block
               do i=1, nx_block
                  if (k /= KMT(i,j,iblock)) then
                     FY(i,j) = p5*(WORK2(i,j) + RHO(i,j,k,curtime,iblock))
                  else
                     FY(i,j) = RHO(i,j,k,curtime,iblock)
                  endif
               end do
               end do
            endif

            do j=1, ny_block
            do i=1, nx_block
               if (k <= KMT(i,j,iblock)) then
                  WORK2(i,j) = (FX(i,j) - FY(i,j))*dzr(k) ! drho/dz
               else
                  WORK2(i,j) = c0
               endif
            end do
            end do

            call zcurl(k,FX,UVEL(:,:,k,curtime,iblock),   &
                            VVEL(:,:,k,curtime,iblock),this_block)
            WORK3D2(:,:,k) = WORK2*(FX*TAREA_R(:,:,iblock) + &
                                    FCORT(:,:,iblock)) ! PV = pot vorticity
         enddo  !  k
         call update_float_buffer(WORK3D2,float_PV,iblock)
      enddo  !  nblocks_clinic

      deallocate (WORK3D2)

   endif  !  float PV requested

!-----------------------------------------------------------------------
!
!  update tracer ghost cells here outside the block loops (it
!  requires synchronization) - only update T,S ghost cells
!  for use in pressure averaging.
!
!-----------------------------------------------------------------------

   if (lpressure_avg .and. leapfrogts) then

      call POP_HaloUpdate(TRACER(:,:,:,1,newtime,:),          &
                              POP_haloClinic,                 &
                              POP_gridHorzLocCenter,          &
                              POP_fieldKindScalar, errorCode, &
                              fillValue = 0.0_POP_r8)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'baroclinic_driver: error updating halo for PT')
         return
      endif

      call POP_HaloUpdate(TRACER(:,:,:,2,newtime,:),          &
                              POP_haloClinic,                 &
                              POP_gridHorzLocCenter,          &
                              POP_fieldKindScalar, errorCode, &
                              fillValue = 0.0_POP_r8)

      if (errorCode /= POP_Success) then
         call POP_ErrorSet(errorCode, &
            'baroclinic_driver: error updating halo for salinity')
         return
      endif

   endif

!-----------------------------------------------------------------------
!
!  now loop over blocks to do momentum equations
!
!-----------------------------------------------------------------------

   !$OMP PARALLEL DO PRIVATE(iblock,k,km1,kp1,n,this_block, &
   !$OMP                     WUK, FX, FY, WORK1, WORK2)

   do iblock = 1,nblocks_clinic

      this_block = get_block(blocks_clinic(iblock),iblock)  

!-----------------------------------------------------------------------
!
!     initialize arrays for vertical sums.
!
!-----------------------------------------------------------------------

      ZX(:,:,iblock) = c0
      ZY(:,:,iblock) = c0

      do k = 1,POP_km

         kp1 = k+1
         km1 = k-1
         if (k == 1) km1 = 1
         if (k == POP_km) kp1 = POP_km

!-----------------------------------------------------------------------
!
!        if pressure averaging is on and it is a leapfrog time step
!        we need the updated density for the pressure averaging
!
!-----------------------------------------------------------------------

         if (lpressure_avg .and. leapfrogts) then
            call state(k,k,TRACER(:,:,k,1,newtime,iblock), &
                           TRACER(:,:,k,2,newtime,iblock), &
                           this_block, RHOOUT=RHO(:,:,k,newtime,iblock))
         endif

!-----------------------------------------------------------------------
!
!        calculate forcing terms (Fx,Fy) at level k.
!
!-----------------------------------------------------------------------

         call clinic(k, FX, FY, WUK,                &
                        UVEL(:,:,:,curtime,iblock), &
                        VVEL(:,:,:,curtime,iblock), &
                        UVEL(:,:,:,oldtime,iblock), &
                        VVEL(:,:,:,oldtime,iblock), &
                        UVEL(:,:,k,mixtime,iblock), &
                        VVEL(:,:,k,mixtime,iblock), &
                        RHO (:,:,k,oldtime,iblock), &
                        RHO (:,:,k,curtime,iblock), &
                        RHO (:,:,k,newtime,iblock), &
                        SMF (:,:,:,iblock),         &
                        DHU (:,:,iblock),           &
                        this_block, errorCode)

         if (errorCode /= POP_Success) then
            call POP_ErrorSet(errorCode, &
               'baroclinic_driver: error in clinic')
            return
         endif

!-----------------------------------------------------------------------
!
!        store forces temporarily in UVEL(newtime),VVEL(newtime).
!
!-----------------------------------------------------------------------

         if (impcor) then   ! implicit treatment

            WORK1 = c2dtu*beta*FCOR(:,:,iblock)
            WORK2 = c2dtu/(c1 + WORK1**2)
            UVEL(:,:,k,newtime,iblock) = (FX + WORK1*FY)*WORK2 
            VVEL(:,:,k,newtime,iblock) = (FY - WORK1*FX)*WORK2 



         else               ! explicit treatment

            UVEL(:,:,k,newtime,iblock) = c2dtu*FX
            VVEL(:,:,k,newtime,iblock) = c2dtu*FY

         endif

!-----------------------------------------------------------------------
!
!        increment sum for vertically-averaged forcing ([Fx],[Fy]).
!
!-----------------------------------------------------------------------

         if (partial_bottom_cells) then
            ZX(:,:,iblock) = ZX(:,:,iblock) + FX*DZU(:,:,k,iblock)
            ZY(:,:,iblock) = ZY(:,:,iblock) + FY*DZU(:,:,k,iblock)
         else
            ZX(:,:,iblock) = ZX(:,:,iblock) + FX*dz(k)
            ZY(:,:,iblock) = ZY(:,:,iblock) + FY*dz(k)
         endif

      enddo ! vertical (k) loop

!-----------------------------------------------------------------------
!
!     normalize sums for vertical averages ([Fx],[Fy]) by dividing
!     by depth at U points.
!
!-----------------------------------------------------------------------

      ZX(:,:,iblock) = ZX(:,:,iblock)*HUR(:,:,iblock)
      ZY(:,:,iblock) = ZY(:,:,iblock)*HUR(:,:,iblock)

!-----------------------------------------------------------------------
!
!     solve tridiagonal system with implicit treatment of vertical 
!     diffusion of velocity.
!
!-----------------------------------------------------------------------

      if (implicit_vertical_mix)                   &
         call impvmixu(UVEL(:,:,:,newtime,iblock), &
                       VVEL(:,:,:,newtime,iblock), & 
                       this_block)

!-----------------------------------------------------------------------
!
!     calculate unnormalized baroclinic velocities (Upp,Vpp)
!
!-----------------------------------------------------------------------

      UVEL(:,:,:,newtime,iblock) = UVEL(:,:,:,oldtime,iblock) + &
                                   UVEL(:,:,:,newtime,iblock)  ! holds c2dtu*Fx
      VVEL(:,:,:,newtime,iblock) = VVEL(:,:,:,oldtime,iblock) + &
                                   VVEL(:,:,:,newtime,iblock)  ! holds c2dtu*Fy

!-----------------------------------------------------------------------
!
!     find vertical averages ([Upp],[Vpp]).
!
!-----------------------------------------------------------------------

      WORK1 = c0  ! initialize sums
      WORK2 = c0

      if (partial_bottom_cells) then
         do k = 1,POP_km
            WORK1 = WORK1 + UVEL(:,:,k,newtime,iblock)*DZU(:,:,k,iblock)
            WORK2 = WORK2 + VVEL(:,:,k,newtime,iblock)*DZU(:,:,k,iblock)
         enddo
      else
         do k = 1,POP_km
            WORK1 = WORK1 + UVEL(:,:,k,newtime,iblock)*dz(k)
            WORK2 = WORK2 + VVEL(:,:,k,newtime,iblock)*dz(k)
         enddo
      endif

      WORK1 = WORK1*HUR(:,:,iblock)  ! normalize by dividing by depth
      WORK2 = WORK2*HUR(:,:,iblock)

!-----------------------------------------------------------------------
!
!     normalize baroclinic velocites by subtracting vertical mean:
!     (Up,Vp) = (Upp,Vpp) - ([Upp],[Vpp]), zero velocities at land pts.
!
!-----------------------------------------------------------------------

      do k = 1,POP_km
         where (k <= KMU(:,:,iblock))
            UVEL(:,:,k,newtime,iblock) = &
            UVEL(:,:,k,newtime,iblock) - WORK1
            VVEL(:,:,k,newtime,iblock) = &
            VVEL(:,:,k,newtime,iblock) - WORK2
         elsewhere 
            UVEL(:,:,k,newtime,iblock) = c0
            VVEL(:,:,k,newtime,iblock) = c0
         endwhere
      enddo

!-----------------------------------------------------------------------
!
!     note:  at this point UVEL(newtime) and VVEL(newtime) contain only 
!     the baroclinic velocities (Up,Vp) at the new time.  they are later
!     updated to the full velocities in step after barotropic is
!     is called, which calculates the barotropic velocites ([U],[V])
!     at the new time.  UVEL and VVEL at time levels oldtime ond curtime 
!     always contain the full velocites.
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!  end of second block loop over vertical levels.
!
!-----------------------------------------------------------------------

   enddo ! second block loop

   !$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!
!  check cfl limits
!
!-----------------------------------------------------------------------

   if (ldiag_cfl) then
      call cfl_check
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine baroclinic_driver

!***********************************************************************
!BOP
! !IROUTINE: baroclinic_correct_adjust
! !INTERFACE:

 subroutine baroclinic_correct_adjust(errorCode)

! !DESCRIPTION:
!  This subroutine finishes updating tracers by performing
!  adjustment-like physics (convection and ice) and completes
!  the corrector step for tracers using the new surface pressure
!  to update tracers.
!
! !REVISION HISTORY:
!  same as module

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode             ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) ::  & 
      k,                  &! vertical level index
      n,                  &! tracer index
      iblock               ! block index

   real (POP_r8), dimension(nx_block,ny_block,POP_nt) :: &
      RHS1                 ! r.h.s. for impvmix on corrector step

   type (block) ::       &
      this_block           ! block information for current block

!-----------------------------------------------------------------------
!
!  increment tlast_ice outside of threaded region
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   call increment_tlast_ice

!-----------------------------------------------------------------------
!
!  do everything for each sub block
!
!-----------------------------------------------------------------------

   do iblock = 1,nblocks_clinic

      this_block = get_block(blocks_clinic(iblock),iblock)  

!-----------------------------------------------------------------------
!
!     do the corrector step for variable thickness surface layer
!
!-----------------------------------------------------------------------

      if (sfc_layer_type == sfc_layer_varthick) then

         if (implicit_vertical_mix) then

            !*** if implicit vertical mixing and pressure averaging:
            !*** correct new T and S and update remaining passive
            !*** tracers with tridiagonal solves

            if (lpressure_avg .and. leapfrogts) then
        
               do n = 1,2   ! corrector for T,S only
                  where (KMT(:,:,iblock) > 0)  ! corrector RHS at k=1
                     RHS1(:,:,n)=((c2*TRACER(:,:,1,n,curtime,iblock) - &
                                      TRACER(:,:,1,n,oldtime,iblock))  &
                                 *(PSURF(:,:,curtime,iblock) -         &
                                   PSURF(:,:,oldtime,iblock))          &
                                - TRACER(:,:,1,n,newtime,iblock)*      &
                                  (PSURF(:,:,newtime,iblock) -         &
                                   PSURF(:,:,curtime,iblock)))/        &
                                   (grav*dz(1)) 
                  elsewhere
                     RHS1(:,:,n) = c0
                  endwhere
               enddo

               !*** T,S update on corrector step
               call impvmixt_correct(TRACER(:,:,:,:,newtime,iblock), &
                                     PSURF (:,:,    newtime,iblock), &
                                     RHS1, 1, 2, this_block)

               do n = 3,POP_nt
                  !*** surface RHS for passive tracers with press avg
                  where (KMT(:,:,iblock) > 0)
                     TRACER(:,:,1,n,newtime,iblock) =               &
                                   TRACER(:,:,1,n,newtime,iblock) - &
                                   TRACER(:,:,1,n,oldtime,iblock)   & 
                                   *(PSURF(:,:,newtime,iblock) -    &
                                     PSURF(:,:,oldtime,iblock))/    &
                                     (grav*dz(1))
                  endwhere
               enddo

               !*** standard update of all passive tracers
               !***  n=3,POP_nt with PSURF(newtime) on LHS at k=1

               call impvmixt(TRACER(:,:,:,:,newtime,iblock), &
                             TRACER(:,:,:,:,oldtime,iblock), &
                             PSURF (:,:,    newtime,iblock), &
                             3, POP_nt, this_block)

            !*** if implicit vertical mixing but no pressure averaging
            !*** update all tracers with tridiagonal solves

            else  ! no pressure averaging or not leapfrog

               do n = 1,POP_nt
                  !***  surface RHS for tracers with pressure avg
                  where (KMT(:,:,iblock) > 0)
                     TRACER(:,:,1,n,newtime,iblock) = &
                     TRACER(:,:,1,n,newtime,iblock) - &
                     TRACER(:,:,1,n,oldtime,iblock)   & 
                     *(PSURF(:,:,newtime,iblock) -    &
                       PSURF(:,:,mixtime,iblock))/(grav*dz(1))
                  endwhere
               enddo

               !*** standard update all tracers:  n=1,POP_nt
               !*** with PSURF(newtime) on the LHS at k=1
               call impvmixt(TRACER(:,:,:,:,newtime,iblock), &
                             TRACER(:,:,:,:,oldtime,iblock), &
                             PSURF (:,:,    newtime,iblock), &
                             1, POP_nt, this_block)

            endif  ! pressure averaging and leapfrog

         else ! no implicit_vertical_mix

            !*** if explicit vertical mixing and pressure averaging:
            !*** correct new tracers level k=1

            if (lpressure_avg .and. leapfrogts) then

               do n = 1,2
                  !*** correct surface T and S with pressure avg
                  where (KMT(:,:,iblock) > 0)
                     TRACER(:,:,1,n,newtime,iblock) =           &
                       (TRACER(:,:,1,n,newtime,iblock)*         &
                       (dz(1) + PSURF(:,:,curtime,iblock)/grav) &
                       +(c2*TRACER(:,:,1,n,curtime,iblock) -    &
                            TRACER(:,:,1,n,oldtime,iblock))     &
                       *(PSURF(:,:,curtime,iblock) -            &
                         PSURF(:,:,oldtime,iblock))/grav)       &
                       /(dz(1) + PSURF(:,:,newtime,iblock)/grav) 
                  elsewhere
                     TRACER(:,:,1,n,newtime,iblock) = c0 ! zero on land pts
                  endwhere
               enddo

               do n = 3,POP_nt
                  where (KMT(:,:,iblock) > 0)
                     TRACER(:,:,1,n,newtime,iblock) =             &
                        (TRACER(:,:,1,n,oldtime,iblock)*          &
                        (dz(1) + PSURF(:,:,oldtime,iblock)/grav)  &
                        + dz(1)*TRACER(:,:,1,n,newtime,iblock))   &
                        /(dz(1) + PSURF(:,:,newtime,iblock)/grav)
                  elsewhere
                     TRACER(:,:,1,n,newtime,iblock) = c0 ! zero on land pts
                  endwhere
               enddo

            else  ! no pressure avg or leapfrog

               do n = 1,POP_nt

                  !*** exact update of all tracers in surface layer
                  where (KMT(:,:,iblock) > 0)
                     TRACER(:,:,1,n,newtime,iblock) =              &
                        (TRACER(:,:,1,n,oldtime,iblock)*           &
                        (dz(1) + PSURF(:,:,mixtime,iblock)/grav)   &
                        + dz(1)*TRACER(:,:,1,n,newtime,iblock))/   &
                        (dz(1) + PSURF(:,:,newtime,iblock)/grav)
                  elsewhere
                     TRACER(:,:,1,n,newtime,iblock) = c0 ! zero on land pts
                  endwhere

               enddo

            endif  ! pressure avg and leapfrog

         endif ! implicit_vertical_mix

      endif ! variable thickness surface layer

!-----------------------------------------------------------------------
!
!     check for surface temperatures below freezing
!     do not reset if ice formation option is on
!
!-----------------------------------------------------------------------

      if (reset_to_freezing .and. .not. liceform) then
         TRACER(:,:,1,1,newtime,iblock) = &
            max(TRACER(:,:,1,1,newtime,iblock),-c2)
      endif

!-----------------------------------------------------------------------
!
!     convective adjustment of tracers - 
!     convad routine does nothing if convective adjustment not chosen
!     otherwise it performs convective adjustment and recomputes
!     density
!
!-----------------------------------------------------------------------

      call convad(TRACER(:,:,:,:,newtime,iblock), &
                  RHO(:,:,:,newtime,iblock), this_block, iblock)

!-----------------------------------------------------------------------
!
!     compute ice formation and adjust temperature due to ice formation
!     if this option was requested
!
!-----------------------------------------------------------------------

      if (liceform .and. mix_pass /= 1) then
         call ice_formation(TRACER(:,:,:,:,newtime,iblock),          &
                            STF(:,:,1,iblock) + SHF_QSW(:,:,iblock), &
                            iblock,this_block,lfw_as_salt_flx)
      endif

!-----------------------------------------------------------------------
!
!     call passive tracer reset subroutines
!
!-----------------------------------------------------------------------

      if (POP_nt > 2) call reset_passive_tracers(  &
         TRACER(:,:,:,:,newtime,iblock))


!-----------------------------------------------------------------------
!
!     compute new density based on new tracers
!     use GPU acceleration if possible
!
!-----------------------------------------------------------------------
      if (use_gpu_state .and. state_range_iopt == state_range_enforce .and. state_itype == state_type_mwjf) then
        write(stdout,'(a25)') 'Going to run state on GPU'

        call mwjf_state(TRACER(:,:,:,1,newtime,iblock), &
                        TRACER(:,:,:,2,newtime,iblock), &
                        1, POP_km, &
                        RHOOUT=RHO(:,:,:,newtime,iblock))

        write(stdout,'(a21)') 'Finished state on GPU, checking result'

        do k = 1,POP_km  ! recalculate new density
          call state(k,k,TRACER(:,:,k,1,newtime,iblock), &
                         TRACER(:,:,k,2,newtime,iblock), &
                         this_block, RHOOUT=RHOREF(:,:,k))

        enddo

        write(stdout,'(a25)') 'Finished state on CPU'

        call gpu_compare(RHO, RHOREF, nx_block*ny_block*POP_km)

      else
        do k = 1,POP_km  ! recalculate new density
          call state(k,k,TRACER(:,:,k,1,newtime,iblock), &
                         TRACER(:,:,k,2,newtime,iblock), &
                         this_block, RHOOUT=RHO(:,:,k,newtime,iblock))

        enddo
      endif

!-----------------------------------------------------------------------
!
!  end of block loop
!
!-----------------------------------------------------------------------

   end do

!-----------------------------------------------------------------------
!EOC

 end subroutine baroclinic_correct_adjust

!***********************************************************************
!BOP
! !IROUTINE: clinic
! !INTERFACE:

 subroutine clinic(k, FX, FY, WUK, UCUR, VCUR, UOLD, VOLD,      &
                      UMIXK, VMIXK, RHOKOLD, RHOKCUR, RHOKNEW,  &
                      SMF_BLOCK, DHU_BLOCK, this_block,errorCode)

! !DESCRIPTION:
!  Calculates forcing terms on r.h.s. of baroclinic momentum eqns.
!  \begin{eqnarray}
!     F_x &=& -L(u) + fv - \nabla p + D_H(u^{n-1}) + D_V(u^{n-1}) \\
!     F_y &=& -L(v) - fu - \nabla p + D_H(v^{n-1}) + D_V(v^{n-1})
!  \end{eqnarray}
!
!  The above equations are written for the case of explicit
!  treatment of the Coriolis terms.  If these terms are treated
!  implicitly, then the coriolis terms above should be replaced by:
!  \begin{eqnarray}
!       +fv &\rightarrow& +f(\gamma v + (1-\gamma)v^{n-1}) \\
!       -fu &\rightarrow& -f(\gamma u + (1-\gamma)u^{n-1})
!  \end{eqnarray}
!  on leapfrog timesteps and
!  \begin{eqnarray}
!       +fv &\rightarrow& +fv^{n-1} \\
!       -fu &\rightarrow& -fu^{n-1}
!  \end{eqnarray}
!  on Matsuno timesteps, where $\gamma$ is a parameter used to vary 
!  the time-centering of the Coriolis and pressure gradient terms on 
!  leapfrog steps.
!
!  The small metric terms for advection and diffusion of the
!  velocity field are calculated in the advection and horizontal 
!  diffusion routines.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: k    ! depth level index

   real (POP_r8), dimension(nx_block,ny_block,POP_km), intent(in) :: &
      UCUR, VCUR,           &! U,V for block at current time
      UOLD, VOLD             ! U,V for block at old     time

   real (POP_r8), dimension(nx_block,ny_block), intent(in) :: &
      UMIXK, VMIXK,         &! U,V at level k and mix time level
      RHOKOLD,              &! density at level k and mix time level
      RHOKCUR,              &! density at level k and cur time level
      RHOKNEW,              &! density at level k and new time level
      DHU_BLOCK              ! change in surface height at U pts

   real (POP_r8), dimension(nx_block,ny_block,2), intent(in) :: &
      SMF_BLOCK              ! surface momentum forcing for this block

   type (block), intent(in) :: &
      this_block             ! block info for the current block

! !INPUT/OUTPUT PARAMETERS:

   real (POP_r8), dimension(nx_block,ny_block), intent(inout) :: &
      WUK             ! vertical velocity at top of U box

! !OUTPUT PARAMETERS:

   real (POP_r8), dimension(nx_block,ny_block), intent(out) :: &
      FX,       &! sum of terms contributing to Fx at level k 
      FY         ! sum of terms contributing to Fy at level k 

   integer (POP_i4), intent(out) :: &
      errorCode             ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables:
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      bid                ! local block address

   real (POP_r8), dimension(nx_block,ny_block) :: &
      WORKX,WORKY     ! local work space to hold forcing terms

!-----------------------------------------------------------------------
!
!  advection L(U),L(V)
!  set vertical velocity at surface
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   bid = this_block%local_id
 
   if (k == 1) WUK = DHU_BLOCK  ! free surface

   call advu(k, WORKX, WORKY, WUK, UCUR, VCUR, this_block)

   FX =  -WORKX   ! advu returns WORKX = +L(U) 
   FY =  -WORKY   ! advu returns WORKY = +L(V) 

   if (ldiag_global) then
      if (partial_bottom_cells) then
         DIAG_KE_ADV_2D(:,:,bid) = DIAG_KE_ADV_2D(:,:,bid) -           &
                                   DZU(:,:,k,bid)*(UCUR(:,:,k)*WORKX + &
                                                   VCUR(:,:,k)*WORKY)
      else
         DIAG_KE_ADV_2D(:,:,bid) = DIAG_KE_ADV_2D(:,:,bid) -  &
                                   dz(k)*(UCUR(:,:,k)*WORKX + & 
                                          VCUR(:,:,k)*WORKY)
      endif
   endif

!-----------------------------------------------------------------------
!
!  coriolis terms
!
!-----------------------------------------------------------------------

   if (impcor .and. leapfrogts) then          ! implicit, leapfrog

      FX = FX + FCOR(:,:,bid)*(      gamma* VCUR(:,:,k) + &
                               (c1 - gamma)*VOLD(:,:,k))
      FY = FY - FCOR(:,:,bid)*(      gamma* UCUR(:,:,k) + & 
                               (c1 - gamma)*UOLD(:,:,k))

   elseif(.not.impcor .and. leapfrogts) then  ! explicit, leapfrog

      FX = FX + FCOR(:,:,bid)*VCUR(:,:,k)
      FY = FY - FCOR(:,:,bid)*UCUR(:,:,k)

   else                                  ! matsuno or foward euler

      FX = FX + FCOR(:,:,bid)*VOLD(:,:,k)
      FY = FY - FCOR(:,:,bid)*UOLD(:,:,k)

   endif

!-----------------------------------------------------------------------
!
!  hydrostatic pressure gradients
!
!-----------------------------------------------------------------------

   call gradp(k,WORKX, WORKY, RHOKOLD, RHOKCUR, RHOKNEW, this_block)

   FX = FX - WORKX   ! gradp returns WORKX as +Gradx(p)
   FY = FY - WORKY   ! gradp returns WORKY as +Grady(p)

   if (partial_bottom_cells) then
      WORKX =  -DZU(:,:,k,bid)*(UCUR(:,:,k)*WORKX + &
                                VCUR(:,:,k)*WORKY)
   else
      WORKX =  -dz(k)*(UCUR(:,:,k)*WORKX + & 
                       VCUR(:,:,k)*WORKY)
   endif

   if (tavg_requested(tavg_UDP) .and. mix_pass /= 1) then
      call accumulate_tavg_field(WORKX,tavg_UDP,bid,k)
   endif

   if (ldiag_global) then
      DIAG_KE_PRESS_2D(:,:,bid) = DIAG_KE_PRESS_2D(:,:,bid) + WORKX
   endif

!-----------------------------------------------------------------------
!
!  horizontal diffusion HDiff(Ub),HDiff(Vb)
!
!-----------------------------------------------------------------------

   call hdiffu(k, WORKX, WORKY, UMIXK, VMIXK, this_block, errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'clinic: error in hdiffu')
      return
   endif

   FX = FX + WORKX
   FY = FY + WORKY

   if (ldiag_global) then
      if (partial_bottom_cells) then
         DIAG_KE_HMIX_2D(:,:,bid) = DIAG_KE_HMIX_2D(:,:,bid) +        &
                                  DZU(:,:,k,bid)*(UCUR(:,:,k)*WORKX + &
                                                  VCUR(:,:,k)*WORKY)
      else
         DIAG_KE_HMIX_2D(:,:,bid) = DIAG_KE_HMIX_2D(:,:,bid) + &
                                    dz(k)*(UCUR(:,:,k)*WORKX + & 
                                           VCUR(:,:,k)*WORKY)
      endif
   endif

!-----------------------------------------------------------------------
!
!  vertical diffusion VDiff(Ub),VDiff(Vb)
!
!-----------------------------------------------------------------------

   call vdiffu(k, WORKX, WORKY, UOLD, VOLD, SMF_BLOCK, this_block)

   FX = FX + WORKX
   FY = FY + WORKY


   if (ldiag_global) then
      if (partial_bottom_cells) then
         DIAG_KE_VMIX_2D(:,:,bid) = DIAG_KE_VMIX_2D(:,:,bid) + &
                           DZU(:,:,k,bid)*(UCUR(:,:,k)*WORKX + &
                                           VCUR(:,:,k)*WORKY)
      else
         DIAG_KE_VMIX_2D(:,:,bid) = DIAG_KE_VMIX_2D(:,:,bid) + &
                                    dz(k)*(UCUR(:,:,k)*WORKX + & 
                                           VCUR(:,:,k)*WORKY)
      endif
   endif

!-----------------------------------------------------------------------
!
!  zero forces (and hence velocities) at land points
!
!-----------------------------------------------------------------------

   where (k > KMU(:,:,bid))
      FX = c0
      FY = c0
   endwhere

!-----------------------------------------------------------------------
!EOC

 end subroutine clinic

!***********************************************************************
!BOP
! !IROUTINE: tracer_update
! !INTERFACE:

 subroutine tracer_update(k, WTK, TNEW, TOLD, TMIX, TCUR,            &
                             UCUR, VCUR, UMIX, VMIX, RHOCUR,         &
                             STF_IN, TFW_IN, QSW, DH_IN, POLD, PCUR, &
                             this_block, errorCode)

! !DESCRIPTION:
!  Computes explicit forcing for tracer equations:
!  \begin{equation}
!     (T^{n+1}-T^{n-1})/(2 \Delta t) = -L(T) + D_H(T^{n-1}) + 
!                                              D_V(T^{n-1}) + S
!  \end{equation}
!  where $L$ is the advection operator, $D_{H,V}$ are the diffusion
!  operators in the horizontal and vertical and $S$ are source terms.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: k     ! depth level index

   real (POP_r8), dimension(nx_block,ny_block,POP_km,POP_nt), intent(in) :: &
      TCUR,                 &! tracers at current time level
      TOLD,                 &! tracers at old     time level
      TMIX                   ! tracers at mix     time level

   real (POP_r8), dimension(nx_block,ny_block,POP_km), intent(in) :: &
      UCUR, VCUR,           &! U,V for block at current time
      UMIX, VMIX,           &! U,V at mix time level
      RHOCUR                 ! density at current time level

   real (POP_r8), dimension(nx_block,ny_block,POP_nt), intent(in) :: &
      STF_IN,               &! surface tracer fluxes
      TFW_IN                 ! tracer concentration in fresh water

   real (POP_r8), dimension(nx_block,ny_block), intent(in) :: &
      DH_IN,                &! sfc height change at tracer points
      POLD,                 &! sfc pressure at old     time
      PCUR,                 &! sfc pressure at current time
      QSW                    ! short-wave heating

   type (block), intent(in) :: &
      this_block             ! block info for the current block

! !INPUT/OUTPUT PARAMETERS:

   real (POP_r8), dimension(nx_block,ny_block), intent(inout) :: &
      WTK          ! on  input, vertical velocity at top    of T box
                   ! on output, vertical velocity at bottom of T box

   real (POP_r8), dimension(nx_block,ny_block,POP_km,POP_nt), intent(inout) :: &
      TNEW                   ! tracers at new time level

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode             ! returned error code

!EOP
!BOC
!-----------------------------------------------------------------------
!         
!  local variables:
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: &
      n,                 &! dummy tracer index
      bid                 ! local_block id

   real (POP_r8), dimension(nx_block,ny_block,POP_nt) :: &
      FT,                &! sum of terms in dT/dt for the nth tracer
      WORKN               ! work array used for various dT/dt terms 

   real (POP_r8), dimension(nx_block,ny_block) :: &
      WORKSW

!-----------------------------------------------------------------------
!
!  initialize some arrays
!
!-----------------------------------------------------------------------

   errorCode = POP_Success

   FT    = c0

   bid = this_block%local_id

!-----------------------------------------------------------------------
!
!  horizontal diffusion HDiff(T)
!
!-----------------------------------------------------------------------


   call hdifft(k, WORKN, TMIX, UMIX, VMIX, this_block, errorCode)

   if (errorCode /= POP_Success) then
      call POP_ErrorSet(errorCode, &
         'tracer_update: error in hdifft')
      return
   endif

   FT = FT + WORKN

   if (ldiag_global) then

      if (partial_bottom_cells) then
         do n=1,POP_nt
            where (k <= KMT(:,:,bid))            &
               DIAG_TRACER_HDIFF_2D(:,:,n,bid) = &
               DIAG_TRACER_HDIFF_2D(:,:,n,bid) + &
               WORKN(:,:,n)*DZT(:,:,k,bid)
         end do
      else
         do n=1,POP_nt


            where (k <= KMT(:,:,bid))            &
               DIAG_TRACER_HDIFF_2D(:,:,n,bid) = & 
               DIAG_TRACER_HDIFF_2D(:,:,n,bid) + &
               WORKN(:,:,n)*dz(k)
         end do
      endif
   endif

!-----------------------------------------------------------------------
!
!  advection L(T)
!  set vertical velocity at surface
!  w = DH = dh/dt or dh/dt - Fw depending on surface type
!
!-----------------------------------------------------------------------


   if (k == 1) WTK = DH_IN

   if (ldiag_global) then

      if (k == 1) then
         where (CALCT(:,:,bid)) &
            DIAG_PE_2D(:,:,bid) = dzw(0)*grav*RHOCUR(:,:,1)*WTK

         if (sfc_layer_type /= sfc_layer_varthick) then
            do n=1,POP_nt
               where (CALCT(:,:,bid))
                  DIAG_TRACER_ADV_2D (:,:,n,bid) = WTK*TCUR(:,:,1,n)
                  DIAG_TRACER_SFC_FLX(:,:,n,bid) = WTK*TCUR(:,:,1,n)
               elsewhere
                  DIAG_TRACER_ADV_2D (:,:,n,bid) = c0
                  DIAG_TRACER_SFC_FLX(:,:,n,bid) = c0
               end where
            end do
         endif

      else

         !*** For energetic consistency, we use dzw even for 
         !*** partial bottom cell case

         where (k <= KMT(:,:,bid))
            DIAG_PE_2D(:,:,bid) = DIAG_PE_2D(:,:,bid) + &
                                  dzw(k-1)*WTK*grav*p5* &
                                  (RHOCUR(:,:,k-1) + RHOCUR(:,:,k))
         endwhere

      endif

   endif

   call advt(k,WORKN,WTK,TMIX,TCUR,UCUR,VCUR,this_block)

   FT = FT - WORKN   ! advt returns WORKN = +L(T) 

   if (ldiag_global) then
     if (partial_bottom_cells) then
       do n=1,POP_nt
         where (k <= KMT(:,:,bid)) DIAG_TRACER_ADV_2D(:,:,n,bid) = & 
                                   DIAG_TRACER_ADV_2D(:,:,n,bid) - &
                                   WORKN(:,:,n)*DZT(:,:,k,bid)
       end do
     else
       do n=1,POP_nt

         where (k <= KMT(:,:,bid)) DIAG_TRACER_ADV_2D(:,:,n,bid) = &
                                   DIAG_TRACER_ADV_2D(:,:,n,bid) - &
                                   WORKN(:,:,n)*dz(k)
       end do
     endif 
   endif

!-----------------------------------------------------------------------
!
!  vertical diffusion VDiff(T)
!
!-----------------------------------------------------------------------

   call vdifft(k, WORKN, TOLD, STF_IN, this_block)

   FT = FT + WORKN

   if (ldiag_global) then
     if (partial_bottom_cells) then
       do n=1,POP_nt
         where (k <= KMT(:,:,bid)) DIAG_TRACER_VDIFF_2D(:,:,n,bid) = &
                                   DIAG_TRACER_VDIFF_2D(:,:,n,bid) + &
                                   WORKN(:,:,n)*DZT(:,:,k,bid)
       end do
     else
       do n=1,POP_nt
         where (k <= KMT(:,:,bid)) DIAG_TRACER_VDIFF_2D(:,:,n,bid) = &
                                   DIAG_TRACER_VDIFF_2D(:,:,n,bid) + &
                                   WORKN(:,:,n)*dz(k)
       end do
     endif
   endif

!-----------------------------------------------------------------------
!
!  add tracer change in surface layer due to freshwater flux
!  if using variable thickness surface layer
!
!-----------------------------------------------------------------------

   if (k == 1 .and. sfc_layer_type == sfc_layer_varthick) then
      do n = 1,POP_nt
         FT(:,:,n) = FT(:,:,n) + dzr(1)*TFW_IN(:,:,n)
      enddo

      if (ldiag_global) then
         do n = 1,POP_nt
            DIAG_TRACER_SFC_FLX(:,:,n,bid) = TFW_IN(:,:,n)
         enddo
      endif
   endif

!-----------------------------------------------------------------------
!
!  add source terms
!
!-----------------------------------------------------------------------

   WORKN = c0

   call set_pt_interior(k,this_block,WORKN(:,:,1))
   call set_s_interior (k,this_block,WORKN(:,:,2))


   if (POP_nt > 2) call set_interior_passive_tracers(k, this_block, WORKN)
   

!-----------------------------------------------------------------------
!
!  add source terms from KPP and from shortwave solar absorption
!    if necessary.
!  NOTE:  this is here instead of in set_{pt,s}_interior in case
!    KPP and/or shortwave solar absorption are turned on but 
!    bulk restoring is not.
!
!-----------------------------------------------------------------------

   !*** does nothing if kpp not chosen - otherwise adds kpp sources
   call add_kpp_sources(WORKN, k, this_block)

   !*** if sw flux available, add penetrative shortwave
   call add_sw_absorb(WORKN, SHF_QSW(:,:,bid), k, this_block)


   FT = FT + WORKN


   if (ldiag_global) then
     if (partial_bottom_cells) then
       do n=1,POP_nt
         where (k <= KMT(:,:,bid)) DIAG_TRACER_SOURCE_2D(:,:,n,bid) = &
                                   DIAG_TRACER_SOURCE_2D(:,:,n,bid) + &
                                   WORKN(:,:,n)*DZT(:,:,k,bid)
       end do
     else
       do n=1,POP_nt
         where (k <= KMT(:,:,bid)) DIAG_TRACER_SOURCE_2D(:,:,n,bid) = &
                                   DIAG_TRACER_SOURCE_2D(:,:,n,bid) + &
                                   WORKN(:,:,n)*dz(k)
       end do
     endif
   endif


!-----------------------------------------------------------------------
!
!  save the explicit part of the RHS in TRACER(newtime)
!  if there is implicit vertical mixing
!
!  with pressure averaging and variable thickness surface layer, 
!  the RHS contains the surface height contribution for the 
!  predictor step (for T,S at k=1 only)
!
!-----------------------------------------------------------------------

   if (implicit_vertical_mix) then

      if (sfc_layer_type == sfc_layer_varthick .and. k == 1 .and. &
          lpressure_avg .and. leapfrogts) then

         do n = 1,2
            where (KMT(:,:,bid) > 0)  ! RHS for predictor
               TNEW(:,:,1,n) = c2dtt(1)*FT(:,:,n) - c2*TCUR(:,:,1,n)* &
                               (PCUR - POLD)/(grav*dz(1))
            endwhere
         enddo

         do n = 3,POP_nt

            TNEW(:,:,k,n) = merge(c2dtt(k)*FT(:,:,n), &
                                  c0, k <= KMT(:,:,bid))
         enddo

      else

         do n = 1,POP_nt
            TNEW(:,:,k,n) = merge(c2dtt(k)*FT(:,:,n), &
                                  c0, k <= KMT(:,:,bid))
         enddo

      endif

!-----------------------------------------------------------------------
!
!  for variable thickness surface layer, update all but surface
!    layers. 
!  at the surface:
!    if explicit vertical mixing and pressure averaging:
!      predict surface T,S and store RHS for all other tracers
!    otherwise
!      store RHS for all tracers for later update with new Psurf
!
!  if not a variable thickness surface layer, update tracers here
!     
!-----------------------------------------------------------------------

   else ! no implicit_vertical_mix

      if (sfc_layer_type == sfc_layer_varthick .and. k == 1) then
         if (lpressure_avg .and. leapfrogts) then

            !*** predict surface T and S with pressure avg

            do n = 1,2
               where (KMT(:,:,bid) > 0)
                  TNEW(:,:,1,n) = TOLD(:,:,1,n)                        &
                           + (c1/(c1 + PCUR/(grav*dz(1))))*            &
                             (c2dtt(1)*FT(:,:,n) - c2*TCUR(:,:,1,n)*   &
                             (PCUR-POLD)/(grav*dz(1)))
               elsewhere
                  TNEW(:,:,1,n) = c0  ! zero tracers on land pts
               endwhere
            end do

            !*** store RHS for other tracer surface layers

            do n = 3,POP_nt
               TNEW(:,:,k,n) = merge(c2dtt(k)*FT(:,:,n), &
                                     c0, k <= KMT(:,:,bid))
            enddo

         else

            !*** store RHS for all tracer surface layers

            do n = 1,POP_nt
               TNEW(:,:,k,n) = merge(c2dtt(k)*FT(:,:,n), &
                                     c0, k <= KMT(:,:,bid))
            enddo

         endif

      else   !*** update all tracers to new time

         do n = 1,POP_nt
            where (k <= KMT(:,:,bid))
               TNEW(:,:,k,n) = TOLD(:,:,k,n) + c2dtt(k)*FT(:,:,n)
            elsewhere
               TNEW(:,:,k,n) = c0  ! zero tracers at land points
            endwhere
         enddo

      endif
   endif ! implicit_vertical_mix


!-----------------------------------------------------------------------
!EOC

 end subroutine tracer_update

!***********************************************************************

 end module baroclinic

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
