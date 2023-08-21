MODULE dynhpg

   !!  v2:  takes into account the stretching of the vertical grid in hpg_ffr
   !!  v3:  takes into account the stretching of the vertical grid in the calculation of the reference profile
   
   !!======================================================================
   !!                       ***  MODULE  dynhpg  ***
   !! Ocean dynamics:  hydrostatic pressure gradient trend
   !!======================================================================
   !! History :  OPA  !  1987-09  (P. Andrich, M.-A. Foujols)  hpg_zco: Original code
   !!            5.0  !  1991-11  (G. Madec)
   !!            7.0  !  1996-01  (G. Madec)  hpg_sco: Original code for s-coordinates
   !!            8.0  !  1997-05  (G. Madec)  split dynber into dynkeg and dynhpg
   !!            8.5  !  2002-07  (G. Madec)  F90: Free form and module
   !!            8.5  !  2002-08  (A. Bozec)  hpg_zps: Original code
   !!   NEMO     1.0  !  2005-10  (A. Beckmann, B.W. An)  various s-coordinate options
   !!                 !         Original code for hpg_ctl, hpg_hel hpg_wdj, hpg_djc, hpg_rot
   !!             -   !  2005-11  (G. Madec) style & small optimisation
   !!            3.3  !  2010-10  (C. Ethe, G. Madec) reorganisation of initialisation phase
   !!            3.4  !  2011-11  (H. Liu) hpg_prj: Original code for s-coordinates
   !!                 !           (A. Coward) suppression of hel, wdj and rot options
   !!            3.6  !  2014-11  (P. Mathiot) hpg_isf: original code for ice shelf cavity
   !!            4.2  !  2020-12  (M. Bell, A. Young) hpg_djc: revised djc scheme
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dyn_hpg      : update the momentum trend with the now horizontal
   !!                  gradient of the hydrostatic pressure
   !!   dyn_hpg_init : initialisation and control of options
   !!       hpg_zco  : z-coordinate scheme
   !!       hpg_zps  : z-coordinate plus partial steps (interpolation)
   !!       hpg_sco  : s-coordinate (standard jacobian formulation)
   !!       hpg_isf  : s-coordinate (sco formulation) adapted to ice shelf
   !!       hpg_djc  : s-coordinate (Density Jacobian with constrained cubic splines (ccs))
   !!       hpg_prj  : s-coordinate (Pressure Jacobian with Cubic polynomial)
   !!       hpg_djr  : s-coordinate (Density Jacobian with ccs subtracting a reference)
   !!       hpg_ffr  : s-coordinate (Forces on faces subtracting a reference profile)
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE isf_oce , ONLY : risfload  ! ice shelf  (risfload variable)
   USE isfload , ONLY : isf_load  ! ice shelf  (isf_load routine )
   USE sbc_oce         ! surface variable (only for the flag with ice shelf)
   USE dom_oce         ! ocean space and time domain
   USE wet_dry         ! wetting and drying
   USE phycst          ! physical constants
   USE trd_oce         ! trends: ocean variables
   USE trddyn          ! trend manager: dynamics
   USE zpshde          ! partial step: hor. derivative     (zps_hde routine)
   !
   USE in_out_manager  ! I/O manager
   USE prtctl          ! Print control
   USE lbclnk          ! lateral boundary condition 
   USE lib_mpp         ! MPP library
   USE eosbn2          ! compute density
   USE timing          ! Timing
   USE iom

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dyn_hpg        ! routine called by step module
   PUBLIC   dyn_hpg_init   ! routine called by opa module

   !                                !!* Namelist namdyn_hpg : hydrostatic pressure gradient
   LOGICAL, PUBLIC ::   ln_hpg_zco   !: z-coordinate - full steps
   LOGICAL, PUBLIC ::   ln_hpg_zps   !: z-coordinate - partial steps (interpolation)
   LOGICAL, PUBLIC ::   ln_hpg_sco   !: s-coordinate (standard jacobian formulation)
   LOGICAL, PUBLIC ::   ln_hpg_djc   !: s-coordinate (Density Jacobian with Cubic polynomial)
   LOGICAL, PUBLIC ::   ln_hpg_prj   !: s-coordinate (Pressure Jacobian scheme)
   LOGICAL, PUBLIC ::   ln_hpg_isf   !: s-coordinate similar to sco modify for isf
   LOGICAL, PUBLIC ::   ln_hpg_djr   !: s-coordinate (density Jacobian with cubic polynomial, 
                                     !:               subtracting a local reference profile)
   LOGICAL, PUBLIC ::   ln_hpg_ffr   !: s-coordinate (forces on faces with subtraction of a 
                                     !:               local reference profile)

   !                                !! Flag to control the type of hydrostatic pressure gradient
   INTEGER, PARAMETER ::   np_ERROR  =-10   ! error in specification of lateral diffusion
   INTEGER, PARAMETER ::   np_zco    =  0   ! z-coordinate - full steps
   INTEGER, PARAMETER ::   np_zps    =  1   ! z-coordinate - partial steps (interpolation)
   INTEGER, PARAMETER ::   np_sco    =  2   ! s-coordinate (standard jacobian formulation)
   INTEGER, PARAMETER ::   np_djc    =  3   ! s-coordinate (Density Jacobian with Cubic polynomial)
   INTEGER, PARAMETER ::   np_prj    =  4   ! s-coordinate (Pressure Jacobian scheme)
   INTEGER, PARAMETER ::   np_isf    =  5   ! s-coordinate similar to sco modify for isf
   INTEGER, PARAMETER ::   np_djr    =  6   ! s-coordinate (density Jacobian with cubic polynomial, 
                                            !:              subtracting a local reference profile) 
   INTEGER, PARAMETER ::   np_ffr    =  7   ! s-coordinate (forces on faces with subtraction of a 
                                            !:              local reference profile) 

   !
   INTEGER, PUBLIC  ::   nhpg         !: type of pressure gradient scheme used 
                                      !: (deduced from ln_hpg_... flags) (PUBLIC for TAM)
   !

   LOGICAL     ::   ln_hpg_bcvN_rhd_hor, ln_hpg_bcvN_rhd_srf   !: flags to specify constrained cubic 
                                                               !: spline (ccs) bdy conditions for rhd & z
   LOGICAL     ::   ln_hpg_bcvN_rhd_bot, ln_hpg_bcvN_z_hor     !: True implies von Neumann bcs; 
                                                               !: False implies linear extrapolation

   REAL(wp)    ::   aco_bc_rhd_hor, bco_bc_rhd_hor, aco_bc_rhd_srf  !: coefficients for hpg_djc & hpg_djr boundary conditions 
   REAL(wp)    ::   bco_bc_rhd_srf, aco_bc_rhd_bot, bco_bc_rhd_bot  !:    "
   REAL(wp)    ::   aco_bc_z_hor, bco_bc_z_hor, aco_bc_z_srf        !:    " 
   REAL(wp)    ::   bco_bc_z_srf, aco_bc_z_bot, bco_bc_z_bot        !:    "

   LOGICAL     ::   ln_hpg_ref      ! T => reference profile is subtracted; 
                                    ! F => no reference profile subtracted
   LOGICAL     ::   ln_hpg_ref_str  ! T => reference profile calculated taking into
                                    !      account vertical variation in grid spacing
   LOGICAL     ::   ln_hpg_ref_ccs  ! T => constrained cubic, 
                                    ! F => simple cubic used for interpolation of reference 
                                    !      to target profiles   
   LOGICAL     ::   ln_hpg_ref_off  ! only applies if ln_hpg_ref_ccs = T; 
                                    ! T => use off-centred simple cubic at upper & lower boundaries 
   INTEGER     ::   nn_loc_ref_tgt  ! 1 = original Mike's algorithm; 
                                    ! 2 = alternative Diego's algorithm;

   LOGICAL     ::   ln_hpg_ffr_hor_ccs  ! T => use constrained cubic spline to 
                                        !      interpolate z_rhd_pmr along upper faces of cells
   LOGICAL     ::   ln_hpg_ffr_hor_cub  ! T => use simple cubic polynomial to interpolate 
                                        !      z_rhd_pmr along upper faces of cells
   LOGICAL     ::   ln_hpg_ffr_vrt_quad ! T => use quadratic fit to z_rhd_pmr in vertical interpoln & integration; 
                                        !       else standard 2nd order scheme

   LOGICAL     ::   ln_dbg_hpg  ! T => debug outputs generated
   LOGICAL     ::   ln_dbg_ij,  ln_dbg_ik,  ln_dbg_jk
   INTEGER     ::   ki_dbg_min, ki_dbg_max, ki_dbg_cen
   INTEGER     ::   kj_dbg_min, kj_dbg_max, kj_dbg_cen
   INTEGER     ::   kk_dbg_min, kk_dbg_max, kk_dbg_cen 


   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "domzgr_substitute.h90"

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: dynhpg.F90 16316 2023-08-08 09:45:11Z diegobruciaferri $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dyn_hpg( kt, Kmm, puu, pvv, Krhs )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_hpg  ***
      !!
      !! ** Method  :   Call the hydrostatic pressure gradient routine
      !!              using the scheme defined in the namelist
      !!
      !! ** Action : - Update (puu(:,:,:,Krhs),pvv(:,:,:,Krhs)) with the now hydrastatic pressure trend
      !!             - send trends to trd_dyn for futher diagnostics (l_trddyn=T)
      !!----------------------------------------------------------------------
      INTEGER                             , INTENT( in )  ::  kt          ! ocean time-step index
      INTEGER                             , INTENT( in )  ::  Kmm, Krhs   ! ocean time level indices
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpt), INTENT(inout) ::  puu, pvv    ! ocean velocities and RHS of momentum equation
      !
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   ztrdu, ztrdv
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('dyn_hpg')
      !
      IF( l_trddyn ) THEN                    ! Temporary saving of puu(:,:,:,Krhs) and pvv(:,:,:,Krhs) trends (l_trddyn)
         ALLOCATE( ztrdu(jpi,jpj,jpk) , ztrdv(jpi,jpj,jpk) )
         ztrdu(:,:,:) = puu(:,:,:,Krhs)
         ztrdv(:,:,:) = pvv(:,:,:,Krhs)
      ENDIF
      !
      SELECT CASE ( nhpg )      ! Hydrostatic pressure gradient computation
      CASE ( np_zco )   ;   CALL hpg_zco    ( kt, Kmm, puu, pvv, Krhs )  ! z-coordinate
      CASE ( np_zps )   ;   CALL hpg_zps    ( kt, Kmm, puu, pvv, Krhs )  ! z-coordinate plus partial steps (interpolation)
      CASE ( np_sco )   ;   CALL hpg_sco    ( kt, Kmm, puu, pvv, Krhs )  ! s-coordinate (standard jacobian formulation)
      CASE ( np_djc )   ;   CALL hpg_djc    ( kt, Kmm, puu, pvv, Krhs )  ! s-coordinate (Density Jacobian with Cubic polynomial)
      CASE ( np_prj )   ;   CALL hpg_prj    ( kt, Kmm, puu, pvv, Krhs )  ! s-coordinate (Pressure Jacobian scheme)
      CASE ( np_isf )   ;   CALL hpg_isf    ( kt, Kmm, puu, pvv, Krhs )  ! s-coordinate similar to sco modify for ice shelf
      CASE ( np_djr )   ;   CALL hpg_djr    ( kt, Kmm, puu, pvv, Krhs )  ! s-coordinate (Density Jacobian with Cubic polynomial 
                                                                         !               subtracting a local reference profile)
      CASE ( np_ffr )   ;   CALL hpg_ffr    ( kt, Kmm, puu, pvv, Krhs )  ! s-coordinate (forces on faces with subtraction of a 
                                                                         !               local reference profile)
      END SELECT
      !
      IF( l_trddyn ) THEN      ! save the hydrostatic pressure gradient trends for momentum trend diagnostics
         ztrdu(:,:,:) = puu(:,:,:,Krhs) - ztrdu(:,:,:)
         ztrdv(:,:,:) = pvv(:,:,:,Krhs) - ztrdv(:,:,:)
         CALL trd_dyn( ztrdu, ztrdv, jpdyn_hpg, kt, Kmm )
         DEALLOCATE( ztrdu , ztrdv )
      ENDIF
      !
      IF(sn_cfctl%l_prtctl)   CALL prt_ctl( tab3d_1=puu(:,:,:,Krhs), clinfo1=' hpg  - Ua: ', mask1=umask,   &
         &                                  tab3d_2=pvv(:,:,:,Krhs), clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
      !
      IF( ln_timing )   CALL timing_stop('dyn_hpg')
      !
   END SUBROUTINE dyn_hpg


   SUBROUTINE dyn_hpg_init( Kmm )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE dyn_hpg_init  ***
      !!
      !! ** Purpose :   initializations for the hydrostatic pressure gradient
      !!              computation and consistency control
      !!
      !! ** Action  :   Read the namelist namdyn_hpg and check the consistency
      !!      with the type of vertical coordinate used (zco, zps, sco)
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) :: Kmm   ! ocean time level index
      !
      INTEGER ::   ioptio = 0      ! temporary integer
      INTEGER ::   ios             ! Local integer output status for namelist read
      !!
      INTEGER  ::   ji, jj, jk, ikt    ! dummy loop indices      ISF
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::  zts_top, zrhd   ! hypothesys on isf density
      REAL(wp), ALLOCATABLE, DIMENSION(:,:)   ::  zrhdtop_isf    ! density at bottom of ISF
      REAL(wp), ALLOCATABLE, DIMENSION(:,:)   ::  ziceload       ! density at bottom of ISF
      !!

      NAMELIST/namdyn_hpg/ ln_hpg_zco, ln_hpg_zps, ln_hpg_sco, ln_hpg_isf,   &
         &                 ln_hpg_djc, ln_hpg_prj, ln_hpg_djr, ln_hpg_ffr,   &
         &                 ln_hpg_bcvN_rhd_hor,    ln_hpg_bcvN_rhd_srf,      & 
         &                 ln_hpg_bcvN_rhd_bot,    ln_hpg_bcvN_z_hor,        &                 
         &                 ln_hpg_ref,             nn_loc_ref_tgt,           &
         &                 ln_hpg_ref_str,         ln_hpg_ref_ccs,           & 
         &                 ln_hpg_ref_off,                                   &
         &                 ln_hpg_ffr_hor_ccs,     ln_hpg_ffr_hor_cub,       &
         &                 ln_hpg_ffr_vrt_quad,                              &   
         &                 ln_dbg_hpg, ln_dbg_ij,  ln_dbg_ik,  ln_dbg_jk,    &
         &                 ki_dbg_min, ki_dbg_max, ki_dbg_cen,               &  
         &                 kj_dbg_min, kj_dbg_max, kj_dbg_cen,               &
         &                 kk_dbg_min, kk_dbg_max, kk_dbg_cen           

      !!----------------------------------------------------------------------
      !
      READ  ( numnam_ref, namdyn_hpg, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namdyn_hpg in reference namelist' )
      !
      READ  ( numnam_cfg, namdyn_hpg, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namdyn_hpg in configuration namelist' )
      IF(lwm) WRITE ( numond, namdyn_hpg )
      !
      IF(lwp) THEN                   ! Control print
         WRITE(numout,*)
         WRITE(numout,*) 'dyn_hpg_init : hydrostatic pressure gradient initialisation'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namdyn_hpg : choice of hpg scheme'
         WRITE(numout,*) '      z-coord. - full steps                              ln_hpg_zco    = ', ln_hpg_zco
         WRITE(numout,*) '      z-coord. - partial steps (interpolation)           ln_hpg_zps    = ', ln_hpg_zps
         WRITE(numout,*) '      s-coord. (standard jacobian formulation)           ln_hpg_sco    = ', ln_hpg_sco
         WRITE(numout,*) '      s-coord. (standard jacobian formulation) for isf   ln_hpg_isf    = ', ln_hpg_isf
         WRITE(numout,*) '      s-coord. (Density Jacobian: Cubic polynomial)      ln_hpg_djc    = ', ln_hpg_djc
         WRITE(numout,*) '      s-coord. (Pressure Jacobian: Cubic polynomial)     ln_hpg_prj    = ', ln_hpg_prj
         WRITE(numout,*) '      s-coord. (Density Jacobian: Cubic minus reference) ln_hpg_djr    = ', ln_hpg_djr
         WRITE(numout,*) '      s-coord. (Density Jacobian: Cubic minus reference) ln_hpg_ffr    = ', ln_hpg_djr
         WRITE(numout,*) '      s-coord. (forces on faces minus local reference)   ln_hpg_ffr    = ', ln_hpg_ffr
      ENDIF
      !
      IF( .NOT.ln_linssh .AND. (ln_hpg_zco.OR.ln_hpg_zps) )   &
         &   CALL ctl_stop( 'dyn_hpg_init : non-linear free surface incompatible with hpg_zco or hpg_zps' )
      !
      IF( (.NOT.ln_hpg_isf .AND. ln_isfcav) .OR. (ln_hpg_isf .AND. .NOT.ln_isfcav) )                  &
         &   CALL ctl_stop( 'dyn_hpg_init : ln_hpg_isf=T requires ln_isfcav=T and vice versa' )  

      !
#if defined key_qco
      IF( ln_hpg_isf ) THEN
         CALL ctl_stop( 'dyn_hpg_init : key_qco and ln_hpg_isf not yet compatible' )
      ENDIF
#endif

      IF( ln_hpg_djr .AND. nn_hls < 2) THEN
         CALL ctl_stop( 'dyn_hpg_init : nn_hls < 2 and ln_hpg_djr not yet compatible' )
      ENDIF

      !
      !                               ! Set nhpg from ln_hpg_... flags & consistency check
      nhpg   = np_ERROR
      ioptio = 0
      IF( ln_hpg_zco ) THEN   ;   nhpg = np_zco   ;   ioptio = ioptio +1   ;   ENDIF
      IF( ln_hpg_zps ) THEN   ;   nhpg = np_zps   ;   ioptio = ioptio +1   ;   ENDIF
      IF( ln_hpg_sco ) THEN   ;   nhpg = np_sco   ;   ioptio = ioptio +1   ;   ENDIF
      IF( ln_hpg_djc ) THEN   ;   nhpg = np_djc   ;   ioptio = ioptio +1   ;   ENDIF
      IF( ln_hpg_prj ) THEN   ;   nhpg = np_prj   ;   ioptio = ioptio +1   ;   ENDIF
      IF( ln_hpg_isf ) THEN   ;   nhpg = np_isf   ;   ioptio = ioptio +1   ;   ENDIF
      IF( ln_hpg_djr ) THEN   ;   nhpg = np_djr   ;   ioptio = ioptio +1   ;   ENDIF
      IF( ln_hpg_ffr ) THEN   ;   nhpg = np_ffr   ;   ioptio = ioptio +1   ;   ENDIF
      !
      IF( ioptio /= 1 )   CALL ctl_stop( 'NO or several hydrostatic pressure gradient options used' )
      ! 
      IF(lwp) THEN
         WRITE(numout,*)
         SELECT CASE( nhpg )
         CASE( np_zco )   ;   WRITE(numout,*) '   ==>>>   z-coord. - full steps '
         CASE( np_zps )   ;   WRITE(numout,*) '   ==>>>   z-coord. - partial steps (interpolation)'
         CASE( np_sco )   ;   WRITE(numout,*) '   ==>>>   s-coord. (standard jacobian formulation)'
         CASE( np_djc )   ;   WRITE(numout,*) '   ==>>>   s-coord. (Density Jacobian: Cubic polynomial)'
         CASE( np_prj )   ;   WRITE(numout,*) '   ==>>>   s-coord. (Pressure Jacobian: Cubic polynomial)'
         CASE( np_isf )   ;   WRITE(numout,*) '   ==>>>   s-coord. (standard jacobian formulation) for isf'
         CASE( np_djr )   ;   WRITE(numout,*) '   ==>>>   s-coord. (Density Jacobian: Cubic polynomial subtracting a local reference profile)'
         CASE( np_ffr )   ;   WRITE(numout,*) '   ==>>>   s-coord. (forces on faces subtracting a local reference profile)'

         END SELECT
         WRITE(numout,*)
      ENDIF
      !                          
      IF ( ln_hpg_djc .OR. ln_hpg_djr .OR. ln_hpg_ffr ) THEN
         IF (ln_hpg_bcvN_rhd_hor) THEN ! Von Neumann boundary condition
           IF(lwp) WRITE(numout,*) '           ccs rhd horizontal bc: von Neumann '
           aco_bc_rhd_hor = 6.0_wp/5.0_wp
           bco_bc_rhd_hor = 7.0_wp/15.0_wp
         ELSE ! Linear extrapolation
           IF(lwp) WRITE(numout,*) '           ccs rhd horizontal bc: linear extrapolation'
           aco_bc_rhd_hor = 3.0_wp/2.0_wp
           bco_bc_rhd_hor = 1.0_wp/2.0_wp
         END IF
         IF (ln_hpg_bcvN_rhd_srf) THEN ! Von Neumann boundary condition
           IF(lwp) WRITE(numout,*) '           ccs rhd surface bc: von Neumann '
           aco_bc_rhd_srf = 6.0_wp/5.0_wp
           bco_bc_rhd_srf = 7.0_wp/15.0_wp
         ELSE ! Linear extrapolation
           IF(lwp) WRITE(numout,*) '           ccs rhd surface bc: linear extrapolation'
           aco_bc_rhd_srf = 3.0_wp/2.0_wp
           bco_bc_rhd_srf = 1.0_wp/2.0_wp
         END IF
         IF (ln_hpg_bcvN_rhd_bot) THEN ! Von Neumann boundary condition
           IF(lwp) WRITE(numout,*) '           ccs rhd bottom bc: von Neumann '
           aco_bc_rhd_bot = 6.0_wp/5.0_wp
           bco_bc_rhd_bot = 7.0_wp/15.0_wp
         ELSE ! Linear extrapolation
           IF(lwp) WRITE(numout,*) '           ccs rhd bottom bc: linear extrapolation'
           aco_bc_rhd_bot = 3.0_wp/2.0_wp
           bco_bc_rhd_bot = 1.0_wp/2.0_wp
         END IF
         IF (ln_hpg_bcvN_z_hor) THEN ! Von Neumann boundary condition
           IF(lwp) WRITE(numout,*) '           ccs z horizontal bc: von Neumann '
           aco_bc_z_hor = 6.0_wp/5.0_wp
           bco_bc_z_hor = 7.0_wp/15.0_wp
         ELSE ! Linear extrapolation
           IF(lwp) WRITE(numout,*) '           ccs z horizontal bc: linear extrapolation'
           aco_bc_z_hor = 3.0_wp/2.0_wp
           bco_bc_z_hor = 1.0_wp/2.0_wp
         END IF

! linear extrapolation used for z in the vertical (surface and bottom)   
         aco_bc_z_srf = 3.0_wp/2.0_wp
         bco_bc_z_srf = 1.0_wp/2.0_wp
         aco_bc_z_bot = 3.0_wp/2.0_wp
         bco_bc_z_bot = 1.0_wp/2.0_wp

      END IF

      IF ( lwp .AND. ln_hpg_ref ) THEN
         IF ( nn_loc_ref_tgt == 1) THEN
            WRITE(numout,*) '       original Mikes algorithm to identify depths for interpolation of reference profile'
         ELSE IF ( nn_loc_ref_tgt == 2) THEN
            WRITE(numout,*) '       alternative Diegos algorithm to identify depths for interpolation of reference profile' 
         ELSE
            WRITE(numout,*) '       ERROR: the chosen algorithm to identify depths for interpolation of reference profile DOES NOT EXIST'
         END IF
         IF ( ln_hpg_ref_ccs .AND. ln_hpg_ref_off .AND. ln_hpg_ref_str ) THEN 
	    WRITE(numout,*) '           interpolation of reference profile by constrained cubic spline with off-centring at boundaries'                  
         ELSE IF ( ln_hpg_ref_ccs ) THEN
	    WRITE(numout,*) '           interpolation of reference profile by constrained cubic spline using boundary conditions'
	 ELSE 
	    WRITE(numout,*) '           interpolation of reference profile by simple cubic spline with off-centring at boundaries'  
         END IF 
         IF ( ln_hpg_ref_str ) THEN 
	    WRITE(numout,*) '           interpolation of reference profile takes account of variation in vertical grid spacing'                  
	 ELSE 
	    WRITE(numout,*) '           interpolation of reference profile does NOT take account of variation in vertical grid spacing' 
         END IF 
      END IF       

      IF ( lwp .AND. ln_hpg_ffr ) THEN 
         IF ( .NOT. ln_hpg_ref  ) THEN
            WRITE(numout,*) '           no reference profile subtracted '                   
         END IF       
         IF ( ln_hpg_ffr_hor_ccs ) THEN 
	    WRITE(numout,*) '           interpolation of z_rhd_pmr profile in horizontal by constrained cubic spline '    
         ELSE if ( ln_hpg_ffr_hor_cub ) THEN 
	    WRITE(numout,*) '           interpolation of z_rhd_pmr profile in horizontal by simple cubic polynomial '
	 ELSE 
	    WRITE(numout,*) '           simple linear interpolation in horizontal of z_rhd_pmr profile'       
         END IF       
         IF ( ln_hpg_ffr_vrt_quad ) THEN 
	    WRITE(numout,*) '           quadratic fit to z_rhd_pmr profile in vertical for interpolation and integration '    
         ELSE 
	    WRITE(numout,*) '           simple second order accurate vertical integration of z_rhd_pmr profile in vertical'       
         END IF       

      END IF       

      !
      IF ( ln_dbg_hpg .AND. lwp ) THEN 
         WRITE(numout,*) '             dyn_hpg diagnostic settings'  
         WRITE(numout,*) ' ki_dbg_min = ', ki_dbg_min, '; ki_dbg_max = ', ki_dbg_max
         WRITE(numout,*) ' kj_dbg_min = ', kj_dbg_min, '; kj_dbg_max = ', kj_dbg_max
         WRITE(numout,*) ' kk_dbg_min = ', kk_dbg_min, '; kk_dbg_max = ', kk_dbg_max
         WRITE(numout,*) ' ki_dbg_cen = ', ki_dbg_cen, '; kj_dbg_cen = ', kj_dbg_cen
         WRITE(numout,*) ' kk_dbg_cen = ', kk_dbg_cen 
      END IF 

   END SUBROUTINE dyn_hpg_init


   SUBROUTINE hpg_zco( kt, Kmm, puu, pvv, Krhs )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE hpg_zco  ***
      !!
      !! ** Method  :   z-coordinate case, levels are horizontal surfaces.
      !!      The now hydrostatic pressure gradient at a given level, jk,
      !!      is computed by taking the vertical integral of the in-situ
      !!      density gradient along the model level from the suface to that
      !!      level:    zhpi = grav .....
      !!                zhpj = grav .....
      !!      add it to the general momentum trend (puu(:,:,:,Krhs),pvv(:,:,:,Krhs)).
      !!            puu(:,:,:,Krhs) = puu(:,:,:,Krhs) - 1/e1u * zhpi
      !!            pvv(:,:,:,Krhs) = pvv(:,:,:,Krhs) - 1/e2v * zhpj
      !!
      !! ** Action : - Update (puu(:,:,:,Krhs),pvv(:,:,:,Krhs)) with the now hydrastatic pressure trend
      !!----------------------------------------------------------------------
      INTEGER                             , INTENT( in )  ::  kt          ! ocean time-step index
      INTEGER                             , INTENT( in )  ::  Kmm, Krhs   ! ocean time level indices
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpt), INTENT(inout) ::  puu, pvv    ! ocean velocities and RHS of momentum equation
      !
      INTEGER  ::   ji, jj, jk       ! dummy loop indices
      REAL(wp) ::   zcoef0, zcoef1   ! temporary scalars
      REAL(wp), DIMENSION(A2D(nn_hls)) ::  zhpi, zhpj
      !!----------------------------------------------------------------------
      !
      IF( .NOT. l_istiled .OR. ntile == 1 )  THEN                       ! Do only on the first tile
         IF( kt == nit000 ) THEN
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'dyn:hpg_zco : hydrostatic pressure gradient trend'
            IF(lwp) WRITE(numout,*) '~~~~~~~~~~~   z-coordinate case '
         ENDIF
      ENDIF
      !
      zcoef0 = - grav * 0.5_wp            ! Local constant initialization
      !
      DO_2D( 0, 0, 0, 0 )                 ! Surface value
         zcoef1 = zcoef0 * e3w(ji,jj,1,Kmm)
         !                                   ! hydrostatic pressure gradient
         zhpi(ji,jj) = zcoef1 * ( rhd(ji+1,jj,1) - rhd(ji,jj,1) ) * r1_e1u(ji,jj)
         zhpj(ji,jj) = zcoef1 * ( rhd(ji,jj+1,1) - rhd(ji,jj,1) ) * r1_e2v(ji,jj)
         !                                   ! add to the general momentum trend
         puu(ji,jj,1,Krhs) = puu(ji,jj,1,Krhs) + zhpi(ji,jj)
         pvv(ji,jj,1,Krhs) = pvv(ji,jj,1,Krhs) + zhpj(ji,jj)
      END_2D
      !
      DO_3D( 0, 0, 0, 0, 2, jpkm1 )        ! interior value (2=<jk=<jpkm1)
         zcoef1 = zcoef0 * e3w(ji,jj,jk,Kmm)
         !                                   ! hydrostatic pressure gradient
         zhpi(ji,jj) = zhpi(ji,jj) + zcoef1 * (  ( rhd(ji+1,jj,jk)+rhd(ji+1,jj,jk-1) )  &
            &                                  - ( rhd(ji  ,jj,jk)+rhd(ji  ,jj,jk-1) )  ) * r1_e1u(ji,jj)

         zhpj(ji,jj) = zhpj(ji,jj) + zcoef1 * (  ( rhd(ji,jj+1,jk)+rhd(ji,jj+1,jk-1) )  &
            &                                  - ( rhd(ji,jj,  jk)+rhd(ji,jj  ,jk-1) )  ) * r1_e2v(ji,jj)
         !                                   ! add to the general momentum trend
         puu(ji,jj,jk,Krhs) = puu(ji,jj,jk,Krhs) + zhpi(ji,jj)
         pvv(ji,jj,jk,Krhs) = pvv(ji,jj,jk,Krhs) + zhpj(ji,jj)
      END_3D
      !
   END SUBROUTINE hpg_zco


   SUBROUTINE hpg_zps( kt, Kmm, puu, pvv, Krhs )
      !!---------------------------------------------------------------------
      !!                 ***  ROUTINE hpg_zps  ***
      !!
      !! ** Method  :   z-coordinate plus partial steps case.  blahblah...
      !!
      !! ** Action  : - Update (puu(:,:,:,Krhs),pvv(:,:,:,Krhs)) with the now hydrastatic pressure trend
      !!----------------------------------------------------------------------
      INTEGER                             , INTENT( in )  ::  kt          ! ocean time-step index
      INTEGER                             , INTENT( in )  ::  Kmm, Krhs   ! ocean time level indices
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpt), INTENT(inout) ::  puu, pvv    ! ocean velocities and RHS of momentum equation
      !!
      INTEGER  ::   ji, jj, jk                       ! dummy loop indices
      INTEGER  ::   iku, ikv                         ! temporary integers
      REAL(wp) ::   zcoef0, zcoef1, zcoef2, zcoef3   ! temporary scalars
      REAL(wp), DIMENSION(A2D(nn_hls),jpk ) :: zhpi, zhpj
      REAL(wp), DIMENSION(A2D(nn_hls),jpts) :: zgtsu, zgtsv
      REAL(wp), DIMENSION(A2D(nn_hls)     ) :: zgru, zgrv
      !!----------------------------------------------------------------------
      !
      IF( .NOT. l_istiled .OR. ntile == 1 )  THEN                       ! Do only on the first tile
         IF( kt == nit000 ) THEN
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'dyn:hpg_zps : hydrostatic pressure gradient trend'
            IF(lwp) WRITE(numout,*) '~~~~~~~~~~~   z-coordinate with partial steps - vector optimization'
         ENDIF
      ENDIF

      ! Partial steps: Compute NOW horizontal gradient of t, s, rd at the last ocean level
      CALL zps_hde( kt, Kmm, jpts, ts(:,:,:,:,Kmm), zgtsu, zgtsv, rhd, zgru , zgrv )

      ! Local constant initialization
      zcoef0 = - grav * 0.5_wp

      !  Surface value (also valid in partial step case)
      DO_2D( 0, 0, 0, 0 )
         zcoef1 = zcoef0 * e3w(ji,jj,1,Kmm)
         ! hydrostatic pressure gradient
         zhpi(ji,jj,1) = zcoef1 * ( rhd(ji+1,jj  ,1) - rhd(ji,jj,1) ) * r1_e1u(ji,jj)
         zhpj(ji,jj,1) = zcoef1 * ( rhd(ji  ,jj+1,1) - rhd(ji,jj,1) ) * r1_e2v(ji,jj)
         ! add to the general momentum trend
         puu(ji,jj,1,Krhs) = puu(ji,jj,1,Krhs) + zhpi(ji,jj,1)
         pvv(ji,jj,1,Krhs) = pvv(ji,jj,1,Krhs) + zhpj(ji,jj,1)
      END_2D

      ! interior value (2=<jk=<jpkm1)
      DO_3D( 0, 0, 0, 0, 2, jpkm1 )
         zcoef1 = zcoef0 * e3w(ji,jj,jk,Kmm)
         ! hydrostatic pressure gradient
         zhpi(ji,jj,jk) = zhpi(ji,jj,jk-1)   &
            &           + zcoef1 * (  ( rhd(ji+1,jj,jk) + rhd(ji+1,jj,jk-1) )   &
            &                       - ( rhd(ji  ,jj,jk) + rhd(ji  ,jj,jk-1) )  ) * r1_e1u(ji,jj)

         zhpj(ji,jj,jk) = zhpj(ji,jj,jk-1)   &
            &           + zcoef1 * (  ( rhd(ji,jj+1,jk) + rhd(ji,jj+1,jk-1) )   &
            &                       - ( rhd(ji,jj,  jk) + rhd(ji,jj  ,jk-1) )  ) * r1_e2v(ji,jj)
         ! add to the general momentum trend
         puu(ji,jj,jk,Krhs) = puu(ji,jj,jk,Krhs) + zhpi(ji,jj,jk)
         pvv(ji,jj,jk,Krhs) = pvv(ji,jj,jk,Krhs) + zhpj(ji,jj,jk)
      END_3D

      ! partial steps correction at the last level  (use zgru & zgrv computed in zpshde.F90)
      DO_2D( 0, 0, 0, 0 )
         iku = mbku(ji,jj)
         ikv = mbkv(ji,jj)
         zcoef2 = zcoef0 * MIN( e3w(ji,jj,iku,Kmm), e3w(ji+1,jj  ,iku,Kmm) )
         zcoef3 = zcoef0 * MIN( e3w(ji,jj,ikv,Kmm), e3w(ji  ,jj+1,ikv,Kmm) )
         IF( iku > 1 ) THEN            ! on i-direction (level 2 or more)
            puu  (ji,jj,iku,Krhs) = puu(ji,jj,iku,Krhs) - zhpi(ji,jj,iku)         ! subtract old value
            zhpi(ji,jj,iku) = zhpi(ji,jj,iku-1)                   &   ! compute the new one
               &            + zcoef2 * ( rhd(ji+1,jj,iku-1) - rhd(ji,jj,iku-1) + zgru(ji,jj) ) * r1_e1u(ji,jj)
            puu  (ji,jj,iku,Krhs) = puu(ji,jj,iku,Krhs) + zhpi(ji,jj,iku)         ! add the new one to the general momentum trend
         ENDIF
         IF( ikv > 1 ) THEN            ! on j-direction (level 2 or more)
            pvv  (ji,jj,ikv,Krhs) = pvv(ji,jj,ikv,Krhs) - zhpj(ji,jj,ikv)         ! subtract old value
            zhpj(ji,jj,ikv) = zhpj(ji,jj,ikv-1)                   &   ! compute the new one
               &            + zcoef3 * ( rhd(ji,jj+1,ikv-1) - rhd(ji,jj,ikv-1) + zgrv(ji,jj) ) * r1_e2v(ji,jj)
            pvv  (ji,jj,ikv,Krhs) = pvv(ji,jj,ikv,Krhs) + zhpj(ji,jj,ikv)         ! add the new one to the general momentum trend
         ENDIF
      END_2D
      !
   END SUBROUTINE hpg_zps


   SUBROUTINE hpg_sco( kt, Kmm, puu, pvv, Krhs )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE hpg_sco  ***
      !!
      !! ** Method  :   s-coordinate case. Jacobian scheme.
      !!      The now hydrostatic pressure gradient at a given level, jk,
      !!      is computed by taking the vertical integral of the in-situ
      !!      density gradient along the model level from the suface to that
      !!      level. s-coordinates (ln_sco): a corrective term is added
      !!      to the horizontal pressure gradient :
      !!         zhpi = grav .....  + 1/e1u mi(rhd) di[ grav dep3w ]
      !!         zhpj = grav .....  + 1/e2v mj(rhd) dj[ grav dep3w ]
      !!      add it to the general momentum trend (puu(:,:,:,Krhs),pvv(:,:,:,Krhs)).
      !!         puu(:,:,:,Krhs) = puu(:,:,:,Krhs) - 1/e1u * zhpi
      !!         pvv(:,:,:,Krhs) = pvv(:,:,:,Krhs) - 1/e2v * zhpj
      !!
      !! ** Action : - Update (puu(:,:,:,Krhs),pvv(:,:,:,Krhs)) with the now hydrastatic pressure trend
      !!----------------------------------------------------------------------
      INTEGER                             , INTENT( in )  ::  kt          ! ocean time-step index
      INTEGER                             , INTENT( in )  ::  Kmm, Krhs   ! ocean time level indices
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpt), INTENT(inout) ::  puu, pvv    ! ocean velocities and RHS of momentum equation
      !!
      INTEGER  ::   ji, jj, jk, jii, jjj           ! dummy loop indices
      REAL(wp) ::   zcoef0, zuap, zvap, ztmp       ! local scalars
      LOGICAL  ::   ll_tmp1, ll_tmp2               ! local logical variables
      REAL(wp), DIMENSION(A2D(nn_hls),jpk)  ::   zhpi, zhpj
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   zcpx, zcpy   !W/D pressure filter
      !!----------------------------------------------------------------------
      !
      IF( ln_wd_il ) ALLOCATE(zcpx(A2D(nn_hls)), zcpy(A2D(nn_hls)))
      !
      IF( .NOT. l_istiled .OR. ntile == 1 )  THEN                       ! Do only on the first tile
         IF( kt == nit000 ) THEN
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'dyn:hpg_sco : hydrostatic pressure gradient trend'
            IF(lwp) WRITE(numout,*) '~~~~~~~~~~~   s-coordinate case, OCE original scheme used'
         ENDIF
      ENDIF
      !
      zcoef0 = - grav * 0.5_wp
      !
      IF( ln_wd_il ) THEN
        DO_2D( 0, 0, 0, 0 )
          ll_tmp1 = MIN(  ssh(ji,jj,Kmm)               ,  ssh(ji+1,jj,Kmm) ) >                &
               &    MAX( -ht_0(ji,jj)               , -ht_0(ji+1,jj) ) .AND.            &
               &    MAX(  ssh(ji,jj,Kmm) +  ht_0(ji,jj),  ssh(ji+1,jj,Kmm) + ht_0(ji+1,jj) )  &
               &                                                       > rn_wdmin1 + rn_wdmin2
          ll_tmp2 = ( ABS( ssh(ji,jj,Kmm)              -  ssh(ji+1,jj,Kmm) ) > 1.E-12 ) .AND. (       &
               &    MAX(   ssh(ji,jj,Kmm)              ,  ssh(ji+1,jj,Kmm) ) >                &
               &    MAX(  -ht_0(ji,jj)              , -ht_0(ji+1,jj) ) + rn_wdmin1 + rn_wdmin2 )

          IF(ll_tmp1) THEN
            zcpx(ji,jj) = 1.0_wp
          ELSE IF(ll_tmp2) THEN
            ! no worries about  ssh(ji+1,jj,Kmm) -  ssh(ji  ,jj,Kmm) = 0, it won't happen ! here
            zcpx(ji,jj) = ABS( (ssh(ji+1,jj,Kmm) + ht_0(ji+1,jj) - ssh(ji,jj,Kmm) - ht_0(ji,jj)) &
                        &    / (ssh(ji+1,jj,Kmm) - ssh(ji  ,jj,Kmm)) )
          ELSE
            zcpx(ji,jj) = 0._wp
          END IF
   
          ll_tmp1 = MIN(  ssh(ji,jj,Kmm)              ,  ssh(ji,jj+1,Kmm) ) >                &
               &    MAX( -ht_0(ji,jj)              , -ht_0(ji,jj+1) ) .AND.            &
               &    MAX(  ssh(ji,jj,Kmm) + ht_0(ji,jj),  ssh(ji,jj+1,Kmm) + ht_0(ji,jj+1) )  &
               &                                                      > rn_wdmin1 + rn_wdmin2
          ll_tmp2 = ( ABS( ssh(ji,jj,Kmm)             -  ssh(ji,jj+1,Kmm) ) > 1.E-12 ) .AND. (        &
               &    MAX(   ssh(ji,jj,Kmm)             ,  ssh(ji,jj+1,Kmm) ) >                &
               &    MAX(  -ht_0(ji,jj)             , -ht_0(ji,jj+1) ) + rn_wdmin1 + rn_wdmin2 )

          IF(ll_tmp1) THEN
            zcpy(ji,jj) = 1.0_wp
          ELSE IF(ll_tmp2) THEN
            ! no worries about  ssh(ji,jj+1,Kmm) -  ssh(ji,jj  ,Kmm) = 0, it won't happen ! here
            zcpy(ji,jj) = ABS( (ssh(ji,jj+1,Kmm) + ht_0(ji,jj+1) - ssh(ji,jj,Kmm) - ht_0(ji,jj)) &
                        &    / (ssh(ji,jj+1,Kmm) - ssh(ji,jj  ,Kmm)) )
          ELSE
            zcpy(ji,jj) = 0._wp
          END IF
        END_2D
      END IF
      !
      DO_2D( 0, 0, 0, 0 )              ! Surface value
         !                                   ! hydrostatic pressure gradient along s-surfaces
         zhpi(ji,jj,1) = zcoef0 * r1_e1u(ji,jj)                      &
            &          * (  e3w(ji+1,jj  ,1,Kmm) * rhd(ji+1,jj  ,1)  &
            &             - e3w(ji  ,jj  ,1,Kmm) * rhd(ji  ,jj  ,1)  )
         zhpj(ji,jj,1) = zcoef0 * r1_e2v(ji,jj)                      &
            &          * (  e3w(ji  ,jj+1,1,Kmm) * rhd(ji  ,jj+1,1)  &
            &             - e3w(ji  ,jj  ,1,Kmm) * rhd(ji  ,jj  ,1)  )
         !                                   ! s-coordinate pressure gradient correction
         zuap = -zcoef0 * ( rhd    (ji+1,jj,1) + rhd    (ji,jj,1) )   &
            &           * ( gde3w(ji+1,jj,1) - gde3w(ji,jj,1) ) * r1_e1u(ji,jj)
         zvap = -zcoef0 * ( rhd    (ji,jj+1,1) + rhd    (ji,jj,1) )   &
            &           * ( gde3w(ji,jj+1,1) - gde3w(ji,jj,1) ) * r1_e2v(ji,jj)
         !
         IF( ln_wd_il ) THEN
            zhpi(ji,jj,1) = zhpi(ji,jj,1) * zcpx(ji,jj)
            zhpj(ji,jj,1) = zhpj(ji,jj,1) * zcpy(ji,jj) 
            zuap = zuap * zcpx(ji,jj)
            zvap = zvap * zcpy(ji,jj)
         ENDIF
         !                                   ! add to the general momentum trend
         puu(ji,jj,1,Krhs) = puu(ji,jj,1,Krhs) + zhpi(ji,jj,1) + zuap
         pvv(ji,jj,1,Krhs) = pvv(ji,jj,1,Krhs) + zhpj(ji,jj,1) + zvap
      END_2D
      !
      DO_3D( 0, 0, 0, 0, 2, jpkm1 )    ! interior value (2=<jk=<jpkm1)
         !                                   ! hydrostatic pressure gradient along s-surfaces
         zhpi(ji,jj,jk) = zhpi(ji,jj,jk-1) + zcoef0 * r1_e1u(ji,jj)                         &
            &           * (  e3w(ji+1,jj,jk,Kmm) * ( rhd(ji+1,jj,jk) + rhd(ji+1,jj,jk-1) )  &
            &              - e3w(ji  ,jj,jk,Kmm) * ( rhd(ji  ,jj,jk) + rhd(ji  ,jj,jk-1) )  )
         zhpj(ji,jj,jk) = zhpj(ji,jj,jk-1) + zcoef0 * r1_e2v(ji,jj)                         &
            &           * (  e3w(ji,jj+1,jk,Kmm) * ( rhd(ji,jj+1,jk) + rhd(ji,jj+1,jk-1) )  &
            &              - e3w(ji,jj  ,jk,Kmm) * ( rhd(ji,jj,  jk) + rhd(ji,jj  ,jk-1) )  )
         !                                   ! s-coordinate pressure gradient correction
         zuap = -zcoef0 * ( rhd  (ji+1,jj  ,jk) + rhd  (ji,jj,jk) ) &
            &           * ( gde3w(ji+1,jj  ,jk) - gde3w(ji,jj,jk) ) * r1_e1u(ji,jj)
         zvap = -zcoef0 * ( rhd  (ji  ,jj+1,jk) + rhd  (ji,jj,jk) ) &
            &           * ( gde3w(ji  ,jj+1,jk) - gde3w(ji,jj,jk) ) * r1_e2v(ji,jj)
         !
         IF( ln_wd_il ) THEN
            zhpi(ji,jj,jk) = zhpi(ji,jj,jk) * zcpx(ji,jj)
            zhpj(ji,jj,jk) = zhpj(ji,jj,jk) * zcpy(ji,jj) 
            zuap = zuap * zcpx(ji,jj)
            zvap = zvap * zcpy(ji,jj)
         ENDIF
         !
         ! add to the general momentum trend
         puu(ji,jj,jk,Krhs) = puu(ji,jj,jk,Krhs) + zhpi(ji,jj,jk) + zuap
         pvv(ji,jj,jk,Krhs) = pvv(ji,jj,jk,Krhs) + zhpj(ji,jj,jk) + zvap
      END_3D
      !
      IF( ln_wd_il )  DEALLOCATE( zcpx , zcpy )
      !
   END SUBROUTINE hpg_sco


   SUBROUTINE hpg_isf( kt, Kmm, puu, pvv, Krhs )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE hpg_isf  ***
      !!
      !! ** Method  :   s-coordinate case. Jacobian scheme.
      !!      The now hydrostatic pressure gradient at a given level, jk,
      !!      is computed by taking the vertical integral of the in-situ
      !!      density gradient along the model level from the suface to that
      !!      level. s-coordinates (ln_sco): a corrective term is added
      !!      to the horizontal pressure gradient :
      !!         zhpi = grav .....  + 1/e1u mi(rhd) di[ grav dep3w ]
      !!         zhpj = grav .....  + 1/e2v mj(rhd) dj[ grav dep3w ]
      !!      add it to the general momentum trend (puu(:,:,:,Krhs),pvv(:,:,:,Krhs)).
      !!         puu(:,:,:,Krhs) = puu(:,:,:,Krhs) - 1/e1u * zhpi
      !!         pvv(:,:,:,Krhs) = pvv(:,:,:,Krhs) - 1/e2v * zhpj
      !!      iceload is added
      !!      
      !! ** Action : - Update (puu(:,:,:,Krhs),pvv(:,:,:,Krhs)) with the now hydrastatic pressure trend
      !!----------------------------------------------------------------------
      INTEGER                             , INTENT( in )  ::  kt          ! ocean time-step index
      INTEGER                             , INTENT( in )  ::  Kmm, Krhs   ! ocean time level indices
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpt), INTENT(inout) ::  puu, pvv    ! ocean velocities and RHS of momentum equation
      !!
      INTEGER  ::   ji, jj, jk             ! dummy loop indices
      INTEGER  ::   ikt ,  ikti1,  iktj1   ! local integer
      REAL(wp) ::   ze3w, ze3wi1, ze3wj1   ! local scalars
      REAL(wp) ::   zcoef0, zuap, zvap     !   -      -
      REAL(wp), DIMENSION(A2D(nn_hls),jpk ) ::  zhpi, zhpj
      REAL(wp), DIMENSION(A2D(nn_hls),jpts) ::  zts_top
      REAL(wp), DIMENSION(A2D(nn_hls))      ::  zrhdtop_oce
      !!----------------------------------------------------------------------
      !
      zcoef0 = - grav * 0.5_wp   ! Local constant initialization
      !
      !                          ! iniitialised to 0. zhpi zhpi 
      zhpi(:,:,:) = 0._wp   ;   zhpj(:,:,:) = 0._wp

      ! compute rhd at the ice/oce interface (ocean side)
      ! usefull to reduce residual current in the test case ISOMIP with no melting
      DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
         ikt = mikt(ji,jj)
         zts_top(ji,jj,1) = ts(ji,jj,ikt,1,Kmm)
         zts_top(ji,jj,2) = ts(ji,jj,ikt,2,Kmm)
      END_2D
      CALL eos( zts_top, risfdep, zrhdtop_oce )

      !                     !===========================!
      !                     !=====  surface value  =====!
      !                     !===========================!
      DO_2D( 0, 0, 0, 0 )
         ikt   = mikt(ji  ,jj  )   ;   ze3w   = e3w(ji  ,jj  ,ikt  ,Kmm)
         ikti1 = mikt(ji+1,jj  )   ;   ze3wi1 = e3w(ji+1,jj  ,ikti1,Kmm)
         iktj1 = mikt(ji  ,jj+1)   ;   ze3wj1 = e3w(ji  ,jj+1,iktj1,Kmm)
         !                          ! hydrostatic pressure gradient along s-surfaces and ice shelf pressure
         !                          ! we assume ISF is in isostatic equilibrium
         zhpi(ji,jj,1) = zcoef0 * r1_e1u(ji,jj) * (   risfload(ji+1,jj) - risfload(ji,jj)  &
            &                                       + 0.5_wp * ( ze3wi1 * ( rhd(ji+1,jj,ikti1) + zrhdtop_oce(ji+1,jj) )     &
            &                                                  - ze3w   * ( rhd(ji  ,jj,ikt  ) + zrhdtop_oce(ji  ,jj) ) )   )
         zhpj(ji,jj,1) = zcoef0 * r1_e2v(ji,jj) * (   risfload(ji,jj+1) - risfload(ji,jj)  &
            &                                       + 0.5_wp * ( ze3wj1 * ( rhd(ji,jj+1,iktj1) + zrhdtop_oce(ji,jj+1) )      &
            &                                                  - ze3w   * ( rhd(ji,jj  ,ikt  ) + zrhdtop_oce(ji,jj  ) ) )   ) 
         !                          ! s-coordinate pressure gradient correction (=0 if z coordinate)
         zuap = -zcoef0 * ( rhd    (ji+1,jj,1) + rhd    (ji,jj,1) )   &
            &           * ( gde3w(ji+1,jj,1) - gde3w(ji,jj,1) ) * r1_e1u(ji,jj)
         zvap = -zcoef0 * ( rhd    (ji,jj+1,1) + rhd    (ji,jj,1) )   &
            &           * ( gde3w(ji,jj+1,1) - gde3w(ji,jj,1) ) * r1_e2v(ji,jj)
         !                          ! add to the general momentum trend
         puu(ji,jj,1,Krhs) = puu(ji,jj,1,Krhs) + (zhpi(ji,jj,1) + zuap) * umask(ji,jj,1)
         pvv(ji,jj,1,Krhs) = pvv(ji,jj,1,Krhs) + (zhpj(ji,jj,1) + zvap) * vmask(ji,jj,1)
      END_2D
      !   
      !                     !=============================!
      !                     !=====  interior values  =====!
      !                     !=============================!
      DO_3D( 0, 0, 0, 0, 2, jpkm1 )
         ze3w   = e3w(ji  ,jj  ,jk,Kmm)
         ze3wi1 = e3w(ji+1,jj  ,jk,Kmm)
         ze3wj1 = e3w(ji  ,jj+1,jk,Kmm)
         !                          ! hydrostatic pressure gradient along s-surfaces
         zhpi(ji,jj,jk) = zhpi(ji,jj,jk-1) + zcoef0 / e1u(ji,jj)   &
            &           * (  ze3wi1 * ( rhd(ji+1,jj,jk) + rhd(ji+1,jj,jk-1) ) * wmask(ji+1,jj,jk)   &
            &              - ze3w   * ( rhd(ji  ,jj,jk) + rhd(ji  ,jj,jk-1) ) * wmask(ji  ,jj,jk)   )
         zhpj(ji,jj,jk) = zhpj(ji,jj,jk-1) + zcoef0 / e2v(ji,jj)   &
            &           * (  ze3wj1 * ( rhd(ji,jj+1,jk) + rhd(ji,jj+1,jk-1) ) * wmask(ji,jj+1,jk)   &
            &              - ze3w   * ( rhd(ji,jj,  jk) + rhd(ji,jj  ,jk-1) ) * wmask(ji,jj  ,jk)   )
         !                          ! s-coordinate pressure gradient correction
         zuap = -zcoef0 * ( rhd   (ji+1,jj  ,jk) + rhd   (ji,jj,jk) )   &
            &           * ( gde3w(ji+1,jj  ,jk) - gde3w(ji,jj,jk) ) / e1u(ji,jj)
         zvap = -zcoef0 * ( rhd   (ji  ,jj+1,jk) + rhd   (ji,jj,jk) )   &
            &           * ( gde3w(ji  ,jj+1,jk) - gde3w(ji,jj,jk) ) / e2v(ji,jj)
         !                          ! add to the general momentum trend
         puu(ji,jj,jk,Krhs) = puu(ji,jj,jk,Krhs) + (zhpi(ji,jj,jk) + zuap) * umask(ji,jj,jk)
         pvv(ji,jj,jk,Krhs) = pvv(ji,jj,jk,Krhs) + (zhpj(ji,jj,jk) + zvap) * vmask(ji,jj,jk)
      END_3D
      !
   END SUBROUTINE hpg_isf


   SUBROUTINE hpg_djc( kt, Kmm, puu, pvv, Krhs )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE hpg_djc  ***
      !!
      !! ** Method  :   Density Jacobian with Cubic polynomial scheme
      !!
      !! Reference: Shchepetkin and McWilliams, J. Geophys. Res., 108(C3), 3090, 2003
      !!----------------------------------------------------------------------
      INTEGER                             , INTENT( in )  ::  kt          ! ocean time-step index
      INTEGER                             , INTENT( in )  ::  Kmm, Krhs   ! ocean time level indices
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpt), INTENT(inout) ::  puu, pvv    ! ocean velocities and RHS of momentum equation
      !!
      INTEGER  ::   ji, jj, jk          ! dummy loop indices
      INTEGER  ::   iktb, iktt          ! jk indices at tracer points for top and bottom points 
      REAL(wp) ::   zcoef0, zep, cffw   ! temporary scalars
      REAL(wp) ::   z_grav_10, z1_12, z1_cff
      REAL(wp) ::   cffu, cffx          !    "         "
      REAL(wp) ::   cffv, cffy          !    "         "
      LOGICAL  ::   ll_tmp1, ll_tmp2    ! local logical variables
      REAL(wp), DIMENSION(A2D(nn_hls),jpk) ::   zhpi, zhpj

      REAL(wp), DIMENSION(A2D(nn_hls),jpk) ::   zdzx, zdzy, zdzz                          ! Primitive grid differences ('delta_xyz')
      REAL(wp), DIMENSION(A2D(nn_hls),jpk) ::   zdz_i, zdz_j, zdz_k                       ! Harmonic average of primitive grid differences ('d_xyz')
      REAL(wp), DIMENSION(A2D(nn_hls),jpk) ::   zdrhox, zdrhoy, zdrhoz                    ! Primitive rho differences ('delta_rho')
      REAL(wp), DIMENSION(A2D(nn_hls),jpk) ::   zdrho_i, zdrho_j, zdrho_k                 ! Harmonic average of primitive rho differences ('d_rho')
      REAL(wp), DIMENSION(A2D(nn_hls),jpk) ::   z_rho_i, z_rho_j, z_rho_k                 ! Face intergrals
      REAL(wp), DIMENSION(A2D(nn_hls))     ::   zz_dz_i, zz_dz_j, zz_drho_i, zz_drho_j    ! temporary arrays
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   zcpx, zcpy   !W/D pressure filter
      !!----------------------------------------------------------------------
      !
      IF( ln_wd_il ) THEN
         ALLOCATE( zcpx(A2D(nn_hls)) , zcpy(A2D(nn_hls)) )
        DO_2D( 0, 0, 0, 0 )
          ll_tmp1 = MIN(  ssh(ji,jj,Kmm)              ,  ssh(ji+1,jj,Kmm) ) >                &
               &    MAX( -ht_0(ji,jj)              , -ht_0(ji+1,jj) ) .AND.            &
               &    MAX(  ssh(ji,jj,Kmm) + ht_0(ji,jj),  ssh(ji+1,jj,Kmm) + ht_0(ji+1,jj) )  &
               &                                                      > rn_wdmin1 + rn_wdmin2
          ll_tmp2 = ( ABS( ssh(ji,jj,Kmm)             -  ssh(ji+1,jj,Kmm) ) > 1.E-12 ) .AND. (        &
               &    MAX(   ssh(ji,jj,Kmm)             ,  ssh(ji+1,jj,Kmm) ) >                &
               &    MAX(  -ht_0(ji,jj)             , -ht_0(ji+1,jj) ) + rn_wdmin1 + rn_wdmin2 )
          IF(ll_tmp1) THEN
            zcpx(ji,jj) = 1.0_wp
          ELSE IF(ll_tmp2) THEN
            ! no worries about  ssh(ji+1,jj,Kmm) -  ssh(ji  ,jj,Kmm) = 0, it won't happen ! here
            zcpx(ji,jj) = ABS( (ssh(ji+1,jj,Kmm) + ht_0(ji+1,jj) - ssh(ji,jj,Kmm) - ht_0(ji,jj)) &
                        &    / (ssh(ji+1,jj,Kmm) - ssh(ji  ,jj,Kmm)) )
          ELSE
            zcpx(ji,jj) = 0._wp
          END IF
   
          ll_tmp1 = MIN(  ssh(ji,jj,Kmm)              ,  ssh(ji,jj+1,Kmm) ) >                &
               &    MAX( -ht_0(ji,jj)              , -ht_0(ji,jj+1) ) .AND.            &
               &    MAX(  ssh(ji,jj,Kmm) + ht_0(ji,jj),  ssh(ji,jj+1,Kmm) + ht_0(ji,jj+1) )  &
               &                                                      > rn_wdmin1 + rn_wdmin2
          ll_tmp2 = ( ABS( ssh(ji,jj,Kmm)             -  ssh(ji,jj+1,Kmm) ) > 1.E-12 ) .AND. (        &
               &    MAX(   ssh(ji,jj,Kmm)             ,  ssh(ji,jj+1,Kmm) ) >                &
               &    MAX(  -ht_0(ji,jj)             , -ht_0(ji,jj+1) ) + rn_wdmin1 + rn_wdmin2 )

          IF(ll_tmp1) THEN
            zcpy(ji,jj) = 1.0_wp
          ELSE IF(ll_tmp2) THEN
            ! no worries about  ssh(ji,jj+1,Kmm) -  ssh(ji,jj  ,Kmm) = 0, it won't happen ! here
            zcpy(ji,jj) = ABS( (ssh(ji,jj+1,Kmm) + ht_0(ji,jj+1) - ssh(ji,jj,Kmm) - ht_0(ji,jj)) &
                        &    / (ssh(ji,jj+1,Kmm) - ssh(ji,jj  ,Kmm)) )
          ELSE
            zcpy(ji,jj) = 0._wp
          END IF
        END_2D
      END IF

      IF( .NOT. l_istiled .OR. ntile == 1 )  THEN                       ! Do only on the first tile
         IF( kt == nit000 ) THEN
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'dyn:hpg_djc : hydrostatic pressure gradient trend'
            IF(lwp) WRITE(numout,*) '~~~~~~~~~~~   s-coordinate case, density Jacobian with cubic polynomial scheme'
         ENDIF
      ENDIF

      ! Local constant initialization
      zcoef0 = - grav * 0.5_wp
      z_grav_10  = grav / 10._wp
      z1_12  = 1.0_wp / 12._wp

      !----------------------------------------------------------------------------------------
      !  1. compute and store elementary vertical differences in provisional arrays 
      !----------------------------------------------------------------------------------------

!!bug gm   Not a true bug, but... zdzz=e3w  for zdzx, zdzy verify what it is really

      DO_3D( 1, 1, 1, 1, 2, jpkm1 )  
         zdrhoz(ji,jj,jk) =   rhd    (ji  ,jj  ,jk) - rhd    (ji,jj,jk-1)
         zdzz  (ji,jj,jk) = - gde3w(ji  ,jj  ,jk) + gde3w(ji,jj,jk-1)
      END_3D

      !-------------------------------------------------------------------------
      ! 2. compute harmonic averages for vertical differences using eq. 5.18
      !-------------------------------------------------------------------------
      zep = 1.e-15

!! mb zdrho_k, zdz_k, zdrho_i, zdz_i, zdrho_j, zdz_j re-centred about the point (ji,jj,jk) 
      zdrho_k(:,:,:) = 0._wp
      zdz_k  (:,:,:) = 0._wp

      DO_3D( 1, 1, 1, 1, 2, jpk-2 )
         cffw = MAX( 2._wp * zdrhoz(ji,jj,jk) * zdrhoz(ji,jj,jk+1), 0._wp )
         z1_cff = zdrhoz(ji,jj,jk) + zdrhoz(ji,jj,jk+1)
         zdrho_k(ji,jj,jk) = cffw / SIGN( MAX( ABS(z1_cff), zep ), z1_cff )
         zdz_k(ji,jj,jk) = 2._wp *   zdzz(ji,jj,jk) * zdzz(ji,jj,jk+1)   &
            &                  / ( zdzz(ji,jj,jk) + zdzz(ji,jj,jk+1) )
      END_3D

      !----------------------------------------------------------------------------------
      ! 3. apply boundary conditions at top and bottom using 5.36-5.37
      !----------------------------------------------------------------------------------

! mb for sea-ice shelves we will need to re-write this upper boundary condition in the same form as the lower boundary condition
      DO_2D( 1, 1, 1, 1 )
         zdrho_k(ji,jj,1) = aco_bc_rhd_srf * ( rhd  (ji,jj,2) - rhd  (ji,jj,1) ) - bco_bc_rhd_srf * zdrho_k(ji,jj,2)
         zdz_k  (ji,jj,1) = aco_bc_z_srf   * (-gde3w(ji,jj,2) + gde3w(ji,jj,1) ) - bco_bc_z_srf   * zdz_k  (ji,jj,2)
      END_2D

      DO_2D( 1, 1, 1, 1 )
         IF ( mbkt(ji,jj)>1 ) THEN
            iktb = mbkt(ji,jj)
            zdrho_k(ji,jj,iktb) = aco_bc_rhd_bot * (    rhd(ji,jj,iktb) -   rhd(ji,jj,iktb-1) ) - bco_bc_rhd_bot * zdrho_k(ji,jj,iktb-1)
            zdz_k  (ji,jj,iktb) = aco_bc_z_bot   * ( -gde3w(ji,jj,iktb) + gde3w(ji,jj,iktb-1) ) - bco_bc_z_bot   * zdz_k  (ji,jj,iktb-1) 
         END IF
      END_2D

      IF ( ln_dbg_hpg ) CALL dbg_3dr( '3. zdz_k', zdz_k ) 
      IF ( ln_dbg_hpg ) CALL dbg_3dr( '3. zdrho_k', zdrho_k ) 

      !--------------------------------------------------------------
      ! 4. Compute side face integrals
      !-------------------------------------------------------------

!! ssh replaces e3w_n ; gde3w is a depth; the formulae involve heights  
!! rho_k stores grav * FX / rho_0  

      !--------------------------------------------------------------
      ! 4. a) Upper half of top-most grid box, compute and store
      !-------------------------------------------------------------
! *** AY note: ssh(ji,jj,Kmm) + gde3w(ji,jj,1) = e3w(ji,jj,1,Kmm)
      DO_2D( 0, 1, 0, 1)
         z_rho_k(ji,jj,1) =  grav * ( ssh(ji,jj,Kmm) + gde3w(ji,jj,1) )                        & 
            &                     * (  rhd(ji,jj,1)                                        &
            &                     + 0.5_wp * (   rhd    (ji,jj,2) - rhd    (ji,jj,1) ) &
            &                              * (   ssh   (ji,jj,Kmm) + gde3w(ji,jj,1) )          &
            &                              / ( - gde3w(ji,jj,2) + gde3w(ji,jj,1) )  )
      END_2D

      !--------------------------------------------------------------
      ! 4. b) Interior faces, compute and store
      !-------------------------------------------------------------

      DO_3D( 0, 1, 0, 1, 2, jpkm1 )
         z_rho_k(ji,jj,jk) = zcoef0 * (   rhd    (ji,jj,jk) + rhd    (ji,jj,jk-1) )                                   &
            &                       * ( - gde3w(ji,jj,jk) + gde3w(ji,jj,jk-1) )                                               &
            &                       + z_grav_10 * (                                                                           &
            &     (   zdrho_k  (ji,jj,jk) - zdrho_k  (ji,jj,jk-1) )                                                           &
            &   * ( - gde3w(ji,jj,jk) + gde3w(ji,jj,jk-1) - z1_12 * ( zdz_k  (ji,jj,jk) + zdz_k  (ji,jj,jk-1) ) )             &
            &   - ( zdz_k    (ji,jj,jk) - zdz_k    (ji,jj,jk-1) )                                                             &
            &   * ( rhd    (ji,jj,jk) - rhd    (ji,jj,jk-1) - z1_12 * ( zdrho_k(ji,jj,jk) + zdrho_k(ji,jj,jk-1) ) )   &
            &                             )
      END_3D

      IF ( ln_dbg_hpg ) CALL dbg_3dr( '4. z_rho_k', z_rho_k ) 

      !----------------------------------------------------------------------------------------
      !  5. compute and store elementary horizontal differences in provisional arrays 
      !----------------------------------------------------------------------------------------
      zdrhox(:,:,:) = 0._wp
      zdzx  (:,:,:) = 0._wp
      zdrhoy(:,:,:) = 0._wp
      zdzy  (:,:,:) = 0._wp

      DO_3D( nn_hls-1, nn_hls-1, nn_hls-1, nn_hls-1, 1, jpkm1 )
         zdrhox(ji,jj,jk) = rhd  (ji+1,jj  ,jk) - rhd  (ji  ,jj  ,jk)
         zdzx  (ji,jj,jk) = gde3w(ji  ,jj  ,jk) - gde3w(ji+1,jj  ,jk)
         zdrhoy(ji,jj,jk) = rhd  (ji  ,jj+1,jk) - rhd  (ji  ,jj  ,jk)
         zdzy  (ji,jj,jk) = gde3w(ji  ,jj  ,jk) - gde3w(ji  ,jj+1,jk)
      END_3D

      IF( nn_hls == 1 ) CALL lbc_lnk( 'dynhpg', zdrhox, 'U', -1._wp, zdzx, 'U', -1._wp, zdrhoy, 'V', -1._wp, zdzy, 'V', -1._wp )

      !-------------------------------------------------------------------------
      ! 6. compute harmonic averages using eq. 5.18
      !-------------------------------------------------------------------------

      DO_3D( 0, 1, 0, 1, 1, jpkm1 )
         cffu = MAX( 2._wp * zdrhox(ji-1,jj,jk) * zdrhox(ji,jj,jk), 0._wp )
         z1_cff = zdrhox(ji-1,jj,jk) + zdrhox(ji,jj,jk)
         zdrho_i(ji,jj,jk) = cffu / SIGN( MAX( ABS(z1_cff), zep ), z1_cff )

         cffx = MAX( 2._wp * zdzx(ji-1,jj,jk)   * zdzx(ji,jj,jk), 0._wp )
         z1_cff = zdzx(ji-1,jj,jk)   + zdzx(ji,jj,jk)
         zdz_i(ji,jj,jk)   = cffx / SIGN( MAX( ABS(z1_cff), zep ), z1_cff )

         cffv = MAX( 2._wp * zdrhoy(ji,jj-1,jk) * zdrhoy(ji,jj,jk), 0._wp )
         z1_cff = zdrhoy(ji,jj-1,jk) + zdrhoy(ji,jj,jk)
         zdrho_j(ji,jj,jk) = cffv / SIGN( MAX( ABS(z1_cff), zep ), z1_cff )

         cffy = MAX( 2._wp * zdzy(ji,jj-1,jk)   * zdzy(ji,jj,jk), 0._wp )
         z1_cff = zdzy(ji,jj-1,jk)   + zdzy(ji,jj,jk)
         zdz_j(ji,jj,jk)   = cffy / SIGN( MAX( ABS(z1_cff), zep ), z1_cff )
      END_3D
      
!!! Note that zdzx, zdzy, zdzz, zdrhox, zdrhoy and zdrhoz should NOT be used beyond this point      

      !----------------------------------------------------------------------------------
      ! 6B. apply boundary conditions at side boundaries using 5.36-5.37
      !----------------------------------------------------------------------------------

      DO jk = 1, jpkm1
         zz_drho_i(:,:) = zdrho_i(:,:,jk)
         zz_dz_i  (:,:) = zdz_i  (:,:,jk)
         zz_drho_j(:,:) = zdrho_j(:,:,jk)
         zz_dz_j  (:,:) = zdz_j  (:,:,jk)
         ! Walls coming from left: should check from 2 to jpi-1 (and jpj=2-jpj)
         DO_2D( 0, 0, 0, 1 )
            IF ( umask(ji,jj,jk) > 0.5_wp .AND. umask(ji-1,jj,jk) < 0.5_wp .AND. umask(ji+1,jj,jk) > 0.5_wp)  THEN
               zz_drho_i(ji,jj) = aco_bc_rhd_hor * ( rhd    (ji+1,jj,jk) - rhd  (ji,jj,jk) ) - bco_bc_rhd_hor * zdrho_i(ji+1,jj,jk)
               zz_dz_i  (ji,jj) = aco_bc_z_hor   * (-gde3w(ji+1,jj,jk)   + gde3w(ji,jj,jk) ) - bco_bc_z_hor   * zdz_i  (ji+1,jj,jk)
            END IF
         END_2D
         ! Walls coming from right: should check from 3 to jpi (and jpj=2-jpj)
         DO_2D( -1, 1, 0, 1 )
            IF ( umask(ji,jj,jk) < 0.5_wp .AND. umask(ji-1,jj,jk) > 0.5_wp .AND. umask(ji-2,jj,jk) > 0.5_wp) THEN
               zz_drho_i(ji,jj) = aco_bc_rhd_hor * ( rhd  (ji,jj,jk) - rhd  (ji-1,jj,jk) ) - bco_bc_rhd_hor * zdrho_i(ji-1,jj,jk)
               zz_dz_i  (ji,jj) = aco_bc_z_hor   * (-gde3w(ji,jj,jk) + gde3w(ji-1,jj,jk) ) - bco_bc_z_hor   * zdz_i  (ji-1,jj,jk)
            END IF
         END_2D
         ! Walls coming from left: should check from 2 to jpj-1 (and jpi=2-jpi)
         DO_2D( 0, 1, 0, 0 )
            IF ( vmask(ji,jj,jk) > 0.5_wp .AND. vmask(ji,jj-1,jk) < 0.5_wp .AND. vmask(ji,jj+1,jk) > 0.5_wp)  THEN
               zz_drho_j(ji,jj) = aco_bc_rhd_hor * ( rhd  (ji,jj+1,jk) - rhd  (ji,jj,jk) ) - bco_bc_rhd_hor * zdrho_j(ji,jj+1,jk)
               zz_dz_j  (ji,jj) = aco_bc_z_hor   * (-gde3w(ji,jj+1,jk) + gde3w(ji,jj,jk) ) - bco_bc_z_hor   * zdz_j  (ji,jj+1,jk)
            END IF
         END_2D
         ! Walls coming from right: should check from 3 to jpj (and jpi=2-jpi)
         DO_2D( 0, 1, -1, 1 )
            IF ( vmask(ji,jj,jk) < 0.5_wp .AND. vmask(ji,jj-1,jk) > 0.5_wp .AND. vmask(ji,jj-2,jk) > 0.5_wp) THEN
               zz_drho_j(ji,jj) = aco_bc_rhd_hor * ( rhd  (ji,jj,jk) - rhd  (ji,jj-1,jk) ) - bco_bc_rhd_hor * zdrho_j(ji,jj-1,jk)
               zz_dz_j  (ji,jj) = aco_bc_z_hor   * (-gde3w(ji,jj,jk) + gde3w(ji,jj-1,jk) ) - bco_bc_z_hor   * zdz_j  (ji,jj-1,jk)
            END IF
         END_2D
         zdrho_i(:,:,jk) = zz_drho_i(:,:)
         zdz_i  (:,:,jk) = zz_dz_i  (:,:)
         zdrho_j(:,:,jk) = zz_drho_j(:,:)
         zdz_j  (:,:,jk) = zz_dz_j  (:,:)
      END DO ! jk 
      
      IF ( ln_dbg_hpg ) THEN 
         CALL dbg_3dr( '6. zdrho_i', zdrho_i ) 
         CALL dbg_3dr( '6. zdrho_j', zdrho_j ) 
      END IF       

      !--------------------------------------------------------------
      ! 7. Calculate integrals on upper/lower faces  
      !-------------------------------------------------------------

      DO_3D( 0, 0, 0, 0, 1, jpkm1 )
! two -ve signs cancel in next two lines (within zcoef0 and because gde3w is a depth not a height)
         z_rho_i(ji,jj,jk) = zcoef0 * ( rhd    (ji+1,jj,jk) + rhd    (ji,jj,jk) )                                       &
             &                    * ( gde3w(ji+1,jj,jk) - gde3w(ji,jj,jk) )                                    
         IF ( umask(ji-1, jj, jk) > 0.5 .OR. umask(ji+1, jj, jk) > 0.5 ) THEN
            z_rho_i(ji,jj,jk) = z_rho_i(ji,jj,jk) - z_grav_10 * (                                                               &
             &     (   zdrho_i  (ji+1,jj,jk) - zdrho_i  (ji,jj,jk) )                                                            &
             &   * ( - gde3w(ji+1,jj,jk) + gde3w(ji,jj,jk) - z1_12 * ( zdz_i  (ji+1,jj,jk) + zdz_i  (ji,jj,jk) ) )              &
             &   - (   zdz_i    (ji+1,jj,jk) - zdz_i    (ji,jj,jk) )                                                            &
             &   * (   rhd    (ji+1,jj,jk) - rhd    (ji,jj,jk) - z1_12 * ( zdrho_i(ji+1,jj,jk) + zdrho_i(ji,jj,jk) ) )  &
             &                                               )
         END IF
  
         z_rho_j(ji,jj,jk) = zcoef0 * ( rhd    (ji,jj+1,jk) + rhd    (ji,jj,jk) )                                       &
             &                    * ( gde3w(ji,jj+1,jk) - gde3w(ji,jj,jk) )                                  
         IF ( vmask(ji, jj-1, jk) > 0.5 .OR. vmask(ji, jj+1, jk) > 0.5 ) THEN
            z_rho_j(ji,jj,jk) = z_rho_j(ji,jj,jk) - z_grav_10 * (                                                               &
             &     (   zdrho_j  (ji,jj+1,jk) - zdrho_j  (ji,jj,jk) )                                                            &
             &   * ( - gde3w(ji,jj+1,jk) + gde3w(ji,jj,jk) - z1_12 * ( zdz_j  (ji,jj+1,jk) + zdz_j  (ji,jj,jk) ) )              &
             &   - (   zdz_j    (ji,jj+1,jk) - zdz_j    (ji,jj,jk) )                                                            &
             &   * (   rhd    (ji,jj+1,jk) - rhd    (ji,jj,jk) - z1_12 * ( zdrho_j(ji,jj+1,jk) + zdrho_j(ji,jj,jk) ) )  &
             &                                                 )
         END IF
      END_3D

      IF ( ln_dbg_hpg ) THEN 
         CALL dbg_3dr( '7. z_rho_i', z_rho_i ) 
         CALL dbg_3dr( '7. z_rho_j', z_rho_j ) 
      END IF       

      !--------------------------------------------------------------
      ! 8. Integrate in the vertical   
      !-------------------------------------------------------------
      !
      ! ---------------
      !  Surface value
      ! ---------------
      DO_2D( 0, 0, 0, 0 )
         zhpi(ji,jj,1) = ( z_rho_k(ji,jj,1) - z_rho_k(ji+1,jj  ,1) - z_rho_i(ji,jj,1) ) * r1_e1u(ji,jj)
         zhpj(ji,jj,1) = ( z_rho_k(ji,jj,1) - z_rho_k(ji  ,jj+1,1) - z_rho_j(ji,jj,1) ) * r1_e2v(ji,jj)
         IF( ln_wd_il ) THEN
           zhpi(ji,jj,1) = zhpi(ji,jj,1) * zcpx(ji,jj)
           zhpj(ji,jj,1) = zhpj(ji,jj,1) * zcpy(ji,jj) 
         ENDIF
         ! add to the general momentum trend
         puu(ji,jj,1,Krhs) = puu(ji,jj,1,Krhs) + zhpi(ji,jj,1)
         pvv(ji,jj,1,Krhs) = pvv(ji,jj,1,Krhs) + zhpj(ji,jj,1)
      END_2D

      ! ----------------
      !  interior value   (2=<jk=<jpkm1)
      ! ----------------
      DO_3D( 0, 0, 0, 0, 2, jpkm1 )
         ! hydrostatic pressure gradient along s-surfaces
         zhpi(ji,jj,jk) = zhpi(ji,jj,jk-1)                                                     &
            &           + (  ( z_rho_k(ji,jj,jk) - z_rho_k(ji+1,jj,jk  ) )                     &
            &              - ( z_rho_i(ji,jj,jk) - z_rho_i(ji  ,jj,jk-1) )  ) * r1_e1u(ji,jj)
         zhpj(ji,jj,jk) = zhpj(ji,jj,jk-1)                                                     &
            &           + (  ( z_rho_k(ji,jj,jk) - z_rho_k(ji,jj+1,jk  ) )                     &
            &               -( z_rho_j(ji,jj,jk) - z_rho_j(ji,jj  ,jk-1) )  ) * r1_e2v(ji,jj)
         IF( ln_wd_il ) THEN
           zhpi(ji,jj,jk) = zhpi(ji,jj,jk) * zcpx(ji,jj)
           zhpj(ji,jj,jk) = zhpj(ji,jj,jk) * zcpy(ji,jj) 
         ENDIF
         ! add to the general momentum trend
         puu(ji,jj,jk,Krhs) = puu(ji,jj,jk,Krhs) + zhpi(ji,jj,jk)
         pvv(ji,jj,jk,Krhs) = pvv(ji,jj,jk,Krhs) + zhpj(ji,jj,jk)
      END_3D
      !

      IF ( ln_dbg_hpg ) THEN 
         CALL dbg_3dr( '8. zhpi', zhpi ) 
         CALL dbg_3dr( '8. zhpj', zhpj ) 
      END IF       
  
      IF( ln_wd_il )   DEALLOCATE( zcpx, zcpy )
      !
   END SUBROUTINE hpg_djc


   SUBROUTINE hpg_prj( kt, Kmm, puu, pvv, Krhs )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE hpg_prj  ***
      !!
      !! ** Method  :   s-coordinate case.
      !!      A Pressure-Jacobian horizontal pressure gradient method
      !!      based on the constrained cubic-spline interpolation for
      !!      all vertical coordinate systems
      !!
      !! ** Action : - Update (puu(:,:,:,Krhs),pvv(:,:,:,Krhs)) with the now hydrastatic pressure trend
      !!----------------------------------------------------------------------
      INTEGER, PARAMETER  :: polynomial_type = 1    ! 1: cubic spline, 2: linear
      INTEGER                             , INTENT( in )  ::  kt          ! ocean time-step index
      INTEGER                             , INTENT( in )  ::  Kmm, Krhs   ! ocean time level indices
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpt), INTENT(inout) ::  puu, pvv    ! ocean velocities and RHS of momentum equation
      !!
      INTEGER  ::   ji, jj, jk, jkk                 ! dummy loop indices
      REAL(wp) ::   zcoef0, znad                    ! local scalars
      !
      !! The local variables for the correction term
      INTEGER  :: jk1, jis, jid, jjs, jjd
      LOGICAL  :: ll_tmp1, ll_tmp2                  ! local logical variables
      REAL(wp) :: zuijk, zvijk, zpwes, zpwed, zpnss, zpnsd, zdeps
      REAL(wp) :: zrhdt1
      REAL(wp) :: zdpdx1, zdpdx2, zdpdy1, zdpdy2
      REAL(wp), DIMENSION(A2D(nn_hls))     ::   zpgu, zpgv   ! 2D workspace
      REAL(wp), DIMENSION(A2D(nn_hls))     ::   zsshu_n, zsshv_n
      REAL(wp), DIMENSION(A2D(nn_hls),jpk) ::   zdept, zrhh
      REAL(wp), DIMENSION(A2D(nn_hls),jpk) ::   zhpi, zu, zv, fsp, xsp, asp, bsp, csp, dsp
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   zcpx, zcpy   !W/D pressure filter
      !!----------------------------------------------------------------------
      !
      IF( .NOT. l_istiled .OR. ntile == 1 )  THEN                       ! Do only on the first tile
         IF( kt == nit000 ) THEN
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'dyn:hpg_prj : hydrostatic pressure gradient trend'
            IF(lwp) WRITE(numout,*) '~~~~~~~~~~~   s-coordinate case, cubic spline pressure Jacobian'
         ENDIF
      ENDIF

      ! Local constant initialization
      zcoef0 = - grav
      znad = 1._wp
      IF( ln_linssh )   znad = 1._wp
      !
      ! ---------------
      !  Surface pressure gradient to be removed
      ! ---------------
      DO_2D( 0, 0, 0, 0 )
         zpgu(ji,jj) = - grav * ( ssh(ji+1,jj,Kmm) - ssh(ji,jj,Kmm) ) * r1_e1u(ji,jj)
         zpgv(ji,jj) = - grav * ( ssh(ji,jj+1,Kmm) - ssh(ji,jj,Kmm) ) * r1_e2v(ji,jj)
      END_2D
      !
      IF( ln_wd_il ) THEN
         ALLOCATE( zcpx(A2D(nn_hls)) , zcpy(A2D(nn_hls)) )
         DO_2D( 0, 0, 0, 0 )
            ll_tmp1 = MIN(   ssh(ji,jj,Kmm)              ,   ssh(ji+1,jj,Kmm)                 ) >       &
               &      MAX( -ht_0(ji,jj)                  , -ht_0(ji+1,jj)                     ) .AND.   &
               &      MAX(   ssh(ji,jj,Kmm) + ht_0(ji,jj),   ssh(ji+1,jj,Kmm) + ht_0(ji+1,jj) ) >       &
               &      rn_wdmin1 + rn_wdmin2
            ll_tmp2 = ( ABS(   ssh(ji,jj,Kmm) -   ssh(ji+1,jj,Kmm) ) > 1.E-12 ) .AND.                   &
               &      ( MAX(   ssh(ji,jj,Kmm) ,   ssh(ji+1,jj,Kmm) ) >                                  &
               &        MAX( -ht_0(ji,jj)     , -ht_0(ji+1,jj)     ) + rn_wdmin1 + rn_wdmin2 )

            IF(ll_tmp1) THEN
               zcpx(ji,jj) = 1.0_wp
            ELSE IF(ll_tmp2) THEN
               ! no worries about  ssh(ji+1,jj,Kmm) -  ssh(ji  ,jj,Kmm) = 0, it won't happen ! here
               zcpx(ji,jj) = ABS( (ssh(ji+1,jj,Kmm) + ht_0(ji+1,jj) - ssh(ji,jj,Kmm) - ht_0(ji,jj)) &
                           &    / (ssh(ji+1,jj,Kmm) -  ssh(ji  ,jj,Kmm)) )
               zcpx(ji,jj) = MAX(MIN( zcpx(ji,jj) , 1.0_wp),0.0_wp)
            ELSE
               zcpx(ji,jj) = 0._wp
            END IF

            ll_tmp1 = MIN(   ssh(ji,jj,Kmm)              ,   ssh(ji,jj+1,Kmm)                 ) >       &
               &      MAX( -ht_0(ji,jj)                  , -ht_0(ji,jj+1)                     ) .AND.   &
               &      MAX(   ssh(ji,jj,Kmm) + ht_0(ji,jj),   ssh(ji,jj+1,Kmm) + ht_0(ji,jj+1) ) >       &
               &      rn_wdmin1 + rn_wdmin2
            ll_tmp2 = ( ABS(   ssh(ji,jj,Kmm) -   ssh(ji,jj+1,Kmm) ) > 1.E-12 ) .AND.                   &
               &      ( MAX(   ssh(ji,jj,Kmm) ,   ssh(ji,jj+1,Kmm) ) >                                  &
               &        MAX( -ht_0(ji,jj)     , -ht_0(ji,jj+1)     ) + rn_wdmin1 + rn_wdmin2 )

            IF(ll_tmp1) THEN
               zcpy(ji,jj) = 1.0_wp
            ELSE IF(ll_tmp2) THEN
               ! no worries about  ssh(ji,jj+1,Kmm) -  ssh(ji,jj  ,Kmm) = 0, it won't happen ! here
               zcpy(ji,jj) = ABS( (ssh(ji,jj+1,Kmm) + ht_0(ji,jj+1) - ssh(ji,jj,Kmm) - ht_0(ji,jj)) &
                           &    / (ssh(ji,jj+1,Kmm) -  ssh(ji,jj  ,Kmm)) )
               zcpy(ji,jj) = MAX(MIN( zcpy(ji,jj) , 1.0_wp),0.0_wp)
            ELSE
               zcpy(ji,jj) = 0._wp
            ENDIF
         END_2D
      ENDIF

      ! Clean 3-D work arrays
      zhpi(:,:,:) = 0._wp
      zrhh(:,:,:) = rhd(A2D(nn_hls),:)

      ! Preparing vertical density profile "zrhh(:,:,:)" for hybrid-sco coordinate
      DO_2D( 1, 1, 1, 1 )
         jk = mbkt(ji,jj)
         IF(     jk <=  1   ) THEN   ;   zrhh(ji,jj,    :   ) = 0._wp
         ELSEIF( jk ==  2   ) THEN   ;   zrhh(ji,jj,jk+1:jpk) = rhd(ji,jj,jk)
         ELSEIF( jk < jpkm1 ) THEN
            DO jkk = jk+1, jpk
               zrhh(ji,jj,jkk) = interp1(gde3w(ji,jj,jkk  ), gde3w(ji,jj,jkk-1),   &
                  &                      gde3w(ji,jj,jkk-2), zrhh (ji,jj,jkk-1), zrhh(ji,jj,jkk-2))
            END DO
         ENDIF
      END_2D

      ! Transfer the depth of "T(:,:,:)" to vertical coordinate "zdept(:,:,:)"
      DO_2D( 1, 1, 1, 1 )
         zdept(ji,jj,1) = 0.5_wp * e3w(ji,jj,1,Kmm) - ssh(ji,jj,Kmm)
      END_2D

      DO_3D( 1, 1, 1, 1, 2, jpk )
         zdept(ji,jj,jk) = zdept(ji,jj,jk-1) + e3w(ji,jj,jk,Kmm)
      END_3D

      fsp(:,:,:) = zrhh (:,:,:)
      xsp(:,:,:) = zdept(:,:,:)

      ! Construct the vertical density profile with the
      ! constrained cubic spline interpolation
      ! rho(z) = asp + bsp*z + csp*z^2 + dsp*z^3
      CALL cspline( fsp, xsp, asp, bsp, csp, dsp, polynomial_type )

      ! Integrate the hydrostatic pressure "zhpi(:,:,:)" at "T(ji,jj,1)"
      DO_2D( 0, 1, 0, 1 )
         zrhdt1 = zrhh(ji,jj,1) - interp3( zdept(ji,jj,1), asp(ji,jj,1), bsp(ji,jj,1),  &
            &                                              csp(ji,jj,1), dsp(ji,jj,1) ) * 0.25_wp * e3w(ji,jj,1,Kmm)

         ! assuming linear profile across the top half surface layer
         zhpi(ji,jj,1) =  0.5_wp * e3w(ji,jj,1,Kmm) * zrhdt1
      END_2D

      ! Calculate the pressure "zhpi(:,:,:)" at "T(ji,jj,2:jpkm1)"
      DO_3D( 0, 1, 0, 1, 2, jpkm1 )
         zhpi(ji,jj,jk) = zhpi(ji,jj,jk-1) +                                  &
            &             integ_spline( zdept(ji,jj,jk-1), zdept(ji,jj,jk),   &
            &                           asp  (ji,jj,jk-1), bsp  (ji,jj,jk-1), &
            &                           csp  (ji,jj,jk-1), dsp  (ji,jj,jk-1)  )
      END_3D

      ! Z coordinate of U(ji,jj,1:jpkm1) and V(ji,jj,1:jpkm1)

      ! Prepare zsshu_n and zsshv_n
      DO_2D( 0, 0, 0, 0 )
!!gm BUG ?    if it is ssh at u- & v-point then it should be:
!          zsshu_n(ji,jj) = (e1e2t(ji,jj) * ssh(ji,jj,Kmm) + e1e2t(ji+1,jj) * ssh(ji+1,jj,Kmm)) * &
!                         & r1_e1e2u(ji,jj) * umask(ji,jj,1) * 0.5_wp 
!          zsshv_n(ji,jj) = (e1e2t(ji,jj) * ssh(ji,jj,Kmm) + e1e2t(ji,jj+1) * ssh(ji,jj+1,Kmm)) * &
!                         & r1_e1e2v(ji,jj) * vmask(ji,jj,1) * 0.5_wp 
!!gm not this:
         zsshu_n(ji,jj) = (e1e2u(ji,jj) * ssh(ji,jj,Kmm) + e1e2u(ji+1, jj) * ssh(ji+1,jj,Kmm)) * &
                        & r1_e1e2u(ji,jj) * umask(ji,jj,1) * 0.5_wp
         zsshv_n(ji,jj) = (e1e2v(ji,jj) * ssh(ji,jj,Kmm) + e1e2v(ji+1, jj) * ssh(ji,jj+1,Kmm)) * &
                        & r1_e1e2v(ji,jj) * vmask(ji,jj,1) * 0.5_wp
      END_2D

      DO_2D( 0, 0, 0, 0 )
         zu(ji,jj,1) = - ( e3u(ji,jj,1,Kmm) - zsshu_n(ji,jj) )
         zv(ji,jj,1) = - ( e3v(ji,jj,1,Kmm) - zsshv_n(ji,jj) )
      END_2D

      DO_3D( 0, 0, 0, 0, 2, jpkm1 )
         zu(ji,jj,jk) = zu(ji,jj,jk-1) - e3u(ji,jj,jk,Kmm)
         zv(ji,jj,jk) = zv(ji,jj,jk-1) - e3v(ji,jj,jk,Kmm)
      END_3D

      DO_3D( 0, 0, 0, 0, 1, jpkm1 )
         zu(ji,jj,jk) = zu(ji,jj,jk) + 0.5_wp * e3u(ji,jj,jk,Kmm)
         zv(ji,jj,jk) = zv(ji,jj,jk) + 0.5_wp * e3v(ji,jj,jk,Kmm)
      END_3D

      DO_3D( 0, 0, 0, 0, 1, jpkm1 )
         zu(ji,jj,jk) = MIN(  zu(ji,jj,jk) , MAX( -zdept(ji,jj,jk) , -zdept(ji+1,jj,jk) )  )
         zu(ji,jj,jk) = MAX(  zu(ji,jj,jk) , MIN( -zdept(ji,jj,jk) , -zdept(ji+1,jj,jk) )  )
         zv(ji,jj,jk) = MIN(  zv(ji,jj,jk) , MAX( -zdept(ji,jj,jk) , -zdept(ji,jj+1,jk) )  )
         zv(ji,jj,jk) = MAX(  zv(ji,jj,jk) , MIN( -zdept(ji,jj,jk) , -zdept(ji,jj+1,jk) )  )
      END_3D


      DO_3D( 0, 0, 0, 0, 1, jpkm1 )
         zpwes = 0._wp; zpwed = 0._wp
         zpnss = 0._wp; zpnsd = 0._wp
         zuijk = zu(ji,jj,jk)
         zvijk = zv(ji,jj,jk)

         !!!!!     for u equation
         IF( jk <= mbku(ji,jj) ) THEN
            IF( -zdept(ji+1,jj,jk) >= -zdept(ji,jj,jk) ) THEN
              jis = ji + 1; jid = ji
            ELSE
              jis = ji;     jid = ji +1
            ENDIF

            ! integrate the pressure on the shallow side
            jk1 = jk
            DO WHILE ( -zdept(jis,jj,jk1) > zuijk )
               IF( jk1 == mbku(ji,jj) ) THEN
                  zuijk = -zdept(jis,jj,jk1)
                  EXIT
               ENDIF
               zdeps = MIN(zdept(jis,jj,jk1+1), -zuijk)
               zpwes = zpwes +                                      &
                  integ_spline(zdept(jis,jj,jk1), zdeps,            &
                                 asp(jis,jj,jk1), bsp(jis,jj,jk1),  &
                                 csp(jis,jj,jk1), dsp(jis,jj,jk1))
               jk1 = jk1 + 1
            END DO

            ! integrate the pressure on the deep side
            jk1 = jk
            DO WHILE ( -zdept(jid,jj,jk1) < zuijk )
               IF( jk1 == 1 ) THEN
                  zdeps = zdept(jid,jj,1) + MIN(zuijk, ssh(jid,jj,Kmm)*znad)
                  zrhdt1 = zrhh(jid,jj,1) - interp3(zdept(jid,jj,1), asp(jid,jj,1), &
                                                    bsp(jid,jj,1)  , csp(jid,jj,1), &
                                                    dsp(jid,jj,1)) * zdeps
                  zpwed  = zpwed + 0.5_wp * (zrhh(jid,jj,1) + zrhdt1) * zdeps
                  EXIT
               ENDIF
               zdeps = MAX(zdept(jid,jj,jk1-1), -zuijk)
               zpwed = zpwed +                                        &
                  integ_spline(zdeps,             zdept(jid,jj,jk1),  &
                               asp(jid,jj,jk1-1), bsp(jid,jj,jk1-1),  &
                               csp(jid,jj,jk1-1), dsp(jid,jj,jk1-1) )
               jk1 = jk1 - 1
            END DO

            ! update the momentum trends in u direction
            zdpdx1 = zcoef0 * r1_e1u(ji,jj) * ( zhpi(ji+1,jj,jk) - zhpi(ji,jj,jk) )
            IF( .NOT.ln_linssh ) THEN
               zdpdx2 = zcoef0 * r1_e1u(ji,jj) * &
                  &    ( REAL(jis-jid, wp) * (zpwes + zpwed) + (ssh(ji+1,jj,Kmm)-ssh(ji,jj,Kmm)) )
            ELSE
               zdpdx2 = zcoef0 * r1_e1u(ji,jj) * REAL(jis-jid, wp) * (zpwes + zpwed)
            ENDIF
            IF( ln_wd_il ) THEN
               zdpdx1 = zdpdx1 * zcpx(ji,jj) * wdrampu(ji,jj)
               zdpdx2 = zdpdx2 * zcpx(ji,jj) * wdrampu(ji,jj)
            ENDIF
            puu(ji,jj,jk,Krhs) = puu(ji,jj,jk,Krhs) + (zdpdx1 + zdpdx2 - zpgu(ji,jj)) * umask(ji,jj,jk)
         ENDIF

         !!!!!     for v equation
         IF( jk <= mbkv(ji,jj) ) THEN
            IF( -zdept(ji,jj+1,jk) >= -zdept(ji,jj,jk) ) THEN
               jjs = jj + 1; jjd = jj
            ELSE
               jjs = jj    ; jjd = jj + 1
            ENDIF

            ! integrate the pressure on the shallow side
            jk1 = jk
            DO WHILE ( -zdept(ji,jjs,jk1) > zvijk )
               IF( jk1 == mbkv(ji,jj) ) THEN
                  zvijk = -zdept(ji,jjs,jk1)
                  EXIT
               ENDIF
               zdeps = MIN(zdept(ji,jjs,jk1+1), -zvijk)
               zpnss = zpnss +                                       &
                  integ_spline(zdept(ji,jjs,jk1), zdeps,             &
                               asp(ji,jjs,jk1),   bsp(ji,jjs,jk1),   &
                               csp(ji,jjs,jk1),   dsp(ji,jjs,jk1) )
              jk1 = jk1 + 1
            END DO

            ! integrate the pressure on the deep side
            jk1 = jk
            DO WHILE ( -zdept(ji,jjd,jk1) < zvijk )
               IF( jk1 == 1 ) THEN
                  zdeps = zdept(ji,jjd,1) + MIN(zvijk, ssh(ji,jjd,Kmm)*znad)
                  zrhdt1 = zrhh(ji,jjd,1) - interp3(zdept(ji,jjd,1), asp(ji,jjd,1), &
                                                    bsp(ji,jjd,1)  , csp(ji,jjd,1), &
                                                    dsp(ji,jjd,1) ) * zdeps
                  zpnsd  = zpnsd + 0.5_wp * (zrhh(ji,jjd,1) + zrhdt1) * zdeps
                  EXIT
               ENDIF
               zdeps = MAX(zdept(ji,jjd,jk1-1), -zvijk)
               zpnsd = zpnsd +                                        &
                  integ_spline(zdeps,             zdept(ji,jjd,jk1),  &
                               asp(ji,jjd,jk1-1), bsp(ji,jjd,jk1-1),  &
                               csp(ji,jjd,jk1-1), dsp(ji,jjd,jk1-1) )
               jk1 = jk1 - 1
            END DO

            ! update the momentum trends in v direction
            zdpdy1 = zcoef0 * r1_e2v(ji,jj) * ( zhpi(ji,jj+1,jk) - zhpi(ji,jj,jk) )
            IF( .NOT.ln_linssh ) THEN
               zdpdy2 = zcoef0 * r1_e2v(ji,jj) * &
                       ( REAL(jjs-jjd, wp) * (zpnss + zpnsd) + (ssh(ji,jj+1,Kmm)-ssh(ji,jj,Kmm)) )
            ELSE
               zdpdy2 = zcoef0 * r1_e2v(ji,jj) * REAL(jjs-jjd, wp) * (zpnss + zpnsd )
            ENDIF
            IF( ln_wd_il ) THEN
               zdpdy1 = zdpdy1 * zcpy(ji,jj) * wdrampv(ji,jj)
               zdpdy2 = zdpdy2 * zcpy(ji,jj) * wdrampv(ji,jj)
            ENDIF

            pvv(ji,jj,jk,Krhs) = pvv(ji,jj,jk,Krhs) + (zdpdy1 + zdpdy2 - zpgv(ji,jj)) * vmask(ji,jj,jk)
         ENDIF
         !
      END_3D
      !
      IF( ln_wd_il )   DEALLOCATE( zcpx, zcpy )
      !
   END SUBROUTINE hpg_prj


   SUBROUTINE cspline( fsp, xsp, asp, bsp, csp, dsp, polynomial_type )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE cspline  ***
      !!
      !! ** Purpose :   constrained cubic spline interpolation
      !!
      !! ** Method  :   f(x) = asp + bsp*x + csp*x^2 + dsp*x^3
      !!
      !! Reference: CJC Kruger, Constrained Cubic Spline Interpoltation
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(A2D(nn_hls),jpk), INTENT(in   ) ::   fsp, xsp           ! value and coordinate
      REAL(wp), DIMENSION(A2D(nn_hls),jpk), INTENT(  out) ::   asp, bsp, csp, dsp ! coefficients of the interpoated function
      INTEGER                             , INTENT(in   ) ::   polynomial_type    ! 1: cubic spline   ;   2: Linear
      !
      INTEGER  ::   ji, jj, jk                 ! dummy loop indices
      REAL(wp) ::   zdf1, zdf2, zddf1, zddf2, ztmp1, ztmp2, zdxtmp
      REAL(wp) ::   zdxtmp1, zdxtmp2, zalpha
      REAL(wp) ::   zdf(jpk)
      !!----------------------------------------------------------------------
      !
      IF (polynomial_type == 1) THEN     ! Constrained Cubic Spline
         DO_2D( 1, 1, 1, 1 )
            !!Fritsch&Butland's method, 1984 (preferred, but more computation)
            !    DO jk = 2, jpkm1-1
            !       zdxtmp1 = xsp(ji,jj,jk)   - xsp(ji,jj,jk-1)
            !       zdxtmp2 = xsp(ji,jj,jk+1) - xsp(ji,jj,jk)
            !       zdf1    = ( fsp(ji,jj,jk)   - fsp(ji,jj,jk-1) ) / zdxtmp1
            !       zdf2    = ( fsp(ji,jj,jk+1) - fsp(ji,jj,jk)   ) / zdxtmp2
            !
            !       zalpha = ( zdxtmp1 + 2._wp * zdxtmp2 ) / ( zdxtmp1 + zdxtmp2 ) / 3._wp
            !
            !       IF(zdf1 * zdf2 <= 0._wp) THEN
            !           zdf(jk) = 0._wp
            !       ELSE
            !         zdf(jk) = zdf1 * zdf2 / ( ( 1._wp - zalpha ) * zdf1 + zalpha * zdf2 )
            !       ENDIF
            !    END DO

            !!Simply geometric average
            DO jk = 2, jpk-2
               zdf1 = (fsp(ji,jj,jk  ) - fsp(ji,jj,jk-1)) / (xsp(ji,jj,jk  ) - xsp(ji,jj,jk-1))
               zdf2 = (fsp(ji,jj,jk+1) - fsp(ji,jj,jk  )) / (xsp(ji,jj,jk+1) - xsp(ji,jj,jk  ))

               IF(zdf1 * zdf2 <= 0._wp) THEN
                  zdf(jk) = 0._wp
               ELSE
                  zdf(jk) = 2._wp * zdf1 * zdf2 / (zdf1 + zdf2)
               ENDIF
            END DO

            zdf(1)     = 1.5_wp * ( fsp(ji,jj,2) - fsp(ji,jj,1) ) / &
                       &          ( xsp(ji,jj,2) - xsp(ji,jj,1) )           -  0.5_wp * zdf(2)
            zdf(jpkm1) = 1.5_wp * ( fsp(ji,jj,jpkm1) - fsp(ji,jj,jpkm1-1) ) / &
                       &          ( xsp(ji,jj,jpkm1) - xsp(ji,jj,jpkm1-1) ) - 0.5_wp * zdf(jpk - 2)

            DO jk = 1, jpk-2
               zdxtmp = xsp(ji,jj,jk+1) - xsp(ji,jj,jk)
               ztmp1  = (zdf(jk+1) + 2._wp * zdf(jk)) / zdxtmp
               ztmp2  =  6._wp * (fsp(ji,jj,jk+1) - fsp(ji,jj,jk)) / zdxtmp / zdxtmp
               zddf1  = -2._wp * ztmp1 + ztmp2
               ztmp1  = (2._wp * zdf(jk+1) + zdf(jk)) / zdxtmp
               zddf2  =  2._wp * ztmp1 - ztmp2

               dsp(ji,jj,jk) = (zddf2 - zddf1) / 6._wp / zdxtmp
               csp(ji,jj,jk) = ( xsp(ji,jj,jk+1) * zddf1 - xsp(ji,jj,jk)*zddf2 ) / 2._wp / zdxtmp
               bsp(ji,jj,jk) = ( fsp(ji,jj,jk+1) - fsp(ji,jj,jk) ) / zdxtmp - &
                             & csp(ji,jj,jk) * ( xsp(ji,jj,jk+1) + xsp(ji,jj,jk) ) - &
                             & dsp(ji,jj,jk) * ((xsp(ji,jj,jk+1) + xsp(ji,jj,jk))**2 - &
                             &                   xsp(ji,jj,jk+1) * xsp(ji,jj,jk))
               asp(ji,jj,jk) = fsp(ji,jj,jk) - xsp(ji,jj,jk) * (bsp(ji,jj,jk) + &
                             &                (xsp(ji,jj,jk) * (csp(ji,jj,jk) + &
                             &                 dsp(ji,jj,jk) * xsp(ji,jj,jk))))
            END DO
         END_2D

      ELSEIF ( polynomial_type == 2 ) THEN     ! Linear
         DO_3D( 1, 1, 1, 1, 1, jpk-2 )
            zdxtmp =xsp(ji,jj,jk+1) - xsp(ji,jj,jk)
            ztmp1 = fsp(ji,jj,jk+1) - fsp(ji,jj,jk)

            dsp(ji,jj,jk) = 0._wp
            csp(ji,jj,jk) = 0._wp
            bsp(ji,jj,jk) = ztmp1 / zdxtmp
            asp(ji,jj,jk) = fsp(ji,jj,jk) - bsp(ji,jj,jk) * xsp(ji,jj,jk)
         END_3D
         !
      ELSE
         CALL ctl_stop( 'invalid polynomial type in cspline' )
      ENDIF
      !
   END SUBROUTINE cspline


   FUNCTION interp1(x, xl, xr, fl, fr)  RESULT(f)
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE interp1  ***
      !!
      !! ** Purpose :   1-d linear interpolation
      !!
      !! ** Method  :   interpolation is straight forward
      !!                extrapolation is also permitted (no value limit)
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(in) ::  x, xl, xr, fl, fr
      REAL(wp)             ::  f ! result of the interpolation (extrapolation)
      REAL(wp)             ::  zdeltx
      !!----------------------------------------------------------------------
      !
      zdeltx = xr - xl
      IF( abs(zdeltx) <= 10._wp * EPSILON(x) ) THEN
         f = 0.5_wp * (fl + fr)
      ELSE
         f = ( (x - xl ) * fr - ( x - xr ) * fl ) / zdeltx
      ENDIF
      !
   END FUNCTION interp1


   FUNCTION interp2( x, a, b, c, d )  RESULT(f)
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE interp1  ***
      !!
      !! ** Purpose :   1-d constrained cubic spline interpolation
      !!
      !! ** Method  :  cubic spline interpolation
      !!
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(in) ::  x, a, b, c, d
      REAL(wp)             ::  f ! value from the interpolation
      !!----------------------------------------------------------------------
      !
      f = a + x* ( b + x * ( c + d * x ) )
      !
   END FUNCTION interp2


   FUNCTION interp3( x, a, b, c, d )  RESULT(f)
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE interp1  ***
      !!
      !! ** Purpose :   Calculate the first order of derivative of
      !!                a cubic spline function y=a+b*x+c*x^2+d*x^3
      !!
      !! ** Method  :   f=dy/dx=b+2*c*x+3*d*x^2
      !!
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(in) ::  x, a, b, c, d
      REAL(wp)             ::  f ! value from the interpolation
      !!----------------------------------------------------------------------
      !
      f = b + x * ( 2._wp * c + 3._wp * d * x)
      !
   END FUNCTION interp3


   FUNCTION integ_spline( xl, xr, a, b, c, d )  RESULT(f)
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE interp1  ***
      !!
      !! ** Purpose :   1-d constrained cubic spline integration
      !!
      !! ** Method  :  integrate polynomial a+bx+cx^2+dx^3 from xl to xr
      !!
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(in) ::  xl, xr, a, b, c, d
      REAL(wp)             ::  za1, za2, za3
      REAL(wp)             ::  f                   ! integration result
      !!----------------------------------------------------------------------
      !
      za1 = 0.5_wp * b
      za2 = c / 3.0_wp
      za3 = 0.25_wp * d
      !
      f  = xr * ( a + xr * ( za1 + xr * ( za2 + za3 * xr ) ) ) - &
         & xl * ( a + xl * ( za1 + xl * ( za2 + za3 * xl ) ) )
      !
   END FUNCTION integ_spline

   SUBROUTINE hpg_djr( kt, Kmm, puu, pvv, Krhs )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE hpg_djr  ***
      !!
      !! ** Method  :   Density Jacobian with Cubic polynomial scheme subtracting a local reference profile (pmr is profile minus reference) 
      !!                This code assumes a 2-point halo 
      !!
      !! Reference: Shchepetkin and McWilliams, J. Geophys. Res., 108(C3), 3090, 2003
      !!----------------------------------------------------------------------
      INTEGER                             , INTENT( in )  ::  kt          ! ocean time-step index
      INTEGER                             , INTENT( in )  ::  Kmm, Krhs   ! ocean time level indices
      REAL(wp), DIMENSION(A2D(nn_hls),jpk,jpt), INTENT(inout) ::  puu, pvv    ! ocean velocities and RHS of momentum equation
      !!
      INTEGER  ::   ji, jj, jk, jr      ! loop indices
      INTEGER  ::   jn_hor_pts          ! number of points in the horizontal stencil
      INTEGER  ::   j_uv                ! 1 for u-cell; 2 for v-cell 
      INTEGER  ::   iktb, iktt          ! jk indices at tracer points for top and bottom points 
      INTEGER  ::   jia,  jib, jja, jjb ! 
      INTEGER  ::   jir, jjr            ! reference (expand)
      INTEGER  ::   jio, jjo            ! offset (expand)
      
      REAL(wp) ::   z_grav_10, z1_12    ! constants 
      REAL(wp) ::   zhta, zhtb          ! temporary scalars
      REAL(wp) ::   zcoef0, zep, cffw   !    "         "
      REAL(wp) ::   aco, bco            !    "         "
      REAL(wp) ::   cffu, cffx, z1_cff  !    "         "
      REAL(wp) ::   cffv, cffy          !    "         "
      REAL(wp) ::   cff_31, cff_42      !    "         "
      LOGICAL  ::   ll_tmp1, ll_tmp2    ! local logical variables

      REAL(wp), DIMENSION(A2D(nn_hls),jpk)   ::   ztmp, zdz_i, zdz_j, zdz_k   ! Harmonic average of primitive grid differences
      REAL(wp), DIMENSION(A2D(nn_hls),jpk)   ::   zrhd_ref                    ! Reference density 
      REAL(wp), DIMENSION(A2D(nn_hls),jpk)   ::   zz_ref                      ! Reference heights
      REAL(wp), DIMENSION(A2D(nn_hls),jpk)   ::   zdrhd_k_ref                 ! Harmonic average of primitive differences for reference field
      REAL(wp), DIMENSION(A2D(nn_hls),jpk,4) ::   z_rhd_pmr                   ! rhd_prm = rhd - rhd_ref (values on the original grid) 
      REAL(wp), DIMENSION(A2D(nn_hls),jpk)   ::   zdrhd_k_pmr                 ! 
      REAL(wp), DIMENSION(A2D(nn_hls),jpk)   ::   z_lmr_k                     ! left minus right density integrals on vertical faces

      INTEGER,  DIMENSION(A2D(nn_hls))       ::   jk_bot_ref                  ! bottom levels in the reference profile
      REAL(wp), DIMENSION(A2D(nn_hls))       ::   zdzx, zdzy                  ! primitive differences in x and y
      REAL(wp), DIMENSION(A2D(nn_hls))       ::   zz_dz_i, zz_dz_j
      REAL(wp), DIMENSION(A2D(nn_hls))       ::   zdrhd_21, zdrhd_32, zdrhd_43
      REAL(wp), DIMENSION(A2D(nn_hls),2)     ::   zz_drhd_ij, zdrhd_ij
      REAL(wp), DIMENSION(A2D(nn_hls))       ::   z_low_ij, z_upp_ij
      REAL(wp), DIMENSION(A2D(nn_hls))       ::   zhpi, zhpj
      LOGICAL                                ::   ln_dbg                      ! DB for debugging

      !!----------------------------------------------------------------------

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn:hpg_djc : hydrostatic pressure gradient trend'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~   s-coordinate case, density Jacobian with cubic polynomial scheme'
      ENDIF

      ! Local constant initialization
      zcoef0 = - grav * 0.5_wp
      z_grav_10  = grav / 10._wp
      z1_12  = 1.0_wp / 12._wp
      zep = 1.e-15

      jn_hor_pts = 4   ! 4 points in the horizontal stencil 

      !------------------------------------------------------------------------------------------------------
      !  1. calculate harmonic averages of differences for z (grid heights) in i, j and k directions
      !------------------------------------------------------------------------------------------------------

      !------------------------------------------------------------------------------------------------------
      !  1.1 compute and store elementary vertical differences then harmonic averages for z using eqn 5.18 
      !      Full domain covered so that _ref profiles can be taken from zdz_k
      !------------------------------------------------------------------------------------------------------

      DO_3D( 2, 2, 2, 2, 2, jpk )  
         ztmp  (ji,jj,jk) = - gde3w(ji  ,jj  ,jk) + gde3w(ji,jj,jk-1)    
      END_3D

      zdz_k  (:,:,1) = 0._wp  ! jk index changed from : to 1 to make computationally less wasteful 

      DO_3D( 2, 2, 2, 2, 2, jpk-1 ) 
         zdz_k(ji,jj,jk) = 2._wp *   ztmp(ji,jj,jk) * ztmp(ji,jj,jk+1)   &
            &                  / ( ztmp(ji,jj,jk) + ztmp(ji,jj,jk+1) )
      END_3D

      !----------------------------------------------------------------------------------
      ! 1.2 apply boundary conditions at top and bottom using 5.36-5.37
      !----------------------------------------------------------------------------------

! mb for sea-ice shelves we will need to re-write this upper boundary condition in the same form as the lower boundary condition
      zdz_k  (:,:,1) = aco_bc_z_srf * (-gde3w(:,:,2) + gde3w(:,:,1) ) - bco_bc_z_srf * zdz_k  (:,:,2)

      DO_2D( 2, 2, 2, 2 )
         iktb = mbkt(ji,jj)
         IF ( iktb > 1 ) THEN
            zdz_k  (ji,jj,iktb) = aco_bc_z_bot * (-gde3w(ji,jj,iktb) + gde3w(ji,jj,iktb-1) ) - bco_bc_z_bot * zdz_k  (ji,jj,iktb-1) 
         END IF
      END_2D

      !----------------------------------------------------------------------------------------
      !  1.3 compute and store elementary horizontal differences then harmonic averages for z using eqn 5.18
      !----------------------------------------------------------------------------------------

      DO jk = 1, jpkm1 
         DO_2D( 1, 1, 1, 1 )
            zdzx  (ji,jj) = - gde3w(ji+1,jj  ,jk) + gde3w(ji,jj,jk  )
            zdzy  (ji,jj) = - gde3w(ji  ,jj+1,jk) + gde3w(ji,jj,jk  )
         END_2D

! this subroutine requires a 2-point halo      CALL lbc_lnk_multi( 'dynhpg', zdzx, 'U', 1., zdzy, 'V', 1. ) 

         DO_2D( 0, 1, 0, 1 )
            cffx = MAX( 2._wp * zdzx  (ji-1,jj) * zdzx  (ji,jj), 0._wp)
            z1_cff = zdzx(ji-1,jj)   + zdzx(ji,jj)
            zdz_i(ji,jj,jk)   = cffx / SIGN( MAX( ABS(z1_cff), zep ), z1_cff )

            cffy = MAX( 2._wp * zdzy  (ji  ,jj-1) * zdzy  (ji,jj), 0._wp)
            z1_cff = zdzy(ji,jj-1)   + zdzy(ji,jj)
            zdz_j(ji,jj,jk)   = cffy / SIGN( MAX( ABS(z1_cff), zep ), z1_cff )
         END_2D
      
      !----------------------------------------------------------------------------------
      ! 1.4 apply boundary conditions at sides using 5.36-5.37
      !----------------------------------------------------------------------------------

         DO_2D( 0, 1, 0, 1)
            ! Walls coming from left: should check from 2 to jpi-1 (and jpj=2-jpj)
            IF ( umask(ji,jj,jk) > 0.5_wp .AND. umask(ji-1,jj,jk) < 0.5_wp .AND. umask(ji+1,jj,jk) > 0.5_wp)  THEN  
               zdz_i(ji,jj,jk) = aco_bc_z_hor * (-gde3w(ji+1,jj,jk) + gde3w(ji,jj,jk) ) - bco_bc_z_hor * zz_dz_i(ji+1,jj)
            END IF
            ! Walls coming from right: should check from 3 to jpi (and jpj=2-jpj)
            IF ( umask(ji,jj,jk) < 0.5_wp .AND. umask(ji-1,jj,jk) > 0.5_wp .AND. umask(ji-2,jj,jk) > 0.5_wp) THEN
               zdz_i(ji,jj,jk) = aco_bc_z_hor * (-gde3w(ji,jj,jk) + gde3w(ji-1,jj,jk) ) - bco_bc_z_hor * zz_dz_i(ji-1,jj)
            END IF
            ! Walls coming from left: should check from 2 to jpj-1 (and jpi=2-jpi)
            IF ( vmask(ji,jj,jk) > 0.5_wp .AND. vmask(ji,jj-1,jk) < 0.5_wp .AND. vmask(ji,jj+1,jk) > 0.5_wp)  THEN
               zdz_j(ji,jj,jk) = aco_bc_z_hor * (-gde3w(ji,jj+1,jk) + gde3w(ji,jj,jk) ) - bco_bc_z_hor * zz_dz_j(ji,jj+1)
            END IF 
            ! Walls coming from right: should check from 3 to jpj (and jpi=2-jpi)
            IF ( vmask(ji,jj,jk) < 0.5_wp .AND. vmask(ji,jj-1,jk) > 0.5_wp .AND. vmask(ji,jj-2,jk) > 0.5_wp) THEN 
               zdz_j(ji,jj,jk) = aco_bc_z_hor * (-gde3w(ji,jj,jk) + gde3w(ji,jj-1,jk) ) - bco_bc_z_hor * zz_dz_j(ji,jj-1)
            END IF
         END_2D
      END DO  ! k 

      IF ( ln_dbg_hpg ) THEN 
         CALL dbg_3dr( '1.4 gde3w', gde3w ) 
         CALL dbg_3dr( '1.4 zdz_i', zdz_i ) 
         CALL dbg_3dr( '1.4 zdz_j', zdz_j ) 
         CALL dbg_3dr( '1.4 zdz_k', zdz_k ) 
         CALL dbg_2di( 'mbkt', mbkt) 
      END IF 
      !----------------------------------------------------------------------------------------
      ! 2.  Start loop over the u and v components and find the reference profile   
      !     The loop ends in section 5.4
      !----------------------------------------------------------------------------------------

      DO j_uv = 1, 2     ! j_uv = 1 is for u-cell ; j_uv = 2 for v-cell

      !----------------------------------------------------------------------------------------
      !  2.1 find reference profiles zrhd_ref and zz_ref and the bottom level of the reference profile  
      !----------------------------------------------------------------------------------------

         IF ( ln_dbg_hpg ) CALL dbg_2dr( '2.1 ht_0', ht_0 ) 
         !CALL calc_rhd_ref(j_uv, jn_hor_pts, gde3w, zrhd_ref, zz_ref, jk_bot_ref)  ! Uses rhd (IN) to 
                                                                                    ! calculate all other 
                                                                                    ! fields (OUT)
         CALL calc_rhd_ref(j_uv, jn_hor_pts, gdept(:,:,:,Kmm), zrhd_ref, zz_ref, jk_bot_ref)

         IF ( ln_dbg_hpg ) THEN 
            CALL dbg_3dr( '2.1 rhd', rhd ) 
            CALL dbg_3dr( '2.1 zrhd_ref', zrhd_ref ) 
            CALL dbg_3dr( '2.1 zz_ref', zz_ref ) 
            CALL dbg_2di( '2.1 jk_bot_ref', jk_bot_ref ) 
         END IF 

      !--------------------------------------------------------------------------------------------------------
      !  2.2  IF  ln_hpg_ref_ccs compute zdrhd_k_ref then set bcs at top & bottom  
      !                              (bcs not needed for simple cubic off-centred at boundaries)
      !--------------------------------------------------------------------------------------------------------

         IF ( ln_hpg_ref_ccs ) THEN
            CALL calc_drhd_k(zrhd_ref, jk_bot_ref, zdrhd_k_ref) 
            IF ( ln_dbg_hpg ) CALL dbg_3dr( '2.3 zdrhd_k_ref', zdrhd_k_ref )	    
         END IF  ! ln_hpg_ref_ccs 

      !----------------------------------------------------------------------------------------
      !  3. interpolate reference profiles to target profiles and form difference profiles z_rhd_pmr
      !----------------------------------------------------------------------------------------

         DO jr = 1, 4      
            IF ( j_uv == 1 ) THEN 
               jio = jr - 2         ! range of jio is -1 to 2 
               jjo = 0
               ln_dbg = .TRUE. 
            ELSE 
               jio = 0
               jjo = jr - 2
               ln_dbg = .FALSE.
            END IF 


            IF ( ln_hpg_ref ) THEN 
               IF ( ln_hpg_ref_str ) THEN 
                  IF ( ln_hpg_ref_ccs ) THEN 
                     IF ( ln_hpg_ref_off ) THEN 
                        !CALL ref_to_tgt_ccs_str_off ( jio, jjo, gde3w, rhd, zz_ref, zrhd_ref, zdrhd_k_ref, jk_bot_ref, z_rhd_pmr(:,:,:,jr) )
                        CALL ref_to_tgt_ccs_str_off ( jio, jjo, gdept(:,:,:,Kmm), rhd, zz_ref, zrhd_ref, zdrhd_k_ref, &
                         &                            jk_bot_ref, z_rhd_pmr(:,:,:,jr) )
                     ELSE
                        !CALL ref_to_tgt_ccs_str ( jio, jjo, gde3w, rhd, zz_ref, zrhd_ref, zdrhd_k_ref, jk_bot_ref, z_rhd_pmr(:,:,:,jr) )
                        CALL ref_to_tgt_ccs_str ( jio, jjo, gdept(:,:,:,Kmm), rhd, zz_ref, zrhd_ref, zdrhd_k_ref, &
                         &                        jk_bot_ref, z_rhd_pmr(:,:,:,jr) )
                     END IF 
                  ELSE  
                     !CALL ref_to_tgt_cub_str ( jio, jjo, gde3w, rhd, zz_ref, zrhd_ref, jk_bot_ref, z_rhd_pmr(:,:,:,jr) )
                     CALL ref_to_tgt_cub_str( jio, jjo, gdept(:,:,:,Kmm), rhd, zz_ref, zrhd_ref, &
                      &                       jk_bot_ref, z_rhd_pmr(:,:,:,jr) )
                     !CALL ref_to_tgt_cub_str_dbg( ln_dbg, jio, jjo, gdept(:,:,:,Kmm), rhd, zz_ref, zrhd_ref, &
                     ! &                           jk_bot_ref, z_rhd_pmr(:,:,:,jr) )
                  END IF
               ELSE         ! these calls are mainly retained to assist in testing
                  IF ( ln_hpg_ref_ccs ) THEN 
                     CALL ref_to_tgt_ccs ( jio, jjo, gde3w, rhd, zz_ref, zrhd_ref, zdrhd_k_ref, jk_bot_ref, z_rhd_pmr(:,:,:,jr) )
                  ELSE  
                     CALL ref_to_tgt_cub ( jio, jjo, gde3w, rhd, zz_ref, zrhd_ref, jk_bot_ref, z_rhd_pmr(:,:,:,jr) )  
                  END IF
               END IF
            END IF 

            IF ( ln_dbg_hpg ) CALL dbg_3dr( '3. z_rhd_pmr', z_rhd_pmr(:,:,:,jr) ) 

         END DO     

      !----------------------------------------------------------------------------------------
      ! 4. Calculations for side-face integrals
      !----------------------------------------------------------------------------------------
 	 
      !----------------------------------------------------------------------------------------
      !  4.1 compute and store elementary vertical differences then harmonic averages 
      !     based on z_rhd_pmr arrays (zdz_k has already been calculated)  
      !----------------------------------------------------------------------------------------

! start loop over the two side faces jr = 2 "left" face; jr = 3 "right" face 
         DO jr = 2, 3

            IF ( j_uv == 1 ) THEN  
               jio = jr - 2 
               jjo = 0 
            ELSE 
               jio = 0 
               jjo = jr - 2 
            END IF 

            DO_3D( 0, 0, 0, 0, 2, jpk )  
               ztmp(ji,jj,jk) =   z_rhd_pmr (ji,jj,jk,jr) - z_rhd_pmr(ji,jj,jk-1,jr)
            END_3D

            zdrhd_k_pmr(:,:,:) = 0._wp   ! should be unnecessary 

            DO_3D( 0, 0, 0, 0, 2, jpk-1 ) 
               cffw = MAX( 2._wp * ztmp(ji,jj,jk) * ztmp(ji,jj,jk+1), 0._wp )
               z1_cff = ztmp(ji,jj,jk) + ztmp(ji,jj,jk+1)
               zdrhd_k_pmr(ji,jj,jk) = cffw / SIGN( MAX( ABS(z1_cff), zep ), z1_cff )
            END_3D

! apply boundary conditions at top and bottom  
            DO_2D( 0, 0, 0, 0 )
               zdrhd_k_pmr(ji,jj,1) = aco_bc_rhd_srf * ( z_rhd_pmr(ji,jj,2,jr) - z_rhd_pmr(ji,jj,1,jr) ) - bco_bc_rhd_srf * zdrhd_k_pmr(ji,jj,2)
               iktb = mbkt(ji+jio,jj+jjo)
               IF ( iktb > 1 ) THEN
                  zdrhd_k_pmr(ji,jj,iktb) = aco_bc_rhd_bot * (z_rhd_pmr(ji,jj,iktb,jr) - z_rhd_pmr(ji,jj,iktb-1,jr) ) - bco_bc_rhd_bot * zdrhd_k_pmr(ji,jj,iktb-1)
               END IF
            END_2D


      !--------------------------------------------------------------
      ! 4.2 Upper half of top-most grid box, compute and store
      !-------------------------------------------------------------
!! ssh replaces e3w_n ; gde3w is a depth; the formulae involve heights  
!! rho_k stores grav * FX / rho_0  
!! *** AY note: ssh(ji,jj,Kmm) + gde3w(ji,jj,1) = e3w(ji,jj,1,Kmm)
            DO_2D( 0, 0, 0, 0)
               ztmp(ji,jj,1) =  grav * ( ssh(ji+jio,jj+jjo,Kmm) + gde3w(ji+jio,jj+jjo,1) )            & 
            &                           * (  z_rhd_pmr(ji,jj,1,jr)                                    &
            &                     + 0.5_wp * (   z_rhd_pmr(ji,jj,2,jr) - z_rhd_pmr(ji,jj,1,jr) )      &
            &                              * (   ssh   (ji+jio,jj+jjo,Kmm) + gde3w(ji+jio,jj+jjo,1) ) &
            &                              / ( - gde3w(ji+jio,jj+jjo,2) + gde3w(ji+jio,jj+jjo,1) )  )
            END_2D

      !--------------------------------------------------------------
      ! 4.3 Interior faces, compute and store
      !-------------------------------------------------------------

            DO_3D( 0, 0, 0, 0, 2, jpkm1 )
               ztmp(ji,jj,jk) = zcoef0 * (   z_rhd_pmr(ji,jj,jk,jr) + z_rhd_pmr(ji,jj,jk-1,jr) )                             &
            &                             * ( - gde3w(ji+jio,jj+jjo,jk) + gde3w(ji+jio,jj+jjo,jk-1) )                            &
            &                 + z_grav_10 * (                                                                                    &
            &     (   zdrhd_k_pmr(ji,jj,jk)    - zdrhd_k_pmr(ji,jj,jk-1) )                                                          &
            &   * ( - gde3w(ji+jio,jj+jjo,jk)  + gde3w(ji+jio,jj+jjo,jk-1)  - z1_12 * ( zdz_k(ji+jio,jj+jjo,jk) + zdz_k  (ji+jio,jj+jjo,jk-1) ) ) &
            &   - (   zdz_k(ji+jio,jj+jjo,jk)  - zdz_k(ji+jio,jj+jjo,jk-1) )                                                                      &
            &   * (   z_rhd_pmr(ji,jj,jk,jr) - z_rhd_pmr(ji,jj,jk-1,jr) - z1_12 * ( zdrhd_k_pmr(ji,jj,jk) + zdrhd_k_pmr(ji,jj,jk-1) ) )               &
            &                             )
            END_3D

! the force on the right face could be set equal to the average of the right face for this cell and the left face for the cell to the right
! this would require an lbc_lnk call  

! lmr stands for left minus right 

            IF ( jr == 2 ) THEN  
               DO_3D( 0, 0, 0, 0, 1, jpkm1 ) 
                  z_lmr_k(ji,jj,jk) =  ztmp(ji,jj,jk)   ! values on left face; 
               END_3D 
            ELSE 
               DO_3D( 0, 0, 0, 0, 1, jpkm1 ) 
                  z_lmr_k(ji,jj,jk) = z_lmr_k(ji,jj,jk) - ztmp(ji,jj,jk)   ! subtract the values on the right face 
               END_3D 
            END IF 

            IF ( ln_dbg_hpg ) CALL dbg_3dr( '4. z_lmr_k', z_lmr_k ) 

         END DO  ! jr 


      !----------------------------------------------------------------------------------------
      !  5. Calculations for upper and lower faces and the vertical integration 
      !----------------------------------------------------------------------------------------

         z_upp_ij(:,:) = 0._wp
         zhpi(:,:)     = 0._wp
         zhpj(:,:)     = 0._wp

         DO jk = 1, jpk -1 

            IF ( ln_dbg_hpg .AND. lwp ) THEN 
               WRITE(numout,*) 
               WRITE(numout,*) ' jk = ', jk 
            END IF 

      !----------------------------------------------------------------------------------------
      !  5.1 compute and store elementary horizontal differences zfor z_rhd_pmr arrays 
      !----------------------------------------------------------------------------------------

            DO_2D( 0, 0, 0, 0 )
               zdrhd_21(ji,jj) =   z_rhd_pmr(ji,jj,jk,2) - z_rhd_pmr(ji,jj,jk,1)
               zdrhd_32(ji,jj) =   z_rhd_pmr(ji,jj,jk,3) - z_rhd_pmr(ji,jj,jk,2)
               zdrhd_43(ji,jj) =   z_rhd_pmr(ji,jj,jk,4) - z_rhd_pmr(ji,jj,jk,3)
            END_2D

            DO_2D( 0, 0, 0, 0 )
               cff_31 = MAX( 2._wp * zdrhd_21(ji,jj) * zdrhd_32(ji,jj), 0._wp ) 
               z1_cff = zdrhd_21(ji,jj) + zdrhd_32(ji,jj)
               zz_drhd_ij(ji,jj,1) = cff_31 / SIGN( MAX( ABS(z1_cff), zep ), z1_cff )

               cff_42 = MAX( 2._wp * zdrhd_32(ji,jj) * zdrhd_43(ji,jj), 0._wp ) 
               z1_cff = zdrhd_32(ji,jj) + zdrhd_43(ji,jj)
               zz_drhd_ij(ji,jj,2) = cff_42 / SIGN( MAX( ABS(z1_cff), zep ), z1_cff ) 
            END_2D

      !----------------------------------------------------------------------------------
      ! 5.2 apply boundary conditions at side boundaries using 5.36-5.37
      !----------------------------------------------------------------------------------


! need to check this sub-section more carefully 

            zdrhd_ij(:,:,:) = zz_drhd_ij(:,:,:)

            IF ( j_uv == 1 ) THEN 

               DO_2D( 0, 0, 0, 0)       
            ! Walls coming from left: should check from 2 to jpi-1 (and jpj=2-jpj)
                  IF ( umask(ji,jj,jk) > 0.5_wp .AND. umask(ji-1,jj,jk) < 0.5_wp .AND. umask(ji+1,jj,jk) > 0.5_wp)  THEN  
                     zdrhd_ij(ji,jj,1) = aco_bc_rhd_hor * ( z_rhd_pmr(ji,jj,jk,3) - z_rhd_pmr(ji,jj,jk,2) ) - bco_bc_rhd_hor * zz_drhd_ij(ji,jj,2) 
                  END IF
            ! Walls coming from right: should check from 3 to jpi (and jpj=2-jpj)
                  IF ( umask(ji,jj,jk) < 0.5_wp .AND. umask(ji-1,jj,jk) > 0.5_wp .AND. umask(ji-2,jj,jk) > 0.5_wp) THEN
                     zdrhd_ij(ji,jj,2) = aco_bc_rhd_hor * ( z_rhd_pmr(ji,jj,jk,3) - z_rhd_pmr(ji,jj,jk,2) ) - bco_bc_rhd_hor * zz_drhd_ij(ji,jj,1)  
                  END IF
               END_2D 

            ELSE ! j_uv == 2 
	     
               DO_2D( 0, 0, 0, 0)    
            ! Walls coming from left: should check from 2 to jpj-1 (and jpi=2-jpi)
                  IF ( vmask(ji,jj,jk) > 0.5_wp .AND. vmask(ji,jj-1,jk) < 0.5_wp .AND. vmask(ji,jj+1,jk) > 0.5_wp)  THEN
                     zdrhd_ij(ji,jj,1) = aco_bc_rhd_hor * ( z_rhd_pmr(ji,jj,jk,3) - z_rhd_pmr(ji,jj,jk,2) ) - bco_bc_rhd_hor * zz_drhd_ij(ji,jj,2)
                  END IF
            ! Walls coming from right: should check from 3 to jpj (and jpi=2-jpi)
                  IF ( vmask(ji,jj,jk) < 0.5_wp .AND. vmask(ji,jj-1,jk) > 0.5_wp .AND. vmask(ji,jj-2,jk) > 0.5_wp) THEN 
                     zdrhd_ij(ji,jj,2) = aco_bc_rhd_hor * ( z_rhd_pmr(ji,jj,jk,3) - z_rhd_pmr(ji,jj,jk,2) ) - bco_bc_rhd_hor * zz_drhd_ij(ji,jj,1) 
                  END IF
               END_2D

            END IF ! j_uv == 2 

            IF ( ln_dbg_hpg ) THEN 
               CALL dbg_2dr( '5.2 zdrhd_ij(:,:,1)', zdrhd_ij(:,:,1) ) 
               CALL dbg_2dr( '5.2 zdrhd_ij(:,:,2)', zdrhd_ij(:,:,2) ) 
            END IF 	       

      !--------------------------------------------------------------
      ! 5.3 Calculate integrals on lower faces  
      !-------------------------------------------------------------

            IF ( j_uv == 1 ) THEN 

               DO_2D( 0, 0, 0, 0 )
! two -ve signs cancel in next two lines (within zcoef0 and because gde3w is a depth not a height)
                  z_low_ij(ji,jj) = zcoef0 * ( z_rhd_pmr(ji,jj,jk,3) + z_rhd_pmr(ji,jj,jk,2) )                                    &
             &                               * ( gde3w(ji+1,jj,jk) - gde3w(ji,jj,jk) )                                    

                  IF ( umask(ji-1, jj, jk) > 0.5 .OR. umask(ji+1, jj, jk) > 0.5 ) THEN
                     z_low_ij(ji,jj) = z_low_ij(ji,jj) - z_grav_10 * (                                                            &
             &           (   zdrhd_ij(ji,jj,2) - zdrhd_ij(ji,jj,1) )                                                              &
             &         * ( - gde3w(ji+1,jj,jk) + gde3w(ji,jj,jk) - z1_12 * ( zdz_i(ji+1,jj,jk) + zdz_i(ji,jj,jk) ) )              &
             &         - (   zdz_i(ji+1,jj,jk) - zdz_i(ji,jj,jk) )                                                                &
             &         * (   z_rhd_pmr(ji,jj,jk,3) - z_rhd_pmr(ji,jj,jk,2) - z1_12 * ( zdrhd_ij(ji,jj,2) + zdrhd_ij(ji,jj,1) ) )  &
             &                                               )
                  END IF
               END_2D

               IF ( ln_dbg_hpg ) CALL dbg_2dr( '5.3 z_low_ij 1', z_low_ij ) 
	       
	    ELSE ! j_uv == 2     
	         
               DO_2D( 0, 0, 0, 0 )
                  z_low_ij(ji,jj) = zcoef0 * ( z_rhd_pmr(ji,jj,jk,3) + z_rhd_pmr(ji,jj,jk,2) )                                    &
             &                    * ( gde3w(ji,jj+1,jk) - gde3w(ji,jj,jk) )                                  

                  IF ( vmask(ji, jj-1, jk) > 0.5 .OR. vmask(ji, jj+1, jk) > 0.5 ) THEN
                     z_low_ij(ji,jj) = z_low_ij(ji,jj) - z_grav_10 * (                                                            &
             &           (   zdrhd_ij(ji,jj,2) - zdrhd_ij(ji,jj,1) )                                                              &
             &         * ( - gde3w(ji,jj+1,jk) + gde3w(ji,jj,jk) - z1_12 * ( zdz_j(ji,jj+1,jk) + zdz_j(ji,jj,jk) ) )              &
             &         - (   zdz_j(ji,jj+1,jk) - zdz_j(ji,jj,jk) )                                                                &
             &         * (   z_rhd_pmr(ji,jj,jk,3) - z_rhd_pmr(ji,jj,jk,2) - z1_12 * ( zdrhd_ij(ji,jj,2) + zdrhd_ij(ji,jj,1) ) )  &
             &                                                 )
                  END IF
               END_2D

               IF ( ln_dbg_hpg ) CALL dbg_2dr( '5.3 z_low_ij 2', z_low_ij ) 
	       
	    END IF ! j_uv      

      !--------------------------------------------------------------
      ! 5.4 Integrate in the vertical (including contributions from both upper and lower and side-faces)    
      !-------------------------------------------------------------
      !
            IF ( j_uv == 1 ) THEN 

               DO_2D( 0, 0, 0, 0 )
                  zhpi(ji,jj) = zhpi(ji,jj) +                                        &
               &            ( z_lmr_k(ji,jj,jk)  - ( z_low_ij(ji,jj) - z_upp_ij(ji,jj) ) ) * r1_e1u(ji,jj)
         ! add to the general momentum trend
                  puu(ji,jj,jk,Krhs) = puu(ji,jj,jk,Krhs) + zhpi(ji,jj)
               END_2D

               IF ( ln_dbg_hpg ) CALL dbg_2dr( '5.4 zhpi', zhpi ) 

            ELSE ! j_uv == 2   

               DO_2D( 0, 0, 0, 0 )
                  zhpj(ji,jj) = zhpj(ji,jj) +                                         &
               &            ( z_lmr_k(ji,jj,jk) - ( z_low_ij(ji,jj) - z_upp_ij(ji,jj) ) ) * r1_e2v(ji,jj)
         ! add to the general momentum trend
                  pvv(ji,jj,jk,Krhs) = pvv(ji,jj,jk,Krhs) + zhpj(ji,jj)
               END_2D

               IF ( ln_dbg_hpg ) CALL dbg_2dr( '5.4 zhpj', zhpj ) 

            END IF ! j_uv  
      
            DO_2D( 0, 0, 0, 0 )
               z_upp_ij(ji,jj) = z_low_ij(ji,jj)
            END_2D 

         END DO ! k 

         ! temporary output of fields for debugging etc. 
         IF ( j_uv == 1) THEN
            CALL iom_put( "e3w", e3w(:,:,:,Kmm) )
            CALL iom_put( "rhd_hpg", z_rhd_pmr(:,:,:,1) )
            IF ( ln_hpg_ref ) THEN
               CALL iom_put( "rhd_ref", zrhd_ref )
               !CALL iom_put( "zz_ref", zz_ref)
               CALL iom_put( "jk_bot_ref", FLOAT( jk_bot_ref(:,:) ) )
            ENDIF
         END IF

      END DO ! j_uv  

   END SUBROUTINE hpg_djr
   
!-----------------------------------------------------------------------------------------------------------------
   
   SUBROUTINE ref_to_tgt_cub ( ki_off_tgt, kj_off_tgt, p_dep_tgt, p_fld_tgt, p_z_ref, p_fld_ref, kk_bot_ref, p_fld_tgt_ref) 
       
      INTEGER,                                INTENT(in)  :: ki_off_tgt    ! offset of points in target array in i-direction   
      INTEGER,                                INTENT(in)  :: kj_off_tgt    ! offset of points in target array in j-direction   
      REAL(wp), DIMENSION(A2D(nn_hls),jpk),   INTENT(in)  :: p_dep_tgt     ! depths of target profiles       
      REAL(wp), DIMENSION(A2D(nn_hls),jpk),   INTENT(in)  :: p_fld_tgt     ! field values on the target grid 
      REAL(wp), DIMENSION(A2D(nn_hls),jpk),   INTENT(in)  :: p_z_ref       ! heights of reference  profiles       
      REAL(wp), DIMENSION(A2D(nn_hls),jpk),   INTENT(in)  :: p_fld_ref     ! field values to be interpolated (in the vertical) on reference grid
      INTEGER,  DIMENSION(A2D(nn_hls)),       INTENT(in)  :: kk_bot_ref    ! bottom levels in the reference profile
      REAL(wp), DIMENSION(A2D(nn_hls),jpk),   INTENT(OUT) :: p_fld_tgt_ref ! target minus reference on target grid         

!---------------------------------------------------------------------------------------------------------------

      INTEGER,  DIMENSION(A2D(nn_hls),jpk) :: jk_ref_for_tgt   ! reference index for interpolation to target grid; target lies between jk_ref and jk_ref-1.
      LOGICAL,  DIMENSION(A2D(nn_hls),jpk) :: ll_tgt_off_cen   ! T => off-centred interpolation (values not used in this routine) 

!---------------------------------------------------------------------------------------------------------------

      INTEGER  :: ji, jj, jk                   ! loop indices 
      INTEGER  :: jkr
      REAL(wp) :: z_r_6, z_r_24                ! constants

      REAL(wp) :: zf_a, zf_b, zf_c, zf_d
      REAL(wp) :: zz_ref_jkr, zz_ref_jkrm1, zeta, zetasq 
      REAL(wp) :: zf_0, zf_1, zf_2, zf_3
      REAL(wp) :: zave_bc, zave_ad, zdif_cb, zdif_da
      
      REAL(wp) :: zz_tgt_lcl, zfld_ref_on_tgt 
      LOGICAL  :: ll_ccs  

!-----------------------------------------------------------------------------------------------------------------

      ll_ccs = .FALSE.      ! set for the call to loc_ref_tgt 
      z_r_6  = 1._wp / 6._wp 
      z_r_24 = 1._wp / 24._wp 

! find jk_ref_for_tgt (bounding levels on reference grid for each target point
      CALL loc_ref_tgt ( ll_ccs, ki_off_tgt, kj_off_tgt, p_dep_tgt, p_z_ref, kk_bot_ref, jk_ref_for_tgt, ll_tgt_off_cen )

      DO_3D( 0, 0, 0, 0, 1, jpk-1 ) 
         zz_tgt_lcl = - p_dep_tgt( ji+ki_off_tgt, jj+kj_off_tgt, jk )

! it would probably be better computationally for fld_ref to have the jk index first. 

!!! jkr >= 2 and p_fld_ref has jk = 0 available   
 
         jkr  = jk_ref_for_tgt(ji,jj,jk)	 
         zf_a = p_fld_ref(ji,jj,jkr-2)
	 zf_b = p_fld_ref(ji,jj,jkr-1)
	 zf_c = p_fld_ref(ji,jj,jkr  )
	 zf_d = p_fld_ref(ji,jj,jkr+1)

         zz_ref_jkrm1 = p_z_ref( ji, jj, jkr - 1 )   
         zz_ref_jkr = p_z_ref( ji, jj, jkr )
         zeta = ( zz_tgt_lcl - 0.5_wp*(zz_ref_jkr+zz_ref_jkrm1) ) / ( zz_ref_jkr - zz_ref_jkrm1 )  
         zetasq = zeta*zeta

         zave_bc = 0.5_wp*(zf_b+zf_c) 
         zave_ad = 0.5_wp*(zf_a+zf_d)
         zdif_cb = zf_c - zf_b
         zdif_da = zf_d - zf_a

         zf_0 = 1.125_wp*zave_bc - 0.125_wp*zave_ad
         zf_1 = 1.125_wp*zdif_cb - z_r_24*zdif_da 
         zf_2 = 0.5_wp*(zave_ad - zave_bc)             ! corrected 12/09/2021
         zf_3 = z_r_6 * zdif_da - 0.5_wp*zdif_cb       ! corrected 12/09/2021

         zfld_ref_on_tgt = zf_0 + zeta*zf_1 + zetasq*(zf_2 + zeta*zf_3) 
	 
! when zfld_ref_on_tgt is commented out in the next line, the results for hpg_djr should agree with those for hpg_djc.   
         p_fld_tgt_ref(ji, jj, jk) = p_fld_tgt(ji+ki_off_tgt, jj+kj_off_tgt, jk) - zfld_ref_on_tgt

      END_3D   

   END SUBROUTINE ref_to_tgt_cub

!---------------------------------------------------------------------------------------------------------------

   SUBROUTINE ref_to_tgt_cub_str ( ki_off_tgt, kj_off_tgt, p_dep_tgt, p_fld_tgt, p_z_ref, p_fld_ref, kk_bot_ref, p_fld_tgt_ref) 
 
      INTEGER,                                INTENT(in)  :: ki_off_tgt    ! offset of points in target array 
                                                                           ! in i-direction   
      INTEGER,                                INTENT(in)  :: kj_off_tgt    ! offset of points in target array 
                                                                           ! in j-direction   
      REAL(wp), DIMENSION(A2D(nn_hls),jpk),   INTENT(in)  :: p_dep_tgt     ! depths of target profiles       
      REAL(wp), DIMENSION(A2D(nn_hls),jpk),   INTENT(in)  :: p_fld_tgt     ! field values on the target grid 
      REAL(wp), DIMENSION(A2D(nn_hls),jpk),   INTENT(in)  :: p_z_ref       ! heights of reference  profiles       
      REAL(wp), DIMENSION(A2D(nn_hls),jpk),   INTENT(in)  :: p_fld_ref     ! field values to be interpolated 
                                                                           ! (in the vertical) on reference grid
      INTEGER,  DIMENSION(A2D(nn_hls)),       INTENT(in)  :: kk_bot_ref    ! bottom levels in the reference profile
      
      REAL(wp), DIMENSION(A2D(nn_hls),jpk),   INTENT(OUT) :: p_fld_tgt_ref ! target minus reference on target grid         

!---------------------------------------------------------------------------------------------------------------

      INTEGER,  DIMENSION(A2D(nn_hls),jpk) :: jk_ref_for_tgt   ! reference index for interpolation to target grid; 
                                                               ! target lies between jk_ref and jk_ref-1.
      LOGICAL,  DIMENSION(A2D(nn_hls),jpk) :: ll_tgt_off_cen   ! T => off-centred interpolation (values not used in 
                                                               !      this routine) 

!---------------------------------------------------------------------------------------------------------------

      INTEGER  :: ji, jj, jk                   ! loop indices 
      INTEGER  :: jkr

      REAL(wp) :: zf_p, zf_0, zf_m, zf_b   ! plus, zero, minus, below ; note that zeta increases with depth ; 
                                           ! zero is zeta = 1/2
      REAL(wp) :: zz_p, zz_0, zz_m, zz_b  
      REAL(wp) :: zs_p, zs_0, zs_m, zs_b
      REAL(wp) :: zz_p_m, zz_p_b, zz_m_b 
      REAL(wp) :: za_1, za_2, za_3   
      
      REAL(wp) :: zz_tgt_lcl, zz_tgt_sq, zfld_ref_on_tgt 
      LOGICAL  :: ll_ccs  

!-----------------------------------------------------------------------------------------------------------------

      ll_ccs = .FALSE.      ! set for the call to loc_ref_tgt 

      ! find jk_ref_for_tgt (bounding levels on reference grid for each target point
      CALL loc_ref_tgt ( ll_ccs, ki_off_tgt, kj_off_tgt, p_dep_tgt, &
       &                 p_z_ref, kk_bot_ref, jk_ref_for_tgt, ll_tgt_off_cen )

      DO_3D( 0, 0, 0, 0, 1, jpk-1 ) 

         ! it would probably be better computationally for fld_ref to have the jk index first. 

         !!! jkr >= 2 and p_fld_ref has jk = 0 available   
 
         jkr  = jk_ref_for_tgt(ji,jj,jk)
         zf_b = p_fld_ref(ji,jj,jkr-2)
         zf_m = p_fld_ref(ji,jj,jkr-1)
         zf_0 = p_fld_ref(ji,jj,jkr  )
         zf_p = p_fld_ref(ji,jj,jkr+1)

         zz_b = p_z_ref( ji, jj, jkr-2) 
         zz_m = p_z_ref( ji, jj, jkr-1)
         zz_0 = p_z_ref( ji, jj, jkr  )   
         zz_b = zz_b - zz_0 
         zz_m = zz_m - zz_0  
         zz_p = p_z_ref( ji, jj, jkr+1) - zz_0

         zz_p_m = zz_p - zz_m 
         zz_p_b = zz_p - zz_b 
         zz_m_b = zz_m - zz_b 

         zs_0 = - zf_0 / ( zz_p * zz_m   * zz_b   )    

         zs_p =   zf_p / ( zz_p * zz_p_m * zz_p_b ) 
         zs_b =   zf_b / ( zz_b * zz_p_b * zz_m_b ) 
         zs_m = - zf_m / ( zz_m * zz_p_m * zz_m_b ) 

         za_1 =    zz_m*zz_b *zs_p +  zz_p*zz_b *zs_m +  zz_p*zz_m *zs_b + (zz_p*zz_m+zz_m*zz_b+zz_p*zz_b)*zs_0  
         za_2 = - (zz_m+zz_b)*zs_p - (zz_p+zz_b)*zs_m - (zz_p+zz_m)*zs_b - (zz_p+zz_m+zz_b)*zs_0 
         za_3 = zs_p + zs_m + zs_0 + zs_b

         zz_tgt_lcl = - p_dep_tgt( ji+ki_off_tgt, jj+kj_off_tgt, jk ) - zz_0
         zz_tgt_sq = zz_tgt_lcl*zz_tgt_lcl
    
         zfld_ref_on_tgt = zf_0 + za_1*zz_tgt_lcl + za_2*zz_tgt_sq + za_3*zz_tgt_lcl*zz_tgt_sq
    
         ! when zfld_ref_on_tgt is commented out in the next line, the results for hpg_djr 
         ! should agree with those for hpg_djc.
         p_fld_tgt_ref(ji, jj, jk) = p_fld_tgt(ji+ki_off_tgt, jj+kj_off_tgt, jk) - zfld_ref_on_tgt

      END_3D   

   END SUBROUTINE ref_to_tgt_cub_str

!-----------------------------------------------------------------------------------------------------------------

   SUBROUTINE ref_to_tgt_cub_str_dbg ( ln_debug, ki_off_tgt, kj_off_tgt, p_dep_tgt, p_fld_tgt, p_z_ref, p_fld_ref, kk_bot_ref, p_fld_tgt_ref)

      LOGICAL,                                INTENT(in)  :: ln_debug      ! DB for debugging
      INTEGER,                                INTENT(in)  :: ki_off_tgt    ! offset of points in target array 
                                                                           ! in i-direction   
      INTEGER,                                INTENT(in)  :: kj_off_tgt    ! offset of points in target array 
                                                                           ! in j-direction   
      REAL(wp), DIMENSION(A2D(nn_hls),jpk),   INTENT(in)  :: p_dep_tgt     ! depths of target profiles       
      REAL(wp), DIMENSION(A2D(nn_hls),jpk),   INTENT(in)  :: p_fld_tgt     ! field values on the target grid 
      REAL(wp), DIMENSION(A2D(nn_hls),jpk),   INTENT(in)  :: p_z_ref       ! heights of reference  profiles       
      REAL(wp), DIMENSION(A2D(nn_hls),jpk),   INTENT(in)  :: p_fld_ref     ! field values to be interpolated 
                                                                           ! (in the vertical) on reference grid
      INTEGER,  DIMENSION(A2D(nn_hls)),       INTENT(in)  :: kk_bot_ref    ! bottom levels in the reference profile

      REAL(wp), DIMENSION(A2D(nn_hls),jpk),   INTENT(OUT) :: p_fld_tgt_ref ! target minus reference on target grid         

!---------------------------------------------------------------------------------------------------------------

      INTEGER,  DIMENSION(A2D(nn_hls),jpk) :: jk_ref_for_tgt   ! reference index for interpolation to target grid; 
                                                               ! target lies between jk_ref and jk_ref-1.
      LOGICAL,  DIMENSION(A2D(nn_hls),jpk) :: ll_tgt_off_cen   ! T => off-centred interpolation (values not used in 
                                                               !      this routine) 
      REAL(wp), DIMENSION(A2D(nn_hls),jpk) :: zfld_ref_on_tgt    

!---------------------------------------------------------------------------------------------------------------

      INTEGER  :: ji, jj, jk                   ! loop indices 
      INTEGER  :: jkr

      REAL(wp) :: zf_p, zf_0, zf_m, zf_b   ! plus, zero, minus, below ; note that zeta increases with depth ; 
                                           ! zero is zeta = 1/2
      REAL(wp) :: zz_p, zz_0, zz_m, zz_b
      REAL(wp) :: zs_p, zs_0, zs_m, zs_b
      REAL(wp) :: zz_p_m, zz_p_b, zz_m_b
      REAL(wp) :: za_1, za_2, za_3

      REAL(wp) :: zz_tgt_lcl, zz_tgt_sq!, zfld_ref_on_tgt
      LOGICAL  :: ll_ccs

!-----------------------------------------------------------------------------------------------------------------

      ll_ccs = .FALSE.      ! set for the call to loc_ref_tgt 

      ! find jk_ref_for_tgt (bounding levels on reference grid for each target point
      CALL loc_ref_tgt ( ll_ccs, ki_off_tgt, kj_off_tgt, p_dep_tgt, &
       &                 p_z_ref, kk_bot_ref, jk_ref_for_tgt, ll_tgt_off_cen )

      ! DB for debugging
      IF ( ln_debug .AND. ki_off_tgt == 0 ) THEN
         CALL iom_put( "jk_ref_for_tgt", FLOAT( jk_ref_for_tgt(:,:,:) ) )
      END IF

      DO_3D( 0, 0, 0, 0, 1, jpk-1 ) 

         ! it would probably be better computationally for fld_ref to have the jk index first. 

         !!! jkr >= 2 and p_fld_ref has jk = 0 available   

         jkr  = jk_ref_for_tgt(ji,jj,jk)
         zf_b = p_fld_ref(ji,jj,jkr-2)
         zf_m = p_fld_ref(ji,jj,jkr-1)
         zf_0 = p_fld_ref(ji,jj,jkr  )
         zf_p = p_fld_ref(ji,jj,jkr+1)

         zz_b = p_z_ref( ji, jj, jkr-2)
         zz_m = p_z_ref( ji, jj, jkr-1)
         zz_0 = p_z_ref( ji, jj, jkr  )
         zz_b = zz_b - zz_0
         zz_m = zz_m - zz_0
         zz_p = p_z_ref( ji, jj, jkr+1) - zz_0

         zz_p_m = zz_p - zz_m
         zz_p_b = zz_p - zz_b
         zz_m_b = zz_m - zz_b

         zs_0 = - zf_0 / ( zz_p * zz_m   * zz_b   )

         zs_p =   zf_p / ( zz_p * zz_p_m * zz_p_b )
         zs_b =   zf_b / ( zz_b * zz_p_b * zz_m_b )
         zs_m = - zf_m / ( zz_m * zz_p_m * zz_m_b )

         za_1 =    zz_m*zz_b *zs_p +  zz_p*zz_b *zs_m +  zz_p*zz_m *zs_b + (zz_p*zz_m+zz_m*zz_b+zz_p*zz_b)*zs_0
         za_2 = - (zz_m+zz_b)*zs_p - (zz_p+zz_b)*zs_m - (zz_p+zz_m)*zs_b - (zz_p+zz_m+zz_b)*zs_0
         za_3 = zs_p + zs_m + zs_0 + zs_b

         zz_tgt_lcl = - p_dep_tgt( ji+ki_off_tgt, jj+kj_off_tgt, jk ) - zz_0
         zz_tgt_sq = zz_tgt_lcl*zz_tgt_lcl

         !zfld_ref_on_tgt = zf_0 + za_1*zz_tgt_lcl + za_2*zz_tgt_sq + za_3*zz_tgt_lcl*zz_tgt_sq
         zfld_ref_on_tgt(ji, jj, jk) = (zf_0 + za_1*zz_tgt_lcl + za_2*zz_tgt_sq + za_3*zz_tgt_lcl*zz_tgt_sq) * &
           &                           tmask(ji+ki_off_tgt, jj+kj_off_tgt, jk)


         ! when zfld_ref_on_tgt is commented out in the next line, the results for hpg_djr 
         ! should agree with those for hpg_djc.
         p_fld_tgt_ref(ji, jj, jk) = ( p_fld_tgt(ji+ki_off_tgt, jj+kj_off_tgt, jk) - zfld_ref_on_tgt(ji, jj, jk) ) * &
           &                         tmask(ji+ki_off_tgt, jj+kj_off_tgt, jk)

      END_3D

      ! DB for debugging
      IF ( ln_debug .AND. ki_off_tgt == 0 ) THEN
         CALL iom_put( "zfld_ref_on_tgt", zfld_ref_on_tgt(:,:,:) )
         CALL iom_put( "p_fld_tgt", p_fld_tgt(:,:,:) ) 
      END IF  

   END SUBROUTINE ref_to_tgt_cub_str_dbg

!-----------------------------------------------------------------------------------------------------------------

   SUBROUTINE ref_to_tgt_ccs ( ki_off_tgt, kj_off_tgt, pdep_tgt, pfld_tgt, pz_ref, pfld_ref, pdfld_k_ref, kk_bot_ref, pfld_tgt_ref) 
   
      INTEGER,                      INTENT(in)  :: ki_off_tgt    ! offset of points in target array in i-direction   
      INTEGER,                      INTENT(in)  :: kj_off_tgt    ! offset of points in target array in j-direction   
      REAL(wp), DIMENSION(A2D(nn_hls),jpk),   INTENT(in)  :: pdep_tgt      ! depths of target profiles       
      REAL(wp), DIMENSION(A2D(nn_hls),jpk),   INTENT(in)  :: pfld_tgt      ! field values on the target grid 
      REAL(wp), DIMENSION(A2D(nn_hls),jpk),   INTENT(in)  :: pz_ref        ! heights of reference  profiles       
      REAL(wp), DIMENSION(A2D(nn_hls),jpk),   INTENT(in)  :: pfld_ref      ! reference field values
      REAL(wp), DIMENSION(A2D(nn_hls),jpk),   INTENT(in)  :: pdfld_k_ref   ! harmonic averages of vertical differences of reference field
      INTEGER,  DIMENSION(A2D(nn_hls)),       INTENT(in)  :: kk_bot_ref    ! bottom levels in the reference profile
      REAL(wp), DIMENSION(A2D(nn_hls),jpk),   INTENT(OUT) :: pfld_tgt_ref  ! target minus reference on target grid         

!-----------------------------------------------------------------------------------------------------------------
      INTEGER,  DIMENSION(A2D(nn_hls),jpk) :: jk_ref_for_tgt   ! reference index for interpolation to target grid
      LOGICAL,  DIMENSION(A2D(nn_hls),jpk) :: ll_tgt_off_cen   ! T => off-centred interpolation (values not used in this routine) 

      INTEGER  :: ji, jj, jk                 ! DO loop indices  
      INTEGER  :: jkr
      REAL(wp) :: zz_tgt_lcl, zz_ref_jkrm1, zz_ref_jkr, zeta, zetasq
      REAL(wp) :: zd_dif, zd_ave, zf_dif, zf_ave
      REAL(wp) :: zcoef_0, zcoef_1, zcoef_2, zcoef_3    
      REAL(wp) :: zfld_ref_on_tgt
      LOGICAL  :: ll_ccs                     ! set to .TRUE. for the call to loc_ref_tgt  

!-----------------------------------------------------------------------------------------------------------------

      ll_ccs = .TRUE.      ! set for the call to loc_ref_tgt 

      ! find jk_ref_for_tgt (bounding levels on reference grid for each target point
      CALL loc_ref_tgt ( ll_ccs, ki_off_tgt, kj_off_tgt, pdep_tgt, pz_ref, kk_bot_ref, jk_ref_for_tgt, ll_tgt_off_cen )

      DO_3D( 0, 0, 0, 0, 1, jpk-1 ) 
         zz_tgt_lcl =  - pdep_tgt( ji+ki_off_tgt, jj+kj_off_tgt, jk )
         jkr = jk_ref_for_tgt( ji, jj, jk )
         zz_ref_jkrm1 = pz_ref( ji, jj, jkr - 1 )   
         zz_ref_jkr = pz_ref( ji, jj, jkr )
         zeta = ( zz_tgt_lcl - 0.5_wp*(zz_ref_jkr+zz_ref_jkrm1) ) / ( zz_ref_jkr - zz_ref_jkrm1 )  
         zetasq = zeta*zeta

         zd_dif =            pdfld_k_ref(ji,jj,jkr) - pdfld_k_ref(ji,jj,jkr-1)  
         zd_ave = 0.5_wp * ( pdfld_k_ref(ji,jj,jkr) + pdfld_k_ref(ji,jj,jkr-1) )   

         zf_dif =            pfld_ref(ji,jj,jkr) - pfld_ref(ji,jj,jkr-1)  
         zf_ave = 0.5_wp * ( pfld_ref(ji,jj,jkr) + pfld_ref(ji,jj,jkr-1) )   

         zcoef_0 =          zf_ave - 0.125_wp * zd_dif
         zcoef_1 = 1.5_wp * zf_dif - 0.5_wp   * zd_ave
         zcoef_2 = 0.5_wp * zd_dif 
         zcoef_3 = 2.0_wp * zd_ave - 2._wp    * zf_dif   
 
         zfld_ref_on_tgt = zcoef_0 + zeta*zcoef_1 + zetasq * ( zcoef_2 + zeta*zcoef_3 )  

         ! when zfld_ref_on_tgt is commented out in the next line, 
         ! the results for hpg_djr should agree with those for hpg_djc.   
 
         pfld_tgt_ref(ji, jj, jk) = pfld_tgt(ji+ki_off_tgt, jj+kj_off_tgt, jk) - zfld_ref_on_tgt  

!         IF ( ln_dbg_hpg .AND. lwp .AND. ji == ki_dbg_cen .AND. jj == kj_dbg_cen .AND. jk == kk_dbg_cen ) THEN
!            WRITE(numout,*) ' zeta, zd_dif, zd_ave, zf_dif, zf_ave = ', zeta, zd_dif, zd_ave, zf_dif, zf_ave
!            WRITE(numout,*) ' zz_tgt_lcl, jkr, zz_ref_jkr, zz_ref_jkrm1 =', zz_tgt_lcl, jkr, zz_ref_jkr, zz_ref_jkrm1 
!            WRITE(numout,*) ' pfld_ref(ji,jj,jkr), pfld_ref(ji,jj,jkr-1) = ', pfld_ref(ji,jj,jkr), pfld_ref(ji,jj,jkr-1)
!	    WRITE(numout,*) ' zfld_ref_on_tgt = ', zfld_ref_on_tgt 
!	    WRITE(numout,*) ' pfld_tgt(ji+ki_off_tgt, jj+kj_off_tgt, jk) = ', pfld_tgt(ji+ki_off_tgt, jj+kj_off_tgt, jk)
!         END IF             

      END_3D   

      IF ( ln_dbg_hpg ) CALL dbg_3dr( 'ref_to_tgt_ccs: pfld_tgt_ref', pfld_tgt_ref ) 

   END SUBROUTINE ref_to_tgt_ccs

!-----------------------------------------------------------------------------------------------------------------

   SUBROUTINE ref_to_tgt_ccs_str ( ki_off_tgt, kj_off_tgt, pdep_tgt, pfld_tgt, pz_ref, pfld_ref, pdfld_k_ref, kk_bot_ref, pfld_tgt_ref) 
       
      INTEGER,                      INTENT(in)  :: ki_off_tgt    ! offset of points in target array in i-direction   
      INTEGER,                      INTENT(in)  :: kj_off_tgt    ! offset of points in target array in j-direction   
      REAL(wp), DIMENSION(A2D(nn_hls),jpk),   INTENT(in)  :: pdep_tgt      ! depths of target profiles       
      REAL(wp), DIMENSION(A2D(nn_hls),jpk),   INTENT(in)  :: pfld_tgt      ! field values on the target grid 
      REAL(wp), DIMENSION(A2D(nn_hls),jpk),   INTENT(in)  :: pz_ref        ! heights of reference  profiles       
      REAL(wp), DIMENSION(A2D(nn_hls),jpk),   INTENT(in)  :: pfld_ref      ! reference field values
      REAL(wp), DIMENSION(A2D(nn_hls),jpk),   INTENT(in)  :: pdfld_k_ref   ! harmonic averages of vertical differences of reference field
      INTEGER,  DIMENSION(A2D(nn_hls)),       INTENT(in)  :: kk_bot_ref    ! bottom levels in the reference profile
      REAL(wp), DIMENSION(A2D(nn_hls),jpk),   INTENT(OUT) :: pfld_tgt_ref  ! target minus reference on target grid         

!-----------------------------------------------------------------------------------------------------------------
      INTEGER,  DIMENSION(A2D(nn_hls),jpk) :: jk_ref_for_tgt   ! reference index for interpolation to target grid
      LOGICAL,  DIMENSION(A2D(nn_hls),jpk) :: ll_tgt_off_cen   ! T => off-centred interpolation (values not used in this routine) 

      INTEGER  :: ji, jj, jk                 ! DO loop indices  
      INTEGER  :: jkr
      REAL(wp) :: zz_tgt_lcl, zz_tgt_sq

      REAL(wp) :: zeta_p, zeta_0, zeta_m, zeta_b   ! plus, zero, minus, below (note that zeta increases with depth) zero is zeta = 1/2
      REAL(wp) :: zz_p, zz_0, zz_m, zz_b  
      REAL(wp) :: zs_p, zs_0, zs_m, zs_b
      REAL(wp) :: zz_p_m, zz_p_b, zz_m_b 
      REAL(wp) :: za_1, za_2, za_3   
      REAL(wp) :: zeta, zetasq

      REAL(wp) :: zd_dif, zd_ave, zf_dif, zf_ave
      REAL(wp) :: zcoef_0, zcoef_1, zcoef_2, zcoef_3    
      REAL(wp) :: zfld_ref_on_tgt
      LOGICAL  :: ll_ccs                     ! set to .FALSE. for the call to loc_ref_tgt  

!-----------------------------------------------------------------------------------------------------------------

      ll_ccs = .TRUE.      ! set for the call to loc_ref_tgt 

! find jk_ref_for_tgt (bounding levels on reference grid for each target point
      CALL loc_ref_tgt ( ll_ccs, ki_off_tgt, kj_off_tgt, pdep_tgt, pz_ref, kk_bot_ref, jk_ref_for_tgt, ll_tgt_off_cen )

      zeta_p =  1.5_wp
      zeta_0 =  0.5_wp
      zeta_m = -0.5_wp
      zeta_b = -1.5_wp

      DO_3D( 0, 0, 0, 0, 1, jpk-1 ) 

         jkr = jk_ref_for_tgt( ji, jj, jk )

         zz_0 = pz_ref( ji, jj, jkr  )   
         zz_m = pz_ref( ji, jj, jkr-1) - zz_0

         IF ( jkr .NE. 2) THEN 
            zz_b = pz_ref( ji, jj, jkr-2) - zz_0  
         ELSE 
            zz_b = 2._wp * zz_m       
         END IF 

         IF ( ll_tgt_off_cen(ji,jj,jk) ) THEN 
            zz_p = - zz_m            
         ELSE  
            zz_p = pz_ref( ji, jj, jkr+1) - zz_0
         END IF 

         zz_p_m = zz_p - zz_m 
         zz_p_b = zz_p - zz_b 
         zz_m_b = zz_m - zz_b 

         zs_0 = - zeta_0 / ( zz_p * zz_m   * zz_b   )    

         zs_p =   zeta_p / ( zz_p * zz_p_m * zz_p_b ) 
         zs_m = - zeta_m / ( zz_m * zz_p_m * zz_m_b ) 
         zs_b =   zeta_b / ( zz_b * zz_p_b * zz_m_b ) 

         za_1 =     zz_m*zz_b*zs_p +   zz_p*zz_b*zs_m +   zz_p*zz_m*zs_b + (zz_p*zz_m+zz_m*zz_b+zz_p*zz_b)*zs_0  
         za_2 = - (zz_m+zz_b)*zs_p - (zz_p+zz_b)*zs_m - (zz_p+zz_m)*zs_b - (zz_p+zz_m+zz_b)*zs_0 
         za_3 = zs_p + zs_m + zs_0 + zs_b

         zz_tgt_lcl =  - pdep_tgt( ji+ki_off_tgt, jj+kj_off_tgt, jk ) - zz_0
         zz_tgt_sq = zz_tgt_lcl*zz_tgt_lcl
    
         zeta = zeta_0 + za_1*zz_tgt_lcl + za_2*zz_tgt_sq + za_3*zz_tgt_lcl*zz_tgt_sq 
         zetasq = zeta*zeta

         zd_dif =            pdfld_k_ref(ji,jj,jkr) - pdfld_k_ref(ji,jj,jkr-1)  
         zd_ave = 0.5_wp * ( pdfld_k_ref(ji,jj,jkr) + pdfld_k_ref(ji,jj,jkr-1) )   

         zf_dif =            pfld_ref(ji,jj,jkr) - pfld_ref(ji,jj,jkr-1)  
         zf_ave = 0.5_wp * ( pfld_ref(ji,jj,jkr) + pfld_ref(ji,jj,jkr-1) )   

         zcoef_0 =          zf_ave - 0.125_wp * zd_dif
         zcoef_1 = 1.5_wp * zf_dif - 0.5_wp   * zd_ave
         zcoef_2 = 0.5_wp * zd_dif 
         zcoef_3 = 2.0_wp * zd_ave - 2._wp    * zf_dif   
 
         zfld_ref_on_tgt = zcoef_0 + zeta*zcoef_1 + zetasq * ( zcoef_2 + zeta*zcoef_3 )  

         ! when zfld_ref_on_tgt is commented out in the next line, the results for hpg_djr 
         ! should agree with those for hpg_djc.   
 
         pfld_tgt_ref(ji, jj, jk) = pfld_tgt(ji+ki_off_tgt, jj+kj_off_tgt, jk) - zfld_ref_on_tgt  

!         IF ( ln_dbg_hpg .AND. lwp .AND. ji == ki_dbg_cen .AND. jj == kj_dbg_cen .AND. jk == kk_dbg_cen ) THEN
!            WRITE(numout,*) ' zeta, zd_dif, zd_ave, zf_dif, zf_ave = ', zeta, zd_dif, zd_ave, zf_dif, zf_ave
!            WRITE(numout,*) ' zz_tgt_lcl, jkr, zz_ref_jkr, zz_ref_jkrm1 =', zz_tgt_lcl, jkr, zz_ref_jkr, zz_ref_jkrm1 
!            WRITE(numout,*) ' pfld_ref(ji,jj,jkr), pfld_ref(ji,jj,jkr-1) = ', pfld_ref(ji,jj,jkr), pfld_ref(ji,jj,jkr-1)
!	    WRITE(numout,*) ' zfld_ref_on_tgt = ', zfld_ref_on_tgt 
!	    WRITE(numout,*) ' pfld_tgt(ji+ki_off_tgt, jj+kj_off_tgt, jk) = ', pfld_tgt(ji+ki_off_tgt, jj+kj_off_tgt, jk)
!         END IF             
      
      END_3D   

      IF ( ln_dbg_hpg ) CALL dbg_3dr( 'ref_to_tgt_ccs_str: pfld_tgt_ref', pfld_tgt_ref ) 

   END SUBROUTINE ref_to_tgt_ccs_str

!-----------------------------------------------------------------------------------------------------------------

   SUBROUTINE ref_to_tgt_ccs_str_off ( ki_off_tgt, kj_off_tgt, pdep_tgt, pfld_tgt, pz_ref, pfld_ref, pdfld_k_ref, kk_bot_ref, pfld_tgt_ref) 

! Coded for a stretched grid and to perform off-centred pure cubic interpolation at the upper and lower boundaries 
       
      INTEGER,                      INTENT(in)  :: ki_off_tgt    ! offset of points in target array in i-direction   
      INTEGER,                      INTENT(in)  :: kj_off_tgt    ! offset of points in target array in j-direction   
      REAL(wp), DIMENSION(A2D(nn_hls),jpk),   INTENT(in)  :: pdep_tgt      ! depths of target profiles       
      REAL(wp), DIMENSION(A2D(nn_hls),jpk),   INTENT(in)  :: pfld_tgt      ! field values on the target grid 
      REAL(wp), DIMENSION(A2D(nn_hls),jpk),   INTENT(in)  :: pz_ref        ! heights of reference  profiles       
      REAL(wp), DIMENSION(A2D(nn_hls),jpk),   INTENT(in)  :: pfld_ref      ! reference field values
      REAL(wp), DIMENSION(A2D(nn_hls),jpk),   INTENT(in)  :: pdfld_k_ref   ! harmonic averages of vertical differences of reference field
      INTEGER,  DIMENSION(A2D(nn_hls)),       INTENT(in)  :: kk_bot_ref    ! bottom levels in the reference profile
      REAL(wp), DIMENSION(A2D(nn_hls),jpk),   INTENT(OUT) :: pfld_tgt_ref  ! target minus reference on target grid         

!-----------------------------------------------------------------------------------------------------------------
      INTEGER,  DIMENSION(A2D(nn_hls),jpk) :: jk_ref_for_tgt   ! reference index for interpolation to target grid
      LOGICAL,  DIMENSION(A2D(nn_hls),jpk) :: ll_tgt_off_cen   ! T => off-centred interpolation

      INTEGER  :: ji, jj, jk                 ! DO loop indices  
      INTEGER  :: jkr
      REAL(wp) :: zz_tgt_lcl, zz_tgt_sq

      REAL(wp) :: zz_p, zz_0, zz_m, zz_b  
      REAL(wp) :: zeta_p, zeta_0, zeta_m, zeta_b
      REAL(wp) :: zf_p, zf_0, zf_m, zf_b
      REAL(wp) :: zs_p, zs_0, zs_m, zs_b
      REAL(wp) :: zz_p_m, zz_p_b, zz_m_b 
      REAL(wp) :: za_1, za_2,za_3   
      REAL(wp) :: zeta, zetasq

      REAL(wp) :: zd_dif, zd_ave, zf_dif, zf_ave
      REAL(wp) :: zcoef_0, zcoef_1, zcoef_2, zcoef_3    
      REAL(wp) :: zfld_ref_on_tgt
      LOGICAL  :: ll_ccs                     ! set to .FALSE. for the call to loc_ref_tgt 

!-----------------------------------------------------------------------------------------------------------------

      ll_ccs = .FALSE.      ! set for the call to loc_ref_tgt
      
      ! find jk_ref_for_tgt (bounding levels on reference grid for each target point) and ll_tgt_off_cen
      CALL loc_ref_tgt ( ll_ccs, ki_off_tgt, kj_off_tgt, pdep_tgt, pz_ref, kk_bot_ref, jk_ref_for_tgt, ll_tgt_off_cen )

      zeta_p =  1.5_wp
      zeta_0 =  0.5_wp
      zeta_m = -0.5_wp
      zeta_b = -1.5_wp

      DO_3D( 0, 0, 0, 0, 1, jpk-1 )  

         jkr = jk_ref_for_tgt( ji, jj, jk )

         zz_b = pz_ref( ji, jj, jkr-2) 
         zz_m = pz_ref( ji, jj, jkr-1)
         zz_0 = pz_ref( ji, jj, jkr  )   
         zz_b = zz_b - zz_0 
         zz_m = zz_m - zz_0  
         zz_p = pz_ref( ji, jj, jkr+1) - zz_0

         zz_tgt_lcl =  - pdep_tgt( ji+ki_off_tgt, jj+kj_off_tgt, jk ) - zz_0
         zz_tgt_sq = zz_tgt_lcl*zz_tgt_lcl

         IF ( ll_tgt_off_cen(ji,jj,jk) ) THEN 

            ! off centred pure cubic fit 

            zf_b = pfld_ref(ji,jj,jkr-2)
            zf_m = pfld_ref(ji,jj,jkr-1)
            zf_0 = pfld_ref(ji,jj,jkr  )
            zf_p = pfld_ref(ji,jj,jkr+1)

            zz_p_m = zz_p - zz_m 
            zz_p_b = zz_p - zz_b 
            zz_m_b = zz_m - zz_b 

            zs_0 = - zf_0 / ( zz_p * zz_m   * zz_b   )    

            zs_p =   zf_p / ( zz_p * zz_p_m * zz_p_b ) 
            zs_b =   zf_b / ( zz_b * zz_p_b * zz_m_b ) 
            zs_m = - zf_m / ( zz_m * zz_p_m * zz_m_b ) 

            za_1 =    zz_m*zz_b *zs_p +  zz_p*zz_b *zs_m +  zz_p*zz_m *zs_b + (zz_p*zz_m+zz_m*zz_b+zz_p*zz_b)*zs_0  
            za_2 = - (zz_m+zz_b)*zs_p - (zz_p+zz_b)*zs_m - (zz_p+zz_m)*zs_b - (zz_p+zz_m+zz_b)*zs_0 
            za_3 = zs_p + zs_m + zs_0 + zs_b

            zfld_ref_on_tgt = zf_0 + za_1*zz_tgt_lcl + za_2*zz_tgt_sq + za_3*zz_tgt_lcl*zz_tgt_sq 

         ELSE 

            zz_p_m = zz_p - zz_m 
            zz_p_b = zz_p - zz_b 
            zz_m_b = zz_m - zz_b 

            zs_0 = - zeta_0 / ( zz_p * zz_m   * zz_b   )    

            zs_p =   zeta_p / ( zz_p * zz_p_m * zz_p_b ) 
            zs_m = - zeta_m / ( zz_m * zz_p_m * zz_m_b ) 
            zs_b =   zeta_b / ( zz_b * zz_p_b * zz_m_b ) 

            za_1 =     zz_m*zz_b*zs_p +   zz_p*zz_b*zs_m +   zz_p*zz_m*zs_b + (zz_p*zz_m+zz_m*zz_b+zz_p*zz_b)*zs_0  
            za_2 = - (zz_m+zz_b)*zs_p - (zz_p+zz_b)*zs_m - (zz_p+zz_m)*zs_b - (zz_p+zz_m+zz_b)*zs_0 
            za_3 = zs_p + zs_m + zs_0 + zs_b
    
            zeta = zeta_0 + za_1*zz_tgt_lcl + za_2*zz_tgt_sq + za_3*zz_tgt_lcl*zz_tgt_sq 
            zetasq = zeta*zeta

            zd_dif =            pdfld_k_ref(ji,jj,jkr) - pdfld_k_ref(ji,jj,jkr-1)  
            zd_ave = 0.5_wp * ( pdfld_k_ref(ji,jj,jkr) + pdfld_k_ref(ji,jj,jkr-1) )   

            zf_dif =            pfld_ref(ji,jj,jkr) - pfld_ref(ji,jj,jkr-1)  
            zf_ave = 0.5_wp * ( pfld_ref(ji,jj,jkr) + pfld_ref(ji,jj,jkr-1) )   

            zcoef_0 =          zf_ave - 0.125_wp * zd_dif
            zcoef_1 = 1.5_wp * zf_dif - 0.5_wp   * zd_ave
            zcoef_2 = 0.5_wp * zd_dif 
            zcoef_3 = 2.0_wp * zd_ave - 2._wp    * zf_dif   
 
            zfld_ref_on_tgt = zcoef_0 + zeta*zcoef_1 + zetasq * ( zcoef_2 + zeta*zcoef_3 )  

         END IF 

         ! when zfld_ref_on_tgt is commented out in the next line, the results for hpg_djr should agree with those for hpg_djc.   
 
         pfld_tgt_ref(ji, jj, jk) = pfld_tgt(ji+ki_off_tgt, jj+kj_off_tgt, jk) - zfld_ref_on_tgt  
    
!         IF ( ln_dbg_hpg .AND. lwp .AND. ji == ki_dbg_cen .AND. jj == kj_dbg_cen .AND. jk == kk_dbg_cen ) THEN
!            WRITE(numout,*) ' zeta, zd_dif, zd_ave, zf_dif, zf_ave = ', zeta, zd_dif, zd_ave, zf_dif, zf_ave
!            WRITE(numout,*) ' zz_tgt_lcl, jkr, zz_ref_jkr, zz_ref_jkrm1 =', zz_tgt_lcl, jkr, zz_ref_jkr, zz_ref_jkrm1 
!            WRITE(numout,*) ' pfld_ref(ji,jj,jkr), pfld_ref(ji,jj,jkr-1) = ', pfld_ref(ji,jj,jkr), pfld_ref(ji,jj,jkr-1)
!	    WRITE(numout,*) ' zfld_ref_on_tgt = ', zfld_ref_on_tgt 
!	    WRITE(numout,*) ' pfld_tgt(ji+ki_off_tgt, jj+kj_off_tgt, jk) = ', pfld_tgt(ji+ki_off_tgt, jj+kj_off_tgt, jk)
!         END IF             
      
      END_3D   

      IF ( ln_dbg_hpg ) CALL dbg_3dr( 'ref_to_tgt_ccs_str_off: pfld_tgt_ref', pfld_tgt_ref ) 

   END SUBROUTINE ref_to_tgt_ccs_str_off

!-----------------------------------------------------------------------------------------------------------------

   SUBROUTINE loc_ref_tgt (ll_ccs, ki_off_tgt, kj_off_tgt, p_dep_tgt, p_z_ref, kk_bot_ref, kk_ref_for_tgt, ll_tgt_off_cen ) 

      LOGICAL,                              INTENT(in)   :: ll_ccs           ! true => constrained cubic spline; 
                                                                             ! false => pure cubic fit 
      INTEGER,                              INTENT(in)   :: ki_off_tgt       ! offset of points in target array in i-direction   
      INTEGER,                              INTENT(in)   :: kj_off_tgt       ! offset of points in target array in j-direction   
      REAL(wp), DIMENSION(A2D(nn_hls),jpk), INTENT(in)   :: p_dep_tgt        ! depths of target profiles       
      REAL(wp), DIMENSION(A2D(nn_hls),jpk), INTENT(in)   :: p_z_ref          ! heights of reference  profiles       
      INTEGER,  DIMENSION(A2D(nn_hls)),     INTENT(in)   :: kk_bot_ref       ! bottom levels in the reference profile
      INTEGER,  DIMENSION(A2D(nn_hls),jpk), INTENT(OUT)  :: kk_ref_for_tgt   ! reference index for interpolation to target grid 
                                                                             ! (the lower point)
                                                                             ! with jkr = kk_ref_for_tgt(ji,jj,jk) the reference 
                                                                             ! levels are jkr-1 and jkr
      LOGICAL,  DIMENSION(A2D(nn_hls),jpk), INTENT(OUT)  :: ll_tgt_off_cen   ! T => target point is off-centred 
                                                                             ! (only applies when ll_ccs is False) 
       
      !-----------------------------------------------------------------------------------------------------------------

      INTEGER,  DIMENSION(A2D(nn_hls))      :: jk_tgt, jk_ref                ! vertical level being processed on target 
                                                                             ! and reference grids respectively 
      INTEGER,  DIMENSION(A2D(nn_hls))      :: z_mbkt_off, z_smsk_off
      REAL(wp), DIMENSION(jpk)              :: zdep_delta, zdep_ref          ! for algorithm 2
      INTEGER                               :: ji, jj, jk_comb, jk_bot  
      INTEGER                               :: jk, jkr, jcount               ! for debugging only 
      REAL(wp)                              :: zdep_off_tgt, zmsk_off_tgt    ! for algorithm 2
      REAL(wp)                              :: zdep_ref_up, zdep_ref_dw      ! for debugging only algorithm 2
      REAL(wp)                              :: z_tgt, z_below, z_above       ! for debugging only algorithm 1
      INTEGER                               :: jk_init                       ! initial jk value for search 
      REAL(wp)                              :: tol_wp                        ! this should be the working precision tolerance 

      tol_wp = 1.E-4 ! the inferred bottom depth seems to be accurate to about 1.E-6.  

      !-----------------------------------------------------------------------------------------------------------------

      ! 1. Initialisation for search for points on reference grid bounding points on the target grid
      ! the first point on the target grid is assumed to lie between the first two points on the reference grid 

      jk_init = 2  !   for all cases; enforcement of off-centred interpolation is now done at the end of the routine 

      jk_ref(:,:) = jk_init
      kk_ref_for_tgt(:,:,1) = jk_init
      jk_tgt(:,:) = jk_init

      ll_tgt_off_cen(:,:,:) = .FALSE.     

      ! 2. Find kk_ref_for_tgt

      ! 2.1 Number of wet-levels and land-sea mask at the surface 
      !     of (ji+ki_off_tgt, jj+kj_off_tgt) stencil point:
      DO_2D( 0, 0, 0, 0 )
         z_smsk_off(ji,jj) = ssmask(ji+ki_off_tgt, jj+kj_off_tgt)
         z_mbkt_off(ji,jj) = mbkt(ji+ki_off_tgt, jj+kj_off_tgt)
      END_2D

      SELECT CASE (nn_loc_ref_tgt)
         CASE (1)
            ! Original Mike's algorithm
            DO jk_comb = 1, 2*(jpk-1) 
               ! if level jk_tgt in target profile is lower than jk_ref in reference profile add one to jk_ref
               DO_2D( 0, 0, 0, 0 ) 
                  IF ( - p_dep_tgt( ji+ki_off_tgt, jj+kj_off_tgt, jk_tgt(ji,jj) ) <  p_z_ref( ji, jj, jk_ref(ji,jj) ) ) THEN
                     IF ( jk_ref(ji,jj) < jpk-1 ) jk_ref(ji,jj) = jk_ref(ji,jj) + 1 
                  END IF 
               END_2D
               
               ! if level jk_tgt lies above or at same level as jk_ref, store jk_ref for this level and add one to jk_tgt
               DO_2D( 0, 0, 0, 0 ) 
                  IF ( - p_dep_tgt( ji+ki_off_tgt, jj+kj_off_tgt, jk_tgt(ji,jj) ) + tol_wp > p_z_ref( ji, jj, jk_ref(ji,jj) ) ) THEN
                     IF ( jk_tgt(ji,jj) < jpk ) THEN
                        kk_ref_for_tgt( ji, jj, jk_tgt(ji,jj) ) = jk_ref(ji,jj)         
                        jk_tgt(ji,jj) = jk_tgt(ji,jj) + 1 
                     END IF
                  END IF
               END_2D
  
               IF ( lwp .AND. ln_dbg_hpg ) THEN
                  CALL dbg_2di_k( 'jk_ref', jk_ref, jk_comb )
                  CALL dbg_2di_k( 'jk_tgt', jk_tgt, jk_comb )
               END IF
            END DO ! jk_comb 
         CASE (2)
            ! Alternative Diego's algorithm
            ! jk_ref for each jk_tgt
            DO jk = jk_init, jpkm1
               DO_2D( 0, 0, 0, 0 )
                  jk_tgt(ji,jj) = jk
                  zdep_off_tgt  = p_dep_tgt(ji+ki_off_tgt, jj+kj_off_tgt, jk_tgt(ji,jj))
                  zdep_ref(:)   = - p_z_ref(ji,jj,:)
                  zdep_delta(:) = zdep_ref(:) - (zdep_off_tgt - tol_wp)
                  IF ( ALL(zdep_delta(:) < 0._wp) ) THEN
                     IF ( z_smsk_off(ji,jj) > 0 ) THEN
                        WRITE(numout,*) 'ji, jj, ki_off_tgt, kj_off_tgt'
                        WRITE(numout,*)  ji, jj, ki_off_tgt, kj_off_tgt
                        WRITE(numout,*) 'zdep_ref(:)', zdep_ref(:)
                        WRITE(numout,*) 'zdep_off_tgt', zdep_off_tgt
                        WRITE(numout,*) 'zdep_delta(:)', zdep_delta(:)
                        WRITE(numout,*) 'z_smsk_off(ji,jj)', z_smsk_off(ji,jj)
                        CALL ctl_stop( 'loc_ref_tgt : kk_ref_for_tgt will not cover all sea-points' )
                     END IF
                  ELSE
                     jk_ref(ji,jj) = MAX( jk_init, MINLOC( zdep_delta, DIM=1, MASK=(zdep_delta > 0._wp) ) )
                  END IF
                  !
                  kk_ref_for_tgt( ji, jj, jk ) = jk_ref(ji,jj)
                  !
               END_2D
            END DO
      END SELECT

      ! 3. Checks to make sure we do not read or write out of bounds 

      ! 3.1 First test on jk_tgt to check that it reaches the bottom level mbkt
      ! DB. I don't think we need this with nn_loc_ref_tgt = 2 since in the 
      !     previous computation we are sure we cover all the sea-point 
      jcount = 0 
      DO_2D(0,0,0,0)
         IF ( jk_tgt(ji,jj) <  z_mbkt_off(ji,jj) ) jcount = jcount + 1 
      END_2D

      IF ( jcount > 0 ) THEN 
         WRITE( numout,*) 'loc_ref_tgt: stopping because kk_ref_for_tgt will not cover all sea-points; jcount = ', jcount 
         CALL ctl_stop( 'dyn_hpg_djr : kk_ref_for_tgt will not cover all sea-points' ) 
      END IF 

      ! 3.2 kk_ref_for_tgt pointing to invalid level in reference profile
      ! DB. I also think we don't need this with nn_loc_ref_tgt = 2
      jcount = 0 
      DO_2D(0,0,0,0)
         IF ( z_smsk_off(ji,jj) > 0 ) THEN
            IF ( kk_ref_for_tgt(ji,jj,z_mbkt_off(ji,jj)) > kk_bot_ref(ji,jj) ) jcount = jcount + 1
         END IF
      END_2D

      IF ( jcount > 0 ) THEN 
         WRITE( numout,*) 'loc_ref_tgt: stopping because kk_ref_for_tgt points to a level below the bottom of the reference profile; jcount = ', jcount 
         CALL ctl_stop( 'dyn_hpg_djr : kk_ref_for_tgt points to a level below the bottom of the reference profile' ) 
      END IF 
   
      SELECT CASE (nn_loc_ref_tgt)
         CASE (1)
            ! Original Mike's algorithm
            ! 3.3 diagnostic check that the right reference levels are chosen for the target profile   
            IF ( lwp .AND. ln_dbg_hpg ) THEN
               WRITE( numout,*) 'loc_ref_tgt: ki_off_tgt, kj_off_tgt = ', ki_off_tgt, kj_off_tgt
               CALL dbg_3dr( '-p_dep_tgt', -p_dep_tgt ) 
               CALL dbg_3dr( 'p_z_ref', p_z_ref ) 
               CALL dbg_3di( 'kk_ref_for_tgt', kk_ref_for_tgt ) 
  
               jcount = 0 
               DO_3D(0,0,0,0,2,jpk-1) 
                  z_tgt = - p_dep_tgt( ji+ki_off_tgt, jj+kj_off_tgt, jk)
                  jkr = kk_ref_for_tgt(ji,jj,jk)
                  z_below = p_z_ref( ji, jj, jkr )
                  z_above = p_z_ref( ji, jj, jkr-1 )
                  IF (tmask(ji+ki_off_tgt, jj+kj_off_tgt, jk) > 0) THEN
                     IF ( ( z_above < z_tgt .AND. jkr > jk_init )  .OR. z_below > (z_tgt + tol_wp) ) THEN
                        IF ( jcount < 10 ) THEN
                           WRITE(numout,*) 'loc_ref_tgt: ji, jj, jk, ki_off_tgt, kj_off_tgt, jkr, z_tgt, z_below, z_above ' 
                           WRITE(numout,*) ji, jj, jk, ki_off_tgt, kj_off_tgt, jkr, z_tgt, z_below, z_above 
                        END IF
                        jcount = jcount + 1 
                     END IF
                  END IF
               END_3D
               WRITE(numout,*) 'loc_ref_tgt: jcount = ', jcount 
            END IF 
  
            ! 3.4 Same check as in 4.2 but detailed diagnostics not written out.
            jcount = 0 
            DO_3D(0,0,0,0,2,jpk-1) 
               z_tgt = - p_dep_tgt( ji+ki_off_tgt, jj+kj_off_tgt, jk)
               jkr = kk_ref_for_tgt(ji,jj,jk)
               z_below = p_z_ref( ji, jj, jkr )
               z_above = p_z_ref( ji, jj, jkr-1 )
               IF (tmask(ji+ki_off_tgt, jj+kj_off_tgt, jk) > 0) THEN
                  IF ( ( z_above < z_tgt .AND. jkr > jk_init)  .OR. z_below > (z_tgt + tol_wp) ) THEN
                     IF ( jcount < 10 ) THEN
                        WRITE(numout,*) 'loc_ref_tgt: ji, jj, jk, ki_off_tgt, kj_off_tgt, jkr, z_tgt, z_below, z_above '
                        WRITE(numout,*) ji, jj, jk, ki_off_tgt, kj_off_tgt, jkr, z_tgt, z_below, z_above
                     END IF
                     jcount = jcount + 1 
                  END IF
               END IF
            END_3D
         CASE (2)
            ! Alternative Diego's algorithm
            ! 3.3 diagnostic check that the right reference levels are chosen for the target profile    
            IF ( lwp .AND. ln_dbg_hpg ) THEN 
               WRITE( numout,*) 'loc_ref_tgt: ki_off_tgt, kj_off_tgt = ', ki_off_tgt, kj_off_tgt
               CALL dbg_3dr( '-p_dep_tgt', -p_dep_tgt ) 
               CALL dbg_3dr( 'p_z_ref', p_z_ref ) 
               CALL dbg_3di( 'kk_ref_for_tgt', kk_ref_for_tgt ) 

               jcount = 0
               DO_3D(0,0,0,0,2,jpk-1)
                  zdep_off_tgt = p_dep_tgt(ji+ki_off_tgt, jj+kj_off_tgt, jk)
                  zmsk_off_tgt = tmask(ji+ki_off_tgt, jj+kj_off_tgt, jk)
                  jkr = kk_ref_for_tgt(ji,jj,jk)
                  zdep_ref_up  = - p_z_ref(ji,jj,jkr-1)
                  zdep_ref_dw  = - p_z_ref(ji,jj,jkr)
                  IF ( zmsk_off_tgt > 0 ) THEN
                      IF ( (zdep_off_tgt-tol_wp) < zdep_ref_up .OR. ((zdep_off_tgt-tol_wp) > zdep_ref_dw .AND. jkr > jk_init) ) THEN
                         IF ( jcount < 10 ) THEN
                            WRITE(numout,*) 'loc_ref_tgt: ji, jj, jk, ki_off_tgt, kj_off_tgt, jkr, zdep_off_tgt, zdep_ref_up, zdep_ref_dw '
                            WRITE(numout,*) ji, jj, jk, ki_off_tgt, kj_off_tgt, jkr, zdep_off_tgt, zdep_ref_up, zdep_ref_dw
                         END IF
                         jcount = jcount + 1
                      END IF
                  END IF
               END_3D
               WRITE(numout,*) 'loc_ref_tgt: jcount = ', jcount 
            END IF 

            ! 3.4 Same check as in 3.3 but detailed diagnostics not written out. 
            jcount = 0 
            DO_3D(0,0,0,0,2,jpk-1)
               zdep_off_tgt = p_dep_tgt(ji+ki_off_tgt, jj+kj_off_tgt, jk)
               zmsk_off_tgt = tmask(ji+ki_off_tgt, jj+kj_off_tgt, jk)
               jkr = kk_ref_for_tgt(ji,jj,jk)
               zdep_ref_up  = - p_z_ref(ji,jj,jkr-1)
               zdep_ref_dw  = - p_z_ref(ji,jj,jkr)
               IF ( zmsk_off_tgt > 0 ) THEN
                  IF ( (zdep_off_tgt-tol_wp) < zdep_ref_up .OR. ((zdep_off_tgt-tol_wp) > zdep_ref_dw .AND. jkr > jk_init) ) THEN
                     IF ( jcount < 10 ) THEN
                        WRITE(numout,*) 'loc_ref_tgt: ji, jj, jk, ki_off_tgt, kj_off_tgt, jkr, zdep_off_tgt, zdep_ref_up, zdep_ref_dw '
                        WRITE(numout,*) ji, jj, jk, ki_off_tgt, kj_off_tgt, jkr, zdep_off_tgt, zdep_ref_up, zdep_ref_dw
                     END IF 
                     jcount = jcount + 1  
                  END IF
               END IF
            END_3D
      END SELECT

      IF ( jcount > 0 ) THEN 
         WRITE( numout,*) 'loc_ref_tgt: stopping because jcount is non-zero; jcount = ', jcount 
         CALL ctl_stop( 'dyn_hpg_djr : loc_ref_tgt failed to locate target points correctly' )
      END IF     
      
      ! 4. Adjust kk_ref_for_tgt so that interpolation of simple cubic is off-centred at the bottom 
      !    (does not require boundary conditions) 
      !    This assumes that every sea point has at least 4 levels.   

      IF ( ll_ccs ) THEN
         DO_3D(0, 0, 0, 0, 2, jpk-1) 
            IF ( kk_ref_for_tgt(ji,jj,jk) == kk_bot_ref(ji,jj) ) THEN 
               ll_tgt_off_cen(ji,jj,jk) = .TRUE.            
            END IF 
         END_3D
      ELSE 
         DO_3D(0, 0, 0, 0, 1, jpk-1) 
            IF ( kk_ref_for_tgt(ji,jj,jk) == 2 ) THEN 
               ll_tgt_off_cen(ji,jj,jk) = .TRUE.            
               kk_ref_for_tgt(ji,jj,jk) = 3
            END IF 
            IF ( kk_ref_for_tgt(ji,jj,jk) == kk_bot_ref(ji,jj) ) THEN 
               ll_tgt_off_cen(ji,jj,jk) = .TRUE.            
               kk_ref_for_tgt(ji,jj,jk) =  kk_bot_ref(ji,jj) - 1
            END IF 
         END_3D
      END IF 

      IF ( lwp .AND. ln_dbg_hpg ) THEN 
         WRITE( numout,*) 'loc_ref_tgt at end of section 4 ', ki_off_tgt, kj_off_tgt
         CALL dbg_3di( 'kk_ref_for_tgt', kk_ref_for_tgt ) 
      END IF 

   END SUBROUTINE loc_ref_tgt

!----------------------------------------------------------------------------

   SUBROUTINE hpg_ffr( kt, Kmm, puu, pvv, Krhs  )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE hpg_ffr  ***
      !!
      !! ** Method  :   s-coordinate case forces on faces scheme. 
      !!                Local density subtracted using cubic or constrained cubic splines (ccs)  
      !!                Remaining densities interpolated to centre either linearly or with ccs  
      !!                Vertical representation either using quadratic density or classic 2nd order accurate 
      !!
      !! ** Action : - Update puu(..,Krhs) and pvv(..,Krhs)  with the now hydrastatic pressure trend
      !!----------------------------------------------------------------------
      INTEGER                             , INTENT( in )  ::  kt          ! ocean time-step index
      INTEGER                             , INTENT( in )  ::  Kmm, Krhs   ! ocean time level indices
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpt), INTENT(inout) ::  puu, pvv    ! ocean velocities and RHS of 
                                                                          ! momentum equation
      !!
      
      REAL(wp), DIMENSION(A2D(nn_hls),jpk)   :: zrhd_ref, zdrhd_k_ref  ! densities (rhd) of reference profile
      REAL(wp), DIMENSION(A2D(nn_hls),jpk)   :: zz_ref                 ! heights of reference profile 
      REAL(wp), DIMENSION(A2D(nn_hls),jpk,4) :: z_rhd_pmr              ! profiles minus reference   
      REAL(wp), DIMENSION(A2D(nn_hls),jpk,4) :: zrhd_ref_on_tgt        ! densities (rhd) of ref. prof. DB debug

      ! The following fields could probably be made single level or at most 2 level fields 
      REAL(wp), DIMENSION(A2D(nn_hls),jpk,4) :: z_dept_pmr, z_depw_pmr ! corresponding gdept and gdepw pmr profiles 
      REAL(wp), DIMENSION(A2D(nn_hls),jpk)   :: z_rhd_mid              ! rhd profiles interpolated to middle of cell 
                                                                       ! (using ccs or cub)  
      REAL(wp), DIMENSION(A2D(nn_hls),jpk)   :: z_dept_mid, z_depw_mid ! dept and depw similarly interpolated to 
                                                                       ! middle of cell (using ccs or cub)  
      REAL(wp), DIMENSION(A2D(nn_hls),jpk)   :: z_e3t_mid              ! e3t profiles interpolated to middle of cell; 
                                                                       ! only printed as a diagnostic  
      REAL(wp), DIMENSION(A2D(nn_hls),jpk,3) :: z_ddepw_ij             ! constrained spline horizontal differences 
                                                                       ! in gdepw 
      REAL(wp), DIMENSION(A2D(nn_hls),jpk,3) :: z_p_pmr, z_F_pmr       ! pressures and forces on vertical faces  
      REAL(wp), DIMENSION(A2D(nn_hls),jpk)   :: z_F_upp                ! forces on upper faces

      INTEGER,  DIMENSION(A2D(nn_hls))       :: jk_bot_ref             ! bottom level in reference profiles 
      INTEGER,  DIMENSION(A2D(nn_hls),3)     :: j_mbkt                 ! bottom levels for the 3 profiles in the cell  
 
      REAL(wp), DIMENSION(A2D(nn_hls),jpk)   :: zpforces               ! total of forces 

      INTEGER  ::   ji, jj, jk, j_uv, jr        ! dummy loop indices
      INTEGER  ::   jio, jjo                    ! offsets for the 2 point stencil (0 or 1)  
      INTEGER  ::   ji_ro, jj_ro, jr_offset     ! offsets for the reference profile 
      INTEGER  ::   jn_hor_pts                  ! number of profiles required in horizontal 
                                                ! (4 if cubic interpolarion or 2 if not) 
      INTEGER  ::   j_t_levs, j_w_levs          ! indicators that field passed in is valid 
                                                ! on t-levels or w-levels
      
      REAL(wp) ::   z_r_6                       ! 1/6 
      REAL(wp) ::   zr_a, zr_b, zr_c            ! rhd evaluated at i, i+1/2 and i+1  
      REAL(wp) ::   za_0, za_1, za_2            ! polynomial coefficients
      LOGICAL  ::   ln_dbg                      ! DB for debugging
      !!----------------------------------------------------------------------

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn:hpg_Lin : hydrostatic pressure gradient trend'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~   s-coordinate case, forces on faces with reference removed'
      ENDIF
      !
      
      z_r_6 = 1.0_wp/6.0_wp
      
      j_t_levs = 1      ! this is not very neat - it is duplicated in the calling routines
      j_w_levs = 2 

      IF ( ln_hpg_ffr_hor_ccs .OR. ln_hpg_ffr_hor_cub ) THEN  
         jn_hor_pts = 4
         jr_offset  = 2
      ELSE 
         jn_hor_pts = 2
         jr_offset  = 1 
      END IF 

      !WRITE( numout, *) ' hpg_ffr: kt, jn_hor_pts, jr_offset = ', kt, jn_hor_pts, jr_offset  

      DO j_uv = 1, 2 
         
         ! set for the whole of the j_uv loop 
         !(loop ends near the bottom of the routine)
         IF ( j_uv == 1 ) THEN 
            jio = 1 
            jjo = 0
            ln_dbg = .TRUE. 
         ELSE 
            jio = 0
            jjo = 1
            ln_dbg = .FALSE. 
         END IF 

         !-------------------------------------------------------------------
         ! 1. Calculate the reference profile from all points in the stencil

         IF ( ln_hpg_ref ) THEN 
            CALL calc_rhd_ref(j_uv, jn_hor_pts, gdept(:,:,:,Kmm), zrhd_ref, zz_ref, jk_bot_ref)    
            IF ( ln_hpg_ref_ccs ) CALL calc_drhd_k(zrhd_ref, jk_bot_ref, zdrhd_k_ref) 
         END IF 

         IF ( ln_dbg_hpg ) THEN 
            CALL dbg_3dr( '2. rhd', rhd ) 
            CALL dbg_3dr( '2. gdept', gdept(:,:,:,Kmm) ) 
         END IF 

         !-------------------------------------------------------------------
         ! 2. Interpolate reference profile to all points in the stencil and 
         !    calculate z_rhd_pmr (profile minus reference)

         DO jr = 1, jn_hor_pts       

           ! range of jio is 
           !  from -1 to 2 if jn_hor_pts = 4;
           !  from  0 to 1 if jn_hor_pts = 2;
           IF ( j_uv == 1 ) THEN 
             ji_ro = jr - jr_offset 
             jj_ro = 0 
           ELSE 
             ji_ro = 0
             jj_ro = jr - jr_offset 
           END IF 

           ! Subtraction of reference profile if requested
           IF ( ln_hpg_ref ) THEN

              IF ( ln_hpg_ref_str ) THEN 
                 IF ( ln_hpg_ref_ccs ) THEN 
                    IF ( ln_hpg_ref_off ) THEN 
                       !CALL ref_to_tgt_ccs_str_off( ji_ro, jj_ro, gdept, rhd, zz_ref, zrhd_ref, &
                       ! &                           zdrhd_k_ref, jk_bot_ref, z_rhd_pmr(:,:,:,jr) )
                       CALL ref_to_tgt_ccs_str_off( ji_ro, jj_ro, gdept(:,:,:,Kmm), rhd, zz_ref, zrhd_ref, zdrhd_k_ref, &
                        &                           jk_bot_ref, z_rhd_pmr(:,:,:,jr) )
                    ELSE
                       !CALL ref_to_tgt_ccs_str( ji_ro, jj_ro, gdept, rhd, zz_ref, zrhd_ref, &
                       ! &                       zdrhd_k_ref, jk_bot_ref, z_rhd_pmr(:,:,:,jr) )
                       CALL ref_to_tgt_ccs_str( ji_ro, jj_ro, gdept(:,:,:,Kmm), rhd, zz_ref, zrhd_ref, zdrhd_k_ref, &
                        &                       jk_bot_ref, z_rhd_pmr(:,:,:,jr) )
                    END IF 
                 ELSE  
                    !CALL ref_to_tgt_cub_str ( ji_ro, jj_ro, gdept, rhd, zz_ref, zrhd_ref,  &
                    ! &                        jk_bot_ref, z_rhd_pmr(:,:,:,jr) )  
                    CALL ref_to_tgt_cub_str( ji_ro, jj_ro, gdept(:,:,:,Kmm), rhd, zz_ref, zrhd_ref, &
                     &                       jk_bot_ref, z_rhd_pmr(:,:,:,jr) )
                    !CALL ref_to_tgt_cub_str_dbg( ln_dbg, ji_ro, jj_ro, gdept(:,:,:,Kmm), rhd, zz_ref, zrhd_ref, &
                    ! &                           jk_bot_ref, z_rhd_pmr(:,:,:,jr) )
                 END IF
              ELSE       ! these calls are mainly retained to assist in testing 
                 IF ( ln_hpg_ref_ccs ) THEN 
                    CALL ref_to_tgt_ccs( ji_ro, jj_ro, gdept, rhd, zz_ref, zrhd_ref, &
                     &                   zdrhd_k_ref, jk_bot_ref, z_rhd_pmr(:,:,:,jr) )
                 ELSE  
                    CALL ref_to_tgt_cub ( ji_ro, jj_ro, gdept, rhd, zz_ref, zrhd_ref, &
                     &                    jk_bot_ref, z_rhd_pmr(:,:,:,jr) )  
                 END IF
              END IF

           ELSE !  NO subtraction of reference profile 

             DO_3D (0,0,0,0,1,jpk)
               z_rhd_pmr(ji,jj,jk,jr) = rhd(ji+ji_ro,jj+jj_ro,jk)
             END_3D

           END IF ! ln_hpg_ref 

           DO_3D (0,0,0,0,1,jpk)
               z_dept_pmr(ji,jj,jk,jr) = gdept(ji+ji_ro,jj+jj_ro,jk,Kmm)
               z_depw_pmr(ji,jj,jk,jr) = gdepw(ji+ji_ro,jj+jj_ro,jk,Kmm)
           END_3D

           IF ( ln_dbg_hpg ) THEN 
              CALL dbg_3dr( '2. z_rhd_pmr', z_rhd_pmr(:,:,:,jr) ) 
           END IF 

         END DO !jr = 1, jn_hor_pts    

         !-------------------------------------------------------------------
         ! 3. Horizontal interpolation to compute intermediate densities 
         !    (either linear or cubic or constrained cubic). Transfers 
         !    data points from reference stencil (2 or 4 point) to a 3 points 
         !    horizontal stencil

         IF ( ln_hpg_ffr_hor_ccs .OR. ln_hpg_ffr_hor_cub) THEN  

            IF ( ln_hpg_ffr_hor_ccs ) THEN
               CALL calc_mid_ccs(j_uv, j_t_levs, aco_bc_rhd_hor, bco_bc_rhd_hor, z_rhd_pmr,  z_rhd_mid)  
               CALL calc_mid_ccs(j_uv, j_t_levs, aco_bc_z_hor,   bco_bc_z_hor,   z_dept_pmr, z_dept_mid)  
               CALL calc_mid_ccs(j_uv, j_w_levs, aco_bc_z_hor,   bco_bc_z_hor,   z_depw_pmr, z_depw_mid)  
               CALL calc_dz_dij_ccs( j_uv, z_depw_pmr, z_ddepw_ij )
            ELSE ! ln_hpg_ffr_hor_cub
               CALL calc_mid_cub(j_uv, j_t_levs, aco_bc_rhd_hor, bco_bc_rhd_hor, z_rhd_pmr,  z_rhd_mid)  
               CALL calc_mid_cub(j_uv, j_t_levs, aco_bc_z_hor,   bco_bc_z_hor,   z_dept_pmr, z_dept_mid)  
               CALL calc_mid_cub(j_uv, j_w_levs, aco_bc_z_hor,   bco_bc_z_hor,   z_depw_pmr, z_depw_mid)  
               CALL calc_dz_dij_cub( j_uv, z_depw_pmr, z_ddepw_ij )
            END IF 

            DO_3D (0,0,0,0,1,jpk)  
               z_rhd_pmr(ji,jj,jk,1) = z_rhd_pmr(ji,jj,jk,2) 
               z_rhd_pmr(ji,jj,jk,2) = z_rhd_mid(ji,jj,jk)
               z_dept_pmr(ji,jj,jk,1) = z_dept_pmr(ji,jj,jk,2) 
               z_dept_pmr(ji,jj,jk,2) = z_dept_mid(ji,jj,jk)
               z_depw_pmr(ji,jj,jk,1) = z_depw_pmr(ji,jj,jk,2) 
               z_depw_pmr(ji,jj,jk,2) = z_depw_mid(ji,jj,jk)
            END_3D

         ELSE !  simple linear interpolation 

            DO_3D (0,0,0,0,1,jpk)  
               z_rhd_pmr(ji,jj,jk,3) =  z_rhd_pmr(ji,jj,jk,2) 
               z_rhd_pmr(ji,jj,jk,2) =  0.5_wp*( z_rhd_pmr(ji,jj,jk,1)  + z_rhd_pmr(ji,jj,jk,2) )
               z_dept_pmr(ji,jj,jk,3) = z_dept_pmr(ji,jj,jk,2)    
               z_dept_pmr(ji,jj,jk,2) = 0.5_wp*( z_dept_pmr(ji,jj,jk,1) + z_dept_pmr(ji,jj,jk,2)  )   
               z_depw_pmr(ji,jj,jk,3) = z_depw_pmr(ji,jj,jk,2)    
               z_depw_pmr(ji,jj,jk,2) = 0.5_wp*( z_depw_pmr(ji,jj,jk,1) + z_depw_pmr(ji,jj,jk,2)  )   
            END_3D
         END IF 

         DO_2D (0,0,0,0) 
            j_mbkt(ji,jj,1) = mbkt(ji, jj) 
            j_mbkt(ji,jj,2) = MIN( mbkt(ji, jj), mbkt(ji+jio, jj+jjo) ) 
            j_mbkt(ji,jj,3) = mbkt(ji+jio, jj+jjo) 
         END_2D

         IF ( ln_dbg_hpg ) THEN 
            CALL dbg_3dr( '3. z_rhd_pmr: 1', z_rhd_pmr(:,:,:,1) ) 
            CALL dbg_3dr( '3. z_rhd_pmr: 2', z_rhd_pmr(:,:,:,2) ) 
            CALL dbg_3dr( '3. z_rhd_pmr: 3', z_rhd_pmr(:,:,:,3) ) 
            CALL dbg_3dr( '3. z_depw_pmr: 1', z_depw_pmr(:,:,:,1) ) 
            CALL dbg_3dr( '3. z_depw_pmr: 2', z_depw_pmr(:,:,:,2) ) 
            CALL dbg_3dr( '3. z_depw_pmr: 3', z_depw_pmr(:,:,:,3) ) 
            DO jr = 1, 3
               DO_3D (0,0,0,0,1,jpk-1)
                  z_e3t_mid(ji,jj,jk) = z_depw_pmr(ji,jj,jk,jr) - z_depw_pmr(ji,jj,jk+1,jr)
               END_3D
               CALL dbg_3dr( '3. z_e3t_mid jr : ', z_e3t_mid(:,:,:) ) 
            END DO 
         END IF 

         !-------------------------------------------------------------------
         ! 4. Vertical interpolation of densities, calculating pressures and 
         !    forces on vertical faces between w levels 	    

         IF ( ln_hpg_ffr_vrt_quad ) THEN 
            DO jr = 1, 3 
               CALL vrt_int_quad( j_mbkt(:,:,jr)      , z_rhd_pmr(:,:,:,jr), z_dept_pmr(:,:,:,jr), &
                &                 z_depw_pmr(:,:,:,jr), z_p_pmr(:,:,:,jr)  , z_F_pmr(:,:,:,jr) ) 
            END DO ! jr

         ELSE  !  .NOT. ln_hpg_ffr_vrt_quad   (simplest vertical integration scheme)  

            DO jr = 1, 3 
               DO_2D(0,0,0,0) 
                  z_p_pmr(ji,jj,1,jr) = 0._wp 
               END_2D
               DO jk = 1, jpk - 1
                  DO_2D(0,0,0,0) 
                     z_p_pmr(ji,jj,jk+1,jr) = z_p_pmr(ji,jj,jk,jr) + grav*(z_depw_pmr(ji,jj,jk+1,jr) &
                       &                      - z_depw_pmr(ji,jj,jk,jr)) *z_rhd_pmr(ji,jj,jk,jr)
                  END_2D
                  DO_2D(0,0,0,0) 
                     z_F_pmr(ji,jj,jk,jr) = 0.5_wp*(z_depw_pmr(ji,jj,jk+1,jr) - z_depw_pmr(ji,jj,jk,jr)) &
                       &                    *( z_p_pmr(ji,jj,jk,jr) + z_p_pmr(ji,jj,jk+1,jr) ) 
                  END_2D
               END DO ! jk 
            END DO ! jr

         END IF  !  ln_hpg_ffr_vrt_quad   

         IF ( ln_dbg_hpg ) THEN 
            CALL dbg_3dr( '4. z_p_pmr: 1', z_p_pmr(:,:,:,1) ) 
            CALL dbg_3dr( '4. z_p_pmr: 2', z_p_pmr(:,:,:,2) ) 
            CALL dbg_3dr( '4. z_p_pmr: 3', z_p_pmr(:,:,:,3) ) 
            CALL dbg_3dr( '4. z_F_pmr: 1', z_F_pmr(:,:,:,1) ) 
            CALL dbg_3dr( '4. z_F_pmr: 3', z_F_pmr(:,:,:,3) ) 
         END IF 

         !-------------------------------------------------------------------
         ! 5. Calculate forces on the upper faces and hence on the total 
         !    forces on the cells (zpforces)  

         DO_2D(0,0,0,0) 
            z_F_upp(ji,jj,1) = 0._wp
         END_2D 

         IF ( ln_hpg_ffr_hor_ccs .OR. ln_hpg_ffr_hor_cub ) THEN     ! use Simpson's rule 
            DO_3D(0,0,0,0,2,jpk) 
               z_F_upp(ji,jj,jk) = - z_r_6 * (     z_ddepw_ij(ji,jj,jk,1)*z_p_pmr(ji,jj,jk,1)  & 
               &                           + 4._wp*z_ddepw_ij(ji,jj,jk,2)*z_p_pmr(ji,jj,jk,2)  & 
               &                           +       z_ddepw_ij(ji,jj,jk,3)*z_p_pmr(ji,jj,jk,3) )   
            END_3D
         ELSE                               ! use trapezoidal rule 
            DO_3D(0,0,0,0,2,jpk) 
               z_F_upp(ji,jj,jk) = 0.5_wp * ( gdepw(ji,jj,jk,Kmm) - gdepw(ji+jio,jj+jjo,jk,Kmm) ) &
               &                          * ( z_p_pmr(ji,jj,jk,1) + z_p_pmr(ji,jj,jk,3) ) 
            END_3D
         END IF 

         IF ( ln_dbg_hpg ) THEN 
            CALL dbg_3dr( '4. z_F_upp: ', z_F_upp ) 
            CALL dbg_3dr( '4. gdepw: ', gdepw )
            CALL dbg_3dr( '4. z_ddepw_ij 1: ', z_ddepw_ij(:,:,:,1) ) 
            CALL dbg_3dr( '4. z_ddepw_ij 2: ', z_ddepw_ij(:,:,:,2) ) 
            CALL dbg_3dr( '4. z_ddepw_ij 3: ', z_ddepw_ij(:,:,:,3) ) 
         END IF 

         IF ( j_uv == 1 ) THEN 
            DO_3D(0,0,0,0,1,jpk-1)
               zpforces(ji,jj,jk) =   z_F_pmr(ji,jj,jk,1) - z_F_pmr(ji,jj,jk,3) &
                  &                 + z_F_upp(ji,jj,jk)   - z_F_upp(ji,jj,jk+1) 
               puu(ji,jj,jk,Krhs) =   puu(ji,jj,jk,Krhs) + (zpforces(ji,jj,jk) / (e1u(ji,jj)*e3u(ji,jj,jk,Kmm))) &
                  &                 * umask(ji,jj,jk) 
            END_3D
            IF ( ln_dbg_hpg ) CALL dbg_3dr( '5. u zpforces: ', zpforces ) 

         ELSE 
            DO_3D(0,0,0,0,1,jpk-1)
               zpforces(ji,jj,jk) =   z_F_pmr(ji,jj,jk,1) - z_F_pmr(ji,jj,jk,3) & 
                  &                 + z_F_upp(ji,jj,jk) - z_F_upp(ji,jj,jk+1) 
               pvv(ji,jj,jk,Krhs) =   pvv(ji,jj,jk,Krhs) + (zpforces(ji,jj,jk) / (e1v(ji,jj)*e3v(ji,jj,jk,Kmm))) &
                  &                 * vmask(ji,jj,jk)   
            END_3D
            IF ( ln_dbg_hpg ) CALL dbg_3dr( '5. v zpforces: ', zpforces ) 
         END IF 

         ! temporary output of fields for debugging etc. 
         IF ( j_uv == 1) THEN 
            !CALL iom_put( "gdepw_hpg", gdepw(:,:,:,Kmm) ) ! DB: We already have this field since 
                                                           !     we have e3t in grid_T
            CALL iom_put( "e3w", e3w(:,:,:,Kmm) )
            CALL iom_put( "rhd_hpg", z_rhd_pmr(:,:,:,1) )
            IF ( ln_hpg_ref ) THEN
               CALL iom_put( "rhd_ref", zrhd_ref )
               !CALL iom_put( "zz_ref", zz_ref)
               CALL iom_put( "jk_bot_ref", FLOAT( jk_bot_ref(:,:) ) )
            ENDIF
            CALL iom_put( "pressure",  z_p_pmr(:,:,:,1) )
            CALL iom_put( "u_force_west", z_F_pmr(:,:,:,1) )
            CALL iom_put( "u_force_upper", z_F_upp )
            !CALL iom_put( "mid_x_press", z_p_pmr(:,:,:,2) ) 
            !CALL iom_put( "x_mid_slope", z_ddepw_ij(:,:,:,2) ) 
            !CALL iom_put( "lhs_x_press", z_p_pmr(:,:,:,1) ) 
            !CALL iom_put( "x_lhs_slope", z_ddepw_ij(:,:,:,1) ) 
            !CALL iom_put( "rhs_x_press", z_p_pmr(:,:,:,3) ) 
            !CALL iom_put( "x_rhs_slope", z_ddepw_ij(:,:,:,3) ) 
         ELSE 
            CALL iom_put( "v_force_south", z_F_pmr(:,:,:,1) )
            CALL iom_put( "v_force_upper", z_F_upp )
            !CALL iom_put( "mid_y_press", z_p_pmr(:,:,:,2) ) 
            !CALL iom_put( "y_mid_slope", z_ddepw_ij(:,:,:,2) ) 
         END IF 

      END DO  ! j_uv  
      !
   END SUBROUTINE hpg_ffr

!------------------------------------------------------------------------------------------
   
   SUBROUTINE calc_rhd_ref (k_uv, kn_hor, pdep, prhd_ref, pz_ref, kk_bot_ref) 

      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE calc_rhd_ref  ***
      !!
      !! ** Method: Find the deepest cell within the stencil. 
      !!              (Later will extend to producing a reference profile 
      !!               that spans the highest and lowest points in the stencil)   
      !!                 
      !!                 
      !!
      !! ** Action: -  Set prhd_ref, pz_ref, kk_bot_ref
      !!----------------------------------------------------------------------
      INTEGER                   , INTENT(in)  ::  k_uv       ! 1 for u-vel; 2 for v_vel
      INTEGER                   , INTENT(in)  ::  kn_hor     ! 4 for cubic; 2 for linear
      REAL(wp), DIMENSION(:,:,:), INTENT(in)  ::  pdep       ! depths (either gdept or gde3w)
      REAL(wp), DIMENSION(:,:,:), INTENT(out) ::  prhd_ref   ! densities    of reference profile
      REAL(wp), DIMENSION(:,:,:), INTENT(out) ::  pz_ref     ! heights      of   "         " 
      INTEGER,  DIMENSION(:,:),   INTENT(out) ::  kk_bot_ref ! bottom level of   "         " 

      INTEGER,  DIMENSION(A2D(nn_hls))        :: ji_ref, jj_ref      
      REAL(wp), DIMENSION(A2D(nn_hls))        :: max_wetdep

      INTEGER  ji, jj, jk    ! standard loop indices
      INTEGER  jio, jjo      ! offsets 
      INTEGER  jib, jjb      ! second set of indices to deeper points 
      INTEGER  jir, jjr      ! indices of reference profile  
      REAL(wp) zhta, zhtb

   IF (k_uv == 1) THEN
      jio = 1
      jjo = 0 
   ELSE 
      jio = 0 
      jjo = 1
   END IF 

   ! DB using depths of the last wet T-level
   ! instead of ht_0 (which is time invariant) 
   ! because model levels are quasi-eulerian, i.e. 
   ! they oscillate with barotropic motion
   DO_2D( 0, 0, 0, 0 )
      max_wetdep(ji,jj) = pdep(ji,jj,mbkt(ji,jj))
   END_2D
   CALL lbc_lnk ( 'calc_rhd_ref', max_wetdep, 'T', 1._wp )

   DO_2D( 0, 0, 0, 0 )   

      !IF ( ht_0(ji+jio,jj+jjo) >= ht_0(ji,jj) ) THEN  
      IF ( max_wetdep(ji+jio,jj+jjo) > max_wetdep(ji,jj) ) THEN ! DB: for points at the boundaries
         zhta = max_wetdep(ji+jio,jj+jjo) 
         ji_ref(ji,jj) = ji+jio
         jj_ref(ji,jj) = jj+jjo
      ELSE 
         zhta = max_wetdep(ji,jj) 
         ji_ref(ji,jj) = ji
         jj_ref(ji,jj) = jj
      END IF 
      
      IF ( kn_hor == 4 ) THEN
         IF ( max_wetdep(ji-jio,jj-jjo) > max_wetdep(ji+2*jio,jj+2*jjo) ) THEN  ! ??? DB >= or >
            zhtb = max_wetdep(ji-jio,jj-jjo) 
            jib = ji-jio
            jjb = jj-jjo
         ELSE 
            zhtb = max_wetdep(ji+2*jio,jj+2*jjo) 
            jib = ji+2*jio
            jjb = jj+2*jjo
         END IF
         IF ( zhta < zhtb ) THEN  
            ji_ref(ji,jj) = jib 
            jj_ref(ji,jj) = jjb 
         END IF
      END IF 

   END_2D 

   DO_3D( 0, 0, 0, 0, 1, jpk-1)
      jir = ji_ref(ji,jj) 
      jjr = jj_ref(ji,jj) 
      prhd_ref (ji,jj,jk) =    rhd(jir, jjr, jk) 
      pz_ref   (ji,jj,jk) = - pdep(jir, jjr, jk) 
   END_3D 

   DO_2D( 0, 0, 0, 0) 
      jir = ji_ref(ji,jj) 
      jjr = jj_ref(ji,jj) 
      kk_bot_ref(ji,jj) = mbkt(jir,jjr)
   END_2D

   END SUBROUTINE calc_rhd_ref

!------------------------------------------------------------------------------------------
   
   SUBROUTINE calc_drhd_k(p_rhd, kk_bot, p_drhd_k) 

      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE calc_drhd_k  ***
      !!
      !! ** Method  :  Calculate harmonic averages of vertical differences and apply upper and lower boundary conditions 
      !!                  
      !! ** Action : - Set p_drhd_k
      !!----------------------------------------------------------------------

      REAL(wp), DIMENSION(A2D(nn_hls),jpk), INTENT(in)  ::  p_rhd    ! densities    of profile
      INTEGER,  DIMENSION(A2D(nn_hls)),     INTENT(in)  ::  kk_bot   ! bottom level of profile
      REAL(wp), DIMENSION(A2D(nn_hls),jpk), INTENT(out) ::  p_drhd_k ! harmonic mean of vertical differences of profile

      INTEGER  ji, jj, jk    ! standard loop indices
      INTEGER  jio, jjo      ! offsets 
      INTEGER  iktb          ! index of the bottom of ref profile
      REAL(wp), DIMENSION(A2D(nn_hls),jpk)  ::  ztmp     ! primitive vertical differences 

      REAL(wp) cffw, z1_cff, zep

      DO_3D( 0, 0, 0, 0, 2, jpk )  
          ztmp(ji,jj,jk) = p_rhd(ji,jj,jk) - p_rhd(ji,jj,jk-1)
      END_3D

      zep = 1.e-15
      DO_3D( 0, 0, 0, 0, 2, jpkm1 ) 
         cffw = MAX( 2._wp * ztmp(ji,jj,jk) * ztmp(ji,jj,jk+1), 0._wp )
         z1_cff = ztmp(ji,jj,jk) + ztmp(ji,jj,jk+1) 
         p_drhd_k(ji,jj,jk) = cffw / SIGN( MAX( ABS(z1_cff), zep ), z1_cff )
      END_3D

      DO_2D( 0, 0, 0, 0 )
         p_drhd_k(ji,jj,1) = aco_bc_rhd_srf * ( p_rhd(ji,jj,2) - p_rhd(ji,jj,1) ) - bco_bc_rhd_srf * p_drhd_k(ji,jj,2)
         iktb = kk_bot(ji,jj)  
         IF ( iktb > 1 ) THEN
            p_drhd_k(ji,jj,iktb) = aco_bc_rhd_bot * (p_rhd(ji,jj,iktb) - p_rhd(ji,jj,iktb-1) ) - bco_bc_rhd_bot * p_drhd_k(ji,jj,iktb-1)
         END IF
      END_2D

   END SUBROUTINE calc_drhd_k

!----------------------------------------------------------------------------
 
   SUBROUTINE calc_mid_ccs( k_uv, k_lev_type, p_aco_bc_fld_hor, p_bco_bc_fld_hor, p_fld_pmr, p_fld_mid ) 

      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE calc_mid_ccs  ***
      !!
      !! ** Method  :   Use constrained cubic spline to interpolate 4 profiles to their central point  
      !!                  
      !! ** Action : - set p_fld_mid 
      !!----------------------------------------------------------------------

      INTEGER                               , INTENT(in)  :: k_uv             ! 1 for u-vel; 2 for v_vel
      INTEGER                               , INTENT(in)  :: k_lev_type       ! 1 for t-level data; 2 for w-level data; used with check of land/sea mask
      REAL(wp)                              , INTENT(in)  :: p_aco_bc_fld_hor ! a coeff for horizontal bc (von Neumann or linear extrapolation)
      REAL(wp)                              , INTENT(in)  :: p_bco_bc_fld_hor ! b coeff for horizontal bc (von Neumann or linear extrapolation)
      REAL(wp), DIMENSION(A2D(nn_hls),jpk,4), INTENT(in)  :: p_fld_pmr        ! field in pmr form
      REAL(wp), DIMENSION(A2D(nn_hls),jpk)  , INTENT(out) :: p_fld_mid        ! field at mid-point

      REAL(wp), DIMENSION(A2D(nn_hls),jpk,3) :: zz_dfld_ij    ! constrained spline horizontal differences (dimension 3 for consistency with calc_dfld_cen_dij)

      INTEGER  ji, jj, jk    ! standard loop indices

      INTEGER  ::   j_t_levs, jk_msk          ! indicators that field passed in is valid on t-levels or w-levels, jk level to use with masks
       
         j_t_levs =  1

         DO jk = 1, jpk  

            IF ( k_lev_type == j_t_levs ) THEN 
	       jk_msk = jk 
            ELSE 
               IF ( jk == 1 ) THEN 
                  jk_msk = 1
	       ELSE 
	          jk_msk = jk - 1 
               END IF 
            END IF 
	    
            CALL calc_dfld_pmr_ij( k_uv, k_lev_type, jk, p_aco_bc_fld_hor, p_bco_bc_fld_hor, p_fld_pmr, zz_dfld_ij )  


! first contribution is simple centred average; p_rhd_pmr has 4 points; 2 and 3 are the central ones 
            DO_2D( 0, 0, 0, 0)       
                p_fld_mid(ji,jj,jk) =     0.5_wp * ( p_fld_pmr(ji,jj,jk,2) + p_fld_pmr(ji,jj,jk,3) )  
            END_2D 

            IF ( k_uv == 1 ) THEN 
               DO_2D( 0, 0, 0, 0)       
                  IF ( umask(ji-1, jj, jk_msk) > 0.5 .OR. umask(ji+1, jj, jk_msk) > 0.5 ) THEN
                     p_fld_mid(ji,jj,jk) = p_fld_mid(ji,jj,jk) - 0.125_wp * ( zz_dfld_ij(ji,jj,jk,2) - zz_dfld_ij(ji,jj,jk,1) )
                  END IF 	    
               END_2D 
            ELSE !  k_uv == 2 
               DO_2D( 0, 0, 0, 0)       
                  IF ( vmask(ji, jj-1, jk_msk) > 0.5 .OR. vmask(ji, jj+1, jk_msk) > 0.5 ) THEN
                     p_fld_mid(ji,jj,jk) = p_fld_mid(ji,jj,jk) - 0.125_wp * ( zz_dfld_ij(ji,jj,jk,2) - zz_dfld_ij(ji,jj,jk,1) )
                  END IF 	    
               END_2D 
            END IF ! k_uv  

	 END DO  ! jk 

         IF ( ln_dbg_hpg ) THEN 
            CALL dbg_3dr( 'calc_mid_ccs: zz_dfld_ij: 1', zz_dfld_ij(:,:,:,1) ) 
            CALL dbg_3dr( 'calc_mid_ccs: zz_dfld_ij: 2', zz_dfld_ij(:,:,:,2) ) 
         END IF 


   END SUBROUTINE calc_mid_ccs

!----------------------------------------------------------------------------
 
   SUBROUTINE calc_mid_cub( k_uv, k_lev_type, p_aco_bc_fld_hor, p_bco_bc_fld_hor, p_fld_pmr, p_fld_mid ) 

      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE calc_mid_cub  ***
      !!
      !! ** Method  :   Use simple cubic polynomial to interpolate 4 profiles to their central point  
      !!                The coefficients p_aco_bc_fld_hor, p_bco_bc_fld_hor are not currently used.  
      !!                The simplest form of von Neumann conditions horizontal bc is used (one could off-centred polynomials)   
      !! ** Action : - set p_fld_mid 
      !!----------------------------------------------------------------------

      INTEGER                               , INTENT(in)  :: k_uv             ! 1 for u-vel; 2 for v_vel
      INTEGER                               , INTENT(in)  :: k_lev_type       ! 1 for t-level data; 2 for w-level data; used with check of land/sea mask
      REAL(wp)                              , INTENT(in)  :: p_aco_bc_fld_hor ! a coeff for horizontal bc (von Neumann or linear extrapolation)
      REAL(wp)                              , INTENT(in)  :: p_bco_bc_fld_hor ! b coeff for horizontal bc (von Neumann or linear extrapolation)
      REAL(wp), DIMENSION(A2D(nn_hls),jpk,4), INTENT(inout)  :: p_fld_pmr        ! field in pmr form
      REAL(wp), DIMENSION(A2D(nn_hls),jpk)  , INTENT(out) :: p_fld_mid        ! field at mid-point

      REAL(wp), DIMENSION(A2D(nn_hls))   :: zdfld_21, zdfld_32, zdfld_43  ! primitive horizontal differences 
      REAL(wp), DIMENSION(A2D(nn_hls),3) :: zz_dfld_ij                    ! constrained spline horizontal differences (dimension 3 for consistency with calc_dfld_cen_dij)

      INTEGER  ji, jj, jk       ! standard loop indices
      REAL(wp) z_1_16, z_9_16   ! temporary sum and products
      REAL(wp) z_cor            ! correction to the central value  

      INTEGER  ::   j_t_levs, jk_msk          ! indicators that field passed in is valid on t-levels or w-levels, jk level to use with masks
       
         j_t_levs =  1

         z_1_16 = 1._wp / 16._wp   ;    z_9_16 = 9._wp / 16._wp        

         DO jk = 1, jpk 

            IF ( k_lev_type == j_t_levs ) THEN 
	       jk_msk = jk 
            ELSE 
               IF ( jk == 1 ) THEN 
                  jk_msk = 1
	       ELSE 
	          jk_msk = jk - 1 
               END IF 
            END IF 

! first contribution is simple centred average plus 1/16 of it; p_rhd_pmr has 4 points; 2 and 3 are the central ones 
            DO_2D( 0, 0, 0, 0)       
                p_fld_mid(ji,jj,jk) =  z_9_16 * ( p_fld_pmr(ji,jj,jk,2) + p_fld_pmr(ji,jj,jk,3) )  
            END_2D 

            IF ( k_uv == 1 ) THEN 
               DO_2D( 0, 0, 0, 0)       
                  IF ( umask(ji-1, jj, jk_msk) > 0.5  ) THEN
                     z_cor = p_fld_pmr(ji,jj,jk,1)
                  ELSE 	    
                     z_cor = p_fld_pmr(ji,jj,jk,2)
                  END IF 	    
                  IF ( umask(ji+1, jj, jk_msk) > 0.5  ) THEN
                     z_cor = z_cor + p_fld_pmr(ji,jj,jk,4)
                  ELSE 	    
                     z_cor = z_cor + p_fld_pmr(ji,jj,jk,3)
                  END IF 	    
                  p_fld_mid(ji,jj,jk) = p_fld_mid(ji,jj,jk) - z_1_16 * z_cor
               END_2D 
            ELSE !  k_uv == 2 
               DO_2D( 0, 0, 0, 0)       
                  IF ( vmask(ji, jj-1, jk_msk) > 0.5 ) THEN  
                     z_cor = p_fld_pmr(ji,jj,jk,1)
                  ELSE 	    
                     z_cor = p_fld_pmr(ji,jj,jk,2)
                  END IF 	    
                  IF ( vmask(ji, jj+1, jk_msk) > 0.5  ) THEN
                     z_cor = z_cor + p_fld_pmr(ji,jj,jk,4)
                  ELSE 	    
                     z_cor = z_cor + p_fld_pmr(ji,jj,jk,3)
                  END IF 	    
                  p_fld_mid(ji,jj,jk) = p_fld_mid(ji,jj,jk) - z_1_16 * z_cor
               END_2D 
            END IF ! k_uv  

	 END DO  ! jk 

   END SUBROUTINE calc_mid_cub

!----------------------------------------------------------------------------
 
   SUBROUTINE calc_dz_dij_ccs( k_uv, p_z_pmr, p_dz_ij ) 

      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE calc_dz_dij_ccs  ***
      !!
      !! ** Method  :   Use constrained cubic spline to determine horizontal derivatives of z at 3 central points   
      !!                The routine is limited to z derivatives because the output is valid on w-levels      
      !! ** Action : - set p_dz_ij 
      !!----------------------------------------------------------------------

      INTEGER                               , INTENT(in)  :: k_uv             ! 1 for u-vel; 2 for v_vel
      REAL(wp), DIMENSION(A2D(nn_hls),jpk,4), INTENT(in)  :: p_z_pmr          ! z field in pmr form
      REAL(wp), DIMENSION(A2D(nn_hls),jpk,3), INTENT(out) :: p_dz_ij          ! constrained cubic horizontal derivatives at -1/2, 0 and 1/2

      REAL(wp), DIMENSION(A2D(nn_hls))   :: zdfld_21, zdfld_32, zdfld_43  ! primitive horizontal differences 

      INTEGER  ji, jj, jk    ! standard loop indices
      INTEGER  jk_msk        ! jk level to use for umask and vmask (we need z on the upper and lower faces)  
      INTEGER  j_w_levs      ! indicator for data on w-levels   
      
         j_w_levs = 2 
                   
         DO jk = 1, jpk 

            CALL calc_dfld_pmr_ij( k_uv, j_w_levs, jk, aco_bc_z_hor, bco_bc_z_hor, p_z_pmr, p_dz_ij )  

! copy the 2nd horizontal element to the 3rd element of the output array (for an isolated velocity cell p_dz_ij is the same for all 3 elements) 
            DO_2D( 0, 0, 0, 0)       
                p_dz_ij(ji,jj,jk,3) = p_dz_ij(ji,jj,jk,2)      
            END_2D 

! set the central element if it is not an isolated velocity cell 
            IF ( jk == 1 ) THEN 
               jk_msk = 1
	    ELSE 
	       jk_msk = jk - 1 
            END IF 

            IF ( k_uv == 1 ) THEN 
               DO_2D( 0, 0, 0, 0)       
                  IF ( umask(ji-1, jj, jk_msk) > 0.5 .OR. umask(ji+1, jj, jk_msk) > 0.5 ) THEN
                     p_dz_ij(ji,jj,jk,2) = 1.5_wp*( p_z_pmr(ji,jj,jk,3) - p_z_pmr(ji,jj,jk,2) ) - 0.25_wp * ( p_dz_ij(ji,jj,jk,3) + p_dz_ij(ji,jj,jk,1) )
                  END IF 	    
               END_2D 
            ELSE !  k_uv == 2 
               DO_2D( 0, 0, 0, 0)       
                  IF ( vmask(ji, jj-1, jk_msk) > 0.5 .OR. vmask(ji, jj+1, jk_msk) > 0.5 ) THEN
                     p_dz_ij(ji,jj,jk,2) = 1.5_wp*( p_z_pmr(ji,jj,jk,3) - p_z_pmr(ji,jj,jk,2) ) - 0.25_wp * ( p_dz_ij(ji,jj,jk,3) + p_dz_ij(ji,jj,jk,1) )
                  END IF 	    
               END_2D 
            END IF ! k_uv  
         
	 END DO  ! jk 

   END SUBROUTINE calc_dz_dij_ccs

!----------------------------------------------------------------------------
 
   SUBROUTINE calc_dz_dij_cub( k_uv, p_z_pmr, p_dz_ij ) 

      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE calc_dz_dij_cub  ***
      !!
      !! ** Method  :   Use simple cubic polynomial to determine horizontal derivatives at 3 central points   
      !!                based on equations (5.5) and (5.6) of SMcW 2003 
      !!                The routine is limited to z derivatives because the output is valid on w-levels      
      !! ** Action : - set p_dz_ij 
      !!----------------------------------------------------------------------

      INTEGER                               , INTENT(in)  :: k_uv             ! 1 for u-vel; 2 for v_vel
      REAL(wp), DIMENSION(A2D(nn_hls),jpk,4), INTENT(in)  :: p_z_pmr          ! z field in pmr form
      REAL(wp), DIMENSION(A2D(nn_hls),jpk,3), INTENT(out) :: p_dz_ij          ! horizontal derivatives at -1/2, 0 and 1/2

      REAL(wp), DIMENSION(A2D(nn_hls)) :: z_z_a, z_z_d       ! z fields with boundary conditions applied 
      REAL(wp) z_1_24                                        ! 1/24 
      REAL(wp) z_a, z_b, z_c, z_d                            ! values of input field at the four points used (-3/2, -1/2, 1/2, 3/2) 
      REAL(wp) z_c_m_b, z_d_m_a                              ! differences c - b   and d - a  
      REAL(wp) z_co_1, z_co_2, z_co_3                        ! coefficients of polynomial (field = z_co_0 + z_co_1 zeta + .. ) 
      INTEGER  ji, jj, jk                                    ! standard loop indices
      INTEGER  jk_msk                                        ! jk level to use for umask and vmask (we need z on the upper and lower faces)  

         z_1_24 = 1.0_wp / 24.0_wp 
       
         DO jk = 1, jpk 

!------------------------------------------------------
! 1. Use simple von Neumann conditions at the boundaries  
!------------------------------------------------------
            IF ( jk == 1 ) THEN 
               jk_msk = 1
	    ELSE 
	       jk_msk = jk - 1 
            END IF 

            IF ( k_uv == 1 ) THEN 
               DO_2D( 0, 0, 0, 0)       
                  IF ( umask(ji-1, jj, jk_msk) > 0.5  ) THEN
                     z_z_a(ji,jj) = p_z_pmr(ji,jj,jk,1)
                  ELSE 	    
                     z_z_a(ji,jj) = p_z_pmr(ji,jj,jk,2)
                  END IF 	    
                  IF ( umask(ji+1, jj, jk_msk) > 0.5  ) THEN
                     z_z_d(ji,jj) = p_z_pmr(ji,jj,jk,4)
                  ELSE 	    
                     z_z_d(ji,jj) = p_z_pmr(ji,jj,jk,3)
                  END IF 	    
               END_2D 
            ELSE !  k_uv == 2 
               DO_2D( 0, 0, 0, 0)       
                  IF ( vmask(ji, jj-1, jk_msk) > 0.5 ) THEN  
                     z_z_a(ji,jj) = p_z_pmr(ji,jj,jk,1)
                  ELSE 	    
                     z_z_a(ji,jj) = p_z_pmr(ji,jj,jk,2)
                  END IF 	    
                  IF ( vmask(ji, jj+1, jk_msk) > 0.5  ) THEN
                     z_z_d(ji,jj) = p_z_pmr(ji,jj,jk,4)
                  ELSE 	    
                     z_z_d(ji,jj) = p_z_pmr(ji,jj,jk,3)
                  END IF 	    
               END_2D 
            END IF ! k_uv  

!------------------------------------------------------
! 2. calculate the coefficients for the polynomial fit (c_1, c_2 and c_3; we don't need c_0)  
!------------------------------------------------------

            DO_2D( 0, 0, 0, 0)       
                z_a = z_z_a(ji,jj)      
                z_b = p_z_pmr(ji,jj,jk,2)
                z_c = p_z_pmr(ji,jj,jk,3) 
                z_d = z_z_d(ji,jj)      

                z_c_m_b = z_c - z_b    ;   z_d_m_a = z_d - z_a 

! eqn (5.6) of SMcW
                z_co_1 = 1.125_wp *  z_c_m_b  - z_1_24 *  z_d_m_a 
                z_co_2 =  -0.5_wp * (z_c+z_b) + 0.5_wp * (z_d+z_a) 
                z_co_3 =  -3.0_wp *  z_c_m_b  +           z_d_m_a 

! eqn (5.5) of SMcW (first derivative of it) 
                p_dz_ij(ji,jj,jk,1) = z_co_1 - 0.5_wp*z_co_2 + 0.125_wp*z_co_3   ! zeta = -0.5
                p_dz_ij(ji,jj,jk,2) = z_co_1                                     ! zeta =  0.0
                p_dz_ij(ji,jj,jk,3) = z_co_1 + 0.5_wp*z_co_2 + 0.125_wp*z_co_3   ! zeta =  0.5

            END_2D 

	 END DO  ! jk 

   END SUBROUTINE calc_dz_dij_cub

!----------------------------------------------------------------------------
 
   SUBROUTINE calc_dfld_pmr_ij( k_uv, k_lev_type, kk, p_aco_bc_fld_hor, p_bco_bc_fld_hor, p_fld_pmr, p_dfld_ij ) 

      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE calc_dfld_pmr_ij  ***
      !!
      !! ** Method  :   Calculate constrained spline horizontal derivatives   
      !!                  
      !! ** Action : - set p_dfld_ij 
      !!----------------------------------------------------------------------

      INTEGER                               , INTENT(in)  :: k_uv             ! 1 for u-vel; 2 for v_vel
      INTEGER                               , INTENT(in)  :: k_lev_type       ! 1 for t-level data; 2 for w-level data; used with check of land/sea mask
      INTEGER                               , INTENT(in)  :: kk               ! vertical level
      REAL(wp)                              , INTENT(in)  :: p_aco_bc_fld_hor ! a coeff for horizontal bc (von Neumann or linear extrapolation)
      REAL(wp)                              , INTENT(in)  :: p_bco_bc_fld_hor ! b coeff for horizontal bc (von Neumann or linear extrapolation)
      REAL(wp), DIMENSION(A2D(nn_hls),jpk,4), INTENT(in)  :: p_fld_pmr        ! field in pmr form
      REAL(wp), DIMENSION(A2D(nn_hls),jpk,3), INTENT(out) :: p_dfld_ij        ! constrained spline horizontal derivatives of the field; only 1 and 2 are set!  

      REAL(wp), DIMENSION(A2D(nn_hls))   :: zdfld_21, zdfld_32, zdfld_43  ! primitive horizontal differences 

      INTEGER  ji, jj        ! standard loop indices
      INTEGER  jio, jjo      ! offsets 
      INTEGER  iktb          ! index of the bottom of ref profile

      REAL(wp) z1_cff, z_cff_31, z_cff_42  ! temporary sum and products
      REAL(wp) zep

      INTEGER  ::   j_t_levs, jk_msk     ! indicators that field passed in is valid on t-levels or w-levels; jk to use with masks 
       
         j_t_levs =  1
       
      !----------------------------------------------------------------------------------------
      !  1. compute and store elementary horizontal differences zfor z_rhd_pmr arrays 
      !----------------------------------------------------------------------------------------
 
         IF ( k_uv == 1) THEN 
            jio = 1
            jjo = 0 
         ELSE 
            jio = 0 
            jjo = 1
         END IF 

         IF ( k_lev_type == j_t_levs ) THEN 
	    jk_msk = kk 
         ELSE 
            IF ( kk == 1 ) THEN 
               jk_msk = 1
	    ELSE 
	       jk_msk = kk - 1 
            END IF 
         END IF

         DO_2D( 0, 0, 0, 0 )
            zdfld_21(ji,jj) =   p_fld_pmr(ji,jj,kk,2) - p_fld_pmr(ji,jj,kk,1)
            zdfld_32(ji,jj) =   p_fld_pmr(ji,jj,kk,3) - p_fld_pmr(ji,jj,kk,2)
            zdfld_43(ji,jj) =   p_fld_pmr(ji,jj,kk,4) - p_fld_pmr(ji,jj,kk,3)
         END_2D

         zep = 1.e-15
         DO_2D( 0, 0, 0, 0 )
            z_cff_31 = MAX( 2._wp * zdfld_21(ji,jj) * zdfld_32(ji,jj), 0._wp ) 
            z1_cff = zdfld_21(ji,jj) + zdfld_32(ji,jj)
            p_dfld_ij(ji,jj,kk,1) = z_cff_31 / SIGN( MAX( ABS(z1_cff), zep ), z1_cff )

            z_cff_42 = MAX( 2._wp * zdfld_32(ji,jj) * zdfld_43(ji,jj), 0._wp ) 
            z1_cff = zdfld_32(ji,jj) + zdfld_43(ji,jj)
            p_dfld_ij(ji,jj,kk,2) = z_cff_42 / SIGN( MAX( ABS(z1_cff), zep ), z1_cff ) 
         END_2D

      !----------------------------------------------------------------------------------
      ! 2. apply boundary conditions at side boundaries using 5.36-5.37    
      !----------------------------------------------------------------------------------

         IF ( k_uv == 1 ) THEN 
            DO_2D( 0, 0, 0, 0)       
            ! Walls coming from left: should check from 2 to jpi-1 (and jpj=2-jpj)
               IF ( umask(ji,jj,jk_msk) > 0.5_wp .AND. umask(ji-1,jj,jk_msk) < 0.5_wp .AND. umask(ji+1,jj,jk_msk) > 0.5_wp)  THEN  
                  p_dfld_ij(ji,jj,kk,1) = p_aco_bc_fld_hor * ( p_fld_pmr(ji,jj,kk,3) - p_fld_pmr(ji,jj,kk,2) ) - p_bco_bc_fld_hor * p_dfld_ij(ji,jj,kk,2) 
               END IF
            ! Walls coming from right: should check from 3 to jpi (and jpj=2-jpj)
               IF ( umask(ji,jj,jk_msk) < 0.5_wp .AND. umask(ji-1,jj,jk_msk) > 0.5_wp .AND. umask(ji-2,jj,jk_msk) > 0.5_wp) THEN
                  p_dfld_ij(ji,jj,kk,2) = p_aco_bc_fld_hor * ( p_fld_pmr(ji,jj,kk,3) - p_fld_pmr(ji,jj,kk,2) ) - p_bco_bc_fld_hor * p_dfld_ij(ji,jj,kk,1)  
               END IF
            END_2D 

            DO_2D( 0, 0, 0, 0)       
            ! For an isolated velocity point use the central difference (this is important) 
               IF ( umask(ji-1, jj, jk_msk) < 0.5 .AND. umask(ji+1, jj, jk_msk) < 0.5 ) THEN
                     p_dfld_ij(ji,jj,kk,1) = p_fld_pmr(ji,jj,kk,3) - p_fld_pmr(ji,jj,kk,2)
                     p_dfld_ij(ji,jj,kk,2) = p_dfld_ij(ji,jj,kk,1)
               END IF 	    
            END_2D 

         ELSE !  k_uv == 2 

            DO_2D( 0, 0, 0, 0)       
            ! Walls coming from left: should check from 2 to jpj-1 (and jpi=2-jpi)
               IF ( vmask(ji,jj,jk_msk) > 0.5_wp .AND. vmask(ji,jj-1,jk_msk) < 0.5_wp .AND. vmask(ji,jj+1,jk_msk) > 0.5_wp)  THEN
                  p_dfld_ij(ji,jj,kk,1) = p_aco_bc_fld_hor * (p_fld_pmr(ji,jj,kk,3) - p_fld_pmr(ji,jj,kk,2) ) - p_bco_bc_fld_hor * p_dfld_ij(ji,jj,kk,2)
               END IF 
            ! Walls coming from right: should check from 3 to jpj (and jpi=2-jpi)
               IF ( vmask(ji,jj,jk_msk) < 0.5_wp .AND. vmask(ji,jj-1,jk_msk) > 0.5_wp .AND. vmask(ji,jj-2,jk_msk) > 0.5_wp) THEN 
                  p_dfld_ij(ji,jj,kk,2) = p_aco_bc_fld_hor * (p_fld_pmr(ji,jj,kk,3) - p_fld_pmr(ji,jj,kk,2) ) - p_bco_bc_fld_hor * p_dfld_ij(ji,jj,kk,1)
               END IF
            END_2D 

            DO_2D( 0, 0, 0, 0)       
            ! For an isolated velocity point use the central difference (this is important) 
               IF ( vmask(ji, jj-1, jk_msk) < 0.5 .AND. vmask(ji, jj+1, jk_msk) < 0.5 ) THEN
                     p_dfld_ij(ji,jj,kk,1) = p_fld_pmr(ji,jj,kk,3) - p_fld_pmr(ji,jj,kk,2)
                     p_dfld_ij(ji,jj,kk,2) = p_dfld_ij(ji,jj,kk,1)
               END IF 	    
            END_2D 

         END IF ! k_uv  

   END SUBROUTINE calc_dfld_pmr_ij

!------------------------------------------------------------------------------------------

   SUBROUTINE vrt_int_quad(k_mbkt, p_rhd, p_dept, p_depw, p_p, p_F) 

      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE calc_rhd_k  ***
      !!
      !! ** Method  :   Use quadratic reconstruction of p_rhd (density profile) to calculate pressure profile and 
      !!                forces on vertical faces between w levels. Formulated for variable vertical grid spacing.    
      !!                  
      !! ** Action : - set p_p and p_F
      !!----------------------------------------------------------------------

      INTEGER,  DIMENSION(A2D(nn_hls)),     INTENT(in)  :: k_mbkt   ! bottom level
      REAL(wp), DIMENSION(A2D(nn_hls),jpk), INTENT(in)  :: p_rhd    ! densities on tracer levels
      REAL(wp), DIMENSION(A2D(nn_hls),jpk), INTENT(in)  :: p_dept   ! depths of T points 
      REAL(wp), DIMENSION(A2D(nn_hls),jpk), INTENT(in)  :: p_depw   ! depths of w points (top & bottom of T cells)  
      REAL(wp), DIMENSION(A2D(nn_hls),jpk), INTENT(out) :: p_p      ! pressures on w levels
      REAL(wp), DIMENSION(A2D(nn_hls),jpk), INTENT(out) :: p_F      ! force on the vertical face between w levels

      INTEGER  ::   ji, jj, jk                  ! dummy loop indices
      REAL(wp) ::   z_1_3                       ! 1/3  

      INTEGER, DIMENSION(A2D(nn_hls))  :: jklev
      REAL(wp), DIMENSION(A2D(nn_hls)) :: zr_a, zr_0, zr_b   ! rhd values (r) at points a (above), 0 (central), b (below) 
      REAL(wp), DIMENSION(A2D(nn_hls)) :: zz_a, zz_0, zz_b   ! heights at points a (above), 0 (central), b (below) 
      REAL(wp), DIMENSION(A2D(nn_hls)) :: zz_upp, zz_low     ! upper & lower w levels 
      INTEGER  ::   jkt                                 ! level being worked on 
      REAL(wp) ::   zet_a, zet_0, zet_b                 ! zeta values for above (a), central (0) and below (b) 
      REAL(wp) ::   zr_d_a, zr_d_0, zr_d_b              ! reciprocals of denominators for terms involving r_a, r_0 and r_b 
      REAL(wp) ::   za_0, za_1, za_2                    ! quadratic polynomial coefficients for r (rhd)  
      REAL(wp) ::   zet_upp, zet_low                    ! zeta for upper and lower w levels  
      REAL(wp) ::   zet_upp_sq, zet_low_sq              ! squares of above values 
      REAL(wp) ::   zp_upp, zp_low                      ! contributions to the pressure 
      REAL(wp) ::   zf_uml                              ! force "upper minus lower" 
      
      z_1_3  = 1.0_wp/3.0_wp

! Integrate densities in the vertical using quadratic fits to the densities 
! zr_a is r_{-1}  zr_b is r_1 ; zr_c is r_3  

      DO_2D(0,0,0,0) 
         p_p(ji,jj,1) = 0.0_wp
      END_2D

! one-sided quadratic at the upper boundary. pressure is zero at the upper surface
      DO jk = 1, jpk 
      
         IF (jk == 1) THEN        !  perform off-centred calculations for the first level 
            DO_2D(0,0,0,0) 
               jklev(ji,jj)  = jk 
	       zr_a(ji,jj)   =   p_rhd(ji,jj,1) 
	       zr_0(ji,jj)   =   p_rhd(ji,jj,2) 
	       zr_b(ji,jj)   =   p_rhd(ji,jj,3) 
	       zz_a(ji,jj)   = - p_dept(ji,jj,1)
	       zz_0(ji,jj)   = - p_dept(ji,jj,2)
               zz_b(ji,jj)   = - p_dept(ji,jj,3)
               zz_upp(ji,jj) = - p_depw(ji,jj,1) 
	       zz_low(ji,jj) = - p_depw(ji,jj,2) 
            END_2D

         ELSE IF ( jk < jpk ) THEN 
            DO_2D(0,0,0,0) 
               jklev(ji,jj)  = jk 
	       zr_a(ji,jj)   =   p_rhd(ji,jj,jk-1) 
	       zr_0(ji,jj)   =   p_rhd(ji,jj,jk) 
	       zr_b(ji,jj)   =   p_rhd(ji,jj,jk+1) 
	       zz_a(ji,jj)   = - p_dept(ji,jj,jk-1)
	       zz_0(ji,jj)   = - p_dept(ji,jj,jk)
               zz_b(ji,jj)   = - p_dept(ji,jj,jk+1)
               zz_upp(ji,jj) = - p_depw(ji,jj,jk) 
	       zz_low(ji,jj) = - p_depw(ji,jj,jk+1) 
            END_2D

         ELSE IF ( jk == jpk ) THEN    ! repeat calculations at bottom level using off-centred quadratic

            DO_2D(0,0,0,0) 
               jkt          = k_mbkt(ji,jj) 
               jklev(ji,jj) = jkt 
	       zr_a(ji,jj)   =   p_rhd(ji,jj,jkt-2) 
	       zr_0(ji,jj)   =   p_rhd(ji,jj,jkt-1) 
	       zr_b(ji,jj)   =   p_rhd(ji,jj,jkt) 
	       zz_a(ji,jj)   = - p_dept(ji,jj,jkt-2)
	       zz_0(ji,jj)   = - p_dept(ji,jj,jkt-1)
               zz_b(ji,jj)   = - p_dept(ji,jj,jkt)
               zz_upp(ji,jj) = - p_depw(ji,jj,jkt) 
	       zz_low(ji,jj) = - p_depw(ji,jj,jkt+1) 
            END_2D

         END IF  ! jk 
	 
         DO_2D(0,0,0,0) 
            jkt  = jklev(ji,jj) 

            zet_a = - ( zz_a(ji,jj) - zz_0(ji,jj) )  
            zet_0 = 0.0_wp
	    zet_b = - ( zz_b(ji,jj) - zz_0(ji,jj) )

            zr_d_a = 1.0_wp / ( zet_a*(zet_a - zet_b) ) 
            zr_d_0 = 1.0_wp / ( zet_a*zet_b ) 
	    zr_d_b = 1.0_wp / ( zet_b*(zet_b - zet_a) )

            za_0 = zr_0(ji,jj)  
            za_1 = - zr_d_0*zr_0(ji,jj)*(zet_a+zet_b) - zr_d_a*zr_a(ji,jj)*zet_b - zr_d_b*zr_b(ji,jj)*zet_a  
            za_2 =   zr_d_0*zr_0(ji,jj)               + zr_d_a*zr_a(ji,jj)       + zr_d_b*zr_b(ji,jj) 
	    
            zet_upp = - ( zz_upp(ji,jj) - zz_0(ji,jj) ) 
            zet_low = - ( zz_low(ji,jj) - zz_0(ji,jj) ) 

            zet_upp_sq = zet_upp*zet_upp   ;  zet_low_sq = zet_low*zet_low
	    
	    zp_upp = za_0*zet_upp + 0.5_wp*za_1*zet_upp_sq + z_1_3*za_2*zet_upp_sq*zet_upp  
            zp_low = za_0*zet_low + 0.5_wp*za_1*zet_low_sq + z_1_3*za_2*zet_low_sq*zet_low 

            p_p(ji,jj,jkt+1) = p_p(ji,jj,jkt) + grav * ( zp_low - zp_upp ) 

            zf_uml =          za_0 * ( 0.5_wp *(zet_low_sq            - zet_upp_sq)            - zet_upp           *(zet_low-zet_upp) ) &
            &        + 0.5_wp*za_1 * ( z_1_3  *(zet_low_sq*zet_low    - zet_upp_sq*zet_upp)    - zet_upp_sq        *(zet_low-zet_upp) ) &
            &        + z_1_3 *za_2 * ( 0.25_wp*(zet_low_sq*zet_low_sq - zet_upp_sq*zet_upp_sq) - zet_upp_sq*zet_upp*(zet_low-zet_upp) ) 
 
! the force on the face is positive (the opposite sign from my notes) 
            p_F(ji,jj,jkt) = p_p(ji,jj,jkt)*( zz_upp(ji,jj) - zz_low(ji,jj) ) + grav * zf_uml   

         END_2D 

      END DO  ! jk  

   END SUBROUTINE vrt_int_quad
   
!------------------------------------------------------------------------------------------

   SUBROUTINE dbg_2di( cc_array, ki_2d ) 
   CHARACTER*(*),                   INTENT(IN) :: cc_array
   INTEGER, DIMENSION(A2D(nn_hls)), INTENT(IN) :: ki_2d   

   INTEGER ::  ji_prt_low, ji_prt_upp, jj_prt_low, jj_prt_upp 
   INTEGER ::  ji_prt, jj_prt

   IF( lwp ) THEN
      ji_prt_low = MAX( ntsi-nn_hls, ki_dbg_min )     
      ji_prt_upp = MIN( ntei+nn_hls, ki_dbg_max )     
      jj_prt_low = MAX( ntsj-nn_hls, kj_dbg_min )     
      jj_prt_upp = MIN( ntej+nn_hls, kj_dbg_max )     

!     print out a 2D horizontal slice 

      IF ( ji_prt_low <= ji_prt_upp .AND. ln_dbg_ij ) THEN 
         IF ( jj_prt_low <= jj_prt_upp ) THEN 
            WRITE(numout,*)
            WRITE(numout,*) cc_array 
            WRITE(numout,*) ' ji_prt_low, ji_prt_upp = ', ji_prt_low, ji_prt_upp
            WRITE(numout,*) ' row/col ', ( ji_prt, ji_prt = ji_prt_low, ji_prt_upp )
            DO jj_prt = jj_prt_low, jj_prt_upp
	      WRITE(numout,*) jj_prt, ( ki_2d(ji_prt, jj_prt), ji_prt = ji_prt_low, ji_prt_upp ) 
            END DO 
         END IF 
      END IF 

   END IF 

   RETURN    
   END SUBROUTINE dbg_2di

!----------------------------------------------------------------------------

   SUBROUTINE dbg_2di_k( cc_array, ki_2d, kk ) 
   CHARACTER*(*),                   INTENT(IN) :: cc_array
   INTEGER, DIMENSION(A2D(nn_hls)), INTENT(IN) :: ki_2d   
   INTEGER,                         INTENT(IN) :: kk   

   INTEGER ::  ji_prt_low, ji_prt_upp, jj_prt_low, jj_prt_upp 
   INTEGER ::  ji_prt, jj_prt, jk_prt

   IF( lwp ) THEN
      ji_prt_low = MAX( ntsi, ki_dbg_min )     
      ji_prt_upp = MIN( ntei, ki_dbg_max )     
      jj_prt_low = MAX( ntsj, kj_dbg_min )     
      jj_prt_upp = MIN( ntej, kj_dbg_max )     

!     print out a 2D (ji, jk) slice 

      IF ( ji_prt_low <= ji_prt_upp .AND. ln_dbg_ik ) THEN 
         IF ( ntsj <= kj_dbg_cen .AND. kj_dbg_cen <= ntej ) THEN 
            IF ( kk == kk_dbg_min ) THEN  
               WRITE(numout,*)
               WRITE(numout,*) cc_array, ' kj = ', kj_dbg_cen
               WRITE(numout,*) 'ji_prt_low, ji_prt_upp = ', ji_prt_low, ji_prt_upp
               WRITE(numout,*) 'kk_dbg_min, kk_dbg_max = ',  kk_dbg_min, kk_dbg_max 
            END IF 
            IF ( kk_dbg_min <= kk .AND. kk <= kk_dbg_max )   & 
     &         WRITE(numout,*) cc_array, kk, ( ki_2d(ji_prt, kj_dbg_cen), ji_prt = ji_prt_low, ji_prt_upp ) 
         END IF 
      END IF 

!     print out a 2D (jj, jk) slice 

      IF ( jj_prt_low <= jj_prt_upp .AND. ln_dbg_jk ) THEN 
         IF ( ntsi <= ki_dbg_cen .AND. ki_dbg_cen <= ntei ) THEN 
            IF ( kk == kk_dbg_min ) THEN  
               WRITE(numout,*)
               WRITE(numout,*) cc_array, ' ki = ', ki_dbg_cen
               WRITE(numout,*) 'jj_prt_low, jj_prt_upp = ', jj_prt_low, jj_prt_upp
               WRITE(numout,*) 'kk_dbg_min, kk_dbg_max = ',  kk_dbg_min, kk_dbg_max 
            END IF 
            IF ( kk_dbg_min <= kk .AND. kk <= kk_dbg_max )   & 
     &         WRITE(numout,*) cc_array, kk, ( ki_2d(ki_dbg_cen, jj_prt), jj_prt = jj_prt_low, jj_prt_upp ) 
         END IF 
      END IF 

   END IF 

   RETURN    
   END SUBROUTINE dbg_2di_k

!----------------------------------------------------------------------------

   SUBROUTINE dbg_2dr( cc_array, pr_2d ) 
   CHARACTER*(*),                    INTENT(IN) :: cc_array
   REAL(wp), DIMENSION(A2D(nn_hls)), INTENT(IN) :: pr_2d   

   INTEGER ::  ji_prt_low, ji_prt_upp, jj_prt_low, jj_prt_upp 
   INTEGER ::  ji_prt, jj_prt

   IF(lwp) THEN
      ji_prt_low = MAX( ntsi-nn_hls, ki_dbg_min )     
      ji_prt_upp = MIN( ntei+nn_hls, ki_dbg_max )     
      jj_prt_low = MAX( ntsj-nn_hls, kj_dbg_min )     
      jj_prt_upp = MIN( ntej+nn_hls, kj_dbg_max )     

!     print out a 2D horizontal slice 

      IF ( ji_prt_low <= ji_prt_upp .AND. ln_dbg_ij ) THEN 
         IF ( jj_prt_low <= jj_prt_upp ) THEN 
            WRITE(numout,*)
            WRITE(numout,*) cc_array 
            WRITE(numout,*) ' row/col ', ( ji_prt, ji_prt = ji_prt_low, ji_prt_upp )
            DO jj_prt = jj_prt_low, jj_prt_upp
	      WRITE(numout,*) jj_prt, ( pr_2d(ji_prt, jj_prt), ji_prt = ji_prt_low, ji_prt_upp ) 
            END DO 
         END IF 
      END IF 

   END IF 

   RETURN    
   END SUBROUTINE dbg_2dr

!----------------------------------------------------------------------------

   SUBROUTINE dbg_3di( cc_array, ki_3d )
    
   CHARACTER*(*),                       INTENT(IN) :: cc_array
   INTEGER, DIMENSION(A2D(nn_hls),jpk), INTENT(IN) :: ki_3d   

   INTEGER ::  ji_prt_low, ji_prt_upp, jj_prt_low, jj_prt_upp 
   INTEGER ::  ji_prt, jj_prt, jk_prt

! output values 


   IF(lwp) THEN
      ji_prt_low = MAX( ntsi-nn_hls, ki_dbg_min )     
      ji_prt_upp = MIN( ntei+nn_hls, ki_dbg_max )     
      jj_prt_low = MAX( ntsj-nn_hls, kj_dbg_min )     
      jj_prt_upp = MIN( ntej+nn_hls, kj_dbg_max )     

!     print out a 2D horizontal slice 

      IF ( ji_prt_low <= ji_prt_upp .AND. ln_dbg_ij ) THEN 
         IF ( jj_prt_low <= jj_prt_upp ) THEN 
            WRITE(numout,*)
            WRITE(numout,*) cc_array, ' level = ', kk_dbg_cen
            WRITE(numout,*) ' row/col ', ( ji_prt, ji_prt = ji_prt_low, ji_prt_upp )
            DO jj_prt = jj_prt_low, jj_prt_upp
	      WRITE(numout,*) jj_prt, ( ki_3d(ji_prt, jj_prt, kk_dbg_cen), ji_prt = ji_prt_low, ji_prt_upp ) 
            END DO 
         END IF 
      END IF 

!     print out a 2D (ji, jk) slice 

      IF ( ji_prt_low <= ji_prt_upp .AND. ln_dbg_ik ) THEN 
         IF ( ntsj <= kj_dbg_cen .AND. kj_dbg_cen <= ntej ) THEN 
            WRITE(numout,*)
            WRITE(numout,*) cc_array, ' kj = ', kj_dbg_cen
            WRITE(numout,*) ' row/lev ', ( ji_prt, ji_prt = ji_prt_low, ji_prt_upp )
            DO jk_prt = kk_dbg_min, kk_dbg_max
	      WRITE(numout,*) jk_prt, ( ki_3d(ji_prt, kj_dbg_cen, jk_prt), ji_prt = ji_prt_low, ji_prt_upp ) 
            END DO 
         END IF 
      END IF 

!     print out a 2D (jj, jk) slice 

      IF ( jj_prt_low <= jj_prt_upp .AND. ln_dbg_jk ) THEN 
         IF ( ntsi <= ki_dbg_cen .AND. ki_dbg_cen <= ntei ) THEN 
            WRITE(numout,*)
            WRITE(numout,*) cc_array, ' ki = ', ki_dbg_cen
            WRITE(numout,*) ' col/lev ', ( jj_prt, jj_prt = jj_prt_low, jj_prt_upp )
            DO jk_prt = kk_dbg_min, kk_dbg_max
	      WRITE(numout,*) jk_prt, ( ki_3d(ki_dbg_cen, jj_prt, jk_prt), jj_prt = jj_prt_low, jj_prt_upp ) 
            END DO 
         END IF 
      END IF 

   END IF 

   RETURN    
   END SUBROUTINE dbg_3di

!----------------------------------------------------------------------------

   SUBROUTINE dbg_3dr( cc_array, pr_3d )
    
   CHARACTER*(*),                        INTENT(IN) :: cc_array
   REAL(wp), DIMENSION(A2D(nn_hls),jpk), INTENT(IN) :: pr_3d   

   INTEGER ::  ji_prt_low, ji_prt_upp, jj_prt_low, jj_prt_upp 
   INTEGER ::  ji_prt, jj_prt, jk_prt

! output values 


   IF(lwp) THEN
      ji_prt_low = MAX( ntsi-nn_hls, ki_dbg_min )     
      ji_prt_upp = MIN( ntei+nn_hls, ki_dbg_max )     
      jj_prt_low = MAX( ntsj-nn_hls, kj_dbg_min )     
      jj_prt_upp = MIN( ntej+nn_hls, kj_dbg_max )     

!     print out a 2D horizontal slice 

      IF ( ji_prt_low <= ji_prt_upp .AND. ln_dbg_ij ) THEN 
         IF ( jj_prt_low <= jj_prt_upp ) THEN 
            WRITE(numout,*)
            WRITE(numout,*) cc_array, ' level = ', kk_dbg_cen
            WRITE(numout,*) ' row/col ', ( ji_prt, ji_prt = ji_prt_low, ji_prt_upp )
            DO jj_prt = jj_prt_low, jj_prt_upp
	      WRITE(numout,*) jj_prt, ( pr_3d(ji_prt, jj_prt, kk_dbg_cen), ji_prt = ji_prt_low, ji_prt_upp ) 
            END DO 
         END IF 
      END IF 

!     print out a 2D (ji, jk) slice 

      IF ( ji_prt_low <= ji_prt_upp .AND. ln_dbg_ik ) THEN 
         IF ( ntsj <= kj_dbg_cen .AND. kj_dbg_cen <= ntej ) THEN 
            WRITE(numout,*)
            WRITE(numout,*) cc_array, ' kj = ', kj_dbg_cen
            WRITE(numout,*) ' row/lev ', ( ji_prt, ji_prt = ji_prt_low, ji_prt_upp )
            DO jk_prt = kk_dbg_min, kk_dbg_max
	      WRITE(numout,*) jk_prt, ( pr_3d(ji_prt, kj_dbg_cen, jk_prt), ji_prt = ji_prt_low, ji_prt_upp ) 
            END DO 
         END IF 
      END IF 

!     print out a 2D (jj, jk) slice 

      IF ( jj_prt_low <= jj_prt_upp .AND. ln_dbg_jk ) THEN 
         IF ( ntsi <= ki_dbg_cen .AND. ki_dbg_cen <= ntei ) THEN 
            WRITE(numout,*)
            WRITE(numout,*) cc_array, ' ki = ', ki_dbg_cen
            WRITE(numout,*) ' col/lev ', ( jj_prt, jj_prt = jj_prt_low, jj_prt_upp )
            DO jk_prt = kk_dbg_min, kk_dbg_max
	      WRITE(numout,*) jk_prt, ( pr_3d(ki_dbg_cen, jj_prt, jk_prt), jj_prt = jj_prt_low, jj_prt_upp ) 
            END DO 
         END IF 
      END IF 

   END IF 

   RETURN    
   END SUBROUTINE dbg_3dr


   !!======================================================================
END MODULE dynhpg
