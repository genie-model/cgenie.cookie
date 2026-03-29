! ******************************************************************************************************************************** !
! sedgem_data.f90
! DATA LOADING/SAVING ROUTINES
! ******************************************************************************************************************************** !


MODULE sedgem_data


  use genie_control
  USE gem_cmn
  USE gem_util
  USE gem_netcdf
  USE sedgem_lib
  USE sedgem_box
  USE sedgem_data_netCDF
  USE sedgem_nnutils
  IMPLICIT NONE
  SAVE


CONTAINS


  ! ****************************************************************************************************************************** !
  ! DATA LOADING ROUTINES
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  ! LOAD SEDGEM 'goin' FILE OPTIONS
  SUBROUTINE sub_load_goin_sedgem()
    USE genie_util, ONLY: check_unit, check_iostat
    ! local variables
    integer::ios
    ! read data_SEDGEM file
    call check_unit(in,__LINE__,__FILE__)
    open(unit=in,file='data_SEDGEM',status='old',action='read',iostat=ios)
    if (ios /= 0) then
       print*,'ERROR: could not open SEDGEM initialisation namelist file'
       stop
    end if
    ! read in namelist and close data_SEDGEM file
    read(UNIT=in,NML=ini_sedgem_nml,IOSTAT=ios)
    if (ios /= 0) then
       print*,'ERROR: could not read SEDGEM namelist'
       stop
    else
       close(unit=in)
    end if
    ! set and report namelist data
    par_indir_name = trim(par_indir_name)//'/'
    par_outdir_name = trim(par_outdir_name)//'/'
    par_inrstdir_name = trim(par_inrstdir_name)//'/'
    par_outrstdir_name = trim(par_outrstdir_name)//'/'
    if (ctrl_debug_init > 0) then
       ! --- RUN CONTROL --------------------------------------------------------------------------------------------------------- !
       print*,'--- RUN CONTROL ------------------------------------'
       print*,'Continuing run?                                     : ',ctrl_continuing
       print*,'Simulation start year [REAL]                        : ',start_year
       print*,'Simulation run length (yr)                          : ',par_misc_t_runtime
       ! --- PHYSICAL CONFIGURATION ---------------------------------------------------------------------------------------------- !
       print*,'--- PHYSICAL CONFIGURATION -------------------------'
       print*,'Top (well-mixed) sediment layer thickness (cm)      : ',par_sed_top_th
       print*,'Sub-surface detrital porosity (cm3 cm-3)            : ',par_sed_poros_det
       print*,'Sub-surface carbonate porosity (cm3 cm-3)           : ',par_sed_poros_CaCO3
       print*,'Maximum depth of shallow water sediments (m)        : ',par_sed_Dmax_neritic
       print*,'Force reef occurrence?                              : ',ctrl_sed_neritic_reef_force
       print*,'Minimum (basic) number of sedcore layers            : ',par_n_sedcore_tot_min
       print*,'Number of dimensioned sedcore layers ka-1 runtime   : ',par_n_sedcore_tot_perky
       print*,'# sedimentary stack sub-layers                      : ',n_sed_tot
       print*,'# initial sedimentary stack sub-layers filled       : ',n_sed_tot_init
       print*,'# sedimentary stack sub-layers to drop off bottom   : ',n_sed_tot_drop
       ! --- DETRITAL CONFIGURATION ---------------------------------------------------------------------------------------------- !
       print*,'Flux of refractory material (g cm-2 kyr-1)          : ',par_sed_fdet
       print*,'No pelagic (dust) detrital contribution?            : ',ctrl_sed_det_NOdust
       ! --- DIAGENESIS SCHEME: SELECTION ---------------------------------------------------------------------------------------- !
       print*,'--- DIAGENESIS SCHEME: SELECTION -------------------'
       print*,'CaCO3 diagenesis scheme                             : ',par_sed_diagen_CaCO3opt
       print*,'opal diagenesis scheme                              : ',par_sed_diagen_opalopt
       print*,'Corg diagenesis scheme                              : ',par_sed_diagen_Corgopt
       print*,'POM_FeOOH diagenesis option                         : ',par_sed_diagen_POM_FeOOH_opt
       ! --- DIAGENESIS SCHEME: CONTROL ------------------------------------------------------------------------------------------ !
       print*,'--- DIAGENESIS SCHEME: CONTROL ---------------------'
       print*,'Bioturbate sediment stack?                          : ',ctrl_sed_bioturb
       print*,'Use Archer et al. [2002] bioturbation scheme?       : ',ctrl_sed_bioturb_Archer
       print*,'maximum layer depth for bioturbation                : ',par_n_sed_mix
       print*,'Max surface bioturbation mixing rate (cm2 yr-1)     : ',par_sed_mix_k_sur_max
       print*,'Min surface bioturbation mixing rate (cm2 yr-1)     : ',par_sed_mix_k_sur_min
       print*,'Prevent CaCO3 erosion (Fdis > Fsed)?                : ',ctrl_sed_noerosion
       print*,'CaCO3 interface dissolution?                        : ',ctrl_sed_interface
       print*,'CaCO3 red tracer tag fraction                       : ',par_sed_CaCO3_fred
       print*,'CaCO3 blue tracer tag fraction                      : ',par_sed_CaCO3_fblue
       print*,'Tag restart CaCO3?                                  : ',ctrl_sed_dyerestart
       print*,'Taged restart CaCO3 depth in layers (n)             : ',par_sed_dyerestart_n
       print*,'modification of calcite k (>1.0 == reduced sol)     : ',par_sed_diagen_kcalmod
       ! --- DIAGENESIS SCHEME: ORGANIC MATTER ----------------------------------------------------------------------------------- !
       print*,'--- DIAGENESIS SCHEME: ORGANIC MATTER --------------'
       print*,'Prevent frac2 from being remineralzied?             : ',ctrl_sed_diagen_preserve_frac2
       print*,'Fractional POC burial scaling (Dunne scheme)        : ',par_sed_diagen_fracCpres_scale
       print*,'Apply Wallmann [2010] C:P remin parameterization?   : ',ctrl_sed_diagen_fracC2Ppres_wallmann2010
       print*,'C:P remin C/P offset                                : ',par_sed_diagen_fracC2Ppres_off
       print*,'C:P remin [O2] threshold (mol kg-1)                 : ',par_sed_diagen_fracC2Ppres_c0_O2 
       print*,'Return of PO4 to ocean in Dunne 2007 scheme?        : ',ctrl_sed_dunne2007_remin_POP
       print*,'Cap Fe2+ dissolution at POM_FeOOH rain flux?        : ',ctrl_sed_diagen_POM_FeOOH_cap
       print*,'Retain original (Redfield) remin transformation?    : ',ctrl_sed_conv_sed_ocn_old
       print*,'[O2] thresh for switching redox arrays (mol kg-1)   : ',par_sed_diagen_O2thresh
       print*,'[NO3] thresh for switching redox arrays (mol kg-1)  : ',par_sed_diagen_NO3thresh
       print*,'[SO4] thresh for switching redox arrays (mol kg-1)  : ',par_sed_diagen_SO4thresh
       print*,'Use BIOGEM redox-dependent remin transformation     : ',ctrl_sed_conv_sed_ocn_redox
       print*,'Use Bohlen 2012 denitrification remin?              : ',ctrl_sed_conv_sedocn_bohlen2012
       ! --- DIAGENESIS SCHEME: HUELSE 2017 -------------------------------------------------------------------------------------- !
       print*,'--- DIAGENESIS SCHEME: HUELSE 2017 -----------------'
       print*,'Corg rate constant parameterization scheme          : ',par_sed_huelse2017_kscheme
       print*,'Corg degradation rates redox dependent?             : ',par_sed_huelse2017_redox
       print*,'labile Corg degradation rate constant (1/yr)        : ',par_sed_huelse2017_k1
       print*,'refractory Corg degradation rate constant (1/yr)    : ',par_sed_huelse2017_k2
       print*,'refractory Corg deg. rate order compared to labile  : ',par_sed_huelse2017_k2_order
       print*,'anoxic refractory Corg deg. rate constant (1/yr)    : ',par_sed_huelse2017_k2_anoxic
       print*,'Include explicit P-cycle in OMEN-SED?               : ',par_sed_huelse2017_P_cycle
       print*,'Remove implicit Alk associated with buried sulf-OM? : ',par_sed_huelse2017_remove_impl_sulALK
       print*,'Simulate ocean Porg loss with buried sulf-OM?       : ',par_sed_huelse2017_sim_P_loss
       print*,'Simulate ocean Porg loss related to Corg burial?    : ',par_sed_huelse2017_sim_P_loss_pres_fracC
       print*,'Simulate increased P-regeneration under anoxia?     : ',par_sed_huelse2017_sim_P_regeneration
       print*,'Return of PO4 to ocean in HUELSE 2017 scheme?       : ',ctrl_sed_huelse2017_remin_POP
      ! --- DIAGENESIS SCHEME: ARCHER 1991 -------------------------------------------------------------------------------------- !
       print*,'--- DIAGENESIS SCHEME: ARCHER 1991 -----------------'
       print*,'dissolution rate constant, units of 1/s             : ',par_sed_archer1991_dissc
       print*,'dissolution rate constant scaling, (%)              : ',par_sed_archer1991_disscpct
       print*,'dissolution rate order                              : ',par_sed_archer1991_dissn
       print*,'organic degradation rate constant, 1/s              : ',par_sed_archer1991_rc
       print*,'loop limit in <o2org> subroutine                    : ',par_sed_archer1991_iterationmax
       print*,'Use old error-catching scheme?                      : ',ctrl_sed_diagen_error_Archer_OLD
       print*,'Replace Archer model calc with lookup estimate?     : ',ctrl_sed_diagen_error_Archer2lookup
       ! --- DIAGENESIS SCHEME: opal --------------------------------------------------------------------------------------------- !
       print*,'base opal KSi value (yr-1)                          : ',par_sed_opal_KSi0
       ! --- CaCO3 PRODUCTION ---------------------------------------------------------------------------------------------------- !
       print*,'--- CaCO3 PRODUCTION -------------------------------'
       print*,'CaCO3 precip scale factor (abiotic) (mol cm-2 yr-1) : ',par_sed_CaCO3precip_sf
       print*,'CaCO3 precip rate law lower (abiotic)               : ',par_sed_CaCO3precip_exp
       print*,'CaCO3 precip scale factor (corals) (mol cm-2 yr-1)  : ',par_sed_reef_CaCO3precip_sf
       print*,'CaCO3 precip rate law lower (corals)                : ',par_sed_reef_CaCO3precip_exp
       print*,'CaCO3 precipitation as calcite (o/w aragonite)?     : ',par_sed_reef_calcite
       print*,'Min threshold for abiotic CaCO3 precipitation       : ',par_sed_CaCO3_abioticohm_min
       print*,'Reef CaCO3 porosity (cm3 cm-3)                      : ',par_sed_poros_CaCO3_reef
       print*,'prescribed CaCO3 production rate (mol cm-2 yr-1)    : ',par_sed_CaCO3burial
       print*,'prescribed global CaCO3 production rate (mol yr-1)  : ',par_sed_CaCO3burialTOT
       print*,'prescribed SrCO3 recryst rate (mol cm-2 yr-1)       : ',par_sed_SrCO3recryst
       print*,'prescribed global SrCO3 recryst rate (mol yr-1)     : ',par_sed_SrCO3recrystTOT
       print*,'carbonate recrystalization r87Sr                    : ',par_r87Sr_SrCO3recryst
       print*,'carbonate recrystalization d88Sr                    : ',par_d88Sr_SrCO3recryst
       ! --- TRACE METALS -------------------------------------------------------------------------------------------------------- !
       print*,'--- TRACE METALS -----------------------------------'
       print*,'Default CaCO3 Ca:Li ratio                           : ',par_bio_red_CaCO3_LiCO3
       print*,'Partition coefficient (alpha)                       : ',par_bio_red_CaCO3_LiCO3_alpha
       print*,'Partition coefficient (alpha)                       : ',par_bio_red_CaCO3_SrCO3_alpha
       ! --- ISOTOPIC FRACTIONATION ---------------------------------------------------------------------------------------------- !
       print*,'--- ISOTOPIC FRACTIONATION -------------------------'
       print*,'set fixed d13C fractionation of Corg w.r.t. CaCO3?  : ',ctrl_sed_Corgburial_fixedD13C
       print*,'fractionation for intercellular C fixation          : ',par_d13C_DIC_Corg_ef
       print*,'Benthic foram 13C fractionation scheme ID string    : ',trim(opt_sed_foram_b_13C_delta)
       print*,'7/6Li fractionation between Li and LiCO3            : ',par_d7Li_LiCO3_epsilon
       print*,'44/40Ca fractionation between Ca and CaCO3 (corals) : ',par_d44Ca_CaCO3_epsilon
       print*,'abiotic 44/40Ca fractionation between Ca and cal    : ',par_d44Ca_abioticcal_epsilon0
       print*,'abiotic 44/40Ca fractionation between Ca and arg    : ',par_d44Ca_abioticarg_epsilon0
       print*,'T-dependence of cal abiotic 44/40Ca fractionation   : ',par_d44Ca_abioticcal_epsilondT
       print*,'T-dependence of arg abiotic 44/40Ca fractionation   : ',par_d44Ca_abioticarg_epsilondT
       print*,'88/86 fractionation between Sr and SrCO3            : ',par_d88Sr_SrCO3_epsilon
       ! --- HYDROTHERMAL, OCEAN CRUSTAL WEATHERING, & CLAY FORMATION ------------------------------------------------------------ !
       print*,'--- HYDROTHERMAL, WEATHERING, & CLAY FORMATION -----'
       print*,'hydrothermal Li flux (mol yr-1)                     : ',par_sed_hydroip_fLi
       print*,'hydrothermal Li flux d7Li (o/oo)                    : ',par_sed_hydroip_fLi_d7Li
       print*,'Li low-T alteration sink (mol yr-1) (Li/Ca norm)    : ',par_sed_lowTalt_fLi_alpha
       print*,'Li low-T alteration sink 7Li epsilon (o/oo)         : ',par_sed_lowTalt_7Li_epsilon
       print*,'fixed (non T-dependent) low-T 7Li epsilon?          : ',ctrl_sed_lowTalt_7Li_epsilon_fixed
       print*,'T-dependent D7Li sensitivity (o/oo K-1)             : ',par_sed_lowTalt_7Li_epsilon_DT
       print*,'Li clay formation sink (mol yr-1) (Li/Ca norm)      : ',par_sed_clay_fLi_alpha
       print*,'Li clay formation sink 7Li epsilon (o/oo)           : ',par_sed_clay_7Li_epsilon
       print*,'fixed (non T-dependent) MACC 7Li epsilon?           : ',ctrl_sed_clay_7Li_epsilon_fixed
       print*,'T-dependent D7Li sensitivity (o/oo K-1)             : ',par_sed_clay_7Li_epsilon_DT
       print*,'hydrothermal Ca flux (mol yr-1)                     : ',par_sed_hydroip_fCa
       print*,'hydrothermal Ca flux d44Ca (o/oo)                   : ',par_sed_hydroip_fCa_d44Ca
       print*,'hydrothermal Mg flux (mol yr-1)                     : ',par_sed_hydroip_fMg
       print*,'Mg -> Ca flux option                                : ',trim(opt_sed_hydroip_MgtoCa)
       print*,'reference Mg/Ca ratio                               : ',par_sed_hydroip_rMgCaREF
       print*,'reference Mg concentration (mol kg-1)               : ',par_sed_hydroip_concMgREF
       print*,'Ca low-T alteration sink (mol yr-1) (Ca/Mg norm)    : ',par_sed_lowTalt_fCa_alpha
       print*,'Ca low-T alteration sink 44Ca epsilon (o/oo)        : ',par_sed_lowTalt_44Ca_epsilon
       print*,'hydrothermal Sr flux (mol yr-1)                     : ',par_sed_hydroip_fSr
       print*,'hydrothermal Sr flux r87Sr (87/86)                  : ',par_sed_hydroip_fSr_r87Sr
       print*,'hydrothermal Sr flux d88Sr (o/oo)                   : ',par_sed_hydroip_fSr_d88Sr
       print*,'Sr low-T alteration sink (mol yr-1)                 : ',par_sed_lowTalt_fSr_alpha
       print*,'CO2 low-T alteration (weathering!) sink (mol yr-1)  : ',par_sed_lowTalt_fCO2
       print*,'hydrothermal CO2 outgassing (mol yr-1)              : ',par_sed_hydroip_fDIC
       print*,'d13C                                                : ',par_sed_hydroip_fDIC_d13C
       print*,'Make hydrothermal input 2D?                         : ',ctrl_sed_Fhydr2D
       print*,'2D hydrothermal input filename                      : ',trim(par_sed_Fhydr2D_name)
       print*,'Mean ocean floor reference temeprature (C)          : ',par_sed_T0C
       ! --- MISC CONTROLS ------------------------------------------------------------------------------------------------------- !
       print*,'--- MISC CONTROLS ----------------------------------'
       print*,'Impose alt detrital burial flux forcing?            : ',ctrl_sed_Fdet
       print*,'Impose alt CaCO3 burial flux forcing?               : ',ctrl_sed_Fcaco3
       print*,'Impose alt opal burial flux forcing?                : ',ctrl_sed_Fopal
       print*,'Impose alt Corg preservation (burial) flux?         : ',ctrl_sed_Pcorg
       print*,'Impose alt Porg preservation (burial) flux?         : ',ctrl_sed_Pporg
       print*,'Impose alt preservation (burial) rain ratio?        : ',ctrl_sed_Prr
       print*,'Set dissolution flux = rain flux for CaCO3 only?    : ',ctrl_force_sed_closedsystem_CaCO3
       print*,'Set dissolution flux = rain flux for opal only?     : ',ctrl_force_sed_closedsystem_opal
       print*,'Impose alt sedimentation rates to sedcores?         : ',ctrl_sed_Fdet_sedcore
       print*,'Assumed fraction of silicate vs. total weathering   : ',par_sed_diag_fracSiweath
       print*,'Assumed d13C of volcanic emissions                  : ',par_sed_diag_volcanicd13C
       print*,'Assumed implicit P:ALK of weathering                : ',par_sed_diag_P2ALK
       print*,'Assumed implicit C:O2 of kerogen weathering         : ',par_sed_diag_C2O2
       print*,'Update kerogen O2 weathering consumption ratio?     : ',ctrl_sed_diag_balanceO2
       ! --- I/O: DIRECTORY DEFINITIONS ------------------------------------------------------------------------------------------ !
       print*,'--- I/O: DIRECTORY DEFINITIONS ---------------------'
       print*,'(Paleo config) input dir. name                      : ',trim(par_pindir_name)
       print*,'Input dir. name                                     : ',trim(par_indir_name)
       print*,'Output dir. name                                    : ',trim(par_outdir_name)
       print*,'Input restart dir. name                             : ',trim(par_inrstdir_name)
       print*,'Output restart dir. name                            : ',trim(par_outrstdir_name)
       print*,'Filename for restart input                          : ',trim(par_infile_name)
       print*,'Filename for restart output                         : ',trim(par_outfile_name)
       print*,'Sediment water depth grid name                      : ',trim(par_sed_topo_D_name)
       print*,'Shallow water sediment (coral reef) mask name       : ',trim(par_sed_reef_mask_name)
       print*,'Sediment core save mask name                        : ',trim(par_sedcore_save_mask_name)
       print*,'Sediment core save list name                        : ',trim(par_sedcore_save_list_name)
       print*,'Biodiffusion profile name                           : ',trim(par_sed_mix_k_name)
       print*,'File containing output years for 0D data            : ',trim(par_output_years_file_0d)
       print*,'File containing output years for 2D data            : ',trim(par_output_years_file_2d)
       print*,'Save threshold of accumulated time (yr)             : ',par_sed_age_save_dt
       print*,'Alt detrital burial flux forcing filename           : ',trim(par_sed_Fdet_name)
       print*,'Alt CaCO3 burial flux forcing filename              : ',trim(par_sed_Fcaco3_name)
       print*,'Alt opal burial flux forcing filename               : ',trim(par_sed_Fopal_name)
       print*,'Alt Corg preservation (burial) flux filename        : ',trim(par_sed_Pcorg_name)
       print*,'Alt Porg preservation (burial) flux filename        : ',trim(par_sed_Pporg_name)
       print*,'Alt preservation (burial)rain ratio filename        : ',trim(par_sed_Prr_name)
       ! --- I/O: MISC ----------------------------------------------------------------------------------------------------------- !
       print*,'--- I/O: MISC --------------------------------------'
       print*,'save timeseries output                              : ',ctrl_timeseries_output
       print*,'append data to output files on restart              : ',ctrl_append_data
       print*,'Save (sedcore) output in ascii format?              : ',ctrl_data_save_ascii
       print*,'Save sedcorenv output (ascii)?                      : ',ctrl_data_save_sedcorenv
       print*,'Report sediment data as a mass fraction?            : ',ctrl_data_save_wtfrac
       print*,'Debug level #1?                                     : ',ctrl_misc_debug1
       print*,'Debug level #2?                                     : ',ctrl_misc_debug2
       print*,'Debug level #3?                                     : ',ctrl_misc_debug3
       print*,'Debug level #4?                                     : ',ctrl_misc_debug4
       print*,'Report errors?                                      : ',ctrl_misc_report_err
       print*,'i sediment coordinate for debug reporting           : ',par_misc_debug_i
       print*,'j sediment coordinate for debug reporting           : ',par_misc_debug_j
       print*,'Report level #1 debug?                              : ',ctrl_debug_lvl1
       ! --- DATA SAVING: MISC --------------------------------------------------------------------------------------------------- !
       print*,'--- DATA SAVING: MISC ------------------------------'
       print*,'Restart in netCDF format?                           : ',ctrl_ncrst
       print*,'netCDF restart file name                            : ',trim(par_ncrst_name)
       print*,'2D netCDF sediments output file name                : ',trim(par_ncout2d_name)
       print*,'1D netCDF sedcore output file name                  : ',trim(par_ncsedcore_name)
       print*,'time interval for averaging final data over (yr)    : ',par_sed_save_av_dtyr
       print*,'Save diagenesis error details?                      : ',ctrl_sed_diagen_error_save
       ! #### INSERT CODE TO LOAD ADDITIONAL PARAMETERS ########################################################################## !
       !
       ! ######################################################################################################################### !
    end if
    ! set ash event flux (g cm-2 kyr-1) (e.g. set at ca. x10 typical dust/detrital value)
    par_sed_ashevent_fash = 1.8
    ! revise diagenesis options
    if (ctrl_sed_Fcaco3) par_sed_diagen_CaCO3opt = 'ALL'
    if (ctrl_sed_Fopal) par_sed_diagen_opalopt   = 'ALL'
    if (ctrl_sed_conv_sedocn_bohlen2012) ctrl_sed_conv_sed_ocn_redox = .true.

  END SUBROUTINE sub_load_goin_sedgem
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  ! *** LOAD SEDGEM RESTART DATA ************************************************************************************************* !
  ! ****************************************************************************************************************************** !
  SUBROUTINE sub_data_load_rst(dum_sfxsumsed,dum_sfxocn)
    USE sedgem_lib
    use gem_netcdf
    USE genie_util, ONLY:check_unit,check_iostat
    ! -------------------------------------------------------- !
    ! DEFINE DUMMY ARGUMENTS
    ! -------------------------------------------------------- !
    real,dimension(n_sed,n_i,n_j),intent(inout)::dum_sfxsumsed
    real,DIMENSION(n_ocn,n_i,n_j),intent(inout)::dum_sfxocn    ! sediment dissolution flux interface array
    ! -------------------------------------------------------- !
    ! DEFINE LOCAL VARIABLES
    ! -------------------------------------------------------- !
    integer::i,j,k,l,io,is,iv                                           ! local counting variables
    integer::ios                                               !
    integer::loc_ncid                                          !
    CHARACTER(len=255)::loc_filename                           ! filename string
    integer::loc_n_l_sed                                       ! number of selected tracers in the re-start file
    integer,DIMENSION(n_sed)::loc_conv_iselected_is            ! number of selected sediment tracers in restart
    real,dimension(n_i,n_j)::loc_ij                            ! 
    integer::loc_ndims,loc_nvars
    integer,ALLOCATABLE,dimension(:)::loc_dimlen
    integer,ALLOCATABLE,dimension(:,:)::loc_varlen
    integer,ALLOCATABLE,dimension(:)::loc_vdims
    character(20),ALLOCATABLE,dimension(:)::loc_varname
    integer::loc_n_sed_tot                                     !
    real,ALLOCATABLE,dimension(:,:,:)::loc_ijk
    ! -------------------------------------------------------- !
    ! INITIALIZE LOCAL VARIABLES
    ! -------------------------------------------------------- !
    IF (ctrl_misc_debug3) print*, 'INITIALIZE LOCAL VARIABLES'
    ! -------------------------------------------------------- ! set filename
    IF (ctrl_misc_debug4) print*, 'set filename'
    IF (ctrl_ncrst) THEN
       loc_filename = TRIM(par_inrstdir_name)//par_ncrst_name
    else
       loc_filename = TRIM(par_inrstdir_name)//trim(par_infile_name)
    endif
    ! -------------------------------------------------------- ! check file status
    IF (ctrl_misc_debug4) print*, 'check file status'
    call check_unit(in,__LINE__,__FILE__)
    OPEN(unit=in,status='old',file=loc_filename,form='unformatted',action='read',IOSTAT=ios)
    close(unit=in)
    If (ios /= 0) then
       CALL sub_report_error( &
            & 'sedgem_data','sub_data_load_restart', &
            & 'You have requested a CONTINUING run, but restart file <'//trim(loc_filename)//'> does not exist', &
            & 'SKIPPING - using default initial values', &
            & (/const_real_null/),.false. &
            & )
    else
       ! -------------------------------------------------------- !
       ! LOAD RESTART
       ! -------------------------------------------------------- !
       IF (ctrl_misc_debug3) print*, 'LOAD RESTART'
       IF (ctrl_ncrst) THEN
          call sub_openfile(loc_filename,loc_ncid)
          ! -------------------------------------------------------- ! determine number of variables
          IF (ctrl_misc_debug4) print*, 'determine number of variables'
          call sub_inqdims (loc_filename,loc_ncid,loc_ndims,loc_nvars)
          ! -------------------------------------------------------- ! allocate arrays
          IF (ctrl_misc_debug4) print*, 'allocate arrays'
          ALLOCATE(loc_dimlen(loc_ndims),STAT=alloc_error)
          call check_iostat(alloc_error,__LINE__,__FILE__)
          ALLOCATE(loc_varlen(2,loc_nvars),STAT=alloc_error)
          call check_iostat(alloc_error,__LINE__,__FILE__)
          ALLOCATE(loc_vdims(loc_nvars),STAT=alloc_error)
          call check_iostat(alloc_error,__LINE__,__FILE__)
          ALLOCATE(loc_varname(loc_nvars),STAT=alloc_error)
          call check_iostat(alloc_error,__LINE__,__FILE__)
          ! -------------------------------------------------------- ! get variable names
          IF (ctrl_misc_debug4) print*, 'get variable names'
          call sub_inqvars(loc_ncid,loc_ndims,loc_nvars,loc_dimlen,loc_varname,loc_vdims,loc_varlen)
          ! -------------------------------------------------------- ! determine stack height and allocate local array
          IF (ctrl_misc_debug4) print*, 'determine stack height and allocate local array'
          loc_n_sed_tot = loc_dimlen(5)
          If (loc_n_sed_tot > n_sed_tot) then
             CALL sub_report_error( &
                  & 'sedgem_data','sub_data_load_restart', &
                  & 'You have compiled in a smaller sediment stack <n='//fun_conv_num_char_n(4,n_sed_tot)// &
                  & '> than the restart <n='//fun_conv_num_char_n(4,loc_n_sed_tot)//'>', &
                  & 'I am not programmed for such liberal activities ... ENDING ...', &
                  & (/const_real_null/),.true. &
                  & )
          elseif (loc_n_sed_tot < n_sed_tot) then
             CALL sub_report_error( &
                  & 'sedgem_data','sub_data_load_restart', &
                  & 'You have compiled in a larger sediment stack <n='//fun_conv_num_char_n(4,n_sed_tot)// &
                  & '> than the restart <n='//fun_conv_num_char_n(4,loc_n_sed_tot)//'> ' // &
                  & 'You may experience a core hiatus ...', &
                  & 'CONTINUING ...', &
                  & (/const_real_null/),.false. &
                  & )
          end If
          ALLOCATE(loc_ijk(n_i,n_j,loc_n_sed_tot),STAT=alloc_error)
          call check_iostat(alloc_error,__LINE__,__FILE__)
          ! -------------------------------------------------------- ! load and apply sediment tracers that are selected
          IF (ctrl_misc_debug4) print*, 'load and apply sediment tracers that are selected'
          ! NOTE: the k dimension is flipped in sub_getvarijk
          IF (ctrl_debug_init == 1) print*,' * Loading sediment stack restart tracers: '
          DO iv=1,loc_nvars
             DO l=1,n_l_sed
                is = conv_iselected_is(l)
                if ('sed_'//trim(string_sed(is)) == trim(loc_varname(iv))) then
                   IF (ctrl_debug_init == 1) print*,'   ',trim(loc_varname(iv))
                   loc_ijk(:,:,:) = 0.0
                   call sub_getvarijk(loc_ncid,'sed_'//trim(string_sed(is)),n_i,n_j,loc_n_sed_tot,loc_ijk(:,:,:))
                   DO i = 1,n_i
                      DO j = 1,n_j
                         if (sed_mask(i,j)) then
                            sed_top(is,i,j) = loc_ijk(i,j,loc_n_sed_tot)
                            do k = 1,loc_n_sed_tot-1
                               sed(is,i,j,(n_sed_tot-loc_n_sed_tot)+k) = loc_ijk(i,j,k)
                            end do
                         else
                            sed(is,i,j,:) = 0.0
                         end if
                      end DO
                   END DO
                endif
             end do
          end DO
          ! -------------------------------------------------------- ! load and apply sediment stack height
          IF (ctrl_misc_debug4) print*, 'load and apply sediment stack height'
          call sub_getvarij(loc_ncid,'phys_dh',n_i,n_j,loc_ij(:,:))
          DO i = 1,n_i
             DO j = 1,n_j
                if (sed_mask(i,j)) then
                   sed_top_h(i,j) = real(n_sed_tot - 2) + loc_ij(i,j)
                end if
             end DO
          end DO
          ! -------------------------------------------------------- ! load and apply dissolution flux tracers
          IF (ctrl_misc_debug4) print*, 'load and apply dissolution flux tracers'
          IF (ctrl_debug_init == 1) print*,' * Loading dissolution flux restart tracers: '
          DO iv=1,loc_nvars
             DO l=1,n_l_ocn
                io = conv_iselected_io(l)
                if ('fdis_'//trim(string_ocn(io)) == trim(loc_varname(iv))) then
                   IF (ctrl_debug_init == 1) print*,'   ',trim(loc_varname(iv))
                   loc_ij(:,:) = 0.0
                   call sub_getvarij(loc_ncid,'fdis_'//trim(string_ocn(io)),n_i,n_j,loc_ij(:,:))
                   DO i = 1,n_i
                      DO j = 1,n_j
                         dum_sfxocn(io,i,j) = loc_ij(i,j)
                      end DO
                   END DO
                endif
             end do
          end DO
          ! -------------------------------------------------------- ! deallocate arrays
          IF (ctrl_misc_debug4) print*, 'deallocate arrays'
          deALLOCATE(loc_dimlen,STAT=alloc_error)
          call check_iostat(alloc_error,__LINE__,__FILE__)
          deALLOCATE(loc_varlen,STAT=alloc_error)
          call check_iostat(alloc_error,__LINE__,__FILE__)
          deALLOCATE(loc_vdims,STAT=alloc_error)
          call check_iostat(alloc_error,__LINE__,__FILE__)
          deALLOCATE(loc_varname,STAT=alloc_error)
          call check_iostat(alloc_error,__LINE__,__FILE__)
          deALLOCATE(loc_ijk,STAT=alloc_error)
          call check_iostat(alloc_error,__LINE__,__FILE__)
          ! -------------------------------------------------------- ! close file
          IF (ctrl_misc_debug4) print*, 'close file'
          call sub_closefile(loc_ncid)
       else
          OPEN(unit=in,status='old',file=loc_filename,form='unformatted',action='read',IOSTAT=ios)
          read(unit=in,iostat=ios)                                        &
               & loc_n_l_sed,                                             &
               & (loc_conv_iselected_is(l),l=1,loc_n_l_sed),              &
               & (sed(loc_conv_iselected_is(l),:,:,:),l=1,loc_n_l_sed),   &
               & (sed_top(loc_conv_iselected_is(l),:,:),l=1,loc_n_l_sed), &
               & sed_top_h(:,:),                                          &
               & (dum_sfxsumsed(conv_iselected_is(l),:,:),l=1,loc_n_l_sed)
          close(unit=in,iostat=ios)
          call check_iostat(ios,__LINE__,__FILE__)
       endif
    end If
    ! -------------------------------------------------------- !
    ! END
    ! -------------------------------------------------------- !
    IF (ctrl_misc_debug3) print*, 'END'
  end SUBROUTINE sub_data_load_rst
  ! ****************************************************************************************************************************** !

  
  ! ****************************************************************************************************************************** !
  ! INITIALIZE SEDIMENT GRID
  ! NOTE: the grid as set up is specific to the GOLDSTEIN ocean model in an equal-area configuration
  !       so this subroutine needs ot be replaced or revised to make SEDGEM compatible with another ocean model
  ! NOTE: the lat-lon grid information is not critical, but sets the details of the grid associated with the saved data
  !       however, the depth set in this subroutine determinds the hydrostatic pressure on the sediments
  !       (and thus the stability of CaCO3 in the sediments)
  ! NOTE: the value of ips_mix_k0 is not set here (even if if perhaps should be) ...
  !       it is set in sedgem_box/sub_update_sed and serves to record the surface bioturbation rate applied in the Archer scheme
  SUBROUTINE sub_init_phys_sed()
    ! local variables
    INTEGER::i,j
    INTEGER::loc_len
    CHARACTER(len=255)::loc_filename
    real::loc_th0,loc_th1,loc_s0,loc_s1,loc_ds
    real,dimension(0:n_j)::loc_s,loc_sv
    real,DIMENSION(n_i,n_j)::loc_ij                  ! 
    ! set alt dir path string length
    loc_len = LEN_TRIM(par_pindir_name)
    ! zero the grid information and 'physics' array
    loc_ij(:,:)     = 0.0
    phys_sed(:,:,:) = 0.0
    ! initialize masks
    sed_mask(:,:)      = .FALSE.
    sed_mask_reef(:,:) = .FALSE.
    sed_mask_muds(:,:) = .FALSE.
    ! calculate local constants
    loc_th0 = -const_pi/2                            ! 
    loc_th1 = const_pi/2                             ! 
    loc_s0 = sin(loc_th0)                            ! 
    loc_s1 = sin(loc_th1)                            !
    loc_ds = (loc_s1-loc_s0)/real(n_j)               ! 
    DO j=0,n_j
       loc_sv(j) = loc_s0 + real(j)*loc_ds           ! 
       loc_s(j) = loc_sv(j) - 0.5*loc_ds             ! 
    end do
    ! initialize array values
    DO i=1,n_i
       DO j=1,n_j
          phys_sed(ips_lat,i,j)  = (180.0/const_pi)*ASIN(loc_s(j))
          phys_sed(ips_lon,i,j)  = (360.0/real(n_i))*(real(i)-0.5) + par_grid_lon_offset
          phys_sed(ips_dlat,i,j) = (180.0/const_pi)*(ASIN(loc_sv(j)) - ASIN(loc_sv(j-1)))
          phys_sed(ips_dlon,i,j) = (360.0/real(n_i))
          phys_sed(ips_latn,i,j) = (180.0/const_pi)*ASIN(loc_sv(j))
          phys_sed(ips_lone,i,j) = (360.0/n_i)*real(i) + par_grid_lon_offset
          phys_sed(ips_A,i,j)    = 2.0*const_pi*(const_rEarth**2)*(1.0/real(n_i))*(loc_sv(j) - loc_sv(j-1))
          phys_sed(ips_rA,i,j)   = 1.0/phys_sed(ips_A,i,j)
       END DO
    END DO
    ! load sediment bathymetry
    if (loc_len > 0) then
        loc_filename = TRIM(par_pindir_name)//TRIM(par_sed_topo_D_name)
    else
        loc_filename = TRIM(par_indir_name)//TRIM(par_sed_topo_D_name)
    endif
    CALL sub_load_data_ij(loc_filename,n_i,n_j,loc_ij(:,:))
    phys_sed(ips_D,:,:) = loc_ij(:,:)
    ! load reef mask
    if (par_sed_Dmax_neritic > -const_real_nullsmall) then
       if (loc_len > 0) then
            loc_filename = TRIM(par_pindir_name)//TRIM(par_sed_reef_mask_name)
       else
            loc_filename = TRIM(par_indir_name)//TRIM(par_sed_reef_mask_name)            
       endif
       CALL sub_load_data_ij(loc_filename,n_i,n_j,loc_ij(:,:))
    else
       loc_ij(:,:) = 0.0
    endif
    ! define sediment masks - used as an area mulitplying factor
    ! (both in logial and area mulitplying factor (real) representations)
    ! NOTE: subsquently, the masks are updated depending on whether there ia an overlying ocean cell or not.
    !       (hence, these masks are just 'potential' locations here at the outset)
    DO i=1,n_i
       DO j=1,n_j
          if (phys_sed(ips_D,i,j) < const_real_nullsmall) then
             ! land! => no sediments!!!
             phys_sed(ips_mask_sed,i,j) = 0.0
             sed_mask(i,j) = .FALSE.
             phys_sed(ips_mask_sed_reef,i,j) = 0.0
             sed_mask_reef(i,j) = .FALSE.
             phys_sed(ips_mask_sed_muds,i,j) = 0.0
             sed_mask_muds(i,j) = .FALSE.
          else
             ! not land(!), so set sediment mask TRUE
             phys_sed(ips_mask_sed,i,j) = 1.0
             sed_mask(i,j) = .TRUE.
             if (phys_sed(ips_D,i,j) < par_sed_Dmax_neritic) then
                ! water shallower than generic neritic depth limit => either reef or mud!
                ! NOTE: if ctrl_sed_neritic_reef_force is set, then shallow points are forced to be reef
                !       (ctrl_sed_neritic_reef_force is .false. by default)
                if ((loc_ij(i,j) > const_real_nullsmall) .OR. ctrl_sed_neritic_reef_force) then
                   ! mask specified as reef ... therefore reef!
                   phys_sed(ips_mask_sed_reef,i,j) = 1.0
                   sed_mask_reef(i,j) = .TRUE.
                   phys_sed(ips_mask_sed_muds,i,j) = 0.0
                   sed_mask_muds(i,j) = .FALSE.
                else
                   ! mask not specified as reef -- you got mud instead!
                   phys_sed(ips_mask_sed_reef,i,j) = 0.0
                   sed_mask_reef(i,j) = .FALSE.
                   phys_sed(ips_mask_sed_muds,i,j) = 1.0
                   sed_mask_muds(i,j) = .TRUE.
                end if
             elseif (ctrl_sed_neritic_reef_force) then
                ! force reef occurrence regardless of depth (assuming depth greater than prescribed neritic limit)
                if (loc_ij(i,j) > const_real_nullsmall) then
                   phys_sed(ips_mask_sed_reef,i,j) = 1.0
                   sed_mask_reef(i,j) = .TRUE.
                   phys_sed(ips_mask_sed_muds,i,j) = 0.0
                   sed_mask_muds(i,j) = .FALSE.
                endif
             else
                ! otherwise ... no reef or mud!
                phys_sed(ips_mask_sed_reef,i,j) = 0.0
                sed_mask_reef(i,j) = .FALSE.
                phys_sed(ips_mask_sed_muds,i,j) = 0.0
                sed_mask_muds(i,j) = .FALSE.
             end if
          end if
       END DO
    END DO
  END SUBROUTINE sub_init_phys_sed
  ! ****************************************************************************************************************************** !

  
  ! ****************************************************************************************************************************** !
  ! META-OPTION SETUP AND PARAMETER VALUE CONSISTENCY CHECK
  SUBROUTINE sub_check_par_sedgem()
    ! local variables
    LOGICAL::loc_flag
    ! initialize variables
    loc_flag = .FALSE.
    ! check that the i,j debug reporting indices specified in sedgem_config.par are within maxis and maxjs
    If (par_misc_debug_i > n_i .OR. par_misc_debug_i < 1) then
       loc_flag = .TRUE.
       par_misc_debug_i = 1
    end if
    If (par_misc_debug_j > n_j .OR. par_misc_debug_j < 1) then
       loc_flag = .TRUE.
       par_misc_debug_j = 1
    end if
    if (loc_flag) then
       CALL sub_report_error( &
            & 'sedgem_data','sub_check_par_sedgem', &
            & 'the i,j indices for spatially-explicit debugging '// &
            & 'must be within the sediment grid limit specification', &
            & 'SETTING OFFENDING PARAMETER VALUES TO 1; CONTINUING', &
            & (/const_real_null/),.false. &
            & )
       loc_flag = .FALSE.
    end If
    IF ((.NOT. sed_select(is_det)) .OR. (.NOT. sed_select(is_ash))) THEN
       CALL sub_report_error( &
            & 'sedgem_data','sub_check_par_sedgem','Both det and ash tracers must be selected '// &
            & '(FILE: gem_config_sed.par)', &
            & 'STOPPING', &
            & (/const_real_null/),.true. &
            & )
    ENDIF
  end SUBROUTINE sub_check_par_sedgem
  ! ****************************************************************************************************************************** !

  
  ! ****************************************************************************************************************************** !
  ! INITIALIZE SEDIMENT PARAMETERS
  SUBROUTINE sub_init_sed()
    ! local variables
    INTEGER::l,is,n                      ! grid and tracer index counters
    integer::loc_len
    CHARACTER(len=255)::loc_filename
    real,DIMENSION(n_i,n_j)::loc_ij      !                ! 
    ! set alt dir path string length
    loc_len = LEN_TRIM(par_pindir_name)
    ! set default array values
    conv_sed_mol_cm3(:)      = 1.0       ! 
    conv_sed_cm3_mol(:)      = 1.0       ! 
    conv_sed_cm3_g(:)        = 1.0       ! 
    conv_sed_g_cm3(:)        = 1.0       ! 
    conv_sed_mask(:)         = 0.0       ! 
    ! zero flux arrays
    sed_fsed(:,:,:) = 0.0                ! 
    sed_fdis(:,:,:) = 0.0                ! 
    ! set up conversion of mol -> cm3 and cm3 -> g (and reciprocals)
    DO l=1,n_l_sed
       is = conv_iselected_is(l)
       ! criterion for particulate organic matter (POM), elemental components, and particle-reactive scavenged elements
       if ( &
            (is == is_POC) &
            & .OR. &
            & (sed_type(is) == par_sed_type_POM) &
            & ) then
          conv_sed_mol_cm3(is) = conv_POC_mol_cm3
          conv_sed_cm3_g(is)   = conv_POC_cm3_g
          ! criterion for carbonate, elemental components, and particle-reactive scavenged elements
       elseif ( &
            & (is == is_CaCO3) &
            & .OR. &
            & (sed_type(is) == par_sed_type_CaCO3) &
            & ) then
          conv_sed_mol_cm3(is) = conv_cal_mol_cm3
          conv_sed_cm3_g(is)   = conv_cal_cm3_g
          ! criterion for opal, elemental components, and particle-reactive scavenged elements
       elseif ( &
            & (is == is_opal) &
            & .OR. &
            & (sed_type(is) == par_sed_type_opal) &
            & ) then
          conv_sed_mol_cm3(is) = conv_opal_mol_cm3
          conv_sed_cm3_g(is)   = conv_opal_cm3_g
          ! detrital and refractory material
       elseif ( &
            & (is == is_det) &
            & .OR. &
            & (sed_type(is) == par_sed_type_abio) &
            & .OR. &
            & (sed_type(is) == par_sed_type_det) &
            & .OR. &
            & (sed_type(is) == par_sed_type_scavenged) &
            & ) then
          conv_sed_mol_cm3(is) = conv_det_mol_cm3
          conv_sed_cm3_g(is)   = conv_det_cm3_g
          ! 'dependent' components (isotopes and 'age')
       elseif ( &
            & (sed_type(is) > 10) &
            & .OR. &
            & (sed_type(is) == par_sed_type_age) &
            & ) then
          conv_sed_mol_cm3(is) = conv_sed_mol_cm3(sed_dep(is))
          conv_sed_cm3_g(is)   = conv_sed_cm3_g(sed_dep(is))
       end if
       ! reciprocal conversion
       if(conv_sed_mol_cm3(is) > const_real_nullsmall) conv_sed_cm3_mol(is) = 1.0/conv_sed_mol_cm3(is)
       if(conv_sed_cm3_g(is) > const_real_nullsmall)   conv_sed_g_cm3(is)   = 1.0/conv_sed_cm3_g(is)
    end DO
    ! set up the mask for defining which sedimentary components contribute to the actual volume of the sediments
    ! (and which are therefore 'virtual')
    ! => POC, CaCO3, opal, miscellaneous detrital material ('det'), ash, iron oxides (FeO)
    do is=1,n_sed
       SELECT CASE (sed_type(is))
       case (par_sed_type_bio,par_sed_type_abio)
          conv_sed_mask(is) = 1.0
       case default
          conv_sed_mask(is) = 0.0
       end select
    end do
    ! allocate size of look-up tables and load data -- CaCO3
    ! NOTE: check for problems allocating array space
    SELECT CASE (par_sed_diagen_CaCO3opt)
    CASE ('ridgwell2001lookup','ridgwell2001lookupvec','archer1991explicit')
       ALLOCATE(lookup_sed_dis_cal( &
            & lookup_i_D_min:lookup_i_D_max, &
            & lookup_i_dCO3_min:lookup_i_dCO3_max, &
            & lookup_i_frac_min:lookup_i_frac_max, &
            & lookup_i_fCorg_min:lookup_i_fCorg_max &
            & ),STAT=error)
       IF (error /= 0) THEN
          CALL sub_report_error( &
               & 'sedgem_data','sub_init_sed', &
               & 'Could not allocate space for CaCO3 diagenesis look-up table array', &
               & 'STOPPING', &
               & (/const_real_zero/),.TRUE. &
               & )
       ENDIF
       call sub_load_sed_dis_lookup_CaCO3()
    END select
    ! allocate and populate lookup table vectors
    ! NOTE: check for problems allocating array space
    SELECT CASE (par_sed_diagen_CaCO3opt)
    CASE ('ridgwell2001lookupvec')
       ALLOCATE(lookup_vec_D(lookup_i_D_min:lookup_i_D_max),STAT=error)
       ALLOCATE(lookup_vec_dco3(lookup_i_dCO3_min:lookup_i_dCO3_max),STAT=error)
       ALLOCATE(lookup_vec_frac(lookup_i_frac_min:lookup_i_frac_max),STAT=error)
       ALLOCATE(lookup_vec_fCorg(lookup_i_fCorg_min:lookup_i_fCorg_max),STAT=error)
       IF (error /= 0) THEN
          CALL sub_report_error( &
               & 'sedgem_data','sub_init_sed', &
               & 'Could not allocate space for look-up table dimension vectors', &
               & 'STOPPING', &
               & (/const_real_zero/),.TRUE. &
               & )
       ENDIF
       lookup_vec_D     = (lookup_D_max/lookup_i_D_max)*(/ (n,n=lookup_i_D_min,lookup_i_D_max) /)
       lookup_vec_dco3  = (lookup_dCO3_max/lookup_i_dCO3_max)*(/ (n,n=lookup_i_dCO3_min,lookup_i_dCO3_max) /)
       lookup_vec_frac  = (lookup_frac_max/lookup_i_frac_max)*(/ (n,n=lookup_i_frac_min,lookup_i_frac_max) /)
       lookup_vec_fCorg = (lookup_fCorg_max/lookup_i_fCorg_max)*(/ (n,n=lookup_i_fCorg_min,lookup_i_fCorg_max) /)
    end select
    ! allocate size of look-up tables and load data -- CaCO3
    ! NOTE: check for problems allocating array space
    SELECT CASE (par_sed_diagen_opalopt)
    CASE ('ridgwelletal2003lookup')
       ALLOCATE(lookup_sed_dis_opal( &
            & lookup_i_opalpc_min:lookup_i_opalpc_max, &
            & lookup_i_concSi_min:lookup_i_concSi_max, &
            & lookup_i_T_min:lookup_i_T_max, &
            & lookup_i_KSi0_min:lookup_i_KSi0_max, &
            & lookup_i_opaltorefrac_min:lookup_i_opaltorefrac_max &
            & ),STAT=error)
       IF (error /= 0) THEN
          CALL sub_report_error( &
               & 'sedgem_data','sub_init_sed', &
               & 'Could not allocate space for opal diagenesis look-up table array', &
               & 'STOPPING', &
               & (/const_real_zero/),.TRUE. &
               & )
       ENDIF
       call sub_load_sed_dis_lookup_opal()
    END select
    ! load and initialize neutral network
    if (par_sed_diagen_CaCO3opt == 'ridgwell2001nn') then
       call sub_init_neuralnetwork()
    end IF
    ! load alternative detrital flux field
    if (ctrl_sed_Fdet) then
       if (loc_len > 0) then
          loc_filename = TRIM(par_pindir_name)//TRIM(par_sed_Fdet_name)
       else
          loc_filename = TRIM(par_indir_name)//TRIM(par_sed_Fdet_name)
       endif
       CALL sub_load_data_ij(loc_filename,n_i,n_j,loc_ij(:,:))
    else
       loc_ij(:,:) = 0.0
    endif
    sed_Fsed_det = loc_ij
    ! load alternative CaCO3 flux field
    if (ctrl_sed_Fcaco3) then
       if (loc_len > 0) then
          loc_filename = TRIM(par_pindir_name)//TRIM(par_sed_Fcaco3_name)
       else
          loc_filename = TRIM(par_indir_name)//TRIM(par_sed_Fcaco3_name)
       endif
       CALL sub_load_data_ij(loc_filename,n_i,n_j,loc_ij(:,:))
    else
       loc_ij(:,:) = 0.0
    endif
    sed_Fsed_caco3 = loc_ij
    ! load alternative opal flux field
    if (ctrl_sed_Fopal) then
       if (loc_len > 0) then
          loc_filename = TRIM(par_pindir_name)//TRIM(par_sed_Fopal_name)
       else
          loc_filename = TRIM(par_indir_name)//TRIM(par_sed_Fopal_name)
       endif
       CALL sub_load_data_ij(loc_filename,n_i,n_j,loc_ij(:,:))
    else
       loc_ij(:,:) = 0.0
    endif
    sed_Fsed_opal = loc_ij
    ! load alternative hydrothermal input mask
    if (ctrl_sed_Fhydr2D) then
       if (loc_len > 0) then
          loc_filename = TRIM(par_pindir_name)//TRIM(par_sed_Fhydr2D_name)
       else
          loc_filename = TRIM(par_indir_name)//TRIM(par_sed_Fhydr2D_name)
       endif
       CALL sub_load_data_ij(loc_filename,n_i,n_j,loc_ij(:,:))
    else
       loc_ij(:,:) = 0.0
    endif
    sed_mask_hydr = loc_ij
    ! load alternative Corg preservation (burial) field (mol cm-2 yr-1)
    if (ctrl_sed_Pcorg) then
       if (loc_len > 0) then
          loc_filename = TRIM(par_pindir_name)//TRIM(par_sed_Pcorg_name)
       else
          loc_filename = TRIM(par_indir_name)//TRIM(par_sed_Pcorg_name)
       endif
       CALL sub_load_data_ij(loc_filename,n_i,n_j,loc_ij(:,:))
    else
       loc_ij(:,:) = 0.0
    endif
    sed_Psed_corg = loc_ij
    ! load alternative Porg preservation (burial) field (mol cm-2 yr-1)
    if (ctrl_sed_Pporg) then
       if (loc_len > 0) then
          loc_filename = TRIM(par_pindir_name)//TRIM(par_sed_Pporg_name)
       else
          loc_filename = TRIM(par_indir_name)//TRIM(par_sed_Pporg_name)
       endif
       CALL sub_load_data_ij(loc_filename,n_i,n_j,loc_ij(:,:))
    else
       loc_ij(:,:) = 0.0
    endif
    sed_Psed_porg = loc_ij
    ! load alternative Porg preservation (burial) rain ratio (C/P) field
    if (ctrl_sed_Prr) then
       if (loc_len > 0) then
          loc_filename = TRIM(par_pindir_name)//TRIM(par_sed_Prr_name)
       else
          loc_filename = TRIM(par_indir_name)//TRIM(par_sed_Prr_name)
       endif
       CALL sub_load_data_ij(loc_filename,n_i,n_j,loc_ij(:,:))
    else
       loc_ij(:,:) = 0.0
    endif
    sed_Psed_rr = loc_ij
    ! initialize diagnostics data array
    sed_diag(:,:,:) = 0.0
    ! initialize average sediment data arrays
    sed_av_fsed(:,:,:)     = 0.0
    sed_av_fdis(:,:,:)     = 0.0
    sed_av_coretop(:,:,:)  = 0.0
    sed_av_diag_err(:,:,:) = 0.0

  END SUBROUTINE sub_init_sed
  ! ****************************************************************************************************************************** !

  
  ! ****************************************************************************************************************************** !
  ! CONFIGURE AND INITIALIZE SEDIMENT LAYERS
  ! NOTE: configured to initialze sediments with ash in the surface layer and detrital material throughout the stack
  SUBROUTINE sub_init_sed_layers_default()
    ! local variables
    INTEGER::i,j,o
    real::loc_sed_poros
    real::loc_sed_poros_top
    ! zero arrays
    sed(:,:,:,:)          = 0.0
    sed_top(:,:,:)        = 0.0
    sed_top_h(:,:)        = 0.0
    sed_top_INTdth(:,:)   = 0.0
    ! grid loop
    DO i=1,n_i
       DO j=1,n_j
          IF (sed_mask(i,j)) THEN
             ! set sediment porosity
             if (sed_mask_reef(i,j)) then
                loc_sed_poros = par_sed_poros_CaCO3_reef
                loc_sed_poros_top = par_sed_poros_CaCO3_reef
             elseif (sed_mask_muds(i,j)) then
                loc_sed_poros = par_sed_poros_det
                loc_sed_poros_top = fun_calc_sed_poros_nsur(0.0,par_sed_top_th)
             else
                loc_sed_poros = par_sed_poros_det
                loc_sed_poros_top = fun_calc_sed_poros_nsur(0.0,par_sed_top_th)
             endif
             ! set default sediment stack values
             ! NOTE: sediment component volumes are in the units of 
             !       actual volume of solid matter per cm2 area of sub-layer
             ! NOTE: the surface layer is initialized with ash to provide a constant sedimentation rate chronology
             sed_top(:,i,j)      = 0.0
             sed_top(is_ash,i,j) = 0.1*(1.0 - loc_sed_poros_top)*par_sed_top_th
             sed_top(is_det,i,j) = 0.9*(1.0 - loc_sed_poros_top)*par_sed_top_th
             if (sed_select(is_det_age)) sed_top(is_det_age,i,j) = par_misc_t_runtime*sed_top(is_det,i,j)
             DO o = 1,n_sed_tot_init
                sed(:,i,j,o)      = 0.0
                sed(is_det,i,j,o) = (1.0 - loc_sed_poros)*1.0
                if (sed_select(is_det_age)) sed(is_det_age,i,j,o) = par_misc_t_runtime*sed(is_det,i,j,o)
             END DO
             DO o = (n_sed_tot_init + 1),n_sed_tot
                sed(:,i,j,o) = 0.0
             END DO
          END if
       end DO
    END DO
    ! set height of top layer of old sediment
    DO i=1,n_i
       DO j=1,n_j
          IF (sed_mask(i,j)) THEN
             sed_top_h(i,j) = REAL(n_sed_tot_init)
          end IF
       end DO
    END DO

  END SUBROUTINE sub_init_sed_layers_default
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  ! *** INITIALIZE SEDCORES ****************************************************************************************************** !
  ! ****************************************************************************************************************************** !
  ! NOTE: this mask sets the grid point locations where synthetic sediment 'cores' will saved, specified in the mask file by;
  !       1.0 = 'save here'
  !       0.0 = 'don't save here'
  !       (other values are not valid, or rather, could give rather unpredictable results ...)
  SUBROUTINE sub_data_sedcore_init()
    USE genie_util, ONLY: check_unit, check_iostat
    ! -------------------------------------------------------- !
    ! DEFINE LOCAL VARIABLES
    ! -------------------------------------------------------- !
    INTEGER::i,j,n
    integer::loc_len_pindir_name,loc_nmax
    CHARACTER(len=255)::loc_filename
    REAL,DIMENSION(n_i,n_j)::loc_ij
    integer,ALLOCATABLE,DIMENSION(:,:)::loc_vij                ! (i,j) vector
    REAL,ALLOCATABLE,DIMENSION(:)::loc_vd                      ! depth vector
    REAL,ALLOCATABLE,DIMENSION(:)::loc_vmar                    ! MAR vector (g cm-1 kyr-1)
    ! -------------------------------------------------------- !
    ! INITIALIZE LOCAL VARIABLES
    ! -------------------------------------------------------- !
    loc_ij(:,:) = 0.0
    sed_save_mask(:,:) = .FALSE.
    ! determine paleo dir path string length (otherwise zero for original directory structure)
    loc_len_pindir_name = LEN_TRIM(par_pindir_name)
    ! -------------------------------------------------------- !
    ! DETERMINE SEDCORES TO BE SAVED
    ! -------------------------------------------------------- !
    ! -------------------------------------------------------- ! load alt sediment core save data (if filename exists)
    !                                                            otherwise load sediment core save mask
    loc_ij(:,:) = 0.0
    if (LEN_TRIM(par_sedcore_save_list_name) > 0) then
       ! accommodate alt paleo data directory structure
       if (loc_len_pindir_name > 0) then
          loc_filename = TRIM(par_pindir_name)//TRIM(par_sedcore_save_list_name)
       else
          loc_filename = TRIM(par_indir_name)//TRIM(par_sedcore_save_list_name)
       endif
       ! determine number of data elements
       loc_nmax = fun_calc_data_n(loc_filename)
       ! allocate local vectors
       ALLOCATE(loc_vij(1:loc_nmax,2),STAT=alloc_error)
       call check_iostat(alloc_error,__LINE__,__FILE__)
       ALLOCATE(loc_vd(1:loc_nmax),STAT=alloc_error)
       call check_iostat(alloc_error,__LINE__,__FILE__)
       ALLOCATE(loc_vmar(1:loc_nmax),STAT=alloc_error)
       call check_iostat(alloc_error,__LINE__,__FILE__)
       ! read file
       ! NOTE: use extended read function if site-specific detrital fluxes (MAR values) are required
       ! NOTE: the same filename is used regardless
       if (ctrl_sed_Fdet_sedcore) then
          call sub_load_data_nptdmar(loc_filename,loc_nmax,loc_vij,loc_vd,loc_vmar)
       else
          call sub_load_data_nptd(loc_filename,loc_nmax,loc_vij,loc_vd)
       end if
       ! populate 2D sedcore mask file
       DO n=1,loc_nmax
          loc_ij(loc_vij(n,1),loc_vij(n,2)) = 1.0
       end do
       ! modify SEDGEM sediment ocean depth
       DO n=1,loc_nmax
          phys_sed(ips_D,loc_vij(n,1),loc_vij(n,2)) = loc_vd(n)
       end do
       ! modify SEDGEM detrital fluxes
       if (ctrl_sed_Fdet_sedcore) then
          DO n=1,loc_nmax
             sed_Fsed_det(loc_vij(n,1),loc_vij(n,2)) = loc_vmar(n)
          end do
       end if
       ! deallocate local vectors
       DEALLOCATE(loc_vij,STAT=dealloc_error)
       call check_iostat(dealloc_error,__LINE__,__FILE__)
       DEALLOCATE(loc_vd,STAT=dealloc_error)
       call check_iostat(dealloc_error,__LINE__,__FILE__)
       DEALLOCATE(loc_vmar,STAT=dealloc_error)
       call check_iostat(dealloc_error,__LINE__,__FILE__)
    else
       if (LEN_TRIM(par_sedcore_save_mask_name) > 0) then
          if (loc_len_pindir_name > 0) then
             loc_filename = TRIM(par_pindir_name)//TRIM(par_sedcore_save_mask_name)
          else
             loc_filename = TRIM(par_indir_name)//TRIM(par_sedcore_save_mask_name)
          endif
          CALL sub_load_data_ij(loc_filename,n_i,n_j,loc_ij(:,:))
       end if
    end if
    ! -------------------------------------------------------- ! set sediment save mask & count number of sedcores
    nv_sedcore = 0
    DO i=1,n_i
       DO j=1,n_j
          if ((loc_ij(i,j) < const_real_nullsmall) .OR. (.NOT. sed_mask(i,j))) then
             sed_save_mask(i,j) = .FALSE.
          else
             sed_save_mask(i,j) = .TRUE.
             nv_sedcore = nv_sedcore + 1
          end if
       end do
    end do
    ! -------------------------------------------------------- !
    ! ALLOCATED ARRAY SPACE
    ! -------------------------------------------------------- !
    ! NOTE: <sedcore_store> is used to accumulate the excess sed layers not retained (pop-ed off of the stack) in the full sed array
    !                       and as such needs to have a tracer dimension equal to the full sed tracer number
    !       <sedcore> is used to reconstruct the sediment cores
    !                 and only needs to be dimensioned as large as the number of saved tracers
    ALLOCATE(vsedcore_store(1:nv_sedcore),STAT=alloc_error)
    call check_iostat(alloc_error,__LINE__,__FILE__)
    do n=1,nv_sedcore
       allocate(vsedcore_store(n)%top(1:n_sedcore_tracer),STAT=alloc_error)
       call check_iostat(alloc_error,__LINE__,__FILE__)
       allocate(vsedcore_store(n)%lay(1:n_sedcore_tracer,1:n_sedcore_tot),STAT=alloc_error)
       call check_iostat(alloc_error,__LINE__,__FILE__)
    end do
    ! -------------------------------------------------------- !
    ! INITIALIZE SEDCORES
    ! -------------------------------------------------------- !
    if (nv_sedcore > 0) then
       DO n=1,nv_sedcore
          vsedcore_store(n)%ht = 0.0
          vsedcore_store(n)%top(:) = 0.0
          vsedcore_store(n)%lay(:,:) = 0.0
       end do
       ! set sedcore (i,j) locations
       ! NOTE: <n> used as counter to index [vsedcore_store]
       n = 0
       DO i=1,n_i
          DO j=1,n_j
             if (sed_save_mask(i,j)) then
                n = n + 1
                vsedcore_store(n)%i = i
                vsedcore_store(n)%j = j
                vsedcore_store(n)%save = .true.
             end if
          end do
       end do
    end if
    ! -------------------------------------------------------- !
    ! END
    ! -------------------------------------------------------- !
  end SUBROUTINE sub_data_sedcore_init
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  ! LOAD SEDIMENT DIAGENESIS LOOK-UP TABLES - CACO3
  SUBROUTINE sub_load_sed_dis_lookup_CaCO3()
    USE genie_util, ONLY: check_unit, check_iostat
    ! local variables
    INTEGER::a,b,d,e
    CHARACTER(len=255)::loc_filename
    integer::ios ! for file checks
    ! *** read in calcite dissolution look-up data ***
    loc_filename = TRIM(par_indir_name)//'lookup_calcite_4.dat'
    call check_unit(in,__LINE__,__FILE__)
    OPEN(unit=in,file=loc_filename,action='read',iostat=ios)
    call check_iostat(ios,__LINE__,__FILE__)
    ! read in data
    DO a = lookup_i_D_min,lookup_i_D_max,1
       DO b = lookup_i_dCO3_min,lookup_i_dCO3_max,1
          DO d = lookup_i_frac_min,lookup_i_frac_max,1
             DO e = lookup_i_fCorg_min,lookup_i_fCorg_max,1
                READ(unit=in,FMT='(F7.3)',iostat=ios) lookup_sed_dis_cal(a,b,d,e)
                call check_iostat(ios,__LINE__,__FILE__)
             END DO
          END DO
       END DO
    END DO
    ! close file pipe
    CLOSE(unit=in,iostat=ios)
    call check_iostat(ios,__LINE__,__FILE__)
    ! change units from (umol cm-2 yr-1) to (mol cm-2 yr-1)
    lookup_sed_dis_cal(:,:,:,:) = conv_umol_mol*lookup_sed_dis_cal(:,:,:,:)
  END SUBROUTINE sub_load_sed_dis_lookup_CaCO3
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  ! LOAD SEDIMENT DIAGENESIS LOOK-UP TABLES - OPAL
  SUBROUTINE sub_load_sed_dis_lookup_opal()
    USE genie_util, ONLY: check_unit, check_iostat
    ! local variables
    INTEGER::a,b,c,d,e
    CHARACTER(len=255)::loc_filename
    integer::ios  ! for file checks
    ! *** read in opal dissolution look-up data ***
    loc_filename = TRIM(par_indir_name)//'lookup_opal_5.dat'
    call check_unit(in,__LINE__,__FILE__)
    OPEN(unit=in,file=loc_filename,action='read',iostat=ios)
    call check_iostat(ios,__LINE__,__FILE__)
    ! read in data
    DO a = lookup_i_opalpc_min,lookup_i_opalpc_max,1
       DO b = lookup_i_concSi_min,lookup_i_concSi_max,1
          DO c = lookup_i_T_min,lookup_i_T_max,1
             DO d = lookup_i_KSi0_min,lookup_i_KSi0_max,1
                DO e = lookup_i_opaltorefrac_min,lookup_i_opaltorefrac_max,1
                   READ(unit=in,FMT='(F7.3)',iostat=ios) lookup_sed_dis_opal(a,b,c,d,e)
                   call check_iostat(ios,__LINE__,__FILE__)
                END DO
             END DO
          END DO
       END DO
    END DO
    ! close file pipe
    CLOSE(unit=in,iostat=ios)
    call check_iostat(ios,__LINE__,__FILE__)
    ! change units from (umol cm-2 yr-1) to (mol cm-2 yr-1)
    lookup_sed_dis_opal(:,:,:,:,:) = conv_umol_mol*lookup_sed_dis_opal(:,:,:,:,:)
  END SUBROUTINE sub_load_sed_dis_lookup_opal
  ! ****************************************************************************************************************************** !


  ! ********************************************************************************************************************************
  ! CONFIGURE AND INITIALIZE NEURAL NETWORK
  SUBROUTINE sub_init_neuralnetwork()
!!$    real(kind=8) :: bias2
!!$    real(kind=8),dimension(par_nn_neurons) :: bias1
!!$    real(kind=8),dimension(par_nn_target,par_nn_neurons) :: wts2
!!$    real(kind=8),dimension(par_nn_neurons,par_nn_input)  :: wts1
!!$    character(50)  :: loc_name
!!$    INTEGER:: loc_iou, loc_ndims, loc_nvars
!!$    INTEGER,dimension(10)  :: loc_dimlen
!!$    INTEGER,dimension(10)  :: loc_vdims
!!$    INTEGER,dimension(2,20):: loc_varlen
!!$    character(20),dimension(10) :: loc_varname
!!$    call sub_nn_allocate_network()
!!$    loc_name = TRIM(par_indir_name)//'nn_calcite_4.nc'
!!$    call sub_openfile (loc_name, loc_iou)
!!$    call sub_inqdims (loc_name, loc_iou, loc_ndims, loc_nvars)
!!$    call sub_inqvars (loc_iou, loc_ndims, loc_nvars, loc_dimlen, loc_varname, &
!!$         & loc_vdims, loc_varlen)
!!$    call sub_getvar1d (loc_iou, loc_varname(1),loc_dimlen(loc_varlen(1,1)),nn_mint)
!!$    call sub_getvar1d (loc_iou, loc_varname(2),loc_dimlen(loc_varlen(1,2)),nn_maxt)
!!$    call sub_getvar1d (loc_iou, loc_varname(3),loc_dimlen(loc_varlen(1,3)),nn_maxp)
!!$    call sub_getvar1d (loc_iou, loc_varname(4),loc_dimlen(loc_varlen(1,4)),nn_minp)
!!$    call sub_getvar2d (loc_iou, loc_varname(5),loc_dimlen(loc_varlen(1,5)), &
!!$                  & loc_dimlen(loc_varlen(2,5)),w1)
!!$    call sub_getvar2d (loc_iou, loc_varname(6),loc_dimlen(loc_varlen(1,6)), &
!!$                  & loc_dimlen(loc_varlen(2,6)),w2)
!!$    call sub_getvar1d (loc_iou, loc_varname(7),loc_dimlen(loc_varlen(1,7)),b1)
!!$    call sub_getvar1d (loc_iou, loc_varname(8),loc_dimlen(loc_varlen(1,8)),b2)
!!$    call sub_closefile(loc_iou)
  END SUBROUTINE sub_init_neuralnetwork
  ! ********************************************************************************************************************************


!!$  ! ****************************************************************************************************************************** !
!!$  ! INITIALIZE SEDIMENT DATA SAVING
!!$  ! NOTE: this mask sets the grid point locations where synthetic sediment 'cores' will saved, specified in the mask file by;
!!$  !       1.0 = 'save here'
!!$  !       0.0 = 'don't save here'
!!$  !       (other values are not valid, or rather, could give rather unpredictable results ...)
!!$  SUBROUTINE sub_init_sedgem_save_sed_data()
!!$    ! local variables
!!$    INTEGER::i,j
!!$    integer::loc_len
!!$    CHARACTER(len=255)::loc_filename
!!$    REAL,DIMENSION(n_i,n_j)::loc_ij             ! 
!!$    ! set alt dir path string length
!!$    loc_len = LEN_TRIM(par_pindir_name)
!!$    ! initialize variables
!!$    loc_ij(:,:) = 0.0
!!$    sed_save_mask(:,:) = .FALSE.
!!$    ! load sediment sediment save mask
!!$    if (loc_len > 0) then
!!$       loc_filename = TRIM(par_pindir_name)//TRIM(par_sedcore_save_mask_name)
!!$    else
!!$       loc_filename = TRIM(par_indir_name)//TRIM(par_sedcore_save_mask_name)           
!!$    endif
!!$    CALL sub_load_data_ij(loc_filename,n_i,n_j,loc_ij(:,:))
!!$    ! set sediment save mask
!!$    DO i=1,n_i
!!$       DO j=1,n_j
!!$          if (loc_ij(i,j) < const_real_nullsmall) then
!!$             sed_save_mask(i,j) = .FALSE.
!!$          else
!!$             sed_save_mask(i,j) = .TRUE.
!!$          end if
!!$       end do
!!$    end do
!!$  end SUBROUTINE sub_init_sedgem_save_sed_data
!!$  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  ! INIT SAVE SEDCORE
  SUBROUTINE sub_sedgem_init_sedcoresenv()
    USE genie_util, ONLY: check_unit, check_iostat
    ! local variables
    INTEGER::i,j
    CHARACTER(len=255)::loc_filename
    integer::ios ! for file checks
    ! create a file and save header information for specified core locations
    ! NOTE: <dum_sed_fdis> passed to sub_sedgem_save_sedcoreenv has units of units of (mol cm-2)
    !       and is converted to umol cm-2 yr-1 when written out
    DO i = 1,n_i
       DO j = 1,n_j
          if (sed_save_mask(i,j)) then
             loc_filename = TRIM(par_outdir_name)//'timeseries_x_sedcores_'// &
                  & fun_conv_num_char_n(2,i)//fun_conv_num_char_n(2,j)// &
                  & string_results_ext
             call check_unit(out,__LINE__,__FILE__)
             OPEN(unit=out,file=loc_filename,action='write',status='replace',iostat=ios)
             call check_iostat(ios,__LINE__,__FILE__)
             Write(unit=out,fmt=*) '% ========================================'
             write(unit=out,fmt='(A2,40A12)',iostat=ios)          &
                  & ' %',                                          &
                  & '        time',                               &
                  & '         age',                               &
                  & '           i',                               &
                  & '           j',                               &
                  & '         lon',                               &
                  & '         lat',                               &
                  & ' ocean depth',                               &
                  & ' k0 (mixing)',                               &
                  & '        temp',                               &
                  & '         sal',                               &
                  & '       [DIC]',                               &
                  & '       [ALK]',                               &
                  & '       [PO4]',                               &
                  & '       [NO3]',                               &
                  & '        [O2]',                               &
                  & '        [Ca]',                               &
                  & '       [SO4]',                               &
                  & '      [SiO4]',                               &
                  & '   d13C(DIC)',                               &
                  & '          pH',                               &
                  & '       [CO3]',                               &
                  & '    ohm(cal)',                               &
                  & '      d[CO3]',                               &
                  & '    POC conc',                               &
                  & '  CaCO3 conc',                               &
                  & ' d13C(CaCO3)',                               &
                  & '   opal conc',                               &
                  & '    det conc',                               &
                  & '    POC rain',                               &
                  & '   d13C(POC)',                               &
                  & '  CaCO3 rain',                               &
                  & ' d13C(CaCO3)',                               &
                  & '   opal rain',                               &
                  & '    det rain',                               &
                  & '     POC dis',                               &
                  & '   d13C(POC)',                               &
                  & '   CaCO3 dis',                               &
                  & ' d13C(CaCO3)',                               &
                  & '    opal dis',                               &
                  & '     det dis'
             call check_iostat(ios,__LINE__,__FILE__)
             write(unit=out,fmt='(A2,40A12)',iostat=ios)          &
                  & ' %',                                          &
                  & '         kyr',                               &
                  & '       kyrBP',                               &
                  & '           -',                               &
                  & '           -',                               &
                  & '       deg E',                               &
                  & '       deg N',                               &
                  & '           m',                               &
                  & '    cm2 yr-1',                               &
                  & '        degC',                               &
                  & '        o/oo',                               &
                  & '          uM',                               &
                  & '          uM',                               &
                  & '          uM',                               &
                  & '          uM',                               &
                  & '          uM',                               &
                  & '          mM',                               &
                  & '          mM',                               &
                  & '          uM',                               &
                  & '        o/oo',                               &
                  & '         SWS',                               &
                  & '          uM',                               &
                  & '           -',                               &
                  & '          uM',                               &
                  & '         wt%',                               &
                  & '         wt%',                               &
                  & '        o/oo',                               &
                  & '         wt%',                               &
                  & '         wt%',                               &
                  & ' molcm-2yr-1',                               &
                  & '        o/oo',                               &
                  & ' molcm-2yr-1',                               &
                  & '        o/oo',                               &
                  & ' molcm-2yr-1',                               &
                  & '  gcm-2kyr-1',                               &
                  & ' molcm-2yr-1',                               &
                  & '        o/oo',                               &
                  & ' molcm-2yr-1',                               &
                  & '        o/oo',                               &
                  & ' molcm-2yr-1',                               &
                  & '  gcm-2kyr-1'
             call check_iostat(ios,__LINE__,__FILE__)
             Write(unit=out,fmt=*) '% ========================================'
             CLOSE(unit=out,iostat=ios)
             call check_iostat(ios,__LINE__,__FILE__)
          end if
       end DO
    end DO
  end SUBROUTINE sub_sedgem_init_sedcoresenv
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  ! SAVE SEDCORE ENVIRONMENT
  SUBROUTINE sub_sedgem_save_sedcoreenv( &
       & dum_dtyr,                       &
       & dum_i,dum_j,                    &
       & dum_sed,                        &
       & dum_sed_fsed,                   &
       & dum_sed_fdis,                   &
       & dum_ocn,                        &
       & dum_sed_carb                    &
       & )
    USE genie_util, ONLY: check_unit, check_iostat
    ! dummy variables
    REAL,INTENT(in)::dum_dtyr                                    ! time-step (years)
    integer,INTENT(in)::dum_i,dum_j                              ! 
    REAL,INTENT(in),DIMENSION(n_sed)::dum_sed                    ! 
    REAL,INTENT(in),DIMENSION(n_sed)::dum_sed_fsed,dum_sed_fdis  ! 
    real,intent(in),DIMENSION(n_ocn)::dum_ocn                    ! ocean composition
    REAL,INTENT(in),DIMENSION(n_carb)::dum_sed_carb              ! 
    ! local variables
    CHARACTER(len=255)::loc_filename                             !
    real::loc_age                                                !
    real::loc_sed_tot_wt                                         !
    REAL,DIMENSION(n_sed)::loc_sed                               !
    real::loc_ocn_DIC_d13C,loc_sed_CaCO3_d13C                                      !
    real::loc_sed_fsed_POC_d13C,loc_sed_fsed_CaCO3_d13C          !
    real::loc_sed_fdis_POC_d13C,loc_sed_fdis_CaCO3_d13C          !
    integer::ios                                                 ! file checks
    ! calculate sediment comcposition (weight fraction)
    loc_sed_tot_wt = fun_calc_sed_mass(dum_sed(:))
    IF (loc_sed_tot_wt > const_real_nullsmall) THEN
       loc_sed(:) = conv_sed_cm3_g(:)*dum_sed(:)/loc_sed_tot_wt
    end IF
    ! calculate sedimentation age
    IF (dum_sed_fsed(is_CaCO3) > const_real_nullsmall) THEN
       loc_age = dum_sed_fsed(is_CaCO3_age)/dum_sed_fsed(is_CaCO3)
    ELSE
       loc_age = 0.0
    ENDIF
    ! calculate local d13C
    ! NOTE: pass -999.999 rather than NaN values to fun_calc_isotope_delta and hence prevent overflow when writing out ASCII
    loc_ocn_DIC_d13C        = &
         & fun_calc_isotope_delta(dum_ocn(io_DIC),dum_ocn(io_DIC_13C),const_standards(11),.FALSE.,const_nulliso)
    loc_sed_CaCO3_d13C      = &
         & fun_calc_isotope_delta(loc_sed(is_CaCO3),loc_sed(is_CaCO3_13C),const_standards(11),.FALSE.,const_nulliso)
    loc_sed_fsed_POC_d13C   = &
         & fun_calc_isotope_delta(dum_sed_fsed(is_POC),dum_sed_fsed(is_POC_13C),const_standards(11),.FALSE.,const_nulliso)
    loc_sed_fsed_CaCO3_d13C = &
         & fun_calc_isotope_delta(dum_sed_fsed(is_CaCO3),dum_sed_fsed(is_CaCO3_13C),const_standards(11),.FALSE.,const_nulliso)
    loc_sed_fdis_POC_d13C   = &
         & fun_calc_isotope_delta(dum_sed_fdis(is_POC),dum_sed_fdis(is_POC_13C),const_standards(11),.FALSE.,const_nulliso)
    loc_sed_fdis_CaCO3_d13C = &
         & fun_calc_isotope_delta(dum_sed_fdis(is_CaCO3),dum_sed_fdis(is_CaCO3_13C),const_standards(11),.FALSE.,const_nulliso)
    ! re-open file and write (append) data
    ! NOTE: for fdet, convert units from (mol cm-2 (per time-step)) to (g cm-2 kyr-1)
    loc_filename = TRIM(par_outdir_name)//'timeseries_x_sedcores_'// &
         & fun_conv_num_char_n(2,dum_i)//fun_conv_num_char_n(2,dum_j)// &
         & string_results_ext
    call check_unit(out,__LINE__,__FILE__)
    OPEN(unit=out,file=loc_filename,action='write',status='old',position='append',iostat=ios)
    call check_iostat(ios,__LINE__,__FILE__)
    write( &
         & unit=out,                                                                                             &
         & fmt='(2X,2f12.4,2i12,2f12.1,f12.1,f12.3,2f12.3,8f12.3,f12.3,f12.3,f12.3,f12.3,f12.3,5f12.3,12f12.3)', &
         & iostat=ios)                                                                                           &
         & conv_yr_kyr*(sed_time-0.5*dum_dtyr),                       &
         & conv_yr_kyr*loc_age,                                       &
         & dum_i,                                                     &
         & dum_j,                                                     &
         & phys_sed(ips_lon,dum_i,dum_j),                             &
         & phys_sed(ips_lat,dum_i,dum_j),                             &
         & phys_sed(ips_D,dum_i,dum_j),                               &
         & phys_sed(ips_mix_k0,dum_i,dum_j),                          &
         & dum_ocn(io_T) - const_zeroC,                               &
         & dum_ocn(io_S),                                             &
         & 1.0E+06*dum_ocn(io_DIC),                                   &
         & 1.0E+06*dum_ocn(io_ALK),                                   &
         & 1.0E+06*dum_ocn(io_PO4),                                   &
         & 1.0E+06*dum_ocn(io_NO3),                                   &
         & 1.0E+06*dum_ocn(io_O2),                                    &
         & 1.0E+03*dum_ocn(io_Ca),                                    &
         & 1.0E+03*dum_ocn(io_SO4),                                   &
         & 1.0E+06*dum_ocn(io_SiO2),                                  &
         & loc_ocn_DIC_d13C,                                          &
         & -log10(dum_sed_carb(ic_H)),                                &
         & 1.0E+06*dum_sed_carb(ic_conc_CO3),                         &
         & dum_sed_carb(ic_ohm_cal),                                  &
         & 1.0E+06*dum_sed_carb(ic_dCO3_cal),                         &
         & 100.0*loc_sed(is_POC),                                     &
         & 100.0*loc_sed(is_CaCO3),                                   &
         & loc_sed_CaCO3_d13C,                                        &
         & 100.0*loc_sed(is_opal),                                    &
         & 100.0*loc_sed(is_det),                                     &
         & 1.0E+06*dum_sed_fsed(is_POC)/dum_dtyr,                     &
         & loc_sed_fsed_POC_d13C,                                     &
         & 1.0E+06*dum_sed_fsed(is_CaCO3)/dum_dtyr,                   &
         & loc_sed_fsed_CaCO3_d13C,                                   &
         & 1.0E+06*dum_sed_fsed(is_opal)/dum_dtyr,                    &
         & dum_sed_fsed(is_det)/(conv_det_g_mol*(conv_yr_kyr*dum_dtyr)), &         
         & 1.0E+06*dum_sed_fdis(is_POC)/dum_dtyr,                     &
         & loc_sed_fdis_POC_d13C,                                     &
         & 1.0E+06*dum_sed_fdis(is_CaCO3)/dum_dtyr,                   &
         & loc_sed_fdis_CaCO3_d13C,                                   &
         & 1.0E+06*dum_sed_fdis(is_opal)/dum_dtyr,                    &
         & dum_sed_fdis(is_det)/(conv_det_g_mol*(conv_yr_kyr*dum_dtyr))
    call check_iostat(ios,__LINE__,__FILE__)
    CLOSE(unit=out,iostat=ios)
    call check_iostat(ios,__LINE__,__FILE__)
  end SUBROUTINE sub_sedgem_save_sedcoreenv
  ! ********************************************************************************************************************************


  ! ****************************************************************************************************************************** !
  ! LOAD IN SEDIMENT BIOTURBATIONAL MIXING PROFILE
  SUBROUTINE sub_load_sed_mix_k()
    USE genie_util, ONLY: check_unit, check_iostat
    ! local variables
    INTEGER::n
    INTEGER::loc_n_elements,loc_n_start
    CHARACTER(len=255)::loc_filename
    integer::ios  ! for file checks
    ! check file format
    loc_filename = TRIM(par_indir_name)//trim(par_sed_mix_k_name)
    CALL sub_check_fileformat(loc_filename,loc_n_elements,loc_n_start)
    ! set maximum number of sediment layers to bioturbate and therefore array size
    par_n_sed_mix = loc_n_elements - 1
    ALLOCATE(par_sed_mix_k(0:par_n_sed_mix),STAT=error)
    ! open file pipe
    call check_unit(in,__LINE__,__FILE__)
    OPEN(unit=in,file=loc_filename,action='read',iostat=ios)
    call check_iostat(ios,__LINE__,__FILE__)
    ! goto start-of-file tag
    DO n = 1,loc_n_start
       READ(unit=in,fmt='(1X)',iostat=ios)
       call check_iostat(ios,__LINE__,__FILE__)
    END DO
    ! read mixed-layer sediment mixing profile (top down)
    ! NOTE: the mixing rate as measured the depth of the top (incomplete stack layer) surface of the sediment stack
    !       has a value equal to the first bioturbation rate value read in, and corresponds to an index value of 'par_n_sed_mix'
    ! NOTE: mixing rate units from data file are (cm2 yr-1)
    ! NOTE: the bottom value in the mixing rate should be zero to properly terminate the profile
    DO n = par_n_sed_mix,0,-1
       READ(unit=in,FMT=*,iostat=ios) par_sed_mix_k(n)
       call check_iostat(ios,__LINE__,__FILE__)
    END DO
    CLOSE(unit=in,iostat=ios)
    call check_iostat(ios,__LINE__,__FILE__)
    ! convert the normalized mixing profile loaded in to a profile of biodiffusion rates
    ! accodring to the parameter par_sed_mix_kmax set in sedgem_config.par
    par_sed_mix_k(:) = par_sed_mix_k_sur_max*par_sed_mix_k(:)
    ! check that the maximum mixing rate in the profile does not exceed the maximum rate 
    ! that can be accomodated by the mixing algorithm (assuming 1 cm stack layer spacing)
    IF (MAXVAL(par_sed_mix_k(:)) > 0.5) THEN
       CALL sub_report_error( &
            & 'sedgem_data','sub_load_sed_mix_k', &
            & 'mixing time-step weighted sediment mixing rate is too large; '//&
            & 'maximum mixing rate in profile (cm2 yr-1) = ', &
            & 'STOPPING', &
            & (/MAXVAL(par_sed_mix_k(:))/),.true. &
            & )
    ENDIF
  END SUBROUTINE sub_load_sed_mix_k
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  ! SAVE SEDIMENT DIAGNOSTICS DATA
  SUBROUTINE sub_data_save_seddiag_GLOBAL(dum_dtyr,dum_sfcsumocn,dum_SLT)
    USE genie_util, ONLY: check_unit, check_iostat
    ! ---------------------------------------------------------- !
    ! DUMMY ARGUMENTS
    ! ---------------------------------------------------------- !
    real,INTENT(in)::dum_dtyr                                  ! 
    real,DIMENSION(n_ocn,n_i,n_j),intent(in)::dum_sfcsumocn    ! 
    real,INTENT(in)::dum_SLT
    ! ---------------------------------------------------------- !
    ! DEFINE LOCAL VARIABLES
    ! ---------------------------------------------------------- !
    INTEGER::i,j,l,is 
    integer::ios  ! for file checks
    CHARACTER(len=255)::loc_filename
    REAL,DIMENSION(n_sed,n_i,n_j)::loc_sed_coretop
    REAL,DIMENSION(n_sed,n_i,n_j)::loc_sed_preservation
    real::loc_tot1_sedgrid,loc_tot1a_sedgrid,loc_tot1b_sedgrid
    real::loc_tot2_sedgrid,loc_tot2a_sedgrid,loc_tot2b_sedgrid
    real::loc_pres_sedgrid
    real::loc_rain_sedgrid
    REAL,DIMENSION(n_sed,n_i,n_j)::loc_fsed                    ! 
    REAL,DIMENSION(n_sed,n_i,n_j)::loc_fdis                    ! 
    real::loc_mean_sedgrid                                     ! 
    real::loc_tot_mask_area                                    ! 
    real::loc_sed_d13C_mean                                    ! 
    real::loc_dt                                               ! local time-step for data saving (and averaging)
    REAL,DIMENSION(n_i,n_j)::loc_area                          ! local area (cm2)
    REAL,DIMENSION(n_i,n_j)::loc_mask                          ! local sediment (total) mask (copy)
    REAL,DIMENSION(n_i,n_j)::loc_mask_reef,loc_mask_muds       ! local reef, shallow sediment masks (copy)
    REAL,DIMENSION(n_i,n_j)::loc_mask_dsea                     ! local deep-sea sediment mask (derived variable)
    real::loc_tot,loc_frac,loc_standard,loc_sig                ! 
    REAL,DIMENSION(n_i,n_j)::loc_CaCO3_d13C,loc_POC_d13C       !
    real::loc_deep_FCaCO3,loc_mud_FCaCO3,loc_reef_FCaCO3
    real::loc_deep_FCaCO3_d13C,loc_mud_FCaCO3_d13C,loc_reef_FCaCO3_d13C
    real::loc_deep_FPOC,loc_mud_FPOC,loc_reef_FPOC
    real::loc_deep_FPOC_d13C,loc_mud_FPOC_d13C,loc_reef_FPOC_d13C
    real::loc_deep_FPOP,loc_mud_FPOP,loc_reef_FPOP
    real::loc_tot_FCaCO3,loc_tot_FPOC,loc_tot_FPOP
    real::loc_tot_FCaCO3_d13C,loc_tot_FPOC_d13C
    real::loc_gamma,loc_Foutgassing,loc_Fkerogen
    real::loc_FCaCO3_d13C,loc_tot_FO2
    ! -------------------------------------------------------- !
    ! INITIALIZE LOCAL VARIABLES
    ! -------------------------------------------------------- !
    ! averaging time-step
    loc_dt = 2.0*dum_dtyr
    ! area (units: cm2)
    loc_area(:,:) = conv_m2_cm2*phys_sed(ips_A,:,:)
    ! masks
    loc_mask(:,:)      = phys_sed(ips_mask_sed,:,:)
    loc_mask_reef(:,:) = phys_sed(ips_mask_sed_reef,:,:)
    loc_mask_muds(:,:) = phys_sed(ips_mask_sed_muds,:,:)
    loc_mask_dsea(:,:) = phys_sed(ips_mask_sed,:,:)*(1.0 - loc_mask_reef(:,:))*(1.0 - loc_mask_muds(:,:))
    ! calculate core-top sediment composition data
    loc_sed_coretop(:,:,:) = fun_sed_coretop()
    ! mean (last 2 time-step averaged) sediemnt and dissolution
    loc_fsed(:,:,:) = (sed_fsed(:,:,:) + sed_fsed_OLD(:,:,:))/loc_dt
    loc_fdis(:,:,:) = (sed_fdis(:,:,:) + sed_fdis_OLD(:,:,:))/loc_dt
    ! calculate local sediment preservation (normalized fraction)
    DO l=1,n_l_sed
       is = conv_iselected_is(l)
       DO i=1,n_i
          DO j=1,n_j
             IF (loc_fsed(is,i,j) > const_real_nullsmall) THEN
                loc_sed_preservation(is,i,j) = (loc_fsed(is,i,j) - loc_fdis(is,i,j))/loc_fsed(is,i,j)
             else
                loc_sed_preservation(is,i,j) = 0.0
             end if
          end do
       end do
    end do
    ! -------------------------------------------------------- ! calculate d13C
    DO i=1,n_i
       DO j=1,n_j
          if (loc_fsed(is_CaCO3,i,j) > const_real_nullsmall) then
             loc_tot  = loc_fsed(sed_dep(is_CaCO3_13C),i,j)
             loc_frac = loc_fsed(is_CaCO3_13C,i,j)
             loc_standard = const_standards(sed_type(is_CaCO3_13C))
             loc_CaCO3_d13C(i,j) = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard,.FALSE.,const_nulliso)
          else
             loc_CaCO3_d13C(i,j) = 0.0
          end if
          if (loc_fsed(is_POC,i,j) > const_real_nullsmall) then
             loc_tot  = loc_fsed(sed_dep(is_POC_13C),i,j)
             loc_frac = loc_fsed(is_POC_13C,i,j)
             loc_standard = const_standards(sed_type(is_POC_13C))
             loc_POC_d13C(i,j) = fun_calc_isotope_delta(loc_tot,loc_frac,loc_standard,.FALSE.,const_nulliso)
          else
             loc_POC_d13C(i,j) = 0.0
          end if
       end do
    end do
    ! -------------------------------------------------------- ! summary component values
    loc_deep_FCaCO3      = 0.0
    loc_mud_FCaCO3       = 0.0
    loc_reef_FCaCO3      = 0.0
    loc_deep_FCaCO3_d13C = 0.0
    loc_mud_FCaCO3_d13C  = 0.0
    loc_reef_FCaCO3_d13C = 0.0
    loc_deep_FPOC        = 0.0
    loc_mud_FPOC         = 0.0
    loc_reef_FPOC        = 0.0
    loc_deep_FPOC_d13C   = 0.0
    loc_mud_FPOC_d13C    = 0.0
    loc_reef_FPOC_d13C   = 0.0
    loc_deep_FPOP        = 0.0
    loc_mud_FPOP         = 0.0
    loc_reef_FPOP        = 0.0
    ! -------------------------------------------------------- !
    ! *** SAVE GLOBAL SUMMARY DATA ***
    ! -------------------------------------------------------- !
    ! -------------------------------------------------------- ! create output file
    ! set filename
    loc_filename = TRIM(par_outdir_name)//'INFO_sediment_summary_AT_END'//string_results_ext
    ! open file
    call check_unit(out,__LINE__,__FILE__)
    OPEN(out,file=TRIM(loc_filename),action='write',iostat=ios)
    call check_iostat(ios,__LINE__,__FILE__)
    ! -------------------------------------------------------- ! write file header
    Write(unit=out,fmt=*) '================================='
    Write(unit=out,fmt=*) '=== GLOBAL SEDIMENT DIAG DATA ==='
    Write(unit=out,fmt=*) '================================='
    ! -------------------------------------------------------- !
    ! DEEP-SEA SEDIMENT GRID DIAGNOSTICS
    ! -------------------------------------------------------- !
    Write(unit=out,fmt=*) ' '
    Write(unit=out,fmt=*) '--- DEEP-SEA SEDIMENT GRID ------'
    write(unit=out,fmt='(A6,f9.3,A19)',iostat=ios) &
         & ' --- > ',par_sed_Dmax_neritic,' m          -------'
    call check_iostat(ios,__LINE__,__FILE__)
    Write(unit=out,fmt=*) ' '
    ! MISC
    Write(unit=out,fmt=*) '---------------------------------'
    write(unit=out,fmt='(A28,I6)',iostat=ios) &
         & ' Total # deep-sea grid pts :',int(sum(loc_mask_dsea(:,:)))
    call check_iostat(ios,__LINE__,__FILE__)
    write(unit=out,fmt='(A28,e14.6,A3,A6,f6.2,A2)',iostat=ios) &
         & ' Total deep-sea area       :',sum(loc_mask_dsea(:,:)*phys_sed(ips_A,:,:)),' m2', &
         & '   =  ',100.0*sum(loc_mask_dsea(:,:)*phys_sed(ips_A,:,:))/sum(loc_mask(:,:)*phys_sed(ips_A,:,:)),' %'
    call check_iostat(ios,__LINE__,__FILE__)
    Write(unit=out,fmt=*) '---------------------------------'
    ! local variables 
    loc_tot_mask_area = sum(loc_mask_dsea(:,:)*loc_area(:,:))
    ! POC
    loc_tot1_sedgrid = sum(loc_mask_dsea(:,:)*loc_area(:,:)*loc_fsed(is_POC,:,:))
    loc_tot2_sedgrid = sum(loc_mask_dsea(:,:)*loc_area(:,:)*loc_fdis(is_POC,:,:))
    if (abs(loc_tot1_sedgrid) > const_real_nullsmall) then 
       loc_pres_sedgrid = 100.0*(loc_tot1_sedgrid - loc_tot2_sedgrid)/loc_tot1_sedgrid
    else
       loc_pres_sedgrid = 0.0
    end if
    if (loc_tot_mask_area > const_real_nullsmall) then 
       loc_mean_sedgrid = sum(loc_mask_dsea(:,:)*loc_sed_coretop(is_POC,:,:)*loc_area(:,:))/loc_tot_mask_area
    else
       loc_mean_sedgrid = 0.0
    end if
    write(unit=out,fmt='(A28,e14.6,A12,f7.3,A9)',iostat=ios) &
         & ' POC rain                  :',loc_tot1_sedgrid,' mol yr-1 = ',1.0E-12*conv_C_mol_kg*loc_tot1_sedgrid,' PgC yr-1'
    call check_iostat(ios,__LINE__,__FILE__)
    write(unit=out,fmt='(A28,e14.6,A12,f7.3,A9)',iostat=ios) &
         & ' POC diss                  :',loc_tot2_sedgrid,' mol yr-1 = ',1.0E-12*conv_C_mol_kg*loc_tot2_sedgrid,' PgC yr-1'
    call check_iostat(ios,__LINE__,__FILE__)
    write(unit=out,fmt='(A28,e14.6,A12,f7.3,A9,A6,f6.2,A2)',iostat=ios) &
         & ' Total POC pres            :',loc_tot1_sedgrid - loc_tot2_sedgrid,' mol yr-1 = ', &
         & 1.0E-12*conv_C_mol_kg*(loc_tot1_sedgrid - loc_tot2_sedgrid),' PgC yr-1','   =  ',loc_pres_sedgrid,' %'
    call check_iostat(ios,__LINE__,__FILE__)
    write(unit=out,fmt='(A28,f6.2,A2)',iostat=ios) &
         & ' Mean wt% POC              :',loc_mean_sedgrid,' %'
    call check_iostat(ios,__LINE__,__FILE__)
    Write(unit=out,fmt=*) '---------------------------------'
    ! SAVE !
    loc_deep_FPOC = loc_tot1_sedgrid - loc_tot2_sedgrid
    ! POC d13C (weighted by area and POC sedimentation rate)
    if (sum(loc_mask_dsea(:,:)*loc_area(:,:)*loc_fsed(is_POC,:,:)) > const_real_nullsmall) then
       loc_sed_d13C_mean = &
            & sum(loc_mask_dsea(:,:)*loc_area(:,:)*loc_fsed(is_POC,:,:)*loc_POC_d13C(:,:))/ &
            & sum(loc_mask_dsea(:,:)*loc_area(:,:)*loc_fsed(is_POC,:,:))
    else
       loc_sed_d13C_mean = -999.9
    end if
    write(unit=out,fmt='(A28,f6.2,A5)',iostat=ios) &
         & ' Mean weighted d13C POC    :',loc_sed_d13C_mean,'o/oo'
    Write(unit=out,fmt=*) '---------------------------------'
    ! SAVE !
    loc_deep_FPOC_d13C = loc_sed_d13C_mean
    ! POP
    loc_tot1_sedgrid = sum(loc_mask_dsea(:,:)*loc_area(:,:)*loc_fsed(is_POP,:,:))
    loc_tot2_sedgrid = sum(loc_mask_dsea(:,:)*loc_area(:,:)*loc_fdis(is_POP,:,:))
    if (abs(loc_tot1_sedgrid) > const_real_nullsmall) then 
       loc_pres_sedgrid = 100.0*(loc_tot1_sedgrid - loc_tot2_sedgrid)/loc_tot1_sedgrid
    else
       loc_pres_sedgrid = 0.0
    end if
    if (loc_pres_sedgrid > const_real_nullsmall) then 
       write(unit=out,fmt='(A28,e14.6,A9)',iostat=ios) &
            & ' POP rain                  :',loc_tot1_sedgrid,' mol yr-1'
       call check_iostat(ios,__LINE__,__FILE__)
       write(unit=out,fmt='(A28,e14.6,A9)',iostat=ios) &
            & ' POP diss                  :',loc_tot2_sedgrid,' mol yr-1'
       call check_iostat(ios,__LINE__,__FILE__)
       write(unit=out,fmt='(A28,e14.6,A12,A22,f6.2,A2)',iostat=ios) &
            & ' Total POP pres            :',loc_tot1_sedgrid - loc_tot2_sedgrid,' mol yr-1   ', &
            & '                 =  ',loc_pres_sedgrid,' %'
       call check_iostat(ios  ,__LINE__,__FILE__)
       Write(unit=out,fmt=*) '---------------------------------'
       ! SAVE ! 
       loc_deep_FPOP = loc_tot1_sedgrid - loc_tot2_sedgrid
       ! C/P
       if (loc_deep_FPOP > const_real_nullsmall) then
          write(unit=out,fmt='(A28,f6.2)',iostat=ios) &
               & ' C/P of Corg burial        :',loc_deep_FPOC/loc_deep_FPOP
          call check_iostat(ios,__LINE__,__FILE__)
          Write(unit=out,fmt=*) '---------------------------------'
       end if
    end if
    ! CaCO3
    loc_tot1_sedgrid = sum(loc_mask_dsea(:,:)*loc_area(:,:)*loc_fsed(is_CaCO3,:,:))
    loc_tot2_sedgrid = sum(loc_mask_dsea(:,:)*loc_area(:,:)*loc_fdis(is_CaCO3,:,:))
    if (abs(loc_tot1_sedgrid) > const_real_nullsmall) then 
       loc_pres_sedgrid = 100.0*(loc_tot1_sedgrid - loc_tot2_sedgrid)/loc_tot1_sedgrid
    else
       loc_pres_sedgrid = 0.0
    end if
    if (loc_tot_mask_area > const_real_nullsmall) then 
       loc_mean_sedgrid = sum(loc_mask_dsea(:,:)*loc_sed_coretop(is_CaCO3,:,:)*loc_area(:,:))/loc_tot_mask_area
    else
       loc_mean_sedgrid = 0.0
    end if
    write(unit=out,fmt='(A28,e14.6,A12,f7.3,A9)',iostat=ios) &
         & ' CaCO3 rain                :',loc_tot1_sedgrid,' mol yr-1 = ',1.0E-12*conv_CaCO3_mol_kgC*loc_tot1_sedgrid,' PgC yr-1'
    call check_iostat(ios,__LINE__,__FILE__)
    write(unit=out,fmt='(A28,e14.6,A12,f7.3,A9)',iostat=ios) &
         & ' CaCO3 diss                :',loc_tot2_sedgrid,' mol yr-1 = ',1.0E-12*conv_CaCO3_mol_kgC*loc_tot2_sedgrid,' PgC yr-1'
    call check_iostat(ios,__LINE__,__FILE__)
    write(unit=out,fmt='(A28,e14.6,A12,f7.3,A9,A6,f6.2,A2)',iostat=ios) &
         & ' Total CaCO3 pres          :',loc_tot1_sedgrid - loc_tot2_sedgrid,' mol yr-1 = ', &
         & 1.0E-12*conv_CaCO3_mol_kgC*(loc_tot1_sedgrid - loc_tot2_sedgrid),' PgC yr-1','   =  ',loc_pres_sedgrid,' %'
    call check_iostat(ios,__LINE__,__FILE__)
    write(unit=out,fmt='(A28,f6.2,A2)',iostat=ios) &
         & ' Mean wt% CaCO3            :',loc_mean_sedgrid,' %'
    call check_iostat(ios,__LINE__,__FILE__)
    Write(unit=out,fmt=*) '---------------------------------'
    ! SAVE !
    loc_deep_FCaCO3 = loc_tot1_sedgrid - loc_tot2_sedgrid 
    ! CaCO3 d13C (weighted by area and CaCO3 sedimentation rate)
    if (sum(loc_mask_dsea(:,:)*loc_area(:,:)*loc_fsed(is_CaCO3,:,:)) > const_real_nullsmall) then
       loc_sed_d13C_mean = &
            & sum(loc_mask_dsea(:,:)*loc_area(:,:)*loc_fsed(is_CaCO3,:,:)*loc_CaCO3_d13C(:,:))/ &
            & sum(loc_mask_dsea(:,:)*loc_area(:,:)*loc_fsed(is_CaCO3,:,:))
    else
       loc_sed_d13C_mean = -999.9
    end if
    write(unit=out,fmt='(A28,f6.2,A5)',iostat=ios) &
         & ' Mean weighted d13C CaCO3  :',loc_sed_d13C_mean,'o/oo'
    Write(unit=out,fmt=*) '---------------------------------'
    ! SAVE !
    loc_deep_FCaCO3_d13C = loc_sed_d13C_mean
    ! CaCO3:POC
    loc_tot1_sedgrid = SUM(loc_mask_dsea(:,:)*loc_fsed(is_POC,:,:))
    loc_tot2_sedgrid = SUM(loc_mask_dsea(:,:)*loc_fsed(is_CaCO3,:,:))
    if (abs(loc_tot1_sedgrid) > const_real_nullsmall) then 
       loc_rain_sedgrid = loc_tot2_sedgrid/loc_tot1_sedgrid
    else
       loc_rain_sedgrid = 0.0
    end if
    write(unit=out,fmt='(A28,f7.3)',iostat=ios) &
         & ' CaCO3/POC rain ratio      :',loc_rain_sedgrid
    call check_iostat(ios,__LINE__,__FILE__)
    Write(unit=out,fmt=*) '---------------------------------'
    ! opal
    if (sed_select(is_opal)) then
       loc_tot1_sedgrid = sum(loc_mask_dsea(:,:)*loc_area(:,:)*loc_fsed(is_opal,:,:))
       loc_tot2_sedgrid = sum(loc_mask_dsea(:,:)*loc_area(:,:)*loc_fdis(is_opal,:,:))
       if (abs(loc_tot1_sedgrid) > 0.0) then 
          loc_pres_sedgrid = 100.0*(loc_tot1_sedgrid - loc_tot2_sedgrid)/loc_tot1_sedgrid
       else
          loc_pres_sedgrid = 0.0
       end if
       if (loc_tot_mask_area > const_real_nullsmall) then 
          loc_mean_sedgrid = sum(loc_mask_dsea(:,:)*loc_sed_coretop(is_opal,:,:)*loc_area(:,:))/loc_tot_mask_area
       else
          loc_mean_sedgrid = 0.0
       end if
       Write(unit=out,fmt=*) '---------------------------------'
       write(unit=out,fmt='(A28,e14.6,A12,f7.3,A9)',iostat=ios) &
            & ' opal rain                 :',loc_tot1_sedgrid,' mol yr-1'
       call check_iostat(ios,__LINE__,__FILE__)
       write(unit=out,fmt='(A28,e14.6,A12,f7.3,A9)',iostat=ios) &
            & ' opal diss                 :',loc_tot2_sedgrid,' mol yr-1'
       call check_iostat(ios,__LINE__,__FILE__)
       write(unit=out,fmt='(A28,e14.6,A12,f6.2,A2)',iostat=ios) &
            & ' Total opal pres           :',loc_tot1_sedgrid - loc_tot2_sedgrid,' mol yr-1 = ',loc_pres_sedgrid,' %'
       call check_iostat(ios,__LINE__,__FILE__)
       write(unit=out,fmt='(A28,f6.2,A2)',iostat=ios) &
            & ' Mean wt% opal             :',loc_mean_sedgrid,' %'
       call check_iostat(ios,__LINE__,__FILE__)
       Write(unit=out,fmt=*) '---------------------------------'
    end if
    ! Li
    if (sed_select(is_LiCO3) .AND. sed_select(is_detLi)) then
       Write(unit=out,fmt=*) '---------------------------------'
       loc_tot1_sedgrid = sum(loc_mask_dsea(:,:)*loc_area(:,:)*loc_fsed(is_LiCO3,:,:))
       loc_tot2_sedgrid = sum(loc_mask_dsea(:,:)*loc_area(:,:)*loc_fdis(is_LiCO3,:,:))
       write(unit=out,fmt='(A28,e14.6,A12,f7.3,A9)',iostat=ios) &
            & ' Li CaCO3 sink             :',loc_tot1_sedgrid,' mol yr-1'
       call check_iostat(ios,__LINE__,__FILE__)
       write(unit=out,fmt='(A28,e14.6,A12,f7.3,A9)',iostat=ios) &
            & ' Li CaCO3 source           :',loc_tot2_sedgrid,' mol yr-1'
       call check_iostat(ios,__LINE__,__FILE__)
       loc_tot1_sedgrid = sum(loc_mask_dsea(:,:)*loc_area(:,:)*loc_fsed(is_detLi,:,:))
       loc_tot2_sedgrid = sum(loc_mask_dsea(:,:)*loc_area(:,:)*loc_fdis(is_detLi,:,:))
       write(unit=out,fmt='(A28,e14.6,A12,f7.3,A9)',iostat=ios) &
            & ' Li detrital sink          :',loc_tot1_sedgrid,' mol yr-1'
       call check_iostat(ios,__LINE__,__FILE__)
       write(unit=out,fmt='(A28,e14.6,A12,f7.3,A9)',iostat=ios) &
            & ' Li detrital source        :',loc_tot2_sedgrid,' mol yr-1'
       call check_iostat(ios,__LINE__,__FILE__)
       Write(unit=out,fmt=*) '---------------------------------'
    end if
    ! -------------------------------------------------------- !
    ! SHALLOW WATER SEDIMENT GRID DIAGNOSTICS
    ! -------------------------------------------------------- !
    Write(unit=out,fmt=*) ' '
    Write(unit=out,fmt=*) '--- SHALLOW SEDIMENT GRID -------'
    write(unit=out,fmt='(A6,f9.3,A19)',iostat=ios) &
         & ' --- < ',par_sed_Dmax_neritic,' m          -------'
    call check_iostat(ios,__LINE__,__FILE__)
    Write(unit=out,fmt=*) ' '
    ! MISC
    Write(unit=out,fmt=*) '---------------------------------'
    write(unit=out,fmt='(A28,I6)',iostat=ios) &
         & ' Total # grid pts          :',int(sum(loc_mask_muds(:,:)))
    call check_iostat(ios,__LINE__,__FILE__)
    write(unit=out,fmt='(A28,e14.6,A3,A6,f6.2,A2)',iostat=ios) &
         & ' Total area                :',sum(loc_mask_muds(:,:)*phys_sed(ips_A,:,:)),' m2', &
         & '   =  ',100.0*sum(loc_mask_muds(:,:)*phys_sed(ips_A,:,:))/sum(loc_mask(:,:)*phys_sed(ips_A,:,:)),' %'
    call check_iostat(ios,__LINE__,__FILE__)
    Write(unit=out,fmt=*) '---------------------------------'
    ! local variables 
    loc_tot_mask_area = sum(loc_mask_muds(:,:)*loc_area(:,:))
    ! POC
    loc_tot1_sedgrid = sum(loc_mask_muds(:,:)*loc_area(:,:)*loc_fsed(is_POC,:,:))
    loc_tot2_sedgrid = sum(loc_mask_muds(:,:)*loc_area(:,:)*loc_fdis(is_POC,:,:))
    if (abs(loc_tot1_sedgrid) > const_real_nullsmall) then 
       loc_pres_sedgrid = 100.0*(loc_tot1_sedgrid - loc_tot2_sedgrid)/loc_tot1_sedgrid
    else
       loc_pres_sedgrid = 0.0
    end if
    if (loc_tot_mask_area > const_real_nullsmall) then 
       loc_mean_sedgrid = sum(loc_mask_muds(:,:)*loc_sed_coretop(is_POC,:,:)*loc_area(:,:))/loc_tot_mask_area
    else
       loc_mean_sedgrid = 0.0
    end if
    write(unit=out,fmt='(A28,e14.6,A12,f7.3,A9)',iostat=ios) &
         & ' POC rain                  :',loc_tot1_sedgrid,' mol yr-1 = ',1.0E-12*conv_C_mol_kg*loc_tot1_sedgrid,' PgC yr-1'
    call check_iostat(ios,__LINE__,__FILE__)
    write(unit=out,fmt='(A28,e14.6,A12,f7.3,A9)',iostat=ios) &
         & ' POC diss                  :',loc_tot2_sedgrid,' mol yr-1 = ',1.0E-12*conv_C_mol_kg*loc_tot2_sedgrid,' PgC yr-1'
    call check_iostat(ios,__LINE__,__FILE__)
    write(unit=out,fmt='(A28,e14.6,A12,f7.3,A9,A6,f6.2,A2)',iostat=ios) &
         & ' Total POC pres            :',loc_tot1_sedgrid - loc_tot2_sedgrid,' mol yr-1 = ', &
         & 1.0E-12*conv_C_mol_kg*(loc_tot1_sedgrid - loc_tot2_sedgrid),' PgC yr-1','   =  ',loc_pres_sedgrid,' %'
    call check_iostat(ios,__LINE__,__FILE__)
    write(unit=out,fmt='(A28,f6.2,A2)',iostat=ios) &
         & ' Mean wt% POC              :',loc_mean_sedgrid,' %'
    call check_iostat(ios,__LINE__,__FILE__)
    Write(unit=out,fmt=*) '---------------------------------'
    ! SAVE !
    loc_mud_FPOC = loc_tot1_sedgrid - loc_tot2_sedgrid
    ! POC d13C (weighted by area and POC sedimentation rate)
    if (sum(loc_mask_muds(:,:)*loc_area(:,:)*loc_fsed(is_POC,:,:)) > const_real_nullsmall) then
       loc_sed_d13C_mean = &
            & sum(loc_mask_muds(:,:)*loc_area(:,:)*loc_fsed(is_POC,:,:)*loc_POC_d13C(:,:))/ &
            & sum(loc_mask_muds(:,:)*loc_area(:,:)*loc_fsed(is_POC,:,:))
    else
       loc_sed_d13C_mean = -999.9
    end if
    write(unit=out,fmt='(A28,f6.2,A5)',iostat=ios) &
         & ' Mean weighted d13C POC    :',loc_sed_d13C_mean,'o/oo'
    Write(unit=out,fmt=*) '---------------------------------' 
    ! SAVE !
    loc_mud_FPOC_d13C = loc_sed_d13C_mean
    ! POP
    loc_tot1_sedgrid = sum(loc_mask_muds(:,:)*loc_area(:,:)*loc_fsed(is_POP,:,:))
    loc_tot2_sedgrid = sum(loc_mask_muds(:,:)*loc_area(:,:)*loc_fdis(is_POP,:,:))
    if (abs(loc_tot1_sedgrid) > const_real_nullsmall) then 
       loc_pres_sedgrid = 100.0*(loc_tot1_sedgrid - loc_tot2_sedgrid)/loc_tot1_sedgrid
    else
       loc_pres_sedgrid = 0.0
    end if
    if (loc_pres_sedgrid > const_real_nullsmall) then 
       write(unit=out,fmt='(A28,e14.6,A9)',iostat=ios) &
            & ' POP rain                  :',loc_tot1_sedgrid,' mol yr-1'
       call check_iostat(ios,__LINE__,__FILE__)
       write(unit=out,fmt='(A28,e14.6,A9)',iostat=ios) &
            & ' POP diss                  :',loc_tot2_sedgrid,' mol yr-1'
       call check_iostat(ios,__LINE__,__FILE__)
       write(unit=out,fmt='(A28,e14.6,A12,A22,f6.2,A2)',iostat=ios) &
            & ' Total POP pres            :',loc_tot1_sedgrid - loc_tot2_sedgrid,' mol yr-1   ', &
            & '                   =  ',loc_pres_sedgrid,' %'
       call check_iostat(ios,__LINE__,__FILE__)
       Write(unit=out,fmt=*) '---------------------------------'  
       ! SAVE !
       loc_mud_FPOP = loc_tot1_sedgrid - loc_tot2_sedgrid  
       ! C/P
       if (loc_mud_FPOP > const_real_nullsmall) then
          write(unit=out,fmt='(A28,f6.2)',iostat=ios) &
               & ' C/P of Corg burial        :',loc_mud_FPOC/loc_mud_FPOP
          call check_iostat(ios,__LINE__,__FILE__)
          Write(unit=out,fmt=*) '---------------------------------'
       end if
    end if
    ! CaCO3
    loc_tot1_sedgrid = sum(loc_mask_muds(:,:)*loc_area(:,:)*loc_fsed(is_CaCO3,:,:))
    loc_tot2_sedgrid = sum(loc_mask_muds(:,:)*loc_area(:,:)*loc_fdis(is_CaCO3,:,:))
    if (abs(loc_tot1_sedgrid) > const_real_nullsmall) then 
       loc_pres_sedgrid = 100.0*(loc_tot1_sedgrid - loc_tot2_sedgrid)/loc_tot1_sedgrid
    else
       loc_pres_sedgrid = 0.0
    end if
    if (loc_tot_mask_area > const_real_nullsmall) then 
       loc_mean_sedgrid = sum(loc_mask_muds(:,:)*loc_sed_coretop(is_CaCO3,:,:)*loc_area(:,:))/loc_tot_mask_area
    else
       loc_mean_sedgrid = 0.0
    end if
    write(unit=out,fmt='(A28,e14.6,A12,f7.3,A9)',iostat=ios) &
         & ' CaCO3 rain                :',loc_tot1_sedgrid,' mol yr-1 = ',1.0E-12*conv_CaCO3_mol_kgC*loc_tot1_sedgrid,' GtC yr-1'
    call check_iostat(ios,__LINE__,__FILE__)
    write(unit=out,fmt='(A28,e14.6,A12,f7.3,A9)',iostat=ios) &
         & ' CaCO3 diss                :',loc_tot2_sedgrid,' mol yr-1 = ',1.0E-12*conv_CaCO3_mol_kgC*loc_tot2_sedgrid,' GtC yr-1'
    call check_iostat(ios,__LINE__,__FILE__)
    write(unit=out,fmt='(A28,e14.6,A12,f7.3,A9,A6,f6.2,A2)',iostat=ios) &
         & ' Total CaCO3 pres          :',loc_tot1_sedgrid - loc_tot2_sedgrid,' mol yr-1 = ', &
         & 1.0E-12*conv_CaCO3_mol_kgC*(loc_tot1_sedgrid - loc_tot2_sedgrid),' PgC yr-1','   =  ',loc_pres_sedgrid,' %'
    call check_iostat(ios,__LINE__,__FILE__)
    write(unit=out,fmt='(A28,f6.2,A2)',iostat=ios) &
         & ' Mean wt% CaCO3            :',loc_mean_sedgrid,' %'
    call check_iostat(ios,__LINE__,__FILE__)
    Write(unit=out,fmt=*) '---------------------------------'
    ! SAVE !
    loc_mud_FCaCO3 = loc_tot1_sedgrid - loc_tot2_sedgrid
    ! CaCO3:POC
    loc_tot1_sedgrid = SUM(loc_mask_muds(:,:)*(sed_fsed(is_POC,:,:) + sed_fsed_OLD(is_POC,:,:)))
    loc_tot2_sedgrid = SUM(loc_mask_muds(:,:)*(sed_fsed(is_CaCO3,:,:) + sed_fsed_OLD(is_CaCO3,:,:)))
    if (abs(loc_tot1_sedgrid) > const_real_nullsmall) then 
       loc_rain_sedgrid = loc_tot2_sedgrid/loc_tot1_sedgrid
    else
       loc_rain_sedgrid = 0.0
    end if
    write(unit=out,fmt='(A28,f7.3)',iostat=ios) &
         & ' CaCO3/POC rain ratio      :', &
         & loc_rain_sedgrid
    call check_iostat(ios,__LINE__,__FILE__)
    Write(unit=out,fmt=*) '---------------------------------'
    ! opal
    if (sed_select(is_opal)) then
       loc_tot1_sedgrid = sum(loc_mask_muds(:,:)*loc_area(:,:)*loc_fsed(is_opal,:,:))
       loc_tot2_sedgrid = sum(loc_mask_muds(:,:)*loc_area(:,:)*loc_fdis(is_opal,:,:))
       if (abs(loc_tot1_sedgrid) > 0.0) then 
          loc_pres_sedgrid = 100.0*(loc_tot1_sedgrid - loc_tot2_sedgrid)/loc_tot1_sedgrid
       else
          loc_pres_sedgrid = 0.0
       end if
       if (loc_tot_mask_area > const_real_nullsmall) then 
          loc_mean_sedgrid = sum(loc_mask_muds(:,:)*loc_sed_coretop(is_opal,:,:)*loc_area(:,:))/loc_tot_mask_area
       else
          loc_mean_sedgrid = 0.0
       end if
       write(unit=out,fmt='(A28,e14.6,A12,f7.3,A9)',iostat=ios) &
            & ' opal rain                 :',loc_tot1_sedgrid,' mol yr-1'
       call check_iostat(ios,__LINE__,__FILE__)
       write(unit=out,fmt='(A28,e14.6,A12,f7.3,A9)',iostat=ios) &
            & ' opal diss                 :',loc_tot2_sedgrid,' mol yr-1'
       call check_iostat(ios,__LINE__,__FILE__)
       write(unit=out,fmt='(A28,e14.6,A12,f6.2,A2)',iostat=ios) &
            & ' Total opal pres           :',loc_tot1_sedgrid - loc_tot2_sedgrid,' mol yr-1 = ',loc_pres_sedgrid,' %'
       call check_iostat(ios,__LINE__,__FILE__)
       write(unit=out,fmt='(A28,f6.2,A2)',iostat=ios) &
            & ' Mean wt% opal             :',loc_mean_sedgrid,' %'
       call check_iostat(ios,__LINE__,__FILE__)
       Write(unit=out,fmt=*) '---------------------------------'
    end if
    ! Li
    if (sed_select(is_LiCO3) .AND. sed_select(is_detLi)) then
       loc_tot1_sedgrid = sum(loc_mask_muds(:,:)*loc_area(:,:)*loc_fsed(is_LiCO3,:,:))
       loc_tot2_sedgrid = sum(loc_mask_muds(:,:)*loc_area(:,:)*loc_fdis(is_LiCO3,:,:))
       write(unit=out,fmt='(A28,e14.6,A12,f7.3,A9)',iostat=ios) &
            & ' Li CaCO3 sink             :',loc_tot1_sedgrid,' mol yr-1'
       call check_iostat(ios,__LINE__,__FILE__)
       write(unit=out,fmt='(A28,e14.6,A12,f7.3,A9)',iostat=ios) &
            & ' Li CaCO3 source           :',loc_tot2_sedgrid,' mol yr-1'
       call check_iostat(ios,__LINE__,__FILE__)
       loc_tot1_sedgrid = sum(loc_mask_muds(:,:)*loc_area(:,:)*loc_fsed(is_detLi,:,:))
       loc_tot2_sedgrid = sum(loc_mask_muds(:,:)*loc_area(:,:)*loc_fdis(is_detLi,:,:))
       write(unit=out,fmt='(A28,e14.6,A12,f7.3,A9)',iostat=ios) &
            & ' Li detrital sink          :',loc_tot1_sedgrid,' mol yr-1'
       call check_iostat(ios,__LINE__,__FILE__)
       write(unit=out,fmt='(A28,e14.6,A12,f7.3,A9)',iostat=ios) &
            & ' Li detrital source        :',loc_tot2_sedgrid,' mol yr-1'
       call check_iostat(ios,__LINE__,__FILE__)
       Write(unit=out,fmt=*) '---------------------------------'
    end if
    ! -------------------------------------------------------- !
    ! CARBONATE REEF GRID DIAGNOSTICS
    ! -------------------------------------------------------- !
    Write(unit=out,fmt=*) ' '
    Write(unit=out,fmt=*) '--- REEF SEDIMENT GRID ----------'
    Write(unit=out,fmt=*) ' '
    ! MISC
    Write(unit=out,fmt=*) '---------------------------------'
    write(unit=out,fmt='(A28,I6)',iostat=ios) &
         & ' Total # reef grid pts     :',int(sum(loc_mask_reef(:,:)))
    call check_iostat(ios,__LINE__,__FILE__)
    write(unit=out,fmt='(A28,e14.6,A3,A6,f6.2,A2)',iostat=ios) &
         & ' Total reef area           :',sum(loc_mask_reef(:,:)*phys_sed(ips_A,:,:)),' m2', &
         & '   =  ',100.0*sum(loc_mask_reef(:,:)*phys_sed(ips_A,:,:))/sum(loc_mask(:,:)*phys_sed(ips_A,:,:)),' %'
    call check_iostat(ios,__LINE__,__FILE__)
    Write(unit=out,fmt=*) '---------------------------------'
    ! local variables 
    loc_tot_mask_area = sum(loc_mask_reef(:,:)*loc_area(:,:))
    ! CaCO3
    loc_tot1_sedgrid = sum(loc_mask_reef(:,:)*loc_area(:,:)*loc_fsed(is_CaCO3,:,:))
    loc_tot2_sedgrid = sum(loc_mask_reef(:,:)*loc_area(:,:)*loc_fdis(is_CaCO3,:,:))
    if (abs(loc_tot1_sedgrid) > const_real_nullsmall) then 
       loc_pres_sedgrid = 100.0*(loc_tot1_sedgrid - loc_tot2_sedgrid)/loc_tot1_sedgrid
    else
       loc_pres_sedgrid = 0.0
    end if
    if (loc_tot_mask_area > const_real_nullsmall) then 
       loc_mean_sedgrid = sum(loc_mask_reef(:,:)*loc_sed_coretop(is_CaCO3,:,:)*loc_area(:,:))/loc_tot_mask_area
    else
       loc_mean_sedgrid = 0.0
    end if
    write(unit=out,fmt='(A28,e14.6,A12,f7.3,A9)',iostat=ios) &
         & ' CaCO3 production          :',loc_tot1_sedgrid,' mol yr-1 = ',1.0E-12*conv_CaCO3_mol_kgC*loc_tot1_sedgrid,' GtC yr-1'
    call check_iostat(ios,__LINE__,__FILE__)
    write(unit=out,fmt='(A28,f6.2,A2)',iostat=ios) &
         & ' Mean wt% CaCO3            :',loc_mean_sedgrid,' %'
    call check_iostat(ios,__LINE__,__FILE__)
    Write(unit=out,fmt=*) '---------------------------------'
    ! SAVE !
    loc_reef_FCaCO3 = loc_tot1_sedgrid
    ! d13C (weighted by area and CaCO3 sedimentation rate)
    if (sum(loc_mask_reef(:,:)*loc_area(:,:)*loc_fsed(is_CaCO3,:,:)) > const_real_nullsmall) then
       loc_sed_d13C_mean = &
            & sum(loc_mask_reef(:,:)*loc_area(:,:)*loc_fsed(is_CaCO3,:,:)*loc_CaCO3_d13C(:,:))/ &
            & sum(loc_mask_reef(:,:)*loc_area(:,:)*loc_fsed(is_CaCO3,:,:))
    else
       loc_sed_d13C_mean = 0.0
    end if
    write(unit=out,fmt='(A28,f6.2,A5)',iostat=ios) &
         & ' Mean weighted d13C CaCO3  :',loc_sed_d13C_mean,'o/oo'
    Write(unit=out,fmt=*) '---------------------------------'
    ! SAVE !
    loc_reef_FCaCO3_d13C = loc_sed_d13C_mean
    ! Li
    if (sed_select(is_LiCO3) .AND. sed_select(is_detLi)) then
       loc_tot1_sedgrid = sum(loc_mask_reef(:,:)*loc_area(:,:)*loc_fsed(is_LiCO3,:,:))
       loc_tot2_sedgrid = sum(loc_mask_reef(:,:)*loc_area(:,:)*loc_fdis(is_LiCO3,:,:))
       write(unit=out,fmt='(A28,e14.6,A12,f7.3,A9)',iostat=ios) &
            & ' Li CaCO3 sink             :',loc_tot1_sedgrid,' mol yr-1'
       call check_iostat(ios,__LINE__,__FILE__)
       write(unit=out,fmt='(A28,e14.6,A12,f7.3,A9)',iostat=ios) &
            & ' Li CaCO3 source           :',loc_tot2_sedgrid,' mol yr-1'
       call check_iostat(ios,__LINE__,__FILE__)
       loc_tot1_sedgrid = sum(loc_mask_reef(:,:)*loc_area(:,:)*loc_fsed(is_detLi,:,:))
       loc_tot2_sedgrid = sum(loc_mask_reef(:,:)*loc_area(:,:)*loc_fdis(is_detLi,:,:))
       write(unit=out,fmt='(A28,e14.6,A12,f7.3,A9)',iostat=ios) &
            & ' Li detrital sink          :',loc_tot1_sedgrid,' mol yr-1'
       call check_iostat(ios,__LINE__,__FILE__)
       write(unit=out,fmt='(A28,e14.6,A12,f7.3,A9)',iostat=ios) &
            & ' Li detrital source        :',loc_tot2_sedgrid,' mol yr-1'
       call check_iostat(ios,__LINE__,__FILE__)
       Write(unit=out,fmt=*) '---------------------------------'
    end if
    ! Sr
    IF (sed_select(is_SrCO3_87Sr) .AND. sed_select(is_SrCO3_88Sr)) THEN
       ! local variables
       loc_tot1_sedgrid  = sum(loc_mask_reef(:,:)*loc_area(:,:)*loc_fsed(is_SrCO3,:,:))
       loc_tot1a_sedgrid = sum(loc_mask_reef(:,:)*loc_area(:,:)*loc_fsed(is_SrCO3_87Sr,:,:))
       loc_tot1b_sedgrid = sum(loc_mask_reef(:,:)*loc_area(:,:)*loc_fsed(is_SrCO3_88Sr,:,:))
       loc_tot2_sedgrid  = sum(loc_mask_reef(:,:)*loc_area(:,:)*loc_fdis(is_SrCO3,:,:))
       loc_tot2a_sedgrid = sum(loc_mask_reef(:,:)*loc_area(:,:)*loc_fdis(is_SrCO3_87Sr,:,:))
       loc_tot2b_sedgrid = sum(loc_mask_reef(:,:)*loc_area(:,:)*loc_fdis(is_SrCO3_88Sr,:,:))
       ! bulk Sr
       write(unit=out,fmt='(A28,e14.6,A12,f7.3,A9)',iostat=ios) &
            & ' SrCO3 sink                :',loc_tot1_sedgrid,' mol yr-1'
       write(unit=out,fmt='(A28,e14.6,A12,f7.3,A9)',iostat=ios) &
            & ' Sr source                 :',loc_tot2_sedgrid,' mol yr-1'
       ! 87Sr
       loc_tot = loc_tot1_sedgrid-loc_tot1a_sedgrid-loc_tot1b_sedgrid
       if (loc_tot > const_real_nullsmall) then
          loc_sig = loc_tot1a_sedgrid/loc_tot
       else
          loc_sig = 0.0
       end if
       write(unit=out,fmt='(A28,f10.6)',iostat=ios) &
            & ' SrCO3 sink -- 87Sr        :',loc_sig
       loc_tot = loc_tot2_sedgrid-loc_tot2a_sedgrid-loc_tot2b_sedgrid
       if (loc_tot > const_real_nullsmall) then
          loc_sig = loc_tot2a_sedgrid/loc_tot
       else
          loc_sig = 0.0
       end if
       write(unit=out,fmt='(A28,f10.6)',iostat=ios) &
            & ' Sr source -- 87Sr         :',loc_sig
       ! 88Sr
       loc_tot = loc_tot1_sedgrid-loc_tot1a_sedgrid-loc_tot1b_sedgrid
       if (loc_tot > const_real_nullsmall) then
          loc_frac     = loc_tot1b_sedgrid
          loc_standard = const_standardsR(ocn_type(io_Sr_88Sr))
          loc_sig      = fun_calc_isotope_deltaR(loc_tot,loc_frac,loc_standard,const_real_null)
       else
          loc_sig = -999.9
       end if
       write(unit=out,fmt='(A28,f10.2,A5)',iostat=ios) &
            & ' SrCO3 sink -- 88Sr        :',loc_sig,' o/oo'
       loc_tot  = loc_tot2_sedgrid-loc_tot2a_sedgrid-loc_tot2b_sedgrid
       if (abs(loc_tot) > const_real_nullsmall) then
          loc_frac     = loc_tot2b_sedgrid
          loc_standard = const_standardsR(ocn_type(io_Sr_88Sr))
          loc_sig      = fun_calc_isotope_deltaR(loc_tot,loc_frac,loc_standard,const_real_null)
       else
          loc_sig = -999.9
       end if
       write(unit=out,fmt='(A28,f10.2,A5)',iostat=ios) &
            & ' Sr source -- 88Sr         :',loc_sig,' o/oo'
       Write(unit=out,fmt=*) '---------------------------------'
    end if
    ! -------------------------------------------------------- !
    ! GLOBAL GRID DIAGNOSTICS
    ! -------------------------------------------------------- !
    Write(unit=out,fmt=*) ' '
    Write(unit=out,fmt=*) '--- TOTAL SEDIMENT GRID ---------'
    Write(unit=out,fmt=*) '--- (as seen by BIOGEM) ---------'
    Write(unit=out,fmt=*) ' '
    ! MISC
    Write(unit=out,fmt=*) '---------------------------------'
    write(unit=out,fmt='(A28,I6)',iostat=ios) &
         & ' Total # sediment grid pts :',int(sum(loc_mask(:,:)))
    call check_iostat(ios,__LINE__,__FILE__)
    write(unit=out,fmt='(A28,e14.6,A3)',iostat=ios) &
         & ' Total sediment area       :',sum(loc_mask(:,:)*phys_sed(ips_A,:,:)),'m2'
    call check_iostat(ios,__LINE__,__FILE__)
    ! local variables 
    loc_tot_mask_area = sum(loc_mask(:,:)*loc_area(:,:))
    if (loc_tot_mask_area > const_real_nullsmall) then 
       loc_mean_sedgrid = sum(loc_mask(:,:)*loc_sed_coretop(is_CaCO3,:,:)*loc_area(:,:))/loc_tot_mask_area
    else
       loc_mean_sedgrid = 0.0
    end if
    Write(unit=out,fmt=*) '---------------------------------'
    write(unit=out,fmt='(A28,f6.2,A2)',iostat=ios) &
         & ' Mean wt% CaCO3            :',loc_mean_sedgrid,' %'
    call check_iostat(ios,__LINE__,__FILE__)
    Write(unit=out,fmt=*) '---------------------------------'
    ! write out mean global detrital flux in units of g cm-2 kyr-1 as a check 
    ! NOTE: local fluxes (loc_fsed) in units of mol cm-2 yr-1
    loc_mean_sedgrid = sum(loc_mask(:,:)*loc_area(:,:)*loc_fsed(is_det,:,:))/loc_tot_mask_area
    loc_mean_sedgrid = loc_mean_sedgrid/(conv_det_g_mol*conv_yr_kyr)
    write(unit=out,fmt='(A28,f7.3,A13)',iostat=ios) &
         & ' Mean detrital flux        :',loc_mean_sedgrid,' g cm-2 kyr-1'
    call check_iostat(ios,__LINE__,__FILE__)
    Write(unit=out,fmt=*) '---------------------------------'
    ! -------------------------------------------------------- !
    ! WEATHERING PARAMETER CALCULATIONS
    ! -------------------------------------------------------- !
    ! NOTE: originally assumed are:
    !       (1) CaSiO3:CaCO3 weathering is assumed to be in a 2:3 proportion
    !       (2) CO2 outgassing d13C is assumed to be -6 o/oo
    !       (3) P:ALK is -16.0 (but really it should be zero)
    !       this is becasue the adjustable parameters are held by ROKGEM (or BIOGEM) and not accessible to SEDGEM
    !       NOW, 2 SEDGEM parameters allow different choices to be made in the open system flux balance calculation:
    !       (1) par_sed_diag_Siweatheringfrac
    !       (2) par_sed_diag_volcanicd13C
    !       and then one further one to complete the ALK mass balance
    !       (3) par_sed_diag_P2ALK
    ! sum global total burial fluxes
    loc_tot_FCaCO3 = loc_deep_FCaCO3+loc_mud_FCaCO3+loc_reef_FCaCO3
    loc_tot_FPOC   = loc_deep_FPOC+loc_mud_FPOC+loc_reef_FPOC
    loc_tot_FPOP   = loc_deep_FPOP+loc_mud_FPOP+loc_reef_FPOP
    ! calculate mean CaCO3 burial d13C
    if (loc_tot_FCaCO3 > const_real_nullsmall) then
       loc_tot_FCaCO3_d13C = &
            & (loc_deep_FCaCO3_d13C*loc_deep_FCaCO3+loc_mud_FCaCO3_d13C*loc_mud_FCaCO3+loc_reef_FCaCO3_d13C*loc_reef_FCaCO3)/ &
            & loc_tot_FCaCO3
    else
       loc_tot_FCaCO3_d13C = 0.0
    end if
    ! calculate mean POC burial d13C
    if (loc_tot_FPOC > const_real_nullsmall) then
       loc_tot_FPOC_d13C = &
            & (loc_deep_FPOC_d13C*loc_deep_FPOC+loc_mud_FPOC_d13C*loc_mud_FPOC+loc_reef_FPOC_d13C*loc_reef_FPOC)/ &
            & loc_tot_FPOC
    else
       loc_tot_FPOC_d13C = 0.0
    end if
    ! calculate gamma
    if ( (abs((loc_tot_FPOC_d13C-loc_tot_FCaCO3_d13C)) > const_real_nullsmall) .AND. (loc_tot_FPOC > const_real_nullsmall) ) then
       loc_gamma = (par_sed_diag_volcanicd13C-loc_tot_FCaCO3_d13C)/(loc_tot_FPOC_d13C-loc_tot_FCaCO3_d13C)
    else
       loc_gamma = 0.0
    end if
    ! outgassing and ketogen carbon fluxes
    loc_Foutgassing = par_sed_diag_fracSiweath*loc_tot_FCaCO3/(1.0-loc_gamma)
    loc_Fkerogen    = loc_tot_FPOC - loc_gamma*par_sed_diag_fracSiweath*loc_tot_FCaCO3/(1.0-loc_gamma)
    ! create simple mass balance if no organic carbon cycle is active
    ! NOTE: loc_FCaCO3_d13C*(1.0-par_sed_diag_fracSiweath)*loc_tot_FCaCO3 + par_sed_diag_volcanicd13C*loc_Foutgassing = 
    !       loc_tot_FCaCO3_d13C*loc_tot_FCaCO3
    if (loc_tot_FPOC > const_real_nullsmall) then
       loc_FCaCO3_d13C = loc_tot_FCaCO3_d13C
    elseif ((1.0-par_sed_diag_fracSiweath)*loc_tot_FCaCO3 < const_real_nullsmall) then
       loc_FCaCO3_d13C = 0.0
    else
       loc_FCaCO3_d13C = (loc_tot_FCaCO3_d13C*loc_tot_FCaCO3 - par_sed_diag_volcanicd13C*loc_Foutgassing)/ &
            & ((1.0-par_sed_diag_fracSiweath)*loc_tot_FCaCO3)
    end if
    ! for updating kerogen weathering O2:C -- calculate O2 gain (mol yr-1) associated with the bural of reduced species
    ! NOTE: the value of sg_par_sed_diag_red_POC_O2 excludes the contribution from P (and N etc.)
    loc_tot_FO2 = (-par_sed_diag_C2O2)*loc_tot_FPOC + 2.0*loc_tot_FPOP
    ! -------------------------------------------------------- !
    ! WRITE WEATHERING PARAMETERS
    ! -------------------------------------------------------- !
    Write(unit=out,fmt=*) ' '
    Write(unit=out,fmt=*) '--- DERIVED PARAMETER VALUES ----'
    Write(unit=out,fmt=*) ' '
    Write(unit=out,fmt=*) '---------------------------------'
    write(unit=out,fmt='(A28,e14.6,A9)',iostat=ios) &
         & ' Global CaCO3 burial       :',loc_tot_FCaCO3,' mol yr-1'
    write(unit=out,fmt='(A28,e14.6,A9,A16,f8.6,A29)',iostat=ios) &
         & ' silicate weathering       :',par_sed_diag_fracSiweath*loc_tot_FCaCO3,' mol yr-1', &
         & ' for an assumed ',par_sed_diag_fracSiweath,' silicate weathering fraction'
    write(unit=out,fmt='(A28,e14.6,A9,A16,f8.6,A30)',iostat=ios) &
         & ' carbonate weathering      :',(1.0-par_sed_diag_fracSiweath)*loc_tot_FCaCO3,' mol yr-1', &
         & ' for an assumed ',(1.0-par_sed_diag_fracSiweath),' carbonate weathering fraction'
    write(unit=out,fmt='(A28,e14.6,A9)',iostat=ios) &
         & ' Global POC burial         :',loc_tot_FPOC,' mol yr-1'
    write(unit=out,fmt='(A28,e14.6,A9)',iostat=ios) &
         & ' Global POP burial         :',loc_tot_FPOP,' mol yr-1'
    Write(unit=out,fmt=*) '---------------------------------'
    write(unit=out,fmt='(A28,e14.6,A5)',iostat=ios) &
         & ' Global CaCO3 d13C         :',loc_FCaCO3_d13C,' o/oo'
    write(unit=out,fmt='(A28,e14.6,A5)',iostat=ios) &
         & ' Global POC d13C           :',loc_tot_FPOC_d13C,' o/oo'
    write(unit=out,fmt='(A28,e14.6,A14,f7.3,A16)',iostat=ios) &
         & ' gamma                     :',loc_gamma,' @ an assumed ',par_sed_diag_volcanicd13C,' o/oo outgassing'
    Write(unit=out,fmt=*) '---------------------------------'
    write(unit=out,fmt='(A28,e14.6,A9)',iostat=ios) &
         & ' CO2 outgassing            :',loc_Foutgassing,' mol yr-1'
    write(unit=out,fmt='(A28,e14.6,A9)',iostat=ios) &
         & ' kerogen weathering        :',loc_Fkerogen,' mol yr-1'
    Write(unit=out,fmt=*) '---------------------------------'
    if (loc_tot_FCaCO3 > const_real_nullsmall .AND. par_sed_diag_fracSiweath > const_real_nullsmall) then
       write(unit=out,fmt='(A28,e14.6)',iostat=ios) &
            & ' kerogen C/silicate ratio  =',loc_Fkerogen/(par_sed_diag_fracSiweath*loc_tot_FCaCO3)
    else
       write(unit=out,fmt='(A32)',iostat=ios) ' kerogen C/silicate ratio  = NaN'
    end if
    if (loc_tot_FCaCO3 > const_real_nullsmall .AND. par_sed_diag_fracSiweath > const_real_nullsmall) then
       write(unit=out,fmt='(A28,e14.6)',iostat=ios) &
            & ' kerogen P/silicate ratio  =',loc_tot_FPOP/(par_sed_diag_fracSiweath*loc_tot_FCaCO3)
    else
       write(unit=out,fmt='(A32)',iostat=ios) ' kerogen P/silicate ratio  = NaN'
    end if
    if (loc_Fkerogen > const_real_nullsmall) then
       write(unit=out,fmt='(A28,e14.6)',iostat=ios) &
            & ' kerogen P/C ratio         =',loc_tot_FPOP/loc_Fkerogen
    else
       write(unit=out,fmt='(A32)',iostat=ios) ' kerogen P/C ratio         = NaN'
    end if
    if (loc_Fkerogen > const_real_nullsmall) then
       write(unit=out,fmt='(A28,e14.6)',iostat=ios) &
            & ' kerogen ALK/C ratio       =',par_sed_diag_P2ALK*loc_tot_FPOP/loc_Fkerogen
    else
       write(unit=out,fmt='(A32)',iostat=ios) ' kerogen ALK/C ratio       = NaN'
    end if
    if (loc_Fkerogen > const_real_nullsmall) then
       write(unit=out,fmt='(A28,e14.6)',iostat=ios) &
            & ' O2 consumption ratio      =',-loc_tot_FO2/loc_Fkerogen
    end if
    Write(unit=out,fmt=*) '---------------------------------'
    ! -------------------------------------------------------- ! write file footer
    Write(unit=out,fmt=*) ' '
    Write(unit=out,fmt=*) '================================='
    ! -------------------------------------------------------- ! close file
    Write(unit=out,fmt=*) ' '
    CLOSE(out,iostat=ios)
    call check_iostat(ios,__LINE__,__FILE__)

    ! *** SAVE WEATHERING PARAMETERS ***
    
    ! set filename
    loc_filename = TRIM(par_outdir_name)//'INFO_weathering_parameters_AT_END'//string_results_ext
    
    ! open file
    call check_unit(out,__LINE__,__FILE__)
    OPEN(out,file=TRIM(loc_filename),action='write',iostat=ios)
    call check_iostat(ios,__LINE__,__FILE__)
    
    ! write out section of text for copy-paste into a user-config
    Write(unit=out,fmt=*) '# --- ROKGEM USER-CONFIG --------'
    Write(unit=out,fmt=*) '#'
    Write(unit=out,fmt=*) '# NOTE: automatically generated by SEDGEM'
    Write(unit=out,fmt=*) '# NOTE: this calculation assumes:'
    write(unit=out,fmt=*) '#       silicate weathering fraction (sg_par_sed_diag_fracSiweath)     == ',par_sed_diag_fracSiweath
    write(unit=out,fmt=*) '#       volcanic outgassing d13C (sg_par_sed_diag_volcanicd13C)        == ',par_sed_diag_volcanicd13C
    write(unit=out,fmt=*) '#       implicit P:ALK in OM N transformation (sg_par_sed_diag_P2ALK)  == ',par_sed_diag_P2ALK
    Write(unit=out,fmt=*) '# NOTE: BE CAREFUL -- the values of parameters #2 and #3 are duplicated in other modules ...'
    Write(unit=out,fmt=*) '#       un-comment the following lines to ensure alignment:'
    Write(unit=out,fmt=*) '#rg_par_outgas_CO2_d13C=',par_sed_diag_volcanicd13C
    Write(unit=out,fmt=*) '#bg_par_bio_red_PON_ALK=',par_sed_diag_P2ALK/16.0
    Write(unit=out,fmt=*) '#'
    Write(unit=out,fmt=*) '# set an OPEN system'
    Write(unit=out,fmt=*) 'bg_ctrl_force_sed_closedsystem=.FALSE.'
    Write(unit=out,fmt=*) '# set CaCO3_weathering-temperature feedback'
    Write(unit=out,fmt=*) 'rg_opt_weather_T_Ca=.true.'
    Write(unit=out,fmt=*) '# set CaSiO3_weathering-temperature feedback'
    Write(unit=out,fmt=*) 'rg_opt_weather_T_Si=.true.'
    Write(unit=out,fmt=*) '# weathering reference mean global land air surface temperature (oC)'
    Write(unit=out,fmt=*) '# NOTE: this value can also be obtained from BIOGEM time-series: biogem_series_misc_SLT.res'
    Write(unit=out,fmt=*) 'rg_par_ref_T0=',dum_SLT
    Write(unit=out,fmt=*) '# global carbonate weathering rate (mol Ca2+ yr-1)'
    write(unit=out,fmt=*) 'rg_par_weather_CaCO3=',(1.0-par_sed_diag_fracSiweath)*loc_tot_FCaCO3
    Write(unit=out,fmt=*) '#  global silicate weathering rate (mol Ca2+ yr-1)'
    Write(unit=out,fmt=*) 'rg_par_weather_CaSiO3=',par_sed_diag_fracSiweath*loc_tot_FCaCO3
    Write(unit=out,fmt=*) '# CO2 outgassing rate (mol C yr-1)'
    Write(unit=out,fmt=*) 'rg_par_outgas_CO2=',loc_Foutgassing
    Write(unit=out,fmt=*) '# set isotopic value of CO2 outgassing (assumed) (o/oo)'
    Write(unit=out,fmt=*) 'rg_par_outgas_CO2_d13C=',par_sed_diag_volcanicd13C
    Write(unit=out,fmt=*) '# set isotopic value of carbonate weathering (o/oo)'
    Write(unit=out,fmt=*) 'rg_par_weather_CaCO3_d13C=',loc_FCaCO3_d13C
    Write(unit=out,fmt=*) '# kerogen POC weathering ratio (to silicate Ca2+) and isotopic composiiton'
    if ((par_sed_diag_fracSiweath*loc_tot_FCaCO3) > const_real_nullsmall) then
       Write(unit=out,fmt=*) 'rg_par_weather_CaSiO3_fracC=',loc_Fkerogen/(par_sed_diag_fracSiweath*loc_tot_FCaCO3)
       Write(unit=out,fmt=*) 'rg_par_weather_CaSiO3_fracC_d13C=',loc_tot_FPOC_d13C
    else
       Write(unit=out,fmt=*) '# WARNING: could not be calculated'
       Write(unit=out,fmt=*) 'rg_par_weather_CaSiO3_fracC=0.0'
       Write(unit=out,fmt=*) 'rg_par_weather_CaSiO3_fracC_d13C=0.0'
    end if
    Write(unit=out,fmt=*) '# kerogen POP weathering settings (as ratios with kerogen C)'
    Write(unit=out,fmt=*) '# NOTE: also account for implicitly associated ALK which is (sg_par_sed_diag_P2ALK):'
    Write(unit=out,fmt=*) '#       the PO4 flux * ',par_sed_diag_P2ALK
    Write(unit=out,fmt=*) '#       (but really ths should simply be zero when simulating geologic carbon cycling)'
    if (loc_Fkerogen > const_real_nullsmall) then
       Write(unit=out,fmt=*) 'rg_par_weather_kerogen_fracP=',loc_tot_FPOP/loc_Fkerogen
       Write(unit=out,fmt=*) 'rg_par_weather_kerogen_fracALK=',par_sed_diag_P2ALK*loc_tot_FPOP/loc_Fkerogen
    else
       Write(unit=out,fmt=*) '# WARNING: could not be calculated'
       Write(unit=out,fmt=*) 'rg_par_weather_kerogen_fracP=0.0'
       Write(unit=out,fmt=*) 'rg_par_weather_kerogen_fracALK=0.0'
    end if
    if (ctrl_sed_diag_balanceO2 .AND. (loc_Fkerogen > const_real_nullsmall)) then
       Write(unit=out,fmt=*) '# updated kerogen O2 weathering consumption ratio to balance global (Corg) O2 budget'
       Write(unit=out,fmt=*) 'rg_par_weather_kerogen_fracO2=',-loc_tot_FO2/loc_Fkerogen
    end if
    Write(unit=out,fmt=*) '#'
    Write(unit=out,fmt=*) '# -------------------------------'
    Write(unit=out,fmt=*) ''
    ! optional/additional SEDGEM parameters
    Write(unit=out,fmt=*) '# --- SEDGEM USER-CONFIG --------'
    Write(unit=out,fmt=*) '#'
    Write(unit=out,fmt=*) '# scale factor for reefal precipitation rate (mol cm-2 yr-1) to achieve a global burial rate of: ', &
         & 1.0E-12*par_sed_CaCO3burialTOT,' (Tmol yr-1):'
    Write(unit=out,fmt=*) 'sg_par_sed_reef_CaCO3precip_sf=',par_sed_reef_CaCO3precip_sf
    Write(unit=out,fmt=*) '#'
    Write(unit=out,fmt=*) '# -------------------------------'
    Write(unit=out,fmt=*) ''
    ! close file
    CLOSE(out,iostat=ios)
    call check_iostat(ios,__LINE__,__FILE__)

    ! *** SAVE FULL CORE-TOP DATA IN TEXT FILE FORMAT ***
    
    ! set filename
    loc_filename = TRIM(par_outdir_name)//'INFO_sediment_locations_AT_END'//string_results_ext
    
    ! open file
    call check_unit(out,__LINE__,__FILE__)
    OPEN(out,file=TRIM(loc_filename),action='write',iostat=ios)
    call check_iostat(ios,__LINE__,__FILE__)

    ! write data
    Write(unit=out,fmt=*) '% ========================================'
    Write(unit=out,fmt=*) '% Sediment diagnostics data'
    Write(unit=out,fmt=*) '% ========================================'
             write(unit=out,fmt='(A2,32A12)',iostat=ios)          &
                  & ' %',                                          &
                  & '           i',                               &
                  & '           j',                               &
                  & '         lon',                               &
                  & '         lat',                               &
                  & ' ocean depth',                               &
                  & '        temp',                               &
                  & '         sal',                               &
                  & '       [DIC]',                               &
                  & '       [ALK]',                               &
                  & '        [O2]',                               &
                  & '        [Ca]',                               &
                  & '      [SiO4]',                               &
                  & '       [CO3]',                               &
                  & '    ohm(cal)',                               &
                  & '      d[CO3]',                               &
                  & '    POC conc',                               &
                  & '  CaCO3 conc',                               &
                  & '   opal conc',                               &
                  & '    det conc',                               &
                  & '    ash conc',                               &
                  & '    POC rain',                               &
                  & '  CaCO3 rain',                               &
                  & '   opal rain',                               &
                  & '    det rain',                               &
                  & '     POC dis',                               &
                  & '   CaCO3 dis',                               &
                  & '    opal dis',                               &
                  & '     det dis',                               &
                  & '  POC burial',                               &
                  & '   CaCO3 bur',                               &
                  & ' opal burial',                               &
                  & '  det burial'
             call check_iostat(ios,__LINE__,__FILE__)
             write(unit=out,fmt='(A2,32A12)',iostat=ios)          &
                  & ' %',                                          &
                  & '           -',                               &
                  & '           -',                               &
                  & '       deg E',                               &
                  & '       deg N',                               &
                  & '           m',                               &
                  & '        degC',                               &
                  & '       o/oo',                                &
                  & '          uM',                               &
                  & '          uM',                               &
                  & '          uM',                               &
                  & '          mM',                               &
                  & '          uM',                               &
                  & '          uM',                               &
                  & '           -',                               &
                  & '          uM',                               &
                  & '         wt%',                               &
                  & '         wt%',                               &
                  & '         wt%',                               &
                  & '         wt%',                               &
                  & '         wt%',                               &
                  & ' molcm-2yr-1',                               &
                  & ' molcm-2yr-1',                               &
                  & ' molcm-2yr-1',                               &
                  & '  gcm-2kyr-1',                               &
                  & ' molcm-2yr-1',                               &
                  & ' molcm-2yr-1',                               &
                  & ' molcm-2yr-1',                               &
                  & '  gcm-2kyr-1',                               &
                  & ' molcm-2yr-1',                               &
                  & ' molcm-2yr-1',                               &
                  & ' molcm-2yr-1',                               &
                  & '  gcm-2kyr-1'
    call check_iostat(ios,__LINE__,__FILE__)
    Write(unit=out,fmt=*) '% ----------------------------------------'
    Write(unit=out,fmt=*) '% (1) Sediment core locations'
    Write(unit=out,fmt=*) '% ----------------------------------------'
    DO i=1,n_i
       DO j=1,n_j
          IF (sed_mask(i,j) .AND. sed_save_mask(i,j)) THEN
             write(unit=out,fmt='(2X,2I12,2f12.1,3f12.1,6f12.1,2f12.3,5f12.1,12f12.3)',iostat=ios) &
                  & i, &
                  & j,                                                     &
                  & phys_sed(ips_lon,i,j),                                   &
                  & phys_sed(ips_lat,i,j),                                   &
                  & phys_sed(ips_D,i,j),                                     &
                  & dum_sfcsumocn(io_T,i,j) - const_zeroC,                                 &
                  & dum_sfcsumocn(io_S,i,j),                                 &
                  & 1.0E+06*dum_sfcsumocn(io_DIC,i,j),                       &
                  & 1.0E+06*dum_sfcsumocn(io_ALK,i,j),                       &
                  & 1.0E+06*dum_sfcsumocn(io_O2,i,j),                        &
                  & 1.0E+03*dum_sfcsumocn(io_Ca,i,j),                        &
                  & 1.0E+06*dum_sfcsumocn(io_SiO2,i,j),                      &
                  & 1.0E+06*sed_carb(ic_conc_CO3,i,j),                       &
                  & sed_carb(ic_ohm_cal,i,j),                                &
                  & 1.0E+06*sed_carb(ic_dCO3_cal,i,j),                       &
                  & loc_sed_coretop(is_POC,i,j),                             &
                  & loc_sed_coretop(is_CaCO3,i,j),                           &
                  & loc_sed_coretop(is_opal,i,j),                            &
                  & loc_sed_coretop(is_det,i,j),                             &
                  & loc_sed_coretop(is_ash,i,j),                             &
                  & 1.0E+06*sed_fsed(is_POC,i,j)/dum_dtyr,                   &
                  & 1.0E+06*sed_fsed(is_CaCO3,i,j)/dum_dtyr,                 &
                  & 1.0E+06*sed_fsed(is_opal,i,j)/dum_dtyr,                  &
                  & sed_fsed(is_det,i,j)/(conv_det_g_mol*(conv_yr_kyr*dum_dtyr)),                   &
                  & 1.0E+06*sed_fdis(is_POC,i,j)/dum_dtyr,                   &
                  & 1.0E+06*sed_fdis(is_CaCO3,i,j)/dum_dtyr,                 &
                  & 1.0E+06*sed_fdis(is_opal,i,j)/dum_dtyr,                  &
                  & sed_fdis(is_det,i,j)/(conv_det_g_mol*(conv_yr_kyr*dum_dtyr)), &
                  & 1.0E+06*(sed_fsed(is_POC,i,j)-sed_fdis(is_POC,i,j))/dum_dtyr,                   &
                  & 1.0E+06*(sed_fsed(is_CaCO3,i,j)-sed_fdis(is_CaCO3,i,j))/dum_dtyr,                 &
                  & 1.0E+06*(sed_fsed(is_opal,i,j)-sed_fdis(is_opal,i,j))/dum_dtyr,                  &
                  & (sed_fsed(is_det,i,j)-sed_fdis(is_det,i,j))/(conv_det_g_mol*(conv_yr_kyr*dum_dtyr))   
             call check_iostat(ios,__LINE__,__FILE__)
          end if
       end do
    end do
    Write(unit=out,fmt=*) '% ----------------------------------------'
    Write(unit=out,fmt=*) '% (2) All other sediemnt locations (excluding sedcores)'
    Write(unit=out,fmt=*) '% ----------------------------------------'
    DO i=1,n_i
       DO j=1,n_j
          IF (sed_mask(i,j) .AND. (.NOT. sed_save_mask(i,j))) THEN
             write(unit=out,fmt='(2X,2I12,2f12.1,3f12.1,6f12.1,2f12.3,5f12.1,12f12.3)',iostat=ios) &
                  & i, &
                  & j,                                                     &
                  & phys_sed(ips_lon,i,j),                                   &
                  & phys_sed(ips_lat,i,j),                                   &
                  & phys_sed(ips_D,i,j),                                     &
                  & dum_sfcsumocn(io_T,i,j) - const_zeroC,                                 &
                  & dum_sfcsumocn(io_S,i,j),                                 &
                  & 1.0E+06*dum_sfcsumocn(io_DIC,i,j),                       &
                  & 1.0E+06*dum_sfcsumocn(io_ALK,i,j),                       &
                  & 1.0E+06*dum_sfcsumocn(io_O2,i,j),                        &
                  & 1.0E+03*dum_sfcsumocn(io_Ca,i,j),                        &
                  & 1.0E+06*dum_sfcsumocn(io_SiO2,i,j),                      &
                  & 1.0E+06*sed_carb(ic_conc_CO3,i,j),                       &
                  & sed_carb(ic_ohm_cal,i,j),                                &
                  & 1.0E+06*sed_carb(ic_dCO3_cal,i,j),                       &
                  & loc_sed_coretop(is_POC,i,j),                             &
                  & loc_sed_coretop(is_CaCO3,i,j),                           &
                  & loc_sed_coretop(is_opal,i,j),                            &
                  & loc_sed_coretop(is_det,i,j),                             &
                  & loc_sed_coretop(is_ash,i,j),                             &
                  & 1.0E+06*sed_fsed(is_POC,i,j)/dum_dtyr,                   &
                  & 1.0E+06*sed_fsed(is_CaCO3,i,j)/dum_dtyr,                 &
                  & 1.0E+06*sed_fsed(is_opal,i,j)/dum_dtyr,                  &
                  & sed_fsed(is_det,i,j)/(conv_det_g_mol*(conv_yr_kyr*dum_dtyr)),                   &
                  & 1.0E+06*sed_fdis(is_POC,i,j)/dum_dtyr,                   &
                  & 1.0E+06*sed_fdis(is_CaCO3,i,j)/dum_dtyr,                 &
                  & 1.0E+06*sed_fdis(is_opal,i,j)/dum_dtyr,                  &
                  & sed_fdis(is_det,i,j)/(conv_det_g_mol*(conv_yr_kyr*dum_dtyr)), &
                  & 1.0E+06*(sed_fsed(is_POC,i,j)-sed_fdis(is_POC,i,j))/dum_dtyr,                   &
                  & 1.0E+06*(sed_fsed(is_CaCO3,i,j)-sed_fdis(is_CaCO3,i,j))/dum_dtyr,                 &
                  & 1.0E+06*(sed_fsed(is_opal,i,j)-sed_fdis(is_opal,i,j))/dum_dtyr,                  &
                  & (sed_fsed(is_det,i,j)-sed_fdis(is_det,i,j))/(conv_det_g_mol*(conv_yr_kyr*dum_dtyr))                   
             call check_iostat(ios,__LINE__,__FILE__)
          end if
       end do
    end do
    Write(unit=out,fmt=*) '% ========================================'
    
    ! close file
    CLOSE(out,iostat=ios)
    call check_iostat(ios,__LINE__,__FILE__)

    ! ---------------------------------------------------------- !
    ! END
    ! ---------------------------------------------------------- !
  end SUBROUTINE sub_data_save_seddiag_GLOBAL
  ! ****************************************************************************************************************************** !


END MODULE sedgem_data

