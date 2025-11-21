# **********************************************************************
# adapted from:	USERCONFIG.PALEO.ECOGEM.PO4.SPIN
# base-config: 	ANY cookie.CBE.*.BASES
# **********************************************************************
#
# *** CLIMATE **********************************************************
#
# SPECIFY WHETHER RADIATIVE FORCING IS DERIVED FROM ATMOSPHERIC CO2	
# NOTE: 'y' requires that CO2 is being calculated in the first place ...
#       'n' can be used for a configuration either with or without explicitly calculated atmopsheric CO2
# e.g., ea_36=y -> raditive forcing is set by the calculated value of atmopsheric CO2
ea_36=y
#
# *** BIOLOGICAL NEW PRODUCTION ****************************************
#
# SPECIFY THE BIOLOGICAL EXPORT SCHEME
# NOTE: with ECOGEM, no BIOGEM biological is selected
bg_par_bio_prodopt="NONE"
# ... if not an implicit biological export scheme, then ...
# SPECIFY THE FILE DEFINING THE ECOSYSTEM STRUCTURE 
# e.g., the full complexity [Ward et al. 2017] ecosystem: eg_par_ecogem_plankton_file='8P8Z.eco'
#       a single (small) phytoplankton: eg_par_ecogem_plankton_file ='NPD.eco'
eg_par_ecogem_plankton_file='8P8Z.eco'
#
# *** EXPORT RATIOS ****************************************************
#
# [dealt with by ECOGEM]
#
# *** REMINERALIZATION *************************************************
#
# SPECIFY PARAMETERS CONTROLLIING WATER-COLUMN DECAY OF DISSOLVED ORGANIC MATTER (DOM)
# DOC lifetime (yrs)
bg_par_bio_remin_DOMlifetime=0.5
# SPECIFY PARAMETERS CONTROLLIING WATER-COLUMN DECAY OF PARTICULATE ORGANIC MATTER (POM)
# initial fractional abundance of POC component #2
bg_par_bio_remin_POC_frac2=0.0557
# depth of remineralization or particulate organic matter
bg_par_bio_remin_POC_eL1=589.9451
# remineralization length #2 for POC
bg_par_bio_remin_POC_eL2=1000000.0
# SPECIFY PARAMETERS CONTROLLIING WATER-COLUMN DECAY OF CaCO3
# initial fractional abundance of CaCO3 component #2
bg_par_bio_remin_CaCO3_frac2=0.45
# depth of remineralization or CaCO3
bg_par_bio_remin_CaCO3_eL1=1.8905e+003
# remineralization length #2 for CaCO3
bg_par_bio_remin_CaCO3_eL2=1000000.0
#
# *** DATA SAVING ******************************************************
#
# (refer to Chapter: 'cookie model output')
# SPECIFY THE SCOPE OF DATA SAVING
# e.g., bg_par_data_save_level=10 -> includes full carbonate chemsitry output
bg_par_data_save_level=10
# SPECIFY THE INTERVAL AT WHICH TIME-SERIES DATA IS SAVED
# e.g., bg_par_data_save_sig_timeinterval=10.0 -> save time-series data every 10.0 years
# NOTE: set bg_par_data_save_sig_timeinterval=0.0 
#       in order to employ a file of save years to determine when time-series data is saved
bg_par_data_save_sig_timeinterval=1.0
# SPECIFY THE INTERVAL AT WHICH TIME-SLICE (SPATIAL) DATA IS SAVED
# e.g., bg_par_data_save_slice_timeinterval=1000.0 -> save time-slice data every 1000.0 years
# NOTE: set bg_par_data_save_slice_timeinterval=0.0 
#       in order to employ a file of save years to determine when time-slice data is saved
# NOTE: additionally set bg_par_infile_slice_name='save_timeslice_NONE.dat'
#       in order to force time-slice save at run end only
bg_par_data_save_slice_timeinterval=0.0
bg_par_infile_slice_name='save_timeslice_NONE.dat'
# SPECIFY USE OF MORE EXPLICIT (EXPERIMENT NAME BASED) netCDF FILENAMES
# e.g., bg_ctrl_ncout_expid_name=.false. -> assign default (short) netCDF filenames
bg_ctrl_ncout_expid_name=.false.
#
# *** FORCINGS *********************************************************
#
# SPEFICY ANY FORCINGS OF OCEAN/ATMOSPHERE/SEDIMENTS
# e.g., bg_par_forcing_name="pyyyyz.RpCO2_Rp13CO2" -> generic forcing of atmopsheric pCO2 and d13C
bg_par_forcing_name="pyyyyz.RpCO2_Rp13CO2"
# SCALE THE FORCING
bg_par_atm_force_scale_val_3=278.0E-06
bg_par_atm_force_scale_val_4=-6.5
#
# *** MISC **********************************************************
#
#
# *** BASIC ECOGEM SETUP ***********************************************
#
# NOTE: these parameters follow Ward et al. [2017]
# SPECIFY WHICH NUTRIENTS (AND Chla) ARE AVAILABLE AND WHETHER A QUOTA MODEL IS USED
# include cellular quota: PO4
eg_usePO4                   =.true.
eg_pquota                   =.true.
# include cellular quota: N
eg_useNO3                   =.false.
eg_nquota                   =.false.
# include cellular quota: Fe
eg_useFe                    =.false.
eg_fquota                   =.false.
# include cellular quota: Chlorophyll
eg_chlquota                 =.true.
# SPECIFY TUNED QUOTA PARAMETERS
# P and Fe quota tuned values
eg_qminP_a  =   3.33e-3
eg_qminFe_a =   1.00e-6
eg_qmaxFe_a =   4.00e-6
# SPECIFY ECOSYSTEM 'LOSS' TERMS
# Prey switching exponent (Holling 2 or 3 [= ns+1])
eg_ns=2
# Carbon respiration rate - intercept
eg_respir_a=0.0
# SPECIFY BULK CARBON EXPORT (AND/OR REDFIELD) RATIOS
# (1) DOM production from POM
eg_beta_mort_a =   0.8 # Maximum fraction to DOM as ESD --> zero
eg_beta_mort_b =   0.4 # Minimum fraction to DOM as ESD --> infinity 
eg_beta_mort_c = 100.0 # Size at 50:50 partition (default = 100 µm^3)
eg_beta_graz_a =   0.8 # Maximum fraction to DOM as ESD --> zero
eg_beta_graz_b =   0.4 # Minimum fraction to DOM as ESD --> infinity
eg_beta_graz_c = 100.0 # Size at 50:50 partition (default = 100 µm^3)
# (2) CaCO3:POC
# exponent for modifier of CaCO3:POC export ratio
# NOTE: set power to zero for invarient ratio
eg_par_bio_red_POC_CaCO3_pP = 0.0
# underlying export CaCO3 as a proportion of particulate organic matter (i.e., CaCO3/POC)
# NOTE: reduced by 50% on the basis that the modern ECOGEM scaling parameter is also ca. 50% lower
eg_par_bio_red_POC_CaCO3 = 0.0100 # underlying export CaCO3 as a proportion of particulate organic matter (i.e., CaCO3/POC)
# (3) d13C fractionation
# d13C fractionaton based on d13C of DIC?
eg_useDIC_13C =.true.
#
# *** ECOGEM BOUNDARY CONDITIONS ***************************************
#
# NOTE: these parameters adjust Ward et al. [2017] for the default ocean physics configuration
# SPEIFY PJHYSICAL ENVIRONMENT
# kraus-turner mixed layer scheme on (1) or off (0)
go_imld = 1
# set mixed layer to be only diagnosed (for ECOGEM)
go_ctrl_diagmld=.true.
# add seaice attenuation of PAR
eg_ctrl_PARseaicelimit=.true.
# ADJUST DOC/POC
# relative partitioning of C into DOM
# NOTE: this improves ocean oxygen by reducing oxidation demand from POC remin.
eg_par_beta_POCtoDOC=0.70
#
# *******************************************************************
# *** END ***********************************************************
# *******************************************************************

