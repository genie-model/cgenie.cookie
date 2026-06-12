# **********************************************************************
# base-config: 	cookie.CB.p_worbe2.BASES
# **********************************************************************
#
# *** CLIMATE **********************************************************
#
# SPECIFY WHETHER RADIATIVE FORCING IS DERIVED FROM ATMOSPHERIC CO2	
# NOTE: 'y' requires that CO2 is being calculated in the first place ...
#       'n' can be used for a configuration either with or without explicitly calculated atmopsheric CO2
# e.g., ea_36=n -> raditive forcing is set by the value of ea_radfor_scl_co2
ea_36=y
#
# *** BIOLOGICAL NEW PRODUCTION ****************************************
#
# SPECIFY THE BIOLOGICAL EXPORT SCHEME
# e.g., bg_par_bio_prodopt='1N1T_PO4MM' -> as used in Ridgwell et al. [2007]
bg_par_bio_prodopt='1N1T_PO4MM'
# SPECIFY BIOLOGICAL EXPORT PARAMETERS
# maximum rate of conversion of dissolved PO4 into organic matter by phytoplankton (mol kg-1 yr-1)
bg_par_bio_k0_PO4=1.9582242E-06
# [PO4] M-M half-sat value (mol kg-1)
bg_par_bio_c0_PO4=2.1989611E-07
#
# *** EXPORT RATIOS ****************************************************
#
# SPECIFY BULK CARBON EXPORT (AND/OR REDFIELD) RATIOS
# production fraction of dissolved organic matter
bg_par_bio_red_DOMfrac=0.66
# fix CaCO3:POC export rain ratio
bg_ctrl_force_CaCO3toPOCrainratio=.true.
bg_par_CaCO3toPOCrainratio_file='CaCO3toPOC.muffin.CB.p_worbe2.BASES.ridgwelletal.SPIN.dat'
#
# *** REMINERALIZATION *************************************************
#
# SPECIFY PARAMETERS CONTROLLIING WATER-COLUMN DECAY OF DISSOLVED ORGANIC MATTER (DOM)
# DOC lifetime (yrs)
bg_par_bio_remin_DOMlifetime=0.5
# SPECIFY PARAMETERS CONTROLLIING WATER-COLUMN DECAY OF PARTICULATE ORGANIC MATTER (POM)
# initial fractional abundance of POC component #2
bg_par_bio_remin_POC_frac2=6.4591110E-02
# depth of remineralization or particulate organic matter
bg_par_bio_remin_POC_eL1=550.5195
# remineralization length #2 for POC
bg_par_bio_remin_POC_eL2=1000000.0
# SPECIFY PARAMETERS CONTROLLIING WATER-COLUMN DECAY OF CaCO3
# initial fractional abundance of CaCO3 component #2
bg_par_bio_remin_CaCO3_frac2=0.4325
# depth of remineralization or CaCO3
bg_par_bio_remin_CaCO3_eL1=1083.486
# remineralization length #2 for CaCO3
bg_par_bio_remin_CaCO3_eL2=1000000.0
#
# *** DATA SAVING ******************************************************
#
# SPECIFY WHAT TO SAVE!
# NOTE: (0) REFER TO THE COOKIE MANUAL FOR DETAILS!
#		(1) all 'false' will leave you with the basic (hidden option!) defaults: 
#           bg_ctrl_save_hidden_grid=.true., bg_ctrl_save_hidden_climate=.true.
#		(2) isotope output is enabled by turning on bg_ctrl_save_basic_proxies
#		(3) not selecting a biologicla scheme (or selection ECOGEM) will turn off bg_ctrl_save_basic_biologicalpump
#       (4) there are 4 'advanced' options corresponding to the 4 'basic' options which add more output
#		(5) turning off a 'basic' option will also turn off the corresponding advanced option
#		(6) bg_ctrl_save_basic_ALL=.true. turns on all basic options; 
#           bg_ctrl_save_advanced_ALL=.true. turns on all advanced options
bg_ctrl_save_basic_reservoirs=.TRUE.
bg_ctrl_save_basic_geochemistry=.TRUE.
bg_ctrl_save_basic_biologicalpump=.TRUE.
bg_ctrl_save_basic_proxies=.false.
bg_ctrl_save_advanced_geochemistry=.true.
# SPECIFY THE INTERVAL AT WHICH TIME-SERIES DATA IS SAVED
# e.g., bg_par_data_save_sig_timeinterval=10.0 -> save time-series data every 10.0 years
# NOTE: to employ a file of save years in determining when time-series data is saved, set:
#       bg_par_data_save_sig_timeinterval=0.0 
bg_par_data_save_sig_timeinterval=10.0
# SPECIFY THE INTERVAL AT WHICH TIME-SLICE (SPATIAL) DATA IS SAVED
# e.g., bg_par_data_save_slice_timeinterval=1000.0 -> save time-slice data every 1000.0 years
# NOTE: in order to employ a file of save years to determine when time-slice data is saved, set:
#       bg_par_data_save_slice_timeinterval=0.0 
#       and then to save time-slice data only for the last year of the experiment, set: 
#       bg_par_infile_slice_name='save_timeslice_NONE.dat'
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
# *** Mg/Ca SETTINGS ************************************************
#
# pH solution criteria (REVISED)
gm_ctrl_carbchem_pH_OLD=.FALSE.
# pH solution tolerance (pH units) (DEFAULT)
gm_par_carbchem_pH_tolerance=0.001
# pH solution tolerance for calculating buffering (DEFAULT)
gm_par_carbchem_pH_tolerance_buffering=0.0001
# factor relaxation in pH solution tolerance for re-try (DEFAULT)
gm_par_carbchem_dpH_tolerance=10.0
# pH solution maximum number of iterations (DEFAULT)
gm_par_carbchem_pH_iterationmax=10
# basic carb constant set (DEFAULT)
gm_par_carbconstset_name='Mehrbach'
# NO Mg/Ca correction
gm_par_adj_carbconst_option='NONE'
# carbchem constant T and S range -- Mehrbach (DEFAULT)
gm_par_carbchem_Tmin =  2.0
gm_par_carbchem_Tmax = 35.0
gm_par_carbchem_Smin = 26.0
gm_par_carbchem_Smax = 43.0
# 'Eocene' [Ca], [Mg]
bg_ocn_init_35=2.0E-02
bg_ocn_init_50=3.0E-02
#
# *******************************************************************
# *** END ***********************************************************
# *******************************************************************
#
