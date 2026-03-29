# **********************************************************************
# adapted from:	cookie.C.p_worjh2.rCARB.SPIN
# base-config: 	cookie.C.p_worjh2.rCARB
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
# e.g., bg_par_bio_prodopt='NONE' -> no biologial export scheme (or ECOGEM will be used)
bg_par_bio_prodopt='NONE'
#
# *** DATA SAVING ******************************************************
#
# (refer to Chapter: 'cookie model output')
# SPECIFY THE SCOPE OF DATA SAVING
# e.g., bg_par_data_save_level=7 -> selects BASIC + biology + tracer + proxy output
bg_par_data_save_level=7
# SPECIFY THE INTERVAL AT WHICH TIME-SERIES DATA IS SAVED
# e.g., bg_par_data_save_sig_timeinterval=10.0 -> save time-series data every 10.0 years
# NOTE: set bg_par_data_save_sig_timeinterval=0.0 
#       in order to employ a file of save years to determine when time-series data is saved
bg_par_data_save_sig_timeinterval=10.0
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
# *** MISC *************************************************************
#
# SPECIFY THAT THE 'RED' NUMERICAL TRACER IS USED TO CALCULATE OCEAN VENTILLATION AGE
# NOTE: the 'red' numerical tracer must be selected in the base-config
bg_ctrl_force_ocn_age1=.true.
#
# **********************************************************************
# *** END **************************************************************
# **********************************************************************

