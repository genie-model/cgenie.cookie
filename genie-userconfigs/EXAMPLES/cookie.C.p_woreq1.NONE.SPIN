# **********************************************************************
# copied from: 	LAB_1.SPIN
# base-config: 	cookie.C.p_woreq1.NONE
# **********************************************************************
#
# *** CLIMATE **********************************************************
#
# solar constant reduced by 5.6% for late Neoproterozoic: 1285.92 W m-2
ma_genie_solar_constant=1285.92
# set NO CO2 climate feedback (climate not responding to changing pCO2)
ea_36=n
# scaling for atmospheric CO2 radiative forcing, relative to 278 ppm
ea_radfor_scl_co2=20.0
#
# *** BIOLOGICAL NEW PRODUCTION ****************************************
#
# No biology
bg_par_bio_prodopt='NONE'
#
# *** DATA SAVING ******************************************************
#
# Preformed diagnostics option == BASIC
bg_par_data_save_level=2
# force time-slice save at run end only
bg_par_infile_slice_name='save_timeslice_NONE.dat'
# Include experiment name in BIOGEM netCDF filenames? 
bg_ctrl_ncout_expid_name=.false.
#
# *** FORCINGS *********************************************************
#
# add a geothermal heat flux (W m-2)
bg_ctrl_force_GOLDSTEInTS=.TRUE.
bg_par_Fgeothermal=100.0E-3
#
# *** MISC *************************************************************
#
#
# **********************************************************************
# *** END **************************************************************
# **********************************************************************

