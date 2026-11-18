# **********************************************************************
# base-config: 	cookie.C.ebP2_yyy.BASESCH4
# **********************************************************************
#
# *** CLIMATE **********************************************************
#
# set NO climate feedback (climate NOT responding to changing pCH4)
ea_36=n
#
# *** ATMOSPHERIC PHOTOCHEMISTRY ***************************************
#
# atmospheric photochemical scheme ID string
ac_par_atm_CH4_photochem="claire06"
#
# *** BIOLOGICAL NEW PRODUCTION ****************************************
#
bg_par_bio_prodopt='1N1T_PO4MM_Tdep'
bg_par_bio_mu1=10
bg_par_bio_c0_PO4=8.9369e-007
#
# *** ORGANIC MATTER EXPORT RATIOS *************************************
#
# production fraction of dissolved organic matter
bg_par_bio_red_DOMfrac=0.66
#
# *** INORGANIC MATTER EXPORT RATIOS ***********************************
#
# SET NO PELAGIC CARBONATE PRODUCTION
bg_par_bio_red_POC_CaCO3=0.0
bg_par_bio_red_POC_CaCO3_pP=0.0
#
# *** REMINERALIZATION *************************************************
#
# DOC lifetime (yrs)
bg_par_bio_remin_DOMlifetime=0.5
# initial fractional abundance of POC component #2
bg_par_bio_remin_POC_frac2=0.0557
# depth of remineralization or particulate organic matter
bg_par_bio_remin_POC_eL1=294.9726
# remineralization length #2 for POC
bg_par_bio_remin_POC_eL2=1000000.0
# initial fractional abundance of CaCO3 component #2
bg_par_bio_remin_CaCO3_frac2=0.45
# depth of remineralization or CaCO3
bg_par_bio_remin_CaCO3_eL1=1.8905e+003
# remineralization length #2 for CaCO3
bg_par_bio_remin_CaCO3_eL2=1000000.0
# sulfide oxidation rate form
bg_opt_bio_remin_oxidize_H2StoSO4="OLD"
# sulfide oxidation rate constant [only used if rate law is linear]
#bg_par_bio_remin_kH2StoSO4=1.0E07
# rate constant for aerobic methanotrophy
bg_par_bio_remin_AER_kAER=0.10
# O2 half-saturation constant for aerobic methanotrophy
bg_par_bio_remin_AER_Km_O2=20.E-06
# rate constant for AOM
bg_par_bio_remin_AOM_kAOM=0.01
# SO4 half-saturation constant for AOM
bg_par_bio_remin_AOM_Km_SO4=500.0E-06
# enable/disable thermodynamic terms
bg_par_bio_remin_AER_thermo="on"
bg_par_bio_remin_AOM_thermo="on"
# biological energy quanta (BEQ) for CH4 oxidation
bg_par_bio_remin_AER_BEQ=15.0
bg_par_bio_remin_AOM_BEQ=15.0
# methanogenesis 13C fractionation
bg_par_d13C_Corg_CH4_epsilon=-35.0
#
bg_ctrl_bio_remin_POC_fixed=.false.
bg_par_bio_remin_sinkingrate=125.0
bg_par_bio_remin_POC_K1=9.0E11
bg_par_bio_remin_POC_Ea1=54000.0
bg_par_bio_remin_POC_K2=1.0E14
bg_par_bio_remin_POC_Ea2=80000.0
bg_par_bio_remin_POC_frac2=0.008
#
# *** DATA SAVING ******************************************************
#
# BASIC + biology + FULL (inc. redox) geochem diagnostics
bg_par_data_save_level=15
# force time-slice save at run end only
bg_par_infile_slice_name='save_timeslice_NONE.dat'
#
# *** FORCINGS *********************************************************
#
# SPEFICY ANY FORCINGS OF OCEAN/ATMOSPHERE/SEDIMENTS
# e.g., bg_par_forcing_name="pyyyyz.RpCO2_Rp13CO2" -> generic forcing of atmopsheric pCO2 and d13C
bg_par_forcing_name="pyyyyz.RpCO2_Rp13CO2"
# SCALE THE FORCING
# NOTE: Crichton et al. [2021] used 280 ppm rather than 278 ppm
bg_par_atm_force_scale_val_3=280.0E-06
bg_par_atm_force_scale_val_4=-6.5
#
# *** MISC *************************************************************
#
# *** optional changes to align with PALEO recommendations **********
# set 'instantaneous' water column remineralziation
bg_par_bio_remin_sinkingrate_physical=9.9E9
bg_par_bio_remin_sinkingrate_reaction=125.0
# maximum time-scale to geochemical reaction completion (days)
bg_par_bio_geochem_tau=45.0
#
# *** S / CH4 cycling settings *****************************************
#
# set remin rate scaling parameters
bg_par_bio_remin_k_O2=1.000
bg_par_bio_remin_k_NO3=0.990
bg_par_bio_remin_k_FeOOH=0.600
bg_par_bio_remin_k_SO4=0.109
bg_par_bio_remin_k_meth=0.050 
# set tuned redox parameters
bg_par_bio_remin_c0_O2=8e-06
bg_par_bio_remin_ci_O2=4e-06       
bg_par_bio_remin_c0_NO3=6e-05  
bg_par_bio_remin_ci_NO3=3.75e-06 
# cookie S-CH4 cycling:
bg_opt_bio_remin_oxidize_H2StoSO4='2ndorder'
bg_par_bio_remin_kH2StoSO4_perM2perhr=0.625E6
bg_par_bio_remin_c0_SO4=5.0E-04
bg_par_bio_remin_ci_SO4=1.0E-03
bg_ctrl_bio_remin_AOM_OLD=.false.
bg_ctrl_bio_remin_redox_save=.true.
#
# *** BOUNDARY CONDITIONS **********************************************
#
# SPECIFY SOLAR CONSTANT (Wm-2)
# NOTE: default (modern) value is 1368.0
ma_genie_solar_constant=1285.92
# SET THE CO2 FORCING
# NOTE: setting pCO2 here will over-ride any previous setting
#       (e.g., the 280.0E-06 default Crichton et al. [2021] value above)
bg_par_atm_force_scale_val_3=278.0E-06
#
# **********************************************************************
# *** END **************************************************************
# **********************************************************************

