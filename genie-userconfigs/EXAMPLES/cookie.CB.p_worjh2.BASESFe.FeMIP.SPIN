# **********************************************************************
# copied from: 	muffin.CB.p_worjh2.BASESFe.FeMIP.SPIN
# base-config: 	cookie.CB.p_worjh2.BASESFe
# **********************************************************************
#
# *** CLIMATE **********************************************************
#
# set climate feedback (climate responding to changing pCO2)
ea_36=y
#
# *** BIOLOGICAL NEW PRODUCTION ****************************************
#
# biological scheme ID string
bg_par_bio_prodopt="bio_PFe"
# [PO4] M-M half-sat value (mol kg-1)
bg_par_bio_c0_PO4=0.10E-06
# [Fe] M-M half-sat value (mol kg-1)
bg_par_bio_c0_Fe=0.10E-09
# biological uptake time-scale
bg_par_bio_tau=63.3827
#
# *** ORGANIC MATTER EXPORT RATIOS *************************************
#
# production fraction of dissolved organic matter
bg_par_bio_red_DOMfrac=0.66
#
# *** INORGANIC MATTER EXPORT RATIOS ***********************************
#
# underlying export CaCO3 as a proportion of particulate organic matter (i.e., CaCO3/POC)
bg_par_bio_red_POC_CaCO3=0.0485
# exponent for modifier of CaCO3:POC export ratio
bg_par_bio_red_POC_CaCO3_pP=0.7440
#
# *** REMINERALIZATION *************************************************
#
# DOC lifetime (yrs)
bg_par_bio_remin_DOMlifetime=1.0
# initial fractional abundance of POC component #2
bg_par_bio_remin_POC_frac2=0.0557
# depth of remineralization or particulate organic matter
bg_par_bio_remin_POC_eL1=589.9451
# remineralization length #2 for POC
bg_par_bio_remin_POC_eL2=1000000.0
# initial fractional abundance of CaCO3 component #2
bg_par_bio_remin_CaCO3_frac2=0.45
# depth of remineralization or CaCO3
bg_par_bio_remin_CaCO3_eL1=1.8905e+003
# remineralization length #2 for CaCO3
bg_par_bio_remin_CaCO3_eL2=1000000.0
#
# *** IRON CYCLING *****************************************************
#
# iron tracer scheme
bg_opt_geochem_Fe='OLD'
# aeolian Fe solubility
bg_par_det_Fe_sol=0.00291468
# modifier of the scavenging rate of dissolved Fe
bg_par_scav_Fe_sf_POC=1.338130
# exponent for aeolian Fe solubility [use 1.0 for uniform solubility]
bg_par_det_Fe_sol_exp=0.500
# no scavenged regeneration
bg_par_scav_fremin=0.0
# return POFe
bg_ctrl_bio_NO_fsedFe=.false.
# Variable Fe:C
bg_ctrl_bio_red_fixedFetoC=.false.
# adjust pK'(FeL)
bg_par_K_FeL_pP=11.0
#(max) C/Fe organic matter ratio
bg_par_bio_red_POFe_POC=250000.0
#[FeT] dependent Fe:C ratio [Ridgwell, 2001] -- power
bg_par_bio_FetoC_pP=-0.4225
#[FeT] dependent Fe:C ratio [Ridgwell, 2001] -- scaling
bg_par_bio_FetoC_K=103684.0
#[FeT] dependent Fe:C ratio [Ridgwell, 2001] -- constant
bg_par_bio_FetoC_C=0.0
#
# *** DATA SAVING ******************************************************
#
# Preformed diagnostics option == BASIC + biology + tracer + proxy
bg_par_data_save_level=11
# force time-slice save at run end only
bg_par_infile_slice_name='save_timeslice_NONE.dat'
# Include experiment name in BIOGEM netCDF filenames? 
bg_ctrl_ncout_expid_name=.false.
#
# *** FORCINGS *********************************************************
#
# specify forcings -- dust field and restoring of atmopsheric pCO2 and d13C
bg_par_forcing_name="p_worjh2.RpCO2_Rp13CO2.Mahowald.0ka"
bg_par_atm_force_scale_val_3=278.0E-06
bg_par_atm_force_scale_val_4=-6.5
#
# *** MISC *************************************************************
#
#
# **********************************************************************
# *** END **************************************************************
# **********************************************************************


