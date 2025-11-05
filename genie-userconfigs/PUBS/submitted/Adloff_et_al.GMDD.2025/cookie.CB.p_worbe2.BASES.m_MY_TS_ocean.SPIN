# *******************************************************************
# *** Mg/Ca PAPER ***************************************************
# *******************************************************************
#
# *** CLIMATE *******************************************************
#
# set climate feedback (climate responding to changing pCO2)
ea_36=y
#
# *** BIOLOGICAL NEW PRODUCTION *************************************
#
# biological scheme ID string
bg_par_bio_prodopt='1N1T_PO4MM'
# maximum rate of conversion of dissolved PO4 into organic matter by phytoplankton (mol kg-1 yr-1)
bg_par_bio_k0_PO4=1.9582242E-06
# [PO4] M-M half-sat value (mol kg-1)
bg_par_bio_c0_PO4=2.1989611E-07
#
# *** ORGANIC MATTER EXPORT RATIOS **********************************
#
# production fraction of dissolved organic matter
bg_par_bio_red_DOMfrac=0.66
#
# *** INORGANIC MATTER EXPORT RATIOS ********************************
#
# fix CaCO3:POC export rain ratio
bg_ctrl_force_CaCO3toPOCrainratio=.true.
bg_par_CaCO3toPOCrainratio_file='CaCO3toPOC.muffin.CB.p_worbe2.BASES.ridgwelletal.SPIN.dat'
#
# *** REMINERALIZATION **********************************************
#
# DOC lifetime (yrs)
bg_par_bio_remin_DOMlifetime=0.5
# initial fractional abundance of POC component #2
bg_par_bio_remin_POC_frac2=6.4591110E-02
# depth of remineralization or particulate organic matter
bg_par_bio_remin_POC_eL1=550.5195
# remineralization length #2 for POC
bg_par_bio_remin_POC_eL2=1000000.0
# initial fractional abundance of CaCO3 component #2
bg_par_bio_remin_CaCO3_frac2=0.4325
# depth of remineralization or CaCO3
bg_par_bio_remin_CaCO3_eL1=1083.486
# remineralization length #2 for CaCO3
bg_par_bio_remin_CaCO3_eL2=1000000.0
#
# *** DATA SAVING ***************************************************
#
# Preformed diagnostics option == OCEAN ACIDIFICATION & FOSSIL FUEL GAMES
bg_par_data_save_level=10
# force time-slice save at run end only
bg_par_infile_slice_name='save_timeslice_NONE.dat'
#
# *** FORCINGS ******************************************************
#
# specify forcings -- generic forcing of atmopsheric pCO2 and d13C
bg_par_forcing_name="pyyyyz.RpCO2_Rp13CO2"
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
# apply MyAmi Mg/Ca correction
gm_par_adj_carbconst_option='Hain2015'
# carbchem constant T and S range -- Mehrbach (DEFAULT)
gm_par_carbchem_Tmin =  2.0
gm_par_carbchem_Tmax = 35.0
gm_par_carbchem_Smin = 26.0
gm_par_carbchem_Smax = 43.0
# modern [Ca], [Mg]
bg_ocn_init_35=1.025E-02
bg_ocn_init_50=5.282E-02
#
# *******************************************************************
# *** END ***********************************************************
# *******************************************************************
#
