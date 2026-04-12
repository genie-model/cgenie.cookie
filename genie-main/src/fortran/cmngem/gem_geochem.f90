! ******************************************************************************************************************************** !
! gem_geochem.f90
! Geochemistry (non carbonate chemistry) Model
! ******************************************************************************************************************************** !


MODULE gem_geochem


  use gem_cmn
  use gem_util
  IMPLICIT NONE
  SAVE


CONTAINS


  ! ****************************************************************************************************************************** !
  ! Fe SPECIATION
  function fun_box_calc_geochem_Fe(dum_FeT,dum_LT,dum_par_K_FeL)
    ! -------------------------------------------------------- !
    ! RESULT VARIABLE
    ! -------------------------------------------------------- !
    real,DIMENSION(1:3)::fun_box_calc_geochem_Fe
    ! -------------------------------------------------------- !
    ! DUMMY ARGUMENTS
    ! -------------------------------------------------------- !
    real,INTENT(in)::dum_FeT,dum_LT
    real,INTENT(in)::dum_par_K_FeL
    ! -------------------------------------------------------- !
    ! DEFINE LOCAL VARIABLES
    ! -------------------------------------------------------- !
    real::loc_Fe,loc_FeL,loc_L
    real,DIMENSION(2)::loc_roots
    ! -------------------------------------------------------- !
    ! INITIALIZE LOCAL VARIABLES
    ! -------------------------------------------------------- ! 
    loc_Fe= 0.0
    loc_FeL = 0.0
    loc_L   = 0.0
    loc_roots(:) = 0.0
    ! -------------------------------------------------------- !
    ! CALCULATE IRON SPECIATION
    ! -------------------------------------------------------- !
    ! -------------------------------------------------------- ! solve Fe speciation equation
    ! K = FeL / (Fe*L) (e.g. see: Parekth et al. [2005])
    ! => FeL = Fe*L*K
    !    conservation relations:
    !    FeL + Fe = FeT => Fe = FeT - FeL
    !    FeL + L  = LT  => L  = LT  - FeL
    !    substitute:
    !    FeL = (FeT - FeL)*(LT - FeL)*K
    !    => FeL/K = FeT*LT + FeL^2 - LT*FeL - FeT*FeL
    !    => FeL/K = FeL^2 - (LT + FeT)*FeL + FeT*LT
    !    => 1.0*FeL^2 - (LT + FeT + 1.0/K)*FeL + FeT*LT = 0.0
    !       solve as: ax2 + bx + c = 0.0
    !                 where x = FeL
    loc_roots(:) = fun_quad_root(1.0,-(dum_LT + dum_FeT + 1.0/dum_par_K_FeL),dum_FeT*dum_LT)
    ! -------------------------------------------------------- ! filter returned roots
    if (maxval(loc_roots(:)) < const_real_nullsmall) then
       IF (ctrl_debug_reportwarnings) THEN
          CALL sub_report_error( &
               & 'biogem_box.f90','sub_calc_geochem_Fe', &
               & 'No REAL root in Fe speciation calculation (or maybe zero ...).'// &
               & ' / Data: loc_FeL(OLD),loc_Fe(OLD),loc_L(OLD),dum_FeT,dum_LT,', &
               & 'SOD THIS FOR A GAME OF SOLDIERS: calculation abondoned ...', &
               & (/loc_FeL,loc_Fe,loc_L,dum_FeT,dum_LT/),.false. &
               & )
          error_stop = .FALSE.
       end IF
    elseif ((minval(loc_roots(:)) > dum_FeT) .AND. (minval(loc_roots(:)) > dum_LT)) then
       IF (ctrl_debug_reportwarnings) THEN
          CALL sub_report_error( &
               & 'biogem_box.f90','sub_calc_geochem_Fe', &
               & 'No solution to Fe speciation calculation possible ... :('// &
               & ' / Data: loc_FeL(OLD),loc_Fe(OLD),loc_L(OLD),dum_FeT,dum_LT,', &
               & 'SOD THIS FOR A GAME OF SOLDIERS: calculation abondoned ...', &
               & (/loc_FeL,loc_Fe,loc_L,dum_FeT,dum_LT/),.false. &
               & )
          error_stop = .FALSE.
       end IF
    else
       if (minval(loc_roots(:)) < const_real_nullsmall) then
          loc_FeL = maxval(loc_roots(:))
       else
          loc_FeL = minval(loc_roots(:))
       end if
       loc_Fe  = dum_FeT - loc_FeL
       loc_L   = dum_LT - loc_FeL
    end if
    ! -------------------------------------------------------- !
    ! RETURN RESULT
    ! -------------------------------------------------------- !
    fun_box_calc_geochem_Fe(1) = loc_Fe
    fun_box_calc_geochem_Fe(2) = loc_FeL
    fun_box_calc_geochem_Fe(3) = loc_L
    ! -------------------------------------------------------- !
    ! END
    ! -------------------------------------------------------- !
  end function fun_box_calc_geochem_Fe
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  ! CALCULATE SOLUBILITY COEFFICIENT
  function fun_calc_solconst(dum_ia,dum_T,dum_S,dum_rho)
    ! -------------------------------------------------------- !
    ! RESULT VARIABLE
    ! -------------------------------------------------------- !
    REAL::fun_calc_solconst
    ! -------------------------------------------------------- !
    ! DUMMY ARGUMENTS
    ! -------------------------------------------------------- !
    integer,INTENT(in)::dum_ia
    real,INTENT(in)::dum_T,dum_S,dum_rho
    ! -------------------------------------------------------- !
    ! DEFINE LOCAL VARIABLES
    ! -------------------------------------------------------- !
    REAL::loc_T,loc_rT,loc_Tr100,loc_S
    ! -------------------------------------------------------- !
    ! INITIALIZE LOCAL VARIABLES
    ! -------------------------------------------------------- ! 
    ! NOTE: pressure in units of (bar) (1 m depth approx = 1 dbar pressure)
    ! NOTE: temperature in K
    ! NOTE: restrict valid T,S range for empirical fit
    ! NOTE: original default was the same as for Mehrbach K1, K2 CONSTANTS:
    !       2.0  <= T <= 35 C
    !       26.0 <= S <= 43 PSU
    ! NOTE: changed from carbchem_* to geochem_*, retaining Mehrbach K1, K2 default range for geochem_*
    if (dum_T <  (const_zeroC +  par_geochem_Tmin))  then
       loc_T = const_zeroC +  par_geochem_Tmin
    elseif (dum_T > (const_zeroC + par_geochem_Tmax)) then
       loc_T = const_zeroC + par_geochem_Tmax
    else
       loc_T = dum_T
    endif
    if (dum_S < par_geochem_Smin) then
       loc_S = par_geochem_Smin
    elseif (dum_S > par_geochem_Smax) then
       loc_S = par_geochem_Smax
    else
       loc_S = dum_S
    endif
    ! -------------------------------------------------------- ! derived local constants
    loc_rT    = 1.0/loc_T
    loc_Tr100 = loc_T/100.0
    ! -------------------------------------------------------- !
    ! CALCULATE Solubility Coefficients (mol/(kg atm)) and return function value
    ! -------------------------------------------------------- !
    ! NOTE: for CO2 and N2O, the soluability coefficient is in units of mol/(kg atm)
    !       rather than as a Bunsen Solubility Coefficient (see Wanninkohf [1992])
    !       => convert units for others
    ! NOTE: for CFC-11 and CFC-12, the soluability coefficient is in units of mol/(kg atm)
    !       rather than as a Bunsen Solubility Coefficient (see Wanninkohf [1992])
    !       (actaully, it is not really this simple and K should be corrected for water vapour pressure and lame things like that)
    SELECT CASE (dum_ia)
    CASE (ia_pCO2,ia_pN2O)
       fun_calc_solconst = EXP( &
            & par_bunsen_coef(1,dum_ia) + par_bunsen_coef(2,dum_ia)*(100*loc_rT) + par_bunsen_coef(3,dum_ia)*LOG(loc_Tr100) + &
            & loc_S* &
            & (par_bunsen_coef(4,dum_ia) + par_bunsen_coef(5,dum_ia)*(loc_Tr100) + par_bunsen_coef(6,dum_ia)*(loc_Tr100)**2) &
            &  )
    CASE (ia_pCFC11,ia_pCFC12)
       fun_calc_solconst = EXP( &
            & par_bunsen_coef(1,dum_ia) + par_bunsen_coef(2,dum_ia)*(100*loc_rT) + par_bunsen_coef(3,dum_ia)*LOG(loc_Tr100) + &
            & loc_S* &
            & (par_bunsen_coef(4,dum_ia) + par_bunsen_coef(5,dum_ia)*(loc_Tr100) + par_bunsen_coef(6,dum_ia)*(loc_Tr100)**2) &
            &  )
    CASE default
       fun_calc_solconst = EXP( &
            & par_bunsen_coef(1,dum_ia) + par_bunsen_coef(2,dum_ia)*(100*loc_rT) + par_bunsen_coef(3,dum_ia)*LOG(loc_Tr100) + &
            & loc_S* &
            & (par_bunsen_coef(4,dum_ia) + par_bunsen_coef(5,dum_ia)*(loc_Tr100) + par_bunsen_coef(6,dum_ia)*(loc_Tr100)**2) &
            &  )/ &
            & (dum_rho*const_V)
    END SELECT
    ! -------------------------------------------------------- !
    ! END
    ! -------------------------------------------------------- !
  end function fun_calc_solconst
  ! ****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  function fun_box_calc_spec_Fe2(dum_Fe2,dum_H2S,dum_par_bio_FeS_abioticohm_cte)
    ! -------------------------------------------------------- !
    ! RESULT VARIABLE
    ! -------------------------------------------------------- !
    real,DIMENSION(1:3)::fun_box_calc_spec_Fe2
    ! -------------------------------------------------------- !
    ! DUMMY ARGUMENTS
    ! -------------------------------------------------------- !
    real,INTENT(in)::dum_Fe2,dum_H2S
    real,INTENT(in)::dum_par_bio_FeS_abioticohm_cte
    ! -------------------------------------------------------- !
    ! DEFINE LOCAL VARIABLES
    ! -------------------------------------------------------- !
    real::loc_Fe2,loc_H2S,loc_FeS
    real,DIMENSION(2)::loc_roots
    ! -------------------------------------------------------- !
    ! INITIALIZE LOCAL VARIABLES
    ! -------------------------------------------------------- ! 
    loc_Fe2 = 0.0
    loc_H2S = 0.0
    loc_FeS = 0.0
    loc_roots(:) = 0.0
    ! -------------------------------------------------------- !
    ! CALCULATE reduced IRON SPECIATION
    ! -------------------------------------------------------- !
    ! -------------------------------------------------------- ! solve reduced Fe speciation equation
    ! 
    ! calculate solution for equilibrium:
    !
    ! Fe2_eq*H2S_eq/FeS_eq = K
    ! Fe2_eq = Fe2_in - x; H2S_eq = H2S_in - x; FeS_eq = FeS_in + x
    ! (Fe2_in - x)(H2S_in - x) = K*FeS_in + K*x
    ! Fe_in*H2S_in - H2S_in*x - Fe_in*x + x^2 = K*FeS_in + K*x
    ! x^2 - (H2S_in + Fe_in + K) + (Fe_in*H2S_in - K*FeS_in) = 0
    !
    ! K = Fe2_eq*H2S_eq/FeS_eq 
    ! => Fe2_eq*H2S_eq = K*FeS_eq
    !    conservation relations:
    !    Fe2_eq + FeS_eq = Fe2_tot => Fe2_eq = Fe2_tot - FeS_eq
    !    H2S_eq + FeS_eq = H2S_tot => H2S_eq = H2S_tot - FeS_eq
    !    substitute:
    !    (Fe2_tot - FeS_eq)(H2S_tot - FeS_eq) = K*FeS_eq
    !    => Fe2_tot*H2S_tot - Fe2_tot*FeS_eq - H2S_tot*FeS_eq + FeS_eq^2 = K*FeS_eq
    !    => 1.0*FeS_eq^2 - (Fe2_tot + H2S_tot + K) + Fe2_tot*H2S_tot = 0.0
    !    solve as: ax2 + bx + c = 0.0
    !                 where x = FeS_eq
    loc_roots(:) = fun_quad_root(1.0,-(dum_Fe2 + dum_H2S + dum_par_bio_FeS_abioticohm_cte),dum_Fe2*dum_H2S)
    ! -------------------------------------------------------- ! filter returned roots
    if (maxval(loc_roots(:)) < const_real_nullsmall) then
       IF (ctrl_debug_reportwarnings) THEN
          CALL sub_report_error( &
               & 'biogem_box.f90','sub_calc_spec_Fe2', &
               & 'No REAL root in Fe2 speciation calculation (or maybe zero ...).'// &
               & ' / Data: loc_Fe2(OLD),loc_H2S(OLD),loc_FeS(OLD),dum_Fe2,dum_H2S,', &
               & 'SOD THIS FOR A GAME OF SOLDIERS: calculation abondoned ...', &
               & (/loc_Fe2,loc_H2S,loc_FeS,dum_Fe2,dum_H2S/),.false. &
               & )
          error_stop = .FALSE.
       end IF
    elseif ((minval(loc_roots(:)) > dum_Fe2) .AND. (minval(loc_roots(:)) > dum_H2S)) then
       IF (ctrl_debug_reportwarnings) THEN
          CALL sub_report_error( &
               & 'biogem_box.f90','sub_calc_spec_Fe2', &
               & 'No solution to Fe2 speciation calculation possible ... :('// &
               & ' / Data: loc_Fe2(OLD),loc_H2S(OLD),loc_FeS(OLD),dum_Fe2,dum_H2S,', &
               & 'SOD THIS FOR A GAME OF SOLDIERS: calculation abondoned ...', &
               & (/loc_Fe2,loc_H2S,loc_FeS,dum_Fe2,dum_H2S/),.false. &
               & )
          error_stop = .FALSE.
       end IF
    else
       if (minval(loc_roots(:)) < const_real_nullsmall) then
          loc_FeS = maxval(loc_roots(:))
       else
          loc_FeS = minval(loc_roots(:))
       end if
       if (loc_FeS > MIN(dum_Fe2,dum_H2S)) then
           loc_FeS = MIN(dum_Fe2,dum_H2S)
       end if
       loc_Fe2  = dum_Fe2 - loc_FeS
       loc_H2S  = dum_H2S - loc_FeS
    end if
    ! -------------------------------------------------------- !
    ! RETURN RESULT
    ! -------------------------------------------------------- !
    fun_box_calc_spec_Fe2(1) = loc_Fe2 
    fun_box_calc_spec_Fe2(2) = loc_H2S
    fun_box_calc_spec_Fe2(3) = loc_FeS
    ! -------------------------------------------------------- !
    ! END
    ! -------------------------------------------------------- !
  end function fun_box_calc_spec_Fe2
  !****************************************************************************************************************************** !


  ! ****************************************************************************************************************************** !
  ! NOTE: this is an initial generic version without isotopc corrections
  ! NOTE: also remove hard threshold and kenetic options
  ! NOTE: required module or global parameters:
  !       conv_ls_lo_c0_O2,    conv_ls_lo_ci_O2,    conv_ls_lo_k_O2
  !       conv_ls_lo_c0_NO3,   conv_ls_lo_ci_NO3,   conv_ls_lo_k_NO3
  !       conv_ls_lo_c0_FeOOH, conv_ls_lo_ci_FeOOH, conv_ls_lo_k_FeOOH
  !       conv_ls_lo_c0_SO4,   conv_ls_lo_ci_SO4,   conv_ls_lo_k_SO4
  !                                                 conv_ls_lo_k_meth
  function fun_conv_ls_lo_remin(dum_ocn)
    ! -------------------------------------------------------- !
    ! RESULT VARIABLE
    ! -------------------------------------------------------- !
    real,dimension(1:n_l_ocn,1:n_l_sed)::fun_conv_ls_lo_remin 
    ! -------------------------------------------------------- !
    ! DUMMY ARGUMENTS
    ! -------------------------------------------------------- !
    real,dimension(1:n_l_ocn),INTENT(in)::dum_ocn              !
    ! -------------------------------------------------------- !
    ! DEFINE LOCAL VARIABLES
    ! -------------------------------------------------------- !
    real::loc_k
    real::loc_O2,loc_NO3,loc_FeOOH,loc_SO4
    real::loc_kO2,loc_kNO3,loc_kFeOOH,loc_kSO4,loc_kmeth
    real::loc_kiO2,loc_kiNO3,loc_kiFeOOH,loc_kiSO4
    real::loc_r15N,loc_r56Fe,loc_r34S
    real::loc_alpha,loc_R
    real::loc_ls_lo_epsilon_NO3,loc_ls_lo_epsilon_SO4,loc_ls_lo_epsilon_FeOOH
    real,dimension(1:n_l_ocn,1:n_l_sed)::loc_conv_ls_lo 
    ! ---------------------------------------------------------- ! initialize local variables
    loc_k = 0.0
    loc_ls_lo_epsilon_NO3   = 0.0
    loc_ls_lo_epsilon_SO4   = 0.0
    loc_ls_lo_epsilon_FeOOH = 0.0
    ! ---------------------------------------------------------- !
    ! CREATE REMIN ARRAY
    ! ---------------------------------------------------------- !
    ! ---------------------------------------------------------- ! set MM-type rate limitations
    ! NOTE: equation form follows Arndt et al. [2013] (ESR) and Boudreau [1997] (book)
    ! NOTE: truncate loc concentrations at zero to avoid negative values being propagated ...
    ! NOTE: catch a local oxidation value of const_real_nullsmall (or less),
    !       as this means that the oxidant was probably negative in the first place
    !       => disable that particular redox remin pathway by setting the kinetic parameter to ZERO
    ! NOTE: to disable inhibition, set the inhibition parameter to a very large number, then
    !       e.g. conv_ls_lo_ci_O2/(conv_ls_lo_ci_O2 + loc_O2) approaches a value of 1.0
    if (ocn_select(io_O2)) then
       loc_O2 = dum_ocn(io2l(io_O2))
       if (loc_O2 <= const_real_nullsmall) then
          loc_O2   = 0.0
          loc_kO2  = 0.0
          loc_kiO2 = 1.0
       else
          loc_kO2 = loc_O2/(loc_O2 + conv_ls_lo_c0_O2)
          loc_kiO2 = conv_ls_lo_ci_O2/(conv_ls_lo_ci_O2 + loc_O2)
       end if
       loc_k    = loc_k + conv_ls_lo_k_O2*loc_kO2
    else
       loc_O2   = 0.0
       loc_kO2  = 0.0
       loc_kiO2 = 1.0
    end if
    if (ocn_select(io_NO3)) then
       loc_NO3 = dum_ocn(io2l(io_NO3))
       if (loc_NO3 <= const_real_nullsmall) then
          loc_NO3   = 0.0
          loc_kNO3  = 0.0
          loc_kiNO3 = 1.0
       else
          loc_kNO3 = loc_NO3/(loc_NO3 + conv_ls_lo_c0_NO3)
          loc_kiNO3 = conv_ls_lo_ci_NO3/(conv_ls_lo_ci_NO3 + loc_NO3)
       end if
       loc_k     = loc_k + conv_ls_lo_k_NO3*loc_kNO3*loc_kiO2
    else
       loc_NO3   = 0.0
       loc_kNO3  = 0.0
       loc_kiNO3 = 1.0
    end if
    if (ocn_select(io_FeOOH)) then
       loc_FeOOH = dum_ocn(io2l(io_FeOOH))
       if (loc_FeOOH <= const_real_nullsmall) then
          loc_FeOOH   = 0.0
          loc_kFeOOH  = 0.0
          loc_kiFeOOH = 1.0
       else
          loc_kFeOOH = loc_FeOOH/(loc_FeOOH + conv_ls_lo_c0_FeOOH)
          loc_kiFeOOH = conv_ls_lo_ci_FeOOH/(conv_ls_lo_ci_FeOOH + loc_FeOOH)
       end if
       loc_k     = loc_k + conv_ls_lo_k_FeOOH*loc_kFeOOH*loc_kiNO3*loc_kiO2
    else
       loc_FeOOH   = 0.0
       loc_kFeOOH  = 0.0
       loc_kiFeOOH = 1.0
    end if
    if (ocn_select(io_SO4)) then
       loc_SO4 = dum_ocn(io2l(io_SO4))
       if (loc_SO4 <= const_real_nullsmall) then
          loc_SO4   = 0.0
          loc_kSO4  = 0.0
          loc_kiSO4 = 1.0
       else
          loc_kSO4  = loc_SO4/(loc_SO4 + conv_ls_lo_c0_SO4)
          loc_kiSO4 = conv_ls_lo_ci_SO4/(conv_ls_lo_ci_SO4 + loc_SO4)
       end if
       loc_k     = loc_k + conv_ls_lo_k_SO4*loc_kSO4*loc_kiFeOOH*loc_kiNO3*loc_kiO2
    else
       loc_SO4   = 0.0
       loc_kSO4  = 0.0
       loc_kiSO4 = 1.0
    end if
    if (ocn_select(io_CH4)) then
       loc_kmeth = 1.0
       loc_k     = loc_k + conv_ls_lo_k_meth*loc_kmeth*loc_kiSO4*loc_kiFeOOH*loc_kiNO3*loc_kiO2
    else
       loc_kmeth = 0.0
    end if
    ! ---------------------------------------------------------- ! check *some* remin occurs
    ! NOTE: parameters adjusted in factor (conv_ls_lo_k_O2*loc_kO2/loc_k) to make this unity
    !       no remin would otherwise occur
    if (loc_k < const_real_nullsmall) then
       loc_kO2 = 1.0
       loc_k   = 1.0/conv_ls_lo_k_O2
    end if
    ! ---------------------------------------------------------- ! calculate weighted remin array
    ! NOTE: normalize to 1.0
    if (ocn_select(io_O2)) then
       if (ocn_select(io_NO3)) then
          if (ocn_select(io_FeOOH)) then
             if (ocn_select(io_SO4)) then
                if (ocn_select(io_CH4)) then
                   ! O2 + NO3 + FeOOH + SO4 + CH4
                   loc_conv_ls_lo(:,:) = &
                        & (conv_ls_lo_k_O2*loc_kO2/loc_k)*conv_ls_lo_O(:,:) + &
                        & (conv_ls_lo_k_NO3*loc_kNO3*loc_kiO2/loc_k)*conv_ls_lo_N(:,:) + &
                        & (conv_ls_lo_k_FeOOH*loc_kFeOOH*loc_kiNO3*loc_kiO2/loc_k)*conv_ls_lo_Fe(:,:) + &
                        & (conv_ls_lo_k_SO4*loc_kSO4*loc_kiFeOOH*loc_kiNO3*loc_kiO2/loc_k)*conv_ls_lo_S(:,:) + &
                        & (conv_ls_lo_k_meth*loc_kmeth*loc_kiSO4*loc_kiFeOOH*loc_kiNO3*loc_kiO2/loc_k)*conv_ls_lo_meth(:,:)
                else
                   ! O2 + NO3 + FeOOH + SO4
                   loc_conv_ls_lo(:,:) = &
                        & (conv_ls_lo_k_O2*loc_kO2/loc_k)*conv_ls_lo_O(:,:) + &
                        & (conv_ls_lo_k_NO3*loc_kNO3*loc_kiO2/loc_k)*conv_ls_lo_N(:,:) + &
                        & (conv_ls_lo_k_FeOOH*loc_kFeOOH*loc_kiNO3*loc_kiO2/loc_k)*conv_ls_lo_Fe(:,:) + &
                        & (conv_ls_lo_k_SO4*loc_kSO4*loc_kiFeOOH*loc_kiNO3*loc_kiO2/loc_k)*conv_ls_lo_S(:,:)
                end if
             else
                ! O2 + NO3 + FeOOH
                loc_conv_ls_lo(:,:) = &
                     & (conv_ls_lo_k_O2*loc_kO2/loc_k)*conv_ls_lo_O(:,:) + &
                     & (conv_ls_lo_k_NO3*loc_kNO3*loc_kiO2/loc_k)*conv_ls_lo_N(:,:) + &
                     & (conv_ls_lo_k_FeOOH*loc_kFeOOH*loc_kiNO3*loc_kiO2/loc_k)*conv_ls_lo_Fe(:,:)
             end if
          else
             if (ocn_select(io_SO4)) then
                if (ocn_select(io_CH4)) then
                   ! O2 + NO3 + SO4 + CH4
                   loc_conv_ls_lo(:,:) = &
                        & (conv_ls_lo_k_O2*loc_kO2/loc_k)*conv_ls_lo_O(:,:) + &
                        & (conv_ls_lo_k_NO3*loc_kNO3*loc_kiO2/loc_k)*conv_ls_lo_N(:,:) + &
                        & (conv_ls_lo_k_SO4*loc_kSO4*loc_kiNO3*loc_kiO2/loc_k)*conv_ls_lo_S(:,:) + &
                        & (conv_ls_lo_k_meth*loc_kmeth*loc_kiSO4*loc_kiNO3*loc_kiO2/loc_k)*conv_ls_lo_meth(:,:)
                else
                   ! O2 + NO3 + SO4
                   loc_conv_ls_lo(:,:) = &
                        & (conv_ls_lo_k_O2*loc_kO2/loc_k)*conv_ls_lo_O(:,:) + &
                        & (conv_ls_lo_k_NO3*loc_kNO3*loc_kiO2/loc_k)*conv_ls_lo_N(:,:) + &
                        & (conv_ls_lo_k_SO4*loc_kSO4*loc_kiNO3*loc_kiO2/loc_k)*conv_ls_lo_S(:,:)
                end if
             else
                ! O2 + NO3
                loc_conv_ls_lo(:,:) = &
                     & (conv_ls_lo_k_O2*loc_kO2/loc_k)*conv_ls_lo_O(:,:) + &
                     & (conv_ls_lo_k_NO3*loc_kNO3*loc_kiO2/loc_k)*conv_ls_lo_N(:,:)
             end if
          endif
       elseif (ocn_select(io_FeOOH)) then
          if (ocn_select(io_SO4)) then
             if (ocn_select(io_CH4)) then
                ! O2 + FeOOH + SO4 + CH4
                loc_conv_ls_lo(:,:) = &
                     & (conv_ls_lo_k_O2*loc_kO2/loc_k)*conv_ls_lo_O(:,:) + &
                     & (conv_ls_lo_k_FeOOH*loc_kFeOOH*loc_kiO2/loc_k)*conv_ls_lo_Fe(:,:) + &
                     & (conv_ls_lo_k_SO4*loc_kSO4*loc_kiFeOOH*loc_kiO2/loc_k)*conv_ls_lo_S(:,:) + &
                     & (conv_ls_lo_k_meth*loc_kmeth*loc_kiSO4*loc_kiFeOOH*loc_kiO2/loc_k)*conv_ls_lo_meth(:,:)
             else
                ! O2 + FeOOH + SO4
                loc_conv_ls_lo(:,:) = &
                     & (conv_ls_lo_k_O2*loc_kO2/loc_k)*conv_ls_lo_O(:,:) + &
                     & (conv_ls_lo_k_FeOOH*loc_kFeOOH*loc_kiO2/loc_k)*conv_ls_lo_Fe(:,:) + &
                     & (conv_ls_lo_k_SO4*loc_kSO4*loc_kiFeOOH*loc_kiO2/loc_k)*conv_ls_lo_S(:,:)
             end if
          else
             ! O2 + FeOOH
             loc_conv_ls_lo(:,:) = &
                  & (conv_ls_lo_k_O2*loc_kO2/loc_k)*conv_ls_lo_O(:,:) + &
                  & (conv_ls_lo_k_FeOOH*loc_kFeOOH*loc_kiO2/loc_k)*conv_ls_lo_Fe(:,:)
          end if
       elseif (ocn_select(io_SO4)) then
          if (ocn_select(io_CH4)) then
             ! O2 + SO4 + CH4
             loc_conv_ls_lo(:,:) = &
                  & (conv_ls_lo_k_O2*loc_kO2/loc_k)*conv_ls_lo_O(:,:) + &
                  & (conv_ls_lo_k_SO4*loc_kSO4*loc_kiO2/loc_k)*conv_ls_lo_S(:,:) + &
                  & (conv_ls_lo_k_meth*loc_kmeth*loc_kiSO4*loc_kiO2/loc_k)*conv_ls_lo_meth(:,:)
          else
             ! O2 + SO4
             loc_conv_ls_lo(:,:) = &
                  & (conv_ls_lo_k_O2*loc_kO2/loc_k)*conv_ls_lo_O(:,:) + &
                  & (conv_ls_lo_k_SO4*loc_kSO4*loc_kiO2/loc_k)*conv_ls_lo_S(:,:)
          end if
       else
          if (ocn_select(io_CH4)) then
             ! O2 + CH4
             loc_conv_ls_lo(:,:) = &
                  & (conv_ls_lo_k_O2*loc_kO2/loc_k)*conv_ls_lo_O(:,:) + &
                  & (conv_ls_lo_k_meth*loc_kmeth*loc_kiO2/loc_k)*conv_ls_lo_meth(:,:)
          else
             ! O2
             loc_conv_ls_lo(:,:) = &
                  & (conv_ls_lo_k_O2*loc_kO2/loc_k)*conv_ls_lo_O(:,:)
          end if
       end if
    else
       loc_conv_ls_lo(:,:) = 0.0
    end if
    ! ---------------------------------------------------------- !
    ! CORRECT FOR ISOTOPES (in oxidants -- NO3, SO4, FeOOH)
    ! ---------------------------------------------------------- !
    ! NOTE: use fractionation factors: loc_ls_lo_epsilon_NO3, loc_ls_lo_epsilon_SO4, loc_ls_lo_epsilon_FeOOH
    !       for now these are local and zero ... but this need not be the case (and could be propegated from BIOGEM)
    !       (we could employ an array to store fractionation factor for all ocean tracers ... or define an oxidants-only array) 
    ! NOTE: for NO3, possible products fomr NO3 reduction are NO2, N2O, N2, NH4 (hence the 4 logical queries)
    ! test for N isotopes selected
    ! NOTE: value of loc_NO3 already determined
    if (ocn_select(io_NO3) .AND. ocn_select(io_NO3_15N)) then
       if (loc_NO3 > const_real_nullsmall) then
          loc_r15N  = dum_ocn(io2l(io_NO3_15N))/loc_NO3
       else
          loc_r15N  = 0.0
       end if
       loc_alpha = 1.0 + loc_ls_lo_epsilon_NO3/1000.0
       loc_R     = loc_r15N/(1.0 - loc_r15N)
       loc_conv_ls_lo(io2l(io_NO3_15N),:) = loc_alpha*loc_R/(1.0 + loc_alpha*loc_R)*loc_conv_ls_lo(io2l(io_NO3),:)
       if (ocn_select(io_NO2) .AND. ocn_select(io_NO2_15N)) then
          loc_conv_ls_lo(io2l(io_NO2_15N),:) = loc_alpha*loc_R/(1.0 + loc_alpha*loc_R)*loc_conv_ls_lo(io2l(io_NO2),:)
       elseif (ocn_select(io_N2O) .AND. ocn_select(io_N2O_15N)) then
          loc_conv_ls_lo(io2l(io_N2O_15N),:) = loc_alpha*loc_R/(1.0 + loc_alpha*loc_R)*loc_conv_ls_lo(io2l(io_N2O),:)
       elseif (ocn_select(io_N2) .AND. ocn_select(io_N2_15N)) then
          loc_conv_ls_lo(io2l(io_N2_15N),:)  = loc_alpha*loc_R/(1.0 + loc_alpha*loc_R)*loc_conv_ls_lo(io2l(io_N2),:)
       elseif (ocn_select(io_NH4) .AND. ocn_select(io_NH4_15N)) then
          loc_conv_ls_lo(io2l(io_NH4_15N),:) = loc_alpha*loc_R/(1.0 + loc_alpha*loc_R)*loc_conv_ls_lo(io2l(io_NH4),:)
       else
          ! [DEFAULT, oxic remin relationship]
       endif
    end if
    ! test for S isotopes selected
    ! NOTE: value of loc_SO4 already determined
    if (ocn_select(io_SO4_34S) .AND. ocn_select(io_H2S_34S)) then
       if (loc_SO4 > const_real_nullsmall) then
          loc_r34S  = dum_ocn(io2l(io_SO4_34S)) / loc_SO4
       else
          loc_r34S  = 0.0
       end if
       loc_alpha = 1.0 + loc_ls_lo_epsilon_SO4/1000.0
       loc_R     = loc_r34S/(1.0 - loc_r34S)
       loc_conv_ls_lo(io2l(io_SO4_34S),:) = loc_alpha*loc_R/(1.0 + loc_alpha*loc_R)*loc_conv_ls_lo(io2l(io_SO4),:)
       loc_conv_ls_lo(io2l(io_H2S_34S),:) = loc_alpha*loc_R/(1.0 + loc_alpha*loc_R)*loc_conv_ls_lo(io2l(io_H2S),:)
    end if
    ! test for Fe isotopes selected
    ! NOTE: value of loc_FeOOH already determined
    if (ocn_select(io_FeOOH_56Fe) .AND. ocn_select(io_Fe2_56Fe)) then
       if (loc_FeOOH > const_real_nullsmall) then
          loc_r56Fe  = dum_ocn(io2l(io_FeOOH_56Fe)) / loc_FeOOH
       else
          loc_r56Fe  = 0.0
       end if
       loc_alpha = 1.0 + loc_ls_lo_epsilon_FeOOH/1000.0
       loc_R     = loc_r56Fe/(1.0 - loc_r56Fe)
       loc_conv_ls_lo(io2l(io_FeOOH_56Fe),:) = loc_alpha*loc_R/(1.0 + loc_alpha*loc_R)*loc_conv_ls_lo(io2l(io_FeOOH_56Fe),:)
       loc_conv_ls_lo(io2l(io_Fe2_56Fe),:)   = loc_alpha*loc_R/(1.0 + loc_alpha*loc_R)*loc_conv_ls_lo(io2l(io_Fe2_56Fe),:)
    end if
    ! -------------------------------------------------------- !
    ! RETURN RESULT
    ! -------------------------------------------------------- !
    fun_conv_ls_lo_remin = loc_conv_ls_lo 
    ! -------------------------------------------------------- !
    ! END
    ! -------------------------------------------------------- !
  end function fun_conv_ls_lo_remin
  ! ****************************************************************************************************************************** !

END MODULE gem_geochem

