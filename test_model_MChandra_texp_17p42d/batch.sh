#!/bin/tcsh

#**********************************************************************
#   Define directories (CAPITAL LETTERS)
#**********************************************************************

# ATOMIC is atomic data directory

setenv OMP_NUM_THREADS 8
setenv ATOMIC       ~jdh/ATOMIC
setenv CMFGEN_PROG  ~kdw25/cur_cmf_gammaray_v2/exe/cmfgen_dev.exe

#***********************************************************************
#    Atomic data soft links. All other shell commands should be
#    place after these links. These links are model dependent. Ther
#    are required by CMFGEN, CMF_FLUX, and DISPGEN
#
#    When N_F=N_S, no F_TO_S link is needed.
#***********************************************************************

# Generic
# -------
 #*****************************************************************************

 ln -sf $ATOMIC/misc/two_phot_data.dat                     TWO_PHOT_DATA
 ln -sf $ATOMIC/misc/xray_phot_fits.dat                    XRAY_PHOT_FITS
 ln -sf $ATOMIC/misc/rs_xray_fluxes_sol.dat                RS_XRAY_FLUXES
 ln -sf $ATOMIC/misc/chg_exch_data.7dec00                  CHG_EXCH_DATA
 ln -sf $ATOMIC/HYD/I/5dec96/hyd_l_data.dat                HYD_L_DATA
 ln -sf $ATOMIC/HYD/I/5dec96/gbf_n_data.dat                GBF_N_DATA

#*****************************************************************************
#    Hydrogen
#   ----------
#*****************************************************************************

 ln -sf  $ATOMIC/HYD/I/5dec96/hiphot.dat                   PHOTHI_A
 ln -sf  $ATOMIC/HYD/I/5dec96/hi_osc.dat                   HI_F_OSCDAT
 ln -sf  $ATOMIC/HYD/I/5dec96/hi_f_to_s_15.dat             HI_F_TO_S
 ln -sf  $ATOMIC/HYD/I/5dec96/hicol.dat                    HI_COL_DATA

#*****************************************************************************
#    Helium
#   --------
#*****************************************************************************

 ln -sf  $ATOMIC/HE/II/5dec96/he2phot.dat                  PHOTHe2_A
 ln -sf  $ATOMIC/HE/II/5dec96/he2_osc.dat                  He2_F_OSCDAT
 ln -sf  $ATOMIC/HE/II/5dec96/he2_f_to_s_10.dat            He2_F_TO_S
 ln -sf  $ATOMIC/HE/II/5dec96/he2col.dat                   He2_COL_DATA

 ln -sf $ATOMIC/HE/I/11may07/heiphot_a7.dat                PHOTHeI_A
 ln -sf $ATOMIC/HE/I/11may07/heioscdat_a7.dat              HeI_F_OSCDAT
 ln -sf $ATOMIC/HE/I/11may07/hei_f_to_s_a7.dat             HeI_F_TO_S
 ln -sf $ATOMIC/HE/I/11may07/heicol.dat                    HeI_COL_DATA

#*****************************************************************************
#    Carbon
#   --------
#*****************************************************************************

 ln -sf  $ATOMIC/CARB/I/12dec04/phot_smooth_50                  PHOTCI_A
 ln -sf  $ATOMIC/CARB/I/12dec04/ci_split_osc                    CI_F_OSCDAT
 ln -sf  $ATOMIC/CARB/I/12dec04/cicol.dat                       CI_COL_DATA
 ln -sf  $ATOMIC/CARB/I/12dec04/f_to_s_122.dat                  CI_F_TO_S

 ln -sf  $ATOMIC/CARB/II/30oct12/phot_sm_0.dat                  PHOTC2_A
 ln -sf  $ATOMIC/CARB/II/30oct12/phot_data_B.dat                PHOTC2_B
 ln -sf  $ATOMIC/CARB/II/30oct12/c2osc_rev.dat                  C2_F_OSCDAT
 ln -sf  $ATOMIC/CARB/II/30oct12/f_to_s_104.dat                 C2_F_TO_S
 ln -sf  $ATOMIC/CARB/II/30oct12/c2col.dat                      C2_COL_DATA
 ln -sf  $ATOMIC/CARB/II/30oct12/c2_auto.dat                    C2_AUTO_DATA

 ln -sf $ATOMIC/CARB/III/23dec04/ciiiphot_sm_a_500.dat          PHOTCIII_A
 ln -sf $ATOMIC/CARB/III/23dec04/ciiiphot_sm_b_500.dat          PHOTCIII_B
 ln -sf $ATOMIC/CARB/III/23dec04/ciii_f_to_s_split_big.dat      CIII_F_TO_S
 ln -sf $ATOMIC/CARB/III/23dec04/ciiiosc_st_split_big.dat       CIII_F_OSCDAT
 ln -sf $ATOMIC/CARB/III/23dec04/ciiicol.dat                    CIII_COL_DATA
 ln -sf $ATOMIC/CARB/III/23dec04/dieciii_ic.dat                 DIECIII

 ln -sf  $ATOMIC/CARB/IV/30oct12/civphot_a12.dat                PHOTCIV_A
 ln -sf  $ATOMIC/CARB/IV/30oct12/civosc_a12_split.dat           CIV_F_OSCDAT
 ln -sf  $ATOMIC/CARB/IV/30oct12/civ_f_to_s_split.dat           CIV_F_TO_S
 ln -sf  $ATOMIC/CARB/IV/30oct12/civcol.dat                     CIV_COL_DATA

 ln -sf  $ATOMIC/CARB/V/20oct02/phot_smooth_3000                PHOTCV_A
 ln -sf  $ATOMIC/CARB/V/20oct02/fin_osc                         CV_F_OSCDAT
 ln -sf  $ATOMIC/CARB/V/20oct02/f_to_s_nosplit_43               CV_F_TO_S
 ln -sf  $ATOMIC/CARB/V/20oct02/col_guess.dat                   CV_COL_DATA

 ln -sf  $ATOMIC/CARB/VI/9may02/phot_smooth_3000                PHOTCSIX_A
 ln -sf  $ATOMIC/CARB/VI/9may02/fin_osc                         CSIX_F_OSCDAT
 ln -sf  $ATOMIC/CARB/VI/9may02/f_to_s_30.dat                   CSIX_F_TO_S
 ln -sf  $ATOMIC/CARB/VI/9may02/col_guess.dat                   CSIX_COL_DATA

#*****************************************************************************
#   Nitrogen
#   --------
#*****************************************************************************

 ln -sf $ATOMIC/NIT/I/12sep12/niphot_a.dat                 PHOTNI_A
 ln -sf $ATOMIC/NIT/I/12sep12/niphot_b.dat                 PHOTNI_B
 ln -sf $ATOMIC/NIT/I/12sep12/niphot_c.dat                 PHOTNI_C
 ln -sf $ATOMIC/NIT/I/12sep12/ni_osc                       NI_F_OSCDAT
 ln -sf $ATOMIC/NIT/I/12sep12/f_to_s_term.dat              NI_F_TO_S

 ln -sf  $ATOMIC/NIT/II/23jan06/fin_osc                    N2_F_OSCDAT
 ln -sf  $ATOMIC/NIT/II/23jan06/phot_nosm_A                PHOTN2_A
 ln -sf  $ATOMIC/NIT/II/23jan06/phot_sm_B                  PHOTN2_B
 ln -sf  $ATOMIC/NIT/II/23jan06/n2col.dat                  N2_COL_DATA
 ln -sf  $ATOMIC/NIT/II/23jan06/f_to_s_157                 N2_F_TO_S
 ln -sf  $ATOMIC/NIT/II/23jan06/n2auto.dat                 N2_AUTO_DATA

 ln -sf  $ATOMIC/NIT/III/24mar07/phot_sm_0_A.dat           PHOTNIII_A
 ln -sf  $ATOMIC/NIT/III/24mar07/phot_sm_0_B.dat           PHOTNIII_B
 ln -sf  $ATOMIC/NIT/III/24mar07/niiiosc_rev.dat           NIII_F_OSCDAT
 ln -sf  $ATOMIC/NIT/III/24mar07/f_to_s_156                NIII_F_TO_S
 ln -sf  $ATOMIC/NIT/III/24mar07/niii_col_data_new.dat     NIII_COL_DATA

 ln -sf  $ATOMIC/NIT/IV/4nov10/PHOTNIV_A_op                PHOTNIV_A
 ln -sf  $ATOMIC/NIT/IV/4nov10/PHOTNIV_B                   PHOTNIV_B
 ln -sf  $ATOMIC/NIT/IV/4nov10/nivosc                      NIV_F_OSCDAT
 ln -sf  $ATOMIC/NIT/IV/4nov10/f_to_s_139                  NIV_F_TO_S
 ln -sf  $ATOMIC/NIT/IV/4nov10/nivcol                      NIV_COL_DATA
 ln -sf  $ATOMIC/NIT/IV/4nov10/nivdie                      DIENIV

 ln -sf  $ATOMIC/NIT/V/30oct12/nvphot_a12.dat              PHOTNV_A
 ln -sf  $ATOMIC/NIT/V/30oct12/nvosc_a12_split.dat         NV_F_OSCDAT
 ln -sf  $ATOMIC/NIT/V/30oct12/nv_f_to_s_split_sm.dat      NV_F_TO_S
 ln -sf  $ATOMIC/NIT/V/30oct12/nvcol.dat                   NV_COL_DATA

#*****************************************************************************
#    Oxygen
#   --------
#*****************************************************************************

 ln -sf $ATOMIC/OXY/I/20sep11/oi_osc_mchf                         OI_F_OSCDAT
 ln -sf $ATOMIC/OXY/I/20sep11/f_to_s_ls                           OI_F_TO_S
 ln -sf $ATOMIC/OXY/I/20sep11/phot_nosm_A                         PHOTOI_A
 ln -sf $ATOMIC/OXY/I/20sep11/phot_nosm_B                         PHOTOI_B
 ln -sf $ATOMIC/OXY/I/20sep11/oi_col                              OI_COL_DATA

 ln -sf $ATOMIC/OXY/II/23mar05/phot_sm_3000.dat                   PHOTO2_A
 ln -sf $ATOMIC/OXY/II/23mar05/o2osc_fin.dat                      O2_F_OSCDAT
 ln -sf $ATOMIC/OXY/II/23mar05/f_to_s_ls.dat                      O2_F_TO_S
 ln -sf $ATOMIC/OXY/II/23mar05/o2col.dat                          O2_COL_DATA

 ln -sf  $ATOMIC/OXY/III/15mar08/phot_sm_0                        PHOTOIII_A
 ln -sf  $ATOMIC/OXY/III/15mar08/oiiiosc                          OIII_F_OSCDAT
 ln -sf  $ATOMIC/OXY/III/15mar08/f_to_s_165                       OIII_F_TO_S
 ln -sf  $ATOMIC/OXY/III/15mar08/col_data_oiii_butler_2012.dat    OIII_COL_DATA
 ln -sf  $ATOMIC/OXY/III/5dec96/oiiidie.dat                       DIEOIII

 ln -sf  $ATOMIC/OXY/IV/28dec02/phot_sm_3000_A                    PHOTOIV_A
 ln -sf  $ATOMIC/OXY/IV/28dec02/phot_sm_3000_B                    PHOTOIV_B
 ln -sf  $ATOMIC/OXY/IV/28dec02/fin_osc                           OIV_F_OSCDAT
 ln -sf  $ATOMIC/OXY/IV/28dec02/f_to_s_ls_242                     OIV_F_TO_S
 ln -sf  $ATOMIC/OXY/IV/28dec02/oivcol.dat                        OIV_COL_DATA
 ln -sf  $ATOMIC/OXY/IV/5dec96/oivdoubdie_103.dat                 DIEOIV

 ln -sf  $ATOMIC/OXY/V/5dec96/ovphot_a.dat                        PHOTOV_A
 ln -sf  $ATOMIC/OXY/V/5dec96/ovphot_b.dat                        PHOTOV_B
 ln -sf  $ATOMIC/OXY/V/5dec96/ovosc_ns_split.dat                  OV_F_OSCDAT
 ln -sf  $ATOMIC/OXY/V/5dec96/ov_f_to_s_split_ext.dat             OV_F_TO_S
 ln -sf  $ATOMIC/OXY/V/5dec96/ovcol.dat                           OV_COL_DATA
 ln -sf  $ATOMIC/OXY/V/5dec96/ovdie.dat                           DIEOV

 ln -sf  $ATOMIC/OXY/VI/30oct12/osixphot_a12.dat                  PHOTOSIX_A
 ln -sf  $ATOMIC/OXY/VI/30oct12/osixosc_a12_split.dat             OSIX_F_OSCDAT
 ln -sf  $ATOMIC/OXY/VI/30oct12/osix_f_to_s_split.dat             OSIX_F_TO_S
 ln -sf  $ATOMIC/OXY/VI/30oct12/osixcol.dat                       OSIX_COL_DATA

#*****************************************************************************
# Neon
# ----
#*****************************************************************************

 ln -sf  $ATOMIC/NEON/I/9sep11/fin_phot                PHOTNeI_A
 ln -sf  $ATOMIC/NEON/I/9sep11/fin_osc                 NeI_F_OSCDAT
 ln -sf  $ATOMIC/NEON/I/9sep11/f_to_s_ls               NeI_F_TO_S
 ln -sf  $ATOMIC/NEON/I/9sep11/col_guess               NeI_COL_DATA

 ln -sf  $ATOMIC/NEON/II/19nov07/phot_nosm                 PHOTNe2_A
 ln -sf  $ATOMIC/NEON/II/19nov07/fin_osc                   Ne2_F_OSCDAT
 ln -sf  $ATOMIC/NEON/II/19nov07/f_to_s_42                 Ne2_F_TO_S
 ln -sf  $ATOMIC/NEON/II/19nov07/col_neii                  Ne2_COL_DATA

 ln -sf  $ATOMIC/NEON/III/19nov07/phot_nosm                PHOTNeIII_A
 ln -sf  $ATOMIC/NEON/III/19nov07/fin_osc                  NeIII_F_OSCDAT
 ln -sf  $ATOMIC/NEON/III/19nov07/f_to_s_57                NeIII_F_TO_S
 ln -sf  $ATOMIC/NEON/III/19nov07/col_neiii               NeIII_COL_DATA

 ln -sf  $ATOMIC/NEON/IV/1dec99/phot_sm_3000.dat           PHOTNeIV_A
 ln -sf  $ATOMIC/NEON/IV/1dec99/fin_osc.dat                NeIV_F_OSCDAT
 ln -sf  $ATOMIC/NEON/IV/1dec99/f_to_s_53.dat              NeIV_F_TO_S
 ln -sf  $ATOMIC/NEON/IV/1dec99/col_data.dat               NeIV_COL_DATA

 ln -sf  $ATOMIC/NEON/V/20jun01/phot_sm_3000.dat             PHOTNeV_A
 ln -sf  $ATOMIC/NEON/V/20jun01/nevosc_rev.dat               NeV_F_OSCDAT
 ln -sf  $ATOMIC/NEON/V/20jun01/f_to_s_60.dat                NeV_F_TO_S
 ln -sf  $ATOMIC/NEON/V/20jun01/col_data.dat                 NeV_COL_DATA

 ln -sf  $ATOMIC/NEON/VI/20jun01/phot_sm_3000.dat           PHOTNeSIX_A
 ln -sf  $ATOMIC/NEON/VI/20jun01/neviosc_rev.dat            NeSIX_F_OSCDAT
 ln -sf  $ATOMIC/NEON/VI/20jun01/f_to_s_43.dat              NeSIX_F_TO_S
 ln -sf  $ATOMIC/NEON/VI/20jun01/col_data.dat               NeSIX_COL_DATA

 ln -sf  $ATOMIC/NEON/VII/20jun01/phot_sm_3000.dat          PHOTNeSEV_A
 ln -sf  $ATOMIC/NEON/VII/20jun01/neviiosc_rev.dat          NeSEV_F_OSCDAT
 ln -sf  $ATOMIC/NEON/VII/20jun01/f_to_s_46.dat             NeSEV_F_TO_S
 ln -sf  $ATOMIC/NEON/VII/20jun01/col_data.dat              NeSEV_COL_DATA

 ln -sf  $ATOMIC/NEON/VIII/20jun01/phot_sm_3000.dat         PHOTNeVIII_A
 ln -sf  $ATOMIC/NEON/VIII/20jun01/neviiiosc_rev.dat        NeVIII_F_OSCDAT
 ln -sf  $ATOMIC/NEON/VIII/20jun01/f_to_s_29.dat            NeVIII_F_TO_S
 ln -sf  $ATOMIC/NEON/VIII/20jun01/col_guess.dat            NeVIII_COL_DATA

#*****************************************************************************
#
# Sodium
#--------
#*****************************************************************************

 ln -sf $ATOMIC/NA/I/5aug97/nai_osc_split.dat             NaI_F_OSCDAT
 ln -sf $ATOMIC/NA/I/5aug97/nai_phot_a.dat                PHOTNaI_A
 ln -sf $ATOMIC/NA/I/5aug97/col_guess.dat                 NaI_COL_DATA
 ln -sf $ATOMIC/NA/I/5aug97/nai_f_to_s_sm.dat             NaI_F_TO_S

 ln -sf  $ATOMIC/NA/II/15feb01/phot_data.dat              PHOTNa2_A
 ln -sf  $ATOMIC/NA/II/15feb01/na2_osc.dat                Na2_F_OSCDAT
 ln -sf  $ATOMIC/NA/II/15feb01/f_to_s_ls_21.dat           Na2_F_TO_S
 ln -sf  $ATOMIC/NA/II/15feb01/col_guess.dat              Na2_COL_DATA

 ln -sf  $ATOMIC/NA/III/15feb01/phot_sm_3000.dat          PHOTNaIII_A
 ln -sf  $ATOMIC/NA/III/15feb01/naiiiosc_rev.dat          NaIII_F_OSCDAT
 ln -sf  $ATOMIC/NA/III/15feb01/f_to_s_26.dat             NaIII_F_TO_S
 ln -sf  $ATOMIC/NA/III/15feb01/col_guess.dat             NaIII_COL_DATA

 ln -sf  $ATOMIC/NA/IV/15feb01/phot_sm_3000.dat           PHOTNaIV_A
 ln -sf  $ATOMIC/NA/IV/15feb01/naivosc_rev.dat            NaIV_F_OSCDAT
 ln -sf  $ATOMIC/NA/IV/15feb01/f_to_s_41.dat              NaIV_F_TO_S
 ln -sf  $ATOMIC/NA/IV/15feb01/col_data.dat               NaIV_COL_DATA

 ln -sf  $ATOMIC/NA/V/20jun01/phot_sm_3000.dat            PHOTNaV_A
 ln -sf  $ATOMIC/NA/V/20jun01/navosc_rev.dat              NaV_F_OSCDAT
 ln -sf  $ATOMIC/NA/V/20jun01/f_to_s_43.dat               NaV_F_TO_S
 ln -sf  $ATOMIC/NA/V/20jun01/col_data.dat                NaV_COL_DATA

 ln -sf  $ATOMIC/NA/VI/20jun01/phot_sm_3000.dat           PHOTNaSIX_A
 ln -sf  $ATOMIC/NA/VI/20jun01/naviosc_rev.dat            NaSIX_F_OSCDAT
 ln -sf  $ATOMIC/NA/VI/20jun01/f_to_s_52.dat              NaSIX_F_TO_S
 ln -sf  $ATOMIC/NA/VI/20jun01/col_data.dat               NaSIX_COL_DATA

 ln -sf  $ATOMIC/NA/VII/20jun01/phot_sm_3000.dat          PHOTNaSEV_A
 ln -sf  $ATOMIC/NA/VII/20jun01/naviiosc_rev.dat          NaSEV_F_OSCDAT
 ln -sf  $ATOMIC/NA/VII/20jun01/f_to_s_37.dat             NaSEV_F_TO_S
 ln -sf  $ATOMIC/NA/VII/20jun01/col_data.dat              NaSEV_COL_DATA

#*****************************************************************************
# Magnesium
#-----------
#*****************************************************************************

 ln -sf $ATOMIC/MG/I/5aug97/mgi_osc_split.dat             MgI_F_OSCDAT
 ln -sf $ATOMIC/MG/I/5aug97/mgi_phot_a.dat                PHOTMgI_A
 ln -sf $ATOMIC/MG/I/5aug97/mgicol.dat                    MgI_COL_DATA
 ln -sf $ATOMIC/MG/I/5aug97/mgi_f_to_s_sm.dat             MgI_F_TO_S

 ln -sf $ATOMIC/MG/II/30oct12/mg2_osc_split.dat           Mg2_F_OSCDAT
 ln -sf $ATOMIC/MG/II/30oct12/mg2_phot_a.dat              PHOTMg2_A
 ln -sf $ATOMIC/MG/II/30oct12/mg2col.dat                  Mg2_COL_DATA
 ln -sf $ATOMIC/MG/II/30oct12/mg2_f_to_s_31               Mg2_F_TO_S

 ln -sf  $ATOMIC/MG/III/20jun01/phot_sm_3000.dat          PHOTMgIII_A
 ln -sf  $ATOMIC/MG/III/20jun01/mgiiiosc_rev.dat          MgIII_F_OSCDAT
 ln -sf  $ATOMIC/MG/III/20jun01/f_to_s_41.dat             MgIII_F_TO_S
 ln -sf  $ATOMIC/MG/III/20jun01/col_guess.dat             MgIII_COL_DATA

 ln -sf  $ATOMIC/MG/IV/20jun01/phot_sm_3000.dat           PHOTMgIV_A
 ln -sf  $ATOMIC/MG/IV/20jun01/mgivosc_rev.dat            MgIV_F_OSCDAT
 ln -sf  $ATOMIC/MG/IV/20jun01/f_to_s_27.dat              MgIV_F_TO_S
 ln -sf  $ATOMIC/MG/IV/20jun01/col_guess.dat              MgIV_COL_DATA
#
 ln -sf  $ATOMIC/MG/V/20jun01/phot_sm_3000.dat            PHOTMgV_A
 ln -sf  $ATOMIC/MG/V/20jun01/mgvosc_rev.dat              MgV_F_OSCDAT
 ln -sf  $ATOMIC/MG/V/20jun01/f_to_s_43.dat               MgV_F_TO_S
 ln -sf  $ATOMIC/MG/V/20jun01/col_data.dat                MgV_COL_DATA

 ln -sf  $ATOMIC/MG/VI/20jun01/phot_sm_3000.dat           PHOTMgSIX_A
 ln -sf  $ATOMIC/MG/VI/20jun01/mgviosc_rev.dat            MgSIX_F_OSCDAT
 ln -sf  $ATOMIC/MG/VI/20jun01/f_to_s_46.dat              MgSIX_F_TO_S
 ln -sf  $ATOMIC/MG/VI/20jun01/col_data.dat               MgSIX_COL_DATA

 ln -sf  $ATOMIC/MG/VII/20jun01/phot_sm_3000.dat          PHOTMgSEV_A
 ln -sf  $ATOMIC/MG/VII/20jun01/mgviiosc_rev.dat          MgSEV_F_OSCDAT
 ln -sf  $ATOMIC/MG/VII/20jun01/f_to_s_54.dat             MgSEV_F_TO_S
 ln -sf  $ATOMIC/MG/VII/20jun01/col_data.dat              MgSEV_COL_DATA

#*****************************************************************************
#   Aluminium
#   --------
#*****************************************************************************

 ln -sf $ATOMIC/AL/I/29jul10/fin_osc                      AlI_F_OSCDAT
 ln -sf $ATOMIC/AL/I/29jul10/phot_smooth_0                PHOTAlI_A
 ln -sf $ATOMIC/AL/I/29jul10/col_data                     AlI_COL_DATA
 ln -sf $ATOMIC/AL/I/29jul10/f_to_s_133                   AlI_F_TO_S

 ln -sf $ATOMIC/AL/II/5aug97/al2_osc_split.dat            Al2_F_OSCDAT
 ln -sf $ATOMIC/AL/II/5aug97/al2_phot_a.dat               PHOTAl2_A
 ln -sf $ATOMIC/AL/II/5aug97/al2col.dat                   Al2_COL_DATA
 ln -sf $ATOMIC/AL/II/5aug97/al2_f_to_s_sm.dat            Al2_F_TO_S

 ln -sf  $ATOMIC/AL/III/30oct12/aliii_phot_a.dat          PHOTAlIII_A
 ln -sf  $ATOMIC/AL/III/30oct12/aliii_osc_split.dat       AlIII_F_OSCDAT
 ln -sf  $ATOMIC/AL/III/30oct12/aliii_f_to_s_31           AlIII_F_TO_S
 ln -sf  $ATOMIC/AL/III/30oct12/aliii_col_data.dat        AlIII_COL_DATA

#*****************************************************************************
#    Silicon
#   --------
#*****************************************************************************

 ln -sf $ATOMIC/SIL/I/23nov11/SiI_OSC                     SkI_F_OSCDAT
 ln -sf $ATOMIC/SIL/I/23nov11/SiI_PHOT_DATA               PHOTSkI_A
 ln -sf $ATOMIC/SIL/I/23nov11/col_data                    SkI_COL_DATA
 ln -sf $ATOMIC/SIL/I/23nov11/f_to_s_ls                   SkI_F_TO_S

 ln -sf $ATOMIC/SIL/II/30oct12/phot_op.dat                PHOTSk2_A
 ln -sf $ATOMIC/SIL/II/30oct12/si2_osc_nahar              Sk2_F_OSCDAT
 ln -sf $ATOMIC/SIL/II/30oct12/f_to_s_ls.dat              Sk2_F_TO_S
 ln -sf $ATOMIC/SIL/II/30oct12/si2_col                    Sk2_COL_DATA

 ln -sf $ATOMIC/SIL/III/5dec96b/osc_op_split_rev.dat      SkIII_F_OSCDAT
 ln -sf $ATOMIC/SIL/III/5dec96b/phot_op.dat               PHOTSkIII_A
 ln -sf $ATOMIC/SIL/III/5dec96b/col_data                  SkIII_COL_DATA
 ln -sf $ATOMIC/SIL/III/5dec96b/f_to_s_split.dat          SkIII_F_TO_S

 ln -sf $ATOMIC/SIL/IV/30oct12/phot_op.dat                PHOTSkIV_A
 ln -sf $ATOMIC/SIL/IV/30oct12/osc_op_split.dat           SkIV_F_OSCDAT
 ln -sf $ATOMIC/SIL/IV/30oct12/f_to_s_split.dat           SkIV_F_TO_S
 ln -sf $ATOMIC/SIL/IV/30oct12/col_data.dat               SkIV_COL_DATA

 ln -sf  $ATOMIC/SIL/V/30mar02/phot_sm_3000.dat           PHOTSkV_A
 ln -sf  $ATOMIC/SIL/V/30mar02/skv_osc.dat                SkV_F_OSCDAT
 ln -sf  $ATOMIC/SIL/V/30mar02/f_to_s_52.dat              SkV_F_TO_S
 ln -sf  $ATOMIC/SIL/V/30mar02/col_guess.dat              SkV_COL_DATA

 ln -sf  $ATOMIC/SIL/VI/30mar02/phot_sm_3000.dat          PHOTSkSIX_A
 ln -sf  $ATOMIC/SIL/VI/30mar02/skvi_osc.dat              SkSIX_F_OSCDAT
 ln -sf  $ATOMIC/SIL/VI/30mar02/f_to_s_66.dat             SkSIX_F_TO_S
 ln -sf  $ATOMIC/SIL/VI/30mar02/col_guess.dat             SkSIX_COL_DATA

#*****************************************************************************
# Phosphorous
# -------
# -------
#*****************************************************************************

 ln -sf  $ATOMIC/PHOS/III/15feb01/phot_data.dat             PHOTPIII_A
 ln -sf  $ATOMIC/PHOS/III/15feb01/osc_op.dat                PIII_F_OSCDAT
 ln -sf  $ATOMIC/PHOS/III/15feb01/f_to_s_36.dat             PIII_F_TO_S
 ln -sf  $ATOMIC/PHOS/III/15feb01/col_guess.dat             PIII_COL_DATA

 ln -sf  $ATOMIC/PHOS/IV/15feb01/phot_data_a.dat            PHOTPIV_A
 ln -sf  $ATOMIC/PHOS/IV/15feb01/phot_data_b.dat            PHOTPIV_B
 ln -sf  $ATOMIC/PHOS/IV/15feb01/pivosc_rev.dat             PIV_F_OSCDAT
 ln -sf  $ATOMIC/PHOS/IV/15feb01/f_to_s_36.dat              PIV_F_TO_S
 ln -sf  $ATOMIC/PHOS/IV/15feb01/col_guess.dat              PIV_COL_DATA

 ln -sf  $ATOMIC/PHOS/V/15feb01/phot_data.dat               PHOTPV_A
 ln -sf  $ATOMIC/PHOS/V/15feb01/pvosc_rev.dat               PV_F_OSCDAT
 ln -sf  $ATOMIC/PHOS/V/15feb01/f_to_s_16.dat               PV_F_TO_S
 ln -sf  $ATOMIC/PHOS/V/15feb01/col_guess.dat               PV_COL_DATA

 ln -sf  $ATOMIC/PHOS/VI/15feb01/phot_op.dat                PHOTPSIX_A
 ln -sf  $ATOMIC/PHOS/VI/15feb01/osc_op.dat                 PSIX_F_OSCDAT
 ln -sf  $ATOMIC/PHOS/VI/15feb01/f_to_s_25.dat              PSIX_F_TO_S
 ln -sf  $ATOMIC/PHOS/VI/15feb01/col_guess.dat              PSIX_COL_DATA

#*****************************************************************************
# Sulpher
# -------
#*****************************************************************************

 ln -sf $ATOMIC/SUL/I/24nov11/SI_OSC                       SI_F_OSCDAT
 ln -sf $ATOMIC/SUL/I/24nov11/f_to_s_ls                    SI_F_TO_S
 ln -sf $ATOMIC/SUL/I/24nov11/SI_PHOT_DATA                 PHOTSI_A
 ln -sf $ATOMIC/SUL/I/24nov11/col_data                     SI_COL_DATA

 ln -sf  $ATOMIC/SUL/II/30oct12/phot_nosm                  PHOTS2_A
 ln -sf  $ATOMIC/SUL/II/30oct12/s2_osc                     S2_F_OSCDAT
 ln -sf  $ATOMIC/SUL/II/30oct12/f_to_s_56.dat              S2_F_TO_S
 ln -sf  $ATOMIC/SUL/II/30oct12/s2_col                     S2_COL_DATA

 ln -sf  $ATOMIC/SUL/III/30oct12/phot_nosm                 PHOTSIII_A
 ln -sf  $ATOMIC/SUL/III/30oct12/siiiosc_fin               SIII_F_OSCDAT
 ln -sf  $ATOMIC/SUL/III/30oct12/f_to_s_127                SIII_F_TO_S
 ln -sf  $ATOMIC/SUL/III/30oct12/col_siii                  SIII_COL_DATA

 ln -sf  $ATOMIC/SUL/IV/19nov07/phot_nosm                  PHOTSIV_A
 ln -sf  $ATOMIC/SUL/IV/19nov07/sivosc_fin                 SIV_F_OSCDAT
 ln -sf  $ATOMIC/SUL/IV/19nov07/f_to_s_69                  SIV_F_TO_S
 ln -sf  $ATOMIC/SUL/IV/19nov07/col_siv                    SIV_COL_DATA

 ln -sf  $ATOMIC/SUL/V/3oct00/phot_sm_3000.dat             PHOTSV_A
 ln -sf  $ATOMIC/SUL/V/3oct00/svosc_fin.dat                SV_F_OSCDAT
 ln -sf  $ATOMIC/SUL/V/3oct00/f_to_s_50.dat                SV_F_TO_S
 ln -sf  $ATOMIC/SUL/V/3oct00/col_sv.dat                   SV_COL_DATA

 ln -sf  $ATOMIC/SUL/VI/3oct00/phot_sm_3000.dat            PHOTSSIX_A
 ln -sf  $ATOMIC/SUL/VI/3oct00/sviosc_fin.dat              SSIX_F_OSCDAT
 ln -sf  $ATOMIC/SUL/VI/3oct00/f_to_s_33.dat               SSIX_F_TO_S
 ln -sf  $ATOMIC/SUL/VI/3oct00/col_guess.dat               SSIX_COL_DATA

#*****************************************************************************
# Chlorine
# --------
#*****************************************************************************

 ln -sf  $ATOMIC/CHL/IV/15feb01/phot_data.dat              PHOTClIV_A
 ln -sf  $ATOMIC/CHL/IV/15feb01/f_to_s_58.dat              ClIV_F_TO_S
 ln -sf  $ATOMIC/CHL/IV/15feb01/clivosc_fin.dat            ClIV_F_OSCDAT
 ln -sf  $ATOMIC/CHL/IV/15feb01/col_data.dat               ClIV_COL_DATA

 ln -sf  $ATOMIC/CHL/V/15feb01/phot_data.dat               PHOTClV_A
 ln -sf  $ATOMIC/CHL/V/15feb01/f_to_s_41.dat               ClV_F_TO_S
 ln -sf  $ATOMIC/CHL/V/15feb01/clvosc_fin.dat              ClV_F_OSCDAT
 ln -sf  $ATOMIC/CHL/V/15feb01/col_data.dat                ClV_COL_DATA

 ln -sf  $ATOMIC/CHL/VI/15feb01/phot_data.dat              PHOTClSIX_A
 ln -sf  $ATOMIC/CHL/VI/15feb01/f_to_s_32.dat              ClSIX_F_TO_S
 ln -sf  $ATOMIC/CHL/VI/15feb01/clviosc_rev.dat            ClSIX_F_OSCDAT
 ln -sf  $ATOMIC/CHL/VI/15feb01/col_guess.dat              ClSIX_COL_DATA

 ln -sf  $ATOMIC/CHL/VII/15feb01/phot_data.dat             PHOTClSEV_A
 ln -sf  $ATOMIC/CHL/VII/15feb01/f_to_s_33.dat             ClSEV_F_TO_S
 ln -sf  $ATOMIC/CHL/VII/15feb01/clviiosc_rev.dat          ClSEV_F_OSCDAT
 ln -sf  $ATOMIC/CHL/VII/15feb01/col_guess.dat             ClSEV_COL_DATA

##*****************************************************************************
# Argon
# -----
#*****************************************************************************

 ln -sf  $ATOMIC/ARG/I/9sep11/phot_nosm                    PHOTArI_A
 ln -sf  $ATOMIC/ARG/I/9sep11/fin_osc                      ArI_F_OSCDAT
 ln -sf  $ATOMIC/ARG/I/9sep11/f_to_s_75                    ArI_F_TO_S
 ln -sf  $ATOMIC/ARG/I/9sep11/col_guess                    ArI_COL_DATA

 ln -sf  $ATOMIC/ARG/II/9sep11/phot_nosm                   PHOTAr2_A
 ln -sf  $ATOMIC/ARG/II/9sep11/fin_osc                     Ar2_F_OSCDAT
 ln -sf  $ATOMIC/ARG/II/9sep11/f_to_s_379                  Ar2_F_TO_S
 ln -sf  $ATOMIC/ARG/II/9sep11/col_data                    Ar2_COL_DATA

 ln -sf  $ATOMIC/ARG/III/19nov07/phot_nosm                 PHOTArIII_A
 ln -sf  $ATOMIC/ARG/III/19nov07/fin_osc                   ArIII_F_OSCDAT
 ln -sf  $ATOMIC/ARG/III/19nov07/f_to_s_32                 ArIII_F_TO_S
 ln -sf  $ATOMIC/ARG/III/19nov07/col_ariii                 ArIII_COL_DATA

 ln -sf  $ATOMIC/ARG/IV/1dec99/phot_sm_3000.dat            PHOTArIV_A
 ln -sf  $ATOMIC/ARG/IV/1dec99/fin_osc.dat                 ArIV_F_OSCDAT
 ln -sf  $ATOMIC/ARG/IV/1dec99/f_to_s_50.dat               ArIV_F_TO_S
 ln -sf  $ATOMIC/ARG/IV/1dec99/col_data.dat                ArIV_COL_DATA
#
 ln -sf  $ATOMIC/ARG/V/1dec99/phot_sm_3000.dat             PHOTArV_A
 ln -sf  $ATOMIC/ARG/V/1dec99/fin_osc.dat                  ArV_F_OSCDAT
 ln -sf  $ATOMIC/ARG/V/1dec99/f_to_s_64.dat                ArV_F_TO_S
 ln -sf  $ATOMIC/ARG/V/1dec99/col_data.dat                 ArV_COL_DATA

 ln -sf  $ATOMIC/ARG/VI/15feb01/phot_sm_3000.dat           PHOTArSIX_A
 ln -sf  $ATOMIC/ARG/VI/15feb01/arviosc_rev.dat            ArSIX_F_OSCDAT
 ln -sf  $ATOMIC/ARG/VI/15feb01/f_to_s_30.dat              ArSIX_F_TO_S
 ln -sf  $ATOMIC/ARG/VI/15feb01/col_data.dat               ArSIX_COL_DATA

 ln -sf  $ATOMIC/ARG/VII/15feb01/phot_sm_3000.dat          PHOTArSEV_A
 ln -sf  $ATOMIC/ARG/VII/15feb01/arviiosc_rev.dat          ArSEV_F_OSCDAT
 ln -sf  $ATOMIC/ARG/VII/15feb01/f_to_s_46.dat             ArSEV_F_TO_S
 ln -sf  $ATOMIC/ARG/VII/15feb01/col_data.dat              ArSEV_COL_DATA

 ln -sf  $ATOMIC/ARG/VIII/15feb01/phot_sm_3000.dat         PHOTArVIII_A
 ln -sf  $ATOMIC/ARG/VIII/15feb01/arviiiosc_rev.dat        ArVIII_F_OSCDAT
 ln -sf  $ATOMIC/ARG/VIII/15feb01/f_to_s_33.dat            ArVIII_F_TO_S
 ln -sf  $ATOMIC/ARG/VIII/15feb01/col_guess.dat            ArVIII_COL_DATA

#*****************************************************************************
# Potasium
# -------
#*****************************************************************************

 ln -sf  $ATOMIC/POT/I/4mar12/phot_ki                       PHOTKI_A
 ln -sf  $ATOMIC/POT/I/4mar12/fin_osc                       KI_F_OSCDAT
 ln -sf  $ATOMIC/POT/I/4mar12/f_to_s_46                     KI_F_TO_S
 ln -sf  $ATOMIC/POT/I/4mar12/COL_DATA                      KI_COL_DATA

 ln -sf  $ATOMIC/POT/II/4mar12/phot_k2                      PHOTK2_A
 ln -sf  $ATOMIC/POT/II/4mar12/fin_osc                      K2_F_OSCDAT
 ln -sf  $ATOMIC/POT/II/4mar12/f_to_s_82                    K2_F_TO_S
 ln -sf  $ATOMIC/POT/II/4mar12/COL_DATA                     K2_COL_DATA

 ln -sf  $ATOMIC/POT/III/15feb01/phot_data.dat              PHOTKIII_A
 ln -sf  $ATOMIC/POT/III/15feb01/kiiiosc_rev.dat            KIII_F_OSCDAT
 ln -sf  $ATOMIC/POT/III/15feb01/f_to_s_ls_20.dat           KIII_F_TO_S
 ln -sf  $ATOMIC/POT/III/15feb01/col_guess.dat              KIII_COL_DATA

 ln -sf  $ATOMIC/POT/IV/15feb01/phot_data.dat               PHOTKIV_A
 ln -sf  $ATOMIC/POT/IV/15feb01/fin_osc.dat                 KIV_F_OSCDAT
 ln -sf  $ATOMIC/POT/IV/15feb01/f_to_s_40.dat               KIV_F_TO_S
 ln -sf  $ATOMIC/POT/IV/15feb01/col_data.dat                KIV_COL_DATA

 ln -sf  $ATOMIC/POT/V/15feb01/phot_data.dat                PHOTKV_A
 ln -sf  $ATOMIC/POT/V/15feb01/fin_osc.dat                  KV_F_OSCDAT
 ln -sf  $ATOMIC/POT/V/15feb01/f_to_s_36.dat                KV_F_TO_S
 ln -sf  $ATOMIC/POT/V/15feb01/col_data.dat                 KV_COL_DATA

 ln -sf  $ATOMIC/POT/VI/15feb01/phot_data.dat               PHOTKSIX_A
 ln -sf  $ATOMIC/POT/VI/15feb01/kviosc_rev.dat              KSIX_F_OSCDAT
 ln -sf  $ATOMIC/POT/VI/15feb01/f_to_s_54.dat               KSIX_F_TO_S
 ln -sf  $ATOMIC/POT/VI/15feb01/col_guess.dat               KSIX_COL_DATA

 ln -sf  $ATOMIC/POT/VII/15feb01/phot_smooth.dat            PHOTKSEV_A
 ln -sf  $ATOMIC/POT/VII/15feb01/osc_op_sp.dat              KSEV_F_OSCDAT
 ln -sf  $ATOMIC/POT/VII/15feb01/f_to_s.dat                 KSEV_F_TO_S
 ln -sf  $ATOMIC/POT/VII/15feb01/col_guess.dat              KSEV_COL_DATA

#*****************************************************************************
# Calcium
# -------
#*****************************************************************************

 ln -sf $ATOMIC/CA/I/5aug97/cai_osc_split.dat              CaI_F_OSCDAT
 ln -sf $ATOMIC/CA/I/5aug97/cai_phot_a.dat                 PHOTCaI_A
 ln -sf $ATOMIC/CA/I/5aug97/caicol.dat                     CaI_COL_DATA
 ln -sf $ATOMIC/CA/I/5aug97/caI_f_to_s_term.dat            CaI_F_TO_S

 ln -sf  $ATOMIC/CA/II/30oct12/ca2_phot_a.dat              PHOTCa2_A
 ln -sf  $ATOMIC/CA/II/30oct12/ca2_osc_split.dat           Ca2_F_OSCDAT
 ln -sf  $ATOMIC/CA/II/30oct12/ca2_f_to_s_sm.dat           Ca2_F_TO_S
 ln -sf  $ATOMIC/CA/II/30oct12/ca2col.dat                  Ca2_COL_DATA

 ln -sf  $ATOMIC/CA/III/10apr99/phot_smooth.dat            PHOTCaIII_A
 ln -sf  $ATOMIC/CA/III/10apr99/osc_op_sp.dat              CaIII_F_OSCDAT
 ln -sf  $ATOMIC/CA/III/10apr99/f_to_s.dat                 CaIII_F_TO_S
 ln -sf  $ATOMIC/CA/III/10apr99/col_guess.dat              CaIII_COL_DATA

 ln -sf  $ATOMIC/CA/IV/10apr99/phot_smooth.dat             PHOTCaIV_A
 ln -sf  $ATOMIC/CA/IV/10apr99/osc_op_sp.dat               CaIV_F_OSCDAT
 ln -sf  $ATOMIC/CA/IV/10apr99/f_to_s.dat                  CaIV_F_TO_S
 ln -sf  $ATOMIC/CA/IV/10apr99/col_guess.dat               CaIV_COL_DATA

 ln -sf  $ATOMIC/CA/V/10apr99/phot_smooth.dat              PHOTCaV_A
 ln -sf  $ATOMIC/CA/V/10apr99/osc_op_sp.dat                CaV_F_OSCDAT
 ln -sf  $ATOMIC/CA/V/10apr99/f_to_s.dat                   CaV_F_TO_S
 ln -sf  $ATOMIC/CA/V/10apr99/col_guess.dat                CaV_COL_DATA

 ln -sf  $ATOMIC/CA/VI/10apr99/phot_smooth.dat             PHOTCaSIX_A
 ln -sf  $ATOMIC/CA/VI/10apr99/osc_op_sp.dat               CaSIX_F_OSCDAT
 ln -sf  $ATOMIC/CA/VI/10apr99/f_to_s.dat                  CaSIX_F_TO_S
 ln -sf  $ATOMIC/CA/VI/10apr99/col_guess.dat               CaSIX_COL_DATA

 ln -sf  $ATOMIC/CA/VII/10apr99/phot_smooth.dat            PHOTCaSEV_A
 ln -sf  $ATOMIC/CA/VII/10apr99/osc_op_sp.dat              CaSEV_F_OSCDAT
 ln -sf  $ATOMIC/CA/VII/10apr99/f_to_s_55.dat              CaSEV_F_TO_S
 ln -sf  $ATOMIC/CA/VII/10apr99/col_guess.dat              CaSEV_COL_DATA

 ln -sf  $ATOMIC/CA/VIII/11jun01/phot_data.dat             PHOTCaVIII_A
 ln -sf  $ATOMIC/CA/VIII/11jun01/f_to_s_54.dat             CaVIII_F_TO_S
 ln -sf  $ATOMIC/CA/VIII/11jun01/caviii_osc.dat            CaVIII_F_OSCDAT
 ln -sf  $ATOMIC/CA/VIII/11jun01/col_guess.dat             CaVIII_COL_DATA

 ln -sf  $ATOMIC/CA/IX/11jun01/phot_data.dat               PHOTCaIX_A
 ln -sf  $ATOMIC/CA/IX/11jun01/f_to_s_50.dat               CaIX_F_TO_S
 ln -sf  $ATOMIC/CA/IX/11jun01/caix_osc.dat                CaIX_F_OSCDAT
 ln -sf  $ATOMIC/CA/IX/11jun01/col_guess.dat               CaIX_COL_DATA

 ln -sf  $ATOMIC/CA/X/30mar02/phot_sm_3000.dat             PHOTCaX_A
 ln -sf  $ATOMIC/CA/X/30mar02/f_to_s_31.dat                CaX_F_TO_S
 ln -sf  $ATOMIC/CA/X/30mar02/cax_osc.dat                  CaX_F_OSCDAT
 ln -sf  $ATOMIC/CA/X/30mar02/col_guess.dat                CaX_COL_DATA


# Scandium
#---------

 ln -sf  $SN_ATOMIC/Sc/I/4nov09/fin_osc                 ScI_F_OSCDAT
 ln -sf  $SN_ATOMIC/Sc/I/4nov09/phot_nosm               PHOTScI_A
 ln -sf  $SN_ATOMIC/Sc/I/4nov09/col_guess               ScI_COL_DATA
 ln -sf  $SN_ATOMIC/Sc/I/4nov09/f_to_s_71               ScI_F_TO_S
 ln -sf  $SN_ATOMIC/Sc/I/4nov09/die_sci                 DIEScI

 ln -sf  $SN_ATOMIC/Sc/II/4nov09/fin_osc                Sc2_F_OSCDAT
 ln -sf  $SN_ATOMIC/Sc/II/4nov09/phot_nosm              PHOTSc2_A
 ln -sf  $SN_ATOMIC/Sc/II/4nov09/col_guess              Sc2_COL_DATA
 ln -sf  $SN_ATOMIC/Sc/II/4nov09/f_to_s_110             Sc2_F_TO_S
 ln -sf  $SN_ATOMIC/Sc/II/4nov09/die_sc2                DIESc2

 ln -sf  $SN_ATOMIC/Sc/III/4nov09/fin_osc               ScIII_F_OSCDAT
 ln -sf  $SN_ATOMIC/Sc/III/4nov09/phot_nosm             PHOTScIII_A
 ln -sf  $SN_ATOMIC/Sc/III/4nov09/col_guess             ScIII_COL_DATA
 ln -sf  $SN_ATOMIC/Sc/III/4nov09/f_to_s_35             ScIII_F_TO_S


#Titanium
#--------

ln -sf $ATOMIC/TIT/II/18oct00/tkii_osc.dat              Tk2_F_OSCDAT
ln -sf $ATOMIC/TIT/II/18oct00/f_to_s_61.dat             Tk2_F_TO_S
ln -sf $ATOMIC/TIT/II/18oct00/phot_data.dat             PHOTTk2_A
ln -sf $ATOMIC/TIT/II/18oct00/col_guess.dat             Tk2_COL_DATA

ln -sf $ATOMIC/TIT/III/18oct00/tkiii_osc.dat            TkIII_F_OSCDAT
ln -sf $ATOMIC/TIT/III/18oct00/f_to_s_49.dat            TkIII_F_TO_S
ln -sf $ATOMIC/TIT/III/18oct00/phot_data.dat            PHOTTkIII_A
ln -sf $ATOMIC/TIT/III/18oct00/col_guess.dat            TkIII_COL_DATA

ln -sf $ATOMIC/TIT/IV/18oct00/col_guess.dat             TkIV_COL_DATA
ln -sf $ATOMIC/TIT/IV/18oct00/f_to_s_18.dat             TkIV_F_TO_S
ln -sf $ATOMIC/TIT/IV/18oct00/tkiv_osc.dat              TkIV_F_OSCDAT
ln -sf $ATOMIC/TIT/IV/18oct00/phot_data.dat             PHOTTkIV_A

# Vanadium
#----------

 ln -sf  $ATOMIC/VAN/I/27may10/col_guess.dat            VI_COL_DATA
 ln -sf  $ATOMIC/VAN/I/27may10/vi_osc                   VI_F_OSCDAT
 ln -sf  $ATOMIC/VAN/I/27may10/vi_phot                  PHOTVI_A

#*****************************************************************************
# Chromium
# --------
#*****************************************************************************

 ln -sf $ATOMIC/CHRO/I/10aug12/cri_osc.dat                  CrI_F_OSCDAT
 ln -sf $ATOMIC/CHRO/I/10aug12/f_to_s_101                   CrI_F_TO_S
 ln -sf $ATOMIC/CHRO/I/10aug12/phot_data.dat                PHOTCrI_A
 ln -sf $ATOMIC/CHRO/I/10aug12/col_guess.dat                CrI_COL_DATA

 ln -sf $ATOMIC/CHRO/II/15aug12/crii_osc.dat                Cr2_F_OSCDAT
 ln -sf $ATOMIC/CHRO/II/15aug12/f_to_s_84.dat               Cr2_F_TO_S
 ln -sf $ATOMIC/CHRO/II/15aug12/phot_data.dat               PHOTCr2_A
 ln -sf $ATOMIC/CHRO/II/15aug12/col_data.dat                Cr2_COL_DATA

 ln -sf $ATOMIC/CHRO/III/18oct00/criii_osc.dat              CrIII_F_OSCDAT
 ln -sf $ATOMIC/CHRO/III/18oct00/f_to_s_68.dat              CrIII_F_TO_S
 ln -sf $ATOMIC/CHRO/III/18oct00/phot_data.dat              PHOTCrIII_A
 ln -sf $ATOMIC/CHRO/III/18oct00/col_guess.dat              CrIII_COL_DATA

 ln -sf  $ATOMIC/CHRO/IV/18oct00/phot_data.dat              PHOTCrIV_A
 ln -sf  $ATOMIC/CHRO/IV/18oct00/f_to_s_48.dat              CrIV_F_TO_S
 ln -sf  $ATOMIC/CHRO/IV/18oct00/criv_osc.dat               CrIV_F_OSCDAT
 ln -sf  $ATOMIC/CHRO/IV/18oct00/col_guess.dat              CrIV_COL_DATA

 ln -sf  $ATOMIC/CHRO/V/18oct00/phot_data.dat               PHOTCrV_A
 ln -sf  $ATOMIC/CHRO/V/18oct00/f_to_s_51.dat               CrV_F_TO_S
 ln -sf  $ATOMIC/CHRO/V/18oct00/crv_osc.dat                 CrV_F_OSCDAT
 ln -sf  $ATOMIC/CHRO/V/18oct00/col_guess.dat               CrV_COL_DATA

 ln -sf  $ATOMIC/CHRO/VI/18oct00/phot_data.dat              PHOTCrSIX_A
 ln -sf  $ATOMIC/CHRO/VI/18oct00/f_to_s_30.dat              CrSIX_F_TO_S
 ln -sf  $ATOMIC/CHRO/VI/18oct00/crvi_osc.dat               CrSIX_F_OSCDAT
 ln -sf  $ATOMIC/CHRO/VI/18oct00/col_guess.dat              CrSIX_COL_DATA

#*****************************************************************************
# Maganese
# --------
#*****************************************************************************

 ln -sf  $ATOMIC/MAN/II/18oct00/mnii_osc.dat                Mn2_F_OSCDAT
 ln -sf  $ATOMIC/MAN/II/18oct00/f_to_s_92.dat               Mn2_F_TO_S
 ln -sf  $ATOMIC/MAN/II/18oct00/phot_data.dat               PHOTMn2_A
 ln -sf  $ATOMIC/MAN/II/18oct00/col_guess.dat               Mn2_COL_DATA

 ln -sf  $ATOMIC/MAN/III/18oct00/mniii_osc.dat              MnIII_F_OSCDAT
 ln -sf  $ATOMIC/MAN/III/18oct00/f_to_s_65.dat              MnIII_F_TO_S
 ln -sf  $ATOMIC/MAN/III/18oct00/phot_data.dat              PHOTMnIII_A
 ln -sf  $ATOMIC/MAN/III/18oct00/col_guess.dat              MnIII_COL_DATA

 ln -sf  $ATOMIC/MAN/IV/18oct00/phot_data.dat               PHOTMnIV_A
 ln -sf  $ATOMIC/MAN/IV/18oct00/mniv_osc.dat                MnIV_F_OSCDAT
 ln -sf  $ATOMIC/MAN/IV/18oct00/f_to_s_47.dat               MnIV_F_TO_S
 ln -sf  $ATOMIC/MAN/IV/18oct00/col_guess.dat               MnIV_COL_DATA

 ln -sf  $ATOMIC/MAN/V/18oct00/phot_data.dat                PHOTMnV_A
 ln -sf  $ATOMIC/MAN/V/18oct00/mnv_osc.dat                  MnV_F_OSCDAT
 ln -sf  $ATOMIC/MAN/V/18oct00/f_to_s_36.dat                MnV_F_TO_S
 ln -sf  $ATOMIC/MAN/V/18oct00/col_guess.dat                MnV_COL_DATA

 ln -sf  $ATOMIC/MAN/VI/18oct00/phot_data.dat               PHOTMnSIX_A
 ln -sf  $ATOMIC/MAN/VI/18oct00/mnvi_osc.dat                MnSIX_F_OSCDAT
 ln -sf  $ATOMIC/MAN/VI/18oct00/f_to_s_41.dat               MnSIX_F_TO_S
 ln -sf  $ATOMIC/MAN/VI/18oct00/col_guess.dat               MnSIX_COL_DATA

 ln -sf  $ATOMIC/MAN/VII/18oct00/phot_data.dat              PHOTMnSEV_A
 ln -sf  $ATOMIC/MAN/VII/18oct00/mnvii_osc.dat              MnSEV_F_OSCDAT
 ln -sf  $ATOMIC/MAN/VII/18oct00/f_to_s_21.dat              MnSEV_F_TO_S
 ln -sf  $ATOMIC/MAN/VII/18oct00/col_guess.dat              MnSEV_COL_DATA

#*****************************************************************************
#  Iron
#  ----
#*****************************************************************************

 ln -sf  $ATOMIC/FE/I/29apr04/phot_smooth_0                 PHOTFeI_A
 ln -sf  $ATOMIC/FE/I/29apr04/fei_osc                       FeI_F_OSCDAT
 ln -sf  $ATOMIC/FE/I/29apr04/fei_f_to_s.dat                FeI_F_TO_S
 ln -sf  $ATOMIC/FE/I/29apr04/col_data                      FeI_COL_DATA

 ln -sf  $ATOMIC/FE/II/24may96/phot_op.dat                  PHOTFe2_A
 ln -sf  $ATOMIC/FE/II/16nov98/fe2osc_nahar_kurucz.dat      Fe2_F_OSCDAT
 ln -sf  $ATOMIC/FE/II/16nov98/f_to_s_sn                    Fe2_F_TO_S
 ln -sf  $ATOMIC/FE/II/16nov98/fe2_col.dat                  Fe2_COL_DATA

 ln -sf  $ATOMIC/FE/III/30oct12/phot_sm_3000.dat            PHOTFeIII_A
 ln -sf  $ATOMIC/FE/III/30oct12/FeIII_OSC                   FeIII_F_OSCDAT
 ln -sf  $ATOMIC/FE/III/30oct12/f_to_s_105                  FeIII_F_TO_S
 ln -sf  $ATOMIC/FE/III/30oct12/col_data.dat                FeIII_COL_DATA

 ln -sf  $ATOMIC/FE/IV/18oct00/phot_sm_3000.dat             PHOTFeIV_A
 ln -sf  $ATOMIC/FE/IV/18oct00/feiv_osc.dat                 FeIV_F_OSCDAT
 ln -sf  $ATOMIC/FE/IV/18oct00/f_to_s_100.dat               FeIV_F_TO_S
 ln -sf  $ATOMIC/FE/IV/18oct00/col_data.dat                 FeIV_COL_DATA

 ln -sf  $ATOMIC/FE/V/18oct00/phot_sm_3000.dat              PHOTFeV_A
 ln -sf  $ATOMIC/FE/V/18oct00/fev_osc.dat                   FeV_F_OSCDAT
 ln -sf  $ATOMIC/FE/V/18oct00/f_to_s_139.dat                FeV_F_TO_S
 ln -sf  $ATOMIC/FE/V/18oct00/col_guess.dat                 FeV_COL_DATA

 ln -sf  $ATOMIC/FE/VI/18oct00/phot_sm_3000.dat             PHOTFeSIX_A
 ln -sf  $ATOMIC/FE/VI/18oct00/fevi_osc.dat                 FeSIX_F_OSCDAT
 ln -sf  $ATOMIC/FE/VI/18oct00/f_to_s_67.dat                FeSIX_F_TO_S
 ln -sf  $ATOMIC/FE/VI/18oct00/col_guess.dat                FeSIX_COL_DATA

 ln -sf  $ATOMIC/FE/VII/18oct00/phot_sm_3000.dat            PHOTFeSEV_A
 ln -sf  $ATOMIC/FE/VII/18oct00/fevii_osc.dat               FeSEV_F_OSCDAT
 ln -sf  $ATOMIC/FE/VII/18oct00/f_to_s_69.dat               FeSEV_F_TO_S
 ln -sf  $ATOMIC/FE/VII/18oct00/col_guess.dat               FeSEV_COL_DATA

 ln -sf  $ATOMIC/FE/VIII/8may97/feviii_phot_op.dat          PHOTFeVIII_A
 ln -sf  $ATOMIC/FE/VIII/8may97/feviii_f_to_s_53.dat        FeVIII_F_TO_S
 ln -sf  $ATOMIC/FE/VIII/8may97/feviii_osc_kb_rk.dat        FeVIII_F_OSCDAT
 ln -sf  $ATOMIC/FE/VIII/8may97/col_guess.dat               FeVIII_COL_DATA

 ln -sf  $ATOMIC/FE/IX/5may02/f_to_s_52.dat                 FeIX_F_TO_S
 ln -sf  $ATOMIC/FE/IX/5may02/osc_split_op                  FeIX_F_OSCDAT
 ln -sf  $ATOMIC/FE/IX/22apr02/phot_smooth_3000             PHOTFeIX_A
 ln -sf  $ATOMIC/FE/IX/22apr02/col_guess.dat                FeIX_COL_DATA

 ln -sf  $ATOMIC/FE/X/5may02/f_to_s_70.dat                  FeX_F_TO_S
 ln -sf  $ATOMIC/FE/X/5may02/osc_split_op                   FeX_F_OSCDAT
 ln -sf  $ATOMIC/FE/X/5may02/phot_smooth_3000               PHOTFeX_A
 ln -sf  $ATOMIC/FE/X/5may02/col_guess.dat                  FeX_COL_DATA

#*****************************************************************************
#    Cobalt
#   ------
#*****************************************************************************

# Co II

 ln -sf $ATOMIC/COB/II/15nov11/fin_osc_bound           Co2_F_OSCDAT
 ln -sf $ATOMIC/COB/II/15nov11/f_to_s_136              Co2_F_TO_S
 ln -sf $ATOMIC/COB/II/15nov11/phot_nosm_rev           PHOTCo2_A
 ln -sf $ATOMIC/COB/II/15nov11/Co2_COL_DATA            Co2_COL_DATA

# CoIII

 ln -sf $ATOMIC/COB/III/30oct12/coiii_osc.dat          CoIII_F_OSCDAT
 ln -sf $ATOMIC/COB/III/30oct12/f_to_s_124             CoIII_F_TO_S
 ln -sf $ATOMIC/COB/III/30oct12/phot_nosm              PHOTCoIII_A
 ln -sf $ATOMIC/COB/III/30oct12/col_data.dat           CoIII_COL_DATA

# CoIV

 ln -sf $ATOMIC/COB/IV/4jan12/coiv_osc.dat             CoIV_F_OSCDAT
 ln -sf $ATOMIC/COB/IV/4jan12/f_to_s_56.dat            CoIV_F_TO_S
 ln -sf $ATOMIC/COB/IV/4jan12/phot_data                PHOTCoIV_A
 ln -sf $ATOMIC/COB/IV/4jan12/col_data.dat             CoIV_COL_DATA

# CoV

 ln -sf $ATOMIC/COB/V/18oct00/cov_osc.dat              CoV_F_OSCDAT
 ln -sf $ATOMIC/COB/V/18oct00/f_to_s_43.dat            CoV_F_TO_S
 ln -sf $ATOMIC/COB/V/18oct00/phot_data.dat            PHOTCoV_A
 ln -sf $ATOMIC/COB/V/18oct00/col_guess.dat            CoV_COL_DATA

# CoSIX

 ln -sf $ATOMIC/COB/VI/18oct00/covi_osc.dat            CoSIX_F_OSCDAT
 ln -sf $ATOMIC/COB/VI/18oct00/f_to_s_41.dat           CoSIX_F_TO_S
 ln -sf $ATOMIC/COB/VI/18oct00/phot_data.dat           PHOTCoSIX_A
 ln -sf $ATOMIC/COB/VI/18oct00/col_guess.dat           CoSIX_COL_DATA

# CoSEV

 ln -sf $ATOMIC/COB/VII/18oct00/covii_osc.dat          CoSEV_F_OSCDAT
 ln -sf $ATOMIC/COB/VII/18oct00/f_to_s_45.dat          CoSEV_F_TO_S
 ln -sf $ATOMIC/COB/VII/18oct00/phot_data.dat          PHOTCoSEV_A
 ln -sf $ATOMIC/COB/VII/18oct00/col_guess.dat          CoSEV_COL_DATA

#*****************************************************************************
#  Nickel
#  ------
#*****************************************************************************

 ln -sf  $ATOMIC/NICK/II/30oct12/nkii_osc.dat                 Nk2_F_OSCDAT
 ln -sf  $ATOMIC/NICK/II/30oct12/f_to_s_59.dat                Nk2_F_TO_S
 ln -sf  $ATOMIC/NICK/II/30oct12/phot_data                    PHOTNk2_A
 ln -sf  $ATOMIC/NICK/II/30oct12/col_data_bautista            Nk2_COL_DATA

 ln -sf  $ATOMIC/NICK/III/27aug12/phot_data.dat               PHOTNkIII_A
 ln -sf  $ATOMIC/NICK/III/27aug12/f_to_s_47.dat               NkIII_F_TO_S
 ln -sf  $ATOMIC/NICK/III/27aug12/nkiii_osc.dat               NkIII_F_OSCDAT
 ln -sf  $ATOMIC/NICK/III/27aug12/col_data.dat                NkIII_COL_DATA

 ln -sf  $ATOMIC/NICK/IV/18oct00/phot_data.dat                PHOTNkIV_A
 ln -sf  $ATOMIC/NICK/IV/18oct00/f_to_s_54.dat                NkIV_F_TO_S
 ln -sf  $ATOMIC/NICK/IV/18oct00/nkiv_osc.dat                 NkIV_F_OSCDAT
 ln -sf  $ATOMIC/NICK/IV/18oct00/col_guess.dat                NkIV_COL_DATA

 ln -sf  $ATOMIC/NICK/V/18oct00/phot_data.dat                 PHOTNkV_A
 ln -sf  $ATOMIC/NICK/V/18oct00/f_to_s_152.dat                NkV_F_TO_S
 ln -sf  $ATOMIC/NICK/V/18oct00/nkv_osc.dat                   NkV_F_OSCDAT
 ln -sf  $ATOMIC/NICK/V/18oct00/col_guess.dat                 NkV_COL_DATA

 ln -sf  $ATOMIC/NICK/VI/18oct00/phot_data.dat                PHOTNkSIX_A
 ln -sf  $ATOMIC/NICK/VI/18oct00/f_to_s_62.dat                NkSIX_F_TO_S
 ln -sf  $ATOMIC/NICK/VI/18oct00/nkvi_osc.dat                 NkSIX_F_OSCDAT
 ln -sf  $ATOMIC/NICK/VI/18oct00/col_guess.dat                NkSIX_COL_DATA

 ln -sf  $ATOMIC/NICK/VII/18oct00/phot_data.dat               PHOTNkSEV_A
 ln -sf  $ATOMIC/NICK/VII/18oct00/f_to_s_61.dat               NkSEV_F_TO_S
 ln -sf  $ATOMIC/NICK/VII/18oct00/nkvii_osc.dat               NkSEV_F_OSCDAT
 ln -sf  $ATOMIC/NICK/VII/18oct00/col_guess.dat               NkSEV_COL_DATA

 ln -sf  $ATOMIC/NICK/VIII/11jun01/phot_data.dat              PHOTNkVIII_A
 ln -sf  $ATOMIC/NICK/VIII/11jun01/f_to_s_48.dat              NkVIII_F_TO_S
 ln -sf  $ATOMIC/NICK/VIII/11jun01/nkviii_osc.dat             NkVIII_F_OSCDAT
 ln -sf  $ATOMIC/NICK/VIII/11jun01/col_guess.dat              NkVIII_COL_DATA

 ln -sf  $ATOMIC/NICK/IX/11jun01/phot_data.dat                PHOTNkIX_A
 ln -sf  $ATOMIC/NICK/IX/11jun01/f_to_s_48.dat                NkIX_F_TO_S
 ln -sf  $ATOMIC/NICK/IX/11jun01/nkix_osc.dat                 NkIX_F_OSCDAT
 ln -sf  $ATOMIC/NICK/IX/11jun01/col_guess.dat                NkIX_COL_DATA

#*****************************************************************************
#
# END OF ATOMIC DATA SOFT LINKS
#
#*****************************************************************************

#If we pass a paremeter to the BATCH file, we assume that we just want to
#ln -sf  data files. In this case we exit.

if ( $1 == '') then

else
    echo "Data files assigned"
    exit
endif


#***********************************************************************
#    Input files: Only needed is different from defaults
#***********************************************************************

# ln -sf  $MODEL/vadat.dat                                  VADAT
# ln -sf  $MODEL/input.dat                                  IN_ITS
# ln -sf  $MODEL/cfdat_in.dat                               CFDAT

#***********************************************************************
#    Departure coeficient input files. T_IN is required for a new
#    model with GRID=FALSE. Other links not needed unless diferent
#    from default.
#***********************************************************************

ln -sf  C2_IN                 T_IN

#***********************************************************************
#    Output files: If disk space is at a preimu, these two files might
#     be better stored on an alternative disk.
#***********************************************************************

#ln -sf /usr/limey/jdh/BAMAT            BAMAT
#ln -sf /usr/limey/jdh/BAION            BAION

#***********************************************************************
#***********************************************************************
#    Run program
#***********************************************************************
#***********************************************************************

#Putting time stamp etc on program

rm -f batch.log

echo " " > 'batch.log'
echo "Model started at:" `date` >> 'batch.log'
echo "Machine name is" `uname -n`  >> 'batch.log'
echo " " >> 'batch.log'
echo "PID of batch.sh is:" $$ >> 'batch.log'

(exec nice -19 $CMFGEN_PROG  >>& 'batch.log')&
echo " " >> 'batch.log'
echo "PID of cmfgen.exe is:" $! >> 'batch.log'
wait

echo " " >> 'batch.log'
echo "Program finished on:" `date` >> 'batch.log'

#***********************************************************************
# Execute the command in next_job.sh. This will allow another job to
# be started.
#***********************************************************************
