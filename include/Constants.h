#ifndef CONSTANTS
#define CONSTANTS
#include <sstream>
#include <vector>
#include <complex>

using std::string;
using std::vector;

namespace constants{
inline	double twopc1Dptmin=0.0;
inline	double twopc1Dptmax=2.0;
inline	double twopc1DetaRange=2.4;
inline	double twopc1Dptmin_N=0.4;
inline	double twopc1Dptmax_N=1000.0;
inline	double twopc1DetaRange_N=2.4;
inline  double twopc1D_Nmax=1000.;
inline  double twopc1D_Nmin=100.;
inline	double twopc2DetaRange=2.4;
inline	double twopc2Dptmin_N=0.4;
inline	double twopc2Dptmax_N=1000.0;
inline	double twopc2DetaRange_N=2.4;
inline	double trig_ptmin=0.5;
inline	double trig_ptmax=3.5;
inline	double assoc_ptmin=0.5;
inline	double assoc_ptmax=3.5;
inline	const double minpt_Rt=5.0;
inline  const double maxPhi_RtTrans=(2.0/3.0)*M_PI;
inline  const double minPhi_RtTrans=(1.0/3.0)*M_PI;
inline	const double etaRange_Rt=0.8;
inline	const double etaRangeMultiplicity=0.5;
inline	const double etaA_3sub = -0.4;
inline	const double etaB_3sub = 0.4;
inline	const double DeltaEtaFULL=4.8;
inline	const double DeltaEtaNS=2.0;
inline	const double DeltaPhiOUT=2.0;
inline	const double DeltaPhiNS=0.85;
inline	const double etaGap=1.0;
inline const double ptmin_cumulantmulti=0.2;
inline const double ptmax_cumulantmulti=3.0;
inline const double ptmin_cumulantpt=0.2;
inline const double ptmax_cumulantpt=5.0;
inline const double ptmin_cumulanteta=0.0;
inline const double ptmax_cumulanteta=3.0;
inline const double ptmin_cumulantmulti_Nch=0.2;
inline const double ptmax_cumulantmulti_Nch=3.0;
inline const double etaRange_cumulantmulti=0.8;
inline const double etaRange_cumulantpt=0.8;
inline const double etaRange_cumulantmulti_Nch=0.8;

inline const double hbarc = 0.197327;// natural unit in [GeV*fm]=1
inline const double TL_tau_00 = 0.10/hbarc;//GeV^-1
inline const double TL_tau_0 = 0.30/hbarc;//GeV^-1
inline const double TL_tau_switch = 0.30/hbarc;//GeV^-1
inline const double TL_tau_switchFM=0.30;//fm
inline const double TL_dtau1=0.01/hbarc;
inline const double TL_dtau1FM=0.010;
inline const double TL_dtau2=0.3/hbarc;
inline const std::string ext_nameTLxy="hydro_profiles_eta_at0__tau_";
inline const std::string ext_nameTLxeta="hydro_profiles_y_at0__tau_";


#ifdef DNDETA_PROTON
	const std::string MODE = "dndeta_proton";
	const string default_out_directory_name="dndeta_proton";
	const double x_max = 20.0;
	const double x_min = -20.0;
	const double d_x = 0.5;
	const int x_cell_capa=1000;
	const string default_ext="/hadronised.txt";
	const double y_max = 0.0;
	const double y_max_HI = 0.0;
	const double y_min = 0.0;
	const double d_y = 0.0;
	const double d_y_HI = 0.0;
	const int y_cell_capa=0;
#elif RT_SPECTRA
	const std::string MODE = "Rt_spectra";
	const string default_out_directory_name="Rt_spectra";
	const double x_max = 10.0;
	const double x_min = 0.0;
	const double x_max_HI = 10.0;
	const double d_x = 0.5;
	const double d_x_HI = 0.5;
	const int x_cell_capa=1000;
	const string default_ext="/hadronised.txt";
	const double y_max = 0.0;
	const double y_max_HI = 0.0;
	const double y_min = 0.0;
	const double d_y = 0.0;
	const double d_y_HI = 0.0;
	const int y_cell_capa=0;
#elif RT_YIELD
	const std::string MODE = "Rt_yield";
	const string default_out_directory_name="Rt_yield";
	const double x_max = 10.0;
	const double x_min = 0.0;
	const double x_max_HI = 10.0;
	const double d_x = 0.5;
	const double d_x_HI = 0.5;
	const int x_cell_capa=1000;
	const string default_ext="/hadronised.txt";
	const double y_max = 0.0;
	const double y_max_HI = 0.0;
	const double y_min = 0.0;
	const double d_y = 0.0;
	const double d_y_HI = 0.0;
	const int y_cell_capa=0;
#elif defined CUMULANTMULTI
	const std::string MODE = "cumulant_multi";
	const string default_out_directory_name="cumulant_multi";
	const double x_max = 150.0;
	const double x_max_HI = 1500.0;
	const double x_min = 0.0;
	const double d_x = 20.0;
	const double d_x_HI = 50.0;
	const int x_cell_capa=1000;
	const string default_ext="/hadronised.txt";
	const double y_max = 0.0;
	const double y_max_HI = 0.0;
	const double y_min = 0.0;
	const double d_y = 0.0;
	const double d_y_HI = 0.0;
	const int y_cell_capa=0;
#elif defined CUMULANTPT
	const std::string MODE = "cumulant_pt";
	const string default_out_directory_name="cumulant_pt";
	const double x_max = 20.0;
	const double x_max_HI = 20.0;
	const double x_min = 0.0;
	const double d_x = 0.5;
	const double d_x_HI = 0.5;
	const int x_cell_capa=1000;
	const string default_ext="/hadronised.txt";
	const double y_max = 0.0;
	const double y_max_HI = 0.0;
	const double y_min = 0.0;
	const double d_y = 0.0;
	const double d_y_HI = 0.0;
	const int y_cell_capa=0;
#elif defined CUMULANTETA
	const std::string MODE = "cumulant_eta";
	const string default_out_directory_name="cumulant_eta";
	const double x_max = 10.0;
	const double x_max_HI = 10.0;
	const double x_min = -10.0;
	const double d_x = 1.5;
	const double d_x_HI = 1.5;
	const int x_cell_capa=1000;
	const string default_ext="/hadronised.txt";
	const double y_max = 0.0;
	const double y_max_HI = 0.0;
	const double y_min = 0.0;
	const double d_y = 0.0;
	const double d_y_HI = 0.0;
	const int y_cell_capa=0;
#elif defined MTSCALING
	const std::string MODE = "mtscaling";
	const string default_out_directory_name="mtscaling";
	const double x_max = 5.0; //Needs to be changed depends on what you wanna see.
	const double x_max_HI = 5.0;
	const double x_min = 0.0;
	const double d_x = 0.05;
	const double d_x_HI = 0.05;
	const int x_cell_capa=1000;
	const string default_ext="/hadronised.txt";
	const double y_max = 0.0;
	const double y_max_HI = 0.0;
	const double y_min = 0.0;
	const double d_y = 0.0;
	const double d_y_HI = 0.0;
	const int y_cell_capa=0;
#elif defined TWOPC2D
	const std::string MODE = "twopc2D";
	const string default_out_directory_name="twopc2D";
	const double x_max = 10;
	const double x_max_HI = 10;
	const double x_min = -10.0;
	const double y_max = (3.0/2.0)*M_PI;
	const double y_max_HI = (3.0/2.0)*M_PI;
	const double y_min = -(1.0/2.0)*M_PI;
	const double d_x = 0.5;
	const double d_y = (2.0*M_PI)/25.0;
	const double d_x_HI = 0.5;
	const double d_y_HI = (2.0*M_PI)/25.0;
	const int x_cell_capa=1000;
	const int y_cell_capa=1000;
	const string default_ext="/hadronised.txt";
#elif defined TWOPC1D
	const std::string MODE = "twopc1D";
	const string default_out_directory_name="twopc1D";
	const double x_max = (3.0/2.0)*M_PI;
	const double x_max_HI = (3.0/2.0)*M_PI;
	const double x_min = -(1.0/2.0)*M_PI;
	const double d_x = (2.0*M_PI)/25.0;
	const double d_x_HI = (2.0*M_PI)/25.0;
	const int x_cell_capa=1000;
	const string default_ext="/hadronised.txt";
	const double y_max = 0.0;
	const double y_max_HI = 0.0;
	const double y_min = 0.0;
	const double d_y = 0.0;
	const double d_y_HI = 0.0;
	const int y_cell_capa=0;
#elif defined TWOPCINTEG
	const std::string MODE = "twopcInteg";
	const string default_out_directory_name="twopcInteg";
	const double x_max = 30.0;
	const double x_max_HI = 1500.0;
	const double x_min = 0.0;
	const double d_x = 5.0;
	const double d_x_HI = 10.0;
	const int x_cell_capa=1000;
	const string default_ext="/hadronised.txt";
	const double y_max = 0.0;
	const double y_max_HI = 0.0;
	const double y_min = 0.0;
	const double d_y = 0.0;
	const double d_y_HI = 0.0;
	const int y_cell_capa=0;
#elif defined MEANMT
	const std::string MODE = "meanmt";
	const string default_out_directory_name="meanmt";
	const double x_max = 150;
	const double x_max_HI = 1500;
	const double x_min = 0.0;
	const double d_x = 2.5;
	const double d_x_HI = 10.0;
	const int x_cell_capa=1000;
	const string default_ext="/hadronised.txt";
	const double y_max = 0.0;
	const double y_max_HI = 0.0;
	const double y_min = 0.0;
	const double d_y = 0.0;
	const double d_y_HI = 0.0;
	const int y_cell_capa=0;
#elif defined MEANPT
	const std::string MODE = "meanpt";
	const string default_out_directory_name="meanpt";
	const double x_max = 150.0;
	const double x_max_HI = 1500.0;
	const double x_min = 0.0;
	const double d_x = 5.0;
	const double d_x_HI = 10.0;
	const int x_cell_capa=1000;
	const string default_ext="/hadronised.txt";
	const double y_max = 0.0;
	const double y_max_HI = 0.0;
	const double y_min = 0.0;
	const double d_y = 0.0;
	const double d_y_HI = 0.0;
	const int y_cell_capa=0;
#elif defined MEANPTPID
	const std::string MODE = "MeanptPID";
	const string default_out_directory_name="MeanptPID";
	const double x_max = 150.0;
	const double x_max_HI = 1500.0;
	const double x_min = 0.0;
	const double d_x = 5.0;
	const double d_x_HI = 10.0;
	const int x_cell_capa=1000;
	const string default_ext="/hadronised.txt";
	const double y_max = 0.0;
	const double y_max_HI = 0.0;
	const double y_min = 0.0;
	const double d_y = 0.0;
	const double d_y_HI = 0.0;
	const int y_cell_capa=0;
#elif defined TIMELAPSE 
	const std::string MODE = "timelapse";
	const string default_out_directory_name="timelapse";
	const double x_max = 15.0/hbarc;
	const double x_max_HI = 15.0/hbarc;
	const double x_min = 0.0;
	const double d_x = 0.1/hbarc;
	const double d_x_HI = 0.1/hbarc;
	const int x_cell_capa=1000;
	const string default_ext="/hadronised.txt";
	const double y_max = 0.0;
	const double y_max_HI = 0.0;
	const double y_min = 0.0;
	const double d_y = 0.0;
	const double d_y_HI = 0.0;
	const int y_cell_capa=0;
#elif defined VERTICES
	const std::string MODE = "vertices";
	const string default_out_directory_name="vertices";
	const double x_max = 20.0;
	const double x_min = -20.0;
	const double d_x = 0.5;
	const int x_cell_capa=1000;
	const string default_ext="/hadronised.txt";
	const double y_max = 0.0;
	const double y_max_HI = 0.0;
	const double y_min = 0.0;
	const double d_y = 0.0;
	const double d_y_HI = 0.0;
	const int y_cell_capa=0;
#elif defined DNDPT
	const std::string MODE = "dndpt";
	const string default_out_directory_name="dndpt";
	const double x_max = 2000.0;
	const double x_min = 0.0;
	const double d_x = 1.0;
	const int x_cell_capa=1000;
	const string default_ext="/hadronised.txt";
	const double y_max = 0.0;
	const double y_max_HI = 0.0;
	const double y_min = 0.0;
	const double d_y = 0.0;
	const double d_y_HI = 0.0;
	const int y_cell_capa=0;
#elif defined JET_PRAC
	const std::string MODE = "JET_PRAC";
	const string default_out_directory_name="JET_PRAC";
	const double x_max = 1.0;
	const double x_max_HI = 1.0;
	const double x_min = 0.0;
	const double d_x = 0.1;
	const double d_x_HI = 0.1;
	const int x_cell_capa=1000;
	const string default_ext="/jetinfo.txt";
	const double y_max = 0.0;
	const double y_max_HI = 0.0;
	const double y_min = 0.0;
	const double d_y = 0.0;
	const double d_y_HI = 0.0;
	const int y_cell_capa=0;
#elif defined DEF
	const std::string MODE = "default";
	const string default_out_directory_name="output";
	const double x_max = 20.0;
	const double x_min = -20.0;
	const double d_x = 0.5;
	const int x_cell_capa=1000;
	const string default_ext="/hadronised.txt";
#else
	const std::string MODE = "default";
	const string default_out_directory_name="output";
	const double x_max = 20.0;
	const double x_min = -20.0;
	const double d_x = 0.5;
	const int x_cell_capa=1000;
	const string default_ext="/hadronised.txt";
#endif





inline	const int id_nucleus=1000000000;
inline	const int id_gluon=21;
inline	const int id_ch_pion=211;
inline	const int id_ch_kaon=321;
inline	const int id_proton=2212;
inline	const int id_lambda=3122;
inline	const int id_cascade=3312;
inline	const int id_omega=3334;
inline	const int id_phi=333;
inline	const int id_K0S=310;

inline	const double w_eta_multiplicity=0.5;
inline	const double w_eta_multiplicity_INEL_lg_0=1.0;
inline	const double w_eta_ATLAS_cut=2.5;
inline	const double V0M_fwd1=2.8;
inline	const double V0M_fwd2=5.1;
inline	const double V0M_bkw1=-3.7;
inline	const double V0M_bkw2=-1.7;
inline	const double outer_SPD=1.4;
inline  const double etaRangeCMSRidge=2.4;
inline  const double MomentumMinCMSRidge=0.4;

inline	const double delta_phi_SMALL = 0.5;
inline	const double delta_eta = 0.5;


inline	const double pPb_mid_rapidity__fwd=0.50;
inline	const double pPb_mid_rapidity__bkw=0.0;
inline	const double pPb_rap_shift_from_lab_to_cm=0.465;
inline	const double default_midy_pm = 0.5;
inline	const double LARGE=1.0e+16;
inline	const int LARGEint=10000;
inline	const double SMALL=1.0e-5;
inline	const double multip_cut_more_than=-LARGE;
inline	const double multip_cut_less_than=LARGE;

inline	const double wBin_pp_1[7]={61.049758987000004, 0.92521, 0.052516, 0.010366, 0.00092042, 0.00019170, 0.000036893};//2.76TeV
inline	const double wBin_pp_2[7]={67.78102915, 3.3182, 0.23192, 0.053267, 0.0057856, 0.0014325, 0.00036575};//7TeV
inline	const double wBin_pp_1tot = 62.039;
inline	const double wBin_pp_2tot = 71.392;
inline	const int nBin=7;
inline  const int n_NtrkClass_pp=13;
inline	const double val_cent_pp[10]={0.0095,0.047,0.095,0.14,0.19,0.28,0.38,0.48,0.68,1.0};
inline	const double val_NtrkClass_pp[n_NtrkClass_pp]={10, 20, 30, 40, 60, 85, 95, 105, 115, 125, 135, 150, 100000};
inline	const double val_cent_original_narrow[10]={0.0025,0.01,0.05,0.10,0.20,0.30,0.40,0.50,0.70,1.0};
inline	const double val_cent_original[4]={0.1,0.3,0.5,1.0};
inline	const double val_cent_pPb[7]={0.05,0.1,0.2,0.4,0.6,0.8,1.0};
inline	const double val_cent_PbPb[11]={0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};
inline  const double val_cent_PbPb_wide[7]={0.05,0.1,0.2,0.4,0.6,0.8,1.0};
inline	const std::string name_cent_pp[10]={"0-0.95pct","0.95-4.7pct", "4.7-9.5pct", "9.5-14pct","14-19pct","19-28pct","28-38pct", "38-48pct", "48-68pct", "68-100pct"};
inline	const std::string name_NtrkClass_pp[13]={"0-10Ntrk", "10-20Ntrk", "20-30Ntrk", "30-40Ntrk", "40-60Ntrk", "60-85Ntrk", "85-95Ntrk", "95-105Ntrk", "105-115Ntrk", "115-125Ntrk", "125-135Ntrk", "135-150Ntrk", "150-infNtrk"};
inline	const std::string name_cent_original_narrow[10]={"0-0.25pct","0.25-1pct", "1-5pct", "5-10pct","10-20pct","20-30pct","30-40pct", "40-50pct", "50-70pct", "70-100pct"};
inline	const std::string name_cent_original[4]={"0-10pct","10-30pct", "30-50pct", "50-100pct"};
inline	const std::string name_cent_pPb[7]={"0-5pct", "5-10pct", "10-20pct", "20-40pct", "40-60pct", "60-80pct", "80-100pct"};
inline	const std::string name_cent_PbPb[11]={"0-5pct","5-10pct", "10-20pct","20-30pct", "30-40pct","40-50pct","50-60pct","60-70pct","70-80pct","80-90pct","90-100pct"};
inline  const std::string name_cent_PbPb_wide[7]={"0-5pct","5-10pct", "10-20pct","20-40pct","40-60pct","60-80pct","80-100pct"};
inline	const double RtBins[5]={0.5, 1.0, 2.0, 3.0, LARGE};
	//const double RtBins[9]={0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.5, 3.5, LARGE};
inline	const double ptMax = 50.0;
inline	const double etaMax = 10.0;
inline	const double yMax = 10.0;

inline	const std::string core_tag="H";
inline	const std::string corona_tag="S";
inline	int num_of_Species_Rt = 7;


inline	const int margin = 10;
inline  const int N_RndomEv=20;
inline  const int Nsample=5;

inline	const std::string data_directory="data";
inline	const double binSize_small=d_x_HI;
inline	const double binSize_log=2.0;
inline	const double switchBin_x = 0.0;

inline	const double minNchHI=20.0;
inline	const double maxNchPP=60.0;

inline	const string default_inputfname = "input";
inline	const string default_out_fname="hist.txt";
inline  const string seed_outputfname="SEED.txt";
inline	const int default_nfiles=10000;
inline	const vector<string> save_settings_fname={"include/Constants.h"};


inline	const std::string settings_outputfname="settings_yuuka_analysis.txt";
inline	const std::string default_subdir_name="centrality_cut";
inline	const std::string fname_eByeInfo="eByeInfo.txt";

inline	const double dummy = -100.0;
inline  const int PrintCounter=1000;
inline  const int PrintCounterTL=100;

inline	std::complex<double> initialval_comp (0.0, 0.0);

	//DEFINITION OF i
	//----------------
inline	const   std::complex<double> i_img(0.0,1.0);  
inline	const double WARNING_IMAGINARY_MAX = 100.0;



}
#endif
