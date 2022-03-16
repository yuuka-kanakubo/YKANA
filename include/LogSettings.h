#ifndef LOGSETTINGS
#define LOGSETTINGS
#include <iostream>
#include <fstream>
#include "Constants.h"

using namespace std;

class LogSettings{


 private:

	 vector<string> BinSettings_str;
	 bool centrality_cut;




 public:

	 vector<string> options;

	 int get_BinSettings_size(){return (int)BinSettings_str.size();}
	 void set_centrality_cut(bool flag){centrality_cut=flag;}


void save_BinSettings(const vector<double> xMin, const vector<double> xMax){


	    {
		    BinSettings_str.push_back("//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-");
		    BinSettings_str.push_back("//Bins..");
		    for(int i=0; i<(int)xMin.size(); i++){
			    ostringstream str;
			    str << xMin[i] << "  " << xMax[i];
			    BinSettings_str.push_back(str.str());
		    }
		    BinSettings_str.push_back("");
	    }

cout << "Saved bin settings " << (int)BinSettings_str.size() << endl;

}



bool archive_settings(const string output_directory_name){


	{
		ifstream ifs(output_directory_name+"/"+constants::settings_outputfname);
		if(ifs.is_open() && this->centrality_cut) {
			return true;
		}
	} 


    vector<string> templine_set;


    templine_set.push_back("//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-");
    templine_set.push_back("");
    for(int i=0; i<(int)options.size(); ++i){
	    templine_set.push_back(options[i]);
    }

    templine_set.push_back("");
    templine_set.push_back("//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-");
    templine_set.push_back("//Constants");
    templine_set.push_back("");

    this->get_constants(templine_set);


    //Output
    //--------
    ofstream ofs;
    ofs.open(output_directory_name+"/"+constants::settings_outputfname,ios::in|ios::app);
    for(int i=0; i<(int)templine_set.size(); ++i)
	    ofs << templine_set[i] << endl;
    cout << __FILE__ << "log " << get_BinSettings_size() << endl;
    for(int i=0; i<(int)BinSettings_str.size(); ++i){
	    ofs << BinSettings_str[i] << endl;
    }

    ofs.close();
    return true;
  }


private:

void get_constants(vector<string> &temp){

//All constants..
//=================
ostringstream oss;
oss <<
"constants::twopc1Dptmin = " << constants::twopc1Dptmin << endl << 
"constants::twopc1Dptmax = " << constants::twopc1Dptmax << endl <<
"constants::twopc1DetaRange = " << constants::twopc1DetaRange << endl <<
"constants::twopc1Dptmin_N = " << constants::twopc1Dptmin_N << endl <<
"constants::twopc1Dptmax_N = " << constants::twopc1Dptmax_N << endl <<
"constants::twopc1DetaRange_N = " << constants::twopc1DetaRange_N << endl <<
"constants::twopc1D_Nmax = " << constants::twopc1D_Nmax << endl <<
"constants::twopc1D_Nmin = " << constants::twopc1D_Nmin << endl <<
"constants::twopc2DetaRange = " << constants::twopc2DetaRange << endl <<
"constants::twopc2Dptmin_N = " << constants::twopc2Dptmin_N << endl <<
"constants::twopc2Dptmax_N = " << constants::twopc2Dptmax_N << endl <<
"constants::twopc2DetaRange_N = " << constants::twopc2DetaRange_N << endl <<
"constants::trig_ptmin = " << constants::trig_ptmin << endl <<
"constants::trig_ptmax = " << constants::trig_ptmax << endl <<
"constants::assoc_ptmin = " << constants::assoc_ptmin << endl <<
"constants::assoc_ptmax = " << constants::assoc_ptmax << endl <<
"constants::minpt_Rt = " << constants::minpt_Rt << endl <<
"constants::maxPhi_RtTrans = " << constants::maxPhi_RtTrans << endl <<
"constants::minPhi_RtTrans = " << constants::minPhi_RtTrans << endl <<
"constants::etaRange_Rt = " << constants::etaRange_Rt << endl <<
"constants::etaRangeMultiplicity = " << constants::etaRangeMultiplicity << endl <<
"constants::etaA_3sub  = " << constants::etaA_3sub  << endl <<
"constants::etaB_3sub  = " << constants::etaB_3sub  << endl <<
"constants::DeltaEtaFULL = " << constants::DeltaEtaFULL << endl <<
"constants::DeltaEtaNS = " << constants::DeltaEtaNS << endl <<
"constants::DeltaPhiOUT = " << constants::DeltaPhiOUT << endl <<
"constants::DeltaPhiNS = " << constants::DeltaPhiNS << endl <<
"constants::etaGap = " << constants::etaGap << endl <<
"constants::ptmin_cumulantmulti = " << constants::ptmin_cumulantmulti << endl <<
"constants::ptmax_cumulantmulti = " << constants::ptmax_cumulantmulti << endl <<
"constants::ptmin_cumulantpt = " << constants::ptmin_cumulantpt << endl <<
"constants::ptmax_cumulantpt = " << constants::ptmax_cumulantpt << endl <<
"constants::ptmin_cumulanteta = " << constants::ptmin_cumulanteta << endl <<
"constants::ptmax_cumulanteta = " << constants::ptmax_cumulanteta << endl <<
"constants::ptmin_cumulantmulti_Nch = " << constants::ptmin_cumulantmulti_Nch << endl <<
"constants::ptmax_cumulantmulti_Nch = " << constants::ptmax_cumulantmulti_Nch << endl <<
"constants::etaRange_cumulantmulti = " << constants::etaRange_cumulantmulti << endl <<
"constants::etaRange_cumulantpt = " << constants::etaRange_cumulantpt << endl <<
"constants::etaRange_cumulantmulti_Nch = " << constants::etaRange_cumulantmulti_Nch << endl <<
"constants::hbarc  = " << constants::hbarc  << endl <<
"constants::TL_tau_00  = " << constants::TL_tau_00  << endl <<
"constants::TL_tau_0  = " << constants::TL_tau_0  << endl <<
"constants::TL_tau_switch  = " << constants::TL_tau_switch  << endl <<
"constants::TL_tau_switchFM = " << constants::TL_tau_switchFM << endl <<
"constants::TL_dtau1 = " << constants::TL_dtau1 << endl <<
"constants::TL_dtau1FM = " << constants::TL_dtau1FM << endl <<
"constants::TL_dtau2 = " << constants::TL_dtau2 << endl <<
"constants::ext_nameTLxy = " << constants::ext_nameTLxy << endl <<
"constants::ext_nameTLxeta = " << constants::ext_nameTLxeta << endl <<
"constants::MODE  = " << constants::MODE  << endl <<
"constants::default_out_directory_name = " << constants::default_out_directory_name << endl <<
"constants::x_max  = " << constants::x_max  << endl <<
"constants::x_min  = " << constants::x_min  << endl <<
"constants::d_x  = " << constants::d_x  << endl <<
"constants::x_cell_capa = " << constants::x_cell_capa << endl <<
"constants::default_ext = " << constants::default_ext << endl <<
"constants::y_max  = " << constants::y_max  << endl <<
"constants::y_max_HI  = " << constants::y_max_HI  << endl <<
"constants::y_min  = " << constants::y_min  << endl <<
"constants::d_y  = " << constants::d_y  << endl <<
"constants::d_y_HI  = " << constants::d_y_HI  << endl <<
"constants::y_cell_capa = " << constants::y_cell_capa << endl <<
"constants::id_nucleus = " << constants::id_nucleus << endl <<
"constants::id_gluon = " << constants::id_gluon << endl <<
"constants::id_ch_pion = " << constants::id_ch_pion << endl <<
"constants::id_ch_kaon = " << constants::id_ch_kaon << endl <<
"constants::id_proton = " << constants::id_proton << endl <<
"constants::id_lambda = " << constants::id_lambda << endl <<
"constants::id_cascade = " << constants::id_cascade << endl <<
"constants::id_omega = " << constants::id_omega << endl <<
"constants::id_phi = " << constants::id_phi << endl <<
"constants::id_K0S = " << constants::id_K0S << endl <<
"constants::w_eta_multiplicity = " << constants::w_eta_multiplicity << endl <<
"constants::w_eta_multiplicity_INEL_lg_0 = " << constants::w_eta_multiplicity_INEL_lg_0 << endl <<
"constants::w_eta_ATLAS_cut = " << constants::w_eta_ATLAS_cut << endl <<
"constants::V0M_fwd1 = " << constants::V0M_fwd1 << endl <<
"constants::V0M_fwd2 = " << constants::V0M_fwd2 << endl <<
"constants::V0M_bkw1 = " << constants::V0M_bkw1 << endl <<
"constants::V0M_bkw2 = " << constants::V0M_bkw2 << endl <<
"constants::outer_SPD = " << constants::outer_SPD << endl <<
"constants::etaRangeCMSRidge = " << constants::etaRangeCMSRidge << endl <<
"constants::MomentumMinCMSRidge = " << constants::MomentumMinCMSRidge << endl <<
"constants::delta_phi_SMALL  = " << constants::delta_phi_SMALL  << endl <<
"constants::delta_eta  = " << constants::delta_eta  << endl <<
"constants::pPb_mid_rapidity__fwd = " << constants::pPb_mid_rapidity__fwd << endl <<
"constants::pPb_mid_rapidity__bkw = " << constants::pPb_mid_rapidity__bkw << endl <<
"constants::pPb_rap_shift_from_lab_to_cm = " << constants::pPb_rap_shift_from_lab_to_cm << endl <<
"constants::default_midy_pm  = " << constants::default_midy_pm  << endl <<
"constants::LARGE = " << constants::LARGE << endl <<
"constants::LARGEint = " << constants::LARGEint << endl <<
"constants::SMALL = " << constants::SMALL << endl <<
"constants::multip_cut_more_than = " << constants::multip_cut_more_than << endl <<
"constants::multip_cut_less_than = " << constants::multip_cut_less_than << endl <<
"constants::ptMax  = " << constants::ptMax  << endl <<
"constants::etaMax  = " << constants::etaMax  << endl <<
"constants::yMax  = " << constants::yMax  << endl <<
"constants::core_tag = " << constants::core_tag << endl <<
"constants::corona_tag = " << constants::corona_tag << endl <<
"constants::num_of_Species_Rt  = " << constants::num_of_Species_Rt  << endl <<
"constants::margin  = " << constants::margin  << endl <<
"constants::N_RndomEv = " << constants::N_RndomEv << endl <<
"constants::Nsample = " << constants::Nsample << endl <<
"constants::data_directory = " << constants::data_directory << endl <<
"constants::binSize_small = " << constants::binSize_small << endl <<
"constants::binSize_log = " << constants::binSize_log << endl <<
"constants::switchBin_x  = " << constants::switchBin_x  << endl <<
"constants::minNchHI = " << constants::minNchHI << endl <<
"constants::maxNchPP = " << constants::maxNchPP << endl <<
"constants::default_inputfname  = " << constants::default_inputfname  << endl <<
"constants::default_out_fname = " << constants::default_out_fname << endl <<
"constants::seed_outputfname = " << constants::seed_outputfname << endl <<
"constants::default_nfiles = " << constants::default_nfiles << endl <<
"constants::settings_outputfname = " << constants::settings_outputfname << endl <<
"constants::default_subdir_name = " << constants::default_subdir_name << endl <<
"constants::fname_eByeInfo = " << constants::fname_eByeInfo << endl <<
"constants::dummy  = " << constants::dummy  << endl <<
"constants::PrintCounter = " << constants::PrintCounter << endl <<
"constants::PrintCounterTL = " << constants::PrintCounterTL << endl <<
"constants::WARNING_IMAGINARY_MAX = " << constants::WARNING_IMAGINARY_MAX << endl;

temp.push_back(oss.str());
}

};
#endif
