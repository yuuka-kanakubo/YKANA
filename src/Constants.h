#ifndef CONSTANTS
#define CONSTANTS
#include <sstream>
#include <vector>

using std::string;
using std::vector;

namespace constants{
  
  const double x_max = 20.0;
  const double x_min = -20.0;
  const double d_x = 0.5;
  const int x_cell_capa=1000;
  const int margin = 10;
  const std::string data_directory="data";

  const string default_inputfname = "input";
  const string default_ext="/hadronised.txt";
  const string default_out_fname="hist.txt";
  const int default_nfiles=10000;
  const string settings_outputfname="settings.txt";
  const vector<string> save_settings_fname={"src/Constants.h"};


#ifdef WORK5
  const std::string MODE = "WORK5";
  const string default_out_directory_name="WORK5";
#elif defined WORK7
  const std::string MODE = "WORK7";
  const string default_out_directory_name="WORK7";
#elif defined WORK3
  const std::string MODE = "WORK3";
  const string default_out_directory_name="output";
#elif defined WORK4
  const std::string MODE = "WORK4";
  const string default_out_directory_name="output";
#elif defined WORK5
  const std::string MODE = "WORK5";
  const string default_out_directory_name="output";
#else
  const std::string MODE = "default";
  const string default_out_directory_name="WORK4";
#endif





  const int id_nucleus=1000000000;
  const int id_gluon=21;
  const int id_ch_pion=211;
  const int id_ch_kaon=321;
  const int id_proton=2212;

 const double w_eta_multiplicity_INEL_lg_0=1.0;
 const double w_eta_ATLAS_cut=2.5;
 const double V0M_fwd1=2.8;
 const double V0M_fwd2=5.1;
 const double V0M_bkw1=-3.7;
 const double V0M_bkw2=-1.7;
 const double outer_SPD=1.4;
 


}
#endif
