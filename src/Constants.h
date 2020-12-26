#ifndef CONSTANTS
#define CONSTANTS
#include <sstream>
#include <vector>

using std::string;
using std::vector;

namespace constants{
  
  const int margin = 10;
  const std::string data_directory="data";

  const string default_inputfname = "input";
  const string default_out_fname="hist.txt";
  const int default_nfiles=10000;
  const string settings_outputfname="settings.txt";
  const vector<string> save_settings_fname={"src/Constants.h"};


#ifdef DNDETA_PROTON
  const std::string MODE = "dndeta_proton";
  const string default_out_directory_name="dndeta_proton";
  const double x_max = 20.0;
  const double x_min = -20.0;
  const double d_x = 0.5;
  const int x_cell_capa=1000;
  const string default_ext="/hadronised.txt";
#elif defined VERTICES
  const std::string MODE = "vertices";
  const string default_out_directory_name="vertices";
  const double x_max = 20.0;
  const double x_min = -20.0;
  const double d_x = 0.5;
  const int x_cell_capa=1000;
  const string default_ext="/hadronised.txt";
#elif defined DNDPT
  const std::string MODE = "dndpt";
  const string default_out_directory_name="dndpt";
  const double x_max = 2000.0;
  const double x_min = 0.0;
  const double d_x = 1.0;
  const int x_cell_capa=1000;
  const string default_ext="/hadronised.txt";
#elif defined JET_PRAC
  const std::string MODE = "JET_PRAC";
  const string default_out_directory_name="JET_PRAC";
  const double x_max = 10.0;
  const double x_min = -10.0;
  const double d_x = 0.5;
  const int x_cell_capa=1000;
  const string default_ext="/jetinfo.txt";
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
 
 const double delta_phi_SMALL = 0.5;
 const double delta_eta = 0.5;


}
#endif
