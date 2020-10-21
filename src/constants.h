#ifndef CONSTANTS
#define CONSTANTS
#include <sstream>
#include <vector>

using std::string;
using std::vector;

namespace constants{
  
  const double x_max = 100.0;
  const double x_min = 0.0;
  const double d_x = 5.0;
  const int x_cell_capa=1000;
  const int margin = 10;
  const std::string data_directory="data";

  const string default_inputfname = "input";
  const string default_out_directory_name="output";
  const string default_ext=".dat";
  const string default_out_fname="hist";
  const int default_nfiles=0;
  const string settings_outputfname="settings.txt";
  const vector<string> save_settings_fname={"src/constants.h"};

}
#endif
