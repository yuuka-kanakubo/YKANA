#ifndef CONSTANTS
#define CONSTANTS
#include <sstream>

using std::string;

namespace constants{
  
  const double x_max = 5000.0;
  const double x_min = 0.0;
  const double d_x = 500.0;
  const int x_cell_capa=1000;
  const int margin = 10;
  const std::string data_directory="data";

  const string default_inputfname = "input";
  const string default_out_directory_name="output";
  const string default_ext=".dat";
  const string default_out_fname="hist";
  const int default_nfiles=0;

}
#endif
