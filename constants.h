#ifndef CONSTANTS
#define CONSTANTS
#include <sstream>

using std::string;

namespace constants{
  
  const double x_max = 2.0;
  const double x_min = 0.0;
  const double d_x = 0.01;
  const int x_cell_capa=1500;
  const std::string data_directory="data";  

  const string default_inputfname = "input";
  const string default_out_directory_name="output";
  const string default_ext=".txt";
  const string default_out_fname="data";
  const int default_nfiles=0;

}
#endif
