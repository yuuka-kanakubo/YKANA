// -*- mode:c++ -*-
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <sstream>
#include <fstream>
#include "constants.h"
#include "util_func.h"
#include "container.h"
#include "message.h"

using namespace std;

util_func uf;
container ct;
message* ms;  


bool read(const string& fname, vector <double> &vec){

  ifstream in;
  in.open(fname.c_str(),ios::in);
  if(!in){ ms->open(fname); return false;}
  
  {
    string templine;
    while(getline(in,templine)) {
      if(templine.find('#')!=std::string::npos) {
      } else if(templine.find('%')!=std::string::npos){
	istringstream iss(templine);
	string pct;
	int event_num;
	int num_jet;
	iss >> pct >> event_num >> num_jet;
	ct.number_of_jet+=num_jet;
      }else{
	istringstream is(templine);
	int jet_id;
	double e, px, py, pz;
	is >> jet_id >> e >> px >> py >> pz;
        vec.push_back(e);
      }
    }
    in.close();
  }
  return true;
}


void analysis(const vector <double> &vec){
  for(int i=0; i<(int)vec.size(); ++i){
    double x_val=vec[i];

    //count.
    //-----------------------------------
    if(x_val<constants::x_min || x_val>constants::x_max) continue;
    int nx=(int) (x_val/constants::d_x)+(fabs(constants::x_min)/constants::d_x);
    ct.Hist[nx]++;
    if(ct.max_nx<nx) ct.max_nx=nx;
  
  }

    //take average and devide by cell width
    //-------------------------------------
    for(int i=0; i<ct.max_nx; ++i){
      ct.Hist[i]/=ct.number_of_jet;
      ct.Hist[i]/=constants::d_x;
    }
 
}
  
bool write(const string& fname){
  ofstream ofs;
  ofs.open(fname.c_str());
  if(!ofs){ms->open(fname); return false;}
  for(int i=0; i<ct.max_nx; ++i){
    double x_axis = ((constants::x_min+(constants::d_x*i))+(constants::x_min+(constants::d_x*(i+1))))/2.0;
    ofs << setw(16) << fixed << setprecision(8) << x_axis
	<< setw(16) << ct.Hist[i] << endl;
  }
  ofs << endl;
  return true;
}
 
 
int main(int argc, char* argv[]){

  if(!ms->enough_argument(argc)) return 1;

  string inputfname = constants::default_inputfname;
  string out_directory_name=constants::default_out_directory_name;
  string ext=constants::default_ext;
  string out_fname=constants::default_out_fname;
  int nfiles=constants::default_nfiles;
  for(int i=1; i<argc; i++) {
    if(!strcmp(argv[i],"-n")) {nfiles = atoi(argv[i+1]);}
    if(!strcmp(argv[i],"-output_dir")){out_directory_name= argv[i+1];}
    if(!strcmp(argv[i],"--input_path")) {inputfname= argv[i+1];}
    if(!strcmp(argv[i],"--ext")) {ext= argv[i+1];}
  }
      
  vector <double> vec;
  for(int i=0; i<nfiles; ++i){
    if(!read(inputfname+uf.generateS(i)+ext, vec)) return 1;
  }
  
  analysis(vec);
    
  uf.make_output_direcrory(out_directory_name);
  string generated_directory_name=uf.get_name_directory();
  if(!write(generated_directory_name+"/"+out_fname+ext)) return 1;

  ms->finish();
  return 0;
}
