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


bool read(const string& fname, vector <double> &vec1, vector<double> &vec2){

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
	ct.Nev_tot+=1;
      }else{
	istringstream is(templine);
	int event_num;
	double Nch, weight;
	is >> event_num >> Nch >> weight;
        vec1.push_back(Nch);
        vec2.push_back(weight);
      }
    }
    in.close();
  }
  return true;
}


void analysis(const vector <double> &vec1, vector<double> &vec2){


  for(int i=0; i<(int)vec1.size(); ++i){
    double x_val=vec1[i];

    //count.
    //-----------------------------------
    if(x_val<constants::x_min || x_val>constants::x_max) continue;
    int nx=(int)((x_val/constants::d_x)+(fabs(constants::x_min)/constants::d_x));
    ct.Hist[nx]+=vec2[i];
    if(ct.max_nx<nx) ct.max_nx=nx;

    ct.sum_weight+=vec2[i];
  
  }


    //take average and devide by cell width
    //-------------------------------------
    for(int i=0; i<ct.max_nx+1; ++i){
      //ct.Hist[i]/=ct.Nev_tot;
      ct.Hist[i]/=ct.sum_weight;
      ct.Hist[i]/=constants::d_x;
    }
 
}
  
bool write(const string& fname){
  ofstream ofs;
  ofs.open(fname.c_str());
  if(!ofs){ms->open(fname); return false;}

  ct.max_nx+=constants::margin;
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
    if(!strcmp(argv[i],"-input_path")) {inputfname= argv[i+1];}
    if(!strcmp(argv[i],"--ext")) {ext= argv[i+1];}
  }
      
  vector <double> vec1, vec2;
  for(int i=0; i<nfiles; ++i){
    if(!(i%1000)) ms->read(i);
    //if(!read(inputfname+uf.generateS(i)+ext, vec1, vec2)) return 1;
    if(!read(inputfname, vec1, vec2)) return 1;
  }
  
  analysis(vec1, vec2);
    
  uf.make_output_directory(out_directory_name);
  string generated_directory_name=uf.get_name_directory();
  if(!write(generated_directory_name+"/"+out_fname+ext)) return 1;

  ms->finish();
  return 0;
}
