// -*- mode:c++ -*-
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <sstream>
#include <fstream>
#include "Constants.h"
#include "Util_func.h"
#include "Container.h"
#include "Message.h"
#include "LogSettings.h"
#include "Settings.h"

using namespace std;


class Analysis{


Util_func uf;
Container ct;
Message* ms;  
const Settings::Options options;
LogSettings set;

public:

 Analysis(const Settings::Options options_in, LogSettings set_in): options(options_in), set(set_in){
	 this->ana();
         ms = new Message();
 };
 ~Analysis(){};



 class ParticleInfo{
	 public:

   double e;
   double px;
   double py;
   double pz;
   double phi;

 };


bool read(const std::string& fname, vector <double> &vec1, bool& INEL_lg_0){

  ifstream in;
  in.open(fname.c_str(),ios::in);
  if(!in){ ms->open(fname); return false;}



bool Multiplicity_INEL_lg_0=false;
  
  {
    std::string templine;
    while(getline(in,templine)) {
      if(templine.find('#')!=std::string::npos) {
      } else if(templine.find('%')!=std::string::npos){
	istringstream iss(templine);
	std::string pct;
	int event_num;
	int num_jet;
	iss >> pct >> event_num >> num_jet;
      }else{
	istringstream is(templine);
	int data1, data2, ID, col, acol;
	double m,e,px,py,pz,x,y,z,t,ft, rap;
	std::string TAG;
	is >> data1 >> data2 >> col >> acol >> ID >> m >> e >> px >> py >> pz >> rap >> x >> y >> z >> t >> ft >> TAG;


		double P_squared=px*px+py*py+pz*pz;
		double P=(P_squared)>0.0 ? sqrt(P_squared):0.0;
		double pt_squared=px*px+py*py;
		double pt=(pt_squared)>0.0 ? sqrt(pt_squared):0.0;
	
		////
		double eta;
		{
			const double LARGE = 10.0;
			if(P==0.0 && pz==0.0) eta=0.0;
			else if(P==-pz) eta=-LARGE;
			else if(P==pz) eta=LARGE;
			else eta=log((P + pz)/(P - pz)) / 2.0;
		}
		////



		if(constants::MODE.find("WORK5")!=string::npos){
			if(ID==constants::id_proton) { 

				vec1.push_back(eta);
			}
		}else if(constants::MODE.find("WORK8")!=string::npos){
			if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon) && fabs(eta)<constants::delta_eta ) { 

				vec1.push_back(pt);

                        }

		}else{
			if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon) ) { 

				vec1.push_back(eta);
			}
		}


	//INEL > 0?
	//-----------
	if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon) && fabs(eta)<constants::w_eta_multiplicity_INEL_lg_0) { 
		Multiplicity_INEL_lg_0=true;
	}
	


      }
    }
    in.close();
  }
  if(Multiplicity_INEL_lg_0) {
	  INEL_lg_0=true;
	  ct.Nev_tot+=1;
  }

  return true;
}




bool read_jetinfo(const std::string& fname, vector <double> &vec1, bool& INEL_lg_0){



//Not consider INEL in jet analysis.
	INEL_lg_0=true;


  ifstream in;
  in.open(fname.c_str(),ios::in);
  if(!in){ ms->open(fname); return false;}

  vector <ParticleInfo> Parents, Jets;
  
  {
    std::string templine;
    while(getline(in,templine)) {
      if(templine.find('G')!=std::string::npos) {
      } else if(templine.find('P')!=std::string::npos){
	istringstream iss(templine);
	double e,px,py,pz;
	std::string TAG;
	iss >> TAG >> e >> px >> py >> pz;
	ParticleInfo info;
	info.e=e;
	info.px=px;
	info.py=py;
	info.pz=pz;
	info.phi=atan2(py, px);
        Parents.push_back(info);
      }else{
	istringstream is(templine);
	double e,px,py,pz;
	int num;
	is >> num >> e >> px >> py >> pz;
	ParticleInfo info;
	info.e=e;
	info.px=px;
	info.py=py;
	info.pz=pz;
	info.phi=atan2(py, px);
        Jets.push_back(info);
      }
    }
    in.close();
  }


 for(int i =0; i<(int)Jets.size(); i++){
	 double phi_J = Jets[i].phi;
	 for(int j =0; j<(int)Parents.size(); j++){
             if(fabs(phi_J-Parents[j].phi)<constants::delta_phi_SMALL){

               double frac = (Jets[i].e - Parents[j].e)/Parents[j].e;
               vec1.push_back(frac);
	       ct.Nev_tot+=1;

             }
  
	 }
 }

  return true;
}


void stat(){

    //take average    
    //-------------------------------------
    for(int i=0; i<ct.max_nx+1; ++i){
      ct.Hist[i]/=ct.Nev_tot;
      //ct.Hist[i]/=ct.sum_weight;
      ct.HistHist[i]/=ct.Nev_tot;
      //ct.HistHist[i]/=ct.sum_weight;

    }

    cout << "Event: " << ct.Nev_tot << endl;

    //Get standard error
    //-------------------------------------
    for(int i=0; i<ct.max_nx+1; ++i){
      double var=ct.HistHist[i]-pow(ct.Hist[i],2.0);
      ct.HistErr[i]=sqrt(var/ct.Nev_tot);
    }


 
}





void fill(const vector <double> &vec1, const bool INEL_lg_0){


if(!INEL_lg_0) return;


  for(int i=0; i<(int)vec1.size(); ++i){
    double x_val=vec1[i];

    //count.
    //-----------------------------------
    if(x_val<constants::x_min || x_val>constants::x_max) continue;
    int nx=(int)((x_val/constants::d_x)+(fabs(constants::x_min)/constants::d_x));
    ct.Hist[nx]+=(1.0/constants::d_x);
    ct.Hist_1ev[nx]+=(1.0/constants::d_x);
    if(ct.max_nx<nx) ct.max_nx=nx;

  }


    for(int i=0; i<ct.max_nx+1; ++i){
      ct.HistHist[i]+=pow(ct.Hist_1ev[i],2.0);
      
     //Clear.
     //---------
     ct.Hist_1ev[i]=0.0;


    }




}




  
bool write(const std::string& fname){
  ofstream ofs;
  ofs.open(fname.c_str());
  if(!ofs){ms->open(fname); return false;}

  ct.max_nx+=constants::margin;
  for(int i=0; i<ct.max_nx; ++i){
    double x_axis = ((constants::x_min+(constants::d_x*i))+(constants::x_min+(constants::d_x*(i+1))))/2.0;
    ofs << setw(16) << fixed << setprecision(8) << x_axis << "  "
	<< setw(16) << ct.Hist[i] << "  "
	<< setw(16) << ct.HistErr[i] << endl;
  }
  ofs << endl;
  return true;
}
 
 
int ana(){
      
  for(int i=0; i<options.nfiles; ++i){
    if(!(i%1000)) ms->read(i);

    std::stringstream ss;
    ss << options.inputfname << "/ev" << setw(9) << setfill('0') << i << options.ext;
    vector <double> vec1;
    bool INEL_lg_0=false;
    if(constants::MODE.find("WORK9")!=std::string::npos){
	    read_jetinfo(ss.str(), vec1, INEL_lg_0);
    }else{
	    read(ss.str(), vec1, INEL_lg_0);
    }
    this->fill(vec1, INEL_lg_0);
  }
  
  this->stat();
    
  uf.make_output_directory(options.out_directory_name);
  std::string generated_directory_name=uf.get_name_directory();
  if(!write(generated_directory_name+"/"+options.out_fname)) return 1;

  set.archive_settings(generated_directory_name);

  return 0;
}





};
