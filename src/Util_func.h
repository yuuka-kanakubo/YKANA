#ifndef UTIL_FUNC
#define UTIL_FUNC
#include <time.h>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <iomanip>
#include "Constants.h"
#include "EbyeInfo.h"

class Util_func{

 private:

    bool isQuark(const int pid) {
      if(abs(pid)<9||this->is_diquark(pid)) return true;
      else return false;
    }


    bool is_diquark(const int pid){
      bool is_diquark=false;
      if(abs(pid)!=2212 && abs(pid)!=2112 && abs(pid)<10000 && 999<abs(pid)) {
	int PID= pid/10;
	if(PID%10==0) is_diquark=true;
      }
      return is_diquark;
    }

  void get_EbyeInfo(const string fname, EbyeInfo& ebye, double rap_shift, bool VZEROAND_trigger, bool parton_level){
    ifstream in;
    in.open(fname.c_str(),ios::in);
    if(!in) {
      cerr << "Error unable to open file " << fname << endl;
      return ;
    } 
    string templine;
    double Multiplicity=0.0;
    double Multiplicity_V0M=0.0;
    double Multiplicity_V0A=0.0;
    bool Multiplicity_INEL_lg_0=false;
    bool V0M_FWD=false;
    bool V0M_BKW=false;
    bool OUTER_SPD=false;
    bool ATLAS_cut=false;
    double Weight=0.0;
    while(getline(in,templine)) {
      if(templine.find('#')!=std::string::npos) {
      }else if(templine.find('%')!=std::string::npos){
	      istringstream is(templine);
	      int iev, nv;
	      string com;
	      double weight_in, tau;
	      is >> com >> iev >> nv >> tau >> weight_in;
              Weight = weight_in; 
      }else{
	istringstream is(templine);
	int ID, data1, data2, col, acol;
	double mass, energy, px, py, pz, x,y,z,t,rap;
	is >> data1 >> data2 >> col >> acol >> ID >> mass >> energy >> px >> py >> pz >> rap >> x >> y >> z >> t;
	if(!parton_level && ID==constants::id_gluon){cout << ":o WARNING: This is partonEvent." << endl; return ;}


	if(rap_shift!=0.0){
		double m_T=energy/cosh(rap);
		double rap__= rap-rap_shift;
		energy=m_T*cosh(rap__);
		pz=m_T*sinh(rap__);
	}


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
	
	if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon) && fabs(eta)<constants::w_eta_multiplicity) { 
		Multiplicity++;
	}
	if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon) && fabs(eta)<constants::w_eta_multiplicity_INEL_lg_0) { 
		Multiplicity_INEL_lg_0=true;
	}
	if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon) && this->in_V0M_detector(eta)) { 
		Multiplicity_V0M++;
	}
	if(abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon){
            if(this->in_V0M_fwd(eta)) { V0M_FWD=true; Multiplicity_V0A++;} 
            if(this->in_V0M_bkw(eta)) { V0M_BKW=true;}
            if(this->in_outerSPD(eta)){ OUTER_SPD=true;}
	}

	if(parton_level){
		if((ID==constants::id_gluon || this->isQuark(ID)) && fabs(eta)<constants::w_eta_multiplicity) { 
			Multiplicity++;
		}
		if((ID==constants::id_gluon || this->isQuark(ID)) && fabs(eta)<constants::w_eta_multiplicity_INEL_lg_0) { 
			Multiplicity_INEL_lg_0=true;
		}
		if((ID==constants::id_gluon || this->isQuark(ID)) && this->in_V0M_detector(eta)) { 
			Multiplicity_V0M++;
		}
		if(ID==constants::id_gluon || this->isQuark(ID)){
			if(this->in_V0M_fwd(eta)) { V0M_FWD=true; Multiplicity_V0A++;} 
			if(this->in_V0M_bkw(eta)) { V0M_BKW=true;}
			if(this->in_outerSPD(eta)){ OUTER_SPD=true;}
		}
	}

	if(abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon){
		if(pt>0.5 && fabs(eta)<constants::w_eta_ATLAS_cut) ATLAS_cut=true;
	}

      }
    }
    in.close();
    ebye.multiplicity=Multiplicity;
    ebye.weight=Weight;
    ebye.multiplicity_INEL_lg_0=Multiplicity_INEL_lg_0;
    ebye.multiplicity_V0M=(VZEROAND_trigger)? Multiplicity_V0A : Multiplicity_V0M;
    ebye.valid=true;
    if(V0M_FWD && V0M_BKW && OUTER_SPD) {ebye.trig_3outof3=true;}
    if(V0M_FWD && V0M_BKW) {ebye.trig_VZEROAND=true;}
    if((V0M_FWD && V0M_BKW) || (V0M_BKW && OUTER_SPD) || (OUTER_SPD && V0M_FWD) ) {ebye.trig_2outof3=true;}
    if(ATLAS_cut) ebye.ATLAS_cut=true;
  }

  bool in_V0M_detector(const double eta){

	  if(this->in_V0M_fwd(eta) || this->in_V0M_bkw(eta)) return true;
	  else return false;

  }


  bool in_V0M_fwd(const double eta){
	  if(eta<constants::V0M_fwd2 && eta>constants::V0M_fwd1) return true;
          else return false;
  }

 bool in_V0M_bkw(const double eta){
	 if(eta>constants::V0M_bkw1 && eta<constants::V0M_bkw2) return true;
	 else return false;
 }
 
 bool in_outerSPD(const double eta){
    if(fabs(eta)<constants::outer_SPD) return true;
    else return false; 
 }
 


 struct stat st;
 std::string name_directory;

 const std::string generateS2(int n){
   std::ostringstream name;
   name  << std::setw(2) << std::setfill('0') << n ;
   return name.str();
 };
 
  const std::string get_date(){
    time_t now = time(NULL);
    struct tm *pnow = localtime(&now);  
    int Y=pnow->tm_year+1900;
    int M=pnow->tm_mon + 1;
    int D=pnow->tm_mday;
    std::string date=generateS(Y)+generateS2(M)+generateS2(D);    
    return date;
  };




  const void make_data_directory(){
    if(stat(constants::data_directory.c_str(),&st) !=0) mkdir(constants::data_directory.c_str(),0775);
    else {};
  };




 public:
  Util_func(){};
  ~Util_func(){};
  
  std::string generateS(double n){
    std::ostringstream oss;
    oss <<  n ; 
    return oss.str();
  };
 
  std::string get_name_directory(){
    return name_directory;
  };





  void get_EbyeInfo_(const string fname, EbyeInfo& ebye, double rap_shift, bool VZEROAND_trigger, bool parton_level){
    this->get_EbyeInfo(fname, ebye, rap_shift, VZEROAND_trigger, parton_level);
  }



 
  
  const void make_output_directory(const std::string& directory_name){
    
    this->make_data_directory();
    
    std::ostringstream name; 
    name  << this->get_date()+"_"+directory_name;
    name_directory=constants::data_directory+"/"+name.str();
    if(stat(name_directory.c_str(),&st)!=0) mkdir(name_directory.c_str(),0775);
    else {};

  };


};
#endif
