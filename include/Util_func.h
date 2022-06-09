#ifndef UTIL_FUNC
#define UTIL_FUNC
#include <time.h>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <iomanip>
#include <memory>
#include "Constants.h"
#include "Rndom.h"
#include "Container.h"
#include "EbyeInfo.h"

using std::string;
using std::cout;
using std::endl;

class Util_func{


 private:

	 std::string DATE;
	 std::shared_ptr<Rndom>& rndom;

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



  void get_EbyeInfo(const string fname, EbyeInfo& ebye, double rap_shift, bool VZEROAND_trigger, bool parton_level, int collision_type){
    std::ifstream in;
    in.open(fname.c_str(),std::ios::in);
    if(!in) {
      std::cerr << "Error unable to open file " << fname << endl;
      return ;
    } 
    string templine;


    double Multiplicity=0.0;
    double N_charge=0.0;
    double Multiplicity_V0M=0.0;
    double Multiplicity_V0A=0.0;
    double N_trk_offline_=0.0;
    bool Multiplicity_INEL_lg_0=false;
    bool V0M_FWD=false;
    bool V0M_BKW=false;
    bool OUTER_SPD=false;
    bool ATLAS_cut=false;
    double Weight=0.0;
    vector<Container::ParticleInfo> sampleSet;

    while(getline(in,templine)) {
      if(templine.find('#')!=std::string::npos) {
      }else if(templine.find('%')!=std::string::npos){
	      std::istringstream is(templine);
	      int iev, nv;
	      string com;
	      double weight_in, tau;
	      is >> com >> iev >> nv >> tau >> weight_in;
              Weight = weight_in; 
      }else{
	std::istringstream is(templine);
	int ID, data1, data2, col, acol;
	double mass, energy, px, py, pz, x,y,z,t,rap, ft;
	string TAG;
	is >> data1 >> data2 >> col >> acol >> ID >> mass >> energy >> px >> py >> pz >> rap >> x >> y >> z >> t >> ft >> TAG;
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
	double phi = atan2(py, px);

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
	
	//if(options.get_flag_SB_CMS()){
	if(pt>constants::assoc_ptmin && pt<constants::assoc_ptmax && fabs(eta)<constants::etaRangeCMSRidge){
		if((flag_only_corona && TAG == constants::corona_tag) || (flag_only_core && TAG == constants::core_tag)|| (!flag_only_core && !flag_only_corona )){
			if((flag_only_core_associates && TAG==constants::core_tag) || (flag_only_corona_associates && TAG==constants::corona_tag) || (!flag_only_corona_associates && !flag_only_core_associates)){
				Container::ParticleInfo part_in;
				part_in.pt=pt;
				part_in.id=ID;
				part_in.eta=eta;
				part_in.phi=phi;
				part_in.TAG=TAG;
				sampleSet.push_back(part_in);
			}
		}
	}

	
	if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon) && fabs(eta)<constants::w_eta_multiplicity) { 
		Multiplicity++;
	}
	if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon) && fabs(eta)<constants::etaRange_Nch && (pt<constants::ptmax_cumulantmulti_Nch && pt>constants::ptmin_cumulantmulti_Nch)) { 
		N_charge++;
	}
	if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon) && fabs(eta)<constants::w_eta_multiplicity_INEL_lg_0) { 
		Multiplicity_INEL_lg_0=true;
	}
	if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon) && this->in_V0M_detector(eta)) { 
		Multiplicity_V0M++;
	}
	if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon) && fabs(eta)<constants::etaRangeCMSRidge && pt > constants::MomentumMinCMSRidge) { 
		N_trk_offline_++;
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

    if((int)sampleSet.size()>0){
	    std::uniform_int_distribution<> rndomSamp(0, (int)sampleSet.size()-1);
	    for(int is=0; is<constants::Nsample; is++){
		    int iSamp=rndomSamp(rndom->generatorSamp);
		    //ebye.sample_part=sampleSet[iSamp];
		    ebye.sample_partSet.push_back(sampleSet[iSamp]);
	    }
	    ebye.valid_assoc=true;
    }
    vector<Container::ParticleInfo>().swap(sampleSet);

    in.close();
    ebye.multiplicity=Multiplicity;
    ebye.Nch=N_charge;
    ebye.weight=Weight;
    ebye.multiplicity_INEL_lg_0=Multiplicity_INEL_lg_0;
    ebye.multiplicity_V0M=(VZEROAND_trigger)? Multiplicity_V0A : Multiplicity_V0M;
    ebye.valid=true;
    ebye.N_trk_offline=N_trk_offline_;
    if(collision_type==101)
	    ebye.set_V0M_class(this->get_NtrkClass(N_trk_offline_));
    if(V0M_FWD && V0M_BKW && OUTER_SPD) {ebye.trig_3outof3=true;}
    if(V0M_FWD && V0M_BKW) {ebye.trig_VZEROAND=true;}
    if((V0M_FWD && V0M_BKW) || (V0M_BKW && OUTER_SPD) || (OUTER_SPD && V0M_FWD) ) {ebye.trig_2outof3=true;}
    if(ATLAS_cut) ebye.ATLAS_cut=true;
  }



int get_NtrkClass(const double val){

	for(int i=0; i<constants::n_NtrkClass_pp; i++){
		if(val>=0.0 && val<constants::val_NtrkClass_pp[0]) {
			return 0;
		}
		if(i>=1 && val>=constants::val_NtrkClass_pp[i-1] && val<constants::val_NtrkClass_pp[i]) {
			return i;
		}
	} 

	    cout << ":( ERROR Somthing wring in " << __LINE__ << " " << __FILE__ << endl;
	    exit(1);
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

 const std::string generateS2(int n){
   std::ostringstream name;
   name  << std::setw(2) << std::setfill('0') << n ;
   return name.str();
 };
 
  const std::string get_date(){
	  return this->DATE;
  };


  void get_nowdata(){
	  time_t now = time(NULL);
	  struct tm *pnow = localtime(&now);  
	  int Y=pnow->tm_year+1900;
	  int M=pnow->tm_mon + 1;
	  int D=pnow->tm_mday;
	  std::string date=generateS(Y)+generateS2(M)+generateS2(D);    
	  this->DATE=date;
  }


 public:

	 bool flag_only_core;
	 bool flag_only_corona;
	 bool flag_only_core_associates;
	 bool flag_only_corona_associates;

  Util_func(std::shared_ptr<Rndom>& rndom_in):DATE(""), rndom(rndom_in), flag_only_core(false), flag_only_corona(false), 
	flag_only_core_associates(false), flag_only_corona_associates(false){
		this->get_nowdata(); 
 };
  ~Util_func(){};
  
  std::string generateS(double n){
    std::ostringstream oss;
    oss <<  n ; 
    return oss.str();
  };
 


  bool checkMassOnShell(const double m, const double e, const double px, const double py, const double pz){

	  double P_squared=px*px+py*py+pz*pz;
	  double E_squared=e*e;

	  bool mos = fabs(m*m - E_squared + P_squared)<constants::SMALL;

	  if(mos){
		  //std::cout << ":) This is mass on shell " << std::fixed << std::setprecision(10) << "E_squared: " << E_squared << "   P_squared: " << P_squared<< "   M_squared: " << m*m << std::endl;
		  return true;
	  }else{
		  std::cout << ":O This is NOT mass on shell " << std::fixed << std::setprecision(10) << "E: " << e << "   P : " << sqrt(P_squared) << "   M: " << m << std::endl;
		  std::cout << "                             " << std::fixed << std::setprecision(10) << "E_squared: " << E_squared << "   P_squared: " << P_squared<< "   M_squared: " << m*m << "    mm - ee + pp " <<  fabs(m*m - E_squared + P_squared)  << std::endl;
		  return false;
	  }
  }



  void get_EbyeInfo_(const string fname, EbyeInfo& ebye, double rap_shift, bool VZEROAND_trigger, bool parton_level, int collision_type){
    this->get_EbyeInfo(fname, ebye, rap_shift, VZEROAND_trigger, parton_level, collision_type);
  }



 
  
  void make_output_directory(const std::string name_){
    
    if(stat(name_.c_str(),&st)!=0) mkdir(name_.c_str(), 0775);
    else {

    };
  };

  std::string get_output_directory(const std::string directory_name){
    
    std::ostringstream name; 
    name  << this->get_date()+"_"+directory_name;
    std::string name_=constants::data_directory+"/"+name.str();
    return name_.c_str();

  };


};
#endif
