#ifndef READ
#define READ

class ReadIn{

private:
		shared_ptr<Message> ms;  
		Options options;
		int ncall_readTimeLapse;
		int nline;

		class print_coodTL{
			public:
				double x_print, y_print, eta_print, tau_print;
		};

		void stockInfo(double& Multiplicity, double& N_charge, double& Multiplicity_V0M, double &Multiplicity_V0A, double& N_trk_offline_, bool& Multiplicity_INEL_lg_0, bool &V0M_FWD, bool &V0M_BKW, bool& OUTER_SPD, bool& ATLAS_cut, Container::ParticleInfo &p){


			double eta=p.rap;

			//Assuming EKRT input which is partons. 
			//So far I do not care about hadron species.
			//=============================================
			//if(is_pikp(p.ID)){
			if(fabs(eta)<constants::w_eta_multiplicity) { 
				Multiplicity++;
			}
			if(fabs(eta)<constants::etaRange_Nch && (p.pt<constants::ptmax_cumulantmulti_Nch && p.pt>constants::ptmin_cumulantmulti_Nch)) { 
				N_charge++;
			}
			if(abs(eta)<constants::w_eta_multiplicity_INEL_lg_0) { 
				Multiplicity_INEL_lg_0=true;
			}
			if(this->in_V0M_detector(eta)) { 
				Multiplicity_V0M++;
			}
			if(fabs(eta)<constants::etaRangeCMSRidge && p.pt > constants::MomentumMinCMSRidge) { 
				N_trk_offline_++;
			}
			if(this->in_V0M_fwd(eta)) { V0M_FWD=true; Multiplicity_V0A++;} 
			if(this->in_V0M_bkw(eta)) { V0M_BKW=true;}
			if(this->in_outerSPD(eta)){ OUTER_SPD=true;}
			if(p.pt>0.5 && fabs(eta)<constants::w_eta_ATLAS_cut) ATLAS_cut=true;
			//}

			return;
		}

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

		bool is_pikp(int ID){
			if (abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon) 
				return true;
			else 
				return false;
		}


public:

ReadIn(shared_ptr<Message> ms_in, Options options_in);
~ReadIn();

bool read(const std::string& fname, Container::EventInfo& oneEventInfo, EbyeInfo& ebye_);
bool readEKRT(const std::string& fname, Container::EventInfo& oneEventInfo);
bool readEKRTbinary(std::vector <Container::EventInfo>& nEventInfo, std::vector<EbyeInfo>& eBye);
bool readXY(const std::string& fname, Container::EventInfo& oneEventInfo);
bool read_jetinfo(const std::string& fname, Container::EventInfo& oneEventInfo);
bool readTimeLapse(const std::string& fname, Container::EventInfo& oneEventInfo, const double weight);
bool get_nline_to_see(int &nline, const std::string fname);
void get_oneline_xeta(istringstream& is, Container::StepInfo& onestep);
void get_oneline_xy(istringstream& is, Container::StepInfo& onestep);
};
#endif
