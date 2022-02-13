// -*- mode:c++ -*-
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <sstream>
#include <fstream>
#include <complex>
#include <math.h>
#include "Constants.h"
#include "Util_func.h"
#include "Container.h"
#include "Message.h"
#include "LogSettings.h"
#include "Settings.h"
#include "CentralityCut.h"
#include "ReadIn.h"

using namespace std;

ReadIn::ReadIn(shared_ptr<Message> ms_in, Settings::Options options_in):ms(ms_in), options(options_in), ncall_readTimeLapse(0){};
ReadIn::~ReadIn(){};

bool ReadIn::read(const std::string& fname, shared_ptr<Container>& ct){

			ifstream in;
			in.open(fname.c_str(),ios::in);
			if(!in){ ms->open(fname); return false;}

			//Variables.
			//part_1ev is a vector containing particle lists
			//info_1ev stores info of one single event.
			//-----------
			Container::EventInfo info_1ev;
			vector<Container::ParticleInfo> part_1ev;
			int Nch=0;
			double weight=1.0;


			{
				std::string templine;
				while(getline(in,templine)) {
					if(templine.find('#')!=std::string::npos) {
					} else if(templine.find('%')!=std::string::npos){
						istringstream iss(templine);
						std::string pct;
						int iev, nv;
						string com;
						double weight_in, tau;
						iss >> com >> iev >> nv >> tau >> weight_in;
						weight=weight_in;
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
						double mt_squared=pt*pt+m*m;
						double mt=(mt_squared)>0.0 ? sqrt(mt_squared):0.0;
						double phi = atan2(py, px);

						////
						double eta;
						{
							const double LARGE = 10.0;
							if(P==0.0 && pz==0.0) eta=0.0;
							else if(P==-pz) eta=-LARGE;
							else if(P==pz) eta=LARGE;
							else eta=std::log((P + pz)/(P - pz)) / 2.0;
						}
						///



						if(constants::MODE.find("dndeta_proton")!=string::npos){
							if(ID==constants::id_proton) { 

								Container::ParticleInfo part_in;
								part_in.eta=eta;
								part_1ev.push_back(part_in);

							}
						}else if(constants::MODE.find("dndpt")!=string::npos || constants::MODE.find("vertices")!=string::npos){
							if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon) && std::fabs(eta)<constants::delta_eta ) { 

								if(constants::MODE.find("dndpt")!=string::npos) {
									Container::ParticleInfo part_in;
									part_in.pt=pt;
									part_1ev.push_back(part_in);
								};
								double r_squared = x*x + y*y;
								double r=(r_squared>0.0)? sqrt(r_squared):0.0;
								if(constants::MODE.find("vertices")!=string::npos){
									Container::ParticleInfo part_in;
									part_in.r=r;
									part_1ev.push_back(part_in);
								}

							}

						}else if(constants::MODE.find("mtscaling")!=string::npos){
							if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon||abs(ID)==constants::id_phi||abs(ID)==constants::id_lambda ||abs(ID)==constants::id_cascade || abs(ID)==constants::id_omega ) && std::fabs(eta)<0.5 ) { 

								if((options.get_flag_only_corona() && TAG == constants::corona_tag) || (options.get_flag_only_core() && TAG == constants::core_tag) || (!options.get_flag_only_core() && !options.get_flag_only_corona() )){



									Container::ParticleInfo part_in;
									part_in.id=ID;
									part_in.mt=mt;
									part_in.e=e;
									part_in.px=px;
									part_in.py=py;
									part_in.pz=pz;
									part_in.m=m;
									part_1ev.push_back(part_in);
								}

							}
							if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon) && std::fabs(eta)<0.3 && pt<10.0 ) Nch++;


						}else if(constants::MODE.find("meanmt")!=string::npos){

							cout << "ERROR:( -DMEANMT is still under construction." << endl;
							exit(1);
							if(abs(ID)==constants::id_cascade && std::fabs(eta)<0.3 && pt>0.15 && pt<10.0 ) { 

								Container::ParticleInfo part_in;
								part_in.mt=mt;
								part_in.id=ID;
								part_1ev.push_back(part_in);

							}
							if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon) && std::fabs(eta)<0.3 && pt<10.0 ) Nch++;


						}else if(constants::MODE.find("meanpt")!=string::npos){
							if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon) && std::fabs(eta)<0.3 && pt>0.15 && pt<10.0 ) { 
								if((options.get_flag_only_corona() && TAG == constants::corona_tag) || (options.get_flag_only_core() && TAG == constants::core_tag)|| (!options.get_flag_only_core() && !options.get_flag_only_corona() )){

								Container::ParticleInfo part_in;
								part_in.pt=pt;
								part_1ev.push_back(part_in);

								}

							}
							if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon) && std::fabs(eta)<0.3 ) Nch++;


						}else if(constants::MODE.find("MeanptPID")!=string::npos){
							if(ID==constants::id_proton && std::fabs(eta)<0.5) { 
								if((options.get_flag_only_corona() && TAG == constants::corona_tag) || (options.get_flag_only_core() && TAG == constants::core_tag)|| (!options.get_flag_only_core() && !options.get_flag_only_corona() )){

								Container::ParticleInfo part_in;
								part_in.pt=pt;
								part_1ev.push_back(part_in);

								}

							}
							if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon) && std::fabs(eta)<0.3 ) Nch++;


						}else if(constants::MODE.find("Rt_spectra")!=string::npos){

							if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon ) && fabs(eta)<constants::etaRange_Rt && pt>0.15){

								if((options.get_flag_only_corona() && TAG == constants::corona_tag) || (options.get_flag_only_core() && TAG == constants::core_tag)|| (!options.get_flag_only_core() && !options.get_flag_only_corona() )){
									Container::ParticleInfo part_in;
									part_in.pt=pt;
									part_in.eta=eta;
									part_in.phi=phi;
									part_in.TAG=TAG;
									part_1ev.push_back(part_in);
								}
							}
							if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon) && fabs(eta)<constants::etaRangeMultiplicity) Nch++;


						}else if(constants::MODE.find("Rt_yield")!=string::npos){

							if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon ||abs(ID)==constants::id_phi||abs(ID)==constants::id_lambda ||abs(ID)==constants::id_cascade || abs(ID)==constants::id_omega) && fabs(eta)<constants::etaRange_Rt && pt>0.15){
								if((options.get_flag_only_corona() && TAG == constants::corona_tag) || (options.get_flag_only_core() && TAG == constants::core_tag)|| (!options.get_flag_only_core() && !options.get_flag_only_corona() )){
									Container::ParticleInfo part_in;
									part_in.id=ID;
									part_in.px=px;
									part_in.py=py;
									part_in.pt=pt;
									part_in.eta=eta;
									part_in.phi=phi;
									part_in.TAG=TAG;
									part_1ev.push_back(part_in);
								}
							}
							if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon) && fabs(eta)<constants::etaRangeMultiplicity) Nch++;


						}else if(constants::MODE.find("cumulant_multi")!=string::npos){

							if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon ) && fabs(eta)<constants::etaRange_cumulantmulti && constants::ptmin_cumulantmulti < pt && pt<constants::ptmax_cumulantmulti ){
								if((options.get_flag_only_corona() && TAG == constants::corona_tag) || (options.get_flag_only_core() && TAG == constants::core_tag)|| (!options.get_flag_only_core() && !options.get_flag_only_corona() )){
									Container::ParticleInfo part_in;
									part_in.pt=pt;
									part_in.eta=eta;
									part_in.phi=phi;
									part_1ev.push_back(part_in);
								}
							}
							if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon) && fabs(eta)<constants::etaRange_cumulantmulti_Nch && constants::ptmin_cumulantmulti_Nch < pt && constants::ptmax_cumulantmulti_Nch > pt) Nch++;


						}else if(constants::MODE.find("cumulant_pt")!=string::npos){

							if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon ) && constants::ptmin_cumulantpt<pt && pt < constants::ptmax_cumulantpt && fabs(eta)<constants::etaRange_cumulantpt){
								if((options.get_flag_only_corona() && TAG == constants::corona_tag) || (options.get_flag_only_core() && TAG == constants::core_tag)|| (!options.get_flag_only_core() && !options.get_flag_only_corona() )){
									Container::ParticleInfo part_in;
									part_in.pt=pt;
									part_in.eta=eta;
									part_in.phi=phi;
									part_1ev.push_back(part_in);
								}
							}


						}else if(constants::MODE.find("cumulant_eta")!=string::npos){

							if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon ) && constants::ptmin_cumulanteta < pt && pt<constants::ptmax_cumulanteta ){
								if((options.get_flag_only_corona() && TAG == constants::corona_tag) || (options.get_flag_only_core() && TAG == constants::core_tag)|| (!options.get_flag_only_core() && !options.get_flag_only_corona() )){
									Container::ParticleInfo part_in;
									part_in.pt=pt;
									part_in.eta=eta;
									part_in.phi=phi;
									part_1ev.push_back(part_in);
								}
							}
							if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon) && fabs(eta)<1.0 && pt<3.0 ) Nch++;


						}else if(constants::MODE.find("twopc2D")!=string::npos){

							if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon ) && pt<constants::twopc2Dptmax && pt>constants::twopc2Dptmin && fabs(eta)<constants::twopc2DetaRange){
								if((options.get_flag_only_corona() && TAG == constants::corona_tag) || (options.get_flag_only_core() && TAG == constants::core_tag)|| (!options.get_flag_only_core() && !options.get_flag_only_corona() )){
									Container::ParticleInfo part_in;
									part_in.id=ID;
									part_in.pt=pt;
									part_in.eta=eta;
									part_in.phi=phi;
									part_in.TAG=TAG;
									part_1ev.push_back(part_in);
								}
							}
							if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon) && fabs(eta)<constants::twopc2DetaRange_N && pt>constants::twopc2Dptmin_N && pt<constants::twopc2Dptmax_N ) Nch++;


						}else if(constants::MODE.find("twopc1D")!=string::npos || constants::MODE.find("twopcInteg")!=string::npos){

							if(((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon ) && fabs(eta)< constants::twopc1DetaRange && pt > constants::twopc1Dptmin && pt<constants::twopc1Dptmax)){

								if((options.get_flag_only_corona() && TAG == constants::corona_tag) || (options.get_flag_only_core() && TAG == constants::core_tag)|| (!options.get_flag_only_core() && !options.get_flag_only_corona() )){
									Container::ParticleInfo part_in;
									part_in.pt=pt;
									part_in.id=ID;
									part_in.eta=eta;
									part_in.phi=phi;
									part_in.TAG=TAG;
									part_1ev.push_back(part_in);
								}
							}
							if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon) && fabs(eta)<constants::twopc1DetaRange_N && pt>constants::twopc1Dptmin_N && pt<constants::twopc1Dptmax_N ) Nch++;


						}else{
							if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon) ) { 

								Container::ParticleInfo part_in;
								part_in.eta=eta;
								part_1ev.push_back(part_in);
							}
						}



						}
					}
					in.close();
				}

				//Store
				//part_1ev is a vector containing particle lists
				//info_1ev stores info of one single event.
				//-------
				info_1ev.weight(weight);
				info_1ev.Nch(Nch);
				info_1ev.part=part_1ev;
				ct->EVENTINFO=info_1ev;
				vector<Container::ParticleInfo>().swap(part_1ev);

				return true;
			}



			
				bool ReadIn::read_jetinfo(const std::string& fname, shared_ptr<Container>& ct){

					ifstream in;
					in.open(fname.c_str(),ios::in);
					if(!in){ ms->open(fname); return false;}
			
					Container::EventInfo info_1ev;
					vector <Container::ParticleInfo> Jets, Gamma;
					{
						std::string templine;
						while(getline(in,templine)) {
							if(templine.find('G')!=std::string::npos) {
								istringstream is(templine);
								double e,px,py,pz;
								int num;
								is >> num >> e >> px >> py >> pz;
								Container::ParticleInfo info;
								info.e=e;
								info.px=px;
								info.py=py;
								info.pz=pz;
								info.pt=sqrt(px*px + py*py);
								info.phi=atan2(py, px);
								Gamma.push_back(info);
							}else if(templine.find('P')!=std::string::npos){
							}else if(templine.find('#')!=std::string::npos){
							}else{
								istringstream is(templine);
								double e,px,py,pz;
								int num;
								is >> num >> e >> px >> py >> pz;
								Container::ParticleInfo info;
								info.e=e;
								info.px=px;
								info.py=py;
								info.pz=pz;
								info.pt=sqrt(px*px + py*py);
								info.phi=atan2(py, px);
								Jets.push_back(info);
							}
						}
						in.close();
					}
			
			
					if((int)Jets.size()<2){ 
						vector<Container::ParticleInfo>().swap(Jets);
						vector<Container::ParticleInfo>().swap(Gamma);
						return false;
					}
					for(int i =0; i<1; i++){
						if((Jets[i].pt+Jets[i+1].pt)<=0.0) {
							vector<Container::ParticleInfo>().swap(Jets);
							vector<Container::ParticleInfo>().swap(Gamma);
							return false;
						}
						info_1ev.Aj(fabs((Jets[i].pt-Jets[i+1].pt)/(Jets[i].pt+Jets[i+1].pt)));
					}
			
					ct->EVENTINFO=info_1ev;

					vector<Container::ParticleInfo>().swap(Jets);
					vector<Container::ParticleInfo>().swap(Gamma);
					return true;
				}



bool ReadIn::readTimeLapse(const std::string& fname, shared_ptr<Container>& ct, const double weight){


	//- at first event, get line num of the point to see.
	int nline=-1;
	if (ncall_readTimeLapse==0)
		this->get_nline_to_see(nline, fname);
	if(nline<0){cout << ":(ERROR " << __FILE__ << __LINE__ << endl; exit(1);}


	//Variables.
	//step_1ev is a vector containing time step info lists
	//info_1ev stores info of one single event.
	//-----------
	Container::EventInfo info_1ev;
	vector<Container::StepInfo> step_1ev;

	//TimeStep
	//---------
	double tau=constants::TL_tau_00;
	for(int step=0; step < constants::LARGEint ; step++){
		double dtau =(tau<= constants::TL_tau_switch)? constants::TL_dtau1:constants::TL_dtau2;
		tau+=dtau;


		//Read
		//----
		ostringstream os;
		os << fname << options.get_ext_nameTL() << fixed << setprecision(2) << tau*constants::hbarc << ".txt";
		ifstream in;
		in.open(os.str().c_str(),ios::in);
		if(!in){return false;}

		std::string templine;
		int nline_current=0;
		while(getline(in,templine)) {
			if(templine.find('#')!=std::string::npos) {
			} else if(templine.find('%')!=std::string::npos){
			}else{
				if(nline_current==nline){
					istringstream is(templine);
					Container::StepInfo step_oneline;
					if(options.get_modeTL()==0) this->get_oneline_xy(is, step_oneline);
					else if(options.get_modeTL()==1)this->get_oneline_xeta(is, step_oneline);
					else{cout << ":( ERROR Something wrong. " << __FILE__ << __LINE__ << endl; exit(1);}
					cout << "get " << step_oneline.tau << "  " << step_oneline.x << "  " << step_oneline.eta << "  " << step_oneline.temp << endl;
					step_1ev.push_back(step_oneline);
				}
			}
			nline_current++;
		}//while

		in.close();

		//Step End
		//---------
		tau+=dtau;
		if(tau*constants::hbarc>100.0) {
			cout << ":oWARNING It seems that something is wrong in hydro profile. The time step loop continues to 100fm." << endl;
			break;
		}

	}//step


	//Store
	//-------
	info_1ev.weight(weight);
	info_1ev.step=step_1ev;
	ct->EVENTINFO=info_1ev;
	vector<Container::StepInfo>().swap(step_1ev);

	ncall_readTimeLapse++;
	return true;

}

bool ReadIn::get_nline_to_see(int &nline, const std::string fname){


	double tau=constants::TL_tau_00;
	tau+=constants::TL_dtau1;
	ostringstream os;
	os << fname << options.get_ext_nameTL() << fixed << setprecision(2) << tau*constants::hbarc << ".txt";
	ifstream in;
	in.open(os.str().c_str(),ios::in);
	if(!in){cout << os.str() << endl;return false;}
	std::string templine;
	double delta=constants::LARGE;
	double xTL = options.get_at_xTL();
	double yTL = options.get_at_yTL();
	double etaTL = options.get_at_etaTL();
	double x_current, y_current, eta_current;
	int nline_=0;
	while(getline(in,templine)) {
		if(templine.find('#')!=std::string::npos) {
		} else if(templine.find('%')!=std::string::npos){
		}else{
			istringstream is(templine);
			Container::StepInfo step_oneline;
			if(options.get_modeTL()==0)this->get_oneline_xy(is, step_oneline);
			else if(options.get_modeTL()==1)this->get_oneline_xeta(is, step_oneline);
			else{cout << ":( ERROR Something wrong. " << __FILE__ << __LINE__ << endl; exit(1);}
			if(options.get_modeTL()==0){
				if(delta>(pow(fabs(xTL-step_oneline.x),2)+pow(fabs(yTL-step_oneline.y),2))){
					x_current=step_oneline.x;
					y_current=step_oneline.y;
					nline=nline_;
					delta=(pow(fabs(xTL-step_oneline.x),2)+pow(fabs(yTL-step_oneline.y),2));
				}
				//if(delta_y>fabs(yTL-step_oneline.y)) y_current=step_online.y;
			}else if(options.get_modeTL()==1){
				if(delta>(pow(fabs(xTL-step_oneline.x),2)+pow(fabs(etaTL-step_oneline.eta),2))){
					x_current=step_oneline.x;
					eta_current=step_oneline.eta;
					nline=nline_;
					delta=(pow(fabs(xTL-step_oneline.x),2)+pow(fabs(etaTL-step_oneline.eta),2));
				}
			}

		}
		nline_++;
	}//while

	cout << "Printing nearest points: x " << x_current << "  y " << y_current << "  eta " << eta_current << endl;
	cout << "                         ==>  " << nline << endl;
	cout << "                 inputs: x " << options.get_at_xTL() << "  y " << options.get_at_yTL() << "  eta " << options.get_at_etaTL()<< endl;

	in.close();
return true;

}


void ReadIn::get_oneline_xeta(istringstream& is, Container::StepInfo& onestep){
			is >> onestep.tau //fm 
				>> onestep.x //fm
				>> onestep.eta //[1]
				>> onestep.e //[GeV^4]
				>> onestep.temp //[GeV]
				>> onestep.s //[GeV^3]
				>> onestep.p //[GeV]
				>> onestep.n //[GeV^3]
				>> onestep.n5 
				>> onestep.mu 
				>> onestep.mu5
				>> onestep.vx  //[1]
				>> onestep.veta  //[1]
				>> onestep.Ueta
				>> onestep.Vtilde
				>> onestep.UABS
				>> onestep.Ux
				>> onestep.Uy
				>> onestep.Ueta
				>> onestep.eEtildedx
				>> onestep.eEtildedy
				>> onestep.eBtildedx
				>> onestep.eBtildedy
				>> onestep.eEtilde_dot_eBtilde
				>> onestep.U4;
}

void ReadIn::get_oneline_xy(istringstream& is, Container::StepInfo& onestep){
			is >> onestep.tau //1
				>> onestep.x //2
				>> onestep.y //3
				>> onestep.e //4
				>> onestep.temp //5
				>> onestep.s //6
				>> onestep.p //7
				>> onestep.n //8
				>> onestep.n5 //9
				>> onestep.mu //10
				>> onestep.mu5 //11
				>> onestep.Vtilde //12
				>> onestep.vx //13
				>> onestep.vy //14
				>> onestep.veta //15	   
				>> onestep.U4 //16
				>> onestep.U_R4 //17
				>> onestep.U_L4 //18	
				>> onestep.eEtildedx  //19
				>> onestep.eEtildedy  //20
				>> onestep.eBtildedx  //21
				>> onestep.eBtildedy  //22
				>> onestep.eEtilde_dot_eBtilde  //23
				>> onestep.Ux //24
				>> onestep.Uy //25	
				>> onestep.Ueta; //26	

}
