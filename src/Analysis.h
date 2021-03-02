// -*- mode:c++ -*-
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <sstream>
#include <fstream>
#include <math.h>
#include "Constants.h"
#include "Util_func.h"
#include "Container.h"
#include "Message.h"
#include "LogSettings.h"
#include "Settings.h"
#include "CentralityCut.h"

using namespace std;


class Analysis{


	shared_ptr<Util_func> uf;
	shared_ptr<Message> ms;  

	Settings::Options options;
	LogSettings log;

	public:

	Analysis(const Settings::Options options_in, LogSettings log_in): options(options_in), log(log_in){
		this->ana();
		ms = make_shared<Message>();
		uf = make_shared<Util_func>();
	};
	~Analysis(){};





	bool read(const std::string& fname, shared_ptr<Container>& ct){

		ifstream in;
		in.open(fname.c_str(),ios::in);
		if(!in){ ms->open(fname); return false;}

		//Variables.
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

					////
					double eta;
					{
						const double LARGE = 10.0;
						if(P==0.0 && pz==0.0) eta=0.0;
						else if(P==-pz) eta=-LARGE;
						else if(P==pz) eta=LARGE;
						else eta=std::log((P + pz)/(P - pz)) / 2.0;
					}
					////



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

					}else if(constants::MODE.find("meanmt")!=string::npos){
						if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon||abs(ID)==constants::id_lambda ||abs(ID)==constants::id_cascade ) && std::fabs(eta)<0.5 ) { 

							Container::ParticleInfo part_in;
							part_in.mt=mt;
							part_in.m=m;
							part_1ev.push_back(part_in);

						}
						if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon) && std::fabs(eta)<0.3 && pt<10.0 ) Nch++;


					}else if(constants::MODE.find("meanpt")!=string::npos){
						if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon) && std::fabs(eta)<0.3 && pt>0.15 && pt<10.0 ) { 

							Container::ParticleInfo part_in;
							part_in.pt=pt;
							part_1ev.push_back(part_in);

						}
						if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon) && std::fabs(eta)<0.3 && pt<10.0 ) Nch++;


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
		//-------
		info_1ev.weight(weight);
		info_1ev.Nch(Nch);
		info_1ev.part=part_1ev;
		ct->EVENTINFO=info_1ev;

		return true;
	}



	//
	//	bool read_jetinfo(const std::string& fname, shared_ptr<Container>& ct){
	//
	//		ifstream in;
	//		in.open(fname.c_str(),ios::in);
	//		if(!in){ ms->open(fname); return false;}
	//
	//		vector <Container::ParticleInfo> Parents, Jets;
	//
	//		{
	//			std::string templine;
	//			while(getline(in,templine)) {
	//				if(templine.find('G')!=std::string::npos) {
	//				} else if(templine.find('P')!=std::string::npos){
	//					istringstream iss(templine);
	//					double e,px,py,pz;
	//					std::string TAG;
	//					iss >> TAG >> e >> px >> py >> pz;
	//					Container::ParticleInfo info;
	//					info.e=e;
	//					info.px=px;
	//					info.py=py;
	//					info.pz=pz;
	//					info.phi=atan2(py, px);
	//					Parents.push_back(info);
	//				}else{
	//					istringstream is(templine);
	//					double e,px,py,pz;
	//					int num;
	//					is >> num >> e >> px >> py >> pz;
	//					Container::ParticleInfo info;
	//					info.e=e;
	//					info.px=px;
	//					info.py=py;
	//					info.pz=pz;
	//					info.phi=atan2(py, px);
	//					Jets.push_back(info);
	//				}
	//			}
	//			in.close();
	//		}
	//
	//
	//		for(int i =0; i<(int)Jets.size(); i++){
	//			double phi_J = Jets[i].phi;
	//			for(int j =0; j<(int)Parents.size(); j++){
	//				if(std::fabs(phi_J-Parents[j].phi)<constants::delta_phi_SMALL){
	//
	//					double frac = (Jets[i].e - Parents[j].e)/Parents[j].e;
	//					ct->value.push_back(frac);
	//					ct->Nev_tot+=1;
	//
	//				}
	//
	//			}
	//		}
	//
	//		return true;
	//	}
	//

	void stat(shared_ptr<Container>& ct){

		//take average    
		//-------------------------------------
		for(int i=0; i<ct->max_nx+1; ++i){
			ct->Hist[i]/=ct->Hist_weight[i];
			ct->HistHist[i]/=ct->Hist_weight[i];

			// devide by cell width 
			//-------------------------------------
			//ct->Hist[i]/=constants::d_x;
		}


		//Get standard error
		//-------------------------------------
		for(int i=0; i<ct->max_nx+1; ++i){
			double var=ct->HistHist[i]-pow(ct->Hist[i],2.0);
			ct->HistErr[i]=sqrt(var/ct->Hist_weight[i]);
		}


	}





	void fill(shared_ptr<Container>& ct){


		Container::EventInfo& EVENT= ct->EVENTINFO;

			//Determine xbin
			//---------------
			double x_val=0.0;
			if(constants::MODE.find("meanpt")!=string::npos){
				x_val=EVENT.Nch();
			}
			if(x_val<constants::x_min || x_val>constants::x_max) return;
			int nx=(int)((x_val/constants::d_x)+(std::fabs(constants::x_min)/constants::d_x));


			//Count particle by particle.
			//----------------------------
			for(int j=0; j<(int)EVENT.part.size(); ++j){
				double y_val=0.0;
				if(constants::MODE.find("meanpt")!=string::npos){
					y_val = EVENT.part[j].pt;
				}else if(constants::MODE.find("meanmt")!=string::npos){
					x_val=EVENT.part[j].m;
					if(x_val<constants::x_min || x_val>constants::x_max) continue;
					nx=(int)((x_val/constants::d_x)+(std::fabs(constants::x_min)/constants::d_x));
					y_val = EVENT.part[j].mt - EVENT.part[j].m;
				}
				ct->Hist[nx]+=y_val*EVENT.weight();
				ct->HistHist[nx]+=pow(y_val,2)*EVENT.weight();
				ct->Hist_weight[nx]+=EVENT.weight();
				if(ct->max_nx<nx) ct->max_nx=nx;

				ct->SumWeight+=EVENT.weight();

			}


	}





	bool write(const std::string& fname, const shared_ptr<Container>& ct){
		ofstream ofs;
		ofs.open((fname+"/"+constants::default_out_fname).c_str());
		if(!ofs){ms->open(fname); return false;}

		ct->max_nx+=constants::margin;
		for(int i=0; i<ct->max_nx; ++i){
			double x_axis = ((constants::x_min+(constants::d_x*i))+(constants::x_min+(constants::d_x*(i+1))))/2.0;
			ofs << setw(16) << fixed << setprecision(8) << x_axis << "  "
				<< setw(16) << ct->Hist[i] << "  "
				<< setw(16) << ct->HistErr[i] << endl;
		}
		ofs << endl;
		return true;
	}


	int ana(){

		//Start Centrality Cut.
		//-----------------------
		int nCent=1;
		vector<EbyeInfo> eBye_CentCut;
		if(options.get_flag_CentralityCut()){
			CentralityCut CentCut(eBye_CentCut, options);
			nCent=(int)options.name_cent.size();
		}


		//Centrality Loop
		for(int iCent=0; iCent<nCent; iCent++){

			if(options.get_flag_CentralityCut())cout << ":)Start analyzing centrality " << options.name_cent[iCent] << "." << endl;



			auto ct = make_shared<Container>();

			//Event Loop
			for(int i=0; i<options.get_nfile(); ++i){
				if(!(i%1000)) ms->read(i);
				std::stringstream ss;
				ss << options.get_dir_name() << "/" << options.get_f_name() << setw(9) << setfill('0') << i << "/" << options.get_ext_name();
				std::string inputpath = ss.str();


				//Centrality cut.
				//---------------
				if(options.get_flag_CentralityCut()){
					if(eBye_CentCut[i].get_V0M_class()!=iCent) continue;
				}


				if(!options.get_flag_CentralityCut()){
					if(options.get_flag_INEL_lg_0()){
						eByeInSettings ebe;
						eByeInSettings::eByeMulti Multi(options, inputpath);
						if(!Multi.ebye.multiplicity_INEL_lg_0) continue;
						ebe.print_eByeInfo(i,Multi);
					}else if(options.get_flag_3outof3_trigger()){
						eByeInSettings ebe;
						eByeInSettings::eByeMulti Multi(options, inputpath);
						if(!Multi.ebye.trig_3outof3) continue;
						ebe.print_eByeInfo(i,Multi);
					}else if(options.get_flag_2outof3_trigger()){
						eByeInSettings ebe;
						eByeInSettings::eByeMulti Multi(options, inputpath);
						if(!Multi.ebye.trig_2outof3) continue;
						ebe.print_eByeInfo(i,Multi);
					}else if(options.get_flag_ATLAS_cut()){
						eByeInSettings ebe;
						eByeInSettings::eByeMulti Multi(options, inputpath);
						if(!Multi.ebye.ATLAS_cut) continue;
						ebe.print_eByeInfo(i,Multi);
					}else if(options.get_flag_VZEROAND_trigger()){
						eByeInSettings ebe;
						eByeInSettings::eByeMulti Multi(options, inputpath);
						if(!Multi.ebye.trig_VZEROAND) continue;
						ebe.print_eByeInfo(i,Multi);
					}else{
						eByeInSettings ebe;
						eByeInSettings::eByeMulti Multi(options, inputpath);
						ebe.print_eByeInfo(i,Multi);
					}
				}



				//Read events
				//---------------
				if(constants::MODE.find("JET_PRAC")!=std::string::npos){
					//read_jetinfo(inputpath, ct);
				}else{
					read(inputpath, ct);
				}
				this->fill(ct);


			}//Event loop

			this->stat(ct);


			//Making output directory namE
			//-----------------------------
			std::string generated_directory_name;
			if(options.get_flag_CentralityCut()){
				generated_directory_name=uf->get_output_directory(options.get_out_directory_name()+="_"+options.name_cent[iCent]);
			}else generated_directory_name=uf->get_output_directory(options.get_out_directory_name());
			uf->make_output_directory(generated_directory_name);


			//Output
			//------
			if(!write(generated_directory_name, ct)) return 1;
			log.archive_settings(generated_directory_name);

		}//Centrality loop


		return 0;
	}





};
