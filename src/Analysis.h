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

using namespace std;



class Analysis{

	private:


		shared_ptr<Util_func> uf;
		shared_ptr<Message> ms;  

		Settings::Options options;
		LogSettings log;


		//Maximum value for histgram
		//-------------------------
		double x_max;
		double d_x;

	public:

		Analysis(const Settings::Options options_in, LogSettings log_in): options(options_in), log(log_in), x_max(constants::x_max), d_x(constants::d_x){
			ms = make_shared<Message>();
			uf = make_shared<Util_func>();
			if(options.get_flag_HI()){
				cout << ":) HI option is called.\n  Maximum value of xaxis is adjusted to HI data." << endl;
				x_max=constants::x_max_HI;
				d_x=constants::d_x_HI;
			}
			this->ana();
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
							if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon||abs(ID)==constants::id_phi||abs(ID)==constants::id_lambda ||abs(ID)==constants::id_cascade ) && std::fabs(eta)<0.5 ) { 

								if((options.get_flag_only_corona() && TAG == constants::corona_tag) || (options.get_flag_only_core() && TAG == constants::core_tag) || (!options.get_flag_only_core() && !options.get_flag_only_corona() )){


									Container::ParticleInfo part_in;
									part_in.id=ID;
									part_in.mt=mt;
									part_in.m=m;
									part_1ev.push_back(part_in);
								}

							}
							if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon) && std::fabs(eta)<0.3 && pt<10.0 ) Nch++;


						}else if(constants::MODE.find("MtNch")!=string::npos){
							//if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon) && std::fabs(eta)<0.3 && pt>0.15 && pt<10.0 ) { 
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
							if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon) && std::fabs(eta)<0.3 && pt<10.0 ) Nch++;


						}else if(constants::MODE.find("cumulant_multi")!=string::npos){

							if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon ) && fabs(eta)<1.0 && 0.2 < pt && pt<3.0 ){
								if((options.get_flag_only_corona() && TAG == constants::corona_tag) || (options.get_flag_only_core() && TAG == constants::core_tag)|| (!options.get_flag_only_core() && !options.get_flag_only_corona() )){
									Container::ParticleInfo part_in;
									part_in.pt=pt;
									part_in.eta=eta;
									part_in.phi=phi;
									part_1ev.push_back(part_in);
								}
							}
							if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon) && fabs(eta)<1.0 && 0.2 < pt && pt<3.0 ) Nch++;


						}else if(constants::MODE.find("cumulant_pt")!=string::npos){

							if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon ) && fabs(eta)<1.0 ){
								if((options.get_flag_only_corona() && TAG == constants::corona_tag) || (options.get_flag_only_core() && TAG == constants::core_tag)|| (!options.get_flag_only_core() && !options.get_flag_only_corona() )){
									Container::ParticleInfo part_in;
									part_in.pt=pt;
									part_in.eta=eta;
									part_in.phi=phi;
									part_1ev.push_back(part_in);
								}
							}
							if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon) && fabs(eta)<1.0 && 0.2 < pt && pt<3.0 ) Nch++;


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
					ct->Hist_x[i]/=ct->Hist_weight[i];
					ct->HistHist[i]/=ct->Hist_weight[i];
					ct->Hist_img_Qvec[i]/=ct->Hist_weight[i];

					// devide by cell width 
					//-------------------------------------
					//ct->Hist[i]/=this->d_x;

					if((constants::MODE.find("cumulant_multi")!=string::npos || constants::MODE.find("cumulant_pt")!=string::npos)){
						if(options.get_flag_vn()){

							//Obtain vn{2} = sqrt(cn{2})
							//--------------------------- 
							ct->Hist[i]=sqrt(ct->Hist[i]);
							ct->HistHist[i]=sqrt(ct->HistHist[i]);

						}

						if(fabs(ct->Hist_img_Qvec[i])>constants::WARNING_IMAGINARY_MAX){
							ms->WARNING_LARGE_IMAGINARYPART(ct->Hist_img_Qvec[i]);
						}

					}
				}


				//Get standard error
				//-------------------------------------
				for(int i=0; i<ct->max_nx+1; ++i){
					double var=ct->HistHist[i]-pow(ct->Hist[i],2.0);
					double error=constants::dummy;
					if((constants::MODE.find("cumulant_multi")!=string::npos || constants::MODE.find("cumulant_pt")!=string::npos) && options.get_flag_vn()){
						double error_=sqrt(var/ct->HistHit[i]);
						error=(1.0/2.0)*(error_/sqrt(ct->Hist[i]));
					}else{
						error=sqrt(var/ct->HistHit[i]);
					}

					ct->HistErr[i]=error;
				}


			}





			void fill(shared_ptr<Container>& ct){


				Container::EventInfo& EVENT= ct->EVENTINFO;

				//Determine xbin
				//---------------
				double x_val=constants::dummy;
				if(constants::MODE.find("meanpt")!=string::npos || constants::MODE.find("MtNch")!=string::npos){
					x_val=EVENT.Nch();
					if(x_val<constants::x_min || x_val>this->x_max) return;
				}
				int nx=(!options.get_flag_HI())? (int)((x_val/this->d_x)+(std::fabs(constants::x_min)/this->d_x)):this->get_cell_index_logplot(x_val);


				//Count particle by particle.
				//----------------------------
				for(int j=0; j<(int)EVENT.part.size(); ++j){
					double y_val=constants::dummy;
					if(constants::MODE.find("meanpt")!=string::npos){
						y_val = EVENT.part[j].pt;
					}else if(constants::MODE.find("MtNch")!=string::npos){
						y_val = EVENT.part[j].mt - EVENT.part[j].m;
					}else if(constants::MODE.find("meanmt")!=string::npos){
						x_val=EVENT.part[j].m;
						if(x_val<constants::x_min || x_val>this->x_max) continue;
						nx=(int)((x_val/this->d_x)+(std::fabs(constants::x_min)/this->d_x));
						this->fix_ax(EVENT.part[j].id, nx, x_val);
						y_val = EVENT.part[j].mt - EVENT.part[j].m;
					}
					ct->Hist[nx]+=y_val*EVENT.weight();
					ct->Hist_x[nx]+=x_val*EVENT.weight();
					ct->HistHit[nx]++;
					ct->HistHist[nx]+=pow(y_val,2)*EVENT.weight();
					ct->Hist_weight[nx]+=EVENT.weight();
					if(ct->max_nx<nx) ct->max_nx=nx;

					ct->SumWeight+=EVENT.weight();

				}


			}


			void fix_ax(const int id, int &nx, double m){

				if(nx!=2 && nx!=9 && nx!=18){

					if(id==constants::id_ch_pion) nx=2;
					else if(id==constants::id_ch_kaon) nx=9;
					else if(id==constants::id_proton) nx=18;
					else if(id==constants::id_phi) nx=20;
					else if(id==constants::id_lambda) nx=22;
					else if(id==constants::id_cascade) nx=26;


				}else return;

			}



			int get_cell_index_logplot(const double x_val_){

				double x_val=x_val_;
				int ncell=0;


				if(x_val<constants::switchBin_x){
					int ncell=(int) floor((x_val-constants::x_min)/constants::binSize_small);
					return ncell;
				}else{

					int ncell_log=1;
					int ncell_start=floor(constants::switchBin_x/constants::binSize_small)-1;
					while(ncell_log<constants::x_cell_capa){

						if(x_val<pow(constants::binSize_log,ncell_log)+constants::switchBin_x){

							ncell=ncell_start+ncell_log;
							break;

						}
						ncell_log++;

					}
					return ncell;
				}
			}



			void fill_vnmulti(shared_ptr<Container>& ct){


				Container::EventInfo& EVENT= ct->EVENTINFO;

				if((int)EVENT.part.size()<2) return;

				//Determine xbin
				//---------------
				double x_val=EVENT.Nch();
				if(x_val<constants::x_min || x_val>this->x_max) return;
				int nx=(!options.get_flag_HI())? (int)((x_val/this->d_x)+(fabs(constants::x_min)/this->d_x)):this->get_cell_index_logplot(x_val);

				//Count particle by particle.
				//----------------------------
				std::complex<double> Qvec_tot=constants::initialval_comp;
				std::complex<double> n_coeff (2.0, 0.0);
				for(int j=0; j<(int)EVENT.part.size(); ++j){
					std::complex<double> phi_ (EVENT.part[j].phi,0.0);
					std::complex<double> Qvec=exp(constants::i_img*n_coeff*phi_);
					Qvec_tot += Qvec;
				}
				double squared_Qvec = real(Qvec_tot * conj(Qvec_tot));
				double squared_Qvec_img = imag(Qvec_tot * conj(Qvec_tot));
				double totN = (double)EVENT.part.size();
				double corr = (squared_Qvec-totN)/(totN*(totN-1.0));

				ct->Hist[nx]+=corr*EVENT.weight();
				ct->Hist_x[nx]+=x_val*EVENT.weight();
				ct->HistHit[nx]++;
				ct->HistHist[nx]+=corr*corr*EVENT.weight();
				ct->Hist_weight[nx]+=EVENT.weight();
				ct->Hist_img_Qvec[nx]+=squared_Qvec_img*EVENT.weight();
				if(ct->max_nx<nx) ct->max_nx=nx;

				ct->SumWeight+=EVENT.weight();

			}





			void fill_vnpt(shared_ptr<Container>& ct){


				Container::EventInfo& EVENT= ct->EVENTINFO;



				//Count particle by particle.
				//----------------------------
				std::complex<double> Qvec_tot[constants::x_cell_capa]={};
				double hit[constants::x_cell_capa]={};
				for(int i=0; i<constants::x_cell_capa; i++){
					Qvec_tot[i]=constants::initialval_comp;
					hit[i]=0.0;
				}
				std::complex<double> n_coeff (2.0, 0.0);
				for(int j=0; j<(int)EVENT.part.size(); ++j){

					//Determine xbin
					//---------------
					double x_val=EVENT.part[j].pt;
					if(x_val<constants::x_min || x_val>this->x_max) continue;
					int nx=(int)((x_val/this->d_x)+(fabs(constants::x_min)/this->d_x));

					std::complex<double> phi_ (EVENT.part[j].phi,0.0);
					std::complex<double> Qvec=exp(constants::i_img*n_coeff*phi_);
					Qvec_tot[nx] += Qvec;
					//      cout << "Qvec " << Qvec << endl;
					hit[nx]++;

					if(ct->max_nx<nx) ct->max_nx=nx;

				}

				for(int nx=0; nx<ct->max_nx+constants::margin; nx++){
					double squared_Qvec = real(Qvec_tot[nx] * conj(Qvec_tot[nx]));
					double squared_Qvec_img = imag(Qvec_tot[nx] * conj(Qvec_tot[nx]));
					double totN = hit[nx];
					if(totN<2) continue;
					double corr = (squared_Qvec-totN)/(totN*(totN-1.0));

					ct->Hist[nx]+=corr*EVENT.weight();
					ct->HistHit[nx]++;
					ct->HistHist[nx]+=corr*corr*EVENT.weight();
					ct->Hist_weight[nx]+=EVENT.weight();
					ct->Hist_img_Qvec[nx]+=squared_Qvec_img*EVENT.weight();

					ct->SumWeight+=EVENT.weight();
				}


			}




			//-------2sub


			void fill_vnpt_sub(shared_ptr<Container>& ct){


				Container::EventInfo& EVENT= ct->EVENTINFO;

				if((int)EVENT.part.size()<2) return;


				//Count particle by particle.
				//----------------------------
				std::complex<double> Qvec_tot_A[constants::x_cell_capa]={};
				std::complex<double> Qvec_tot_B[constants::x_cell_capa]={};
				double hit_A[constants::x_cell_capa]={};
				double hit_B[constants::x_cell_capa]={};
				for(int i=0; i<constants::x_cell_capa; i++){
					Qvec_tot_A[i]=constants::initialval_comp;
					Qvec_tot_B[i]=constants::initialval_comp;
					hit_A[i]=0.0;
					hit_B[i]=0.0;
				}
				std::complex<double> n_coeff (2.0, 0.0);
				for(int j=0; j<(int)EVENT.part.size(); ++j){

					//Determine xbin
					//---------------
					double x_val=EVENT.part[j].pt;
					if(x_val<constants::x_min || x_val>this->x_max) continue;
					int nx=(int)((x_val/this->d_x)+(fabs(constants::x_min)/this->d_x));

					std::complex<double> phi_ (EVENT.part[j].phi,0.0);
					std::complex<double> Qvec=exp(constants::i_img*n_coeff*phi_);
					if(EVENT.part[j].eta<-0.0){
						Qvec_tot_A[nx] += Qvec;
						hit_A[nx]++;
					}else if(EVENT.part[j].eta>0.0){
						Qvec_tot_B[nx] += Qvec;
						hit_B[nx]++;
					}

					if(ct->max_nx<nx) ct->max_nx=nx;

				}

				for(int nx=0; nx<ct->max_nx+constants::margin; nx++){
					double squared_Qvec = real(Qvec_tot_A[nx] * conj(Qvec_tot_B[nx]));
					double squared_Qvec_img = imag(Qvec_tot_A[nx] * conj(Qvec_tot_B[nx]));
					if(hit_A[nx]==0.0 || hit_B[nx]==0.0) continue;
					double corr = (squared_Qvec)/(hit_A[nx]*hit_B[nx]);

					ct->Hist[nx]+=corr*EVENT.weight();
					ct->HistHit[nx]++;
					ct->HistHist[nx]+=corr*corr*EVENT.weight();
					ct->Hist_weight[nx]+=EVENT.weight();
					ct->Hist_img_Qvec[nx]+=squared_Qvec_img*EVENT.weight();

					ct->SumWeight+=EVENT.weight();
				}


			}





			void fill_vnmulti_sub(shared_ptr<Container>& ct){


				Container::EventInfo& EVENT= ct->EVENTINFO;


				//Determine xbin
				//---------------
				double x_val=EVENT.Nch();
				if(x_val<constants::x_min || x_val>this->x_max) return;
				int nx=(!options.get_flag_HI())? (int)((x_val/this->d_x)+(fabs(constants::x_min)/this->d_x)):this->get_cell_index_logplot(x_val);

				//Count particle by particle.
				//----------------------------
				std::complex<double> Qvec_tot_A=constants::initialval_comp;
				std::complex<double> Qvec_tot_B=constants::initialval_comp;
				std::complex<double> n_coeff (2.0, 0.0);
				double hit_A=0.0, hit_B=0.0;
				for(int j=0; j<(int)EVENT.part.size(); ++j){
					std::complex<double> phi_ (EVENT.part[j].phi,0.0);
					std::complex<double> Qvec=exp(constants::i_img*n_coeff*phi_);
					if(EVENT.part[j].eta<-0.0){
						Qvec_tot_A += Qvec;
						hit_A++;
					}else if(EVENT.part[j].eta>0.0){
						Qvec_tot_B += Qvec;
						hit_B++;
					}
				}
				double squared_Qvec = real(Qvec_tot_A * conj(Qvec_tot_B));
				double squared_Qvec_img = imag(Qvec_tot_A * conj(Qvec_tot_B));

				if(hit_A==0.0 || hit_B==0.0) return;

				double corr = (squared_Qvec)/(hit_A*hit_B);

				ct->Hist[nx]+=corr*EVENT.weight();
				ct->Hist_x[nx]+=x_val*EVENT.weight();
				ct->HistHit[nx]++;
				ct->HistHist[nx]+=corr*corr*EVENT.weight();
				ct->Hist_weight[nx]+=EVENT.weight();
				ct->Hist_img_Qvec[nx]+=squared_Qvec_img*EVENT.weight();
				if(ct->max_nx<nx) ct->max_nx=nx;

				ct->SumWeight+=EVENT.weight();

			}


			//----------2sub











			bool write(const std::string& fname, const shared_ptr<Container>& ct){
				ofstream ofs;
				ofs.open((fname+"/"+constants::default_out_fname).c_str());
				cout <<"write" << endl;
				if(!ofs){ms->open(fname+"/"+constants::default_out_fname); return false;}

				ct->max_nx+=constants::margin;
				for(int i=0; i<ct->max_nx; ++i){

					if(ct->HistHit[i]==0) continue;

					double x_axis =(constants::MODE.find("cumulant_pt")!=string::npos)? ((constants::x_min+(this->d_x*i))+(constants::x_min+(this->d_x*(i+1))))/2.0: ct->Hist_x[i];
					ofs << setw(16) << fixed << setprecision(8) << x_axis << "  "
						<< setw(16) << ct->Hist[i] << "  "
						<< setw(16) << ct->HistErr[i] << "  "
						<< setw(16) << ct->HistHit[i] << endl;
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
							if(!read(inputpath, ct)) continue;
						}
						if(constants::MODE.find("cumulant_multi")!=std::string::npos){
							if(options.get_flag_2subevent()) this->fill_vnmulti_sub(ct); 
							else this->fill_vnmulti(ct);
						}else if(constants::MODE.find("cumulant_pt")!=std::string::npos){
							if(options.get_flag_2subevent()) this->fill_vnpt_sub(ct); 
							else this->fill_vnpt(ct);
						}else{
							this->fill(ct);
						}

					}//Event loop

					this->stat(ct);


					//Making output directory name
					//-----------------------------
					std::string generated_directory_name;
					if(options.get_flag_CentralityCut()){
						if(iCent==0){
							generated_directory_name=uf->get_output_directory(options.get_out_directory_name());
							uf->make_output_directory(generated_directory_name);
						}
						generated_directory_name=uf->get_output_directory(options.get_out_directory_name()+="/Cent"+options.name_cent[iCent]);
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
