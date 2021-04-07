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

							if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon ) && fabs(eta)<constants::etaRange_cumulantmulti && 0.2 < pt && pt<3.0 ){
								if((options.get_flag_only_corona() && TAG == constants::corona_tag) || (options.get_flag_only_core() && TAG == constants::core_tag)|| (!options.get_flag_only_core() && !options.get_flag_only_corona() )){
									Container::ParticleInfo part_in;
									part_in.pt=pt;
									part_in.eta=eta;
									part_in.phi=phi;
									part_1ev.push_back(part_in);
								}
							}
							if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon) && fabs(eta)<constants::etaRange_cumulantmulti && 0.2 < pt && pt<3.0 ) Nch++;


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
					ct->FinalHist[i]=ct->Hist[i]/ct->Hist_weight[i];
					ct->Hist_x[i]/=ct->Hist_weight[i];
					double meanxx  = ct->HistHist[i]/ct->Hist_weight[i];
					double meanx = ct->FinalHist[i];

					// devide by cell width 
					//-------------------------------------
					//ct->Hist[i]/=this->d_x;


					//Get standard error
					//-------------------------------------
					double var=meanxx-pow(meanx,2.0);
					double error=sqrt(var/ct->HistHit[i]);
					ct->HistErr[i]=error;
				}


			}





			void stat_flow(shared_ptr<Container>& ct){

				//take average    
				//-------------------------------------
				for(int i=0; i<ct->max_nx+1; ++i){
					ct->Hist[i]/=ct->Hist_weight[i];
					ct->Hist_sub[i]/=ct->Hist_weight[i];
					ct->Hist_subsub[i]/=ct->Hist_weight[i];
					ct->Hist_x[i]/=ct->Hist_weight[i];
					ct->HistHist[i]/=ct->Hist_weight[i];
					ct->HistHist_sub[i]/=ct->Hist_weight[i];
					ct->HistHist_subsub[i]/=ct->Hist_weight[i];
					ct->Hist_img_Qvec[i]/=ct->Hist_weight[i];

					// devide by cell width 
					//-------------------------------------
					//ct->Hist[i]/=this->d_x;

					if(fabs(ct->Hist_img_Qvec[i])>constants::WARNING_IMAGINARY_MAX){
						ms->WARNING_LARGE_IMAGINARYPART(ct->Hist_img_Qvec[i]);
					}
				}

				//CUMULANT
				//-------------

					if(!options.get_flag__4particle()){

						//c2{2}
						//---------

						//Get standard error
						//c2{2} = <<2>>, so no need to worry about error propagation.
						//-------------------------------------
						for(int i=0; i<ct->max_nx+1; ++i){
							double var=ct->HistHist[i]-pow(ct->Hist[i],2.0);
							ct->HistErr[i]=sqrt(var/ct->HistHit[i]);
							ct->FinalHist[i]=ct->Hist[i];
						}


					}else{


						if(!options.get_flag_3subevent()){

							//c2{4} (standard., 2subevent)
							//---------
							for(int i=0; i<ct->max_nx+1; ++i){
								//Get c2{4} = <<4>>-2*<<2>>^2
								//------------------------
								double c24=ct->Hist[i] - 2.0 * pow(ct->Hist_sub[i],2);

								//Prepare standard error of <<2>> and <<4>>.
								//----------------------------------------
								double var2part=ct->HistHist_sub[i]-pow(ct->Hist_sub[i],2.0);
								double std2part =sqrt(var2part/ct->HistHit[i]);
								double var4part=ct->HistHist[i]-pow(ct->Hist[i],2.0);
								double std4part =sqrt(var4part/ct->HistHit[i]);

								//Get error  delta c2{4}
								//-------------------------
								double err = sqrt(std4part*std4part + 16.0*ct->Hist_sub[i]*ct->Hist_sub[i]*std2part*std2part);
								ct->HistErr[i]=err;
								ct->FinalHist[i]=c24;
							}

						}else{

							//c2{4} (3subevent)
							//---------
							for(int i=0; i<ct->max_nx+1; ++i){
								//Get c2{4} = <<4>>aa|bc - 2 * <<2>>a|b * <<2>>a|c
								//----------------------------------------------------
								double c24=ct->Hist[i] - 2.0 * ct->Hist_sub[i]*ct->Hist_subsub[i];

								//Prepare standard error of <<2>>a|b, <<2>>a|c, and <<4>>aa|bb.
								//----------------------------------------
								double var2partAB=ct->HistHist_sub[i]-pow(ct->Hist_sub[i],2.0);
								double std2partAB =sqrt(var2partAB/ct->HistHit[i]);
								double var2partAC=ct->HistHist_subsub[i]-pow(ct->Hist_subsub[i],2.0);
								double std2partAC =sqrt(var2partAC/ct->HistHit[i]);
								double var4part=ct->HistHist[i]-pow(ct->Hist[i],2.0);
								double std4part =sqrt(var4part/ct->HistHit[i]);

								//Get error  delta c2{4}
								//-------------------------
								double err = sqrt(std4part*std4part 
										+ 4.0*pow(ct->Hist_subsub[i],2)*pow(std2partAB,2)
										+ 4.0*pow(ct->Hist_sub[i],2)*pow(std2partAC,2));
								ct->HistErr[i]=err;
								ct->FinalHist[i]=c24;
							}


						}


					}


					//FOURIER COEFFICIENT
					//----------------------

					if(!options.get_flag__4particle()){

						//v2{2}
						//---------

						for(int i=0; i<ct->max_nx+1; ++i){
							//Obtain vn{2} = sqrt(cn{2})
							//--------------------------- 
							double v22 =sqrt(ct->Hist[i]);

							//Prepare standard error of <<2>>.
							//----------------------------------------
							double var2part=ct->HistHist[i]-pow(ct->Hist[i],2.0);
							double std2part =sqrt(var2part/ct->HistHit[i]);

							//Get error  delta v2{2}
							//-------------------------
							double err=(1.0/2.0)*(std2part/sqrt(v22));
							ct->HistErr_vn[i]=err;
							ct->FinalHist_vn[i]=v22;
						}




					}else{


						if(!options.get_flag_3subevent()){

							//v2{4} (standard., 2-subevent)
							//---------
							for(int i=0; i<ct->max_nx+1; ++i){
								//Get c2{4} = <<4>>-2*<<2>>
								//------------------------
								double c24=ct->Hist[i] - 2.0 * pow(ct->Hist_sub[i],2);

								//Prepare standard error of <<2>> and <<4>>.
								//----------------------------------------
								double var2part=ct->HistHist_sub[i]-pow(ct->Hist_sub[i],2.0);
								double std2part =sqrt(var2part/ct->HistHit[i]);
								double var4part=ct->HistHist[i]-pow(ct->Hist[i],2.0);
								double std4part =sqrt(var4part/ct->HistHit[i]);

								//Get error  delta c2{4}
								//-------------------------
								double errc24 = sqrt(std4part*std4part + 16.0*ct->Hist_sub[i]*ct->Hist_sub[i]*std2part*std2part);

								//Get v2{4} = (-c2{4})**(1/4)
								//---------------------------
								double v24 = pow(-c24, 1.0/4.0);

								//Get error delta v2{4}.
								//----------------------------
								double err =  (1.0/4.0)*errc24*pow(c24, -3.0/4.0);
								ct->FinalHist_vn[i]=v24;
								ct->HistErr_vn[i]=err;
							}

						}else{


							//v2{4} (3subevent)
							//---------
							for(int i=0; i<ct->max_nx+1; ++i){

								//Get c2{4} = <<4>>aa|bc - 2 * <<2>>a|b * <<2>>a|c
								//----------------------------------------------------
								double c24=ct->Hist[i] - 2.0 * ct->Hist_sub[i]*ct->Hist_subsub[i];

								//Prepare standard error of <<2>>a|b, <<2>>a|c, and <<4>>aa|bb.
								//----------------------------------------
								double var2partAB=ct->HistHist_sub[i]-pow(ct->Hist_sub[i],2.0);
								double std2partAB =sqrt(var2partAB/ct->HistHit[i]);
								double var2partAC=ct->HistHist_subsub[i]-pow(ct->Hist_subsub[i],2.0);
								double std2partAC =sqrt(var2partAC/ct->HistHit[i]);
								double var4part=ct->HistHist[i]-pow(ct->Hist[i],2.0);
								double std4part =sqrt(var4part/ct->HistHit[i]);

								//Get error  delta c2{4}
								//-------------------------
								double errc24 = sqrt(std4part*std4part 
										+ 4.0*pow(ct->Hist_subsub[i],2)*pow(std2partAB,2)
										+ 4.0*pow(ct->Hist_sub[i],2)*pow(std2partAC,2));

								//Get v2{4} = (-c2{4})**(1/4)
								//---------------------------
								double v24 = pow(-c24, 1.0/4.0);

								//Get error delta v2{4}.
								//----------------------------
								double err =  (1.0/4.0)*errc24*pow(c24, -3.0/4.0);
								ct->FinalHist_vn[i]=v24;
								ct->HistErr_vn[i]=err;
							}



						}
					}




			}





			void fill(shared_ptr<Container>& ct){


				Container::EventInfo& EVENT= ct->EVENTINFO;

				//Determine xbin
				//---------------
				double x_val=constants::dummy;
				if(constants::MODE.find("meanpt")!=string::npos || constants::MODE.find("meanmt")!=string::npos){
					x_val=EVENT.Nch();
					if(x_val<constants::x_min || x_val>this->x_max) return;
				}
				int nx=(!options.get_flag_HI())? this->get_cell_index(x_val):this->get_cell_index_logplot(x_val);


				//Count particle by particle.
				//----------------------------
				for(int j=0; j<(int)EVENT.part.size(); ++j){
					double y_val=constants::dummy;
					if(constants::MODE.find("meanpt")!=string::npos){
						y_val = EVENT.part[j].pt;
					}else if(constants::MODE.find("meanmt")!=string::npos){
						y_val = EVENT.part[j].mt - EVENT.part[j].m;
					}else if(constants::MODE.find("mtscaling")!=string::npos){
						x_val=EVENT.part[j].m;
						if(x_val<constants::x_min || x_val>this->x_max) continue;
						nx=(int)((x_val/this->d_x)+(std::fabs(constants::x_min)/this->d_x));

			 			if(EVENT.part[j].id==constants::id_phi){
							uf->checkMassOnShell(EVENT.part[j].m, EVENT.part[j].e, EVENT.part[j].px, EVENT.part[j].py, EVENT.part[j].pz);
						}

						if(!this->fix_ax(EVENT.part[j].id, nx, x_val)) continue;
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


			bool fix_ax(const int id, int &nx, double m){

					//Current 
					//----------
					if(abs(id)==constants::id_ch_pion) nx=0;
					else if(abs(id)==constants::id_ch_kaon) nx=1;
					else if(abs(id)==constants::id_proton) nx=2;
					else if(abs(id)==constants::id_phi) nx=3;
					else if(abs(id)==constants::id_lambda) nx=4;
					else if(abs(id)==constants::id_cascade) nx=5;
					else if(abs(id)==constants::id_omega) nx=6;
					else {
						cout << "continue " << endl;
						return false;
					}
					return true;
			}


			int get_cell_index(const double x_val_){
				int ncell = (int)((x_val_/this->d_x)+(fabs(constants::x_min)/this->d_x));
				if(constants::MODE.find("cumulant_multi")!=string::npos){
					if(x_val_>constants::maxNchPP) ncell=6;
				}
				return  ncell;
			}



			int get_cell_index_logplot(const double x_val_){

				double x_val=x_val_;
				int ncell=0;


				if(x_val<constants::switchBin_x){
					int ncell=(int) floor((x_val-constants::x_min)/constants::binSize_small);

					//Put small multiplicity events into one bin.
					//---------------------------------------------
					if(constants::MODE.find("cumulant_multi")!=string::npos && x_val<constants::minNchHI) ncell=0;

					return ncell;
				}else{
					//Put small multiplicity events into one bin.
					//---------------------------------------------
					if(constants::MODE.find("cumulant_multi")!=string::npos && x_val<constants::minNchHI){
						ncell=0;
						return ncell;
					}

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



			void fill_vn4multi(shared_ptr<Container>& ct){


				Container::EventInfo& EVENT= ct->EVENTINFO;

				if((int)EVENT.part.size()<4) return;

				//Determine xbin
				//---------------
				double x_val=EVENT.Nch();
				if(x_val<constants::x_min || x_val>this->x_max) return;
				int nx=(!options.get_flag_HI())? this->get_cell_index(x_val):this->get_cell_index_logplot(x_val);

				//Count particle by particle.
				//----------------------------
				std::complex<double> Qvec_tot=constants::initialval_comp;
				std::complex<double> Qvec2_tot=constants::initialval_comp;
				std::complex<double> n_coeff (2.0, 0.0);
				std::complex<double> n2_coeff (2.0*2.0, 0.0);
				for(int j=0; j<(int)EVENT.part.size(); ++j){
					std::complex<double> phi_ (EVENT.part[j].phi,0.0);
					std::complex<double> Qvec=exp(constants::i_img*n_coeff*phi_);
					std::complex<double> Qvec2=exp(constants::i_img*n2_coeff*phi_);
					Qvec_tot += Qvec;
					Qvec2_tot += Qvec2;
				}

				double totN = (double)EVENT.part.size();
				double corr4 = real(Qvec_tot*Qvec_tot*conj(Qvec_tot)*conj(Qvec_tot) 
						- Qvec_tot * conj(Qvec_tot)*(4.0*totN-8.0)
						-conj(Qvec2_tot)* pow(Qvec_tot,2) - Qvec2_tot* conj(pow(Qvec_tot,2)) 
						+  Qvec2_tot*conj(Qvec2_tot) - 6.0*totN +2.0*pow(totN,2))/(pow(totN,4)-6.0*pow(totN,3)+11.0*pow(totN,2)-6.0*totN);

				double corr4_img = imag(Qvec_tot*Qvec_tot*conj(Qvec_tot)*conj(Qvec_tot) 
						- Qvec_tot * conj(Qvec_tot)*(4.0*totN-8.0)
						-conj(Qvec2_tot)* pow(Qvec_tot,2) - Qvec2_tot* conj(pow(Qvec_tot,2)) 
						+  Qvec2_tot*conj(Qvec2_tot) - 6.0*totN +2.0*pow(totN,2))/(pow(totN,4)-6.0*pow(totN,3)+11.0*pow(totN,2)-6.0*totN);

				//Obtain 2particle correlation
				//-------------------------------
				double squared_Qvec = real(Qvec_tot * conj(Qvec_tot));
				double corr2 = (squared_Qvec-totN)/(totN*(totN-1.0));


				ct->Hist[nx]+=corr4*EVENT.weight();
				ct->Hist_sub[nx]+=corr2*EVENT.weight();
				ct->Hist_x[nx]+=x_val*EVENT.weight();
				ct->HistHit[nx]++;
				ct->HistHist[nx]+=corr4*corr4*EVENT.weight();
				ct->HistHist_sub[nx]+=corr2*corr2*EVENT.weight();
				ct->Hist_weight[nx]+=EVENT.weight();
				ct->Hist_img_Qvec[nx]+=corr4_img*EVENT.weight();
				if(ct->max_nx<nx) ct->max_nx=nx;

				ct->SumWeight+=EVENT.weight();

			}



			void fill_vnmulti(shared_ptr<Container>& ct){


				Container::EventInfo& EVENT= ct->EVENTINFO;

				if((int)EVENT.part.size()<2) return;

				//Determine xbin
				//---------------
				double x_val=EVENT.Nch();
				if(x_val<constants::x_min || x_val>this->x_max) return;
				int nx=(!options.get_flag_HI())? this->get_cell_index(x_val):this->get_cell_index_logplot(x_val);

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


			void fill_vnpt_2sub(shared_ptr<Container>& ct){


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





			void fill_vnmulti_2sub(shared_ptr<Container>& ct){


				Container::EventInfo& EVENT= ct->EVENTINFO;


				//Determine xbin
				//---------------
				double x_val=EVENT.Nch();
				if(x_val<constants::x_min || x_val>this->x_max) return;
				int nx=(!options.get_flag_HI())? this->get_cell_index(x_val):this->get_cell_index_logplot(x_val);

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


			void fill_vn4multi_2sub(shared_ptr<Container>& ct){


				Container::EventInfo& EVENT= ct->EVENTINFO;


				//Determine xbin
				//---------------
				double x_val=EVENT.Nch();
				if(x_val<constants::x_min || x_val>this->x_max) return;
				int nx=(!options.get_flag_HI())? this->get_cell_index(x_val):this->get_cell_index_logplot(x_val);

				//Count particle by particle.
				//----------------------------
				std::complex<double> Qvec_tot_A=constants::initialval_comp;
				std::complex<double> Qvec_tot_B=constants::initialval_comp;
				std::complex<double> Qvec2_tot_A=constants::initialval_comp;
				std::complex<double> Qvec2_tot_B=constants::initialval_comp;
				std::complex<double> n_coeff (2.0, 0.0);
				std::complex<double> n2_coeff (2.0*2.0, 0.0);
				double hit_A=0.0, hit_B=0.0;
				for(int j=0; j<(int)EVENT.part.size(); ++j){
					std::complex<double> phi_ (EVENT.part[j].phi,0.0);
					std::complex<double> Qvec=exp(constants::i_img*n_coeff*phi_);
					std::complex<double> Qvec2=exp(constants::i_img*n2_coeff*phi_);
					if(EVENT.part[j].eta<-0.0){
						Qvec_tot_A += Qvec;
						Qvec2_tot_A += Qvec2;
						hit_A++;
					}else if(EVENT.part[j].eta>0.0){
						Qvec_tot_B += Qvec;
						Qvec2_tot_B += Qvec2;
						hit_B++;
					}
				}
				double squared_Qvec = real((pow(Qvec_tot_A,2) -  Qvec2_tot_A)* conj(pow(Qvec_tot_B,2) - Qvec2_tot_B));
				double squared_Qvec_img = imag((pow(Qvec_tot_A,2) -  Qvec2_tot_A)* conj(pow(Qvec_tot_B,2) - Qvec2_tot_B));

				if(hit_A<=1.0 || hit_B<=1.0) return;

				double corr4 = (squared_Qvec)/(hit_A*(hit_A-1)*hit_B*(hit_B-1));

				//Obtain 2particle correlation
				//------------------------------
				double squared_Qvec_2part = real(Qvec_tot_A * conj(Qvec_tot_B));
				double corr2 = (squared_Qvec_2part)/(hit_A*hit_B);


				ct->Hist[nx]+=corr4*EVENT.weight();
				ct->Hist_sub[nx]+=corr2*EVENT.weight();
				ct->Hist_x[nx]+=x_val*EVENT.weight();
				ct->HistHit[nx]++;
				ct->HistHist[nx]+=corr4*corr4*EVENT.weight();
				ct->HistHist_sub[nx]+=corr2*corr2*EVENT.weight();
				ct->Hist_weight[nx]+=EVENT.weight();
				ct->Hist_img_Qvec[nx]+=squared_Qvec_img*EVENT.weight();
				if(ct->max_nx<nx) ct->max_nx=nx;

				ct->SumWeight+=EVENT.weight();

			}



			//----------3sub

			void fill_vn4multi_3sub(shared_ptr<Container>& ct){


				Container::EventInfo& EVENT= ct->EVENTINFO;


				//Determine xbin
				//---------------
				double x_val=EVENT.Nch();
				if(x_val<constants::x_min || x_val>this->x_max) return;
				int nx=(!options.get_flag_HI())? this->get_cell_index(x_val):this->get_cell_index_logplot(x_val);

				//Count particle by particle.
				//----------------------------
				std::complex<double> Qvec_tot_A=constants::initialval_comp;
				std::complex<double> Qvec_tot_B=constants::initialval_comp;
				std::complex<double> Qvec_tot_C=constants::initialval_comp;
				std::complex<double> Qvec2_tot_A=constants::initialval_comp;
				std::complex<double> Qvec2_tot_B=constants::initialval_comp;
				std::complex<double> Qvec2_tot_C=constants::initialval_comp;
				std::complex<double> n_coeff (2.0, 0.0);
				std::complex<double> n2_coeff (2.0*2.0, 0.0);
				double hit_A=0.0, hit_B=0.0, hit_C=0.0;
				for(int j=0; j<(int)EVENT.part.size(); ++j){
					std::complex<double> phi_ (EVENT.part[j].phi,0.0);
					std::complex<double> Qvec=exp(constants::i_img*n_coeff*phi_);
					std::complex<double> Qvec2=exp(constants::i_img*n2_coeff*phi_);
					if(EVENT.part[j].eta<constants::etaA_3sub){
						Qvec_tot_B += Qvec;
						Qvec2_tot_B += Qvec2;
						hit_B++;
					}else if(EVENT.part[j].eta<=constants::etaB_3sub){
						Qvec_tot_A += Qvec;
						Qvec2_tot_A += Qvec2;
						hit_A++;
					}else{
						Qvec_tot_C += Qvec;
						Qvec2_tot_C += Qvec2;
						hit_C++;
					}
				}
				double squared_Qvec = real((pow(Qvec_tot_A,2)-Qvec2_tot_A)*conj(Qvec_tot_B)*conj(Qvec_tot_C));
				double squared_Qvec_img = imag((pow(Qvec_tot_A,2)-Qvec2_tot_A)*conj(Qvec_tot_B)*conj(Qvec_tot_C));

				if(hit_A<=1.0 || hit_B<=0.0|| hit_C<=0.0) return;

				double corr4 = (squared_Qvec)/(hit_A*(hit_A-1)*hit_B*hit_C);

				//Obtain 2particle correlation btw A and B.
				//----------------------------------------
				double squared_Qvec_2partAB = real(Qvec_tot_A * conj(Qvec_tot_B));
				double corr2_AB = (squared_Qvec_2partAB)/(hit_A*hit_B);


				//Obtain 2particle correlation btw A and C.
				//----------------------------------------
				double squared_Qvec_2partAC = real(Qvec_tot_A * conj(Qvec_tot_C));
				double corr2_AC = (squared_Qvec_2partAC)/(hit_A*hit_C);


				ct->Hist[nx]+=corr4*EVENT.weight();
				ct->Hist_sub[nx]+=corr2_AB*EVENT.weight();
				ct->Hist_subsub[nx]+=corr2_AC*EVENT.weight();
				ct->Hist_x[nx]+=x_val*EVENT.weight();
				ct->HistHit[nx]++;
				ct->HistHist[nx]+=corr4*corr4*EVENT.weight();
				ct->HistHist_sub[nx]+=corr2_AB*corr2_AB*EVENT.weight();
				ct->HistHist_subsub[nx]+=corr2_AC*corr2_AC*EVENT.weight();
				ct->Hist_weight[nx]+=EVENT.weight();
				ct->Hist_img_Qvec[nx]+=squared_Qvec_img*EVENT.weight();
				if(ct->max_nx<nx) ct->max_nx=nx;

				ct->SumWeight+=EVENT.weight();

			}












			bool write(const std::string& fname, const shared_ptr<Container>& ct){
				ofstream ofs;
				ofs.open((fname+"/"+constants::default_out_fname).c_str());
				cout <<"write" << endl;
				if(!ofs){ms->open(fname+"/"+constants::default_out_fname); return false;}

				ct->max_nx+=constants::margin;

				if(constants::MODE.find("cumulant_multi")!=string::npos) {

					for(int i=0; i<ct->max_nx; ++i){

						if(ct->HistHit[i]==0) continue;

						double x_axis =ct->Hist_x[i];
						ofs << setw(16) << fixed << setprecision(8) << x_axis << "  "
							<< setw(16) << ct->FinalHist[i] << "  "
							<< setw(16) << ct->HistErr[i] << "  "
							<< setw(16) << ct->FinalHist_vn[i] << "  "
							<< setw(16) << ct->HistErr_vn[i] << "  "
							<< setw(16) << ct->HistHit[i] << endl;
					}

				}else{
					for(int i=0; i<ct->max_nx; ++i){

						if(ct->HistHit[i]==0) continue;

						double x_axis =(constants::MODE.find("cumulant_pt")!=string::npos)? ((constants::x_min+(this->d_x*i))+(constants::x_min+(this->d_x*(i+1))))/2.0: ct->Hist_x[i];
						ofs << setw(16) << fixed << setprecision(8) << x_axis << "  "
							<< setw(16) << ct->FinalHist[i] << "  "
							<< setw(16) << ct->HistErr[i] << "  "
							<< setw(16) << ct->HistHit[i] << endl;
					}
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
							if(options.get_flag_2subevent()){ 
								if(options.get_flag__4particle()){
									this->fill_vn4multi_2sub(ct); 
								}else{		
									this->fill_vnmulti_2sub(ct); 
								}              
							}else if(options.get_flag_3subevent()){
								if(options.get_flag__4particle()){
									this->fill_vn4multi_3sub(ct); 
								}else{		
									cout << "ERROR :( Option is wrong. --threesub should be used with --4particle." << endl;
									exit(1);
								}              
							}else{
								if(options.get_flag__4particle()){
									this->fill_vn4multi(ct);
								}else{
									this->fill_vnmulti(ct);
								}
							}
						}else if(constants::MODE.find("cumulant_pt")!=std::string::npos){
							if(options.get_flag_2subevent()) this->fill_vnpt_2sub(ct); 
							else this->fill_vnpt(ct);
						}else{
							this->fill(ct);
						}

					}//Event loop

					if(constants::MODE.find("cumulant_multi")!=string::npos || constants::MODE.find("cumulant_pt")!=string::npos) this->stat_flow(ct);
					else this->stat(ct);


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
