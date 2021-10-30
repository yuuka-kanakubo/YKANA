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
#include "InfoHist.h"
#include "Stat.h"

using namespace std;


Stat::Stat(shared_ptr<Message> ms_in, Settings::Options options_in, shared_ptr<InfoHist> infohist_in, shared_ptr<Util_func> uf_in):ms(ms_in), options(options_in), infohist(infohist_in), uf(uf_in){};
Stat::~Stat(){};

			void Stat::stat_twopc(shared_ptr<Container>& ct){

				//take average    
				//-------------------------------------
				for(int i=0; i<ct->max_nx+1; ++i){
					for(int j=0; j<ct->max_ny+1; ++j){
						ct->Final2DHist[i][j]=ct->Hist2D[i][j]/ct->SumPair;
						ct->Hist2D_x[i][j]/=ct->Hist2DPartHit[i][j];
						ct->Hist2D_y[i][j]/=ct->Hist2DPartHit[i][j];

						//Devide by bin width
						//---------------------------
						ct->Final2DHist[i][j]/=this->infohist->d_x;
						ct->Final2DHist[i][j]/=this->infohist->d_y;
					}
				}

			}

			void Stat::stat_twopc1D(shared_ptr<Container>& ct){

				//take average    
				//-------------------------------------
				for(int i=0; i<ct->max_nx+1; ++i){
					if(options.get_flag_tagged()){
						if(constants::MODE.find("twopcInteg")!=string::npos) ct->FinalHist[i]=ct->Hist[i]/ct->HistHit[i];
						else ct->FinalHist[i]=ct->Hist[i]/ct->SumTrig;
						if(options.get_flag__2PCfull()){
							ct->FinalHist[i]/=2.0*constants::DeltaEtaFULL;
							if(constants::MODE.find("twopcInteg")!=string::npos) ct->FinalHist[i]/=2.0*M_PI;
						}else if(options.get_flag__2PCnearside()){
							ct->FinalHist[i]/=2.0*constants::DeltaEtaNS;
							if(constants::MODE.find("twopcInteg")!=string::npos) ct->FinalHist[i]/=2.0*constants::DeltaPhiNS;
						}else if(options.get_flag__2PCout()){
							ct->FinalHist[i]/=2.0*fabs(constants::DeltaEtaFULL-constants::DeltaEtaNS);
							if(constants::MODE.find("twopcInteg")!=string::npos) ct->FinalHist[i]/=fabs(constants::DeltaPhiOUT-constants::DeltaPhiNS);
						}else{
							cout << ":( ERROR Something wrong with the flags. " << endl;
							exit(1);
						}
					}else{
						//ct->FinalHist[i]=ct->Hist[i]/ct->SumPair;
						ct->FinalHist[i]=ct->Hist[i]/ct->SumTrig;
						ct->FinalHist[i]/=this->infohist->d_x;
					}
					ct->Hist_x[i]/=ct->Hist_weight[i];

				}

			}

			void Stat::stat(shared_ptr<Container>& ct){

				//take average    
				//-------------------------------------
				for(int i=0; i<ct->max_nx+1; ++i){
					ct->FinalHist[i]=ct->Hist[i]/ct->Hist_weight[i];
					ct->Hist_x[i]/=ct->Hist_weight[i];
					double meanxx  = ct->HistHist[i]/ct->Hist_weight[i];
					double meanx = ct->FinalHist[i];

					// devide by cell width 
					//-------------------------------------
					//ct->Hist[i]/=this->infohist->d_x;


					//Get standard error
					//-------------------------------------
					double var=meanxx-pow(meanx,2.0);
					double error=sqrt(var/ct->HistHit[i]);
					ct->HistErr[i]=error;
				}


			}

			void Stat::stat_Rt(shared_ptr<Container>& ct){

				//Get <Nt>
				//-------------------------------------
				double Nt_tot=0.0;
				for(int i=0; i<(int)ct->Nt_eBye.size(); ++i){
					Nt_tot+=(double)ct->Nt_eBye[i]*ct->weight_eBye[i];
				}
				double meanNt = Nt_tot/ct->SumWeight;

				for(int i=0; i<(int)ct->Nt_eBye.size(); ++i){


					//Rt
					//----
					double x_val=((double)ct->Nt_eBye[i])/meanNt;
					if(x_val<constants::x_min || x_val>this->infohist->x_max) continue;
					int nx=(int)((x_val/this->infohist->d_x)+(std::fabs(constants::x_min)/this->infohist->d_x));

					//cout << ct->dNdeta_eBye[i]  << "   " << x_val << endl;

					//If the event is weighted with w_i where \Sum_i^ev w_i = w_tot, 
					//The event should be counted as (w_i/w_tot) event.
					//---------------------------------------------------------- 
					ct->Hist[nx]+=1.0*(ct->weight_eBye[i]/ct->SumWeight);
					ct->Hist_x[nx]+=x_val*ct->weight_eBye[i];
					ct->Hist_weight[nx]+=ct->weight_eBye[i];
					ct->HistHit[nx]++;
					if(ct->max_nx<nx) ct->max_nx=nx;

				}

				for(int i=0; i<ct->max_nx+1; ++i){
					ct->FinalHist[i]=ct->Hist[i];
					ct->Hist_x[i]/=ct->Hist_weight[i];

					// devide by cell width 
					//-------------------------------------
					//ct->Hist[i]/=this->infohist->d_x;

					ct->HistErr[i]=0.0;
				}


			}


			int Stat::get_xaxis_RtClass(double xval){
				if(xval<constants::RtBins[0]) return 0;
				else if(xval<constants::RtBins[1]) return 1;
				else if(xval<constants::RtBins[2]) return 2;
				else if(xval<constants::RtBins[3]) return 3;
				else if(xval<constants::RtBins[4]) return 4;
				else if(xval<constants::RtBins[5]) return 5;
				else if(xval<constants::RtBins[6]) return 6;
				else if(xval<constants::RtBins[7]) return 7;
				else if(xval<constants::RtBins[8]) return 8;
				else{
					cout << "ERROR:( Something is wrong! Rt:" << xval << endl; 
					exit(1);
				}
			}


			void Stat::stat_RtYield(shared_ptr<Container>& ct){

				//Get <Nt>
				//-------------------------------------
				double Nt_tot=0.0;
				double Ntmin_tot=0.0;
				double Ntmax_tot=0.0;
				for(int i=0; i<(int)ct->Nt_eBye.size(); ++i){
					Nt_tot+=(double)ct->Nt_eBye[i]*ct->weight_eBye[i];
					Ntmin_tot+=(double)ct->Ntmin_eBye[i]*ct->weight_eBye[i];
					Ntmax_tot+=(double)ct->Ntmax_eBye[i]*ct->weight_eBye[i];
				}
				double meanNt = Nt_tot/ct->SumWeight;
				//double meanNtmin = Ntmin_tot/ct->SumWeight;
				//double meanNtmax = Ntmax_tot/ct->SumWeight;
				ct->meanNt=meanNt;
				if((int)ct->Nt_eBye.size()!=(int)ct->TowardYield_eBye.size() || (int)ct->Nt_eBye.size()!=(int)ct->TransYield_eBye.size()){
					cout << "ERROR:( Something wrong in stat_RtYield." << endl;
					exit(1);
				}

				for(int i=0; i<(int)ct->Nt_eBye.size(); ++i){

					//Rt
					//----
					double x_val=((double)ct->Nt_eBye[i])/meanNt;
					//double Rtmin=((double)ct->Ntmin_eBye[i])/meanNtmin;
					//double Rtmax=((double)ct->Ntmax_eBye[i])/meanNtmax;
					if(x_val<constants::x_min || x_val>this->infohist->x_max) continue;
					//int nx=(int)((x_val/this->infohist->d_x)+(std::fabs(constants::x_min)/this->infohist->d_x));
					int nx=this->get_xaxis_RtClass(x_val);

					//dndeta, RcoreT, RcoreN, Rt, Rtmin, Rtmax 
					//------------------------------------------
					//cout << fixed << setprecision(8) << ct->TagEventNum[i] << "  "  << ct->dNdeta_eBye[i]  << "  " << ct->CoreT_eBye[i] << "   " <<ct->CoreN_eBye[i] << "   " << x_val << "  " << Rtmin << "  " << Rtmax << endl;

					for(int sp=0; sp<constants::num_of_Species_Rt; sp++){
						double y_val_trans = ct->TransYield_eBye[i].get_sp(sp);
						//double y_val_trans = ct->CoreT_eBye[i];
						ct->RtHist_RtTrans_yield[sp][nx] += y_val_trans*ct->weight_eBye[i];
						ct->RtHist_RtTrans_yieldyield[sp][nx] += y_val_trans*y_val_trans*ct->weight_eBye[i];
						double y_val_toward = ct->TowardYield_eBye[i].get_sp(sp);
						//double y_val_toward = ct->CoreN_eBye[i];
						ct->RtHist_RtToward_yield[sp][nx] += y_val_toward*ct->weight_eBye[i];
						ct->RtHist_RtToward_yieldyield[sp][nx] += y_val_toward*y_val_toward*ct->weight_eBye[i];
						ct->HistHit_Rt[sp][nx]++;
					}
					ct->Hist_x[nx]+=x_val*ct->weight_eBye[i];
					ct->Hist_weight[nx]+=ct->weight_eBye[i];
					ct->HistHit[nx]++;
					if(ct->max_nx<nx) ct->max_nx=nx;

				}

				for(int i=0; i<ct->max_nx+1; ++i){
					for(int sp=0; sp<constants::num_of_Species_Rt; sp++){
						//TRANS
						//--------
						ct->RtHist_RtTrans_yield[sp][i]/=ct->Hist_weight[i];

						double meanxx_trans  = ct->RtHist_RtTrans_yieldyield[sp][i]/ct->Hist_weight[i];
						double meanx_trans = ct->RtHist_RtTrans_yield[sp][i];

						double var_trans=meanxx_trans-pow(meanx_trans,2.0);
						double error_trans=sqrt(var_trans/ct->HistHit_Rt[sp][i]);

						ct->HistErrTrans_Rt[sp][i]=error_trans;

						//TOWARD
						//--------
						ct->RtHist_RtToward_yield[sp][i]/=ct->Hist_weight[i];

						double meanxx_toward  = ct->RtHist_RtToward_yieldyield[sp][i]/ct->Hist_weight[i];
						double meanx_toward = ct->RtHist_RtToward_yield[sp][i];

						double var_toward=meanxx_toward-pow(meanx_toward,2.0);
						double error_toward=sqrt(var_toward/ct->HistHit_Rt[sp][i]);
						ct->HistErrToward_Rt[sp][i]=error_toward;
					}
					ct->Hist_x[i]/=ct->Hist_weight[i];


					// devide by cell width 
					//-------------------------------------
					//ct->Hist[i]/=this->infohist->d_x;


				}


			}

			void Stat::stat_jets(shared_ptr<Container>& ct){


				//take average    
				//-------------------------------------
				for(int i=0; i<ct->max_nx+1; ++i){
					ct->FinalHist[i]=ct->Hist[i]/ct->SumWeight;
					ct->Hist_x[i]/=ct->Hist_weight[i];

					// devide by cell width 
					//-------------------------------------
					ct->FinalHist[i]/=this->infohist->d_x;


					//Get standard error
					//-------------------------------------
					ct->HistErr[i]=0.0;
				}


			}





			void Stat::stat_flow(shared_ptr<Container>& ct){

				StatCumulant stc;

				//take average    
				//-------------------------------------
				for(int i=0; i<ct->max_nx+1; ++i){
					//Take average of numerator and denominator for each first.
					//Hist2 is dealing with denominator.
					//==========================================================
					ct->Hist[i]/=ct->Hist_weight[i];
					ct->Hist_sub[i]/=ct->Hist_weight[i];
					ct->Hist_subsub[i]/=ct->Hist_weight[i];
					ct->Hist2[i]/=ct->Hist_weight[i];
					ct->Hist2_sub[i]/=ct->Hist_weight[i];
					ct->Hist2_subsub[i]/=ct->Hist_weight[i];
					if(constants::MODE.find("cumulant_eta")!=string::npos || constants::MODE.find("cumulant_pt")!=string::npos){
						//When the xaxis is info of 1 particle.
						//-----------------------------------
						ct->Hist_x[i]/=ct->HistPartHit[i];
					}else{
						//When the xaxis is info of 1 event.
						//-----------------------------------
						ct->Hist_x[i]/=ct->Hist_weight[i];
					}
					ct->HistHist[i]/=ct->Hist_weight[i];
					ct->HistHist_sub[i]/=ct->Hist_weight[i];
					ct->HistHist_subsub[i]/=ct->Hist_weight[i];
					ct->HistHist2[i]/=ct->Hist_weight[i];
					ct->HistHist2_sub[i]/=ct->Hist_weight[i];
					ct->HistHist2_subsub[i]/=ct->Hist_weight[i];
					ct->Hist_img_Qvec[i]/=ct->Hist_weight[i];

					//Need to add for higher correlation.
					stc.qvec2corr=ct->Hist[i];
					stc.numpair2corr=ct->Hist2[i];


					// devide by cell width 
					//-------------------------------------
					//ct->Hist[i]/=this->infohist->d_x;

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
						//
						//1. Obtain error of num and denom.
						//2. Calcualte error of <<2>> from 1.
						//=========================================
						for(int i=0; i<ct->max_nx+1; ++i){
							//Numerator. Q*Q + delta Q*Q
							//-----------------------------
							double var=ct->HistHist[i]-pow(ct->Hist[i],2.0);
							stc.qvec2corr_err=sqrt(var/ct->HistHit[i]);
							//Denominator. Npair + delta Npair
							//-----------------------------
							double var2=ct->HistHist2[i]-pow(ct->Hist2[i],2.0);
							stc.numpair2corr_err=sqrt(var2/ct->HistHit[i]);

							//Obtain <<2>> + delta <<2>>!
							//----------------------------
							ct->FinalHist[i]=stc.qvec2corr/stc.numpair2corr;
							ct->HistErr[i]=(stc.qvec2corr/stc.numpair2corr)
								*sqrt(pow(stc.qvec2corr_err/stc.qvec2corr, 2.0)
										+pow(stc.numpair2corr_err/stc.numpair2corr, 2.0));
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




