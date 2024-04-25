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


Stat::Stat(shared_ptr<Message> ms_in, Options options_in, shared_ptr<InfoHist> infohist_in, shared_ptr<Util_func> uf_in):ms(ms_in), options(options_in), infohist(infohist_in), uf(uf_in){};
Stat::~Stat(){};

			void Stat::stat_twopc(shared_ptr<Container>& ct){


				//take average    
				//-------------------------------------
				for(int i=0; i<ct->max_nx+1; ++i){
					for(int j=0; j<ct->max_ny+1; ++j){
						if(!options.get_flag_SB_CMS()){
							ct->Final2DHist[i][j]=ct->Hist2D[i][j]/ct->Hist2DPartHit[i][j]/ct->SumWeight;
							ct->Final2DHistSub[i][j]=ct->HistSub2D[i][j]/ct->Hist2DPartHit[i][j]/ct->SumWeight;
						}else{
							ct->Hist2D[i][j]=ct->Hist2D[i][j]/ct->SumWeight;
						}
						ct->Hist2D_x[i][j]/=ct->Hist2DPartHit[i][j];
						ct->Hist2D_y[i][j]/=ct->Hist2DPartHit[i][j];

						ct->Final2DHit[i][j]=ct->Hist2DPartHit[i][j]/ct->SumWeight;

						//Devide by bin width
						//---------------------------
						if(!options.get_flag_SB_CMS()){
							ct->Final2DHit[i][j]/=this->infohist->d_x;
							ct->Final2DHit[i][j]/=this->infohist->d_y;
						}else{
							ct->Hist2D[i][j]/=this->infohist->d_x;
							ct->Hist2D[i][j]/=this->infohist->d_y;
						}
					}
				}


				if(options.get_flag_SB_CMS()){

					vector<double> X, Y, Z;
					for(int i=0; i<ct->max_nx+1; ++i){
						for(int j=0; j<ct->max_ny+1; ++j){
							ct->HistSub2D[i][j]=ct->HistSub2D[i][j]/ct->SumWeight;
							ct->HistSub2D_x[i][j]/=ct->HistSub2DPartHit[i][j];
							ct->HistSub2D_y[i][j]/=ct->HistSub2DPartHit[i][j];

							//Devide by bin width
							//---------------------------
							ct->HistSub2D[i][j]/=this->infohist->d_x;
							ct->HistSub2D[i][j]/=this->infohist->d_y;

							X.push_back(ct->HistSub2D_x[i][j]);
							Y.push_back(ct->HistSub2D_y[i][j]);
							Z.push_back(ct->HistSub2D[i][j]);
						}
					}


					//Obtain S/B*(B(0,0))
					//============
					double B00 = this->get_B00(X, Y, Z);
					cout << "B00 " << B00 << endl;
					ct->B00 = B00;
					for(int i=0; i<ct->max_nx+1; ++i){
						for(int j=0; j<ct->max_ny+1; ++j){

							ct->Final2DHist[i][j]=(ct->HistSub2D[i][j]>constants::SMALL)? B00*(ct->Hist2D[i][j]/ct->HistSub2D[i][j]):0.0;

						}
					}

				}

return;
			}


double Stat::get_B00(const vector<double>& X, const vector<double>& Y, const vector<double>& Z){

	double dl=constants::LARGE;
	double B00=0.0;
	for(int i=0; i<(int)X.size(); i++){

			if(dl>pow(fabs(X[i]),2)+pow(fabs(Y[i]),2)){
				dl=pow(fabs(X[i]), 2) + pow(fabs(Y[i]), 2);
				B00=Z[i];
			}

	}

return B00;

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
					ct->FinalHist[i]=ct->Hist[i]/ct->SumWeight;
					ct->Hist_x[i]/=ct->Hist_weight[i];

					double meanxx  = ct->HistHist[i]/ct->SumWeight;
					double meanx = ct->FinalHist[i];

					//Get standard error
					//-------------------------------------
					double var=meanxx-pow(meanx,2.0);
					double error=sqrt(var/ct->SumWeight);
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

				auto stc=make_shared<StatCumulant>();

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
					if(!options.get_flag__4particle()){
						stc->qvec2corr[i]=ct->Hist[i];
						stc->numpair2corr[i]=ct->Hist2[i];
					}else{
						if(!options.get_flag_3subevent()){
							stc->qvec2corr[i]=ct->Hist_sub[i];
							stc->numpair2corr[i]=ct->Hist2_sub[i];
							stc->qvec4corr[i]=ct->Hist[i];
							stc->numpair4corr[i]=ct->Hist2[i];
						}else{
							stc->qvec2ABcorr[i]=ct->Hist_sub[i];
							stc->numpair2ABcorr[i]=ct->Hist2_sub[i];
							stc->qvec2ACcorr[i]=ct->Hist_subsub[i];
							stc->numpair2ACcorr[i]=ct->Hist2_subsub[i];
							stc->qvec4corr[i]=ct->Hist[i];
							stc->numpair4corr[i]=ct->Hist2[i];
						}
					}


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
							stc->qvec2corr_err[i]=sqrt(var/ct->HistHit[i]);
							//Denominator. Npair + delta Npair
							//-----------------------------
							double var2=ct->HistHist2[i]-pow(ct->Hist2[i],2.0);
							stc->numpair2corr_err[i]=sqrt(var2/ct->HistHit[i]);

							//Obtain <<2>> + delta <<2>>!
							//----------------------------
							ct->FinalHist[i]=stc->qvec2corr[i]/stc->numpair2corr[i];
							ct->HistErr[i]=(stc->qvec2corr[i]/stc->numpair2corr[i])
								*sqrt(pow(stc->qvec2corr_err[i]/stc->qvec2corr[i], 2.0)
										+pow(stc->numpair2corr_err[i]/stc->numpair2corr[i], 2.0));
							//Store info
							//-------------
							stc->Corr2[i]=ct->FinalHist[i];
							stc->Corr2Err[i]=ct->HistErr[i];
						}


					}else{


						if(!options.get_flag_3subevent()){

							//c2{4} (standard., 2subevent)
							//
							//1. Obtain <<2>>=Q*Q/Npair, <<4>>=Q*Q/Npair and their errors.
							//2. Obtain C2{4} and its error using 1.
							//==================================================================
							for(int i=0; i<ct->max_nx+1; ++i){
								//1. Get <<4>> and <<2>>
								//------------------------
								stc->Corr4[i]=stc->qvec4corr[i]/stc->numpair4corr[i];
								stc->Corr2[i]=stc->qvec2corr[i]/stc->numpair2corr[i];

                                                                //Prepare standard error of Q*Q and Npair of each <2> and <4>
								double var2part=ct->HistHist_sub[i]-pow(ct->Hist_sub[i],2.0);
								stc->qvec2corr_err[i]=sqrt(var2part/ct->HistHit[i]);
								double var2part_2=ct->HistHist2_sub[i]-pow(ct->Hist2_sub[i],2.0);
								stc->numpair2corr_err[i]=sqrt(var2part_2/ct->HistHit[i]);
								double var4part=ct->HistHist[i]-pow(ct->Hist[i],2.0);
								stc->qvec4corr_err[i]=sqrt(var4part/ct->HistHit[i]);
								double var4part_2=ct->HistHist2[i]-pow(ct->Hist2[i],2.0);
								stc->numpair4corr_err[i]=sqrt(var4part_2/ct->HistHit[i]);

								//Obtain <<2>> + delta <<2>>, <<4>> + delta <<4>>!
								stc->Corr2Err[i]=(stc->qvec2corr[i]/stc->numpair2corr[i])
									*sqrt(pow(stc->qvec2corr_err[i]/stc->qvec2corr[i], 2.0)
											+pow(stc->numpair2corr_err[i]/stc->numpair2corr[i], 2.0));
								stc->Corr4Err[i]=(stc->qvec4corr[i]/stc->numpair4corr[i])
									*sqrt(pow(stc->qvec4corr_err[i]/stc->qvec4corr[i], 2.0)
											+pow(stc->numpair4corr_err[i]/stc->numpair4corr[i], 2.0));

								//2. Get c2{4} + delta c2{4}
								//----------------------------
								double c24 = stc->Corr4[i]-2.0*stc->Corr2[i]*stc->Corr2[i];
								double c24_err = sqrt(stc->Corr4Err[i]*stc->Corr4Err[i]
										-8.0*stc->Corr2[i]*stc->Corr2Err[i]*stc->Corr4Err[i]
										+16.0*stc->Corr2[i]*stc->Corr2[i]*stc->Corr2Err[i]*stc->Corr2Err[i]);
								ct->FinalHist[i]=c24;
								ct->HistErr[i]=c24_err;

								//Store info
								//-----------
								stc->Cumu4[i]=ct->FinalHist[i];
								stc->Cumu4Err[i]=ct->HistErr[i];
							}

						}else{

							//c2{4} (3subevent)
							//
							//1. Obtain <<2>>a|b=Q*Q/Npair, <<2>>a|c=Q*Q/Npair, <<4>>a|b,a|c=Q*Q/Npair and their errors.
							//2. Obtain c2{4} and its error using 1.
							//==================================================================
							for(int i=0; i<ct->max_nx+1; ++i){
								//1. <<2>>a|b+delta<<2>>a|b, <<2>>a|c+delta<<2>>a|c, <<4>>+delta<<4>> 
								//-----------------------------------------------------------------------------
								stc->Corr4[i]=stc->qvec4corr[i]/stc->numpair4corr[i];
								stc->Corr2AB[i]=stc->qvec2ABcorr[i]/stc->numpair2ABcorr[i];
								stc->Corr2AC[i]=stc->qvec2ACcorr[i]/stc->numpair2ACcorr[i];

                                                                //Prepare standard error of Q*Q and Npair of each <2>, <2> and <4>
								double var2ABpart=ct->HistHist_sub[i]-pow(ct->Hist_sub[i],2.0);
								stc->qvec2ABcorr_err[i]=sqrt(var2ABpart/ct->HistHit[i]);
								double var2ABpart_2=ct->HistHist2_sub[i]-pow(ct->Hist2_sub[i],2.0);
								stc->numpair2ABcorr_err[i]=sqrt(var2ABpart_2/ct->HistHit[i]);

								double var2ACpart=ct->HistHist_subsub[i]-pow(ct->Hist_subsub[i],2.0);
								stc->qvec2ACcorr_err[i]=sqrt(var2ACpart/ct->HistHit[i]);
								double var2ACpart_2=ct->HistHist2_subsub[i]-pow(ct->Hist2_subsub[i],2.0);
								stc->numpair2ACcorr_err[i]=sqrt(var2ACpart_2/ct->HistHit[i]);

								double var4part=ct->HistHist[i]-pow(ct->Hist[i],2.0);
								stc->qvec4corr_err[i]=sqrt(var4part/ct->HistHit[i]);
								double var4part_2=ct->HistHist2[i]-pow(ct->Hist2[i],2.0);
								stc->numpair4corr_err[i]=sqrt(var4part_2/ct->HistHit[i]);

								//Obtain <<2>> + delta <<2>>, <<4>> + delta <<4>>!
								stc->Corr2ABErr[i]=(stc->qvec2ABcorr[i]/stc->numpair2ABcorr[i])
									*sqrt(pow(stc->qvec2ABcorr_err[i]/stc->qvec2ABcorr[i], 2.0)
											+pow(stc->numpair2ABcorr_err[i]/stc->numpair2ABcorr[i], 2.0));
								stc->Corr2ACErr[i]=(stc->qvec2ACcorr[i]/stc->numpair2ACcorr[i])
									*sqrt(pow(stc->qvec2ACcorr_err[i]/stc->qvec2ACcorr[i], 2.0)
											+pow(stc->numpair2ACcorr_err[i]/stc->numpair2ACcorr[i], 2.0));
								stc->Corr4Err[i]=(stc->qvec4corr[i]/stc->numpair4corr[i])
									*sqrt(pow(stc->qvec4corr_err[i]/stc->qvec4corr[i], 2.0)
											+pow(stc->numpair4corr_err[i]/stc->numpair4corr[i], 2.0));


								//2. Get c2{4} + delta c2{4}
								//----------------------------
								double c24 = stc->Corr4[i]-2.0*stc->Corr2AB[i]*stc->Corr2AC[i];
								double c24_err = sqrt(stc->Corr4Err[i]*stc->Corr4Err[i]
										+4.0*pow(stc->Corr2AB[i],2.0)*pow(stc->Corr2ACErr[i],2.0)
										+4.0*pow(stc->Corr2AC[i],2.0)*pow(stc->Corr2ABErr[i],2.0)
										-4.0*stc->Corr4Err[i]*stc->Corr2ABErr[i]*stc->Corr2AC[i]
										-4.0*stc->Corr4Err[i]*stc->Corr2ACErr[i]*stc->Corr2AB[i]
										+8.0*stc->Corr2ABErr[i]*stc->Corr2ACErr[i]*stc->Corr2AB[i]*stc->Corr2AC[i]);
								ct->FinalHist[i]=c24;
								ct->HistErr[i]=c24_err;

								//Store info
								//-----------
								stc->Cumu4[i]=ct->FinalHist[i];
								stc->Cumu4Err[i]=ct->HistErr[i];

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
							//*** Corr2 = Cumu2***
							//--------------------------- 
							double v22 =sqrt(stc->Corr2[i]);

							//Prepare standard error of <<2>>.
							//----------------------------------------
							double std2part=stc->Corr2Err[i];

							//Get error  delta v2{2}
							//-------------------------
							double err=(1.0/2.0)*(std2part/sqrt(v22));
							ct->FinalHist_vn[i]=v22;
							ct->HistErr_vn[i]=err;
						}




					}else{

							//v2{4} (standard., 2-subevent, 3-subevent)
							//---------
							for(int i=0; i<ct->max_nx+1; ++i){
								//Get v2{4} = (-c2{4})**(1/4)
								//---------------------------
								double v24 = pow(-1.0*stc->Cumu4[i], 1.0/4.0);
								double err =  (-1.0/4.0)*pow(-1.0*stc->Cumu4[i], -3.0/4.0)*stc->Cumu4Err[i];
								ct->FinalHist_vn[i]=v24;
								ct->HistErr_vn[i]=err;
							}

					}




			}




