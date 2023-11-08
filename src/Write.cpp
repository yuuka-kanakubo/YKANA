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
#include "Write.h"

using namespace std;


Write::Write(shared_ptr<Message> ms_in, Settings::Options options_in, shared_ptr<InfoHist> infohist_in, shared_ptr<Util_func> uf_in):ms(ms_in), options(options_in), infohist(infohist_in), uf(uf_in){};
Write::~Write(){};


		int Write::getMapEdgeX(const double maxval){
			int n=(int)((maxval/this->infohist->d_x)+(std::fabs(constants::x_min)/this->infohist->d_x));
			return n;
		}
		int Write::getMapEdgeY(const double maxval){
			int n=(int)((maxval/this->infohist->d_y)+(std::fabs(constants::y_min)/this->infohist->d_y));
			return n;
		}


			bool Write::write(const std::string& fname, const shared_ptr<Container>& ct){
				ofstream ofs;
				ofs.open((fname+"/"+constants::default_out_fname).c_str());
				cout <<"write" << endl;
				if(!ofs){ms->open(fname+"/"+constants::default_out_fname); return false;}

				ct->max_nx+=constants::margin;

				if(constants::MODE.find("cumulant_pt")!=string::npos || constants::MODE.find("cumulant_eta")!=string::npos || constants::MODE.find("cumulant_multi")!=string::npos) {

					for(int i=0; i<ct->max_nx+1; ++i){

						if(ct->HistHit[i]==0) continue;

						double x_axis =ct->Hist_x[i];
						ofs << setw(16) << fixed << setprecision(8) << x_axis << "  "
							<< setw(16) << ct->FinalHist[i] << "  "
							<< setw(16) << ct->HistErr[i] << "  "
							<< setw(16) << ct->FinalHist_vn[i] << "  "
							<< setw(16) << ct->HistErr_vn[i] << "  "
							<< setw(16) << ct->HistHit[i] << endl;
					}

				}else if(constants::MODE.find("twopc2D")!=string::npos 
						|| constants::MODE.find("twodm")!=string::npos) {

					if(!options.get_flag_SB_CMS()){
						for(int i=0; i<this->getMapEdgeX(this->infohist->x_max); ++i){
							for(int j=0; j<this->getMapEdgeY(this->infohist->y_max); ++j){

								double xaxis=((constants::x_min+(this->infohist->d_x*i))+(constants::x_min+(this->infohist->d_x*(i+1))))/2.0;
								double yaxis=((constants::y_min+(this->infohist->d_y*j))+(constants::y_min+(this->infohist->d_y*(j+1))))/2.0;
								//ct->Hist2D_x[i][j];
								//ct->Hist2D_y[i][j];
								ofs << fixed << setprecision(8) 
									<< setw(16) << xaxis << "  "
									<< setw(16) << yaxis << "  "
									//<< setw(16) << ct->Hist2D_x[i][j] << "  "
									//<< setw(16) << ct->Hist2D_x[i][j] << "  "
									<< setw(16) << ct->Final2DHist[i][j] << "  "
									<< setw(16) << ct->Final2DHistSub[i][j] << "  "
									<< setw(16) << ct->Final2DHit[i][j] << endl;
							}
							ofs << endl;
						}
					}else{
						for(int i=0; i<this->getMapEdgeX(this->infohist->x_max); ++i){
							for(int j=0; j<this->getMapEdgeY(this->infohist->y_max); ++j){

								double xaxis=((constants::x_min+(this->infohist->d_x*i))+(constants::x_min+(this->infohist->d_x*(i+1))))/2.0;
								double yaxis=((constants::y_min+(this->infohist->d_y*j))+(constants::y_min+(this->infohist->d_y*(j+1))))/2.0;
								//ct->Hist2D_x[i][j];
								//ct->Hist2D_y[i][j];
								ofs << setw(16) << fixed << setprecision(8) << xaxis << "  "
									<< setw(16) << yaxis << "  "
									<< setw(16) << ct->Final2DHist[i][j] << "  "
									<< setw(16) << ct->Hist2D[i][j] << "  "
									<< setw(16) << ct->HistSub2D[i][j] << "  "
									<< setw(16) << ct->Hist2DPartHit[i][j] << endl;
							}
							ofs << endl;
						}
					}

				}else{
					for(int i=0; i<this->getMapEdgeX(this->infohist->x_max); ++i){

						if(ct->HistHit[i]==0) continue;
						double x_axis =ct->Hist_x[i];
 	  cout << "(:3 = )3 ? " << __FILE__ << " (" << __LINE__ << ")" << endl;
						ofs << setw(16) << fixed << setprecision(8) << x_axis << "  "
							<< setw(16) << ct->FinalHist[i] << "  "
							<< setw(16) << ct->HistErr[i] << "  "
							<< setw(16) << ct->HistHit[i] << endl;
					}
				}
				if(constants::MODE.find("Rt_spectra")!=string::npos){
					ofs << "%Mean Nt:" << ct->meanNt << endl;
				}
				ofs << endl;
				return true;
			}





			bool Write::write_RtYield(const std::string& fname, const shared_ptr<Container>& ct){
				cout <<"writeRt" << endl;
				Container::yield spname;
				for(int sp=0; sp<constants::num_of_Species_Rt; sp++){
					//Make output file
					//-------------------
					ofstream ofsTrans;
					ofsTrans.open((fname+"/"+"Trans_"+spname.get_particleName(sp)+constants::default_out_fname).c_str());
					if(!ofsTrans){ms->open(fname+"/"+"Trans_"+spname.yield::get_particleName(sp)+constants::default_out_fname); return false;}

					ct->max_nx+=constants::margin;

					for(int i=0; i<ct->max_nx+1; ++i){
						if(ct->HistHit_Rt[sp][i]==0) continue;
						double x_axis =ct->Hist_x[i];
						ofsTrans << setw(16) << fixed << setprecision(8) << x_axis << "  "
							<< setw(16) << ct->RtHist_RtTrans_yield[sp][i] << "  "
							<< setw(16) << ct->HistErrTrans_Rt[sp][i] << "  "
							<< setw(16) << ct->HistHit_Rt[sp][i] << endl;
					}
					ofsTrans << "%Mean Nt:" << ct->meanNt << endl;
					ofsTrans.close();

					ofstream ofsToward;
					ofsToward.open((fname+"/"+"Toward_"+spname.get_particleName(sp)+constants::default_out_fname).c_str());
					if(!ofsToward){ms->open(fname+"/"+"Toward_"+spname.yield::get_particleName(sp)+constants::default_out_fname); return false;}

					ct->max_nx+=constants::margin;

					for(int i=0; i<ct->max_nx+1; ++i){
						if(ct->HistHit_Rt[sp][i]==0) continue;
						double x_axis =ct->Hist_x[i];
						ofsToward << setw(16) << fixed << setprecision(8) << x_axis << "  "
							<< setw(16) << ct->RtHist_RtToward_yield[sp][i] << "  "
							<< setw(16) << ct->HistErrToward_Rt[sp][i] << "  "
							<< setw(16) << ct->HistHit_Rt[sp][i] << endl;
					}
					ofsToward << "%Mean Nt:" << ct->meanNt << endl;
					ofsToward.close();
					
				}
				return true;
			}
