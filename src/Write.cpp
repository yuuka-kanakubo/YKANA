#include "Write.h"


Write::Write(shared_ptr<Message> ms_in, Options options_in, shared_ptr<Util_func> uf_in):ms(ms_in), options(options_in), uf(uf_in){};
Write::~Write(){};


int Write::getMapEdgeX(const double maxval){
	double half_dx = options.ih.d_x / 2.0;
	double binZeroCenter = 0; 
        int n = std::round((maxval - binZeroCenter + std::fabs(options.ih.x_edge_min) - half_dx) / options.ih.d_x);
	if(n<0 || n>constants::x_cell_capa){
		std::cout << "ERROR:( n is beyond the capacity.  n:" << n << " out of constants::x_cell_capa " << constants::x_cell_capa << std::endl;
		exit(EXIT_FAILURE);
	}
	return n;
}

int Write::getMapEdgeY(const double maxval){
	double half_dy = options.ih.d_y / 2.0;
	double binZeroCenter = 0; 
        int n = std::round((maxval - binZeroCenter + std::fabs(options.ih.y_edge_min) - half_dy) / options.ih.d_y);
	if(n<0 || n>constants::y_cell_capa){
		std::cout << "ERROR:( n is beyond the capacity.  n:" << n << " out of constants::y_cell_capa " << constants::y_cell_capa << std::endl;
		exit(EXIT_FAILURE);
	}
	return n;
}


bool Write::write_BSTR(const std::string& fname, const Container& ct){
	std::ofstream ofs;
	ofs.open((fname+"/"+constants::default_out_fname).c_str());
	if(!ofs){ms->open(fname+"/"+constants::default_out_fname); return false;}

	for(int i=0; i<ct.max_nx; ++i){
		ofs << setw(16) << fixed << setprecision(8) << ct.Hist_x[i] << "  "
			<< setw(16) << ct.FinalHist[i] << "  "
			<< setw(16) << ct.HistErr[i] << "  "
			<< setw(16) << ct.HistHit[i] << std::endl;
	}

	ofs << std::endl;
	return true;
}



			bool Write::write(const std::string& fname, const shared_ptr<Container>& ct){
				std::ofstream ofs;
				ofs.open((fname+"/"+constants::default_out_fname).c_str());
				if(!ofs){ms->open(fname+"/"+constants::default_out_fname); return false;}

				ct->max_nx+=constants::margin;

				if(constants::MODE.find("cumulant_pt")!=std::string::npos || constants::MODE.find("cumulant_eta")!=std::string::npos || constants::MODE.find("cumulant_multi")!=std::string::npos) {

					for(int i=0; i<ct->max_nx+1; ++i){

						if(ct->HistHit[i]==0) continue;
						double x_axis= options.ih.x_edge_min + (i + 0.5) * options.ih.d_x;
						if(options.get_xaxis_type()==1){
							x_axis =ct->Hist_x[i];
						}
						ofs << setw(16) << fixed << setprecision(8) << x_axis << "  "
							<< setw(16) << ct->FinalHist[i] << "  "
							<< setw(16) << ct->HistErr[i] << "  "
							<< setw(16) << ct->FinalHist_vn[i] << "  "
							<< setw(16) << ct->HistErr_vn[i] << "  "
							<< setw(16) << ct->HistHit[i] << std::endl;
					}

				}else if(constants::MODE.find("twopc2D")!=std::string::npos 
						|| constants::MODE.find("twodm")!=std::string::npos) {

					if(!options.get_flag_SB_CMS()){
						for(int i=0; i<this->getMapEdgeX(options.ih.x_edge_max); ++i){
							for(int j=0; j<this->getMapEdgeY(options.ih.y_edge_max); ++j){

								double xaxis= options.ih.x_edge_min + (i + 0.5) * options.ih.d_x;
								double yaxis= options.ih.y_edge_min + (j + 0.5) * options.ih.d_y;
								//ct->Hist2D_x[i][j];
								//ct->Hist2D_y[i][j];
								ofs << fixed << setprecision(8) 
									<< setw(16) << xaxis << "  "
									<< setw(16) << yaxis << "  "
									//<< setw(16) << ct->Hist2D_x[i][j] << "  "
									//<< setw(16) << ct->Hist2D_x[i][j] << "  "
									<< setw(16) << ct->Final2DHist[i][j] << "  "
									<< setw(16) << ct->Final2DHistSub[i][j] << "  "
									<< setw(16) << ct->Final2DHit[i][j] << std::endl;
							}
							ofs << std::endl;
						}
					}else{
						for(int i=0; i<this->getMapEdgeX(options.ih.x_edge_max); ++i){
							for(int j=0; j<this->getMapEdgeY(options.ih.y_edge_max); ++j){

								double xaxis= options.ih.x_edge_min + (i + 0.5) * options.ih.d_x;
								double yaxis= options.ih.y_edge_min + (j + 0.5) * options.ih.d_y;
								//ct->Hist2D_x[i][j];
								//ct->Hist2D_y[i][j];
								ofs << setw(16) << fixed << setprecision(8) << xaxis << "  "
									<< setw(16) << yaxis << "  "
									<< setw(16) << ct->Final2DHist[i][j] << "  "
									<< setw(16) << ct->Hist2D[i][j] << "  "
									<< setw(16) << ct->HistSub2D[i][j] << "  "
									<< setw(16) << ct->Hist2DPartHit[i][j] << std::endl;
							}
							ofs << std::endl;
						}
					}

				}else{
					for(int i=0; i<this->getMapEdgeX(options.ih.x_edge_max); ++i){

						double x_axis= options.ih.x_edge_min + (i + 0.5) * options.ih.d_x;
						if(options.get_xaxis_type()==1){
							x_axis =ct->Hist_x[i];
						}
						ofs << setw(16) << fixed << setprecision(8) << x_axis << "  "
							<< setw(16) << ct->FinalHist[i] << "  "
							<< setw(16) << ct->HistErr[i] << "  "
							<< setw(16) << ct->HistHit[i] << std::endl;
					}
				}
				if(constants::MODE.find("Rt_spectra")!=string::npos){
					ofs << "%Mean Nt:" << ct->meanNt << std::endl;
				}
				ofs << std::endl;
				return true;
			}





			bool Write::write_RtYield(const std::string& fname, const shared_ptr<Container>& ct){
				Container::yield spname;
				for(int sp=0; sp<constants::num_of_Species_Rt; sp++){
					//Make output file
					//-------------------
					std::ofstream ofsTrans;
					ofsTrans.open((fname+"/"+"Trans_"+spname.get_particleName(sp)+constants::default_out_fname).c_str());
					if(!ofsTrans){ms->open(fname+"/"+"Trans_"+spname.yield::get_particleName(sp)+constants::default_out_fname); return false;}

					ct->max_nx+=constants::margin;

					for(int i=0; i<ct->max_nx+1; ++i){
						if(ct->HistHit_Rt[sp][i]==0) continue;
						double x_axis =ct->Hist_x[i];
						ofsTrans << setw(16) << fixed << setprecision(8) << x_axis << "  "
							<< setw(16) << ct->RtHist_RtTrans_yield[sp][i] << "  "
							<< setw(16) << ct->HistErrTrans_Rt[sp][i] << "  "
							<< setw(16) << ct->HistHit_Rt[sp][i] << std::endl;
					}
					ofsTrans << "%Mean Nt:" << ct->meanNt << std::endl;
					ofsTrans.close();

					std::ofstream ofsToward;
					ofsToward.open((fname+"/"+"Toward_"+spname.get_particleName(sp)+constants::default_out_fname).c_str());
					if(!ofsToward){ms->open(fname+"/"+"Toward_"+spname.yield::get_particleName(sp)+constants::default_out_fname); return false;}

					ct->max_nx+=constants::margin;

					for(int i=0; i<ct->max_nx+1; ++i){
						if(ct->HistHit_Rt[sp][i]==0) continue;
						double x_axis =ct->Hist_x[i];
						ofsToward << setw(16) << fixed << setprecision(8) << x_axis << "  "
							<< setw(16) << ct->RtHist_RtToward_yield[sp][i] << "  "
							<< setw(16) << ct->HistErrToward_Rt[sp][i] << "  "
							<< setw(16) << ct->HistHit_Rt[sp][i] << std::endl;
					}
					ofsToward << "%Mean Nt:" << ct->meanNt << std::endl;
					ofsToward.close();
					
				}
				return true;
			}
