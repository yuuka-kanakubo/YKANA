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
#include "Fill.h"

using namespace std;


Fill::Fill(shared_ptr<Message> ms_in, Options options_in, shared_ptr<Util_func> uf_in):N_iCentEv(0), ms(ms_in), options(options_in), uf(uf_in), rndom(nullptr){

	if(options.get_xaxis_type()==3){
		this->SetCustomBin();
	}

};
Fill::~Fill(){};

void Fill::nextCent(shared_ptr<Rndom> rndom_in){
        this->rndom=rndom_in;
	this->N_iCentEv=this->rndom->get_iEv_Cent().size();
	if(this->N_iCentEv<constants::N_RndomEv){
		cout << ":o Warning. This centrality contains only  " << this->N_iCentEv << " events. The background might be not random enough." << endl;
	}
	//ms->print(this->N_iCentEv);
}

void Fill::fill_jets(shared_ptr<Container>& ct){


				Container::EventInfo& EVENT= ct->EVENTINFO;

				//Determine xbin
				//---------------
				double x_val=EVENT.Aj();
				double y_val=1.0;
				if(x_val<options.ih.x_edge_min || x_val>options.ih.x_edge_max) return;
				int nx=this->get_cell_index(x_val);

				ct->Hist[nx]+=y_val;
				ct->Hist_x[nx]+=x_val;
				ct->Hist_weight[nx]+=1.0;
				if(ct->max_nx<nx) ct->max_nx=nx;

				ct->SumWeight+=1.0;

			}





			void Fill::fill_2D(shared_ptr<Container>& ct){


				Container::EventInfo& EVENT= ct->EVENTINFO;

				//Count particle by particle.
				//----------------------------
				double **Hit1ev, **zVal1ev, **zSubVal1ev;
				Hit1ev = new double *[constants::x_cell_capa];
				zVal1ev = new double *[constants::x_cell_capa];
				zSubVal1ev = new double *[constants::x_cell_capa];
				for(int i_cell=0; i_cell<constants::x_cell_capa; i_cell++){
					Hit1ev[i_cell] = new double [constants::y_cell_capa];
					zVal1ev[i_cell] = new double [constants::y_cell_capa];
					zSubVal1ev[i_cell] = new double [constants::y_cell_capa];
				}
				for(int i=0; i<constants::x_cell_capa; i++){
					for(int j=0; j<constants::y_cell_capa; j++){
						Hit1ev[i][j]=0.0;
						zVal1ev[i][j]=0.0;
						zSubVal1ev[i][j]=0.0;
					}
				}
				int max_nx = 0, max_ny = 0;
				for(int j=0; j<(int)EVENT.part.size(); ++j){

					double x_val=EVENT.part[j].x;
					if(x_val<options.ih.x_edge_min || x_val>options.ih.x_edge_max) continue;
					int nx=this->get_cell_index(x_val);
					double y_val=EVENT.part[j].y;
					if(y_val<options.ih.y_edge_min || y_val>options.ih.y_edge_max) continue;
					int ny=this->get_cell_index_y(y_val);

					if(max_nx<nx) max_nx = nx;
					if(max_ny<ny) max_ny = ny;

					zVal1ev[nx][ny]+=EVENT.part[j].tata;
					zSubVal1ev[nx][ny]+=EVENT.part[j].mt;
					Hit1ev[nx][ny]+=1;
					ct->Hist2D_x[nx][ny]+=x_val*EVENT.weight();
					ct->Hist2D_y[nx][ny]+=y_val*EVENT.weight();
					if(ct->max_nx<nx) ct->max_nx=nx;
					if(ct->max_ny<ny) ct->max_ny=ny;

				}

				for(int nx = 0; nx<max_nx+1; nx++){
					for(int ny = 0; ny<max_ny+1; ny++){
						if(fabs(Hit1ev[nx][ny])<constants::SMALL){
							ct->Hist2D[nx][ny]+=0.0;
						}else{
							ct->Hist2D[nx][ny]+=zVal1ev[nx][ny]*EVENT.weight();
							ct->HistSub2D[nx][ny]+=zSubVal1ev[nx][ny]*EVENT.weight();
						}
						ct->Hist2DPartHit[nx][ny]+=Hit1ev[nx][ny]*EVENT.weight();
					}
				}

				ct->SumWeight+=EVENT.weight();

				for(int i = 0; i < constants::x_cell_capa; i++) {
					delete[] Hit1ev[i];
					delete[] zVal1ev[i];
					delete[] zSubVal1ev[i];
				}
				delete[] Hit1ev;
				delete[] zVal1ev;
				delete[] zSubVal1ev;

				return;
			}



			void Fill::fill_twopc(shared_ptr<Container>& ct){


				Container::EventInfo& EVENT= ct->EVENTINFO;

				//Count particle by particle.
				//----------------------------
				double **Hit1ev;
				Hit1ev = new double *[constants::x_cell_capa];
				for(int i_cell=0; i_cell<constants::x_cell_capa; i_cell++){
					Hit1ev[i_cell] = new double [constants::y_cell_capa];
				}
				for(int i=0; i<constants::x_cell_capa; i++){
					for(int j=0; j<constants::y_cell_capa; j++){
						Hit1ev[i][j]=0.0;
					}
				}
				int max_nx = 0, max_ny = 0;
				int NumPair=0;
				int NumTrig=0;
				for(int i=0; i<(int)EVENT.part.size(); ++i){

					if(!(EVENT.part[i].pt>constants::trig_ptmin && EVENT.part[i].pt<constants::trig_ptmax)) continue;
					if(options.get_flag_only_core_triggers() && EVENT.part[i].TAG==constants::corona_tag) continue;
					if(options.get_flag_only_corona_triggers() && EVENT.part[i].TAG==constants::core_tag) continue;
					NumTrig++;

					for(int j=0; j<(int)EVENT.part.size(); ++j){

						if (i==j) continue;
						if(!(EVENT.part[j].pt>constants::assoc_ptmin && EVENT.part[j].pt<constants::assoc_ptmax)) continue;
						if(options.get_flag_only_core_associates() && EVENT.part[j].TAG==constants::corona_tag) continue;
						if(options.get_flag_only_corona_associates() && EVENT.part[j].TAG==constants::core_tag) continue;

						double x_val=EVENT.part[i].eta - EVENT.part[j].eta;
						if(x_val<options.ih.x_edge_min || x_val>options.ih.x_edge_max) continue;
						int nx=this->get_cell_index(x_val);

						double y_val=this->getDeltaPhi(EVENT.part[i].phi, EVENT.part[j].phi);
						if(y_val<options.ih.y_edge_min || y_val>options.ih.y_edge_max) continue;
						int ny=this->get_cell_index_y(y_val);

						if(max_nx<nx) max_nx = nx;
						if(max_ny<ny) max_ny = ny;
						Hit1ev[nx][ny]+=1.0;
						ct->Hist2D_x[nx][ny]+=x_val*EVENT.weight();
						ct->Hist2D_y[nx][ny]+=y_val*EVENT.weight();
						if(ct->max_nx<nx) ct->max_nx=nx;
						if(ct->max_ny<ny) ct->max_ny=ny;
						NumPair++;

					}
				}
				//---------------
				
				for(int nx = 0; nx<max_nx+1; nx++){
					for(int ny = 0; ny<max_ny+1; ny++){
						ct->Hist2D[nx][ny]+=Hit1ev[nx][ny]*EVENT.weight()/((double)NumTrig);
						ct->Hist2DPartHit[nx][ny]+=Hit1ev[nx][ny]*EVENT.weight();
					}
				}

				ct->SumWeight+=EVENT.weight();


				for(int i = 0; i < constants::x_cell_capa; i++) {
					delete[] Hit1ev[i];
				}
				delete[] Hit1ev;

			}


void Fill::fill_twopc_B_CMS(shared_ptr<Container>& ct, const vector<EbyeInfo>& eBye_All){

				Container::EventInfo& EVENT= ct->EVENTINFO;

				//Count particle by particle.
				//----------------------------
				double **Hit1ev;
				Hit1ev = new double *[constants::x_cell_capa];
				for(int i_cell=0; i_cell<constants::x_cell_capa; i_cell++){
					Hit1ev[i_cell] = new double [constants::y_cell_capa];
				}
				for(int i=0; i<constants::x_cell_capa; i++){
					for(int j=0; j<constants::y_cell_capa; j++){
						Hit1ev[i][j]=0.0;
					}
				}
				int max_nx = 0, max_ny = 0;
				int NumPair=0;
				int NumTrig=0;
				for(int i=0; i<(int)EVENT.part.size(); ++i){

					if(!(EVENT.part[i].pt>constants::trig_ptmin && EVENT.part[i].pt<constants::trig_ptmax)) continue;
					if(options.get_flag_only_core_triggers() && EVENT.part[i].TAG==constants::corona_tag) continue;
					if(options.get_flag_only_corona_triggers() && EVENT.part[i].TAG==constants::core_tag) continue;
					NumTrig++;

					//Select N_RndomEv
					//==================
					vector <Container::ParticleInfo> Mixedpart = this->select_N_RndomEv(eBye_All);

					for(int j=0; j<(int)Mixedpart.size(); ++j){

						string TAG = Mixedpart[j].TAG;
						if(options.get_flag_only_core_associates() && TAG==constants::corona_tag) continue;
						if(options.get_flag_only_corona_associates() && TAG==constants::core_tag) continue;
						if(!(Mixedpart[j].pt>constants::assoc_ptmin && Mixedpart[j].pt<constants::assoc_ptmax)) continue;

						double x_val=EVENT.part[i].eta - Mixedpart[j].eta;
						if(x_val<options.ih.x_edge_min || x_val>options.ih.x_edge_max) continue;
						int nx=this->get_cell_index(x_val);

						double y_val=this->getDeltaPhi(EVENT.part[i].phi, Mixedpart[j].phi);
						if(y_val<options.ih.y_edge_min || y_val>options.ih.y_edge_max) continue;
						int ny=this->get_cell_index_y(y_val);

						if(max_nx<nx) max_nx = nx;
						if(max_ny<ny) max_ny = ny;
						Hit1ev[nx][ny]+=1.0;
						ct->HistSub2D_x[nx][ny]+=x_val*EVENT.weight();
						ct->HistSub2D_y[nx][ny]+=y_val*EVENT.weight();
						if(ct->max_nx<nx) ct->max_nx=nx;
						if(ct->max_ny<ny) ct->max_ny=ny;
						NumPair++;

					}
				}
				//---------------
				
				for(int nx = 0; nx<max_nx+1; nx++){
					for(int ny = 0; ny<max_ny+1; ny++){
						ct->HistSub2D[nx][ny]+=Hit1ev[nx][ny]*EVENT.weight()/((double)NumTrig);
						ct->HistSub2DPartHit[nx][ny]+=Hit1ev[nx][ny]*EVENT.weight();
					}
				}

				//These are already calculated in signal calculation, so no need to do here.
				//====================================================================
				//ct->SumWeight+=EVENT.weight();


				for(int i = 0; i < constants::x_cell_capa; i++) {
					delete[] Hit1ev[i];
				}
				delete[] Hit1ev;

}



			void Fill::fill_Rt(shared_ptr<Container>& ct){

				Container::EventInfo& EVENT= ct->EVENTINFO;

				//Find max pt particle.
				//----------------------------
				double max_pt=-1.0;
				int itrig=-1;
				for(int i=0; i<(int)EVENT.part.size(); ++i){
					if(fabs(EVENT.part[i].pt)<constants::minpt_Rt) continue;
					int ID = (int)EVENT.part[i].ID;
					if(!(abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon)) continue;
					if(max_pt<EVENT.part[i].pt) {
						max_pt = EVENT.part[i].pt;
						itrig=i;
					}
				}
				if(itrig<0) {
					return;
				}

				//Get Nt.
				//-------------------------------
				int Nt=0;
				for(int j=0; j<(int)EVENT.part.size(); ++j){
					//Count Nt
					//-----------
					int ID = (int)EVENT.part[j].ID;
					if(!(abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon)) continue;
					if(fabs(EVENT.part[j].phi-EVENT.part[itrig].phi)<constants::maxPhi_RtTrans 
							&& fabs(EVENT.part[j].phi-EVENT.part[itrig].phi)>constants::minPhi_RtTrans){
						Nt++;
					}
				}
				//---------------
				//Store 
				//--------
				ct->Nt_eBye.push_back(Nt);
				ct->dNdeta_eBye.push_back(EVENT.Nch());
				ct->weight_eBye.push_back(EVENT.weight());
				ct->SumWeight+=EVENT.weight();
				return;
			}

			void Fill::fill_RtYield(shared_ptr<Container>& ct){

				Container::EventInfo& EVENT= ct->EVENTINFO;

				//Find max pt particle.
				//----------------------------
				double max_pt=-1.0;
				int itrig=-1;
				for(int i=0; i<(int)EVENT.part.size(); ++i){
					int ID = (int)EVENT.part[i].ID;
					if(!(abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon)) continue;
					if(fabs(EVENT.part[i].pt)<constants::minpt_Rt) continue;
					if(max_pt<EVENT.part[i].pt) {
						max_pt = EVENT.part[i].pt;
						itrig=i;
					}
				}
				if(itrig<0) {
					return;
				}

				//Get Nt.
				//-------------------------------
				int Nt=0;
				int Ntmin=0;
				int Ntmax=0;
				int Ntcore=0;
				int Ntcorona=0;
				int Nncore=0;
				int Nncorona=0;
				Container::yield TransYield;
				Container::yield TowardYield;
				for(int j=0; j<(int)EVENT.part.size(); ++j){
					//Count Nt (multiplicity in transverse region)
					//----------------------------------------------
					int ID = (int)EVENT.part[j].ID;
					if(fabs(EVENT.part[j].phi-EVENT.part[itrig].phi)<constants::maxPhi_RtTrans && fabs(EVENT.part[j].phi-EVENT.part[itrig].phi)>constants::minPhi_RtTrans && (abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon)){
						Nt++;
						if(EVENT.part[j].phi-EVENT.part[itrig].phi>0.0) Ntmin++;
						else Ntmax++;
					}
					if(Ntmin>Ntmax){
						int Nttemp=Ntmin;
						Ntmin=Ntmax;
						Ntmax=Nttemp;
					}

					//Count Nth(yield in transverse region)
					//----------------------------------------------
					if(fabs(EVENT.part[j].phi-EVENT.part[itrig].phi)<constants::maxPhi_RtTrans && fabs(EVENT.part[j].phi-EVENT.part[itrig].phi)>constants::minPhi_RtTrans){
						int ID = (int)EVENT.part[j].ID;
						string TAG = EVENT.part[j].TAG;
						if(abs(ID)==constants::id_proton) TransYield.add_ppbar(1.0);
							if(abs(ID)==constants::id_ch_pion)TransYield.add_chpi(1.0); 
							if(abs(ID)==constants::id_ch_kaon)TransYield.add_chkaon(1.0); 
							if(abs(ID)==constants::id_phi)TransYield.add_phi(1.0); 
							if(abs(ID)==constants::id_lambda)TransYield.add_lambda(1.0); 
							if(abs(ID)==constants::id_cascade)TransYield.add_cascade(1.0); 
							if(abs(ID)==constants::id_omega)TransYield.add_omega(1.0); 

							if(TAG == constants::core_tag)Ntcore++;
							else if(TAG == constants::corona_tag)Ntcorona++;

					}
					//Count Nth(yield in towards region)
					//----------------------------------------------
					else if(fabs(EVENT.part[j].phi-EVENT.part[itrig].phi)<constants::minPhi_RtTrans){
						int ID = (int)EVENT.part[j].ID;
						string TAG = EVENT.part[j].TAG;
						if(abs(ID)==constants::id_proton) TowardYield.add_ppbar(1.0);
							if(abs(ID)==constants::id_ch_pion)TowardYield.add_chpi(1.0); 
							if(abs(ID)==constants::id_ch_kaon)TowardYield.add_chkaon(1.0); 
							if(abs(ID)==constants::id_phi)TowardYield.add_phi(1.0); 
							if(abs(ID)==constants::id_lambda)TowardYield.add_lambda(1.0); 
							if(abs(ID)==constants::id_cascade)TowardYield.add_cascade(1.0); 
							if(abs(ID)==constants::id_omega)TowardYield.add_omega(1.0); 

							if(TAG == constants::core_tag)Nncore++;
							else if(TAG == constants::corona_tag)Nncorona++;
					}else{

					}
				}
				//---------------
				//Store 
				//--------
				ct->TransYield_eBye.push_back(TransYield);
				ct->TowardYield_eBye.push_back(TowardYield);
				ct->Nt_eBye.push_back(Nt);
				ct->dNdeta_eBye.push_back(EVENT.Nch());
				ct->CoreT_eBye.push_back(Ntcorona);
				ct->CoreN_eBye.push_back(Nncorona);
				ct->Ntmin_eBye.push_back(Ntmin);
				ct->Ntmax_eBye.push_back(Ntmax);
				ct->weight_eBye.push_back(EVENT.weight());
				ct->SumWeight+=EVENT.weight();
				ct->TagEventNum.push_back(ct->CountEv-1);
				return;
			}

			void Fill::fill_twopc1D_taggedInteg(shared_ptr<Container>& ct){

				Container::EventInfo& EVENT= ct->EVENTINFO;

				//Count particle by particle.
				//----------------------------
				int NumPair=0;
				int NumTrig=0;
				double x_val=EVENT.Nch();
				if(x_val<options.ih.x_edge_min || x_val>options.ih.x_edge_max) return;
				int nx=this->get_cell_index(x_val);
				for(int i=0; i<(int)EVENT.part.size(); ++i){
					if(EVENT.part[i].pt<constants::trig_ptmin) continue;
					NumTrig++;

					//Seeing associates.
					//-------------------------------
					for(int j=0; j<(int)EVENT.part.size(); ++j){
						if(i==j) continue;
						string TAG = EVENT.part[j].TAG;
						if(options.get_flag_only_core_associates() && TAG==constants::corona_tag) continue;
						if(options.get_flag_only_corona_associates() && TAG==constants::core_tag) continue;

						double DeltaPhi=this->getDeltaPhi_twopc1D(EVENT.part[i].phi, EVENT.part[j].phi);
						double DeltaEta=fabs(EVENT.part[i].eta - EVENT.part[j].eta);
						if(options.get_flag__2PCfull() && DeltaEta> constants::DeltaEtaFULL ) continue;
						else if(options.get_flag__2PCnearside() && (DeltaEta>constants::DeltaEtaNS || fabs(DeltaPhi)>constants::DeltaPhiNS)) continue;
						else if(options.get_flag__2PCout() && (DeltaEta<constants::DeltaEtaNS || DeltaEta>constants::DeltaEtaFULL || DeltaPhi<constants::DeltaPhiNS || DeltaPhi>constants::DeltaPhiOUT)) continue;


						ct->Hist[nx]+=1.0*EVENT.weight();
						NumPair++;

					}
				}
				//---------------
				ct->HistHit[nx]+=((double)NumTrig)*EVENT.weight();
				ct->Hist_x[nx]+=x_val*EVENT.weight();
				ct->Hist_weight[nx]+=EVENT.weight();
				ct->SumWeight+=EVENT.weight();
				if(ct->max_nx<nx) ct->max_nx=nx;
				return;
			}

			void Fill::fill_twopc1D_tagged(shared_ptr<Container>& ct){

				Container::EventInfo& EVENT= ct->EVENTINFO;

				//Count particle by particle.
				//----------------------------
				int NumPair=0;
				int NumTrig=0;
				for(int i=0; i<(int)EVENT.part.size(); ++i){
					if(EVENT.part[i].pt<constants::trig_ptmin) continue;
					NumTrig++;

					//Seeing associates.
					//-------------------------------
					for(int j=0; j<(int)EVENT.part.size(); ++j){
						if(i==j) continue;
						string TAG = EVENT.part[j].TAG;
						if(options.get_flag_only_core_associates() && TAG==constants::corona_tag) continue;
						if(options.get_flag_only_corona_associates() && TAG==constants::core_tag) continue;

						double x_val=this->getDeltaPhi_twopc1D(EVENT.part[i].phi, EVENT.part[j].phi);
						double DeltaEta=fabs(EVENT.part[i].eta - EVENT.part[j].eta);
						if(options.get_flag__2PCfull() && DeltaEta> constants::DeltaEtaFULL) continue;
						else if(options.get_flag__2PCnearside() && DeltaEta>constants::DeltaEtaNS) continue;
						else if(options.get_flag__2PCout() && (DeltaEta<constants::DeltaEtaNS || DeltaEta>constants::DeltaEtaFULL)) continue;

						if(x_val<options.ih.x_edge_min || x_val>options.ih.x_edge_max) continue;
						int nx=this->get_cell_index(x_val);

						ct->Hist[nx]+=1.0*EVENT.weight();
						ct->Hist_x[nx]+=x_val*EVENT.weight();
						ct->Hist_weight[nx]+=EVENT.weight();
						ct->HistHit[nx]++;
						if(ct->max_nx<nx) ct->max_nx=nx;
						NumPair++;

					}
				}
				//---------------
				ct->SumWeight+=EVENT.weight();
				return;
			}


			void Fill::fill_twopc1D_tagged_1particle(shared_ptr<Container>& ct){

				Container::EventInfo& EVENT= ct->EVENTINFO;

				//Count particle by particle.
				//----------------------------
				int NumPair=0;
				int NumTrig=0;
				int itrig=-1;
				double max_pt=-1.0;
				for(int i=0; i<(int)EVENT.part.size(); ++i){
					if(abs((int)EVENT.part[i].ID)==constants::id_K0S) continue;
					if(max_pt<EVENT.part[i].pt) {
						max_pt = EVENT.part[i].pt;
						itrig=i;
					}
				}
				if(itrig<0) {
					NumTrig++;
				}

				//Seeing K0S.
				//-------------------------------
				for(int j=0; j<(int)EVENT.part.size(); ++j){
					if(itrig<0) {
						break;
					}
					if(abs((int)EVENT.part[j].ID)!=constants::id_K0S) continue;

					double x_val=this->getDeltaPhi_twopc1D(EVENT.part[itrig].phi, EVENT.part[j].phi);
					double DeltaEta=fabs(EVENT.part[itrig].eta - EVENT.part[j].eta);
					if(options.get_flag__2PCfull() && DeltaEta> constants::DeltaEtaFULL) continue;
					else if(options.get_flag__2PCnearside() && DeltaEta>constants::DeltaEtaNS) continue;
					else if(options.get_flag__2PCout() && (DeltaEta<constants::DeltaEtaNS || DeltaEta>constants::DeltaEtaFULL)) continue;

					if(x_val<options.ih.x_edge_min || x_val>options.ih.x_edge_max) continue;
					int nx=this->get_cell_index(x_val);

					ct->Hist[nx]+=1.0*EVENT.weight();
					ct->Hist_x[nx]+=x_val*EVENT.weight();
					ct->Hist_weight[nx]+=EVENT.weight();
					ct->HistHit[nx]++;
					if(ct->max_nx<nx) ct->max_nx=nx;
					NumPair++;

				}
				//---------------
				ct->SumWeight+=EVENT.weight();
				return;
			}

			void Fill::fill_twopc1D(shared_ptr<Container>& ct){

				Container::EventInfo& EVENT= ct->EVENTINFO;
				double N=EVENT.Nch();
				if(N<constants::twopc1D_Nmin || N>constants::twopc1D_Nmax) return;

				//Count particle by particle.
				//----------------------------
				int NumPair=0;
				int NumTrig=0;
				for(int i=0; i<(int)EVENT.part.size(); ++i){

					for(int j=0; j<(int)EVENT.part.size(); ++j){
						if(i==j) continue;

						double x_val=this->getDeltaPhi_twopc1D(EVENT.part[i].phi, EVENT.part[j].phi);
						double DeltaEta=fabs(EVENT.part[i].eta - EVENT.part[j].eta);
						if(options.get_flag__2PCfull() && DeltaEta>constants::DeltaEtaFULL) continue;
						else if(options.get_flag__2PCnearside() && DeltaEta>constants::DeltaEtaNS) continue;
						else if(options.get_flag__2PCout() && (DeltaEta<constants::DeltaEtaNS || DeltaEta>constants::DeltaEtaFULL)) continue;

						if(x_val<options.ih.x_edge_min || x_val>options.ih.x_edge_max) continue;
						int nx=this->get_cell_index(x_val);

						ct->Hist[nx]+=1.0*EVENT.weight();
						ct->Hist_x[nx]+=x_val*EVENT.weight();
						ct->Hist_weight[nx]+=EVENT.weight();
						ct->HistHit[nx]++;
						if(ct->max_nx<nx) ct->max_nx=nx;
						NumPair++;

					}
					NumTrig++;
				}
				//---------------

				ct->SumWeight+=EVENT.weight();
				return;
			}


			double Fill::getDeltaPhi(const double phi1, const double phi2){
				double deltaPhi = phi1 - phi2;
				if(deltaPhi<0.0) deltaPhi += 2.0*M_PI;

				if(deltaPhi>=options.ih.y_edge_min && deltaPhi<=options.ih.y_edge_max){
					return deltaPhi;
				}else if(deltaPhi<options.ih.y_edge_min){
					return deltaPhi + 2.0*M_PI;
				}else if(deltaPhi>options.ih.y_edge_max){
					return deltaPhi-2.0*M_PI;
				}else{
					cout << ":( Error. Something wrong in double getDeltaPhi." << deltaPhi << endl;
					exit(1);
				}

			}

			double Fill::getDeltaPhi_twopc1D(const double phi1, const double phi2){
				//1. Fix range from -2pi<phi<2pi to 0<phi<2pi
				//----------------------------
				double deltaPhi = phi1-phi2;
				if(deltaPhi<0.0) deltaPhi += 2.0*M_PI;

				if(deltaPhi>=options.ih.x_edge_min && deltaPhi<=options.ih.x_edge_max){
					return deltaPhi;
				}else if(deltaPhi<options.ih.x_edge_min){
					return deltaPhi + 2.0*M_PI;
				}else if(deltaPhi>options.ih.x_edge_max){
					return deltaPhi-2.0*M_PI;
				}else{
					cout << ":( Error. Something wrong in double getDeltaPhi." << deltaPhi << endl;
					exit(1);
				}

			}


			void Fill::fill(shared_ptr<Container>& ct, const double dNdeta, const int bin){


				Container::EventInfo& EVENT= ct->EVENTINFO;

				//Determine xbin
				//---------------
				double x_val=constants::dummy;
				int nx=-1;
				if(constants::MODE.find("meanpt")!=string::npos || constants::MODE.find("meanmt")!=string::npos || constants::MODE.find("MeanptPID")!=string::npos){
					x_val=EVENT.Nch();
					nx=(!options.get_flag_HI())? this->get_cell_index(x_val):this->get_cell_index_logplot(x_val);
					if(x_val<options.ih.x_edge_min || x_val>options.ih.x_edge_max) return;
				}

				//For mean pt PID vs. dNdeta
				//--------------------------------------------------
				if(options.get_flag_vs_Multi()){
					nx=bin;
					x_val=dNdeta;
				}



				//Count particle by particle.
				//----------------------------
				for(int j=0; j<(int)EVENT.part.size(); ++j){
					double y_val=constants::dummy;
					if(constants::MODE.find("meanpt")!=string::npos || constants::MODE.find("MeanptPID")!=string::npos){
						y_val = EVENT.part[j].pt;
					}else if(constants::MODE.find("meanmt")!=string::npos){
						y_val = EVENT.part[j].mt - EVENT.part[j].m;
					}else if(constants::MODE.find("mtscaling")!=string::npos){
						x_val=EVENT.part[j].m;
						if(x_val<options.ih.x_edge_min || x_val>options.ih.x_edge_max) continue;
						nx=this->get_cell_index(x_val);

			 			if(EVENT.part[j].ID==constants::id_phi){
							uf->checkMassOnShell(EVENT.part[j].m, EVENT.part[j].e, EVENT.part[j].px, EVENT.part[j].py, EVENT.part[j].pz);
						}

						if(!this->fix_ax(EVENT.part[j].ID, nx, x_val)) continue;
						y_val = EVENT.part[j].mt - EVENT.part[j].m;
					}else{
						//Default filling.
						//You can put whatever you want in x and y.
						//==========================================
						x_val=EVENT.part[j].rap;
						y_val = EVENT.part[j].e/(options.ih.d_x);

						//DNDPT
						//=====
						if(options.get_obs_type().find("dndpt")!=string::npos){
							//if(fabs(EVENT.part[j].rap)>0.5) continue;
							x_val = EVENT.part[j].pt;
							y_val = 1.0/(options.ih.d_x);
						}

						//Find bin
						//========
						nx=this->get_cell_index(x_val);
					}

					ct->Hist[nx]+=y_val*EVENT.weight();
					ct->Hist_x[nx]+=x_val*EVENT.weight();
					ct->HistHit[nx]++;

					//Temporal archive of Njet in one event for each bin
					//=================================================
					ct->Hist_1ev[nx]+=y_val;

					ct->Hist_weight[nx]+=EVENT.weight();
					if(ct->max_nx<nx) ct->max_nx=nx;


				}//particle loop

				ct->SumWeight+=EVENT.weight();

				//loop over bins. This should be done after particle loop
				//======================================================
				for(int i=0; i<ct->max_nx+1; ++i){

					//I want to get sum_ev Nminijet*Nminijet 
					//in [nx]
					//===================================
					ct->HistHist[i]+=pow(ct->Hist_1ev[i], 2);
					ct->Hist_1ev[i]=0.0;
				}

			}



int Fill::get_cell_index_cstm(const double val){

	if(val<xMin_cstm[0] && fabs(xMin_cstm[0])<constants::SMALL){
		return 0;
	}
	if(val>=xMax_cstm[(int)xMax_cstm.size()-1]){
		return -1;
	}
	for(int i=0; i<(int)xMin_cstm.size(); i++){
		if(val>=xMin_cstm[i] && val<xMax_cstm[i]) {
			if(fabs(xMin_cstm[0])<constants::SMALL) return i+1;//#0 is 0.0-xMin[0]
			else return i;
		}
	} 

	cout << __FILE__ << " line " << __LINE__ << " Something is wrong " << endl;
	exit(1);

}


			bool Fill::fix_ax(const int id, int &nx, double m){

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


			int Fill::get_cell_index(const double x_val){
				double half_dx = options.ih.d_x / 2.0;
				double binZeroCenter = 0;
				int ncell = std::round((x_val - binZeroCenter + std::fabs(options.ih.x_edge_min) - half_dx) / options.ih.d_x);
				if(options.get_xaxis_type()==3){
					return this->get_cell_index_cstm(x_val);
				}
				if(ncell<0){
					std::cout << "ERROR:( something is wrong.  ncell:" << ncell << " out of constants::x_cell_capa " << constants::x_cell_capa << std::endl;
					std::cout << "                                       x_val:" << x_val << " options.ih.x_edge_min " << options.ih.x_edge_min << std::endl;
					exit(EXIT_FAILURE);
				}
				if(ncell>constants::x_cell_capa){
					std::cout << "WARNING:O ncell is over the capacity.  ncell:" << ncell << " out of constants::x_cell_capa " << constants::x_cell_capa << std::endl;
					std::cout << "                                       x_val:" << x_val << " options.ih.x_edge_min " << options.ih.x_edge_min << std::endl;
                                        ncell=constants::x_cell_capa-1;//Put everything into the last cell.
				}
				return  ncell;
			}
			int Fill::get_cell_index_y(const double y_val){
				double half_dy = options.ih.d_y / 2.0;
				double binZeroCenter = 0;
				int ncell = std::round((y_val - binZeroCenter + std::fabs(options.ih.y_edge_min) - half_dy) / options.ih.d_y);
				if(ncell<0 || ncell>constants::y_cell_capa){
					std::cout << "ERROR:( ncell is beyond the capacity.  ncell:" << ncell << " out of constants::y_cell_capa " << constants::y_cell_capa << std::endl;
					std::cout << "                                       y_val:" << y_val << " options.ih.y_edge_min " << options.ih.y_edge_min << std::endl;
					exit(EXIT_FAILURE);
				}
				return  ncell;
			}



			int Fill::get_cell_index_logplot(const double x_val_){

				double x_val=x_val_;
				int ncell=0;


				if(x_val<constants::switchBin_x){
					int ncell=(int) floor((x_val-options.ih.x_edge_min)/constants::binSize_small);

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



			void Fill::fill_vn4multi(shared_ptr<Container>& ct, const double Nch){


				Container::EventInfo& EVENT= ct->EVENTINFO;

				if((int)EVENT.part.size()<4) return;

				//Determine xbin
				//---------------
				double x_val=EVENT.Nch();
				if(x_val<options.ih.x_edge_min || x_val>options.ih.x_edge_max) return;
				int nx=(!options.get_flag_HI())? this->get_cell_index(x_val):this->get_cell_index_logplot(x_val);
				//If xval is calculated from different particle list.
				//======================================================
				if(options.get_flag_vs_Multi()){
					x_val=Nch;
					nx=(!options.get_flag_HI())? this->get_cell_index(x_val):this->get_cell_index_logplot(x_val);
				}


				//Count particle by particle.
				//----------------------------
				std::complex<double> Qvec_tot=constants::initialval_comp;
				std::complex<double> Qvec2_tot=constants::initialval_comp;
				std::complex<double> n_coeff (options.ih.N_coeff, 0.0);
				std::complex<double> n2_coeff (options.ih.N_coeff*2.0, 0.0);
				for(int j=0; j<(int)EVENT.part.size(); ++j){
					std::complex<double> phi_ (EVENT.part[j].phi,0.0);
					std::complex<double> Qvec=exp(constants::i_img*n_coeff*phi_);
					std::complex<double> Qvec2=exp(constants::i_img*n2_coeff*phi_);
					Qvec_tot += Qvec;
					Qvec2_tot += Qvec2;
				}

				double totN = (double)EVENT.part.size();
				double corr4_num = real(Qvec_tot*Qvec_tot*conj(Qvec_tot)*conj(Qvec_tot) 
						+  Qvec2_tot*conj(Qvec2_tot)
						-2.0*real(Qvec2_tot*conj(Qvec_tot)*conj(Qvec_tot))
						- Qvec_tot*conj(Qvec_tot)*(4.0*totN-8.0)
						 - 6.0*totN +2.0*pow(totN,2));
				double corr4_denom = pow(totN,4)-6.0*pow(totN,3)+11.0*pow(totN,2)-6.0*totN;


				//Obtain 2particle correlation
				//-------------------------------
				double squared_Qvec = real(Qvec_tot * conj(Qvec_tot));
				double corr2_num = squared_Qvec-totN;
				double corr2_denom = totN*(totN-1.0);


				ct->Hist[nx]+=corr4_num*EVENT.weight();
				ct->Hist2[nx]+=corr4_denom*EVENT.weight();
				ct->Hist_sub[nx]+=corr2_num*EVENT.weight();
				ct->Hist2_sub[nx]+=corr2_denom*EVENT.weight();
				ct->Hist_x[nx]+=x_val*EVENT.weight();
				ct->HistHit[nx]++;
				ct->HistHist[nx]+=corr4_num*corr4_num*EVENT.weight();
				ct->HistHist2[nx]+=corr4_denom*corr4_denom*EVENT.weight();
				ct->HistHist_sub[nx]+=corr2_num*corr2_num*EVENT.weight();
				ct->HistHist2_sub[nx]+=corr2_denom*corr2_denom*EVENT.weight();
				ct->Hist_weight[nx]+=EVENT.weight();
				if(ct->max_nx<nx) ct->max_nx=nx;

				ct->SumWeight+=EVENT.weight();

			}



			void Fill::fill_vnmulti(shared_ptr<Container>& ct){


				Container::EventInfo& EVENT= ct->EVENTINFO;

				if((int)EVENT.part.size()<2) return;

				//Determine xbin
				//---------------
				double x_val=EVENT.Nch();
				if(x_val<options.ih.x_edge_min || x_val>options.ih.x_edge_max) return;
				int nx=(!options.get_flag_HI())? this->get_cell_index(x_val):this->get_cell_index_logplot(x_val);

				//Count particle by particle.
				//----------------------------
				std::complex<double> Qvec_tot=constants::initialval_comp;
				std::complex<double> n_coeff (options.ih.N_coeff, 0.0);
				for(int j=0; j<(int)EVENT.part.size(); ++j){
					std::complex<double> phi_ (EVENT.part[j].phi,0.0);
					std::complex<double> Qvec=exp(constants::i_img*n_coeff*phi_);
					Qvec_tot += Qvec;
				}
				double squared_Qvec = real(Qvec_tot * conj(Qvec_tot));
				double squared_Qvec_img = imag(Qvec_tot * conj(Qvec_tot));
				double totN = (double)EVENT.part.size();
				double corr_num = squared_Qvec-totN;
				double corr_denom = totN*(totN-1.0);

				ct->Hist[nx]+=corr_num*EVENT.weight();
				ct->Hist2[nx]+=corr_denom*EVENT.weight();
				ct->Hist_x[nx]+=x_val*EVENT.weight();
				ct->HistHit[nx]++;
				ct->HistHist[nx]+=corr_num*corr_num*EVENT.weight();
				ct->HistHist2[nx]+=corr_denom*corr_denom*EVENT.weight();
				ct->Hist_weight[nx]+=EVENT.weight();
				ct->Hist_img_Qvec[nx]+=squared_Qvec_img*EVENT.weight();
				if(ct->max_nx<nx) ct->max_nx=nx;

				ct->SumWeight+=EVENT.weight();

			}





			void Fill::fill_vneta(shared_ptr<Container>& ct){


				Container::EventInfo& EVENT= ct->EVENTINFO;

				//Count particle by particle.
				//----------------------------
				std::complex<double> Qvec_tot[constants::x_cell_capa]={};
				double hit[constants::x_cell_capa]={};
				for(int i=0; i<constants::x_cell_capa; i++){
					Qvec_tot[i]=constants::initialval_comp;
					hit[i]=0.0;
				}
				std::complex<double> n_coeff (options.ih.N_coeff, 0.0);
				for(int j=0; j<(int)EVENT.part.size(); ++j){

					//Determine xbin
					//---------------
					double x_val=EVENT.part[j].eta;
					if(x_val<options.ih.x_edge_min || x_val>options.ih.x_edge_max) continue;
					int nx=this->get_cell_index(x_val);

					std::complex<double> phi_ (EVENT.part[j].phi,0.0);
					std::complex<double> Qvec=exp(constants::i_img*n_coeff*phi_);
					Qvec_tot[nx] += Qvec;
					//      cout << "Qvec " << Qvec << endl;
					hit[nx]++;
					ct->Hist_x[nx]+=x_val*EVENT.weight();
					ct->HistPartHit[nx]+=EVENT.weight();

					if(ct->max_nx<nx) ct->max_nx=nx;

				}

				for(int nx=0; nx<ct->max_nx+constants::margin; nx++){
					double squared_Qvec = real(Qvec_tot[nx] * conj(Qvec_tot[nx]));
					double squared_Qvec_img = imag(Qvec_tot[nx] * conj(Qvec_tot[nx]));
					double totN = hit[nx];
					if(totN<2) continue;
					double corr_num = squared_Qvec-totN;
					double corr_denom = totN*(totN-1.0);

					ct->Hist[nx]+=corr_num*EVENT.weight();
					ct->Hist2[nx]+=corr_denom*EVENT.weight();
					ct->HistHit[nx]++;
					ct->HistHist[nx]+=corr_num*corr_num*EVENT.weight();
					ct->HistHist2[nx]+=corr_denom*corr_denom*EVENT.weight();
					ct->Hist_weight[nx]+=EVENT.weight();
					ct->Hist_img_Qvec[nx]+=squared_Qvec_img*EVENT.weight();

					ct->SumWeight+=EVENT.weight();
				}


			}


			void Fill::fill_vnpt(shared_ptr<Container>& ct){


				Container::EventInfo& EVENT= ct->EVENTINFO;



				//Count particle by particle.
				//----------------------------
				std::complex<double> Qvec_tot[constants::x_cell_capa]={};
				double hit[constants::x_cell_capa]={};
				for(int i=0; i<constants::x_cell_capa; i++){
					Qvec_tot[i]=constants::initialval_comp;
					hit[i]=0.0;
				}
				std::complex<double> n_coeff (options.ih.N_coeff, 0.0);
				for(int j=0; j<(int)EVENT.part.size(); ++j){

					//Determine xbin
					//---------------
					double x_val=EVENT.part[j].pt;
					if(x_val<options.ih.x_edge_min || x_val>options.ih.x_edge_max) continue;
					int nx= this->get_cell_index(x_val);

					std::complex<double> phi_ (EVENT.part[j].phi,0.0);
					std::complex<double> Qvec=exp(constants::i_img*n_coeff*phi_);
					Qvec_tot[nx] += Qvec;
					//      cout << "Qvec " << Qvec << endl;
					hit[nx]++;
					ct->Hist_x[nx]+=x_val*EVENT.weight();
					ct->HistPartHit[nx]+=EVENT.weight();

					if(ct->max_nx<nx) ct->max_nx=nx;

				}

				for(int nx=0; nx<ct->max_nx+constants::margin; nx++){
					double squared_Qvec = real(Qvec_tot[nx] * conj(Qvec_tot[nx]));
					double squared_Qvec_img = imag(Qvec_tot[nx] * conj(Qvec_tot[nx]));
					double totN = hit[nx];
					if(totN<2) continue;
					double corr_num = squared_Qvec-totN;
					double corr_denom = totN*(totN-1.0);

					ct->Hist[nx]+=corr_num*EVENT.weight();
					ct->Hist2[nx]+=corr_denom*EVENT.weight();
					ct->HistHit[nx]++;
					ct->HistHist[nx]+=corr_num*corr_num*EVENT.weight();
					ct->HistHist2[nx]+=corr_denom*corr_denom*EVENT.weight();
					ct->Hist_weight[nx]+=EVENT.weight();
					ct->Hist_img_Qvec[nx]+=squared_Qvec_img*EVENT.weight();

					ct->SumWeight+=EVENT.weight();
				}


			}




			//-------2sub


			void Fill::fill_vnpt_2sub(shared_ptr<Container>& ct){


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
				std::complex<double> n_coeff (options.ih.N_coeff, 0.0);
				for(int j=0; j<(int)EVENT.part.size(); ++j){

					//Determine xbin
					//---------------
					double x_val=EVENT.part[j].pt;
					if(x_val<options.ih.x_edge_min || x_val>options.ih.x_edge_max) continue;
					int nx= this->get_cell_index(x_val);

					std::complex<double> phi_ (EVENT.part[j].phi,0.0);
					std::complex<double> Qvec=exp(constants::i_img*n_coeff*phi_);
					if(EVENT.part[j].eta<-(constants::etaGap/2.0)){
						Qvec_tot_A[nx] += Qvec;
						hit_A[nx]++;
					}else if(EVENT.part[j].eta>(constants::etaGap/2.0)){
						Qvec_tot_B[nx] += Qvec;
						hit_B[nx]++;
					}

					if(ct->max_nx<nx) ct->max_nx=nx;
					ct->Hist_x[nx]+=x_val*EVENT.weight();
					ct->HistPartHit[nx]+=EVENT.weight();

				}

				for(int nx=0; nx<ct->max_nx+constants::margin; nx++){
					double squared_Qvec = real(Qvec_tot_A[nx] * conj(Qvec_tot_B[nx]));
					double squared_Qvec_img = imag(Qvec_tot_A[nx] * conj(Qvec_tot_B[nx]));
					if(hit_A[nx]==0.0 || hit_B[nx]==0.0) continue;
					double corr_num = squared_Qvec;
					double corr_denom = hit_A[nx]*hit_B[nx];

					ct->Hist[nx]+=corr_num*EVENT.weight();
					ct->Hist2[nx]+=corr_denom*EVENT.weight();
					ct->HistHit[nx]++;
					ct->HistHist[nx]+=corr_num*corr_num*EVENT.weight();
					ct->HistHist2[nx]+=corr_denom*corr_denom*EVENT.weight();
					ct->Hist_weight[nx]+=EVENT.weight();
					ct->Hist_img_Qvec[nx]+=squared_Qvec_img*EVENT.weight();

					ct->SumWeight+=EVENT.weight();
				}


			}





			void Fill::fill_vnmulti_2sub(shared_ptr<Container>& ct){


				Container::EventInfo& EVENT= ct->EVENTINFO;


				//Determine xbin
				//---------------
				double x_val=EVENT.Nch();
				if(x_val<options.ih.x_edge_min || x_val>options.ih.x_edge_max) return;
				int nx=(!options.get_flag_HI())? this->get_cell_index(x_val):this->get_cell_index_logplot(x_val);

				//Count particle by particle.
				//----------------------------
				std::complex<double> Qvec_tot_A=constants::initialval_comp;
				std::complex<double> Qvec_tot_B=constants::initialval_comp;
				std::complex<double> n_coeff (options.ih.N_coeff, 0.0);
				double hit_A=0.0, hit_B=0.0;
				for(int j=0; j<(int)EVENT.part.size(); ++j){
					std::complex<double> phi_ (EVENT.part[j].phi,0.0);
					std::complex<double> Qvec=exp(constants::i_img*n_coeff*phi_);
					if(EVENT.part[j].eta<-(constants::etaGap/2.0)){
						Qvec_tot_A += Qvec;
						hit_A++;
					}else if(EVENT.part[j].eta>(constants::etaGap/2.0)){
						Qvec_tot_B += Qvec;
						hit_B++;
					}
				}
				double squared_Qvec = real(Qvec_tot_A * conj(Qvec_tot_B));
				double squared_Qvec_img = imag(Qvec_tot_A * conj(Qvec_tot_B));

				if(hit_A==0.0 || hit_B==0.0) return;

				double corr_num = squared_Qvec;
				double corr_denom = hit_A*hit_B;

				ct->Hist[nx]+=corr_num*EVENT.weight();
				ct->Hist2[nx]+=corr_denom*EVENT.weight();
				ct->Hist_x[nx]+=x_val*EVENT.weight();
				ct->HistHit[nx]++;
				ct->HistHist[nx]+=corr_num*corr_num*EVENT.weight();
				ct->HistHist2[nx]+=corr_denom*corr_denom*EVENT.weight();
				ct->Hist_weight[nx]+=EVENT.weight();
				ct->Hist_img_Qvec[nx]+=squared_Qvec_img*EVENT.weight();
				if(ct->max_nx<nx) ct->max_nx=nx;

				ct->SumWeight+=EVENT.weight();

			}


			void Fill::fill_vn4multi_2sub(shared_ptr<Container>& ct){


				Container::EventInfo& EVENT= ct->EVENTINFO;


				//Determine xbin
				//---------------
				double x_val=EVENT.Nch();
				if(x_val<options.ih.x_edge_min || x_val>options.ih.x_edge_max) return;
				int nx=(!options.get_flag_HI())? this->get_cell_index(x_val):this->get_cell_index_logplot(x_val);

				//Count particle by particle.
				//----------------------------
				std::complex<double> Qvec_tot_A=constants::initialval_comp;
				std::complex<double> Qvec_tot_B=constants::initialval_comp;
				std::complex<double> Qvec2_tot_A=constants::initialval_comp;
				std::complex<double> Qvec2_tot_B=constants::initialval_comp;
				std::complex<double> n_coeff (options.ih.N_coeff, 0.0);
				std::complex<double> n2_coeff (options.ih.N_coeff*2.0, 0.0);
				double hit_A=0.0, hit_B=0.0;
				for(int j=0; j<(int)EVENT.part.size(); ++j){
					std::complex<double> phi_ (EVENT.part[j].phi,0.0);
					std::complex<double> Qvec=exp(constants::i_img*n_coeff*phi_);
					std::complex<double> Qvec2=exp(constants::i_img*n2_coeff*phi_);
					if(EVENT.part[j].eta<-(constants::etaGap/2.0)){
						Qvec_tot_A += Qvec;
						Qvec2_tot_A += Qvec2;
						hit_A++;
					}else if(EVENT.part[j].eta>(constants::etaGap/2.0)){
						Qvec_tot_B += Qvec;
						Qvec2_tot_B += Qvec2;
						hit_B++;
					}
				}
				double squared_Qvec = real((pow(Qvec_tot_A,2) -  Qvec2_tot_A)* conj(pow(Qvec_tot_B,2) - Qvec2_tot_B));
				double squared_Qvec_img = imag((pow(Qvec_tot_A,2) -  Qvec2_tot_A)* conj(pow(Qvec_tot_B,2) - Qvec2_tot_B));

				if(hit_A<=1.0 || hit_B<=1.0) return;

				double corr4_num = squared_Qvec;
				double corr4_denom = hit_A*(hit_A-1)*hit_B*(hit_B-1);

				//Obtain 2particle correlation
				//------------------------------
				double squared_Qvec_2part = real(Qvec_tot_A * conj(Qvec_tot_B));
				double corr2_num = squared_Qvec_2part;
				double corr2_denom = hit_A*hit_B;


				ct->Hist[nx]+=corr4_num*EVENT.weight();
				ct->Hist_sub[nx]+=corr2_num*EVENT.weight();
				ct->Hist2[nx]+=corr4_denom*EVENT.weight();
				ct->Hist2_sub[nx]+=corr2_denom*EVENT.weight();
				ct->Hist_x[nx]+=x_val*EVENT.weight();
				ct->HistHit[nx]++;
				ct->HistHist[nx]+=corr4_num*corr4_num*EVENT.weight();
				ct->HistHist_sub[nx]+=corr2_num*corr2_num*EVENT.weight();
				ct->HistHist2[nx]+=corr4_denom*corr4_denom*EVENT.weight();
				ct->HistHist2_sub[nx]+=corr2_denom*corr2_denom*EVENT.weight();
				ct->Hist_weight[nx]+=EVENT.weight();
				ct->Hist_img_Qvec[nx]+=squared_Qvec_img*EVENT.weight();
				if(ct->max_nx<nx) ct->max_nx=nx;

				ct->SumWeight+=EVENT.weight();

			}



			//----------3sub

			void Fill::fill_vn4multi_3sub(shared_ptr<Container>& ct){


				Container::EventInfo& EVENT= ct->EVENTINFO;


				//Determine xbin
				//---------------
				double x_val=EVENT.Nch();
				if(x_val<options.ih.x_edge_min || x_val>options.ih.x_edge_max) return;
				int nx=(!options.get_flag_HI())? this->get_cell_index(x_val):this->get_cell_index_logplot(x_val);

				//Count particle by particle.
				//----------------------------
				std::complex<double> Qvec_tot_A=constants::initialval_comp;
				std::complex<double> Qvec_tot_B=constants::initialval_comp;
				std::complex<double> Qvec_tot_C=constants::initialval_comp;
				std::complex<double> Qvec2_tot_A=constants::initialval_comp;
				std::complex<double> Qvec2_tot_B=constants::initialval_comp;
				std::complex<double> Qvec2_tot_C=constants::initialval_comp;
				std::complex<double> n_coeff (options.ih.N_coeff, 0.0);
				std::complex<double> n2_coeff (options.ih.N_coeff*2.0, 0.0);
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

				double corr4_num = squared_Qvec;
				double corr4_denom = hit_A*(hit_A-1)*hit_B*hit_C;

				//Obtain 2particle correlation btw A and B.
				//----------------------------------------
				double squared_Qvec_2partAB = real(Qvec_tot_A * conj(Qvec_tot_B));
				double corr2_AB_num = squared_Qvec_2partAB;
				double corr2_AB_denom = hit_A*hit_B;


				//Obtain 2particle correlation btw A and C.
				//----------------------------------------
				double squared_Qvec_2partAC = real(Qvec_tot_A * conj(Qvec_tot_C));
				double corr2_AC_num = squared_Qvec_2partAC;
				double corr2_AC_denom = hit_A*hit_C;


				ct->Hist[nx]+=corr4_num*EVENT.weight();
				ct->Hist2[nx]+=corr4_denom*EVENT.weight();
				ct->Hist_sub[nx]+=corr2_AB_num*EVENT.weight();
				ct->Hist2_sub[nx]+=corr2_AB_denom*EVENT.weight();
				ct->Hist_subsub[nx]+=corr2_AC_num*EVENT.weight();
				ct->Hist2_subsub[nx]+=corr2_AC_denom*EVENT.weight();
				ct->Hist_x[nx]+=x_val*EVENT.weight();
				ct->HistHit[nx]++;
				ct->HistHist[nx]+=corr4_num*corr4_num*EVENT.weight();
				ct->HistHist2[nx]+=corr4_denom*corr4_denom*EVENT.weight();
				ct->HistHist_sub[nx]+=corr2_AB_num*corr2_AB_num*EVENT.weight();
				ct->HistHist2_sub[nx]+=corr2_AB_denom*corr2_AB_denom*EVENT.weight();
				ct->HistHist_subsub[nx]+=corr2_AC_num*corr2_AC_num*EVENT.weight();
				ct->HistHist2_subsub[nx]+=corr2_AC_denom*corr2_AC_denom*EVENT.weight();
				ct->Hist_weight[nx]+=EVENT.weight();
				ct->Hist_img_Qvec[nx]+=squared_Qvec_img*EVENT.weight();
				if(ct->max_nx<nx) ct->max_nx=nx;

				ct->SumWeight+=EVENT.weight();

			}





			void Fill::fill_TimeLapse(shared_ptr<Container>& ct){


				Container::EventInfo& EVENT= ct->EVENTINFO;

				//Step
				//----------------------------
				for(int j=0; j<(int)EVENT.step.size(); ++j){

					//Fill temp for each tau.
					//========================

					//Determine xbin
					//---------------
					int nx=EVENT.step[j].nstep;
					double x_val=EVENT.step[j].tau;
					double y_val=0.0;
					if(options.get_valTL().find("temp")!=string::npos){
						y_val=EVENT.step[j].temp;
					}else if(options.get_valTL().find("entropy")!=string::npos){
						y_val=EVENT.step[j].s;
					}else if(options.get_valTL().find("pressure")!=string::npos){
						y_val=EVENT.step[j].p;
					}else if(options.get_valTL().find("energy")!=string::npos){
						y_val=EVENT.step[j].e;
					}else{
						cout << ":(ERROR Set proper value with -valTL" << endl;
						exit(1);
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



vector<Container::ParticleInfo> Fill::select_N_RndomEv (const vector<EbyeInfo>& eBye_All){

        vector <Container::ParticleInfo> Mixedpart;
    
	std::uniform_int_distribution<> rndomEv(0,this->N_iCentEv-1);
	std::uniform_int_distribution<> rndomSamp_i(0,constants::Nsample-1);
	for(int i=0; i<constants::N_RndomEv; i++){
		int i_=rndomEv(rndom->generator);
		int iSamp_=rndomSamp_i(rndom->generatorSamp_i);
		int iEv = rndom->get_iEv_Cent()[i_];
                //Mixedpart.push_back(eBye_All[iEv].sample_part);
                Mixedpart.push_back(eBye_All[iEv].sample_partSet[iSamp_]);
	}

return Mixedpart;
}

