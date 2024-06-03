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
#include "Options.h"
#include "Util_func.h"
#include "Container.h"
#include "Message.h"
#include "LogSettings.h"
#include "Settings.h"
#include "CentralityCut.h"
#include "ReadIn.h"
#include "Fill.h"
#include "Stat.h"
#include "Write.h"
#include "BSTR.h"


class Analysis{

	private:

		int PrintCounter;
		std::shared_ptr<Util_func> uf;
		std::shared_ptr<Message> ms;  
		std::shared_ptr<ReadIn> read;
		std::shared_ptr<Fill> fill;
		std::shared_ptr<Stat> stat;
		std::shared_ptr<Write> write;
		std::shared_ptr<Rndom> rndom;
		std::shared_ptr<BSTR> bstr;
		std::vector<Container> ct_ALL;
		std::vector<Container::EventInfo> nEventInfo_for_BSTR;
		std::vector<EbyeInfo> eBye_CentCut_for_BSTR;

		Options& options;
		LogSettings& log;



		void shuffling_or_takeout_nevents(vector<Container::EventInfo>& nEventInfo, vector<EbyeInfo>& eBye){

			if(options.get_flag_shuffle() || options.get_flag_BSTR()){

				//Shuffle!
				//=======
				cout << ":D Shuffling the reading events!  nEventInfo.size():" << (int) nEventInfo.size() 
					<< ",  eBye.size():" << (int)eBye.size() 
					<< endl;
				std::random_device rnd_device;
				std::shuffle(std::begin(nEventInfo), std::end(nEventInfo), mt19937(rnd_device()));

				//I want to pick up the first n events from the shuffled nEventInfo vector.
				//Then I want to delete corresponding elements of eBye to the deleted ones.
				//==================================================
				int Ndiscard = (int)nEventInfo.size() - options.get_nfile();
				nEventInfo.erase(nEventInfo.begin(), nEventInfo.begin()+Ndiscard);
				if(options.get_flag_CentralityCut()){
					vector<EbyeInfo> eBye_picked;
					for (int k = 0; k<(int)nEventInfo.size(); k++){
						eBye_picked.push_back(eBye[nEventInfo[k].order_reading()]);
						//order_reading() is only used here for shuffling. This should be the nth events which is read by this event and archived as nth component of eBye.
				    }
				    eBye=eBye_picked;
				    std::vector<EbyeInfo>().swap(eBye_picked);
				}
				cout << ":D After Shuffling the reading events!  nEventInfo.size():" << (int) nEventInfo.size() 
				    << ", eBye.size():" << (int)eBye.size() 
				    << endl;

		    }else if((int)nEventInfo.size()>(int)options.get_nfile()){
			    cout << "Input has " << (int)nEventInfo.size() << " events but I am going to analyse " << (int)options.get_nfile() << " events." << endl;
			    std::vector<Container::EventInfo> nEventInfo_picked(nEventInfo.begin(), nEventInfo.begin()+options.get_nfile());
			    if(options.get_flag_CentralityCut()){
				    vector<EbyeInfo> eBye_picked;
				    for (int k = 0; k<(int)nEventInfo_picked.size(); k++){
					    eBye_picked.push_back(eBye[nEventInfo_picked[k].order_reading()]);
					    //order_reading() is only used here for shuffling. This should be the nth events which is read by this event and archived as nth component of eBye.
				    }
				    eBye=eBye_picked;
				    std::vector<EbyeInfo>().swap(eBye_picked);
			    }
			    nEventInfo=nEventInfo_picked;
			    std::vector<Container::EventInfo>().swap(nEventInfo_picked);
		    }

			return;
		}



	public:

		Analysis(Options options_in, LogSettings log_in): PrintCounter(0),options(options_in), log(log_in){

			//Random number for thermal parton sampling.
			//-------------------------------------------
			srand((unsigned) time(NULL));
			rndom = std::make_shared<Rndom>(111);
			//rndom = make_shared<Rndom>(rand());

			ms = std::make_shared<Message>();
			uf = std::make_shared<Util_func>(this->rndom);
			if(options.get_flag_set_Ncoeff()){
				options.ih.N_coeff = options.get_Ncoeff();
			}
			read = std::make_shared<ReadIn>(this->ms, this->options);
			fill = std::make_shared<Fill>(this->ms, this->options,  this->uf);
			stat = std::make_shared<Stat>(this->ms, this->options,  this->uf);
			write = std::make_shared<Write>(this->ms, this->options, this->uf);
			if(constants::MODE.find("timelapse")!=std::string::npos) this->PrintCounter=constants::PrintCounterTL;
			else this->PrintCounter=constants::PrintCounter;
			if(options.get_flag_BSTR()) {
				bstr = std::make_shared<BSTR>(this->options);
				for(int iCent=0; iCent<(int)options.name_cent.size(); iCent++){
							ct_ALL.emplace_back(this->options.get_flag_SB_CMS());
				}
			}

			this->ana();

		};
		~Analysis(){};






		int ana(){


			//To get bootstrap error
			//======================
			for(int iBSTR=0; iBSTR<options.get_nBSTR(); iBSTR++){
				options.set_current_iBSTR(iBSTR);
				if(options.get_flag_BSTR()){
					std::cout << "BSTRBSTRBSTRBSTRBSTRBSTRBSTRBSTRBSTR" << iBSTR << std::endl;
					std::cout << ":D iBSTR -- " << iBSTR << std::endl;
				}

				//Start Centrality Cut.
				//-----------------------
				int nCent=1;
				std::vector<EbyeInfo> eBye_CentCut;
				std::vector <Container::EventInfo> nEventInfo;//Archive all event info
				if(options.get_flag_CentralityCut() || options.get_flag_vs_Multi()){
					CentralityCut CentCut(eBye_CentCut, nEventInfo, options, this->rndom);
					if(options.get_flag_BSTR()){
						if(options.get_current_iBSTR()==0){
							//Archiving data that I initially read.
							//For BSTR, I pick up nevents from (int)nEventInfo_for_BSTR.size().
							//I need to archieve the information because nEventInfo size will be reduced to nevents
							//for analysis.
							//===============================================================
							this->nEventInfo_for_BSTR=CentCut.get_nEventInfo();
							this->eBye_CentCut_for_BSTR=CentCut.get_eBye_CentCut();
						}else{
							nEventInfo = this->nEventInfo_for_BSTR;
							eBye_CentCut = this->eBye_CentCut_for_BSTR;
						}
					}

					this->shuffling_or_takeout_nevents(nEventInfo, eBye_CentCut);
					CentCut.ClassifyCentrality(eBye_CentCut);

					if(!options.get_flag_vs_Multi()) nCent=(int)options.name_cent.size();
				}else if(!options.get_flag_CentralityCut() && this->options.get_flag_SB_CMS()){
					CentralityCut CentCut(eBye_CentCut, nEventInfo, options, this->rndom);
				}else{

					//Read binary filne and archive all
					//minijets info from all events.
					//No centrality cut but since it is one file input, 
					//I need to read all events before event loop.
					//================================
					if(options.get_flag_EKRTformat() && options.get_flag_EKRTbinary()){
						std::vector<EbyeInfo> dmmy;
						if(!options.get_flag_BSTR() || (options.get_flag_BSTR() && options.get_current_iBSTR()==0)){
							read->readEKRTbinary(nEventInfo, dmmy);
							//copy nEVentInfo for BSTR
							if(options.get_flag_BSTR()){
								if(options.get_current_iBSTR()==0){
									this->nEventInfo_for_BSTR = nEventInfo;
								}else{
									nEventInfo = this->nEventInfo_for_BSTR;
								}
							}
						}

						this->shuffling_or_takeout_nevents(nEventInfo, dmmy);

					}//TODO: implement shuffling for nonEKRTbinary format as well.

				}


				//Centrality Loop
				for(int iCent=0; iCent<nCent; iCent++){

					if(options.get_flag_CentralityCut()){
						cout << ":)Start analyzing centrality " << options.name_cent[iCent] << "." << endl;
					}

					//Archive event numbering.
					//========================
					if(this->options.get_flag_SB_CMS()){
						this->rndom->set_flag_CentCut(options.get_flag_CentralityCut());
						this->rndom->Archive_iEv_Cent(iCent, eBye_CentCut);
						this->fill->nextCent(this->rndom);
					}


					auto ct = std::make_shared<Container>(this->options.get_flag_SB_CMS());



					//Event Loop
					//==========
					int EV_Count=0;
					int pct=0;
					for(int i=options.get_beginfile(); i<options.get_nfile(); ++i){
						if(!(i%this->PrintCounter)) ms->read(i);
						if(fabs((double)i/(double)options.get_nfile()-pct*0.1)<constants::SMALL){
							cout << ":D Analysing... " << pct*10 << "\% " << endl;
							pct++;
						}
						ct->CountEv++;
						std::stringstream ss;
						ss << options.get_dir_name() << "/" << options.get_f_name() << setw(9) << setfill('0') << i << "/" << options.get_ext_name();
						std::string inputpath = ss.str();


						//Centrality cut.
						//---------------
						if(options.get_flag_CentralityCut()){
							if(eBye_CentCut[EV_Count].get_V0M_class()!=iCent) {
								EV_Count++;
								continue;
							}
						}

						//Obtain event info for TimeLapse
						//-------------------------------
						double weight_TL=(options.get_flag_CentralityCut() && constants::MODE.find("timelapse")!=std::string::npos)? eBye_CentCut[EV_Count].weight:1.0;
						double dNdeta_fill=(options.get_flag_vs_Multi())? eBye_CentCut[EV_Count].multiplicity:0.0;
						double Nch_fill=(options.get_flag_vs_Multi())? eBye_CentCut[EV_Count].Nch:0.0;
						int bin = (options.get_flag_vs_Multi())?eBye_CentCut[EV_Count].get_V0M_class():0.0;
						EV_Count++;


						if(options.get_flag_vs_Multi() && bin<0) continue;
						if(!options.get_flag_CentralityCut() && !options.get_flag_vs_Multi()){
							if(options.get_flag_INEL_lg_0()){
								eByeInSettings ebe;
								eByeInSettings::eByeMulti Multi(options, inputpath, this->rndom);
								weight_TL=Multi.ebye.weight;
								if(!Multi.ebye.multiplicity_INEL_lg_0) continue;
								ebe.print_eByeInfo(i,Multi);
							}else if(options.get_flag_3outof3_trigger()){
								eByeInSettings ebe;
								eByeInSettings::eByeMulti Multi(options, inputpath, this->rndom);
								weight_TL=Multi.ebye.weight;
								if(!Multi.ebye.trig_3outof3) continue;
								ebe.print_eByeInfo(i,Multi);
							}else if(options.get_flag_2outof3_trigger()){
								eByeInSettings ebe;
								eByeInSettings::eByeMulti Multi(options, inputpath, this->rndom);
								weight_TL=Multi.ebye.weight;
								if(!Multi.ebye.trig_2outof3) continue;
								ebe.print_eByeInfo(i,Multi);
							}else if(options.get_flag_ATLAS_cut()){
								eByeInSettings ebe;
								eByeInSettings::eByeMulti Multi(options, inputpath, this->rndom);
								weight_TL=Multi.ebye.weight;
								if(!Multi.ebye.ATLAS_cut) continue;
								ebe.print_eByeInfo(i,Multi);
							}else if(options.get_flag_VZEROAND_trigger()){
								eByeInSettings ebe;
								eByeInSettings::eByeMulti Multi(options, inputpath, this->rndom);
								weight_TL=Multi.ebye.weight;
								if(!Multi.ebye.trig_VZEROAND) continue;
								ebe.print_eByeInfo(i,Multi);
							}else{
								//eByeInSettings ebe;
								//eByeInSettings::eByeMulti Multi(options, inputpath, this->rndom);
								//weight_TL=Multi.ebye.weight;
								//ebe.print_eByeInfo(i,Multi);
							}
						}





						//Read events
						//===========
						if(!options.get_flag_CentralityCut() && !options.get_flag_EKRTbinary()){
							if(constants::MODE.find("JET_PRAC")!=std::string::npos){
								if(!read->read_jetinfo(inputpath, ct->EVENTINFO)) continue;
							}else{
								if(constants::MODE.find("timelapse")!=std::string::npos){
									std::stringstream ssTL;
									ssTL << options.get_dir_name() << "/" << options.get_f_name() << setw(9) << setfill('0') << i << "/";
									if(!read->readTimeLapse(ssTL.str(), ct->EVENTINFO, weight_TL)) continue;
								}else if(constants::MODE.find("readXY")!=std::string::npos){
									if(!read->readXY(inputpath, ct->EVENTINFO)) continue;
								}else if(options.get_flag_EKRTformat()){
									if(!read->readEKRT(inputpath, ct->EVENTINFO)) continue;
								}else{
									EbyeInfo dmmy;
									if(!read->read(inputpath, ct->EVENTINFO, dmmy)) continue;
								}
							}
						}else{
							ct->EVENTINFO = nEventInfo[i];
						}


						//Filling
						//=======
						if(constants::MODE.find("cumulant_multi")!=std::string::npos){
							if(options.get_flag_2subevent()){ 
								if(options.get_flag__4particle()){
									fill->fill_vn4multi_2sub(ct); 
								}else{		
									fill->fill_vnmulti_2sub(ct); 
								}              
							}else if(options.get_flag_3subevent()){
								if(options.get_flag__4particle()){
									fill->fill_vn4multi_3sub(ct); 
								}else{		
									cout << "ERROR :( Option is wrong. --threesub should be used with --4particle." << endl;
									exit(1);
								}              
							}else{
								if(options.get_flag__4particle()){
									//TODO: Currently Nch_fill is added like this.. but need to integrate information into ct.
									fill->fill_vn4multi(ct, Nch_fill);
								}else{
									fill->fill_vnmulti(ct);
								}
							}
						}else if(constants::MODE.find("cumulant_pt")!=std::string::npos){
							if(options.get_flag_2subevent()) fill->fill_vnpt_2sub(ct); 
							else fill->fill_vnpt(ct);
						}else if(constants::MODE.find("cumulant_eta")!=std::string::npos){
							fill->fill_vneta(ct);
						}else if(constants::MODE.find("twopc2D")!=std::string::npos){
							if(options.get_flag_SB_CMS()){
								fill->fill_twopc_B_CMS(ct, eBye_CentCut);
								fill->fill_twopc(ct);
							}else{
								fill->fill_twopc(ct);
							}
						}else if(constants::MODE.find("twopc1D")!=std::string::npos){
							if(options.get_flag_tagged()) fill->fill_twopc1D_tagged(ct);
							else fill->fill_twopc1D(ct);
						}else if(constants::MODE.find("twopcInteg")!=string::npos){
							fill->fill_twopc1D_taggedInteg(ct);
						}else if(constants::MODE.find("JET_PRAC")!=std::string::npos){
							fill->fill_jets(ct);
						}else if(constants::MODE.find("Rt_spectra")!=string::npos){
							fill->fill_Rt(ct);
						}else if(constants::MODE.find("Rt_yield")!=string::npos){
							fill->fill_RtYield(ct);
						}else if(constants::MODE.find("timelapse")!=string::npos){
							fill->fill_TimeLapse(ct);
						}else if(constants::MODE.find("twodm")!=string::npos){
							fill->fill_2D(ct);
						}else{
							fill->fill(ct, dNdeta_fill, bin);
						}

					}//Event loop


					if(constants::MODE.find("cumulant_eta")!=string::npos || constants::MODE.find("cumulant_multi")!=string::npos || constants::MODE.find("cumulant_pt")!=string::npos){ 
						stat->stat_flow(ct);
					}else if(constants::MODE.find("JET_PRAC")!=string::npos){
						stat->stat_jets(ct);
					}else if(constants::MODE.find("twopc2D")!=string::npos){
						stat->stat_twopc(ct);
					}else if(constants::MODE.find("twopc1D")!=string::npos){
						stat->stat_twopc1D(ct);
					}else if(constants::MODE.find("twodm")!=string::npos){
						stat->stat_twopc(ct);
					}else if(constants::MODE.find("Rt_spectra")!=string::npos){
						stat->stat_Rt(ct);
					}else if(constants::MODE.find("Rt_yield")!=string::npos){
						stat->stat_RtYield(ct);
					}else{
						stat->stat(ct);
						//Error bar with bootstrap is available only here.
						if(options.get_flag_BSTR()){
							bstr->fill_iBSTR(iCent, ct, this->ct_ALL[iCent]);
						}
					}


					//Making output directory name
					//-----------------------------
					if(!options.get_flag_BSTR()){
						std::string generated_directory_name;
						if(options.get_flag_CentralityCut()){
							if(iCent==0){
								generated_directory_name=uf->get_output_directory(options.get_out_directory_name());
								uf->make_output_directory(generated_directory_name);
							}
							generated_directory_name=uf->get_output_directory(options.get_out_directory_name()+="/Cent"+options.name_cent[iCent]);
						}else generated_directory_name=uf->get_output_directory(options.get_out_directory_name()); uf->make_output_directory(generated_directory_name);


						//Output
						//------
						if(constants::MODE.find("Rt_yield")!=string::npos){
							if(!write->write_RtYield(generated_directory_name, ct)) return 1;
						}else {
							if(!write->write(generated_directory_name, ct)) return 1;
						}
						if(!log.archive_settings(generated_directory_name)) return 1;
					}

				std::cout << std::endl;
				}//Centrality loop

				std::cout << std::endl;
			}//Bootstrap loop

			if(options.get_flag_BSTR()){
				for(int iCent=0; iCent<(int)options.name_cent.size(); iCent++){
					bstr->stat_iBSTR(this->ct_ALL[iCent]);
					std::string generated_directory_name;
					if(options.get_flag_CentralityCut()){
						if(iCent==0){
							generated_directory_name=uf->get_output_directory(options.get_out_directory_name());
							uf->make_output_directory(generated_directory_name);
						}
						generated_directory_name=uf->get_output_directory(options.get_out_directory_name()+="/Cent"+options.name_cent[iCent]);
					}else generated_directory_name=uf->get_output_directory(options.get_out_directory_name()); uf->make_output_directory(generated_directory_name);

					//Output
					//------
					if(!write->write_BSTR(generated_directory_name, this->ct_ALL[iCent])) return 1;
					if(!log.archive_settings(generated_directory_name)) return 1;
				}

			}

				return 0;
			}



		};
