#ifndef SETTINGS
#define SETTINGS
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <math.h> 
#include <vector>
#include <sstream>
#include <fstream>
#include <memory>
#include "Constants.h"
#include "Message.h"
#include "LogSettings.h"
#include "EbyeInfo.h"
#include "Util_func.h"
#include "Options.h"

using std::cout;
using std::endl;

class Settings{

	private:

		shared_ptr<Message> ms;


	public:



		Options options;
		LogSettings log;

		Settings(int argc, char* argv[]){
			cout << __FILE__ << " Settings constructor." << endl;

			ms = std::make_shared<Message>();
			init(argc, argv);
			consistency_check();
			if(options.get_xaxis_type()==3){
				this->log.set_centrality_cut(this->options.get_flag_CentralityCut());
				this->log.save_BinSettings(this->options.xMin_cstm, this->options.xMax_cstm);
			}
		};
		~Settings(){};


		void init(int argc, char* argv[]){

			//if(!ms->enough_argument(argc)) exit (1);
			cout << "argc: " << argc << endl;

			for(int i=1; i<argc; i++) {
				log.options.push_back(argv[i-1]);
				log.options.push_back(argv[i]);
				cout << __FILE__ << "  " << argv[i] << endl;
				if(!strcmp(argv[i],"-dir")){options.set_dir_name ( argv[i+1]); i++;}///FOR OUTPUT DIR.NAME.
				else if(!strcmp(argv[i],"-f")){options.set_f_name ( argv[i+1]); i++;}
				else if(!strcmp(argv[i],"-ext")){options.set_ext_name ( argv[i+1]); i++;}
				else if(!strcmp(argv[i],"-n")){options.set_nfile (atoi(argv[i+1]));}
				else if(!strcmp(argv[i],"-from")){options.set_beginfile (atoi(argv[i+1])); options.set_flag_specify_startingfile();}
				else if(!strcmp(argv[i],"-to")){options.set_nfile (atoi(argv[i+1])); options.set_flag_specify_startingfile();}
				else if(!strcmp(argv[i],"-outdir")){options.set_out_directory_name(argv[i+1]);}
				else if(!strcmp(argv[i],"--multip_cut_more_than")){options.set_flag_multiplicity_cut_more_than(atof(argv[i+1]));}
				else if(!strcmp(argv[i],"--multip_cut_less_than")){options.set_flag_multiplicity_cut_less_than(atof(argv[i+1]));}
				else if(!strcmp(argv[i],"-modeTL")){options.set_modeTL(atoi(argv[i+1]));}
				else if(!strcmp(argv[i],"--ID")){options.set_specify_ID(atoi(argv[i+1]));}
				else if(!strcmp(argv[i],"-valTL")){options.set_valTL(argv[i+1]);}
				else if(!strcmp(argv[i],"-look_at_xTL")){options.set_at_xTL(atof(argv[i+1]));}
				else if(!strcmp(argv[i],"-look_at_yTL")){options.set_at_yTL(atof(argv[i+1]));}
				else if(!strcmp(argv[i],"-look_at_etaTL")){options.set_at_etaTL(atof(argv[i+1]));}
				else if(!strcmp(argv[i],"--high_pt_mode")){options.set_flag_high_pt_mode();}
				else if(!strcmp(argv[i],"--shuffle")){options.set_flag_shuffle();}
				else if(!strcmp(argv[i],"--BSTR")){options.set_flag_BSTR(); options.set_nBSTR(atoi(argv[i+1]));}
				else if(!strcmp(argv[i],"--npickup_BSTR")){options.set_npickupBSTR(atoi(argv[i+1]));}
				else if(!strcmp(argv[i],"--CentralityCutsumEt")){options.set_flag_CentralityCutsumEt();}
				else if(!strcmp(argv[i],"--pPb_cm")){options.set_flag_pPb_cm_calculation();}
				else if(!strcmp(argv[i],"--long_range_cut_type")){options.set_mid_rapidity_cut_type(atoi(argv[i+1]));}//0 or 1
				else if(!strcmp(argv[i],"--parton")){options.set_parton_level_hist();}
				else if(!strcmp(argv[i],"--nozeros")){options.set_flag_zerofill();}
				else if(!strcmp(argv[i],"--vs_Multi")){options.set_flag_vs_Multi(atoi(argv[i+1]));}//Need to specify collision system in the case that centrality determination is required.
				else if(!strcmp(argv[i],"--xaxis3_input")){options.set_axis3_input(argv[i+1]);}
				//else if(!strcmp(argv[i],"--range")){options.set_xmax(atof(argv[i+1]));}///FOR x RANGE.
				else if(!strcmp(argv[i],"--yshift")){options.dlty = atof(argv[i+1]); options.set_rapidity_shift_hist();} ///FOR rapidity shift.
				else if(!strcmp(argv[i],"--xaxis")){options.set_xaxis_type(atoi(argv[i+1]));}//For output xaxis
				else if(!strcmp(argv[i],"--INEL_lg_0")){options.set_INEL_lg_0();}
				else if(!strcmp(argv[i],"--trig_3outof3")){options.set_trig_3outof3();}
				else if(!strcmp(argv[i],"--trig_2outof3")){options.set_trig_2outof3();}
				else if(!strcmp(argv[i],"--trig_VZEROAND")){options.set_trig_VZEROAND();}
				else if(!strcmp(argv[i],"--ATLAS_cut")){options.set_ATLAS_cut();}
				else if(!strcmp(argv[i],"--EKRTformat")){options.set_EKRTformat();}
				else if(!strcmp(argv[i],"--EKRTbinary")){options.set_EKRTbinary();}
				else if(!strcmp(argv[i],"--only_core")){options.set_flag_only_core();}
				else if(!strcmp(argv[i],"--only_corona")){options.set_flag_only_corona();}
				else if(!strcmp(argv[i],"--only_core_triggers")){options.set_flag_only_core_triggers();}
				else if(!strcmp(argv[i],"--only_corona_triggers")){options.set_flag_only_corona_triggers();}
				else if(!strcmp(argv[i],"--only_core_associates")){options.set_flag_only_core_associates();}
				else if(!strcmp(argv[i],"--only_corona_associates")){options.set_flag_only_corona_associates();}
				else if(!strcmp(argv[i],"--print_dndmt")){options.set_flag_print_dndmt();}
				else if(!strcmp(argv[i],"--print_eBye")){options.set_flag_print_eBye();}
				else if(!strcmp(argv[i],"--CentralityCut")){options.set_flag_CentralityCut(atoi(argv[i+1]));}
				else if(!strcmp(argv[i],"--CentralityCut_dir")){options.set_flag_Specify_dir_for_CentralityCut(); options.set_dir_name_CentCut(argv[i+1]);}
				else if(!strcmp(argv[i],"--CentralityCut_f")){options.set_flag_Specify_f_for_CentralityCut(); options.set_f_name_CentCut(argv[i+1]);}
				else if(!strcmp(argv[i],"--CentralityCut_ext")){options.set_flag_Specify_ext_for_CentralityCut(); options.set_ext_name_CentCut(argv[i+1]);}
				else if(!strcmp(argv[i],"--twosub")){options.set_flag_2subevent();}
				else if(!strcmp(argv[i],"--threesub")){options.set_flag_3subevent();}
				else if(!strcmp(argv[i],"--HI")){options.set_flag_HI();}
				else if(!strcmp(argv[i],"--4particle")){options.set_flag__4particle(); cout << ";) Four particle cumulant calculation started!" << endl;}
				else if(!strcmp(argv[i],"--SB_CMS")){options.set_flag_SB_CMS(); cout << ";) 2PC correlation in 2D calculation. S/B with CMS definition will be performed." << endl;}
				else if(!strcmp(argv[i],"--2PCfull")){options.set_flag__2PCfull(); cout << ";) 2PC correlation in 1D calculation started with --2PCfull." << endl;}
				else if(!strcmp(argv[i],"--2PCnearside")){options.set_flag__2PCnearside(); cout << ";) 2PC correlation in 1D calculation started with --2PCnearside." << endl;}
				else if(!strcmp(argv[i],"--2PCout")){options.set_flag__2PCout(); cout << ";) 2PC correlation in 1D calculation started with --2PCout." << endl;}
				else if(!strcmp(argv[i],"--tagged")){options.set_flag_tagged(); cout << ";) 2PC TAGGED calculation will be performed." << endl;}
				else if(!strcmp(argv[i],"--setNcoeff")){options.set_flag_set_Ncoeff(); options.Ncoeff=(double)atoi(argv[i+1]); cout << ";) v"<< options.Ncoeff << " calculation is performed. " << endl;}
				else { 
					string opt_in(argv[i]);
					if(opt_in.find('-')==string::npos)continue;
					cout << "ERROR:( There is no such an option: " << opt_in << endl; 
					exit(1);
				}
			}
		}

void consistency_check(){

	if(constants::MODE.find("twopcInteg")!=string::npos || constants::MODE.find("twopc1D")!=std::string::npos){
		if(!options.get_flag__2PCfull() && !options.get_flag__2PCnearside() && !options.get_flag__2PCout()) {
			cout << ":( ERROR. Option --2PCfull, --2PCnearside, or --2PCout is needed. " << endl;
			exit(1);
		}
		if((options.get_flag__2PCfull() && options.get_flag__2PCnearside()) || (options.get_flag__2PCfull() && options.get_flag__2PCout()) || (options.get_flag__2PCnearside() && options.get_flag__2PCout())) {
			cout << ":( ERROR. Use only one option from --2PCfull, --2PCnearside, or --2PCout. " << endl;
			exit(1);
		}
	}else{
		if(options.get_flag_tagged()) cout << ":o WARNING. --tagged option works when TWOPC1D is activated." << endl;
	}

	if(constants::MODE.find("twopcInteg")!=string::npos){
		options.set_flag_tagged();
		cout << ":) Mode TWOPCINTEG is called. Option --tagged is forced to be used. " << endl;
	}

        if(options.get_flag_specify_startingfile() && (options.get_nfile()<options.get_beginfile())){
		cout << ":( ERROR. Specify reading file number with -from and -to. -from " << options.get_beginfile() << " -to " << options.get_nfile() << endl;
		exit(1);
	}

	if(constants::MODE.find("twopc2D")!=string::npos){
		if(options.get_flag_only_core()){
				options.set_flag_only_core_triggers();
				options.set_flag_only_core_associates();
		}else if(options.get_flag_only_corona()){
				options.set_flag_only_corona_triggers();
				options.set_flag_only_corona_associates();
		}
	}

	if (constants::MODE.find("MeanptPID")!=string::npos && options.get_specify_ID() == -10000){
		cout << "ERROR:( Specify ID for meanpt PID calculation." << endl;
		exit(1);
	}

}



};


class eByeInSettings{


	private:

				ofstream ofs_eBye;

	public:

				//Functions for eBye info
				//----------------------------
				class eByeMulti{
					public:
						EbyeInfo ebye;
						eByeMulti(Options option_in, const std::string inputpath, shared_ptr<Rndom>& rndom_in):rndom(rndom_in){
							EbyeInfo ebye_;
							//this->get_EbyeInfo(option_in, inputpath, ebye_);
							ebye=ebye_;
						}
					private:
						shared_ptr<Rndom>& rndom;
						std::vector<Container::EventInfo> nEventInfo;
						void get_EbyeInfo(Options options, const std::string inputpath, EbyeInfo& ebye){
							//auto utl_ = make_shared<Util_func>(this->rndom);
							//utl_->get_EbyeInfo(inputpath, nEventInfo, ebye, options);
						}

				};
				void open_eBye_output(const std::string output){
					std::string output_eBye=output+"/"+constants::fname_eByeInfo;
					ofs_eBye.open(output_eBye.c_str());  
				}
				void print_eByeInfo(const int i, const eByeMulti Multi){
					ofs_eBye << i << "  " << Multi.ebye.multiplicity << "  " << Multi.ebye.multiplicity_V0M << "  " << Multi.ebye.weight << endl;
				}


				eByeInSettings(){};
				~eByeInSettings(){};


};
#endif
