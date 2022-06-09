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

using std::cout;
using std::endl;

class Settings{

	private:

		shared_ptr<Message> ms;


	public:

		//------------------------------------------------------------------------
		class Options{

			public:
				//TODO: MOve this to private
				//-------------------------
				vector<double> val_cent;
				vector<string> name_cent;
				vector<double> xMin_cstm, xMax_cstm;
				double dlty;
				double Ncoeff;
				string axis3_inputf;

			private:

				std::string dir_name;
				std::string f_name;
				std::string ext_name;
				std::string out_directory_name;
				std::string out_fname;
				std::string dir_name_CentCut;
				std::string f_name_CentCut;
				std::string ext_name_CentCut;
				int nfile;
				int beginfile;
				bool specify_startingfile;

				//Centrality Cut
				bool CentralityCut;
				bool Specify_dir_for_CentralityCut;
				bool Specify_f_for_CentralityCut;
				bool Specify_ext_for_CentralityCut;
				bool parton_level;
				bool rapidity_shift;
				bool high_pt_mode;
				bool cut_INEL_lg_0;
				bool trig_3outof3;
				bool trig_2outof3;
				bool trig_VZEROAND;
				bool ATLAS_cut;
				bool only_core;
				bool only_corona;
				bool only_core_triggers;
				bool only_corona_triggers;
				bool only_core_associates;
				bool only_corona_associates;
				bool print_eBye;
				bool print_dndmt;
				bool flag_pPb_cm_calculation;
				bool flag_multiplicity_cut;
				bool flag_vs_Multi;
				bool zerofill;
				bool two_subevent;
				bool three_subevent;
				bool HI;
				bool _4particle;
				bool _SB_CMS;
				bool _2PCfull;
				bool _2PCnearside;
				bool _2PCout;
				bool tagged;
				bool set_Ncoeff;
				unsigned int mid_rapidity_cut_type;
				unsigned int axis_type;
				unsigned int collision_type;
				double long_range_hist_pm;
				double multip_cut_more_than;
				double multip_cut_less_than;
                                int modeTL;//0:xy, 1:xeta
                                std::string valTL;
                                double at_xTL, at_yTL, at_etaTL;
				int ID;

				void GetBinSettings(){

					cout << "Getting Bin Settings..." << endl;

					ifstream in;
					in.open(axis3_inputf.c_str(),ios::in);
					if(!in) {
						cout << __FILE__ << " line:"<< __LINE__ << " WARNING:o unable to open file. " << endl;
						exit(1);
					}
					string templine;
					while(getline(in,templine)) {
						istringstream ist(templine);
						double xmin, xmax;
						ist >> xmin >> xmax;
						this->xMin_cstm.push_back(xmin);
						this->xMax_cstm.push_back(xmax);
					}
				}



			public:

				//Set functions.
				//-----------------------
				void set_dir_name(std::string name){this->dir_name=name;}
				void set_f_name(std::string name){this->f_name=name;}
				void set_ext_name(std::string name){this->ext_name=name;}
				void set_nfile(int n){this->nfile=n;}
				void set_beginfile(int n){this->beginfile=n;}
				void set_flag_specify_startingfile(){this->specify_startingfile=true;}
				void set_dir_name_CentCut(std::string name){this->dir_name_CentCut=name;}
				void set_f_name_CentCut(std::string name){this->f_name_CentCut=name;}
				void set_ext_name_CentCut(std::string name){this->ext_name_CentCut=name;}
				void set_out_directory_name(std::string name){this->out_directory_name=name;}

				void set_parton_level_hist(){parton_level=true;};
				void set_rapidity_shift_hist(){rapidity_shift=true;};
				void set_flag_high_pt_mode(){high_pt_mode=true; cout << ":o HIGH PT MODE is called. Currently PP7TeV events are assumed. Please change a settings if you are analyzing different energy or system." << endl;};
				void set_flag_zerofill(){zerofill=true;};
				void set_axis3_input(const string path){axis3_inputf=path;  this->GetBinSettings();}
				void set_xaxis_type(const int i){this->axis_type=i;}
				void set_INEL_lg_0(){cut_INEL_lg_0=true;};
				void set_trig_3outof3(){trig_3outof3=true;};
				void set_trig_2outof3(){trig_2outof3=true;};
				void set_trig_VZEROAND(){trig_VZEROAND=true;};
				void set_ATLAS_cut(){ ATLAS_cut=true; 
					this->set_longrange_of_hist__plus_minus(constants::w_eta_ATLAS_cut);
					mid_rapidity_cut_type=3;
				};
				void set_flag_multiplicity_cut_more_than(const double multip_cut_more_than_in){
					flag_multiplicity_cut=true;
					multip_cut_more_than=multip_cut_more_than_in;
				};
				void set_flag_multiplicity_cut_less_than(const double multip_cut_less_than_in){
					flag_multiplicity_cut=true;
					multip_cut_less_than=multip_cut_less_than_in;
				};
				void set_longrange_of_hist__plus_minus(const double dy){long_range_hist_pm=dy;};
				void set_mid_rapidity_cut_type(const int cut_type){mid_rapidity_cut_type=cut_type;};
				void set_flag_pPb_cm_calculation(){
					flag_pPb_cm_calculation=true;
					mid_rapidity_cut_type=2;
					cout << ":O mid-rapidity cut is automatically set to be -0.5<y_{cm}<0." << endl;
				}
				void set_flag_vs_Multi(const int i){this->flag_vs_Multi=true; this->collision_type=i;}
				void set_modeTL(const int val){this->modeTL=val;}
				void set_valTL(const std::string val){this->valTL=val;}
				void set_at_xTL(const double val){this->at_xTL=val;}
				void set_at_yTL(const double val){this->at_yTL=val;}
				void set_at_etaTL(const double val){this->at_etaTL=val;}
				void set_flag_only_core(){only_core=true;}
				void set_flag_only_corona(){only_corona=true;}
				void set_flag_only_core_triggers(){only_core_triggers=true;}
				void set_flag_only_corona_triggers(){only_corona_triggers=true;}
				void set_flag_only_core_associates(){only_core_associates=true;}
				void set_flag_only_corona_associates(){only_corona_associates=true;}
				void set_flag_print_eBye(){print_eBye=true;}
				void set_flag_print_dndmt(){this->print_dndmt=true;}
				void set_flag_2subevent(){this->two_subevent=true;}
				void set_flag_3subevent(){this->three_subevent=true;}
				void set_flag_HI(){this->HI=true;}
				void set_flag__4particle(){this->_4particle=true;}
				void set_flag_SB_CMS(){this->_SB_CMS=true;}
				void set_flag__2PCfull(){this->_2PCfull=true;}
				void set_flag__2PCnearside(){this->_2PCnearside=true;}
				void set_flag__2PCout(){this->_2PCout=true;}
				void set_flag_tagged(){this->tagged=true;}
				void set_flag_set_Ncoeff(){this->set_Ncoeff=true;}
				void set_flag_Specify_dir_for_CentralityCut(){this->Specify_dir_for_CentralityCut=true;}
				void set_flag_Specify_f_for_CentralityCut(){this->Specify_f_for_CentralityCut=true;}
				void set_flag_Specify_ext_for_CentralityCut(){this->Specify_ext_for_CentralityCut=true;}
				void set_specify_ID(int i){this->ID=i;}
				void set_flag_CentralityCut(int collision_type_in=-1){
					CentralityCut=true;
					collision_type=collision_type_in;
					if(collision_type==1){
						cout << ":)Centrality cut for pPb." << endl;
					}else if(collision_type==2){
						cout << ":)Centrality cut for PbPb." << endl;
					}else if(collision_type==3){
						cout << ":)Centrality cut for pp." << endl;
					}else if(collision_type==4){
						cout << ":)Centrality cut for PbPb (wide)." << endl;
					}else if(collision_type==8){
						cout << ":)Centrality cut original (narrow)." << endl;
					}else if(collision_type==9){
						cout << ":)Centrality cut original." << endl;
					}else if(collision_type==101){
						cout << ":)Centrality cut CMS Ntrk." << endl;
					}else if(collision_type==10){
						cout << ":)Centrality cut pp 5TeV." << endl;
					}else if(collision_type==11){
						cout << ":)Centrality cut pp 5TeV (Xi)." << endl;
					}else if(collision_type==12){
						cout << ":)Centrality cut pp 5TeV (Omega)." << endl;
					}else{
						cout << "ERROR:( Something wrong with --CentralityCut. Specify appropriate collision type. ex) --CentralityCut 1" << endl;
						cout << "        1: pPb, 2:PbPb, 3:pp, 4:PbPb (wide), 8: original(narrow), 9: original " << endl;
						cout << "        101: pp13TeV (CMS) " << endl;
						exit(1);
					}
				}

				//Get functions.
				//-----------------------
				std::string get_dir_name()const{return this->dir_name;}
				std::string get_f_name()const{return this->f_name;}
				std::string get_ext_name()const{return this->ext_name;}
				std::string get_ext_nameTL()const{
					if(this->modeTL==0) return constants::ext_nameTLxy;
					else if(this->modeTL==1) return constants::ext_nameTLxeta;
					else{

						cout << ":(ERROR Set proper option with -modeTL. 0: xy-plane, 1: xeta-plane." << endl;
						exit(1); 
					}
				}
				std::string get_dir_name_CentCut()const{return this->dir_name_CentCut;}
				std::string get_f_name_CentCut()const{return this->f_name_CentCut;}
				std::string get_ext_name_CentCut()const{return this->ext_name_CentCut;}
				std::string get_out_directory_name()const{return this->out_directory_name;}
				int get_nfile()const{return this->nfile;}
				int get_beginfile()const{return this->beginfile;}
				bool get_flag_specify_startingfile()const{return this->specify_startingfile;}
				bool get_flag_vs_Multi()const{return this->flag_vs_Multi;}
				int get_xaxis_type()const{return axis_type;};
				bool get_hist_parton_level()const{return parton_level;};
				bool get_hist_rapidity_shift()const{return rapidity_shift;};
				bool get_flag_high_pt_mode()const{return high_pt_mode;};
				bool get_flag_zerofill()const{return zerofill;};
				double get_d_longrange_pm()const{return long_range_hist_pm;};
				double get_d_longrange()const{return long_range_hist_pm*2.0;};
				double get_d_longrange_mid_rapidity_cut_type2()const{return std::fabs(constants::pPb_mid_rapidity__bkw-constants::pPb_mid_rapidity__fwd);};
				int get_mid_rapidity_cut_type()const{return mid_rapidity_cut_type;};
				int get_collision_type()const{return collision_type;}
				double get_pPb_mid_rapidity__bkw()const{return constants::pPb_mid_rapidity__bkw;}
				double get_pPb_mid_rapidity__fwd()const{return constants::pPb_mid_rapidity__fwd;}
				bool get_flag_pPb_cm_calculation()const{return flag_pPb_cm_calculation;}
				double get_multiplicity_cut_more_than()const{return multip_cut_more_than;};
				double get_multiplicity_cut_less_than()const{return multip_cut_less_than;};
				double get_Ncoeff()const{return Ncoeff;}
				int get_modeTL(){return this->modeTL;}
				std::string get_valTL(){return this->valTL;}
				double get_at_xTL(){return this->at_xTL;}
				double get_at_yTL(){return this->at_yTL;}
				double get_at_etaTL(){return this->at_etaTL;}
				bool get_flag_multiplicity_cut()const{return flag_multiplicity_cut;};
				bool get_flag_INEL_lg_0()const{return cut_INEL_lg_0;};
				bool get_flag_3outof3_trigger()const{return trig_3outof3;};
				bool get_flag_2outof3_trigger()const{return trig_2outof3;};
				bool get_flag_VZEROAND_trigger()const{return trig_VZEROAND;};
				bool get_flag_ATLAS_cut()const{return ATLAS_cut;};
				bool get_flag_only_core()const{return only_core;} 
				bool get_flag_only_corona()const{return only_corona;}
				bool get_flag_only_core_triggers()const{return only_core_triggers;} 
				bool get_flag_only_corona_triggers()const{return only_corona_triggers;}
				bool get_flag_only_core_associates()const{return only_core_associates;} 
				bool get_flag_only_corona_associates()const{return only_corona_associates;}
				bool get_flag_print_eBye()const{return print_eBye;}
				bool get_flag_CentralityCut() const{return CentralityCut;}
				bool get_flag_print_dndmt()const{return print_dndmt;} 
				bool get_flag_Specify_dir_for_CentralityCut()const{return this->Specify_dir_for_CentralityCut;}
				bool get_flag_Specify_f_for_CentralityCut()const{return this->Specify_f_for_CentralityCut;}
				bool get_flag_Specify_ext_for_CentralityCut()const{return this->Specify_ext_for_CentralityCut;}
				bool get_flag_2subevent()const{return this->two_subevent;}
				bool get_flag_3subevent()const{return this->three_subevent;}
				bool get_flag_HI()const{return this->HI;}
				bool get_flag__4particle()const{return this->_4particle;}
				bool get_flag_SB_CMS()const{return this->_SB_CMS;}
				bool get_flag__2PCfull()const{return this->_2PCfull;}
				bool get_flag__2PCnearside()const{return this->_2PCnearside;}
				bool get_flag__2PCout()const{return this->_2PCout;}
				bool get_flag_tagged()const{return this->tagged;}
				bool get_flag_set_Ncoeff()const{return this->set_Ncoeff;}
				int get_specify_ID(){return this->ID;}




				//Constructor
				//--------------
				Options():
					dlty(0.0),
					dir_name(""),
					f_name("xxx"),
					ext_name("xxx"),
					out_directory_name(constants::default_out_directory_name),
					out_fname(constants::default_out_fname),
					dir_name_CentCut(""),
					f_name_CentCut("xxx"),
					ext_name_CentCut("xxx"),
					nfile(1),
					beginfile(0),
					specify_startingfile(false),
					CentralityCut(false),
					Specify_dir_for_CentralityCut(false),
					Specify_f_for_CentralityCut(false),
					Specify_ext_for_CentralityCut(false),
					parton_level(false),
					rapidity_shift(false),
					high_pt_mode(false),
					cut_INEL_lg_0(false),
					trig_3outof3(false),
					trig_2outof3(false),
					trig_VZEROAND(false),
					ATLAS_cut(false),
					only_core(false),
					only_corona(false),
					only_core_triggers(false),
					only_corona_triggers(false),
					only_core_associates(false),
					only_corona_associates(false),
					print_eBye(false),
					print_dndmt(false),
					flag_pPb_cm_calculation(false),
					flag_multiplicity_cut(false),
					zerofill(true),
					two_subevent(false),
					three_subevent(false),
					HI(false),
					_4particle(false),
					_SB_CMS(false),
					_2PCfull(false),
					_2PCnearside(false),
					_2PCout(false),
					tagged(false),
					set_Ncoeff(false),

					mid_rapidity_cut_type(0),
					axis_type(0),
					collision_type(3),
					long_range_hist_pm(constants::default_midy_pm),
					multip_cut_more_than(constants::multip_cut_more_than),
					multip_cut_less_than(constants::multip_cut_less_than),
					modeTL(0),
					valTL("temp"),
					at_xTL(0.0),
					at_yTL(0.0),
					at_etaTL(0.0),
					ID(-10000)
					{};

		};
		//<--- Class Options


		Options options;
		LogSettings log;

		Settings(int argc, char* argv[]){

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

			for(int i=1; i<argc; i++) {
				log.options.push_back(argv[i]);
				if(!strcmp(argv[i],"-dir")){options.set_dir_name ( argv[i+1]);}///FOR OUTPUT DIR.NAME.
				else if(!strcmp(argv[i],"-f")){options.set_f_name ( argv[i+1]);}
				else if(!strcmp(argv[i],"-ext")){options.set_ext_name ( argv[i+1]);}
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
						eByeMulti(Settings::Options option_in, const std::string inputpath, shared_ptr<Rndom>& rndom_in):rndom(rndom_in){
							EbyeInfo ebye_;
							this->get_EbyeInfo(option_in, inputpath, ebye_);
							ebye=ebye_;
						}
					private:
						shared_ptr<Rndom>& rndom;
						void get_EbyeInfo(Settings::Options options, const std::string inputpath, EbyeInfo& ebye){
							auto utl_ = make_shared<Util_func>(this->rndom);
							double rap_shift=0.0;
							if(options.get_hist_rapidity_shift() || options.get_flag_pPb_cm_calculation()){
								rap_shift=(options.get_flag_pPb_cm_calculation())? constants::pPb_rap_shift_from_lab_to_cm:options.dlty;
							}
							utl_->get_EbyeInfo_(inputpath, ebye, rap_shift, options.get_flag_VZEROAND_trigger(), options.get_hist_parton_level(), options.get_collision_type());
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
