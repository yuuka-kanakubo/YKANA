#ifndef OPTIONS_H
#define OPTIONS_H
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <math.h> 
#include <vector>
#include <sstream>
#include <fstream>
#include <memory>

		class Options{

			private:

				class InfoHist{

						double x_max;
						double y_max;
						double x_min;
						double y_min;

						bool isthisOddInt(const double val){

							int val_floored = std::floor(val);
							double decimal = val - (double) val_floored;
							if(fabs(decimal)<constants::TINY){
								if (std::fabs(val_floored%2)<constants::TINY){
									//Then even
									return false;
								}else{
									return true;
								}
							}else{ 
								return false;
							}
						}


					public:

						InfoHist():N_coeff(2.0){
							this->x_max = constants::x_max;
							this->x_min = constants::x_min;
							this->y_max = constants::y_max;
							this->y_min = constants::y_min;
							this->d_x = constants::d_x;
							this->d_y = constants::d_y;

						};
						~InfoHist(){};



						void tailor_hist(Options& options__){

							if((x_max > x_min) && (y_max > y_min) && (d_x<fabs(x_max - x_min))&& (d_y<fabs(y_max - y_min))) {
								std::cout << "Setting up InfoHist. Everything is sane." << std::endl;
							}else{
								std::cout << "ERROR:( in InfoHist. Check the followings: " 
									<< "   x_max:" << x_max
									<< "   x_min:" << x_min
									<< "   y_max:" << y_max
									<< "   y_min:" << y_min
									<< "   d_x:" << d_x 
									<< "   d_y:" << d_y 
									<< std::endl;
								exit(EXIT_FAILURE);
							}
							//Getting edge.
							this->x_edge_max = this->x_max + this->d_x/2.0;
							this->x_edge_min = -1.0 * this->x_edge_max;
							this->y_edge_max = this->y_max + this->d_y/2.0;
							this->y_edge_min = -1.0 * this->y_edge_max;

							if(options__.get_flag_hist_ZeroCentered()){
								//Check if fabs(x_edge_max)/d_x returns int.
								//if((fabs(x_edge_max - x_edge_min)%d_x)<constants::TINY)
								if(this->isthisOddInt(fabs(x_edge_max - x_edge_min)/d_x)){
									std::cout << "(:3 = )3 ? " << __FILE__ << " (" << __LINE__ << ") this->isthisOddInt : true." << fabs(x_edge_max - x_edge_min)/d_x << std::endl;
									return;//passed the criteria.
								}else{
									//Original cell n 
									//Original x_edge_min, max
									//Find closest odd int.
									double orig_n = fabs(x_edge_max - x_edge_min)/d_x;
									int orig_n_floored = std::floor(orig_n);
									std::cout << "(:3 = )3 ? " << __FILE__ << " (" << __LINE__ << ") orig_n_floored " << orig_n_floored << ",   orig_n " << orig_n << std::endl;
									if (std::fabs(orig_n_floored%2)<constants::TINY){
										orig_n_floored++;
										this->x_edge_max = orig_n_floored * this->d_x/2.0;
										this->x_edge_min = -1.0 * this->x_edge_max;
										std::cout << "(:3 = )3 ? " << __FILE__ << " (" << __LINE__ << ") x_edge_max " << x_edge_max << std::endl;
										std::cout << "(:3 = )3 ? " << __FILE__ << " (" << __LINE__ << ") x_edge_min " << x_edge_min << std::endl;
									}else{
										this->x_edge_max = orig_n_floored * this->d_x/2.0;
										this->x_edge_min = -1.0 * this->x_edge_max;
										std::cout << "(:3 = )3 ? " << __FILE__ << " (" << __LINE__ << ") x_edge_max " << x_edge_max << std::endl;
										std::cout << "(:3 = )3 ? " << __FILE__ << " (" << __LINE__ << ") x_edge_min " << x_edge_min << std::endl;
									}
									return;
								}
							}


							return;
						}

						//Maximum value for histgram
						//-------------------------
						double d_x;
						double d_y;
						double N_coeff;
						double x_edge_min;
						double x_edge_max;
						double y_edge_min;
						double y_edge_max;


				};

			public:
				//TODO: MOve this to private
				//-------------------------
				std::vector<double> val_cent;
				std::vector<std::string> name_cent;
				std::vector<double> xMin_cstm, xMax_cstm;
				double dlty;
				double Ncoeff;
				std::string axis3_inputf;

                                InfoHist ih;

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
				bool EKRTformat;
				bool EKRTbinary;
				bool hist_ZeroCentered;
				bool only_core;
				bool only_corona;
				bool only_core_triggers;
				bool only_corona_triggers;
				bool only_core_associates;
				bool only_corona_associates;
				bool print_eBye;
				bool print_dndmt;
				bool flag_pPb_cm2lab;
				bool flag_pPb_lab2cm;
				bool flag_multiplicity_cut;
				bool flag_vs_Multi;
				bool shuffle;
				bool BSTR;
				int nBSTR;
				int current_iBSTR;
				std::string obs_type;
				bool flag_sortsumEt;
				bool flag_sortV0A;
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

					std::cout << "Getting Bin Settings..." << std::endl;

					std::ifstream in;
					in.open(axis3_inputf.c_str(),std::ios::in);
					if(!in) {
						std::cout << __FILE__ << " line:"<< __LINE__ << " WARNING:o unable to open file. " << std::endl;
						exit(1);
					}
					std::string templine;
					while(getline(in,templine)) {
						std::istringstream ist(templine);
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
				void set_flag_shuffle(){this->shuffle=true;}
				void set_flag_BSTR(){this->BSTR=true;}
				void set_nBSTR(int n){this->nBSTR=n;}
				void set_current_iBSTR(int n){this->current_iBSTR=n;}
				void set_flag_sortsumEt(){this->flag_sortsumEt=true;}
				void set_flag_sortV0A(){this->flag_sortV0A=true;}
				void set_obs_type(std::string str_type){

					if(str_type.find("detdy")!=std::string::npos){
						this->obs_type="detdy";
					}else if(str_type.find("dedy")!=std::string::npos){
						this->obs_type="dedy";
					}else if(str_type.find("dndy")!=std::string::npos){
						this->obs_type="dndy";
					}else if(str_type.find("dndeta")!=std::string::npos){
						this->obs_type="dndeta";
					}else if(str_type.find("dndpt")!=std::string::npos){
						this->obs_type="dndpt";
					}else if(str_type.find("dndcoordx")!=std::string::npos){
						this->obs_type="dndcoordx";
					}else if(str_type.find("dndcoordy")!=std::string::npos){
						this->obs_type="dndcoordy";
					}else if(str_type.find("dndphi")!=std::string::npos){
						this->obs_type="dndphi";
					}else{
						std::cout << __FILE__ << "(" << __LINE__ << ")" << "ERROR :( Unknown option " << str_type << std::endl; 
						exit(EXIT_FAILURE);
					}


				}

				void set_dir_name_CentCut(std::string name){this->dir_name_CentCut=name;}
				void set_f_name_CentCut(std::string name){this->f_name_CentCut=name;}
				void set_ext_name_CentCut(std::string name){this->ext_name_CentCut=name;}
				void set_out_directory_name(std::string name){this->out_directory_name=name;}

				void set_parton_level_hist(){parton_level=true;};
				void set_rapidity_shift_hist(){rapidity_shift=true;};
				void set_flag_high_pt_mode(){high_pt_mode=true; std::cout << ":o HIGH PT MODE is called. Currently PP7TeV events are assumed. Please change a settings if you are analyzing different energy or system." << std::endl;};
				void set_flag_zerofill(){zerofill=true;};
				void set_axis3_input(const std::string path){axis3_inputf=path;  this->GetBinSettings();}
				void set_xaxis_type(const int i){this->axis_type=i;}
				void set_INEL_lg_0(){cut_INEL_lg_0=true;};
				void set_trig_3outof3(){trig_3outof3=true;};
				void set_trig_2outof3(){trig_2outof3=true;};
				void set_trig_VZEROAND(){trig_VZEROAND=true;};
				void set_ATLAS_cut(){ ATLAS_cut=true; 
					this->set_longrange_of_hist__plus_minus(constants::w_eta_ATLAS_cut);
					mid_rapidity_cut_type=3;
				};
				void set_EKRTformat(){ 
					EKRTformat=true;
					std::cout << "Now EKRT format is assumed for input." << std::endl; 
				};
				void set_EKRTbinary(){ 
					EKRTbinary=true;
					std::cout << "Reading EKRT binary files..." << std::endl; 
				};
				void set_hist_ZeroCentered(){ 
					hist_ZeroCentered=true;
					std::cout << "Histgram is going to span [-x, x] with 0 centered. ..." << std::endl; 
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
				void set_flag_pPb_lab2cm(){
					flag_pPb_lab2cm=true;
					mid_rapidity_cut_type=2;
					std::cout << ":O mid-rapidity cut is automatically set to be -0.5<y_{cm}<0." << std::endl;
				}
				void set_flag_pPb_cm2lab(){
					flag_pPb_cm2lab=true;
					mid_rapidity_cut_type=2;
					std::cout << ":O mid-rapidity cut is automatically set to be -0.5<y_{cm}<0." << std::endl;
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
						std::cout << ":)Centrality cut for pPb." << std::endl;
					}else if(collision_type==2){
						std::cout << ":)Centrality cut for PbPb." << std::endl;
					}else if(collision_type==3){
						std::cout << ":)Centrality cut for pp." << std::endl;
					}else if(collision_type==4){
						std::cout << ":)Centrality cut for PbPb (wide)." << std::endl;
					}else if(collision_type==8){
						std::cout << ":)Centrality cut original (narrow)." << std::endl;
					}else if(collision_type==9){
						std::cout << ":)Centrality cut original." << std::endl;
					}else if(collision_type==101){
						std::cout << ":)Centrality cut CMS Ntrk." << std::endl;
					}else if(collision_type==10){
						std::cout << ":)Centrality cut pp 5TeV." << std::endl;
					}else if(collision_type==11){
						std::cout << ":)Centrality cut pp 5TeV (Xi)." << std::endl;
					}else if(collision_type==12){
						std::cout << ":)Centrality cut pp 5TeV (Omega)." << std::endl;
					}else{
						std::cout << "ERROR:( Something wrong with --CentralityCut. Specify appropriate collision type. ex) --CentralityCut 1" << std::endl;
						std::cout << "        1: pPb, 2:PbPb, 3:pp, 4:PbPb (wide), 8: original(narrow), 9: original " << std::endl;
						std::cout << "        101: pp13TeV (CMS) " << std::endl;
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

						std::cout << ":(ERROR Set proper option with -modeTL. 0: xy-plane, 1: xeta-plane." << std::endl;
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
				bool get_flag_shuffle()const{return shuffle;};
				bool get_flag_BSTR()const{return BSTR;};
				int get_nBSTR()const{return nBSTR;};
				int get_current_iBSTR()const{return current_iBSTR;};
				bool get_flag_sortsumEt()const{return flag_sortsumEt;};
				bool get_flag_sortV0A()const{return flag_sortV0A;};
				std::string get_obs_type()const{return obs_type;};
				double get_d_longrange_pm()const{return long_range_hist_pm;};
				double get_d_longrange()const{return long_range_hist_pm*2.0;};
				double get_d_longrange_mid_rapidity_cut_type2()const{return std::fabs(constants::pPb_mid_rapidity__bkw-constants::pPb_mid_rapidity__fwd);};
				int get_mid_rapidity_cut_type()const{return mid_rapidity_cut_type;};
				int get_collision_type()const{return collision_type;}
				double get_pPb_mid_rapidity__bkw()const{return constants::pPb_mid_rapidity__bkw;}
				double get_pPb_mid_rapidity__fwd()const{return constants::pPb_mid_rapidity__fwd;}
				bool get_flag_pPb_cm2lab()const{return flag_pPb_cm2lab;}
				bool get_flag_pPb_lab2cm()const{return flag_pPb_lab2cm;}
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
				bool get_flag_EKRTformat()const{return EKRTformat;};
				bool get_flag_EKRTbinary()const{return EKRTbinary;};
				bool get_flag_hist_ZeroCentered()const{return hist_ZeroCentered;};
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
					EKRTformat(false),
					EKRTbinary(false),
					hist_ZeroCentered(false),
					only_core(false),
					only_corona(false),
					only_core_triggers(false),
					only_corona_triggers(false),
					only_core_associates(false),
					only_corona_associates(false),
					print_eBye(false),
					print_dndmt(false),
					flag_pPb_cm2lab(false),
					flag_pPb_lab2cm(false),
					flag_multiplicity_cut(false),
					shuffle(false),
					BSTR(false),
					nBSTR(1),
					current_iBSTR(0),
					obs_type("dedy"),
					flag_sortsumEt(false),
					flag_sortV0A(false),
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
					{
					};




		};
#endif
