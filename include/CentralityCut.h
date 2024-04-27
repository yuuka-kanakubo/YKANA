#ifndef CENTRALITYCUT
#define CENTRALITYCUT
#include "EbyeInfo.h"
#include "Classification.h"
#include "ReadIn.h"

class CentralityCut{


public:

	CentralityCut(vector<EbyeInfo>& eBye_in, vector<Container::EventInfo>& nEventInfo, Options& options_in, shared_ptr<Rndom>& rndom_in):eBye(eBye_in), options(options_in), rndom(rndom_in){
		cout << ":)Start analysis for centrality cut." << endl;

			if(options.get_collision_type()==1){
				for(int i=0; i<7; ++i){ 
					options.val_cent.push_back(constants::val_cent_pPb[i]); 
					options.name_cent.push_back(constants::name_cent_pPb[i]);
				}
			}else if(options.get_collision_type()==2){
				for(int i=0; i<11; ++i){ 
					options.val_cent.push_back(constants::val_cent_PbPb[i]); 
					options.name_cent.push_back(constants::name_cent_PbPb[i]);
				}
			}else if(options.get_collision_type()==4){
				for(int i=0; i<7; ++i){ 
					options.val_cent.push_back(constants::val_cent_PbPb_wide[i]); 
					options.name_cent.push_back(constants::name_cent_PbPb_wide[i]);
				}
			}else if(options.get_collision_type()==3){
				for(int i=0; i<10; ++i){ 
					options.val_cent.push_back(constants::val_cent_pp[i]); 
					options.name_cent.push_back(constants::name_cent_pp[i]);
				}

			}else if(options.get_collision_type()==8){
				for(int i=0; i<10; ++i){ 
					options.val_cent.push_back(constants::val_cent_original_narrow[i]); 
					options.name_cent.push_back(constants::name_cent_original_narrow[i]);
				}

			}else if(options.get_collision_type()==9){
				for(int i=0; i<4; ++i){ 
					options.val_cent.push_back(constants::val_cent_original[i]); 
					options.name_cent.push_back(constants::name_cent_original[i]);
				}

			}else if(options.get_collision_type()==101){
				for(int i=0; i<13; ++i){ 
					options.val_cent.push_back(constants::val_NtrkClass_pp[i]); 
					options.name_cent.push_back(constants::name_NtrkClass_pp[i]);
				}

			}else if(options.get_collision_type()==10){
				for(int i=0; i<11; ++i){ 
					options.val_cent.push_back(constants::val_cent_pp5[i]); 
					options.name_cent.push_back(constants::name_cent_pp5[i]);
				}

			}else if(options.get_collision_type()==11){
				for(int i=0; i<11; ++i){ 
					options.val_cent.push_back(constants::val_cent_pp5_Xi[i]); 
					options.name_cent.push_back(constants::name_cent_pp5_Xi[i]);
				}

			}else if(options.get_collision_type()==12){
				for(int i=0; i<7; ++i){ 
					options.val_cent.push_back(constants::val_cent_pp5_Omega[i]); 
					options.name_cent.push_back(constants::name_cent_pp5_Omega[i]);
				}

			}else{
				cout << "ERROR:( no such a collision type " << options.get_collision_type() << endl;
				exit(1);
			}


			this->read_events(nEventInfo);


	}

         void ClassifyCentrality(){



			Classification* ev_cl=new Classification(eBye, options);
			if(options.get_collision_type()==101){
				cout << "This is Ntrk classification. No need to sort." << endl;
			}else{
				ev_cl->V0M_classification();
			}
			delete ev_cl;
			cout << ":)Finish centrality classification." << endl;


	 }


         vector<EbyeInfo>& eBye;

private:
	Options& options;
	shared_ptr<Rndom>& rndom;
	shared_ptr<ReadIn> read;
	shared_ptr<Message> ms;



void read_events(vector<Container::EventInfo>& nEventInfo){

	ms = make_shared<Message>();
	read = make_shared<ReadIn>(this->ms, this->options);
	int EV_Counter=0;
	int nfile=options.get_nfile();
	if(options.get_flag_EKRTbinary()) nfile=1;
	for(int i=options.get_beginfile();i<nfile;i++) {

	    if(i%1000==0) cout << ":) Reading " << i << "files." << endl; 
	    ostringstream oss;
	    string dirname = (!options.get_flag_Specify_dir_for_CentralityCut())?  options.get_dir_name(): options.get_dir_name_CentCut();
	    string fname = (!options.get_flag_Specify_f_for_CentralityCut())?  options.get_f_name() : options.get_f_name_CentCut();
	    string extname = (!options.get_flag_Specify_ext_for_CentralityCut())? options.get_ext_name() : options.get_ext_name_CentCut();
	    if(!options.get_flag_zerofill()) oss << dirname << "/" << fname << i << "/" << extname;
	    else if (options.get_flag_EKRTbinary()) oss << options.get_dir_name() << "/" << options.get_ext_name();
	    else oss << dirname << "/" << fname << setw(9) << setfill('0') << i << "/" << extname;
	    string inputpath=oss.str();

	    //This is the case when reading files are ebye.
	    //==========================================
	    if (!options.get_flag_EKRTbinary()) 
	    {
		    EbyeInfo ebye_;
		    Container::EventInfo oneEventInfo;
		    read->read(inputpath, oneEventInfo, ebye_);
		    ebye_.orig_eventNum=EV_Counter;
		    this->eBye.push_back(ebye_);
		    nEventInfo.push_back(oneEventInfo);
		    EV_Counter++;
	    }
	    else// EKRT binary files are input then get vector of eBye[nev] in the following.
	    {
		    read->readEKRTbinary(nEventInfo, eBye);
	    }
	}

return;
}





};
#endif
