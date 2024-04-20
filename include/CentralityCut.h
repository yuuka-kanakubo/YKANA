#ifndef CENTRALITYCUT
#define CENTRALITYCUT
#include "EbyeInfo.h"
#include "Classification.h"

class CentralityCut{

public:

	CentralityCut(vector<EbyeInfo>& eBye_in, Settings::Options& options_in, shared_ptr<Rndom>& rndom_in):eBye(eBye_in), options(options_in), rndom(rndom_in){
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


			this->read_events();


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
	Settings::Options& options;
	shared_ptr<Rndom>& rndom;



void read_events(){

	double rap_shift=0.0;
	if(options.get_hist_rapidity_shift() || options.get_flag_pPb_cm_calculation()){
		rap_shift=(options.get_flag_pPb_cm_calculation())? constants::pPb_rap_shift_from_lab_to_cm:options.dlty;
	}

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

	    auto utl_ = make_shared<Util_func>(this->rndom);

	    //This need to be organized. It would be better to make option class?
	    //===========================
	    if(options.get_flag_only_core())utl_->flag_only_core=true;
	    if(options.get_flag_only_corona())utl_->flag_only_corona=true;
	    if(options.get_flag_only_core_associates())utl_->flag_only_core_associates=true;
	    if(options.get_flag_only_corona_associates())utl_->flag_only_corona_associates=true;

	    //This is the case when reading files are ebye.
	    //==========================================
	    {
		    EbyeInfo ebye_;
		    utl_->get_EbyeInfo_(inputpath, ebye_, rap_shift, options.get_flag_VZEROAND_trigger(), options.get_hist_parton_level(), options.get_collision_type());
		    ebye_.orig_eventNum=EV_Counter;
		    this->eBye.push_back(ebye_);
		    EV_Counter++;
	    }
	    //else EKRT binary files are input then put ebye info below.
	    {
		    utl_->get_EbyeInfo_forAlleventsBinaryEKRTformat(inputpath, this->eBye, rap_shift, options.get_flag_VZEROAND_trigger(), options.get_hist_parton_level(), options.get_collision_type());
	    }
	}

return;
}





};
#endif
