#include "EbyeInfo.h"
#include "Classification.h"

class CentralityCut{

public:

	CentralityCut(vector<EbyeInfo>& eBye_in, Settings::Options& options_in):eBye(eBye_in), options(options_in){
		cout << ":)Start analysis for centrality cut." << endl;
                this->read_events();


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
			}else if(options.get_collision_type()==3){
				for(int i=0; i<10; ++i){ 
					options.val_cent.push_back(constants::val_cent_pp[i]); 
					options.name_cent.push_back(constants::name_cent_pp[i]);
				}

			}else if(options.get_collision_type()==9){
				for(int i=0; i<4; ++i){ 
					options.val_cent.push_back(constants::val_cent_original[i]); 
					options.name_cent.push_back(constants::name_cent_original[i]);
				}

			}else{
				cout << "ERROR:( no such a collision type " << options.get_collision_type() << endl;
				exit(1);
			}



		Classification* ev_cl=new Classification(eBye, options);
                ev_cl->V0M_classification();
                delete ev_cl;
		cout << ":)Finish centrality classification." << endl;
	}


         vector<EbyeInfo>& eBye;

private:
	Settings::Options& options;



void read_events(){

	double rap_shift=0.0;
	if(options.get_hist_rapidity_shift() || options.get_flag_pPb_cm_calculation()){
		rap_shift=(options.get_flag_pPb_cm_calculation())? constants::pPb_rap_shift_from_lab_to_cm:options.dlty;
	}

    for(int i=0;i<options.get_nfile();i++) {

	    if(i%1000==0) cout << ":) Reading " << i << "files." << endl; 
	    ostringstream oss;
	    string dirname = (!options.get_flag_Specify_dir_for_CentralityCut())?  options.get_dir_name(): options.get_dir_name_CentCut();
	    string fname = (!options.get_flag_Specify_f_for_CentralityCut())?  options.get_f_name() : options.get_f_name_CentCut();
	    string extname = (!options.get_flag_Specify_ext_for_CentralityCut())? options.get_ext_name() : options.get_ext_name_CentCut();
	    if(!options.get_flag_zerofill()) oss << dirname << "/" << fname << i << "/" << extname;
	    else oss << dirname << "/" << fname << setw(9) << setfill('0') << i << "/" << extname;
	    string inputpath=oss.str();

	    auto utl_ = make_shared<Util_func>();

	    EbyeInfo ebye_;
	    utl_->get_EbyeInfo_(inputpath, ebye_, rap_shift, options.get_flag_VZEROAND_trigger(), options.get_hist_parton_level());
            ebye_.orig_eventNum=i;
	    this->eBye.push_back(ebye_);

    }


}





};
