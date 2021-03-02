#include "EbyeInfo.h"
#include "classification.h"

class CentralityCut{

public:

	CentralityCut(vector<EbyeInfo>& eBye_in,const shared_ptr<settings> set_in):eBye(eBye_in),set(set_in){
		cout << ":)Start analysis for centrality cut." << endl;
                this->read_events();


			if(set->collision_type==1){
				for(int i=0; i<7; ++i){ 
					set->val_cent.push_back(constants::val_cent_pPb[i]); 
					set->name_cent.push_back(constants::name_cent_pPb[i]);
				}
			}else if(set->collision_type==2){
				for(int i=0; i<11; ++i){ 
					set->val_cent.push_back(constants::val_cent_PbPb[i]); 
					set->name_cent.push_back(constants::name_cent_PbPb[i]);
				}
			}else if(set->collision_type==3){
				for(int i=0; i<10; ++i){ 
					set->val_cent.push_back(constants::val_cent_pp[i]); 
					set->name_cent.push_back(constants::name_cent_pp[i]);
				}

			}else if(set->collision_type==9){
				for(int i=0; i<4; ++i){ 
					set->val_cent.push_back(constants::val_cent_original[i]); 
					set->name_cent.push_back(constants::name_cent_original[i]);
				}

			}else{
				cout << "ERROR:( no such a collision type " << set->collision_type << endl;
				exit(1);
			}



		classification* ev_cl=new classification(eBye, set);
                ev_cl->V0M_classification();
                delete ev_cl;
		cout << ":)Finish centrality classification." << endl;
	}


         vector<EbyeInfo>& eBye;

private:
	shared_ptr<settings> set;



void read_events(){

	double rap_shift=0.0;
	if(set->get_hist_rapidity_shift() || set->get_flag_pPb_cm_calculation()){
		rap_shift=(set->get_flag_pPb_cm_calculation())? constants::pPb_rap_shift_from_lab_to_cm:set->dlty;
	}

    for(int i=0;i<set->get_nfile();i++) {

	    if(i%1000==0) cout << ":) Reading " << i << "files." << endl; 
	    ostringstream oss;
	    string dirname = (!set->get_flag_Specify_dir_for_CentralityCut())?  set->get_dir_name(): set->get_dir_name_CentCut();
	    string fname = (!set->get_flag_Specify_f_for_CentralityCut())?  set->get_f_name() : set->get_f_name_CentCut();
	    string extname = (!set->get_flag_Specify_ext_for_CentralityCut())? set->get_ext_name() : set->get_ext_name_CentCut();
	    if(!set->zerofill) oss << dirname << "/" << fname << i << "/" << extname;
	    else oss << dirname << "/" << fname << setw(9) << setfill('0') << i << "/" << extname;
	    string inputpath=oss.str();

	    auto utl_ = make_shared<utility>();

	    EbyeInfo ebye_;
	    utl_->get_EbyeInfo_(inputpath, ebye_, rap_shift, set->get_flag_VZEROAND_trigger(), set->get_hist_parton_level());
            ebye_.orig_eventNum=i;
	    this->eBye.push_back(ebye_);

    }


}





};
