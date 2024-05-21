#ifndef CENTRALITYCUT
#define CENTRALITYCUT
#include <algorithm>
#include "EbyeInfo.h"
#include "Classification.h"
#include "ReadIn.h"

class CentralityCut{


public:

	CentralityCut(vector<EbyeInfo>& eBye_in, vector<Container::EventInfo>& nEventInfo, Options& options_in, shared_ptr<Rndom>& rndom_in):options(options_in), rndom(rndom_in){
		cout << ":)Start analysis for centrality cut." << endl;

			this->read_events(nEventInfo, eBye_in);

	}

         void ClassifyCentrality(vector<EbyeInfo>& eBye){

			Classification* ev_cl=new Classification(eBye, options);
			if(options.get_collision_type()==101){
				cout << "This is Ntrk classification. No need to sort." << endl;
			}else{
				ev_cl->V0M_classification();
			}
			delete ev_cl;
			cout << ":)Finish centrality classification." << endl;

	 }



	std::vector<Container::EventInfo> get_nEventInfo(){return this->nEventInfo_for_BSTR;};
	std::vector<EbyeInfo> get_eBye_CentCut(){return this->eBye_CentCut_for_BSTR;};


private:
	Options& options;
	shared_ptr<Rndom>& rndom;
	shared_ptr<ReadIn> read;
	shared_ptr<Message> ms;
	std::vector<Container::EventInfo> nEventInfo_for_BSTR;
	std::vector<EbyeInfo> eBye_CentCut_for_BSTR;



void read_events(vector<Container::EventInfo>& nEventInfo, vector<EbyeInfo>& eBye){

	ms = make_shared<Message>();
	read = make_shared<ReadIn>(this->ms, this->options);
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
		    eBye.push_back(ebye_);
		    oneEventInfo.order_reading(i); 
		    nEventInfo.push_back(oneEventInfo);
	    }
	    else// EKRT binary files are input then get vector of eBye[nev] in the following.
	    {
		    if(!options.get_flag_BSTR() || (options.get_flag_BSTR() && options.get_current_iBSTR()==0))
			    read->readEKRTbinary(nEventInfo, eBye);


	    }


	}//event loop

	if(options.get_flag_BSTR() && options.get_current_iBSTR()==0){
		this->nEventInfo_for_BSTR=nEventInfo;
		this->eBye_CentCut_for_BSTR=eBye;
	}


return;
}





};
#endif
