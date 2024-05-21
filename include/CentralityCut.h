#ifndef CENTRALITYCUT
#define CENTRALITYCUT
#include <algorithm>
#include "EbyeInfo.h"
#include "Classification.h"
#include "ReadIn.h"

class CentralityCut{


public:

	CentralityCut(vector<EbyeInfo>& eBye_in, vector<Container::EventInfo>& nEventInfo, Options& options_in, shared_ptr<Rndom>& rndom_in):eBye(eBye_in), options(options_in), rndom(rndom_in){
		cout << ":)Start analysis for centrality cut." << endl;

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
		    this->eBye.push_back(ebye_);
		    oneEventInfo.order_reading(i); 
		    nEventInfo.push_back(oneEventInfo);
	    }
	    else// EKRT binary files are input then get vector of eBye[nev] in the following.
	    {
		    if(!options.get_flag_BSTR() || (options.get_flag_BSTR() && options.get_current_iBSTR()==0))
			    read->readEKRTbinary(nEventInfo, eBye);

		    if(options.get_flag_shuffle() || options.get_flag_BSTR()){

			    //Shuffle!
			    //=======
			    cout << ":D Shuffling the reading events!  nEventInfo.size():" << (int) nEventInfo.size() 
				    << "                                  eBye.size()      :" << (int)eBye.size() 
				    << endl;
                            //unsigned seed = 123;//Seed is intentionally fixed here
                            std::random_device rnd_device;
			    std::shuffle(std::begin(nEventInfo), std::end(nEventInfo), mt19937(rnd_device()));

			    //n_events = options.get_nfile();
			    //I want to pick up the first n events from the shuffled nEventInfo vector.
			    //Then I want to delete corresponding elements of eBye to the deleted ones.
			    //==================================================
			    int nn=0;//Can be changed.
			    std::vector<Container::EventInfo> nEventInfo_picked(nEventInfo.begin()+options.get_nfile()*nn, nEventInfo.begin()+options.get_nfile()*(nn+1));
			    vector<EbyeInfo> eBye_picked;
			    for (int k = 0; k<(int)nEventInfo_picked.size(); k++){
				    eBye_picked.push_back(eBye[nEventInfo_picked[k].order_reading()]);
				    //order_reading() is only used here for shuffling. This should be the nth events which is read by this event and archived as nth component of eBye.
			    }
			    nEventInfo=nEventInfo_picked;
			    eBye=eBye_picked;
			    std::vector<Container::EventInfo>().swap(nEventInfo_picked);
			    std::vector<EbyeInfo>().swap(eBye_picked);
			    cout << ":D After Shuffling the reading events!  nEventInfo.size():" << (int) nEventInfo.size() 
				    << "                                  eBye.size()      :" << (int)eBye.size() 
				    << endl;

		    }else if((int)nEventInfo.size()>(int)options.get_nfile()){
			    cout << "Input has " << (int)nEventInfo.size() << " events but I am going to analyse " << (int)options.get_nfile() << " events." << endl;
			    std::vector<Container::EventInfo> nEventInfo_picked(nEventInfo.begin(), nEventInfo.begin()+options.get_nfile());
			    vector<EbyeInfo> eBye_picked;
			    for (int k = 0; k<(int)nEventInfo_picked.size(); k++){
				    eBye_picked.push_back(eBye[nEventInfo_picked[k].order_reading()]);
				    //order_reading() is only used here for shuffling. This should be the nth events which is read by this event and archived as nth component of eBye.
			    }
			    nEventInfo=nEventInfo_picked;
			    eBye=eBye_picked;
			    std::vector<Container::EventInfo>().swap(nEventInfo_picked);
			    std::vector<EbyeInfo>().swap(eBye_picked);
		    }


	    }
	}


return;
}





};
#endif
