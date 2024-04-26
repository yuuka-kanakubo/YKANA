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
#include "Util_func.h"
#include "Container.h"
#include "Message.h"
#include "LogSettings.h"
#include "Settings.h"
#include "CentralityCut.h"
#include "ReadIn.h"

using namespace std;

ReadIn::ReadIn(shared_ptr<Message> ms_in, Options options_in):ms(ms_in), options(options_in), ncall_readTimeLapse(0), nline(-1){};
ReadIn::~ReadIn(){};




bool ReadIn::readXY(const std::string& fname, shared_ptr<Container>& ct){

			ifstream in;
			in.open(fname.c_str(),ios::in);
			if(!in){ ms->open(fname); return false;}

			//Variables.
			//part_1ev is a vector containing particle lists
			//info_1ev stores info of one single event.
			//-----------
			Container::EventInfo info_1ev;
			vector<Container::ParticleInfo> part_1ev;

			{
				std::string templine;
				while(getline(in,templine)) {
					if(templine.find('#')!=std::string::npos) {
					} else if(templine.find('%')!=std::string::npos){
					}else{
						istringstream is(templine);
						double data1, data2;//I dare to leave them as arbitral values so that this function is available for anything.
						is >> data1 >> data2;


						Container::ParticleInfo part_in;
						part_in.r=data1;
						part_in.vt=data2;
						part_1ev.push_back(part_in);

					}
				}

			}

			info_1ev.part=part_1ev;
			ct->EVENTINFO=info_1ev;
			vector<Container::ParticleInfo>().swap(part_1ev);

			return true;
}



bool ReadIn::readEKRT(const std::string& fname, shared_ptr<Container>& ct){

	ifstream in;
	in.open(fname.c_str(),ios::in);
	if(!in){ ms->open(fname); return false;}

	//Variables.
	//part_1ev is a vector containing particle lists
	//info_1ev stores info of one single event.
	//-----------
	Container::EventInfo info_1ev;
	vector<Container::ParticleInfo> part_1ev;
	int Nch=0;
	double weight=1.0;


	{
		std::string templine;
		while(getline(in,templine)) {
			if(templine.find('#')!=std::string::npos) {
			} else if(templine.find('%')!=std::string::npos){
			}else{
				double x, y, pt, rap, phi, t;
				double tata;
				int ID;
				int momID1;
				int momID2;
				int imomNucleon1;
				int imomNucleon2;
				double xmomNucleon1;
				double xmomNucleon2;
				double ymomNucleon1;
				double ymomNucleon2;
				double zmomNucleon1;
				double zmomNucleon2;
				bool is_mom1Neutron;
				bool is_mom2Neutron;
				istringstream is(templine);
				is >> t 
					>> x 
					>> y 
					>> pt 
					>> rap 
					>> phi 
					>> tata 
					>> ID 
					>> momID1 
					>> momID2 
					>> imomNucleon1 
					>> imomNucleon2 
					>> xmomNucleon1 
					>> xmomNucleon2 
					>> ymomNucleon1 
					>> ymomNucleon2 
					>> zmomNucleon1 
					>> zmomNucleon2 
					>> is_mom1Neutron 
					>> is_mom2Neutron;

				Container::ParticleInfo part_in;

				part_in.x=x;//[fm]
				part_in.y=y;//[fm]
				double rr = pow(x,2)+pow(y,2);
				part_in.r=(rr<constants::SMALL)? 0.0:sqrt(rr);//[fm]
				part_in.t=t;//[fm]
				part_in.phi=phi;
				part_in.rap=rap;
				part_in.pt=pt;//[GeV]
				part_in.tata=tata;
				part_in.ID=ID;
				part_in.momID1=momID1;
				part_in.momID2=momID2;
				part_in.imomNucleon1=imomNucleon1;
				part_in.imomNucleon2=imomNucleon2;
				part_in.xmomNucleon1=xmomNucleon1;
				part_in.xmomNucleon2=xmomNucleon2;
				part_in.ymomNucleon1=ymomNucleon1;
				part_in.ymomNucleon2=ymomNucleon2;
				part_in.zmomNucleon1=zmomNucleon1;
				part_in.zmomNucleon2=zmomNucleon2;
				part_in.is_mom1Neutron=is_mom1Neutron;
				part_in.is_mom2Neutron=is_mom2Neutron;

				double mtsq=pow(part_in.pt,2)+pow(part_in.m,2);
				double mt=(mtsq>0.0) ? sqrt(mtsq):0.0;
				part_in.mt=mt;
				part_in.px=part_in.pt*cos(part_in.phi);
				part_in.py=part_in.pt*sin(part_in.phi);
				part_in.pz=part_in.mt*sinh(part_in.rap);
				part_in.e=part_in.mt*cosh(part_in.rap);
				part_in.tau=part_in.t/cosh(part_in.rap);
				part_in.z=part_in.tau*sinh(part_in.rap);
				double msq=pow(part_in.e,2)
					-pow(part_in.px,2)
					-pow(part_in.py,2)
					-pow(part_in.pz,2);
				double m=(msq>0.)? sqrt(msq):0.;
				part_in.m=m;
				//if(std::fabs(part_in.eta)<constants::delta_eta ) { 
					//if(std::fabs(part_in.x)<constants::delta_xcoord) { 
						part_1ev.push_back(part_in);
					//}
				//}

			}//particle loop if no # %
		}//particle loop
		in.close();
	}//namespace 

	//Store
	//part_1ev is a vector containing particle lists
	//info_1ev stores info of one single event.
	//-------
	info_1ev.weight(weight);
	info_1ev.Nch(Nch);
	info_1ev.part=part_1ev;
	ct->EVENTINFO=info_1ev;
	vector<Container::ParticleInfo>().swap(part_1ev);

	return true;


}



bool ReadIn::readEKRTbinary(std::vector <Container::EventInfo>& nEventInfo){

	std::stringstream ss;
	ss << options.get_dir_name() << "/" << options.get_ext_name();
	std::string inputf = ss.str();
	std::ifstream in(inputf, std::ios::in | std::ios::binary);
	if(!in.is_open()) {
		cout << "ERROR:( " << __FILE__ << " (" << __LINE__ << " )unable to open file " << inputf << endl;
		exit(EXIT_FAILURE);
	} 
	uint64_t n_events, n_jets;
	double t01, t02, x, y, pt, y1, y2, phi, tata;
	int_fast16_t init1, init2, final1, final2;
	uint_fast16_t ia, ib;
	double xa, ya, za, xb, yb, zb;
	bool a_is_neutron, b_is_neutron;
	//double TotPt=0.0;
	//double TotE=0.0;
	//int nEvent=0;
	if (sizeof n_jets != 8 || sizeof t01 != 8 )
	{
		std::cout << "Change the types! n_jets should be 64 bit unsigned integer"
			<< " and the coordinates should all be 64 bit float numbers" << std::endl;
		exit(EXIT_FAILURE);
	}

	in.read(reinterpret_cast<char*>(&n_events), sizeof n_events);
	uint64_t n_jet_total = 0;
	int pct = 0;
	n_events = options.get_nfile();
	std::cout << "Trying to read " << n_events << " events from the file "
              << inputf << " ..." << std::endl;
	for (uint64_t i=0; i<n_events; i++){
		if(fabs((double)i/(double)n_events-pct*0.1)<constants::SMALL){
			cout << ":) " << pct*10 << "\% is done " << endl;
			pct++;
		}

		in.read(reinterpret_cast<char*>(&n_jets), sizeof n_jets);
		//cout << " n_jets " << n_jets << endl;
		//if(i>=1) break;
		//cout << i << "  in  " << n_events << " : " << n_jets<< endl;

		//part_1ev is a vector containing particle lists
		//-----------
		Container::EventInfo oneEventInfo;
		vector<Container::ParticleInfo> part_1ev;
		//This should be each minijet loop
		for (uint64_t ii=0; ii<n_jets; ii++, n_jet_total++){


			in.read(reinterpret_cast<char*>(&t01)  , sizeof t01);  //t01
			in.read(reinterpret_cast<char*>(&t02)  , sizeof t02);  //t02
			in.read(reinterpret_cast<char*>(&x)    , sizeof x);    //x
			in.read(reinterpret_cast<char*>(&y)    , sizeof y);    //y
			in.read(reinterpret_cast<char*>(&pt)   , sizeof pt);   //pT
			in.read(reinterpret_cast<char*>(&y1)   , sizeof y1);   //y1
			in.read(reinterpret_cast<char*>(&y2)   , sizeof y2);   //y2
			in.read(reinterpret_cast<char*>(&phi)  , sizeof phi);  //phi_1
			in.read(reinterpret_cast<char*>(&tata) , sizeof tata); //T_A * T_A
			in.read(reinterpret_cast<char *>(&init1), sizeof init1);   // flavour of incoming 1
			in.read(reinterpret_cast<char *>(&init2), sizeof init2);   // flavour of incoming 2
			in.read(reinterpret_cast<char *>(&final1), sizeof final1); // flavour of outgoing 1
			in.read(reinterpret_cast<char *>(&final2), sizeof final2); // flavour of outgoing 2

			in.read(reinterpret_cast<char *>(&ia), sizeof ia); // index of mother a
			in.read(reinterpret_cast<char *>(&ib), sizeof ib); // index of mother b
			in.read(reinterpret_cast<char *>(&xa), sizeof xa);
			in.read(reinterpret_cast<char *>(&ya), sizeof ya);
			in.read(reinterpret_cast<char *>(&za), sizeof za);
			in.read(reinterpret_cast<char *>(&xb), sizeof xb);
			in.read(reinterpret_cast<char *>(&yb), sizeof yb);
			in.read(reinterpret_cast<char *>(&zb), sizeof zb);
			in.read(reinterpret_cast<char *>(&a_is_neutron), sizeof a_is_neutron);
			in.read(reinterpret_cast<char *>(&b_is_neutron), sizeof b_is_neutron);

			//I will put the information into container
			//=========================================
			//Minijet 1.
			Container::ParticleInfo parton_in1;
			parton_in1.x=x;//[fm]
			parton_in1.y=y;//[fm]
			parton_in1.t=t01*constants::hbarc;//[fm]
			parton_in1.phi=phi;
			parton_in1.rap=y1;
			parton_in1.pt=pt;//[GeV]
			parton_in1.tata=tata;
			parton_in1.ID=final1;
			parton_in1.momID1=init1;
			parton_in1.momID2=init2;
			parton_in1.imomNucleon1=ia;
			parton_in1.imomNucleon2=ib;
			parton_in1.xmomNucleon1=xa;
			parton_in1.xmomNucleon2=xb;
			parton_in1.ymomNucleon1=ya;
			parton_in1.ymomNucleon2=yb;
			parton_in1.zmomNucleon1=za;
			parton_in1.zmomNucleon2=zb;
			parton_in1.is_mom1Neutron=a_is_neutron;
			parton_in1.is_mom2Neutron=b_is_neutron;

			double mtsq=pow(parton_in1.pt,2)+pow(0.0,2);
			double mt=(mtsq>0.0) ? sqrt(mtsq):0.0;
			parton_in1.mt=mt;
			parton_in1.px=parton_in1.pt*cos(parton_in1.phi);
			parton_in1.py=parton_in1.pt*sin(parton_in1.phi);
			parton_in1.pz=parton_in1.mt*sinh(parton_in1.rap);
			parton_in1.e=parton_in1.mt*cosh(parton_in1.rap);
			parton_in1.tau=parton_in1.t/cosh(parton_in1.rap);
			parton_in1.z=parton_in1.tau*sinh(parton_in1.rap);
			double msq=pow(parton_in1.e,2)
				-pow(parton_in1.px,2)
				-pow(parton_in1.py,2)
				-pow(parton_in1.pz,2);
			double m=(msq>0.)? sqrt(msq):0.;
			if(m>constants::MEDSMALL){
				cout << __FILE__ << " Should be mass-less in EKRT. m: " << m << endl;
				exit(1);
			}
			parton_in1.m=m;

			part_1ev.push_back(parton_in1);

			//Minijet 2.
			Container::ParticleInfo parton_in2;
			parton_in2.x=x;//[fm]
			parton_in2.y=y;//[fm]
			parton_in2.t=t02*constants::hbarc;//[fm]
			parton_in2.phi=phi+M_PI;
			parton_in2.rap=y2;
			parton_in2.pt=pt;//[GeV]
			parton_in2.tata=tata;
			parton_in2.ID=final2;
			parton_in2.momID1=init1;
			parton_in2.momID2=init2;
			parton_in2.imomNucleon1=ia;
			parton_in2.imomNucleon2=ib;
			parton_in2.xmomNucleon1=xa;
			parton_in2.xmomNucleon2=xb;
			parton_in2.ymomNucleon1=ya;
			parton_in2.ymomNucleon2=yb;
			parton_in2.zmomNucleon1=za;
			parton_in2.zmomNucleon2=zb;
			parton_in2.is_mom1Neutron=a_is_neutron;
			parton_in2.is_mom2Neutron=b_is_neutron;

			mtsq=pow(parton_in2.pt,2)+pow(0.0,2);
			mt=(mtsq>0.0) ? sqrt(mtsq):0.0;
			parton_in2.mt=mt;
			parton_in2.px=parton_in2.pt*cos(parton_in2.phi);
			parton_in2.py=parton_in2.pt*sin(parton_in2.phi);
			parton_in2.pz=parton_in2.mt*sinh(parton_in2.rap);
			parton_in2.e=parton_in2.mt*cosh(parton_in2.rap);
			parton_in2.tau=parton_in2.t/cosh(parton_in2.rap);
			parton_in2.z=parton_in2.tau*sinh(parton_in2.rap);
			msq=pow(parton_in2.e,2)
				-pow(parton_in2.px,2)
				-pow(parton_in2.py,2)
				-pow(parton_in2.pz,2);
			m=(msq>0.)? sqrt(msq):0.;
			if(m>constants::MEDSMALL){
				cout << __FILE__ << " Should be mass-less in EKRT. m: " << m << endl;
				exit(1);
			}
			parton_in2.m=m;

			part_1ev.push_back(parton_in2);


		}//minijet loop

		//cout << "Finished minijet loop :D " << endl;


		//Once found a last minijet from one single event,
		//put info of one event into 'oneEventInfo'
		//====================================================
		oneEventInfo.part = part_1ev;
		//oneEventInfo. = Nch 
		//oneEventInfo. = weight;
		//oneEventInfo. = ....;
		vector<Container::ParticleInfo>().swap(part_1ev);

		nEventInfo.push_back(oneEventInfo);

	}//Event loop

	in.close();

	//	//Put info of Count Ev
	//	//ct->CountEv++;
	return true;
}




bool ReadIn::read(const std::string& fname, shared_ptr<Container>& ct){

			ifstream in;
			in.open(fname.c_str(),ios::in);
			if(!in){ ms->open(fname); return false;}

			//Variables.
			//part_1ev is a vector containing particle lists
			//info_1ev stores info of one single event.
			//-----------
			Container::EventInfo info_1ev;
			vector<Container::ParticleInfo> part_1ev;
			int Nch=0;
			double weight=1.0;


			{
				std::string templine;
				while(getline(in,templine)) {
					if(templine.find('#')!=std::string::npos) {
					} else if(templine.find('%')!=std::string::npos){
						istringstream iss(templine);
						std::string pct;
						int iev, nv;
						string com;
						double weight_in, tau;
						iss >> com >> iev >> nv >> tau >> weight_in;
						weight=weight_in;
					}else{
						istringstream is(templine);
						int data1, data2, ID, col, acol;
						double m,e,px,py,pz,x,y,z,t,ft, rap;
						std::string TAG;
						is >> data1 >> data2 >> col >> acol >> ID >> m >> e >> px >> py >> pz >> rap >> x >> y >> z >> t >> ft >> TAG;


						double P_squared=px*px+py*py+pz*pz;
						double P=(P_squared)>0.0 ? sqrt(P_squared):0.0;
						double pt_squared=px*px+py*py;
						double pt=(pt_squared)>0.0 ? sqrt(pt_squared):0.0;
						double mt_squared=pt*pt+m*m;
						double mt=(mt_squared)>0.0 ? sqrt(mt_squared):0.0;
						double phi = atan2(py, px);

						////
						double eta;
						{
							const double LARGE = 10.0;
							if(P==0.0 && pz==0.0) eta=0.0;
							else if(P==-pz) eta=-LARGE;
							else if(P==pz) eta=LARGE;
							else eta=std::log((P + pz)/(P - pz)) / 2.0;
						}
						///



						if(constants::MODE.find("dndeta_proton")!=string::npos){
							if(ID==constants::id_proton) { 

								Container::ParticleInfo part_in;
								part_in.eta=eta;
								part_1ev.push_back(part_in);

							}
						}else if(constants::MODE.find("dndpt")!=string::npos || constants::MODE.find("vertices")!=string::npos){
							if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon) && std::fabs(eta)<constants::delta_eta ) { 

								if(constants::MODE.find("dndpt")!=string::npos) {
									Container::ParticleInfo part_in;
									part_in.pt=pt;
									part_1ev.push_back(part_in);
								};
								double r_squared = x*x + y*y;
								double r=(r_squared>0.0)? sqrt(r_squared):0.0;
								if(constants::MODE.find("vertices")!=string::npos){
									Container::ParticleInfo part_in;
									part_in.r=r;
									part_1ev.push_back(part_in);
								}

							}

						}else if(constants::MODE.find("mtscaling")!=string::npos){
							if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon||abs(ID)==constants::id_phi||abs(ID)==constants::id_lambda ||abs(ID)==constants::id_cascade || abs(ID)==constants::id_omega ) && std::fabs(eta)<0.5 ) { 

								if((options.get_flag_only_corona() && TAG == constants::corona_tag) || (options.get_flag_only_core() && TAG == constants::core_tag) || (!options.get_flag_only_core() && !options.get_flag_only_corona() )){



									Container::ParticleInfo part_in;
									part_in.id=ID;
									part_in.mt=mt;
									part_in.e=e;
									part_in.px=px;
									part_in.py=py;
									part_in.pz=pz;
									part_in.m=m;
									part_1ev.push_back(part_in);
								}

							}
							if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon) && std::fabs(eta)<0.3 && pt<10.0 ) Nch++;


						}else if(constants::MODE.find("meanmt")!=string::npos){

							cout << "ERROR:( -DMEANMT is still under construction." << endl;
							exit(1);
							if(abs(ID)==constants::id_cascade && std::fabs(eta)<0.3 && pt>0.15 && pt<10.0 ) { 

								Container::ParticleInfo part_in;
								part_in.mt=mt;
								part_in.id=ID;
								part_1ev.push_back(part_in);

							}
							if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon) && std::fabs(eta)<0.3 && pt<10.0 ) Nch++;


						}else if(constants::MODE.find("meanpt")!=string::npos){
							if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon) && std::fabs(eta)<0.3 && pt>0.15 && pt<10.0 ) { 
								if((options.get_flag_only_corona() && TAG == constants::corona_tag) || (options.get_flag_only_core() && TAG == constants::core_tag)|| (!options.get_flag_only_core() && !options.get_flag_only_corona() )){

								Container::ParticleInfo part_in;
								part_in.pt=pt;
								part_1ev.push_back(part_in);

								}

							}
							if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon) && std::fabs(eta)<0.3 ) Nch++;


						}else if(constants::MODE.find("MeanptPID")!=string::npos){
							if(abs(ID)==options.get_specify_ID() && std::fabs(eta)<constants::default_midy_pm) { 
								if((options.get_flag_only_corona() && TAG == constants::corona_tag) || (options.get_flag_only_core() && TAG == constants::core_tag)|| (!options.get_flag_only_core() && !options.get_flag_only_corona() )){

								Container::ParticleInfo part_in;
								part_in.pt=pt;
								part_1ev.push_back(part_in);

								}

							}
							if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon) && std::fabs(eta)<0.3 ) Nch++;


						}else if(constants::MODE.find("Rt_spectra")!=string::npos){

							if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon ) && fabs(eta)<constants::etaRange_Rt && pt>0.15){

								if((options.get_flag_only_corona() && TAG == constants::corona_tag) || (options.get_flag_only_core() && TAG == constants::core_tag)|| (!options.get_flag_only_core() && !options.get_flag_only_corona() )){
									Container::ParticleInfo part_in;
									part_in.pt=pt;
									part_in.eta=eta;
									part_in.phi=phi;
									part_in.TAG=TAG;
									part_1ev.push_back(part_in);
								}
							}
							if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon) && fabs(eta)<constants::etaRangeMultiplicity) Nch++;


						}else if(constants::MODE.find("Rt_yield")!=string::npos){

							if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon ||abs(ID)==constants::id_phi||abs(ID)==constants::id_lambda ||abs(ID)==constants::id_cascade || abs(ID)==constants::id_omega) && fabs(eta)<constants::etaRange_Rt && pt>0.15){
								if((options.get_flag_only_corona() && TAG == constants::corona_tag) || (options.get_flag_only_core() && TAG == constants::core_tag)|| (!options.get_flag_only_core() && !options.get_flag_only_corona() )){
									Container::ParticleInfo part_in;
									part_in.id=ID;
									part_in.px=px;
									part_in.py=py;
									part_in.pt=pt;
									part_in.eta=eta;
									part_in.phi=phi;
									part_in.TAG=TAG;
									part_1ev.push_back(part_in);
								}
							}
							if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon) && fabs(eta)<constants::etaRangeMultiplicity) Nch++;


						}else if(constants::MODE.find("cumulant_multi")!=string::npos){

							if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon ) && fabs(eta)<constants::etaRange_cumulantmulti && constants::ptmin_cumulantmulti < pt && pt<constants::ptmax_cumulantmulti ){
								if((options.get_flag_only_corona() && TAG == constants::corona_tag) || (options.get_flag_only_core() && TAG == constants::core_tag)|| (!options.get_flag_only_core() && !options.get_flag_only_corona() )){
									Container::ParticleInfo part_in;
									part_in.pt=pt;
									part_in.eta=eta;
									part_in.phi=phi;
									part_1ev.push_back(part_in);
								}
							}
							if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon) && fabs(eta)<constants::etaRange_cumulantmulti_Nch && constants::ptmin_cumulantmulti_Nch < pt && constants::ptmax_cumulantmulti_Nch > pt) Nch++;


						}else if(constants::MODE.find("cumulant_pt")!=string::npos){

							if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon ) && constants::ptmin_cumulantpt<pt && pt < constants::ptmax_cumulantpt && fabs(eta)<constants::etaRange_cumulantpt){
								if((options.get_flag_only_corona() && TAG == constants::corona_tag) || (options.get_flag_only_core() && TAG == constants::core_tag)|| (!options.get_flag_only_core() && !options.get_flag_only_corona() )){
									Container::ParticleInfo part_in;
									part_in.pt=pt;
									part_in.eta=eta;
									part_in.phi=phi;
									part_1ev.push_back(part_in);
								}
							}


						}else if(constants::MODE.find("cumulant_eta")!=string::npos){

							if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon ) && constants::ptmin_cumulanteta < pt && pt<constants::ptmax_cumulanteta ){
								if((options.get_flag_only_corona() && TAG == constants::corona_tag) || (options.get_flag_only_core() && TAG == constants::core_tag)|| (!options.get_flag_only_core() && !options.get_flag_only_corona() )){
									Container::ParticleInfo part_in;
									part_in.pt=pt;
									part_in.eta=eta;
									part_in.phi=phi;
									part_1ev.push_back(part_in);
								}
							}
							if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon) && fabs(eta)<1.0 && pt<3.0 ) Nch++;


						}else if(constants::MODE.find("twopc2D")!=string::npos){

							if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon ) && fabs(eta)<constants::twopc2DetaRange){
								if((options.get_flag_only_corona() && TAG == constants::corona_tag) || (options.get_flag_only_core() && TAG == constants::core_tag)|| (!options.get_flag_only_core() && !options.get_flag_only_corona() )){
									Container::ParticleInfo part_in;
									part_in.id=ID;
									part_in.pt=pt;
									part_in.eta=eta;
									part_in.phi=phi;
									part_in.TAG=TAG;
									part_1ev.push_back(part_in);
								}
							}
							if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon) && fabs(eta)<constants::twopc2DetaRange_N && pt>constants::twopc2Dptmin_N && pt<constants::twopc2Dptmax_N ) Nch++;


						}else if(constants::MODE.find("twopc1D")!=string::npos || constants::MODE.find("twopcInteg")!=string::npos){

							if(((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon ) && fabs(eta)< constants::twopc1DetaRange && pt > constants::twopc1Dptmin && pt<constants::twopc1Dptmax)){

								if((options.get_flag_only_corona() && TAG == constants::corona_tag) || (options.get_flag_only_core() && TAG == constants::core_tag)|| (!options.get_flag_only_core() && !options.get_flag_only_corona() )){
									Container::ParticleInfo part_in;
									part_in.pt=pt;
									part_in.id=ID;
									part_in.eta=eta;
									part_in.phi=phi;
									part_in.TAG=TAG;
									part_1ev.push_back(part_in);
								}
							}
							if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon) && fabs(eta)<constants::twopc1DetaRange_N && pt>constants::twopc1Dptmin_N && pt<constants::twopc1Dptmax_N ) Nch++;


						}else{
							if((abs(ID)==constants::id_proton||abs(ID)==constants::id_ch_pion||abs(ID)==constants::id_ch_kaon) ) { 

								Container::ParticleInfo part_in;
								part_in.eta=eta;
								part_1ev.push_back(part_in);
							}
						}



						}
					}
					in.close();
				}

				//Store
				//part_1ev is a vector containing particle lists
				//info_1ev stores info of one single event.
				//-------
				info_1ev.weight(weight);
				info_1ev.Nch(Nch);
				info_1ev.part=part_1ev;
				ct->EVENTINFO=info_1ev;
				vector<Container::ParticleInfo>().swap(part_1ev);

				return true;
			}





			
				bool ReadIn::read_jetinfo(const std::string& fname, shared_ptr<Container>& ct){

					ifstream in;
					in.open(fname.c_str(),ios::in);
					if(!in){ ms->open(fname); return false;}
			
					Container::EventInfo info_1ev;
					vector <Container::ParticleInfo> Jets, Gamma;
					{
						std::string templine;
						while(getline(in,templine)) {
							if(templine.find('G')!=std::string::npos) {
								istringstream is(templine);
								double e,px,py,pz;
								int num;
								is >> num >> e >> px >> py >> pz;
								Container::ParticleInfo info;
								info.e=e;
								info.px=px;
								info.py=py;
								info.pz=pz;
								info.pt=sqrt(px*px + py*py);
								info.phi=atan2(py, px);
								Gamma.push_back(info);
							}else if(templine.find('P')!=std::string::npos){
							}else if(templine.find('#')!=std::string::npos){
							}else{
								istringstream is(templine);
								double e,px,py,pz;
								int num;
								is >> num >> e >> px >> py >> pz;
								Container::ParticleInfo info;
								info.e=e;
								info.px=px;
								info.py=py;
								info.pz=pz;
								info.pt=sqrt(px*px + py*py);
								info.phi=atan2(py, px);
								Jets.push_back(info);
							}
						}
						in.close();
					}
			
			
					if((int)Jets.size()<2){ 
						vector<Container::ParticleInfo>().swap(Jets);
						vector<Container::ParticleInfo>().swap(Gamma);
						return false;
					}
					for(int i =0; i<1; i++){
						if((Jets[i].pt+Jets[i+1].pt)<=0.0) {
							vector<Container::ParticleInfo>().swap(Jets);
							vector<Container::ParticleInfo>().swap(Gamma);
							return false;
						}
						info_1ev.Aj(fabs((Jets[i].pt-Jets[i+1].pt)/(Jets[i].pt+Jets[i+1].pt)));
					}
			
					ct->EVENTINFO=info_1ev;

					vector<Container::ParticleInfo>().swap(Jets);
					vector<Container::ParticleInfo>().swap(Gamma);
					return true;
				}



bool ReadIn::readTimeLapse(const std::string& fname, shared_ptr<Container>& ct, const double weight){


	//- at first event, get line num of the point to see.
	if (ncall_readTimeLapse==0)
		if(!this->get_nline_to_see(nline, fname)) return false;
	if(nline<0){cout << ":(ERROR " << __FILE__ << __LINE__ << endl; exit(1);}


	//Variables.
	//step_1ev is a vector containing time step info lists
	//info_1ev stores info of one single event.
	//-----------
	Container::EventInfo info_1ev;
	vector<Container::StepInfo> step_1ev;

	//TimeStep
	//---------
	double tau=constants::TL_tau_00;
	for(int step=0; step < constants::LARGEint ; step++){
		double dtau =(tau<= constants::TL_tau_switch)? constants::TL_dtau1:constants::TL_dtau2;
		if(tau*constants::hbarc>constants::TL_tau_switchFM-constants::TL_dtau1FM && tau*constants::hbarc<=constants::TL_tau_switchFM)dtau=constants::TL_dtau2;

		//The following line is for the DCCI with bug in File_fluid_manager setting tau00 as 0.01.
		//--------------------------------
		if(step==0) tau=0.0;
		if(step==1) tau=constants::TL_tau_00+dtau;
		ostringstream os;
		os << fname << options.get_ext_nameTL() << fixed << setprecision(2) << tau*constants::hbarc << ".txt";
		ifstream in;
		in.open(os.str().c_str(),ios::in);
		if(!in){
					Container::StepInfo step_oneline;
					//Set time
					//-----------
					stringstream tmp;
					double tau_=tau*constants::hbarc;
					tmp << setprecision(2) << fixed << tau_;
					double tau__ = stod(tmp.str());
					step_oneline.tau=tau__;
					step_oneline.nstep=step;
					step_1ev.push_back(step_oneline);
		}

		std::string templine;
		int nline_current=0;
		while(getline(in,templine)) {
			if(templine.find('#')!=std::string::npos) {
			} else if(templine.find('%')!=std::string::npos){
			} else if(templine.empty()){
			}else{
				if(nline_current==nline){
					istringstream is(templine);
					Container::StepInfo step_oneline;
					if(options.get_modeTL()==0) this->get_oneline_xy(is, step_oneline);
					else if(options.get_modeTL()==1)this->get_oneline_xeta(is, step_oneline);
					else{cout << ":( ERROR Something wrong. " << __FILE__ << __LINE__ << endl; exit(1);}
					step_oneline.nstep=step;
					step_1ev.push_back(step_oneline);
				}
			}
			nline_current++;
		}//while

		in.close();
		//Step
		//----
		tau+=dtau;

		if(tau>constants::x_max) {
			break;
		}
		if(tau*constants::hbarc>100.0) {
			cout << ":oWARNING It seems that something is wrong in hydro profile. The time step loop continues to 100fm." << endl;
			break;
		}

	}//step


	//Store
	//-------
	info_1ev.weight(weight);
	info_1ev.step=step_1ev;
	ct->EVENTINFO=info_1ev;
	vector<Container::StepInfo>().swap(step_1ev);

	ncall_readTimeLapse++;
	return true;

}

bool ReadIn::get_nline_to_see(int &nline, const std::string fname){


	double tau=constants::TL_tau_00;
	tau+=constants::TL_dtau1;
	ostringstream os;
	os << fname << options.get_ext_nameTL() << fixed << setprecision(2) << tau*constants::hbarc << ".txt";
	ifstream in;
	in.open(os.str().c_str(),ios::in);
	if(!in){cout << os.str() << endl;return false;}
	std::string templine;
	double delta=constants::LARGE;
	double xTL = options.get_at_xTL();
	double yTL = options.get_at_yTL();
	double etaTL = options.get_at_etaTL();
	double x_current, y_current, eta_current;
	int nline_=0;
	while(getline(in,templine)) {
		if(templine.find('#')!=std::string::npos) {
		} else if(templine.find('%')!=std::string::npos){
		} else if(templine.empty()){
		}else{
			istringstream is(templine);
			Container::StepInfo step_oneline;
			if(options.get_modeTL()==0)this->get_oneline_xy(is, step_oneline);
			else if(options.get_modeTL()==1)this->get_oneline_xeta(is, step_oneline);
			else{cout << ":( ERROR Something wrong. " << __FILE__ << __LINE__ << endl; exit(1);}
			if(options.get_modeTL()==0){
				if(delta>(pow(fabs(xTL-step_oneline.x),2)+pow(fabs(yTL-step_oneline.y),2))){
					x_current=step_oneline.x;
					y_current=step_oneline.y;
					nline=nline_;
					delta=(pow(fabs(xTL-step_oneline.x),2)+pow(fabs(yTL-step_oneline.y),2));
				}
				//if(delta_y>fabs(yTL-step_oneline.y)) y_current=step_online.y;
			}else if(options.get_modeTL()==1){
				if(delta>(pow(fabs(xTL-step_oneline.x),2)+pow(fabs(etaTL-step_oneline.eta),2))){
					x_current=step_oneline.x;
					eta_current=step_oneline.eta;
					nline=nline_;
					delta=(pow(fabs(xTL-step_oneline.x),2)+pow(fabs(etaTL-step_oneline.eta),2));
				}
			}

		}
		nline_++;
	}//while

	cout << "Printing nearest points: x " << x_current << "  y " << y_current << "  eta " << eta_current << endl;
	cout << "                         ==>  " << nline << endl;
	cout << "                 inputs: x " << options.get_at_xTL() << "  y " << options.get_at_yTL() << "  eta " << options.get_at_etaTL()<< endl;

	in.close();
return true;

}


void ReadIn::get_oneline_xeta(istringstream& is, Container::StepInfo& onestep){
			is >> onestep.tau //[fm]
				>> onestep.x //fm
				>> onestep.eta //[1]
				>> onestep.e //[GeV^4]
				>> onestep.temp //[GeV]
				>> onestep.s //[GeV^3]
				>> onestep.p //[GeV]
				>> onestep.n //[GeV^3]
				>> onestep.n5 
				>> onestep.mu 
				>> onestep.mu5
				>> onestep.vx  //[1]
				>> onestep.veta  //[1]
				>> onestep.Ueta
				>> onestep.Vtilde
				>> onestep.UABS
				>> onestep.Ux
				>> onestep.Uy
				>> onestep.Ueta
				>> onestep.eEtildedx
				>> onestep.eEtildedy
				>> onestep.eBtildedx
				>> onestep.eBtildedy
				>> onestep.eEtilde_dot_eBtilde
				>> onestep.U4;
}

void ReadIn::get_oneline_xy(istringstream& is, Container::StepInfo& onestep){
			is >> onestep.tau //1
				>> onestep.x //2
				>> onestep.y //3
				>> onestep.e //4
				>> onestep.temp //5
				>> onestep.s //6
				>> onestep.p //7
				>> onestep.n //8
				>> onestep.n5 //9
				>> onestep.mu //10
				>> onestep.mu5 //11
				>> onestep.Vtilde //12
				>> onestep.vx //13
				>> onestep.vy //14
				>> onestep.veta //15	   
				>> onestep.U4 //16
				>> onestep.U_R4 //17
				>> onestep.U_L4 //18	
				>> onestep.eEtildedx  //19
				>> onestep.eEtildedy  //20
				>> onestep.eBtildedx  //21
				>> onestep.eBtildedy  //22
				>> onestep.eEtilde_dot_eBtilde  //23
				>> onestep.Ux //24
				>> onestep.Uy //25	
				>> onestep.Ueta; //26	

}
