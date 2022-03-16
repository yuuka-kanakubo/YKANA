#ifndef RNDOM
#define RNDOM
#include <random>
#include "Constants.h"
#include "EbyeInfo.h"

class Rndom{

private:

	int EVENT_SEED;
	vector<int> iEv_Cent;
	bool flag_CentCut;

	public:

		//This seed shoule be initialize one per execution.
		//=================================================
		std::ranlux48 generator;
		std::ranlux48 generatorSamp;
		std::ranlux48 generatorSamp_i;

		Rndom(const int seed): 
			generator(seed), generatorSamp(seed+1), generatorSamp_i(seed+10){
				EVENT_SEED=seed;
                                flag_CentCut=false;
			};

		void keep_seed_info(const std::string path){

			std::ofstream ofs;
			ofs.open(path+"/"+constants::seed_outputfname);
			std::cout << "open: " <<path+"/"+constants::seed_outputfname << std::endl; 
			std::cout << "SEED: " << this->get_EVENT_SEED() << std::endl; 
			ofs << "# EVENT SEED" << std::endl;
			ofs << this->get_EVENT_SEED() << std::endl;
			ofs.close();

               }

		void set_flag_CentCut(bool flag){this->flag_CentCut=flag;}

		int get_EVENT_SEED(){return this->EVENT_SEED;}

		void Archive_iEv_Cent(const int iCent, const vector<EbyeInfo> eBye_CentCut){

			//Release memory
			//=================
			vector<int>().swap(this->iEv_Cent);

			for(int i=0; i<(int)eBye_CentCut.size(); i++){
				if(!eBye_CentCut[i].valid) continue;
				if(!eBye_CentCut[i].valid_assoc) continue;
				if(this->flag_CentCut){
					if(eBye_CentCut[i].get_V0M_class()==iCent) {
						this->iEv_Cent.push_back(i);
						//cout << "  ==> Pushing back " << i  << "  " << eBye_CentCut[i].N_trk_offline << "  " << eBye_CentCut[i].multiplicity << endl;
					} 
				}else{
					this->iEv_Cent.push_back(i);
				}
			}
		}

                vector<int> get_iEv_Cent(){return this->iEv_Cent;}


};
#endif
