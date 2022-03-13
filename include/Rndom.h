#ifndef RNDOM
#define RNDOM
#include <random>
#include "Constants.h"

class Rndom{

private:

	int EVENT_SEED;

	public:

		//This seed shoule be initialize one per execution.
		//=================================================
		std::ranlux48 generator;
		std::ranlux48 generatorSamp;
		
		Rndom(const int seed): 
			generator(seed), generatorSamp(seed+1){
				EVENT_SEED=seed;
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


		int get_EVENT_SEED(){return this->EVENT_SEED;}



};
#endif
