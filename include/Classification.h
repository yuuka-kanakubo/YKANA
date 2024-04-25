#ifndef CLASSIFICATION
#define CLASSIFICATION
#include <algorithm>
#include <functional>
#include "EbyeInfo.h"
#include "Options.h"

class Classification{

	private:
		double sumN;
                vector<EbyeInfo>& eBye;
		const Options options;

	public:
		Classification(vector<EbyeInfo>& eBye_in, const Options options_in):sumN(0.0), eBye(eBye_in), options(options_in){


		};
		~Classification(){};

       
        void V0M_classification(){
		this->get_sumN();
		this->sort_in_descending_order_of_V0M();
	};


       private:
          void get_sumN(){ 
		  int Nev=(int)eBye.size();
		  for(int i=0; i<Nev; i++){
			  if(eBye[i].valid){
				  if(options.get_flag_INEL_lg_0() && !eBye[i].multiplicity_INEL_lg_0) { continue;}
				  if(options.get_flag_2outof3_trigger() && !eBye[i].trig_2outof3) { continue;}
				  if(options.get_flag_VZEROAND_trigger() && !eBye[i].trig_VZEROAND) { continue;}
				  sumN+=eBye[i].weight;
			  }
		  }
	  }

          void sort_in_descending_order_of_V0M(){



		  //Sorting...
		  //-------------
		  vector<EbyeInfo> eBye_temp =eBye;
		  std::sort(eBye_temp.begin(), eBye_temp.end(), std::greater<EbyeInfo>());

		  //Classifying...
		  //---------------
		  double sum_w=0.0;
		  double sum_w_tot=0.0;
		  int hit=0;
                  int ct=0;
		  for(int i=0; i<(int)eBye_temp.size(); i++){

			  if(!eBye_temp[i].valid) continue;
			  if(options.get_flag_INEL_lg_0() && !eBye_temp[i].multiplicity_INEL_lg_0) continue;
			  if(options.get_flag_2outof3_trigger() && !eBye_temp[i].trig_2outof3) { continue;}
			  if(options.get_flag_VZEROAND_trigger() && !eBye[i].trig_VZEROAND) { continue;}

			  sum_w+=eBye_temp[i].weight;
			  sum_w_tot+=eBye_temp[i].weight;
			  hit++;
			  eBye[eBye_temp[i].orig_eventNum].set_V0M_class(ct);

			  if(sum_w_tot>=options.val_cent[ct]*sumN){
				  ct++;
				  sum_w=0.0;
				  hit=0;
				  if(ct>((int)options.val_cent.size())-1) break;
			  }


		  }


         }


};
#endif
