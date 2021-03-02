#include "Constants.h"

class Container{


	public:

	class ParticleInfo{

		public:

			double e;
			double px;
			double py;
			double pz;
			double pt;
			double x;
			double y;
			double z;
			double t;
			double r;
			double rap;
			double eta;
			double phi;

	};


	//Info of 1 event.
	//---------------
	class EventInfo{

		private:

			double weight_;
			int Nch_;

		public:
			EventInfo(): weight_(1.0), Nch_(0){};

			vector<ParticleInfo> part;
			void weight(double weight_in){this->weight_=weight_in;}
			void Nch(int Nch_in){this->Nch_=Nch_in;}
			double weight()const{return this->weight_;}
			int Nch()const{return this->Nch_;}

	};



       //Store information in this vector while reading files. Vector size --> # of Event
       //----------------------------------
       std::vector <EventInfo> EVENTINFO;


	public:
	double Hist[constants::x_cell_capa]={};
	double Hist_1ev[constants::x_cell_capa]={};
	double Hist_weight[constants::x_cell_capa]={};
	double HistHist[constants::x_cell_capa]={};
	double HistErr[constants::x_cell_capa]={};
	double SumWeight;
	int max_nx=-1;



        Container():SumWeight(0.0){};

};
