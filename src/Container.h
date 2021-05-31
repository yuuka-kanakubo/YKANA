#ifndef CONTAINER
#define CONTAINER
#include "Constants.h"

class Container{


	public:

	class ParticleInfo{

		public:

			int id;
			double e;
			double m;
			double mt;
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
			double Aj_;

		public:
			EventInfo(): weight_(1.0), Nch_(0), Aj_(1.0){};
			~EventInfo(){
				vector<ParticleInfo>().swap(part);
			};

			vector<ParticleInfo> part;
			void weight(double weight_in){this->weight_=weight_in;}
			void Nch(int Nch_in){this->Nch_=Nch_in;}
			void Aj(double Aj_in){this->Aj_=Aj_in;}
			double weight()const{return this->weight_;}
			int Nch()const{return this->Nch_;}
			double Aj()const{return this->Aj_;}

	};



       //Store information of 1event.
       //----------------------------------
       EventInfo EVENTINFO;


	public:

       Container():SumWeight(0.0), SumPair(0.0), SumTrig(0.0){
cout << "Calling Container." << endl;
	       if(constants::MODE.find("twopc")!=string::npos){
		       Hist2D = new double *[constants::x_cell_capa];
		       Hist2D_x= new double *[constants::x_cell_capa];
		       Hist2D_y= new double *[constants::x_cell_capa];
		       Hist2DPartHit= new double *[constants::x_cell_capa];
		       Final2DHist= new double *[constants::x_cell_capa];

		       for(int i=0; i<constants::x_cell_capa; i++){
			       Hist2D[i] = new double[constants::y_cell_capa];
			       Hist2D_x[i]= new double[constants::y_cell_capa];
			       Hist2D_y[i]= new double[constants::y_cell_capa];
			       Hist2DPartHit[i]= new double[constants::y_cell_capa];
			       Final2DHist[i]= new double[constants::y_cell_capa];
		       }
	       }

	       for(int i=0; i<constants::x_cell_capa; i++){
		       for(int j=0; j<constants::y_cell_capa; j++){
			       Hist2D[i][j]=0.0;
			       Hist2D_x[i][j]=0.0;
			       Hist2D_y[i][j]=0.0;
			       Hist2DPartHit[i][j]=0.0;
			       Final2DHist[i][j]=0.0;
		       }
	       }


       };

       ~Container(){
	       cout << "Calling Deconstructore of Container." << endl;
	       //Free each sub-array
	       if(constants::MODE.find("twopc")!=string::npos){
		       for(int i = 0; i < constants::x_cell_capa; i++) {
			       delete[] Hist2D[i];
			       delete[] Hist2D_x[i];
			       delete[] Hist2D_y[i];
			       delete[] Hist2DPartHit[i];
			       delete[] Final2DHist[i];
		       }
		       //Free the array of pointers
		       delete[] Hist2D;
		       delete[] Hist2D_x;
		       delete[] Hist2D_y;
		       delete[] Hist2DPartHit;
		       delete[] Final2DHist;
	       }
       };

	double Hist[constants::x_cell_capa]={};
	double Hist_sub[constants::x_cell_capa]={};
	double Hist_subsub[constants::x_cell_capa]={};
	double Hist_x[constants::x_cell_capa]={};
	double Hist_1ev[constants::x_cell_capa]={};
	double Hist_weight[constants::x_cell_capa]={};
	double HistHist[constants::x_cell_capa]={};
	double HistHist_sub[constants::x_cell_capa]={};
	double HistHist_subsub[constants::x_cell_capa]={};

	//For write out
	//---------------
	double HistErr[constants::x_cell_capa]={};
	double FinalHist[constants::x_cell_capa]={};
	double HistErr_vn[constants::x_cell_capa]={};
	double FinalHist_vn[constants::x_cell_capa]={};

	double Hist_img_Qvec[constants::x_cell_capa]={};
	int HistHit[constants::x_cell_capa]={};
	double HistPartHit[constants::x_cell_capa]={};
	double SumWeight;
	double SumPair;
	double SumTrig;
	int max_nx=-1;
	int max_ny=-1;

	double **Hist2D;
	double **Hist2D_x;
	double **Hist2D_y;
	double **Hist2DPartHit;
	double **Final2DHist;

};
#endif
