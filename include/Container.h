#ifndef CONTAINER
#define CONTAINER
#include <iostream>
#include "Constants.h"


class Container{
 
   private:

	   bool flag_SB_CMS;

	public:

	class ParticleInfo{

		public:

			double e;
			double m;
			double mt;
			double px;
			double py;
			double pz;
			double pt;
			double vt;
			double x;
			double y;
			double z;
			double t;
			double tau;
			double r;
			double rap;
			double eta;
			double phi;
			std::string TAG;
			double tata;
			uint_fast16_t ID;
			uint_fast16_t momID1;
			uint_fast16_t momID2;
			uint_fast16_t imomNucleon1;
			uint_fast16_t imomNucleon2;
			double xmomNucleon1;
			double xmomNucleon2;
			double ymomNucleon1;
			double ymomNucleon2;
			double zmomNucleon1;
			double zmomNucleon2;
			bool is_mom1Neutron;
			bool is_mom2Neutron;

			ParticleInfo():e(0.0), m(0.0), mt(0.0), px(0.0), py(0.0), pz(0.0), pt(0.0), vt(0.0), x(0.0), y(0.0), z(0.0), t(0.0), tau(0.0), r(0.0), rap(0.0), eta(0.0), phi(0.0), TAG("na"), tata(0.0), ID(0), momID1(0), momID2(0), imomNucleon1(0), imomNucleon2(0), xmomNucleon1(0.0), xmomNucleon2(0.0), ymomNucleon1(0.0), ymomNucleon2(0.0), zmomNucleon1(0.0), zmomNucleon2(0.0), is_mom1Neutron(false), is_mom2Neutron(false){}


	};

        class StepInfo{

		public:

			StepInfo():nstep(0), tau(0.0),x(0.0),y(0.0), eta(0.0), e(0.0),temp(0.0){};
			~StepInfo(){};

			int nstep;
			double tau;
			double x;//fm
			double y;//fm
			double eta;
			double e; //[GeV^4]
			double temp; //[GeV]
			double s; //[GeV^3]
			double p; //[GeV]
			double n; //[GeV^3]
			double n5; 
			double mu; 
			double mu5;
			double vx;  //[1]
			double vy;  //[1]
			double veta;  //[1]
			double Vtilde;
			double UABS;
			double U4;
			double U_R4;
			double U_L4;
			double Ux;
			double Uy;
			double Ueta;
			double eEtildedx;
			double eEtildedy;
			double eBtildedx;
			double eBtildedy;
			double eEtilde_dot_eBtilde;

	};


	//Info of 1 event.
	//---------------
	class EventInfo{

		private:

			double weight_;
			int Nt_;
			int Nch_;
			int order_reading_;
			double Aj_;

		public:
			EventInfo(): weight_(1.0), Nt_(0), Nch_(0), order_reading_(0), Aj_(1.0){};
			~EventInfo(){
				std::vector<ParticleInfo>().swap(part);
				std::vector<StepInfo>().swap(step);
			};

			std::vector<ParticleInfo> part;
			std::vector<StepInfo> step;
			void weight(double weight_in){this->weight_=weight_in;}
			void Nt(int Nt_in){this->Nt_=Nt_in;}
			void Nch(int Nch_in){this->Nch_=Nch_in;}
			void Aj(double Aj_in){this->Aj_=Aj_in;}
			void order_reading(int i){this->order_reading_=i;}
			double weight()const{return this->weight_;}
			int Nch()const{return this->Nch_;}
			int order_reading()const{return this->order_reading_;}
			double Aj()const{return this->Aj_;}

	};



       //Store information of 1event.
       //----------------------------------
       EventInfo EVENTINFO;


	public:

       Container(bool SB_CMS):flag_SB_CMS(SB_CMS),SumWeight(0.0), SumPair(0.0), SumTrig(0.0), CountEv(0), B00(0.0), meanNt(-1.0){
	       //std::cout << "Calling Container." << std::endl;
	       if(constants::MODE.find("twopc")!=std::string::npos
			       || constants::MODE.find("twodm")!=std::string::npos){
		       Hist2D = new double *[constants::x_cell_capa];
		       HistSub2D = new double *[constants::x_cell_capa];
		       Hist2D_x= new double *[constants::x_cell_capa];
		       Hist2D_y= new double *[constants::x_cell_capa];
		       Hist2DPartHit= new double *[constants::x_cell_capa];
		       if(flag_SB_CMS){
			       HistSub2D_x= new double *[constants::x_cell_capa];
			       HistSub2D_y= new double *[constants::x_cell_capa];
			       HistSub2DPartHit= new double *[constants::x_cell_capa];
		       }
		       Final2DHist= new double *[constants::x_cell_capa];
		       Final2DHistSub= new double *[constants::x_cell_capa];
		       Final2DHit= new double *[constants::x_cell_capa];

		       for(int i=0; i<constants::x_cell_capa; i++){
			       Hist2D[i] = new double[constants::y_cell_capa];
			       HistSub2D[i] = new double[constants::y_cell_capa];
			       Hist2D_x[i]= new double[constants::y_cell_capa];
			       Hist2D_y[i]= new double[constants::y_cell_capa];
			       Hist2DPartHit[i]= new double[constants::y_cell_capa];
			       if(flag_SB_CMS){
				       HistSub2D_x[i]= new double[constants::y_cell_capa];
				       HistSub2D_y[i]= new double[constants::y_cell_capa];
				       HistSub2DPartHit[i]= new double[constants::y_cell_capa];
			       }
			       Final2DHist[i]= new double[constants::y_cell_capa];
			       Final2DHistSub[i]= new double[constants::y_cell_capa];
			       Final2DHit[i]= new double[constants::y_cell_capa];
		       }

		       for(int i=0; i<constants::x_cell_capa; i++){
			       for(int j=0; j<constants::y_cell_capa; j++){
				       Hist2D[i][j]=0.0;
				       HistSub2D[i][j]=0.0;
				       Hist2D_x[i][j]=0.0;
				       Hist2D_y[i][j]=0.0;
				       Hist2DPartHit[i][j]=0.0;
				       if(flag_SB_CMS){
					       HistSub2D_x[i][j]=0.0;
					       HistSub2D_y[i][j]=0.0;
					       HistSub2DPartHit[i][j]=0.0;
				       }
				       Final2DHist[i][j]=0.0;
				       Final2DHistSub[i][j]=0.0;
				       Final2DHit[i][j]=0.0;
			       }
		       }

	       }else if(constants::MODE.find("Rt_yield")!=std::string::npos){
		       RtHist_RtTrans_yield = new double *[constants::num_of_Species_Rt];
		       RtHist_RtToward_yield = new double *[constants::num_of_Species_Rt];
		       RtHist_RtTrans_yieldyield = new double *[constants::num_of_Species_Rt];
		       RtHist_RtToward_yieldyield = new double *[constants::num_of_Species_Rt];
		       HistHit_Rt = new double *[constants::num_of_Species_Rt];
		       HistErrTrans_Rt = new double *[constants::num_of_Species_Rt];
		       HistErrToward_Rt = new double *[constants::num_of_Species_Rt];
		       for(int i=0; i<constants::num_of_Species_Rt; i++){
			       RtHist_RtTrans_yield[i] = new double [constants::x_cell_capa];
			       RtHist_RtToward_yield[i] = new double [constants::x_cell_capa];
			       RtHist_RtTrans_yieldyield[i] = new double [constants::x_cell_capa];
			       RtHist_RtToward_yieldyield[i] = new double [constants::x_cell_capa];
			       HistHit_Rt[i] = new double [constants::x_cell_capa];
			       HistErrTrans_Rt[i] = new double [constants::x_cell_capa];
			       HistErrToward_Rt[i] = new double [constants::x_cell_capa];
		       }

		       for(int j=0; j<constants::num_of_Species_Rt; j++){
			       for(int i=0; i<constants::x_cell_capa; i++){
				       RtHist_RtTrans_yield[j][i]=0.0;
				       RtHist_RtToward_yield[j][i]=0.0;
				       RtHist_RtTrans_yieldyield[j][i]=0.0;
				       RtHist_RtToward_yieldyield[j][i]=0.0;
				       HistHit_Rt[j][i]=0.0;
				       HistErrTrans_Rt[j][i]=0.0;
				       HistErrToward_Rt[j][i]=0.0;
			       }
		       }

	       }


       };

       ~Container(){
	       //std::cout << "Calling Deconstructore of Container." << std::endl;
	       //Free each sub-array
	       if(constants::MODE.find("twopc")!=std::string::npos
			       || constants::MODE.find("twodm")!=std::string::npos){
		       for(int i = 0; i < constants::x_cell_capa; i++) {
			       delete[] Hist2D[i];
			       delete[] HistSub2D[i];
			       delete[] Hist2D_x[i];
			       delete[] Hist2D_y[i];
			       delete[] Hist2DPartHit[i];
			       if(flag_SB_CMS){
				       delete[] HistSub2D_x[i];
				       delete[] HistSub2D_y[i];
				       delete[] HistSub2DPartHit[i];
			       }
			       delete[] Final2DHist[i];
			       delete[] Final2DHistSub[i];
			       delete[] Final2DHit[i];
		       }
		       //Free the array of pointers
		       delete[] Hist2D;
		       delete[] HistSub2D;
		       delete[] Hist2D_x;
		       delete[] Hist2D_y;
		       delete[] Hist2DPartHit;
		       if(flag_SB_CMS){
			       delete[] HistSub2D_x;
			       delete[] HistSub2D_y;
			       delete[] HistSub2DPartHit;
		       }
		       delete[] Final2DHist;
		       delete[] Final2DHistSub;
		       delete[] Final2DHit;
}else if(constants::MODE.find("Rt_yield")!=std::string::npos){
		       for(int i = 0; i < constants::num_of_Species_Rt; i++) {
			       delete[] RtHist_RtTrans_yield[i];
			       delete[] RtHist_RtToward_yield[i];
			       delete[] RtHist_RtTrans_yieldyield[i];
			       delete[] RtHist_RtToward_yieldyield[i];
			       delete[] HistHit_Rt[i];
			       delete[] HistErrTrans_Rt[i];
			       delete[] HistErrToward_Rt[i];
		       }
		       delete[] RtHist_RtTrans_yield;
		       delete[] RtHist_RtToward_yield;
		       delete[] RtHist_RtTrans_yieldyield;
		       delete[] RtHist_RtToward_yieldyield;
		       delete[] HistHit_Rt;
		       delete[] HistErrTrans_Rt;
		       delete[] HistErrToward_Rt;
	       }
       };

	double Hist[constants::x_cell_capa]={};
	double Hist_sub[constants::x_cell_capa]={};
	double Hist_subsub[constants::x_cell_capa]={};
	double Hist2[constants::x_cell_capa]={};
	double Hist2_sub[constants::x_cell_capa]={};
	double Hist2_subsub[constants::x_cell_capa]={};
	double Hist_x[constants::x_cell_capa]={};
	double Hist_1ev[constants::x_cell_capa]={};
	double Hist_weight[constants::x_cell_capa]={};
	double HistHist[constants::x_cell_capa]={};
	double HistHist_sub[constants::x_cell_capa]={};
	double HistHist_subsub[constants::x_cell_capa]={};
	double HistHist2[constants::x_cell_capa]={};
	double HistHist2_sub[constants::x_cell_capa]={};
	double HistHist2_subsub[constants::x_cell_capa]={};

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
	int CountEv;
	double B00;

	//Special case: Rt analysis needs eBye Nt. So need to store info of all events.
	//------------------------------------------------------------------
	class yield{
		private:
			double _chpi;
			double _chkaon;
			double _ppbar;
			double _cascade;
			double _lambda;
			double _phi;
			double _omega;

		public:
			yield():_chpi(0.0), _chkaon(0.0), _ppbar(0.0), _cascade(0.0), _lambda(0.0), _phi(0.0), _omega(0.0){};

			double chpi(){return _chpi;}
			double chkaon(){return _chkaon;}
			double ppbar(){return _ppbar;}
			double cascade(){return _cascade;}
			double phi(){return _phi;}
			double lambda(){return _lambda;}
			double omega(){return _omega;}

			double get_sp(int sp){
				if(sp==0) return _chpi;
				else if(sp==1) return _chkaon;
				else if(sp==2) return _ppbar;
				else if(sp==3) return _phi;
				else if(sp==4) return _lambda;
				else if(sp==5) return _cascade;
				else if(sp==6) return _omega;
				else {
					std::cout << "ERROR :( out of range in get_sp" << sp << std::endl;
					exit(1);
				}
			}

			std::string get_particleName(int sp){
				if(sp==0) return "chpi";
				else if(sp==1) return "chkaon";
				else if(sp==2) return "ppbar";
				else if(sp==3) return "phi";
				else if(sp==4) return "lambda";
				else if(sp==5) return "cascade";
				else if(sp==6) return "omega";
				else {
					std::cout << "ERROR :( out of range in get_sp" << sp << std::endl;
					exit(1);
				}
			}

			void chpi(double val_in){ _chpi = val_in;}
			void chkaon(double val_in){ _chkaon= val_in;}
			void ppbar(double val_in){ _ppbar= val_in;}
			void cascade(double val_in){ _cascade= val_in;}
			void phi(double val_in){ _phi= val_in;}
			void lambda(double val_in){ _lambda= val_in;}
			void omega(double val_in){ _omega= val_in;}

			void add_chpi(double val_in){ _chpi += val_in;}
			void add_chkaon(double val_in){ _chkaon+= val_in;}
			void add_ppbar(double val_in){ _ppbar+= val_in;}
			void add_cascade(double val_in){ _cascade+= val_in;}
			void add_phi(double val_in){ _phi+= val_in;}
			void add_lambda(double val_in){ _lambda+= val_in;}
			void add_omega(double val_in){ _omega+= val_in;}

	};


	std::vector<int> Nt_eBye;
	std::vector<int> Ntmin_eBye;
	std::vector<int> Ntmax_eBye;
	std::vector<int> TagEventNum;
	std::vector<double> weight_eBye;
	std::vector<double> dNdeta_eBye;
	std::vector<double> CoreT_eBye;
	std::vector<double> CoreN_eBye;
	std::vector<yield> TransYield_eBye;
	std::vector<yield> TowardYield_eBye;
	double meanNt;
	double **RtHist_RtTrans_yield;
	double **RtHist_RtToward_yield;
	double **RtHist_RtTrans_yieldyield;
	double **RtHist_RtToward_yieldyield;
	double **HistHit_Rt;
	double **HistErrTrans_Rt;
	double **HistErrToward_Rt;

	int max_nx=-1;
	int max_ny=-1;

	double **Hist2D;
	double **Hist2D_x;
	double **Hist2D_y;
	double **Hist2DPartHit;
	double **HistSub2D;
	double **HistSub2D_x;
	double **HistSub2D_y;
	double **HistSub2DPartHit;
	double **Final2DHist;
	double **Final2DHistSub;
	double **Final2DHit;

};
#endif
