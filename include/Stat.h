#ifndef STAT
#define STAT

class Stat{

private:

		shared_ptr<Message> ms;  
		Options options;
		shared_ptr<InfoHist> infohist;
		shared_ptr<Util_func> uf;
		int get_xaxis_RtClass(double xval);
		double get_B00(const vector<double>& X, const vector<double>& Y, const vector<double>& Z);

public:

Stat(shared_ptr<Message> ms_in, Options options_in, shared_ptr<InfoHist> info, shared_ptr<Util_func> uf);
~Stat();

			void stat_twopc(shared_ptr<Container>& ct);
			void stat_twopc1D(shared_ptr<Container>& ct);
			void stat(shared_ptr<Container>& ct);
			void stat_Rt(shared_ptr<Container>& ct);
			void stat_RtYield(shared_ptr<Container>& ct);
			void stat_jets(shared_ptr<Container>& ct);
			void stat_flow(shared_ptr<Container>& ct);

class StatCumulant{

public:
 
 StatCumulant(){
   this->init();
 };

 // Qn*Qn (qvec corr), Npart(num pair)
 //-----------------------------------
 double qvec2corr[constants::x_cell_capa]={};
 double numpair2corr[constants::x_cell_capa]={};
 double qvec2ABcorr[constants::x_cell_capa]={};
 double numpair2ABcorr[constants::x_cell_capa]={};
 double qvec2ACcorr[constants::x_cell_capa]={};
 double numpair2ACcorr[constants::x_cell_capa]={};
 double qvec4corr[constants::x_cell_capa]={};
 double numpair4corr[constants::x_cell_capa]={};
 double qvec2corr_err[constants::x_cell_capa]={};
 double numpair2corr_err[constants::x_cell_capa]={};
 double qvec2ABcorr_err[constants::x_cell_capa]={};
 double numpair2ABcorr_err[constants::x_cell_capa]={};
 double qvec2ACcorr_err[constants::x_cell_capa]={};
 double numpair2ACcorr_err[constants::x_cell_capa]={};
 double qvec4corr_err[constants::x_cell_capa]={};
 double numpair4corr_err[constants::x_cell_capa]={};

 //<<nparticle>>
 //---------------
 double Corr2[constants::x_cell_capa]={};
 double Corr2Err[constants::x_cell_capa]={};
 double Corr2AB[constants::x_cell_capa]={};
 double Corr2ABErr[constants::x_cell_capa]={};
 double Corr2AC[constants::x_cell_capa]={};
 double Corr2ACErr[constants::x_cell_capa]={};
 double Corr4[constants::x_cell_capa]={};
 double Corr4Err[constants::x_cell_capa]={};

 //cn{nparticle}
 //--------------------
 double Cumu4[constants::x_cell_capa]={};
 double Cumu4Err[constants::x_cell_capa]={};

 void init(){
	 for(int i=0; i<constants::x_cell_capa; i++){
		 qvec2corr[i]=0.0;
		 numpair2corr[i]=0.0;
		 qvec2ABcorr[i]=0.0;
		 numpair2ABcorr[i]=0.0;
		 qvec2ACcorr[i]=0.0;
		 numpair2ACcorr[i]=0.0;
		 qvec4corr[i]=0.0;
		 numpair4corr[i]=0.0;
		 qvec2corr_err[i]=0.0;
		 numpair2corr_err[i]=0.0;
		 qvec4corr_err[i]=0.0;
		 numpair4corr_err[i]=0.0;
		 Corr2[i]=0.0;
		 Corr2Err[i]=0.0;
		 Corr2AB[i]=0.0;
		 Corr2ABErr[i]=0.0;
		 Corr2AC[i]=0.0;
		 Corr2ACErr[i]=0.0;
		 Corr4[i]=0.0;
		 Corr4Err[i]=0.0;
		 Cumu4[i]=0.0;
		 Cumu4Err[i]=0.0;
	 }
 }


};


};
#endif
