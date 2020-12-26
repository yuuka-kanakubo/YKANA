#include "Constants.h"

class Container{

 class ParticleInfo{
	 public:

   double e;
   double px;
   double py;
   double pz;
   double x;
   double y;
   double z;
   double t;
   double rap;
   double eta;
   double phi;

 };

public:

  class Event{

	  public:
		  vector<ParticleInfo> part;
          
                  vector<double> value;

  };






 public:
  double Hist[constants::x_cell_capa]={};
  double Hist_1ev[constants::x_cell_capa]={};
  double HistHist[constants::x_cell_capa]={};
  double HistErr[constants::x_cell_capa]={};
  int max_nx=-1;
  int Nev_tot=0;
  double sum_weight=0.0;

};
