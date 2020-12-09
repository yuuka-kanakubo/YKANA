#include "Constants.h"

class Container{

 public:
  double Hist[constants::x_cell_capa]={};
  double HistHist[constants::x_cell_capa]={};
  double HistErr[constants::x_cell_capa]={};
  int max_nx=-1;
  int Nev_tot=0;
  double sum_weight=0.0;

};
