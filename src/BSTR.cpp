#include "BSTR.h"

void BSTR::fill_iBSTR(const int iCent, std::shared_ptr<Container>& ct_in, Container& ct_ALL_BSTRarchive){
	ct_ALL_BSTRarchive.max_nx = ct_in->max_nx;
	for(int i=0; i<ct_ALL_BSTRarchive.max_nx+1; ++i){
		//Summing up iBSTR
		//================
		ct_ALL_BSTRarchive.FinalHist[i] += ct_in->FinalHist[i];
		ct_ALL_BSTRarchive.Hist_x[i] += ct_in->Hist_x[i];
		ct_ALL_BSTRarchive.HistHist[i] += pow(ct_in->FinalHist[i], 2);
	}

	return;
}

void BSTR::stat_iBSTR(Container& ct){

	for(int i=0; i<ct.max_nx+1; ++i){
		//Summing up iBSTR
		//================
		ct.FinalHist[i]=ct.FinalHist[i]/this->nBSTR;
		ct.Hist_x[i]=ct.Hist_x[i]/this->nBSTR;

		double meanxx  = ct.HistHist[i]/this->nBSTR;
		double meanx = ct.FinalHist[i];


		double var=meanxx-pow(meanx,2.0);
		double error=sqrt(var/(this->nBSTR-1.0));
		ct.HistErr[i] = error;
	}

	return;

}
