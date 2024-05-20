#include "BSTR.h"

void BSTR::fill_iBSTR(const int iCent, std::shared_ptr<Container>& ct_in){
std::cout << "In BSTR filling " << iCent << std::endl;
	Container& ct = this->ct_ALL[iCent];
	std::cout << __FILE__ << "  Filling to BSTR: " << ct.SumWeight << std::endl;
        ct.max_nx = ct_in->max_nx;
	for(int i=0; i<ct.max_nx+1; ++i){
		//Summing up iBSTR
		//================
		ct.FinalHist[i] += ct_in->FinalHist[i];
		ct.Hist_x[i] += ct_in->Hist_x[i];
		ct.HistHist[i] += pow(ct_in->FinalHist[i], 2);
	}

	return;
}

void BSTR::stat_iBSTR(const int iCent){

	Container& ct = this->ct_ALL[iCent];
	std::cout << __FILE__ << "  Stat with nBSTR: " << this->nBSTR << std::endl;
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
