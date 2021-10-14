#ifndef FILL
#define FILL
#include "Constants.h"
#include "Util_func.h"
#include "Container.h"
#include "Message.h"
#include "LogSettings.h"
#include "Settings.h"
#include "CentralityCut.h"
#include "InfoHist.h"

class Fill{

private:

		shared_ptr<Message> ms;  
		Settings::Options options;
		shared_ptr<InfoHist> infohist;
		shared_ptr<Util_func> uf;

			bool fix_ax(const int id, int &nx, double m);
			int get_cell_index(const double x_val_);
			int get_cell_index_logplot(const double x_val_);
			double getDeltaPhi(const double phi1, const double phi2);
			double getDeltaPhi_twopc1D(const double phi1, const double phi2);
	

public:

Fill(shared_ptr<Message> ms_in, Settings::Options options_in, shared_ptr<InfoHist> info, shared_ptr<Util_func> uf);
~Fill();

			void fill_jets(shared_ptr<Container>& ct);
			void fill_twopc_tagged(shared_ptr<Container>& ct);
			void fill_twopc(shared_ptr<Container>& ct);
			void fill_Rt(shared_ptr<Container>& ct);
			void fill_RtYield(shared_ptr<Container>& ct);
			void fill_twopc1D_taggedInteg(shared_ptr<Container>& ct);
			void fill_twopc1D_tagged(shared_ptr<Container>& ct);
			void fill_twopc1D_tagged_1particle(shared_ptr<Container>& ct);
			void fill_twopc1D(shared_ptr<Container>& ct);
			void fill(shared_ptr<Container>& ct);
			void fill_vn4multi(shared_ptr<Container>& ct);
			void fill_vnmulti(shared_ptr<Container>& ct);
			void fill_vneta(shared_ptr<Container>& ct);
			void fill_vnpt(shared_ptr<Container>& ct);
			void fill_vnpt_2sub(shared_ptr<Container>& ct);
			void fill_vnmulti_2sub(shared_ptr<Container>& ct);
			void fill_vn4multi_2sub(shared_ptr<Container>& ct);
			void fill_vn4multi_3sub(shared_ptr<Container>& ct);

};
#endif
