#ifndef STAT
#define STAT

class Stat{

private:

		shared_ptr<Message> ms;  
		Settings::Options options;
		shared_ptr<InfoHist> infohist;
		shared_ptr<Util_func> uf;
		int get_xaxis_RtClass(double xval);

public:

Stat(shared_ptr<Message> ms_in, Settings::Options options_in, shared_ptr<InfoHist> info, shared_ptr<Util_func> uf);
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
 
 StatCumulant():
	 qvec2corr(0.0),
	 numpair2corr(0.0),
	 qvec4corr(0.0),
	 numpair4corr(0.0),
	 qvec2corr_err(0.0),
	 numpair2corr_err(0.0),
	 qvec4corr_err(0.0),
	 numpair4corr_err(0.0)
	{};

 double qvec2corr;
 double numpair2corr;
 double qvec4corr;
 double numpair4corr;
 double qvec2corr_err;
 double numpair2corr_err;
 double qvec4corr_err;
 double numpair4corr_err;

};


};
#endif
