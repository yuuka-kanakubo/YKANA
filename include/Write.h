#ifndef WRITE
#define WRITE

class Write{

private:

		shared_ptr<Message> ms;  
		Options options;
		shared_ptr<InfoHist> infohist;
		shared_ptr<Util_func> uf;

public:

Write(shared_ptr<Message> ms_in, Options options_in, shared_ptr<InfoHist> info, shared_ptr<Util_func> uf);
~Write();
		int getMapEdgeX(const double maxval);
		int getMapEdgeY(const double maxval);
		bool write(const std::string& fname, const shared_ptr<Container>& ct);
		bool write_RtYield(const std::string& fname, const shared_ptr<Container>& ct);

};
#endif
