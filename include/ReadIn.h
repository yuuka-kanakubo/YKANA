#ifndef READ
#define READ

class ReadIn{

private:
		shared_ptr<Message> ms;  
		Options options;
		int ncall_readTimeLapse;
		int nline;

class print_coodTL{

public:
	double x_print, y_print, eta_print, tau_print;
};


public:

ReadIn(shared_ptr<Message> ms_in, Options options_in);
~ReadIn();

bool read(const std::string& fname, shared_ptr<Container>& ct);
bool readEKRT(const std::string& fname, shared_ptr<Container>& ct);
bool readEKRTbinary(std::vector <Container::EventInfo>& nEventInfo);
bool readXY(const std::string& fname, shared_ptr<Container>& ct);
bool read_jetinfo(const std::string& fname, shared_ptr<Container>& ct);
bool readTimeLapse(const std::string& fname, shared_ptr<Container>& ct, const double weight);
bool get_nline_to_see(int &nline, const std::string fname);
void get_oneline_xeta(istringstream& is, Container::StepInfo& onestep);
void get_oneline_xy(istringstream& is, Container::StepInfo& onestep);
};
#endif
