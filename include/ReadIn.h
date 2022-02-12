#ifndef READ
#define READ

class ReadIn{

private:
		shared_ptr<Message> ms;  
		Settings::Options options;

public:

ReadIn(shared_ptr<Message> ms_in, Settings::Options options_in);
~ReadIn();

bool read(const std::string& fname, shared_ptr<Container>& ct);
bool read_jetinfo(const std::string& fname, shared_ptr<Container>& ct);
bool readTimeLapse(const std::string& fname, shared_ptr<Container>& ct, const double weight);

};
#endif
