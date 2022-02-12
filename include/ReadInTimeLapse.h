#ifndef READTL
#define READTL

class ReadInTimeLapse{

private:
		shared_ptr<Message> ms;  
		Settings::Options options;

public:

ReadInTimeLapse(shared_ptr<Message> ms_in, Settings::Options options_in);
~ReadInTimeLapse();

bool read(const std::string& fname, shared_ptr<Container>& ct);

};
#endif
