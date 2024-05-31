// -*- mode:c++ -*-
#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <sstream>
#include <fstream>
#include <complex>
#include <math.h>
#include "Constants.h"
#include "Util_func.h"
#include "Container.h"
#include "Message.h"
#include "LogSettings.h"
#include "Settings.h"
#include "CentralityCut.h"
#ifndef WRITE
#define WRITE

class Write{

private:

		shared_ptr<Message> ms;  
		Options options;
		shared_ptr<Util_func> uf;

public:

Write(shared_ptr<Message> ms_in, Options options_in, shared_ptr<Util_func> uf);
~Write();
		int getMapEdgeX(const double maxval);
		int getMapEdgeY(const double maxval);
		bool write_BSTR(const std::string& fname, const Container& ct);
		bool write(const std::string& fname, const shared_ptr<Container>& ct);
		bool write_RtYield(const std::string& fname, const shared_ptr<Container>& ct);

};
#endif
