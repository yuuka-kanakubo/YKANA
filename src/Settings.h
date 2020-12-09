#ifndef SETTINGS
#define SETTINGS
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <sstream>
#include <fstream>
#include "Constants.h"
#include "Message.h"
#include "LogSettings.h"

using std::cout;
using std::endl;

class Settings{

private:

 Message* ms;


public:

//------------------------------------------------------------------------
 class Options{

  public:

  std::string inputfname;
  std::string out_directory_name;
  std::string ext;
  std::string out_fname;
  int nfiles;

  Options():
	  inputfname(constants::default_inputfname),
	  out_directory_name(constants::default_out_directory_name),
	  ext(constants::default_ext),
	  out_fname(constants::default_out_fname),
	  nfiles(constants::default_nfiles){};

 };
//------------------------------------------------------------------------


 Options options;
 LogSettings set;

 Settings(int argc, char* argv[]){
 
 ms = new Message();
 init(argc, argv);

 };
 ~Settings(){};


void init(int argc, char* argv[]){

  //if(!ms->enough_argument(argc)) exit (1);

  for(int i=1; i<argc; i++) {
    set.options.push_back(argv[i]);
    if(!strcmp(argv[i],"-n")) {options.nfiles = atoi(argv[i+1]);}
    else if(!strcmp(argv[i],"-outdir")){options.out_directory_name= argv[i+1];}
    else if(!strcmp(argv[i],"-PATH")) {options.inputfname= argv[i+1];}
    else if(!strcmp(argv[i],"-ext")) {options.ext= argv[i+1];}
    else { 
	    std::string opt_in(argv[i]);
	    if(opt_in.find('-')==std::string::npos)continue;
	    cout << "ERROR:( There is no such an option: " << opt_in << endl; 
	    exit(1);
    }
  }
}


};
#endif
