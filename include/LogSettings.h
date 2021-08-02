#ifndef LOGSETTINGS
#define LOGSETTINGS
#include <iostream>
#include <fstream>
#include "Constants.h"

using namespace std;

class LogSettings{

 public:

vector<string> options;

int archive_settings(const string output_directory_name){


	{
		ifstream ifs(output_directory_name+"/"+constants::settings_outputfname);
		if(ifs.is_open()) {
			return 0;
		}
	} 


    vector<string> templine_set;


    templine_set.push_back("//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-");
    templine_set.push_back("");
    for(int i=0; i<(int)options.size(); ++i){
	    templine_set.push_back(options[i]);
    }

    templine_set.push_back("");

    for(int i=0; i<(int)constants::save_settings_fname.size(); ++i){
	    {
		    templine_set.push_back("//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-");
		    templine_set.push_back("");
		    ifstream in;
		    in.open(constants::save_settings_fname[i], ios::in);	
		    if(!in) {
			    cout << "Error unable to open file " << constants::save_settings_fname[i] << endl;
			    return 1;
		    }    
		    string templine;
		    while(getline(in,templine)) {
			    istringstream ist(templine);
			    templine_set.push_back(templine);
		    }
		    templine_set.push_back("");
	    }
    }

    ofstream ofs;
    ofs.open(output_directory_name+"/"+constants::settings_outputfname,ios::in|ios::app);
    for(int i=0; i<(int)templine_set.size(); ++i)
     ofs << templine_set[i] << endl;

    ofs.close();
    return 0;
  }

};
#endif
