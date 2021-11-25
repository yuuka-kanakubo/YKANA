#ifndef LOGSETTINGS
#define LOGSETTINGS
#include <iostream>
#include <fstream>
#include "Constants.h"

using namespace std;

class LogSettings{


 private:

	 vector<string> BinSettings_str;
	 bool centrality_cut;


 public:

	 vector<string> options;

	 int get_BinSettings_size(){return (int)BinSettings_str.size();}
	 void set_centrality_cut(bool flag){centrality_cut=flag;}


void save_BinSettings(string input){


	    {
		    BinSettings_str.push_back("//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-");
		    BinSettings_str.push_back("//Bins..");
		    ifstream in;
		    in.open(input, ios::in);
		    if(!in) {
			    cout << "Error unable to open file " << input << endl;
			    return;
		    }    
		    string str;
		    while(getline(in,str)) {
			    istringstream ist(str);
			    BinSettings_str.push_back(str);
		    }
		    BinSettings_str.push_back("");
	    }

cout << "Saved bin settings " << (int)BinSettings_str.size() << endl;

}



bool archive_settings(const string output_directory_name){


	{
		ifstream ifs(output_directory_name+"/"+constants::settings_outputfname);
		if(ifs.is_open() && this->centrality_cut) {
			return true;
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
			    return false;
		    }    
		    string templine;
		    while(getline(in,templine)) {
			    istringstream ist(templine);
			    templine_set.push_back(templine);
		    }
		    templine_set.push_back("");
	    }
    }




    //Output
    //--------
    ofstream ofs;
    ofs.open(output_directory_name+"/"+constants::settings_outputfname,ios::in|ios::app);
    for(int i=0; i<(int)templine_set.size(); ++i)
	    ofs << templine_set[i] << endl;
    cout << __FILE__ << "log " << get_BinSettings_size() << endl;
    for(int i=0; i<(int)BinSettings_str.size(); ++i){
	    ofs << BinSettings_str[i] << endl;
    }

    ofs.close();
    return true;
  }

};
#endif
