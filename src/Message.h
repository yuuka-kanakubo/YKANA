#ifndef MESSAGE
#define MESSAGE
#include <iostream>
#include <fstream>
#include <sstream>
#include "Constants.h"

using std::cout;
using std::endl;
using std::ifstream;

class Message{

 private:
  
 public:
  
 Message(){};
  ~Message(){};
  

  const bool enough_argument(const int& argc){
    if(argc< 5){
      cout << "ERROR :( Need option. " << endl;
      cout << "         ./analysis -n *  -output  *  " << endl;
      return false;
    }else return true;
  };

  const void read(const int& ifile){
    cout << "Reading ;) " << ifile << " files.." << endl;
  };
  
  const void open(const std::string& fname){
      cout << "ERROR :(  unable to open file. " << fname << endl;
  }

  const void finish(){
    cout << "Done :) " << endl;
  };

  const void TestMode(){
    cout << "MODE: " << constants::MODE << endl;
 }


  const void WARNING_LARGE_IMAGINARYPART(const double IMG){
	  cout << "WARNING :o  You've got large value for imaginary part of Q vector. " << IMG << endl;
  }

};
#endif
