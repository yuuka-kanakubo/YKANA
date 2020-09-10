#include <time.h>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <iomanip>
#include "constants.h"

class util_func{

 public:
  util_func(){};
  ~util_func(){};
  
  std::string generateS(double n){
    std::ostringstream oss;
    oss <<  n ; 
    return oss.str();
  };
 
  std::string get_name_directory(){
    return name_directory;
  };
 
 private: 

 struct stat st;
 std::string name_directory;

 const std::string generateS2(int n){
   std::ostringstream name;
   name  << std::setw(2) << std::setfill('0') << n ;
   return name.str();
 };
 
  const std::string get_date(){
    time_t now = time(NULL);
    struct tm *pnow = localtime(&now);  
    int Y=pnow->tm_year+1900;
    int M=pnow->tm_mon + 1;
    int D=pnow->tm_mday;
    std::string date=generateS(Y)+generateS2(M)+generateS2(D);    
    return date;
  };

  
 public:
  const void make_output_directory(const std::string& directory_name){
    
    this->make_data_directory();
    
    std::ostringstream name; 
    name  << this->get_date()+"_"+directory_name;
    name_directory=constants::data_directory+"/"+name.str();
    if(stat(name_directory.c_str(),&st)!=0) mkdir(name_directory.c_str(),0775);
    else {};

  };

 private:
  const void make_data_directory(){
    if(stat(constants::data_directory.c_str(),&st) !=0) mkdir(constants::data_directory.c_str(),0775);
    else {};
  };


};

