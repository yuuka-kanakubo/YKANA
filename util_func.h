#include <time.h>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <iomanip>

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

  std::string generateS2(int n){
    std::ostringstream name;
    name  << std::setw(2) << std::setfill('0') << n ;
    return name.str();
  };

  std::string get_date(){
    time_t now = time(NULL);
    struct tm *pnow = localtime(&now);  
    int Y=pnow->tm_year+1900;
    int M=pnow->tm_mon + 1;
    int D=pnow->tm_mday;
    std::string date=generateS(Y)+generateS2(M)+generateS2(D);    
    return date;
  };

  
 public:
  void make_output_direcrory(const std::string& directory_name){

    this->make_data_direcrory();

    std::ostringstream name;
    name  << this->get_date()+"_"+directory_name;
    name_directory=name.str();
    if(stat(name_directory.c_str(),&st)!=0) mkdir(name_directory.c_str(),0775);
    else {};

  };

 private:
  void make_data_direcrory(){
    if(stat(constants::data_directory.c_str(),&st) !=0) mkdir(constants::data_directory.c_str(),0775);
    else {};
  };


};

