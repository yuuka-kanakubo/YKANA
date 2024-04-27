#ifndef UTIL_FUNC
#define UTIL_FUNC
#include <time.h>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <iomanip>
#include <memory>
#include "Constants.h"
#include "Rndom.h"
#include "Container.h"
#include "EbyeInfo.h"
#include "Options.h"

using std::string;
using std::cout;
using std::endl;

class Util_func{


 private:

	 std::string DATE;
	 std::shared_ptr<Rndom>& rndom;

 


 struct stat st;

 const std::string generateS2(int n){
   std::ostringstream name;
   name  << std::setw(2) << std::setfill('0') << n ;
   return name.str();
 };
 
  const std::string get_date(){
	  return this->DATE;
  };


  void get_nowdata(){
	  time_t now = time(NULL);
	  struct tm *pnow = localtime(&now);  
	  int Y=pnow->tm_year+1900;
	  int M=pnow->tm_mon + 1;
	  int D=pnow->tm_mday;
	  std::string date=generateS(Y)+generateS2(M)+generateS2(D);    
	  this->DATE=date;
  }


 public:

  Util_func(std::shared_ptr<Rndom>& rndom_in):DATE(""), rndom(rndom_in){
		this->get_nowdata(); 
 };
  ~Util_func(){};
  
  std::string generateS(double n){
    std::ostringstream oss;
    oss <<  n ; 
    return oss.str();
  };
 


  bool checkMassOnShell(const double m, const double e, const double px, const double py, const double pz){

	  double P_squared=px*px+py*py+pz*pz;
	  double E_squared=e*e;

	  bool mos = fabs(m*m - E_squared + P_squared)<constants::SMALL;

	  if(mos){
		  //std::cout << ":) This is mass on shell " << std::fixed << std::setprecision(10) << "E_squared: " << E_squared << "   P_squared: " << P_squared<< "   M_squared: " << m*m << std::endl;
		  return true;
	  }else{
		  std::cout << ":O This is NOT mass on shell " << std::fixed << std::setprecision(10) << "E: " << e << "   P : " << sqrt(P_squared) << "   M: " << m << std::endl;
		  std::cout << "                             " << std::fixed << std::setprecision(10) << "E_squared: " << E_squared << "   P_squared: " << P_squared<< "   M_squared: " << m*m << "    mm - ee + pp " <<  fabs(m*m - E_squared + P_squared)  << std::endl;
		  return false;
	  }
  }



  
  void make_output_directory(const std::string name_){
    
    if(stat(name_.c_str(),&st)!=0) mkdir(name_.c_str(), 0775);
    else {

    };
  };

  std::string get_output_directory(const std::string directory_name){
    
    std::ostringstream name; 
    name  << this->get_date()+"_"+directory_name;
    std::string name_=constants::data_directory+"/"+name.str();
    return name_.c_str();

  };


};
#endif
