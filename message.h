#include <iostream>
#include <fstream>
#include <sstream>

using std::cout;
using std::endl;
using std::ifstream;

class message{

 private:
  
 public:
  
 message(){};
  ~message(){};
  

  const bool enough_argument(const int& argc){
    if(argc< 5){
      cout << "Need option. " << endl;
      cout << "./analysis -n *  -output  *  " << endl;
      return false;
    }else return true;
  };

  
  const void open(const std::string& fname){
      cout << "Error unable to open file. " << fname << endl;
  }

  const void finish(){
    cout << "Done." << endl;
  };


};
