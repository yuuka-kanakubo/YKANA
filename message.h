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
  

  bool enough_argument(const int& argc){
    if(argc< 5){
      cout << "Need option. " << endl;
      cout << "./analysis -n *  -output  *  " << endl;
      return false;
    }else return true;
  };

  
  void open(const std::string fname){
      cout << "Error unable to open file. " << fname << endl;
  }

  void finish(){
    cout << "Done." << endl;
  };


};
