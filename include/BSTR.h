#include <memory>
#include "Container.h"
#include "Options.h"

class BSTR{

private:
	int nBSTR;
	Options& options;

public:
	std::vector<Container> ct_ALL;

	BSTR(Options& options_in):options(options_in){
                this->nBSTR=options.get_nBSTR();
		for(int iCent=0; iCent<(int)options.name_cent.size(); iCent++){
					Container ct_iBSTR(this->options.get_flag_SB_CMS());
					ct_ALL.push_back(ct_iBSTR);
		}
                std::cout << "size of ct_ALL " << (int) ct_ALL.size() << std::endl;
	};
	~BSTR(){};

	void fill_iBSTR(const int iCent, std::shared_ptr<Container>& ct_in);
	void stat_iBSTR(const int iCent);



};
