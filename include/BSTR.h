#include <memory>
#include "Container.h"
#include "Options.h"

class BSTR{

private:
	int nBSTR;
	Options& options;

public:

	BSTR(Options& options_in):options(options_in){
                this->nBSTR=options.get_nBSTR();
	};
	~BSTR(){};

	void fill_iBSTR(const int iCent, std::shared_ptr<Container>& ct_in, Container& ct_ALL_BSTRarchive);
	void stat_iBSTR(Container& ct);



};
