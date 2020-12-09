#include <iostream>
#include "Settings.h"
#include "Analysis.h"
#include "Message.h"

using std::cout;
using std::endl;

int main(int argc, char* argv[]){

	Message* ms = new Message();
	ms->TestMode();

	Settings setting(argc, argv);
        Analysis analysis(setting.options, setting.set);

	ms->finish();

	return 0;
}
