#include <iostream>
#include "Settings.h"
#include "Analysis.h"
#include "Message.h"

int main(int argc, char* argv[]){

	auto ms = make_shared<Message>();
	ms->TestMode();

	auto setting = make_shared<Settings>(argc, argv);
        Analysis analysis(setting->options, setting->log);

	ms->finish();

	return 0;
}
