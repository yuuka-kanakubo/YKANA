#ifndef INFOHIST
#define INFOHIST
class InfoHist{

public:


InfoHist(double x_max_in, double y_max_in, double d_x_in, double d_y_in, double N_coeff_in): x_max(x_max_in),  y_max(y_max_in),  d_x(d_x_in),  d_y(d_y_in),  N_coeff(N_coeff_in){};
InfoHist(){};
~InfoHist(){};

		//Maximum value for histgram
		//-------------------------
		double x_max;
		double y_max;
		double d_x;
		double d_y;
		double N_coeff;


};
#endif
