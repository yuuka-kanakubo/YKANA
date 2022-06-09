#ifndef EBYEINFO
#define EBYEINFO
#include "Container.h"

class EbyeInfo{
 public:

     EbyeInfo():multiplicity(-1.0),multiplicity_V0M(-1.0),N_trk_offline(-1.0),Nch(-1.0), weight(-1.0), multiplicity_INEL_lg_0(false),
	     trig_3outof3(false), trig_2outof3(false), ATLAS_cut(false),trig_VZEROAND(false), valid(false), valid_assoc(false), orig_eventNum(-1), V0M_class(-1){};

     double multiplicity;
     double multiplicity_V0M;
	double N_trk_offline;
	double Nch;
	double weight;

	bool multiplicity_INEL_lg_0;
	bool trig_3outof3;
	bool trig_2outof3;
	bool ATLAS_cut;
	bool trig_VZEROAND;
        bool valid;
        bool valid_assoc;
	Container::ParticleInfo sample_part;
	vector<Container::ParticleInfo> sample_partSet;

	int orig_eventNum;

	void set_V0M_class(const int CLASS){this->V0M_class=CLASS;}
	int get_V0M_class()const{return this->V0M_class;}


	bool operator > (const EbyeInfo& event_info) const {
		return (multiplicity_V0M > event_info.multiplicity_V0M);
	}

 private:
	int V0M_class;

};
#endif
