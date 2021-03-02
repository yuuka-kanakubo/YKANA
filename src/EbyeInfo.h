#ifndef EBYEINFO
#define EBYEINFO
class EbyeInfo{
 public:

     EbyeInfo():multiplicity(-1.0),multiplicity_V0M(-1.0),multiplicity_INEL_lg_0(false),
	     trig_3outof3(false), trig_2outof3(false), ATLAS_cut(false),trig_VZEROAND(false), weight(-1.0), valid(false), orig_eventNum(-1), V0M_class(-1){};

     double multiplicity;
	double multiplicity_V0M;
	bool multiplicity_INEL_lg_0;
	bool trig_3outof3;
	bool trig_2outof3;
	bool ATLAS_cut;
	bool trig_VZEROAND;
	double weight;
        bool valid;

	int orig_eventNum;

	bool set_V0M_class(const int CLASS){this->V0M_class=CLASS;};
	int get_V0M_class(){return this->V0M_class;};


	bool operator > (const EbyeInfo& event_info) const {
		return (multiplicity_V0M > event_info.multiplicity_V0M);
	}

 private:
	int V0M_class;

};
#endif
