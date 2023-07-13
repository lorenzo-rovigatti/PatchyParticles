#ifndef _SUS
#define _SUS

void susConstructor(input_file *input,int ncolloids);
void freeSus();
void updateHistograms(double op,double op0,double k,double temperature);
void saveResetHistograms(llint time);


#endif
