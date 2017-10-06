#ifndef AVB_H_
#define AVB_H_

typedef struct _avbmc {
	PatchyParticle **neighbours;
	int num_neighbours;
	
	double avb_vin, avb_vout;
	double avb_p;
} avbmc;


#endif