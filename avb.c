#include <math.h>
#include "avb.h"


avbmc *avbdata;

void AVBMC_init(input_file *input, System *syst, Output *IO)
{
	avbdata=malloc(sizeof(avbmc));
	
	int numpatches=syst->n_patches;
	
	avbdata->neighbours=malloc(numpatches*sizeof(PatchyParticle*));
	avbdata->num_neighbours=0;
	
	avbdata->avb_vin = syst->n_patches*syst->n_patches * (M_PI*(syst->kf_delta*syst->kf_delta*syst->kf_delta + 3.*SQR(syst->kf_delta) +
			3.*syst->kf_delta) * SQR(1. - syst->kf_cosmax)/3.);
	
	output_log_msg(IO, "Vavb = %lf\n", syst->avb_vin);
	
	syst->avb_vout = syst->V - syst->avb_vin;

	syst->avb_p = 0.5;
	getInputDouble(input, "avb_p", &syst->avb_p, 0);
	
}

void AVBMC_free()
{
	free(avbdata->neighbours);
	free(avbdata);
}
