/*
 * neighs.h
 *
 *  Created on: 08/nov/2011
 *      Author: lorenzo
 */

#ifndef NEIGHS_H_
#define NEIGHS_H_

#include "defs.h"

int kf_would_interact(System *syst, PatchyParticle *p, vector r, matrix orient, int *onp, int *onq);
int kf_interact(System *syst, PatchyParticle *p, PatchyParticle *q, int *onp, int *onq);

#endif /* NEIGHS_H_ */
