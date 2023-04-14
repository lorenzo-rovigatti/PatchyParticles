#ifndef INTERACTION_MAP_H
#define INTERACTION_MAP_H

typedef struct _interactionmap {
	int num;		// number particles interacting
	int *who;		// particles interacting
	int **with;		// interacting partners
	int *howmany;	// number of interacting partners
	int num_bonds;
	int size;		// total number of particles
	double **rij2;
	jvector **rij;
} interactionmap;

interactionmap* createInteractionMap(int max_elements,int max_neighbours);
void freeInteractionMap(interactionmap *i);
void resetInteractionMap(interactionmap *i);
void buildImeFromIm(interactionmap *im,interactionmap *ime);

#endif
