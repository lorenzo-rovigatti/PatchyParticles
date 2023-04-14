#ifndef ORDINATOR_H
#define ORDINATOR_H

typedef struct _ordinator {
	int *keyvector;
	int numkey;
} ordinator;

ordinator *getOrdinator(int maxelements);
void freeOrdinator(ordinator *o);
int ordinatorGetOrder(ordinator *o,int key);
int ordinatorGetKey(ordinator *o,int order);

#endif