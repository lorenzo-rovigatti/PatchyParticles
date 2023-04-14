#include <stdio.h>
#include <stdlib.h>

#include "ordinator.h"


ordinator *getOrdinator(int maxelements)
{
	ordinator *o=malloc(sizeof(ordinator));
	
	o->keyvector=calloc(maxelements,sizeof(int));
	o->numkey=0;
	
	return o;
}

void freeOrdinator(ordinator *o)
{
	free(o->keyvector);
	free(o);
}

int ordinatorGetOrder(ordinator *o,int key)
{
	int i=0;
	
	while ((i<o->numkey) && (key!=o->keyvector[i]))
	{
		i++;
	}
	
	if (i==o->numkey)
	{
		o->keyvector[i]=key;
		o->numkey++;
	}
	
	return i;
}

int ordinatorGetKey(ordinator *o,int order)
{
	if (order<o->numkey)
	{
		return o->keyvector[order];
	}
	else
		return -1;
}

