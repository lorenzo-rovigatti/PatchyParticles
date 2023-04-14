#ifndef BILISTA_H
#define BILISTA_H


typedef struct _bilista {
	int prev;
	int next;
} bilista;


bilista* bilistaGet(int size);

bilista* bilistaCopy(bilista *b,int size);

void saveBilista(bilista *b,int size,FILE *pfile);

void getBilista(bilista *b,int size,FILE *pfile);

int bilistaCompare(bilista *b1,bilista *b2,int size);

void bilistaFree(bilista *b);

int bilistaInsert(bilista *b,int el);

int bilistaIsIn(bilista *b,int el);

int bilistaNext(bilista *b,int *el);

int bilistaRemove(bilista *b,int el);

void bilistaReset(bilista *b,int size);

int bilistaPop(bilista *b,int *el);

#endif
