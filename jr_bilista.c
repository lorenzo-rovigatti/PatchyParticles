#include <stdio.h>
#include <stdlib.h>
#include "jr_bilista.h"

bilista* bilistaGet(int size)
{
	bilista *b=calloc(1+size,sizeof(bilista));
	
	int i;
	
	b[0].prev=-1;
	b[0].next=-1;
	
	for (i=1;i<=size;i++)
	{
		b[i].prev=-2;
		b[i].next=-1;
	}
	
	b+=1;
	
	return b;
}

void bilistaFree(bilista *b)
{
	free(b-1);
}

bilista* bilistaCopy(bilista *b,int size)
{
	bilista *bcopy=calloc(1+size,sizeof(bilista));
	bcopy+=1;

	int i;
	for (i=-1;i<size;i++)
	{
		bcopy[i]=b[i];
	}
	return bcopy;	
}

void saveBilista(bilista *b,int size,FILE *pfile)
{
	fwrite(b-1,sizeof(bilista),1+size,pfile);
}

void getBilista(bilista *b,int size,FILE *pfile)
{
	size_t status;
	status=fread(b-1,sizeof(bilista),1+size,pfile);
}

int bilistaCompare(bilista *b1,bilista *b2,int size)
{
	int state=0;
	int i=-1;
	while (bilistaNext(b1,&i))
	{
		if (!bilistaIsIn(b2,i))
			state=1;
	}
	i=-1;
	while (bilistaNext(b2,&i))
	{
		if (!bilistaIsIn(b1,i))
			state=1;
	}


	return state;
}

int bilistaInsert(bilista *b,int el)
{
	if (b[el].prev!=-2)
	{
		// element already present
		return 0;
	}
	else
	{
		int tail=b[-1].prev;
		
		b[el].prev=tail;
		
		b[tail].next=el;
		
		b[-1].prev=el;
		
		return 1;
	}
}

int bilistaIsIn(bilista *b,int el)
{
	if (b[el].prev!=-2)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

int bilistaNext(bilista *b,int *el)
{
	/* l'uso classico in un ciclo potrebbe essere
	
	int i=-1;
	while (bilistaNext(b,&i))
		fai qualcosa con i
	
	*/
	
	*el=b[*el].next;
	
	if (*el==-1)
		return 0;
	else
		return 1;
	
}

int bilistaRemove(bilista *b,int el)
{
	if (b[el].prev==-2)
	{
		// element not present
		return 0;
	}
	else
	{
		
		int prev=b[el].prev;
		int next=b[el].next;
		
		b[prev].next=next;
		b[next].prev=prev;
		
		b[el].prev=-2;
		b[el].next=-1;
		
		return 1;
	}
}

void bilistaReset(bilista *b,int size)
{
	b-=1;
	
	int i;
	
	b[0].prev=-1;
	b[0].next=-1;
	
	for (i=1;i<=size;i++)
	{
		b[i].prev=-2;
		b[i].next=-1;
	}
	
	b+=1;
	
}


/*
void bilistaHeadInsert(bilista *b,int el)
{
	
}
*/

int bilistaPop(bilista *b,int *el)
{
	*el=b[-1].next;
	
	if (*el==-1)
		return 0;
	else
	{
		// remove
		
		int prev=b[*el].prev;
		int next=b[*el].next;
		
		b[prev].next=next;
		b[next].prev=prev;
		
		b[*el].prev=-2;
		b[*el].next=-1;
		
		return 1;
	}
}

/*
int main()
{
	bilista *b=bilistaGet(10);
	
	bilistaInsert(b,0);
	bilistaInsert(b,9);
	bilistaInsert(b,7);
	
	
	
// 	bilistaRemove(b,7);
//	bilistaRemove(b,0);
//	bilistaRemove(b,9);
	
//	if (bilistaIsIn(b,9))
//		printf("9 in bilista\n");
	
//	int i=-1;
//	while (bilistaNext(b,&i))
//		printf("%d\n",i);
	
	int position=1;
	int i=0;
	int selected=-1;
	
	// selezione saltando 9
// 	while ( (bilistaNext(b,&selected)) && (selected!=9 ? 1 : ++position) && (i++<position) );
// 	
// 	printf("Selected %d\n",selected);
	
	while (bilistaPop(b,&i))
	{
		printf("%d\n",i);
	}
	
	
	bilistaFree(b);
	return 0;
}
*/

