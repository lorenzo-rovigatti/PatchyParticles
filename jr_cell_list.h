#ifndef CELL_LIST_H
#define CELL_LIST_H

typedef struct _listcell {
	int *HoC;                            // Head of Chain for linked list
	int *LinkedList;                     // linked list
	int NumberCells_x;                     // number of cells in one direction
	int NumberCells_y;
	int NumberCells_z;
	double CellSize_x;                     // cell edge size
	double CellSize_y;
	double CellSize_z;
} listcell;


void celllistConstructor(FILE *config_file,jvector *pos,int *ncolloids,jvector *box,double *cutoff,
			 interactionmap *interactioList,double *cellsize,
			 listcell **list,void (**getInteractionMap)(args *argv,int argc),args **interactionMapArguments,int *numinteractionargs);

listcell* getList(jvector box_sides,double cutoff,int num_particles);
void freeList(listcell *l);
void updateList(listcell *l,const jvector *pos,int num);
void fullUpdateList(listcell *l,const jvector *pos,int num,jvector box_sides,double cutoff);
void changeCell(listcell *l,const jvector *oldpos,const jvector *newpos,int num);
void calculateInteractionMap(listcell *l,interactionmap *im);
int calculateSystemInteractionMap(listcell *l,interactionmap *im,jvector *pos,jvector *box,double cutoff);
int calculateSystemInteractionMapOrdered(listcell *l,interactionmap *im,jvector *pos,jvector *box,double cutoff);
void calculateInteractionMapWithCutoffDistance(listcell *l,interactionmap *im,jvector *pos,jvector *box,double cutoff);
void calculateExtendedInteractionMapWithCutoffDistance(listcell *l,interactionmap *im,interactionmap *ime,jvector *pos,jvector *box,double cutoff);
void calculateInteractionMapWithCutoffDistanceOrdered(listcell *l,interactionmap *ime,jvector *pos,jvector *box,double cutoff);
void addToList(listcell *l,const jvector *pos,int num);
void get_perpendicular_versor(jvector v, jvector res);

#endif
