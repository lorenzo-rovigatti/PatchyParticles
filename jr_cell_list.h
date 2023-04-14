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


void celllistConstructor(FILE *config_file,vector *pos,int *ncolloids,vector *box,double *cutoff,
			 interactionmap *interactioList,double *cellsize,
			 listcell **list,void (**getInteractionMap)(args *argv,int argc),args **interactionMapArguments,int *numinteractionargs);

listcell* getList(vector box_sides,double cutoff,int num_particles);
void freeList(listcell *l);
void updateList(listcell *l,const vector *pos,int num);
void fullUpdateList(listcell *l,const vector *pos,int num,vector box_sides,double cutoff);
void changeCell(listcell *l,const vector *oldpos,const vector *newpos,int num);
void calculateInteractionMap(listcell *l,interactionmap *im);
int calculateSystemInteractionMap(listcell *l,interactionmap *im,vector *pos,vector *box,double cutoff);
int calculateSystemInteractionMapOrdered(listcell *l,interactionmap *im,vector *pos,vector *box,double cutoff);
void calculateInteractionMapWithCutoffDistance(listcell *l,interactionmap *im,vector *pos,vector *box,double cutoff);
void calculateExtendedInteractionMapWithCutoffDistance(listcell *l,interactionmap *im,interactionmap *ime,vector *pos,vector *box,double cutoff);
void calculateInteractionMapWithCutoffDistanceOrdered(listcell *l,interactionmap *ime,vector *pos,vector *box,double cutoff);
void addToList(listcell *l,const vector *pos,int num);
void get_perpendicular_versor(vector v, vector res);

#endif
