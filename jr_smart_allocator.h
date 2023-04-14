#ifndef SMART_ALLOCATOR_H
#define SMART_ALLOCATOR_H

/* these functions also initialize to zero the values of the arrays */

#define Vector1D(array,dim,type)                                               \
array=calloc(dim,sizeof(type));


#define Free1D(array)                                                          \
free(array);

#define Matrix2D(array,dim1,dim2,type)                                         \
do {                                                                           \
int j;                                                                         \
array=(type**)malloc((dim1)*sizeof(type*));                                    \
array[0]=(type*)calloc((dim1)*(dim2),sizeof(type));                            \
for (j=0;j<(dim1);j++)                                                         \
array[j]=array[0]+j*(dim2);                                                    \
} while (0);

#define Matrix2DSafe(array,dim1,dim2,buffer_space,type)                                         \
do {                                                                           \
int j;                                                                         \
array=(type**)malloc((dim1)*sizeof(type*));                                    \
array[0]=(type*)calloc((dim1)*(dim2)+buffer_space,sizeof(type));                            \
for (j=0;j<(dim1);j++)                                                         \
array[j]=array[0]+j*(dim2);                                                    \
} while (0);

/* the first element of array[0] is just a buffer */
#define UpperTriangularMatrix(array,dim,type)                                  \
do {                                                                           \
int j;                                                                         \
array=(type**)malloc(((dim)-1)*sizeof(type*));                                 \
array[0]=(type*)calloc(1+(((dim)*(dim-1))/2),sizeof(type));                    \
for (j=0;j<(dim)-1;j++)                                                        \
array[j]=array[0]+1+j*(dim)-((j+1)*(j+2))/2;                                   \
} while (0);


/* the first element of array[0] is just a buffer */
#define DiagonalUpperMatrix(array,dim,type)                                    \
do {                                                                           \
int j;                                                                         \
array=(type**)malloc(dim*sizeof(type*));                                       \
array[0]=(type*)calloc(((dim)*(dim+1))/2,sizeof(type));                        \
int cumulant=0;                                                                \
for (j=0;j<dim;j++)                                                            \
{                                                                              \
cumulant+=j;                                                                    \
array[j]=array[0]+j*dim-cumulant;                                              \
}                                                                              \
} while (0);


#define Free2D(array)                                                          \
free(array[0]);                                                                \
free(array);


#define Matrix3D(array,dim1,dim2,dim3,type)                                    \
do{                                                                            \
int j,k;                                                                       \
type *mem;                                                                     \
mem=(type*)calloc((dim1)*(dim2)*(dim3),sizeof(type));                          \
array=(type***)malloc((dim1)*sizeof(type**));                                  \
for (k=0;k<(dim1);k++)                                                         \
{                                                                              \
array[k]=(type**)malloc((dim2)*sizeof(type*));                                 \
for (j=0;j<(dim2);j++)                                                         \
{                                                                              \
array[k][j]=mem+k*(dim2)*(dim3)+j*(dim3);                                      \
}                                                                              \
}                                                                              \
} while (0);


#define Free3D(array,dim1)                                                     \
do{                                                                            \
int k;                                                                         \
free(array[0][0]);                                                             \
for (k=0;k<(dim1);k++)                                                         \
{                                                                              \
free(array[k]);                                                                \
}                                                                              \
free(array);                                                                   \
} while (0);



#define Matrix4D(array,dim1,dim2,dim3,dim4,type)                               \
do{                                                                            \
int i,j,k;                                                                     \
type *mem;                                                                     \
mem=(type*)calloc((dim1)*(dim2)*(dim3)*(dim4),sizeof(type));                   \
array=(type****)malloc((dim1)*sizeof(type***));                                \
for (k=0;k<(dim1);k++)                                                         \
{                                                                              \
array[k]=(type***)malloc((dim2)*sizeof(type**));                               \
for (j=0;j<(dim2);j++)                                                         \
{                                                                              \
array[k][j]=(type**)malloc((dim3)*sizeof(type*));                              \
for (i=0;i<(dim3);i++)                                                         \
{                                                                              \
array[k][j][i]=mem+k*(dim2)*(dim3)*(dim4)+j*(dim3)*(dim4)+i*(dim4);            \
}                                                                              \
}                                                                              \
}                                                                              \
} while (0);


#define Free4D(array,dim1,dim2)                                                \
do{                                                                            \
int k,j;                                                                       \
free(array[0][0][0]);                                                          \
for (k=0;k<(dim1);k++)                                                         \
{                                                                              \
for (j=0;j<(dim2);j++)                                                         \
{                                                                              \
free(array[k][j]);                                                             \
}                                                                              \
free(array[k]);                                                                \
}                                                                              \
free(array);                                                                   \
} while (0);


#endif

