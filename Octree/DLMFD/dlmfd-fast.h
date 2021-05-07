/** 
# Set of structures and functions to loop through the boundary point stencils 
and compute scalar products in a faster way.
*/

# ifndef DLMBLOCK 
#   define DLMBLOCK 5000
# endif

# if dimension == 2
#   define NDOFSTENCIL 9
# else
#   define NDOFSTENCIL 27
# endif


/** Structure and routines to store pointers and coefficients from the
initialization of the Uzawa algorithm and use them over the iterative process
to compute faster M_u^T*w = <DLM_w, v>_P(t) over the boundary points */
typedef struct {
  int n;
  int nm;
  double** qux;
  double** quy;
  double** dlmwx;
  double** dlmwy;
# if dimension == 3  
    double** quz;
    double** dlmwz;
# endif   
  double* weight; 
} BPFastLoop_LambdaMom;




//----------------------------------------------------------------------------
void initialize_and_allocate_BPFastLoop_LambdaMom( BPFastLoop_LambdaMom* loop, 
	const int nm ) 
//----------------------------------------------------------------------------
{
  loop->qux = (double**) calloc( nm, sizeof(double*) ); 
  loop->quy = (double**) calloc( nm, sizeof(double*) );  
  loop->dlmwx = (double**) calloc( nm, sizeof(double*) ); 
  loop->dlmwy = (double**) calloc( nm, sizeof(double*) ); 
# if dimension == 3  
    loop->quz = (double**) calloc( nm, sizeof(double*) );  
    loop->dlmwz = (double**) calloc( nm, sizeof(double*) );
# endif 
  loop->weight = (double*) calloc( nm, sizeof(double) );
  loop->n = 0;
  loop->nm = nm;
}




//----------------------------------------------------------------------------
void free_BPFastLoop_LambdaMom( BPFastLoop_LambdaMom* loop ) 
//----------------------------------------------------------------------------
{
  free(loop->qux); loop->qux = NULL;
  free(loop->quy); loop->quy = NULL;
  free(loop->dlmwx); loop->dlmwx = NULL;
  free(loop->dlmwy); loop->dlmwy = NULL;
# if dimension == 3  
    free(loop->quz); loop->quz = NULL;
    free(loop->dlmwz); loop->dlmwz = NULL;
# endif 
  free(loop->weight); loop->weight = NULL;
  loop->n = 0;
  loop->nm = 0;
}




//----------------------------------------------------------------------------
void append_BPFastLoop_LambdaMom( BPFastLoop_LambdaMom* loop,
	double* qux_, double* quy_,  double* quz_, 
	double* dlmwx_, double* dlmwy_, double* dlmwz_, double weight_ )
//----------------------------------------------------------------------------
{
  if ( loop->n >= loop->nm )
  {
    loop->nm += DLMBLOCK;
    loop->qux = (double**) realloc( loop->qux, loop->nm*sizeof(double*) ); 
    loop->quy = (double**) realloc( loop->quy, loop->nm*sizeof(double*) ); 
    loop->quz = (double**) realloc( loop->quz, loop->nm*sizeof(double*) );   
    loop->dlmwx = (double**) realloc( loop->dlmwx, loop->nm*sizeof(double*) ); 
    loop->dlmwy = (double**) realloc( loop->dlmwy, loop->nm*sizeof(double*) ); 
    loop->dlmwz = (double**) realloc( loop->dlmwz, loop->nm*sizeof(double*) );
    loop->weight = (double*) realloc( loop->weight, loop->nm*sizeof(double) );
  }

  int n = loop->n;
  loop->qux[n] = qux_;
  loop->quy[n] = quy_;  
  loop->quz[n] = quz_;  
  loop->dlmwx[n] = dlmwx_;
  loop->dlmwy[n] = dlmwy_;  
  loop->dlmwz[n] = dlmwz_;      
  loop->weight[n] = weight_;
  loop->n += 1;    
}




//----------------------------------------------------------------------------
void append_BPFastLoop_LambdaMom_2D( BPFastLoop_LambdaMom* loop,
	double* qux_, double* quy_, double* dlmwx_, double* dlmwy_, 
	double weight_ )
//----------------------------------------------------------------------------
{
  if ( loop->n >= loop->nm )
  {
    loop->nm += DLMBLOCK;
    loop->qux = (double**) realloc( loop->qux, loop->nm*sizeof(double*) ); 
    loop->quy = (double**) realloc( loop->quy, loop->nm*sizeof(double*) );  
    loop->dlmwx = (double**) realloc( loop->dlmwx, loop->nm*sizeof(double*) ); 
    loop->dlmwy = (double**) realloc( loop->dlmwy, loop->nm*sizeof(double*) ); 
    loop->weight = (double*) realloc( loop->weight, loop->nm*sizeof(double) );
  }

  int n = loop->n;
  loop->qux[n] = qux_;
  loop->quy[n] = quy_;  
  loop->dlmwx[n] = dlmwx_;
  loop->dlmwy[n] = dlmwy_;       
  loop->weight[n] = weight_;
  loop->n += 1;    
}








/** Structure and routines to store pointers and coefficients from the
initialization of the Uzawa algorithm and use them over the iterative process
to compute faster M_u*tu = <alpha, tu>_P(t) over the boundary points */
typedef struct {
  int nv;
  int nvm;
  double** dlmvx;
  double** dlmvy;
  int* ndof;
  int ntu;
  int ntum;      
  double** tux;
  double** tuy;
# if dimension == 3 
    double** dlmvz; 
    double** tuz;
# endif 
  double* weight; 
} BPFastLoop_ResU;




//----------------------------------------------------------------------------
void initialize_and_allocate_BPFastLoop_ResU( BPFastLoop_ResU* loop, 
	const int nm ) 
//----------------------------------------------------------------------------
{
  loop->dlmvx = (double**) calloc( nm, sizeof(double*) ); 
  loop->dlmvy = (double**) calloc( nm, sizeof(double*) ); 
  loop->ndof = (int*) calloc( nm, sizeof(int) );
  loop->nv = 0;
  loop->nvm = nm;  
       
  loop->ntum = NDOFSTENCIL*nm;
  loop->tux = (double**) calloc( loop->ntum, sizeof(double*) ); 
  loop->tuy = (double**) calloc( loop->ntum, sizeof(double*) ); 

# if dimension == 3 
    loop->dlmvz = (double**) calloc( nm, sizeof(double*) );
    loop->tuz = (double**) calloc( loop->ntum, sizeof(double*) );
# endif 
    
  loop->weight = (double*) calloc( loop->ntum, sizeof(double) );
  loop->ntu = 0;  
}




//----------------------------------------------------------------------------
void free_BPFastLoop_ResU( BPFastLoop_ResU* loop ) 
//----------------------------------------------------------------------------
{
  free(loop->tux); loop->tux = NULL;
  free(loop->tuy); loop->tuy = NULL;
  free(loop->dlmvx); loop->dlmvx = NULL;
  free(loop->dlmvy); loop->dlmvy = NULL;
# if dimension == 3 
    free(loop->tuz); loop->tuz = NULL;
    free(loop->dlmvz); loop->dlmvz = NULL;
# endif 
  free(loop->weight); loop->weight = NULL;
  free(loop->ndof); loop->ndof = NULL;  
  loop->nv = 0;
  loop->nvm = 0;
  loop->ntu = 0;
  loop->ntum = 0;  
}




//----------------------------------------------------------------------------
void append_BPFastLoop_ResU_tu( BPFastLoop_ResU* loop,
	double* tux_, double* tuy_,  double* tuz_, double weight_ )
//----------------------------------------------------------------------------
{
  if ( loop->ntu >= loop->ntum )
  {
    loop->ntum += NDOFSTENCIL*DLMBLOCK;
    loop->tux = (double**) realloc( loop->tux, loop->ntum*sizeof(double*) ); 
    loop->tuy = (double**) realloc( loop->tuy, loop->ntum*sizeof(double*) ); 
    loop->tuz = (double**) realloc( loop->tuz, loop->ntum*sizeof(double*) );   
    loop->weight = (double*) realloc( loop->weight, loop->ntum*sizeof(double) );
  }

  int n = loop->ntu;
  loop->tux[n] = tux_;
  loop->tuy[n] = tuy_;  
  loop->tuz[n] = tuz_;       
  loop->weight[n] = weight_;
  loop->ntu += 1;    
}




//----------------------------------------------------------------------------
void append_BPFastLoop_ResU_v( BPFastLoop_ResU* loop,
	double* dlmvx_, double* dlmvy_,  double* dlmvz_, int ndof_ )
//----------------------------------------------------------------------------
{
  if ( loop->nv >= loop->nvm )
  {
    loop->nvm += DLMBLOCK;
    loop->dlmvx = (double**) realloc( loop->dlmvx, loop->nvm*sizeof(double*) ); 
    loop->dlmvy = (double**) realloc( loop->dlmvy, loop->nvm*sizeof(double*) ); 
    loop->dlmvz = (double**) realloc( loop->dlmvz, loop->nvm*sizeof(double*) );
    loop->ndof = (int*) realloc( loop->ndof, loop->nvm*sizeof(int) );
  }

  int n = loop->nv;
  loop->dlmvx[n] = dlmvx_;
  loop->dlmvy[n] = dlmvy_;  
  loop->dlmvz[n] = dlmvz_;       
  loop->ndof[n] = ndof_;
  loop->nv += 1;    
}




//----------------------------------------------------------------------------
void append_BPFastLoop_ResU_tu_2D( BPFastLoop_ResU* loop,
	double* tux_, double* tuy_,  double weight_ )
//----------------------------------------------------------------------------
{
  if ( loop->ntu >= loop->ntum )
  {
    loop->ntum += NDOFSTENCIL*DLMBLOCK;
    loop->tux = (double**) realloc( loop->tux, loop->ntum*sizeof(double*) ); 
    loop->tuy = (double**) realloc( loop->tuy, loop->ntum*sizeof(double*) );   
    loop->weight = (double*) realloc( loop->weight, loop->ntum*sizeof(double) );
  }

  int n = loop->ntu;
  loop->tux[n] = tux_;
  loop->tuy[n] = tuy_;      
  loop->weight[n] = weight_;
  loop->ntu += 1;    
}




//----------------------------------------------------------------------------
void append_BPFastLoop_ResU_v_2D( BPFastLoop_ResU* loop,
	double* dlmvx_, double* dlmvy_,  double* dlmvz_, int ndof_ )
//----------------------------------------------------------------------------
{
  if ( loop->nv >= loop->nvm )
  {
    loop->nvm += DLMBLOCK;
    loop->dlmvx = (double**) realloc( loop->dlmvx, loop->nvm*sizeof(double*) );
    loop->dlmvy = (double**) realloc( loop->dlmvy, loop->nvm*sizeof(double*) );
    loop->ndof = (int*) realloc( loop->ndof, loop->nvm*sizeof(int) );
  }

  int n = loop->nv;
  loop->dlmvx[n] = dlmvx_;
  loop->dlmvy[n] = dlmvy_;  
  loop->ndof[n] = ndof_;
  loop->nv += 1;    
}
