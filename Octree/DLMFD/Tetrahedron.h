/** 
# Set of functions for a tetrahedron 
*/


/** Tests whether a point lies inside the tetrahedron */
//----------------------------------------------------------------------------
bool is_in_Tetrahedron( const double x, const double y, const double z, 
	const GeomParameter gp )
//----------------------------------------------------------------------------
{
  return false;
}




/** Computes the number of boundary points on the surface of the tetrahedron */
//----------------------------------------------------------------------------
void compute_nboundary_Tetrahedron( GeomParameter* gcp, int* nb, int* lN ) 
//----------------------------------------------------------------------------
{

}




/** Creates boundary points on the surface of the tetrahedron */
//----------------------------------------------------------------------------
void create_FD_Boundary_Tetrahedron( GeomParameter* gcp, 
	SolidBodyBoundary* dlm_bd, const int m, const int lN, vector pshift ) 
//----------------------------------------------------------------------------
{

}




/** Finds cells lying inside the tetrahedron */
//----------------------------------------------------------------------------
void create_FD_Interior_Tetrahedron( particle * p, vector Index_lambda, 
	vector shift ) 
//----------------------------------------------------------------------------
{

}




/** Reads geometric parameters of the tetrahedron */
//----------------------------------------------------------------------------
void update_Tetrahedron( GeomParameter* gcp ) 
//----------------------------------------------------------------------------
{  

}




/** Frees the geometric parameters of the tetrahedron */
//----------------------------------------------------------------------------
void free_Tetrahedron( GeomParameter* gcp ) 
//----------------------------------------------------------------------------
{  

}
