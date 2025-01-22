/** 
# Set of functions that computes the weights associated to a Lagrange
  multiplier
*/

/** 
## Some definitions
*/
# define sup2deltax ( relnl.x < - 1.9 * delta ) && ( relnl.x > - 2.1 * delta )
# define sup1deltax ( relnl.x < - 0.9 * delta ) && ( relnl.x > - 1.1 * delta )
# define sup0deltax fabs(relnl.x) < 0.1 * delta

# define sup2deltay ( relnl.y < - 1.9 * delta ) && ( relnl.y > - 2.1 * delta )
# define sup1deltay ( relnl.y < - 0.9 * delta ) && ( relnl.y > - 1.1 * delta )
# define sup0deltay  fabs(relnl.y) < 0.1 * delta 

# define sup2deltaz ( relnl.z < - 1.9 * delta ) && ( relnl.z > - 2.1 * delta )
# define sup1deltaz ( relnl.z < - 0.9 * delta ) && ( relnl.z > - 1.1 * delta )
# define sup0deltaz fabs(relnl.z) < 0.1 * delta

# define inf2deltax ( relnl.x < 2.1 * delta ) && ( relnl.x > 1.9 * delta )
# define inf1deltax ( relnl.x < 1.1 * delta ) && ( relnl.x > 0.9 * delta )
  
# define inf2deltay ( relnl.y < 2.1 * delta ) && ( relnl.y > 1.9 * delta )
# define inf1deltay ( relnl.y < 1.1 * delta ) && ( relnl.y > 0.9 * delta )

# define inf2deltaz ( relnl.z < 2.1 * delta ) && ( relnl.z > 1.9 * delta )
# define inf1deltaz ( relnl.z < 1.1 * delta ) && ( relnl.z > 0.9 * delta )

/** Compute the 1D distance for stencils in Periodic conditions */
#define ACROSS_PERIODIC( a, b ) ( a - b > L0/2. || a - b < - L0/2. )
#define PERIODIC_1DIST( a, b ) ( ( a - L0 - b > L0/2. \
	|| a - L0 - b < - L0/2. ) ? a + L0 - b : a - L0 - b )
#define GENERAL_1DIST( a, b ) ( ACROSS_PERIODIC( a, b ) ? \
	PERIODIC_1DIST( a, b ) : a - b )



//----------------------------------------------------------------------------
void assign_dial( coord rel, int* CX ) 
//----------------------------------------------------------------------------
{
# if dimension == 2
    if ( rel.x > 0. )
    {
      if ( rel.y > 0. )
        /* We are on the lower left cell */
        *CX = 1;
      else // rel.y <= 0.
        /* We are on the upper left cell */
        *CX = 4;
    }  
    else // rel.x <= 0.
    {
      if ( rel.y > 0. )
        /* we are on the lower right cell */
        *CX = 2;
      else // rel.y <= 0. ) {
        /* We are on the upper right cell */
        *CX = 3;
    }

# elif dimension == 3 
    if ( rel.z > 0. )
    {
      /* one has to pick the forward cells + the middle ones according 
      to their 2D counterpart */
      if ( rel.x > 0. )
      {
        if ( rel.y > 0. )
	  /* We are on the lower left cell */
	  *CX = 10;
        else // rel.y <= 0. )
	  /* We are on the upper left cell */
	  *CX = 40;
      }
      else // rel.x <= 0. 
      {
        if ( rel.y > 0. )
	  /* we are on the lower right cell */
	  *CX = 20;
        else // rel.y <= 0 )
	  /* We are on the upper right cell */
	  *CX = 30;
      }
    }
    else // rel.z <= 0. 
    {
      /* one has to pick the backward cells + the middle ones according 
      to their 2D counterpart */
      if ( rel.x > 0. )
      {
        if ( rel.y > 0. )
	  /* We are on the lower left cell */
	  *CX = 11;
        else // rel.y <= 0. 
	  /* We are on the upper left cell */
	  *CX = 41;
      }
      else // rel.x <= 0. )
      {
        if ( rel.y > 0. )
	  /* We are on the lower right cell */
	  *CX = 21;
        else // rel.y <= 0. ) 
	  /* We are on the upper right cell */
	  *CX = 31;
      }
    }
# endif 
}




//----------------------------------------------------------------------------
void NCX1_CX1( const coord relnl, const double delta, 
	const int weight_numbers[9], int* weight_id, size_t* goflag ) 
//----------------------------------------------------------------------------
{
  // weight_numbers[2] -> right column top
  // weight_numbers[5] -> right column middle
  // weight_numbers[1] -> right column bottom

  // weight_numbers[6] -> middle column top
  // weight_numbers[8] -> middle colum middle
  // weight_numbers[4] -> middle column bottom

  // weight_numbers[3] -> left column top
  // weight_numbers[7] -> left column middle
  // weight_numbers[0] -> left column bottom

  /* Stencil config, O is the cell containing the Lagrange multiplier */
  /* |_|_|x|x|x| --> |_|_|3|6|2| */
  /* |_|_|x|x|x| --> |_|_|7|8|5| */
  /* |_|_|O|x|x| --> |_|_|0|4|1| */ 
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  
  // Right column, sup2deltax
  if ( relnl.x < - 1.9 * delta && relnl.x > - 2.1 * delta ) 
  {     
    // top row, sup2deltay
    if ( relnl.y < - 1.9 * delta && relnl.y > - 2.1 * delta ) 
    {
      *weight_id = weight_numbers[2]; 
      *goflag = 1;
    }
    // middle row, sup1deltay
    else if ( relnl.y < - 0.9 * delta && relnl.y > - 1.1 * delta ) 
    {
      *weight_id = weight_numbers[5]; 
      *goflag = 1;
    }
    // bottom row, sup0deltay
    else if ( fabs(relnl.y) < 0.1 * delta ) 
    {
      *weight_id = weight_numbers[1]; 
      *goflag = 1;
    }
  }
  
  // Middle column, sup1deltax
  else if ( relnl.x < - 0.9 * delta && relnl.x > - 1.1 * delta ) 
  {
    // top row, sup2deltay
    if ( relnl.y < - 1.9 * delta && relnl.y > - 2.1 * delta ) 
    {
      *weight_id = weight_numbers[6]; 
      *goflag = 1;
    }
    // middle row, sup1deltay
    else if ( relnl.y < - 0.9 * delta && relnl.y > - 1.1 * delta ) 
    {
      *weight_id = weight_numbers[8]; 
      *goflag = 1;
    }
    // bottom row, sup0deltay
    else if ( fabs(relnl.y) < 0.1 * delta ) 
    {
      *weight_id = weight_numbers[4]; 
      *goflag = 1;
    }
  }
   
  // Left column, sup0deltax
  else if ( fabs(relnl.x) < 0.1 * delta ) 
  {    
    // top row, sup2deltay
    if ( relnl.y < - 1.9 * delta && relnl.y > - 2.1 * delta ) 
    {
      *weight_id = weight_numbers[3]; 
      *goflag = 1;
    }
    // middle row, sup1deltay
    if ( relnl.y < - 0.9 * delta && relnl.y > - 1.1 * delta ) 
    {
      *weight_id = weight_numbers[7]; 
      *goflag = 1;
    }
    // bottom row, sup0deltay
    if ( fabs(relnl.y) < 0.1 * delta ) 
    {
      *weight_id = weight_numbers[0]; 
      *goflag = 1;
    }
  }  
}




//----------------------------------------------------------------------------
void NCX1_CX2( const coord relnl, const double delta, 
	const int weight_numbers[9], int* weight_id, size_t* goflag ) 
//----------------------------------------------------------------------------
{
  // weight_numbers[2] -> right column top
  // weight_numbers[5] -> right column middle
  // weight_numbers[1] -> right column bottom

  // weight_numbers[6] -> middle column top
  // weight_numbers[8] -> middle colum middle
  // weight_numbers[4] -> middle column bottom

  // weight_numbers[3] -> left column top
  // weight_numbers[7] -> left column middle
  // weight_numbers[0] -> left column bottom

  /* Stencil config, O is the cell containing the Lagrange multiplier */
  /* |_|x|x|x|_| -->  |_|3|6|2|_| */
  /* |_|x|x|x|_| -->  |_|7|8|5|_| */
  /* |_|x|O|x|_| -->  |_|0|4|1|_| */ 
  /* |_|_|_|_|_| -->  |_|_|_|_|_| */
  /* |_|_|_|_|_| -->  |_|_|_|_|_| */

  // Right column, sup1deltax
  if ( relnl.x < - 0.9 * delta && relnl.x > - 1.1 * delta ) 
  {    
    // top row, sup2deltay
    if ( relnl.y < - 1.9 * delta && relnl.y > - 2.1 * delta ) 
    { 
      *weight_id = weight_numbers[2]; 
      *goflag = 1;
    }
    // middle row, sup1deltay
    else if ( relnl.y < - 0.9 * delta && relnl.y > - 1.1 * delta ) 
    { 
      *weight_id = weight_numbers[5]; 
      *goflag = 1;
    }
    // bottom row, sup0deltay
    else if ( fabs(relnl.y) < 0.1 * delta ) 
    {
      *weight_id = weight_numbers[1]; 
      *goflag = 1;
    }
  }

  // Middle column, sup0deltax
  else if ( fabs(relnl.x) < 0.1 * delta ) 
  {
    // top row, sup2deltay
    if ( relnl.y < - 1.9 * delta && relnl.y > - 2.1 * delta ) 
    {
      *weight_id = weight_numbers[6]; 
      *goflag = 1;
    }
    // middle row, sup1deltay
    else if ( relnl.y < - 0.9 * delta && relnl.y > - 1.1 * delta ) 
    {
      *weight_id = weight_numbers[8]; 
      *goflag = 1;
    }
    // bottom row, sup0deltay
    else if ( fabs(relnl.y) < 0.1 * delta ) 
    {
      *weight_id = weight_numbers[4]; 
      *goflag = 1;
    }
  }

  // Left column
  else if ( relnl.x < 1.1 * delta && relnl.x > 0.9 * delta ) 
  {    
    // top row, sup2deltay
    if ( relnl.y < - 1.9 * delta && relnl.y > - 2.1 * delta ) 
    {
      *weight_id = weight_numbers[3]; 
      *goflag = 1;
    }
    // middle row, sup1deltay
    else if ( relnl.y < - 0.9 * delta && relnl.y > - 1.1 * delta ) 
    {
      *weight_id = weight_numbers[7]; 
      *goflag = 1;
    }
    // bottom row, sup0deltay
    else if ( fabs(relnl.y) < 0.1 * delta ) 
    {
      *weight_id = weight_numbers[0]; 
      *goflag = 1;
    }
  }
}




//----------------------------------------------------------------------------
void NCX_centred( const coord relnl, const double delta, 
	const int weight_numbers[9], int* weight_id, size_t* goflag ) 
//----------------------------------------------------------------------------
{
  /* Stencil config, O is the cell containing the Lagrange multiplier */
  /* |_|_|_|_|_| -->  |_|_|_|_|_| */
  /* |_|x|x|x|_| -->  |_|3|6|2|_| */
  /* |_|x|O|x|_| -->  |_|7|8|5|_| */ 
  /* |_|x|x|x|_| -->  |_|0|4|1|_| */
  /* |_|_|_|_|_| -->  |_|_|_|_|_| */

  
  // Right column, sup1deltax
  if (relnl.x < - 0.9 * delta && relnl.x > - 1.1 * delta) 
  {    
    // top row, sup1deltay
    if ( relnl.y < - 0.9 * delta && relnl.y > - 1.1 * delta ) 
    { 
      *weight_id = weight_numbers[2]; 
      *goflag = 1;
    }
    // middle row, sup0deltay
    else if ( fabs(relnl.y) <  0.1 * delta ) 
    {
      *weight_id = weight_numbers[5]; 
      *goflag = 1;
    }
    // bottom row, inf1deltay
    else if ( relnl.y < 1.1 * delta && relnl.y > 0.9 * delta ) 
    {
      *weight_id = weight_numbers[1]; 
      *goflag = 1;
    }
  }
  
  // Middle column, sup0deltax
  else if ( fabs(relnl.x) < 0.1 * delta ) 
  {    
    // top row, sup1deltay
    if ( relnl.y < - 0.9 * delta && relnl.y > - 1.1 * delta ) 
    {
      *weight_id = weight_numbers[6]; 
      *goflag = 1;
    }
    // middle row, sup0deltay
    else if ( fabs(relnl.y) <  0.1 * delta ) 
    {
      *weight_id = weight_numbers[8]; 
      *goflag = 1;
    }
    // bottom row, inf1deltay
    else if ( relnl.y < 1.1 * delta && relnl.y > 0.9 * delta ) 
    { 
      *weight_id = weight_numbers[4]; 
      *goflag = 1;
    }
  }
  
  // Left column, inf1deltax
  else if ( relnl.x < 1.1 * delta && relnl.x > 0.9 * delta ) 
  {
    // top row, sup1deltay
    if ( relnl.y < - 0.9 * delta && relnl.y > - 1.1 * delta ) 
    {
      *weight_id = weight_numbers[3]; 
      *goflag = 1;
    }
    // middle row, sup0deltay
    else if ( fabs(relnl.y) <  0.1 * delta ) 
    {
      *weight_id = weight_numbers[7]; 
      *goflag = 1; 
    }
    // bottom row, inf1deltay
    else if ( relnl.y < 1.1 * delta && relnl.y > 0.9 * delta ) 
    {
      *weight_id = weight_numbers[0]; 
      *goflag = 1;
    }
  }
}




//----------------------------------------------------------------------------
void NCX1_CX4( const coord relnl, const double delta, 
	const int weight_numbers[9], int* weight_id, size_t* goflag ) 
//----------------------------------------------------------------------------
{
  // weight_numbers[2] -> right column top
  // weight_numbers[5] -> right column middle
  // weight_numbers[1] -> right column bottom

  // weight_numbers[6] -> middle column top
  // weight_numbers[8] -> middle colum middle
  // weight_numbers[4] -> middle column bottom

  // weight_numbers[3] -> left column top
  // weight_numbers[7] -> left column middle
  // weight_numbers[0] -> left column bottom


  /* Stencil config, O is the cell containing the Lagrange multiplier */
  /* |_|_|_|_|_| -->  |_|_|_|_|_| */
  /* |_|_|x|x|x| -->  |_|_|3|6|2| */
  /* |_|_|O|x|x| -->  |_|_|7|8|5| */ 
  /* |_|_|x|x|x| -->  |_|_|0|4|1| */
  /* |_|_|_|_|_| -->  |_|_|_|_|_| */

  
  // Right column, sup2deltax
  if ( relnl.x < -1.9 * delta  &&  relnl.x > - 2.1 * delta ) 
  {
    // top row, sup1deltay
    if ( relnl.y < - 0.9 * delta && relnl.y > - 1.1 * delta ) 
    {
      *weight_id = weight_numbers[2]; 
      *goflag = 1;
    }
    // middle row, sup0deltay
    else if ( fabs(relnl.y) <  0.1 * delta ) 
    {
      *weight_id = weight_numbers[5]; 
      *goflag = 1;
    }
    // bottom row, inf1deltay
    else if ( relnl.y < 1.1 * delta && relnl.y > 0.9 * delta ) 
    {
      *weight_id = weight_numbers[1]; 
      *goflag = 1;
    }
  }
  
  // Middle column, sup1deltax
  else if ( relnl.x < - 0.9 * delta && relnl.x > - 1.1 * delta ) 
  {
    // top row, sup1deltay
    if ( relnl.y < - 0.9 * delta && relnl.y > - 1.1 * delta ) 
    {
      *weight_id = weight_numbers[6]; 
      *goflag = 1;
    }
    // middle row, sup0deltay
    else if ( fabs(relnl.y) <  0.1 * delta ) 
    {
      *weight_id = weight_numbers[8]; 
      *goflag = 1;
    }
    // bottom row, inf1deltay
    else if ( relnl.y < 1.1 * delta && relnl.y > 0.9 * delta ) 
    { 
      *weight_id = weight_numbers[4]; 
      *goflag = 1;
    }
  }

  // Left column, sup0deltax
  else if ( fabs(relnl.x) < 0.1 * delta ) 
  {
    // top row, sup1deltay
    if ( relnl.y < - 0.9 * delta && relnl.y > - 1.1 * delta ) 
    {
      *weight_id = weight_numbers[3]; 
      *goflag = 1;
    }
    // middle row, sup0deltay
    else if( fabs(relnl.y) <  0.1 * delta ) 
    {
      *weight_id = weight_numbers[7]; 
      *goflag = 1;
    }
    // bottom row, inf1deltay
    else if( relnl.y < 1.1 * delta && relnl.y > 0.9 * delta ) 
    {
      *weight_id = weight_numbers[0]; 
      *goflag = 1;
    }
  }
}




//----------------------------------------------------------------------------
void NCX2_CX2( const coord relnl, const double delta, 
	const int weight_numbers[9], int* weight_id, size_t* goflag ) 
//----------------------------------------------------------------------------
{
  // weight_numbers[2] -> right column top
  // weight_numbers[5] -> right column middle
  // weight_numbers[1] -> right column bottom

  // weight_numbers[6] -> middle column top
  // weight_numbers[8] -> middle colum middle
  // weight_numbers[4] -> middle column bottom

  // weight_numbers[3] -> left column top
  // weight_numbers[7] -> left column middle
  // weight_numbers[0] -> left column bottom


  /* Stencil config, O is the cell containing the Lagrange multiplier */
  /* |x|x|x|_|_| -->  |3|6|2|_|_| */
  /* |x|x|x|_|_| -->  |7|8|5|_|_| */
  /* |x|x|O|_|_| -->  |0|4|1|_|_| */ 
  /* |_|_|_|_|_| -->  |_|_|_|_|_| */
  /* |_|_|_|_|_| -->  |_|_|_|_|_| */

  
  // Right column, sup0deltax
  if ( fabs(relnl.x) < 0.1 * delta ) 
  {
    // top row, sup2deltay
    if ( relnl.y < - 1.9 * delta && relnl.y > - 2.1 * delta ) 
    { 
      *weight_id = weight_numbers[2]; 
      *goflag = 1;
    }
    // middle row, sup1deltay
    else if ( relnl.y < - 0.9 * delta && relnl.y > - 1.1 * delta )
     {
      *weight_id = weight_numbers[5]; 
      *goflag = 1;
    }
    // bottom row, sup0deltay
    else if ( fabs(relnl.y) < 0.1 * delta ) 
    {
      *weight_id = weight_numbers[1]; 
      *goflag = 1;
    }
  }
  
  // Middle column, inf1deltax 
  else if ( relnl.x < 1.1 * delta && relnl.x > 0.9 * delta ) 
  {
    // top row, sup2deltay
    if ( relnl.y < - 1.9 * delta && relnl.y > - 2.1 * delta ) 
    {
      *weight_id = weight_numbers[6]; 
      *goflag = 1;
    }
    // middle row, sup1deltay
    else if ( relnl.y < - 0.9 * delta && relnl.y > - 1.1 * delta ) 
    {
      *weight_id = weight_numbers[8]; 
      *goflag = 1;
    }
    // bottom row, sup0deltay
    else if ( fabs(relnl.y) < 0.1 * delta ) 
    {
      *weight_id = weight_numbers[4]; 
      *goflag = 1;
    }
  }
  
  // Left column, inf2deltax
  else if ( relnl.x < 2.1 * delta  &&  relnl.x > 1.9 * delta ) 
  {
    // top row, sup2deltay
    if ( relnl.y < - 1.9 * delta && relnl.y > - 2.1 * delta ) 
    {
      *weight_id = weight_numbers[3]; 
      *goflag = 1;
    }
    // middle row, sup1deltay
    else if ( relnl.y < - 0.9 * delta && relnl.y > - 1.1 * delta ) 
    {
      *weight_id = weight_numbers[7]; 
      *goflag = 1;
    }
    // bottom row, sup0deltay
    else if ( fabs(relnl.y) < 0.1 * delta ) 
    {
      *weight_id = weight_numbers[0]; 
      *goflag = 1;
    }
  }
}




//----------------------------------------------------------------------------
void NCX2_CX3( const coord relnl, const double delta, 
	const int weight_numbers[9], int* weight_id, size_t* goflag ) 
//----------------------------------------------------------------------------
{
  // weight_numbers[2] -> right column top
  // weight_numbers[5] -> right column middle
  // weight_numbers[1] -> right column bottom

  // weight_numbers[6] -> middle column top
  // weight_numbers[8] -> middle colum middle
  // weight_numbers[4] -> middle column bottom

  // weight_numbers[3] -> left column top
  // weight_numbers[7] -> left column middle
  // weight_numbers[0] -> left column bottom

  /* Stencil config, O is the cell containing the Lagrange multiplier */
  /* |_|_|_|_|_| -->  |_|_|_|_|_| */
  /* |x|x|x|_|_| -->  |3|6|2|_|_| */
  /* |x|x|O|_|_| -->  |7|8|5|_|_| */ 
  /* |x|x|x|_|_| -->  |0|4|1|_|_| */
  /* |_|_|_|_|_| -->  |_|_|_|_|_| */

  
  // Right column, sup0deltax
  if ( fabs(relnl.x) < 0.1 * delta ) 
  {
    // top row, sup1deltay
    if ( relnl.y < - 0.9 * delta && relnl.y > - 1.1 * delta ) 
    {
      *weight_id = weight_numbers[2]; 
      *goflag = 1;
    }
    // middle row, sup0deltay
    else if ( fabs(relnl.y) <  0.1 * delta ) 
    {
      *weight_id = weight_numbers[5]; 
      *goflag = 1;
    }
    // bottom row, inf1deltay
    else if ( relnl.y < 1.1 * delta && relnl.y > 0.9 * delta ) 
    {
      *weight_id = weight_numbers[1]; 
      *goflag = 1;
    }
  }
  
  // Middle column, inf1deltax
  else if ( relnl.x < 1.1 * delta && relnl.x > 0.9 * delta ) 
  {
    // top row, sup1deltay
    if ( relnl.y < - 0.9 * delta && relnl.y > - 1.1 * delta ) 
    {
      *weight_id = weight_numbers[6]; 
      *goflag = 1;
    }
    // middle row, sup0deltay
    else if ( fabs(relnl.y) <  0.1 * delta ) 
    {
      *weight_id = weight_numbers[8];
       *goflag = 1;
    }
    // bottom row, inf1deltay
    else if ( relnl.y < 1.1 * delta && relnl.y > 0.9 * delta ) 
    {
      *weight_id = weight_numbers[4]; 
      *goflag = 1;
    }
  }
  
  // Left column, inf2deltax
  else if ( relnl.x < 2.1 * delta  &&  relnl.x > 1.9 * delta ) 
  {
    // top row, sup1deltay
    if ( relnl.y < - 0.9 * delta && relnl.y > - 1.1 * delta ) 
    {
      *weight_id = weight_numbers[3]; 
      *goflag = 1;
    }
    // middle row, sup0deltay
    else if ( fabs(relnl.y) <  0.1 * delta ) 
    {
      *weight_id = weight_numbers[7]; 
      *goflag = 1;
    }
    // bottom row, inf1deltay
    else if ( relnl.y < 1.1 * delta && relnl.y > 0.9 * delta ) 
    {
      *weight_id = weight_numbers[0]; 
      *goflag = 1;
    }	
  }
}




//----------------------------------------------------------------------------
void NCX3_CX3( const coord relnl, const double delta, 
	const int weight_numbers[9], int* weight_id, size_t* goflag ) 
//----------------------------------------------------------------------------
{
  // weight_numbers[2] -> right column top
  // weight_numbers[5] -> right column middle
  // weight_numbers[1] -> right column bottom

  // weight_numbers[6] -> middle column top
  // weight_numbers[8] -> middle colum middle
  // weight_numbers[4] -> middle column bottom

  // weight_numbers[3] -> left column top
  // weight_numbers[7] -> left column middle
  // weight_numbers[0] -> left column bottom

  /* Stencil config, O is the cell containing the Lagrange multiplier */
  /* |_|_|_|_|_| -->  |_|_|_|_|_| */
  /* |_|_|_|_|_| -->  |_|_|_|_|_| */
  /* |x|x|O|_|_| -->  |3|6|2|_|_| */ 
  /* |x|x|x|_|_| -->  |7|8|5|_|_| */
  /* |x|x|x|_|_| -->  |0|4|1|_|_| */
  
  // Right column, sup0deltax
  if ( fabs(relnl.x) < 0.1 * delta ) 
  {      
    // top row, sup0deltay
    if ( fabs(relnl.y) < 0.1 * delta ) 
    {
      *weight_id = weight_numbers[2]; 
      *goflag = 1;
    }
    // middle row, inf1deltay
    else if ( relnl.y < 1.1 * delta && relnl.y > 0.9 * delta ) 
    {
      *weight_id = weight_numbers[5]; 
      *goflag = 1; 
    }
    // bottom row, inf2deltay
    else if ( relnl.y < 2.1 * delta  &&  relnl.y > 1.9 * delta ) 
    {
      *weight_id = weight_numbers[1]; 
      *goflag = 1;
    }
  }
  
  // Middle column, inf1deltax
  else if ( relnl.x < 1.1 * delta && relnl.x > 0.9 * delta ) 
  {
    // top row, sup0deltay
    if ( fabs(relnl.y) < 0.1 * delta ) 
    {
      *weight_id = weight_numbers[6]; 
      *goflag = 1;
    }
    // middle row, inf1deltay
    else if ( relnl.y < 1.1 * delta && relnl.y > 0.9 * delta ) 
    {
      *weight_id = weight_numbers[8]; 
      *goflag = 1; 
    }
    // bottom row, inf2deltay
    else if ( relnl.y < 2.1 * delta  &&  relnl.y > 1.9 * delta ) 
    {
      *weight_id = weight_numbers[4]; 
      *goflag = 1;
    }
  }
  
  // Left column, inf2deltax
  else if ( relnl.x < 2.1 * delta  &&  relnl.x > 1.9 * delta ) 
  {	
    // top row : sup0deltay
    if ( fabs(relnl.y) < 0.1 * delta ) 
    {
      *weight_id = weight_numbers[3]; 
      *goflag = 1;
    }
    // middle row, inf1deltay
    else if ( relnl.y < 1.1 * delta && relnl.y > 0.9 * delta ) 
    {
      *weight_id = weight_numbers[7]; 
      *goflag = 1; 
    }
    // bottom row, inf2deltay
    else if ( relnl.y < 2.1 * delta  &&  relnl.y > 1.9 * delta ) 
    {
      *weight_id = weight_numbers[0]; 
      *goflag = 1;
    }
  }

}




//----------------------------------------------------------------------------
void NCX3_CX4( const coord relnl, const double delta, 
	const int weight_numbers[9], int* weight_id, size_t* goflag ) 
//----------------------------------------------------------------------------
{
  // weight_numbers[2] -> right column top
  // weight_numbers[5] -> right column middle
  // weight_numbers[1] -> right column bottom

  // weight_numbers[6] -> middle column top
  // weight_numbers[8] -> middle colum middle
  // weight_numbers[4] -> middle column bottom

  // weight_numbers[3] -> left column top
  // weight_numbers[7] -> left column middle
  // weight_numbers[0] -> left column bottom

  /* Stencil config, O is the cell containing the Lagrange multiplier */
  /* |_|_|_|_|_| -->  |_|_|_|_|_| */
  /* |_|_|_|_|_| -->  |_|_|_|_|_| */
  /* |_|x|O|x|_| -->  |_|3|6|2|_| */ 
  /* |_|x|x|x|_| -->  |_|7|8|5|_| */
  /* |_|x|x|x|_| -->  |_|0|4|1|_| */
  
  // Right column, sup1deltax
  if ( relnl.x < - 0.9 * delta && relnl.x > - 1.1 * delta ) 
  {
    // top row, sup0deltay
    if ( fabs(relnl.y) < 0.1 * delta ) 
    {
      *weight_id = weight_numbers[2];
      *goflag = 1;
    }
    // middle row, inf1deltay
    else if ( relnl.y < 1.1 * delta && relnl.y > 0.9 * delta ) 
    {
      *weight_id = weight_numbers[5]; 
      *goflag = 1; 
    }
    // bottom row, inf2deltay
    else if ( relnl.y < 2.1 * delta  &&  relnl.y > 1.9 * delta ) 
    {
      *weight_id = weight_numbers[1]; 
      *goflag = 1;
    }
  }

  // Middle column, sup0deltax
  else if ( fabs(relnl.x) < 0.1 * delta ) 
  {
    // top row, sup0deltay
    if ( fabs(relnl.y) < 0.1 * delta ) 
    {
      *weight_id = weight_numbers[6]; 
      *goflag = 1;
    }
    // middle row, inf1deltay
    else if ( relnl.y < 1.1 * delta && relnl.y > 0.9 * delta ) 
    {
      *weight_id = weight_numbers[8]; 
      *goflag = 1; 
    }
    // bottom row, inf2deltay
    else if ( relnl.y < 2.1 * delta  &&  relnl.y > 1.9 * delta ) 
    {
      *weight_id = weight_numbers[4];
      *goflag = 1;
    }
  }
  
  // Left column, inf1deltax
  else if ( relnl.x < 1.1 * delta && relnl.x > 0.9 * delta ) 
  {
    // top row, sup0deltay
    if ( fabs(relnl.y) < 0.1 * delta ) 
    {
      *weight_id = weight_numbers[3]; 
      *goflag = 1;
    }
    // middle row, inf1deltay
    else if ( relnl.y < 1.1 * delta && relnl.y > 0.9 * delta ) 
    {
      *weight_id = weight_numbers[7]; 
      *goflag = 1; 
    }
    // bottom row, inf2deltay
    else if ( relnl.y < 2.1 * delta  &&  relnl.y > 1.9 * delta ) 
    {
      *weight_id = weight_numbers[0]; 
      *goflag = 1;
    }
  }
}




//----------------------------------------------------------------------------
void NCX4_CX4( const coord relnl, const double delta, 
	const int weight_numbers[9], int* weight_id, size_t* goflag ) 
//----------------------------------------------------------------------------
{
  // weight_numbers[2] -> right column top
  // weight_numbers[5] -> right column middle
  // weight_numbers[1] -> right column bottom

  // weight_numbers[6] -> middle column top
  // weight_numbers[8] -> middle colum middle
  // weight_numbers[4] -> middle column bottom

  // weight_numbers[3] -> left column top
  // weight_numbers[7] -> left column middle
  // weight_numbers[0] -> left column bottom

  /* Stencil config, O is the cell containing the Lagrange multiplier */
  /* |_|_|_|_|_| -->  |_|_|_|_|_| */
  /* |_|_|_|_|_| -->  |_|_|_|_|_| */
  /* |_|_|O|x|x| -->  |_|_|3|6|2| */ 
  /* |_|_|x|x|x| -->  |_|_|7|8|5| */
  /* |_|_|x|x|x| -->  |_|_|0|4|1| */
  
  // Right column, sup2deltax
  if ( relnl.x < -1.9 * delta  &&  relnl.x > - 2.1 * delta ) 
  {
    // top row, sup0deltay
    if ( fabs(relnl.y) < 0.1 * delta ) 
    { 
      *weight_id = weight_numbers[2]; 
      *goflag = 1;
    }
    // middle row, inf1deltay
    else if ( relnl.y < 1.1 * delta && relnl.y > 0.9 * delta ) 
    {
      *weight_id = weight_numbers[5]; 
      *goflag = 1; 
    }
    // bottom row, inf2deltay
    else if ( relnl.y < 2.1 * delta  &&  relnl.y > 1.9 * delta ) 
    {
      *weight_id = weight_numbers[1];
      *goflag = 1;
    }
  }
  
  // Middle column, sup1deltax
  else if ( relnl.x < - 0.9 * delta && relnl.x > - 1.1 * delta ) 
  {
    // top row, sup0deltay
    if ( fabs(relnl.y) < 0.1 * delta ) 
    {
      *weight_id = weight_numbers[6]; 
      *goflag = 1;
    }
    // middle row, inf1deltay
    else if ( relnl.y < 1.1 * delta && relnl.y > 0.9 * delta ) 
    {
      *weight_id = weight_numbers[8]; 
      *goflag = 1; 
    }
    // bottom row, inf2deltay
    else if ( relnl.y < 2.1 * delta  &&  relnl.y > 1.9 * delta ) 
    {
      *weight_id = weight_numbers[4]; 
      *goflag = 1;
    }
  }

  // Left column, sup0deltax
  else if ( fabs(relnl.x) < 0.1 * delta ) 
  {    
    // top row, sup0deltay
    if ( fabs(relnl.y) < 0.1 * delta ) 
    {
      *weight_id = weight_numbers[3]; 
      *goflag = 1;
    }
    // middle row, inf1deltay
    else if ( relnl.y < 1.1 * delta && relnl.y > 0.9 * delta ) 
    {
      *weight_id = weight_numbers[7]; 
      *goflag = 1; 
    }
    // bottom row, inf2deltay
    else if ( relnl.y < 2.1 * delta  &&  relnl.y > 1.9 * delta ) 
    {
      *weight_id = weight_numbers[0]; 
      *goflag = 1;
    }
  }
}




//----------------------------------------------------------------------------
void fill_weight_numbers( const int wisz, int* weight_numbers ) 
//----------------------------------------------------------------------------
{
  // In the plane:
  // weight_numbers[2] -> right column top
  // weight_numbers[5] -> right column middle
  // weight_numbers[1] -> right column bottom

  // weight_numbers[6] -> middle column top
  // weight_numbers[8] -> middle colum middle
  // weight_numbers[4] -> middle column bottom

  // weight_numbers[3] -> left column top
  // weight_numbers[7] -> left column middle
  // weight_numbers[0] -> left column bottom

  for (int i = 0; i < 9; i ++)  weight_numbers[i] = -1;

  switch( wisz )
  {    
    // wisz == -1 backward plane
    case -1:
      weight_numbers[3] = 6; weight_numbers[6] = 7; weight_numbers[2] = 8;
      weight_numbers[7] = 5; weight_numbers[8] = 4; weight_numbers[5] = 3;
      weight_numbers[0] = 0; weight_numbers[4] = 1; weight_numbers[1] = 2;
      break;
  
    // wisz == 0 middle plane 
    case 0:
      weight_numbers[3] = 15; weight_numbers[6] = 16; weight_numbers[2] = 17;
      weight_numbers[7] = 14; weight_numbers[8] = 13; weight_numbers[5] = 12;
      weight_numbers[0] = 9;  weight_numbers[4] = 10; weight_numbers[1] = 11;
      break;

    // wisz == 1 forward plane 
    case 1: 
      weight_numbers[3] = 24; weight_numbers[6] = 25; weight_numbers[2] = 26;
      weight_numbers[7] = 23; weight_numbers[8] = 22; weight_numbers[5] = 21;
      weight_numbers[0] = 18; weight_numbers[4] = 19; weight_numbers[1] = 20;
      break;
      
    default:
      printf("Warning: unexpected value of wisz in fill_weight_numbers\n");
  }
}




#if dimension == 3 
// NCX10_CX Bloc, 8 functions
//----------------------------------------------------------------------------
void NCX10_CX10( const coord relnl, const double delta, int* weight_id, 
	size_t* goflag ) 
//----------------------------------------------------------------------------
{
  int weight_numbers[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};
  
  /* Stencil config, O is the cell containing the Lagrange multiplier */
  
  /* z = localplane (plane containing the multiplier) */
  /* |_|_|x|x|x| --> |_|_|6|7|8| */
  /* |_|_|x|x|x| --> |_|_|5|4|3| */
  /* |_|_|O|x|x| --> |_|_|0|1|2| */ 
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  
  // backward z plane: z = zlocal : sup0deltaz
  if ( fabs(relnl.z) < 0.1 * delta ) 
  {
    fill_weight_numbers( -1, &weight_numbers[0] );
    NCX1_CX1( relnl, delta, weight_numbers, weight_id, goflag );
  }

  /* z = localplane + delta */
  /* |_|_|x|x|x| --> |_|_|15|16|18| */
  /* |_|_|x|x|x| --> |_|_|14|13|12| */
  /* |_|_|x|x|x| --> |_|_|_9|10|11| */ 
  /* |_|_|_|_|_| --> |_|_|__|__|__| */
  /* |_|_|_|_|_| --> |_|_|__|__|__| */
  
  // middle z plane: z = zlocal + Delta : sup1deltaz
  else if ( relnl.z < - 0.9 * delta && relnl.z > - 1.1 * delta ) 
  {
    fill_weight_numbers( 0, &weight_numbers[0] );
    NCX1_CX1( relnl, delta, weight_numbers, weight_id, goflag );
  }

  /* z = localplane + 2delta */
  /* |_|_|x|x|x| --> |_|_|24|25|26| */
  /* |_|_|x|x|x| --> |_|_|23|22|21| */
  /* |_|_|x|x|x| --> |_|_|18|19|20| */ 
  /* |_|_|_|_|_| --> |_|_|__|__|__| */
  /* |_|_|_|_|_| --> |_|_|__|__|__| */
  
  // forward z plane: z = zlocal + 2. * Delta : sup2deltaz
  else if ( relnl.z < - 1.9 * delta && relnl.z > - 2.1 * delta ) 
  {
    fill_weight_numbers( 1, &weight_numbers[0] );
    NCX1_CX1( relnl, delta, weight_numbers, weight_id, goflag );
  }
}




//----------------------------------------------------------------------------
void NCX10_CX20( const coord relnl, const double delta, int* weight_id, 
	size_t* goflag ) 
//----------------------------------------------------------------------------
{
  int weight_numbers[9] =  {-1, -1, -1, -1, -1, -1, -1, -1, -1};
  
 /* Stencil config, O is the cell containing the Lagrange multiplier */
  
  /* backward z plane: z = localplane : sup0deltaz (plane containing the 
  multiplier) */
  /* |_|x|x|x|_| --> |_|6|7|8|_| */
  /* |_|x|x|x|_| --> |_|5|4|3|_| */
  /* |_|x|O|x|_| --> |_|0|1|2|_| */ 
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  
  if ( fabs(relnl.z) < 0.1 * delta ) 
  {
    fill_weight_numbers( -1, &weight_numbers[0] );
    NCX1_CX2( relnl, delta, weight_numbers, weight_id, goflag );
  }

  /* middle z plane: z = zlocal + Delta : sup1deltaz */
  /* |_|x|x|x|_| --> |_|15|16|17|_| */
  /* |_|x|x|x|_| --> |_|14|13|12|_| */
  /* |_|x|x|x|_| --> |_|_9|10|11|_| */ 
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  
  else if ( relnl.z < - 0.9 * delta && relnl.z > - 1.1 * delta ) 
  {
    fill_weight_numbers( 0, &weight_numbers[0] );
    NCX1_CX2( relnl, delta, weight_numbers, weight_id, goflag );
  }
  
  /* forward z plane: z = zlocal + 2Delta : sup2deltaz */
  /* |_|x|x|x|_| --> |_|24|25|26|_| */
  /* |_|x|x|x|_| --> |_|23|22|21|_| */
  /* |_|x|x|x|_| --> |_|18|19|20|_| */ 
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  
  else if ( relnl.z < - 1.9 * delta && relnl.z > - 2.1 * delta ) 
  {
    fill_weight_numbers(1, &weight_numbers[0] );
    NCX1_CX2( relnl, delta, weight_numbers, weight_id, goflag );
  }
}




//----------------------------------------------------------------------------
void NCX10_CX30( const coord relnl, const double delta, int* weight_id, 
	size_t* goflag ) 
//----------------------------------------------------------------------------
{
  int weight_numbers[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};

  /* Stencil config, O is the cell containing the Lagrange multiplier */
  
  /* backward z plane: z = localplane : sup0deltaz (plane containing the 
  multiplier) */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|x|x|x|_| --> |_|6|7|8|_| */
  /* |_|x|O|x|_| --> |_|5|4|3|_| */ 
  /* |_|x|x|x|_| --> |_|0|1|2|_| */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */

  if ( fabs(relnl.z) < 0.1 * delta ) 
  {
    fill_weight_numbers( -1, &weight_numbers[0] );
    NCX_centred( relnl, delta, weight_numbers, weight_id, goflag );
  }

  /* middle z plane: z = zlocal + Delta : sup1deltaz */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  /* |_|x|x|x|_| --> |_|15|16|17|_| */
  /* |_|x|x|x|_| --> |_|14|13|12|_| */ 
  /* |_|x|x|x|_| --> |_|_9|10|11|_| */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  
  else if ( relnl.z < - 0.9 * delta && relnl.z > - 1.1 * delta ) 
  {
    fill_weight_numbers( 0, &weight_numbers[0] );
    NCX_centred( relnl, delta, weight_numbers, weight_id, goflag );
  }

  /* middle z plane: z = zlocal + 2Delta : sup2deltaz */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  /* |_|x|x|x|_| --> |_|24|25|26|_| */
  /* |_|x|x|x|_| --> |_|23|22|21|_| */ 
  /* |_|x|x|x|_| --> |_|18|19|20|_| */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  
  else if ( relnl.z < - 1.9 * delta && relnl.z > - 2.1 * delta ) 
  {
    fill_weight_numbers( 1, &weight_numbers[0] );
    NCX_centred( relnl, delta, weight_numbers, weight_id, goflag );
  }
}




//----------------------------------------------------------------------------
void NCX10_CX40( const coord relnl, const double delta, int* weight_id, 
	size_t* goflag ) 
//----------------------------------------------------------------------------
{
  int weight_numbers[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};

  /* Stencil config, O is the cell containing the Lagrange multiplier */
  /* backward z plane: z = localplane : sup0deltaz (plane containing the 
  multiplier) */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|_|O|x|x| --> |_|6|7|8|_| */ 
  /* |_|_|x|x|x| --> |_|5|4|3|_| */
  /* |_|_|x|x|x| --> |_|0|1|2|_| */

  if ( fabs(relnl.z) < 0.1 * delta ) 
  {
    fill_weight_numbers( -1, &weight_numbers[0] );
    NCX1_CX4( relnl, delta, weight_numbers, weight_id, goflag );
  }

  /* middle z plane: z = zlocal + Delta : sup1deltaz */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  /* |_|_|x|x|x| --> |_|15|16|17|_| */ 
  /* |_|_|x|x|x| --> |_|14|13|12|_| */
  /* |_|_|x|x|x| --> |_|_9|10|11|_| */

  else if ( relnl.z < - 0.9 * delta && relnl.z > - 1.1 * delta ) 
  {
    fill_weight_numbers( 0, &weight_numbers[0] );
    NCX1_CX4(relnl, delta, weight_numbers, weight_id, goflag );
  }
  
  /* forward z plane: z = zlocal + 2. * Delta : sup2deltaz */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  /* |_|_|x|x|x| --> |_|24|25|26|_| */ 
  /* |_|_|x|x|x| --> |_|23|22|21|_| */
  /* |_|_|x|x|x| --> |_|18|19|20|_| */
  
  else if ( relnl.z < - 1.9 * delta && relnl.z > - 2.1 * delta ) 
  {
    fill_weight_numbers( 1, &weight_numbers[0] );
    NCX1_CX4( relnl, delta, weight_numbers, weight_id, goflag );
  }
}




//----------------------------------------------------------------------------
void NCX10_CX11( const coord relnl, const double delta, int* weight_id, 
	size_t* goflag ) 
//----------------------------------------------------------------------------
{
  int weight_numbers[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};

  /* Stencil config, O is the cell containing the Lagrange multiplier */
  
  /* backward z plane: z = zlocal - Delta : inf1deltaz*/
  /* |_|_|x|x|x| --> |_|_|6|7|8| */
  /* |_|_|x|x|x| --> |_|_|5|4|3| */
  /* |_|_|x|x|x| --> |_|_|0|1|2| */ 
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  
  if ( relnl.z < 1.1 * delta && relnl.z > 0.9 * delta ) 
  {
    fill_weight_numbers( -1, &weight_numbers[0] );
    NCX1_CX1( relnl, delta, weight_numbers, weight_id, goflag );
  }

  /* middle z plane: z = zlocal : sup0deltaz */
  /* |_|_|x|x|x| --> |_|_|15|16|17| */
  /* |_|_|x|x|x| --> |_|_|14|13|12| */
  /* |_|_|O|x|x| --> |_|_|_9|10|11| */ 
  /* |_|_|_|_|_| --> |_|_|__|__|__| */
  /* |_|_|_|_|_| --> |_|_|__|__|__| */
    
  else if ( fabs(relnl.z) < 0.1 * delta ) 
  {
    fill_weight_numbers( 0, &weight_numbers[0] );
    NCX1_CX1( relnl, delta, weight_numbers, weight_id, goflag );
  }
  
  /* forward z plane: z = zlocal + Delta : sup1deltaz */
  /* |_|_|x|x|x| --> |_|_|24|25|26| */
  /* |_|_|x|x|x| --> |_|_|23|22|21| */
  /* |_|_|x|x|x| --> |_|_|18|19|20| */ 
  /* |_|_|_|_|_| --> |_|_|__|__|__| */
  /* |_|_|_|_|_| --> |_|_|__|__|__| */

  else if ( relnl.z < - 0.9 * delta && relnl.z > - 1.1 * delta ) 
  {
    fill_weight_numbers( 1, &weight_numbers[0] );
    NCX1_CX1( relnl, delta, weight_numbers, weight_id, goflag );
  }
}




//----------------------------------------------------------------------------
void NCX10_CX21( const coord relnl, const double delta, int* weight_id, 
	size_t* goflag ) 
//----------------------------------------------------------------------------
{
  int weight_numbers[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};

  /* Stencil config, O is the cell containing the Lagrange multiplier */
  /* backward z plane: z = zlocal - Delta : inf1deltaz*/
  /* |_|x|x|x|_| -->  |_|6|7|8|_| */
  /* |_|x|x|x|_| -->  |_|5|4|3|_| */
  /* |_|x|x|x|_| -->  |_|0|1|2|_| */ 
  /* |_|_|_|_|_| -->  |_|_|_|_|_| */
  /* |_|_|_|_|_| -->  |_|_|_|_|_| */

  if ( relnl.z < 1.1 * delta && relnl.z > 0.9 * delta ) 
  {
    fill_weight_numbers( -1, &weight_numbers[0] );
    NCX1_CX2( relnl, delta, weight_numbers, weight_id, goflag );
  }
  /* middle z plane: z = zlocal : sup0deltaz*/
  /* |_|x|x|x|_| -->  |_|15|16|17|_| */
  /* |_|x|x|x|_| -->  |_|14|13|12|_| */
  /* |_|x|O|x|_| -->  |_|_9|10|11|_| */ 
  /* |_|_|_|_|_| -->  |_|__|__|__|_| */
  /* |_|_|_|_|_| -->  |_|__|__|__|_| */
  
  else if ( fabs(relnl.z) < 0.1 * delta ) 
  {
    fill_weight_numbers( 0, &weight_numbers[0] );
    NCX1_CX2( relnl, delta, weight_numbers, weight_id, goflag );
  }
  
  /* forward z plane: z = zlocal + Delta : sup1deltaz */
  /* |_|x|x|x|_| -->  |_|24|25|26|_| */
  /* |_|x|x|x|_| -->  |_|23|22|21|_| */
  /* |_|x|x|x|_| -->  |_|18|19|20|_| */ 
  /* |_|_|_|_|_| -->  |_|__|__|__|_| */
  /* |_|_|_|_|_| -->  |_|__|__|__|_| */
  
  else if ( relnl.z < - 0.9 * delta && relnl.z > - 1.1 * delta ) 
  {
    fill_weight_numbers( 1, &weight_numbers[0] );
    NCX1_CX2( relnl, delta, weight_numbers, weight_id, goflag );
  }
}




//----------------------------------------------------------------------------
void NCX10_CX31( const coord relnl, const double delta, int* weight_id, 
	size_t* goflag ) 
//----------------------------------------------------------------------------
{
  int weight_numbers[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};

  /* Stencil config, O is the cell containing the Lagrange multiplier */
  /* backward z plane: z = zlocal - Delta : inf1deltaz*/
  /* |_|_|_|_|_| -->  |_|_|_|_|_| */
  /* |_|x|x|x|_| -->  |_|6|7|8|_| */
  /* |_|x|x|x|_| -->  |_|5|4|3|_| */ 
  /* |_|x|x|x|_| -->  |_|0|1|2|_| */
  /* |_|_|_|_|_| -->  |_|_|_|_|_| */
 
  if ( relnl.z < 1.1 * delta && relnl.z > 0.9 * delta ) 
  {
    fill_weight_numbers( -1, &weight_numbers[0] );
    NCX_centred( relnl, delta, weight_numbers, weight_id, goflag );
  }

  /* middle z plane: z = zlocal : sup0deltaz */
  /* |_|_|_|_|_| -->  |_|__|__|__|_| */
  /* |_|x|x|x|_| -->  |_|15|16|17|_| */
  /* |_|x|O|x|_| -->  |_|14|13|12|_| */ 
  /* |_|x|x|x|_| -->  |_|_9|10|11|_| */
  /* |_|_|_|_|_| -->  |_|__|__|__|_| */
  
  else if ( fabs(relnl.z) < 0.1 * delta ) 
  {
    fill_weight_numbers( 0, &weight_numbers[0] );
    NCX_centred( relnl, delta, weight_numbers, weight_id, goflag );
  }
  
  /* forward z plane: z = zlocal + Delta : sup1deltaz */
  /* |_|_|_|_|_| -->  |_|__|__|__|_| */
  /* |_|x|x|x|_| -->  |_|24|25|26|_| */
  /* |_|x|x|x|_| -->  |_|23|22|21|_| */ 
  /* |_|x|x|x|_| -->  |_|18|19|20|_| */
  /* |_|_|_|_|_| -->  |_|__|__|__|_| */
  
  else if ( relnl.z < - 0.9 * delta && relnl.z > - 1.1 * delta ) 
  {
    fill_weight_numbers( 1, &weight_numbers[0] );
    NCX_centred( relnl, delta, weight_numbers, weight_id, goflag );
  }
}




//----------------------------------------------------------------------------
void NCX10_CX41( const coord relnl, const double delta, int* weight_id, 
	size_t* goflag ) 
//----------------------------------------------------------------------------
{
  int weight_numbers[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};
  
  /* Stencil config, O is the cell containing the Lagrange multiplier */
  /* backward z plane: z = zlocal - Delta : inf1deltaz */
  /* |_|_|_|_|_| -->  |_|_|_|_|_| */
  /* |_|_|x|x|x| -->  |_|_|6|7|8| */
  /* |_|_|x|x|x| -->  |_|_|5|4|3| */ 
  /* |_|_|x|x|x| -->  |_|_|0|1|2| */
  /* |_|_|_|_|_| -->  |_|_|_|_|_| */

  if ( relnl.z < 1.1 * delta && relnl.z > 0.9 * delta ) 
  {
    fill_weight_numbers( -1, &weight_numbers[0] );
    NCX1_CX4( relnl, delta, weight_numbers, weight_id, goflag );
  }

  /* middle z plane: z = zlocal : sup0deltaz */
  /* |_|_|_|_|_| -->  |_|_|__|__|__| */
  /* |_|_|x|x|x| -->  |_|_|15|16|17| */
  /* |_|_|O|x|x| -->  |_|_|14|13|12| */ 
  /* |_|_|x|x|x| -->  |_|_|_9|10|11| */
  /* |_|_|_|_|_| -->  |_|_|__|__|__| */

  else if ( fabs(relnl.z) < 0.1 * delta ) 
  {
    fill_weight_numbers( 0, &weight_numbers[0] );
    NCX1_CX4( relnl, delta, weight_numbers, weight_id, goflag );
  }
  
  // forward z plane: z = zlocal + Delta : sup1deltaz
  /* |_|_|_|_|_| -->  |_|_|__|__|__| */
  /* |_|_|x|x|x| -->  |_|_|24|25|26| */
  /* |_|_|x|x|x| -->  |_|_|23|22|21| */ 
  /* |_|_|x|x|x| -->  |_|_|18|19|20| */
  /* |_|_|_|_|_| -->  |_|_|__|__|__| */
  
  else if ( relnl.z < - 0.9 * delta && relnl.z > - 1.1 * delta ) 
  {
    fill_weight_numbers( 1, &weight_numbers[0] );
    NCX1_CX4( relnl, delta, weight_numbers, weight_id, goflag );
  }
}




// NCX11_CX Bloc, 4 functions

/* NCX11_CX20 = NCX10_CX21 */

/* NCX11_CX30 = NCX10_CX31 */

/* NCX11_CX40 = NCX10_CX41 */

//----------------------------------------------------------------------------
void NCX11_CX11( const coord relnl, const double delta, int* weight_id, 
	size_t* goflag ) 
//----------------------------------------------------------------------------
{
  int weight_numbers[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};
  
  /* Stencil config, O is the cell containing the Lagrange multiplier */
  
  /* backward z plane: z = zlocal - 2Delta : inf2deltaz */
  /* |_|_|x|x|x| --> |_|_|6|7|8| */
  /* |_|_|x|x|x| --> |_|_|5|4|3| */
  /* |_|_|x|x|x| --> |_|_|0|1|2| */ 
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  
  if ( relnl.z < 2.1 * delta && relnl.z > 1.9 * delta ) 
  {
    fill_weight_numbers( -1, &weight_numbers[0] );
    NCX1_CX1( relnl, delta, weight_numbers, weight_id, goflag );
  }

  /* middle z plane: z = zlocal - Delta : inf1deltaz */
  /* |_|_|x|x|x| --> |_|_|15|16|17| */
  /* |_|_|x|x|x| --> |_|_|14|13|12| */
  /* |_|_|x|x|x| --> |_|_|_9|10|11| */ 
  /* |_|_|_|_|_| --> |_|_|__|__|__| */
  /* |_|_|_|_|_| --> |_|_|__|__|__| */
  
  else if ( relnl.z < 1.1 * delta && relnl.z > 0.9 * delta ) 
  {
    fill_weight_numbers( 0, &weight_numbers[0] );
    NCX1_CX1( relnl, delta, weight_numbers, weight_id, goflag );
  }
  
  /* forward z plane: z = zlocal : sup0deltaz */
  /* |_|_|x|x|x| --> |_|_|24|25|26| */
  /* |_|_|x|x|x| --> |_|_|23|22|21| */
  /* |_|_|O|x|x| --> |_|_|18|19|20| */ 
  /* |_|_|_|_|_| --> |_|_|__|__|__| */
  /* |_|_|_|_|_| --> |_|_|__|__|__| */
  
  else if ( fabs(relnl.z) < 0.1 * delta ) 
  {
    fill_weight_numbers( 1, &weight_numbers[0] );
    NCX1_CX1( relnl, delta, weight_numbers, weight_id, goflag );
  }
}




//----------------------------------------------------------------------------
void NCX11_CX21( const coord relnl, const double delta, int* weight_id, 
	size_t* goflag ) 
//----------------------------------------------------------------------------
{
  int weight_numbers[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};
  /* Stencil config, O is the cell containing the Lagrange multiplier */
  
  /* backward z plane: z = zlocal - 2Delta : inf2deltaz */
  /* |_|x|x|x|_| --> |_|6|7|8|_| */
  /* |_|x|x|x|_| --> |_|5|4|3|_| */
  /* |_|x|x|x|_| --> |_|0|1|2|_| */ 
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  
  if ( relnl.z < 2.1 * delta && relnl.z > 1.9 * delta ) 
  {
    fill_weight_numbers( -1, &weight_numbers[0] );
    NCX1_CX2( relnl, delta, weight_numbers, weight_id, goflag );
  }

  /* middle z plane: z = zlocal - Delta : inf1deltaz */
  /* |_|x|x|x|_| --> |_|15|16|17|_| */
  /* |_|x|x|x|_| --> |_|14|13|12|_| */
  /* |_|x|x|x|_| --> |_|_9|10|11|_| */ 
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  
  else if ( relnl.z < 1.1 * delta && relnl.z > 0.9 * delta ) 
  {
    fill_weight_numbers( 0, &weight_numbers[0] );
    NCX1_CX2( relnl, delta, weight_numbers, weight_id, goflag );
  }
  
  /* forward z plane: z = zlocal : sup0deltaz */
  /* |_|x|x|x|_| --> |_|24|25|26|_| */
  /* |_|x|x|x|_| --> |_|23|22|21|_| */
  /* |_|x|O|x|_| --> |_|18|19|20|_| */ 
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  
  else if ( fabs(relnl.z) < 0.1 * delta ) 
  {
    fill_weight_numbers( 1, &weight_numbers[0] );
    NCX1_CX2( relnl, delta, weight_numbers, weight_id, goflag );
  }
}




//----------------------------------------------------------------------------
void NCX11_CX31( const coord relnl, const double delta, int* weight_id, 
	size_t* goflag ) 
//----------------------------------------------------------------------------
{
  int weight_numbers[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};
  /* Stencil config, O is the cell containing the Lagrange multiplier */

  /* backward z plane: z = zlocal - 2Delta : inf2deltaz */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|x|x|x|_| --> |_|6|7|8|_| */
  /* |_|x|x|x|_| --> |_|5|4|3|_| */ 
  /* |_|x|x|x|_| --> |_|0|1|2|_| */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  
  if ( relnl.z < 2.1 * delta && relnl.z > 1.9 * delta ) 
  { 
    fill_weight_numbers( -1, &weight_numbers[0] );
    NCX_centred( relnl, delta, weight_numbers, weight_id, goflag );
  }

  /* middle z plane: z = zlocal - Delta : inf1deltaz */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  /* |_|x|x|x|_| --> |_|15|16|17|_| */
  /* |_|x|x|x|_| --> |_|14|13|12|_| */ 
  /* |_|x|x|x|_| --> |_|_9|10|11|_| */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  
  else if ( relnl.z < 1.1 * delta && relnl.z > 0.9 * delta ) 
  {
    fill_weight_numbers( 0, &weight_numbers[0] );
    NCX_centred( relnl, delta, weight_numbers, weight_id, goflag );
  }
  
  /* forward z plane: z = zlocal : sup0deltaz */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  /* |_|x|x|x|_| --> |_|24|25|26|_| */
  /* |_|x|x|x|_| --> |_|23|22|21|_| */ 
  /* |_|x|x|x|_| --> |_|18|19|20|_| */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  
  else if ( fabs(relnl.z) < 0.1 * delta ) 
  {
    fill_weight_numbers( 1, &weight_numbers[0] );
    NCX_centred( relnl, delta, weight_numbers, weight_id, goflag );
  }
}




//----------------------------------------------------------------------------
void NCX11_CX41( const coord relnl, const double delta, int* weight_id, 
	size_t* goflag ) 
//----------------------------------------------------------------------------
{
  int weight_numbers[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};
  /* Stencil config, O is the cell containing the Lagrange multiplier */

  /* backward z plane: z = zlocal - 2Delta : inf2deltaz */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|_|x|x|x| --> |_|_|6|7|8| */
  /* |_|_|x|x|x| --> |_|_|5|4|3| */ 
  /* |_|_|x|x|x| --> |_|_|0|1|2| */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  
  if ( relnl.z < 2.1 * delta && relnl.z > 1.9 * delta ) 
  {
    fill_weight_numbers( -1, &weight_numbers[0] );
    NCX1_CX4( relnl, delta, weight_numbers, weight_id, goflag );
  }

  /* middle z plane: z = zlocal - Delta : inf1deltaz */
  /* |_|_|_|_|_| --> |_|_|__|__|__| */
  /* |_|_|x|x|x| --> |_|_|15|16|17| */
  /* |_|_|x|x|x| --> |_|_|14|13|12| */ 
  /* |_|_|x|x|x| --> |_|_|_9|10|11| */
  /* |_|_|_|_|_| --> |_|_|__|__|__| */
  
  else if ( relnl.z < 1.1 * delta && relnl.z > 0.9 * delta ) 
  {
    fill_weight_numbers( 0, &weight_numbers[0] );
    NCX1_CX4( relnl, delta, weight_numbers, weight_id, goflag );
  }
  
  /* forward z plane: z = zlocal : sup0deltaz */
  /* |_|_|_|_|_| --> |_|_|__|__|__| */
  /* |_|_|x|x|x| --> |_|_|24|25|26| */
  /* |_|_|x|x|x| --> |_|_|23|22|21| */ 
  /* |_|_|x|x|x| --> |_|_|18|19|20| */
  /* |_|_|_|_|_| --> |_|_|__|__|__| */
  
  else if ( fabs(relnl.z) < 0.1 * delta ) 
  {
    fill_weight_numbers( 1, &weight_numbers[0] );
    NCX1_CX4( relnl, delta, weight_numbers, weight_id, goflag );
  }
}




// NCX20_CX Bloc, 6 functions

//----------------------------------------------------------------------------
void NCX20_CX20( const coord relnl, const double delta, int* weight_id, 
	size_t* goflag ) 
//----------------------------------------------------------------------------
{
  int weight_numbers[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};
  
  /* Stencil config, O is the cell containing the Lagrange multiplier */
  
  /* z = localplane (plane containing the multiplier) */
  /* backward z plane: z = zlocal : sup0deltaz */
  /* |x|x|x|_|_| --> |6|7|8|_|_| */
  /* |x|x|x|_|_| --> |5|4|3|_|_| */
  /* |x|x|O|_|_| --> |0|1|2|_|_| */ 
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  
  if ( fabs(relnl.z) < 0.1 * delta ) 
  {
    fill_weight_numbers( -1, &weight_numbers[0] );
    NCX2_CX2( relnl, delta, weight_numbers, weight_id, goflag );
  }

  /* middle z plane: z = zlocal + Delta : sup1delta*/
  /* |x|x|x|_|_| --> |15|16|17|_|_| */
  /* |x|x|x|_|_| --> |14|13|12|_|_| */
  /* |x|x|x|_|_| --> |_9|10|11|_|_| */ 
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  /* |_|_|_|_|_| --> |__|__|__|_|_| */

  else if ( relnl.z < - 0.9 * delta && relnl.z > - 1.1 * delta ) 
  {
    fill_weight_numbers( 0, &weight_numbers[0] );
    NCX2_CX2( relnl, delta, weight_numbers, weight_id, goflag );
  }
  
  /* forward z plane: z = zlocal + 2. * Delta : sup2delta*/
  /* |x|x|x|_|_| --> |24|25|26|_|_| */
  /* |x|x|x|_|_| --> |23|22|21|_|_| */
  /* |x|x|x|_|_| --> |18|19|20|_|_| */ 
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  
  else if ( relnl.z < - 1.9 * delta && relnl.z > - 2.1 * delta ) 
  {
    fill_weight_numbers( 1, &weight_numbers[0] );
    NCX2_CX2( relnl, delta, weight_numbers, weight_id, goflag );
  }
}




//----------------------------------------------------------------------------
void NCX20_CX30( const coord relnl, const double delta, int* weight_id, 
	size_t* goflag ) 
//----------------------------------------------------------------------------
{
  int weight_numbers[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};
  
  /* Stencil config, O is the cell containing the Lagrange multiplier */
  
  /* z = localplane (plane containing the multiplier) */
  /* backward z plane: z = zlocal : sup0deltaz */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |x|x|x|_|_| --> |6|7|8|_|_| */
  /* |x|x|O|_|_| --> |5|4|3|_|_| */ 
  /* |x|x|x|_|_| --> |0|1|2|_|_| */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  
  if ( fabs(relnl.z) < 0.1 * delta ) 
  {
    fill_weight_numbers( -1, &weight_numbers[0] );
    NCX2_CX3( relnl, delta, weight_numbers, weight_id, goflag );
  }

  /* middle z plane: z = zlocal + Delta : sup1deltaz */
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  /* |x|x|x|_|_| --> |15|16|17|_|_| */
  /* |x|x|x|_|_| --> |14|13|12|_|_| */ 
  /* |x|x|x|_|_| --> |_9|10|11|_|_| */
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  
  else if ( relnl.z < - 0.9 * delta && relnl.z > - 1.1 * delta ) 
  {
    fill_weight_numbers( 0, &weight_numbers[0] );
    NCX2_CX3( relnl, delta, weight_numbers, weight_id, goflag );
  }
  
  /* forward z plane: z = zlocal + 2. * Delta : sup2deltaz */
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  /* |x|x|x|_|_| --> |24|25|26|_|_| */
  /* |x|x|x|_|_| --> |23|22|21|_|_| */ 
  /* |x|x|x|_|_| --> |18|19|20|_|_| */
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  
  else if ( relnl.z < - 1.9 * delta && relnl.z > - 2.1 * delta ) 
  {
    fill_weight_numbers( 1, &weight_numbers[0] );
    NCX2_CX3( relnl, delta, weight_numbers, weight_id, goflag );
  }  
}




//----------------------------------------------------------------------------
void NCX20_CX40( const coord relnl, const double delta, int* weight_id, 
	size_t* goflag ) 
//----------------------------------------------------------------------------
{
  int weight_numbers[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};
  
  /* Stencil config, O is the cell containing the Lagrange multiplier */
  
  /* z = localplane (plane containing the multiplier) */
  /* backward z plane: z = zlocal : sup0deltaz */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|x|x|x|_| --> |_|6|7|8|_| */
  /* |_|x|O|x|_| --> |_|5|4|3|_| */ 
  /* |_|x|x|x|_| --> |_|0|1|2|_| */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  
  if ( fabs(relnl.z) < 0.1 * delta ) 
  {
    fill_weight_numbers( -1, &weight_numbers[0] );
    NCX_centred( relnl, delta, weight_numbers, weight_id, goflag );
  }

  // middle z plane: z = zlocal + Delta : sup1deltaz
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  /* |_|x|x|x|_| --> |_|15|16|17|_| */
  /* |_|x|x|x|_| --> |_|14|13|12|_| */ 
  /* |_|x|x|x|_| --> |_|_9|10|11|_| */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  
  else if ( relnl.z < - 0.9 * delta && relnl.z > - 1.1 * delta ) 
  {
    fill_weight_numbers( 0, &weight_numbers[0] );
    NCX_centred( relnl, delta, weight_numbers, weight_id, goflag );
  }
  
  /* forward z plane: z = zlocal + 2. * Delta : sup2deltaz */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  /* |_|x|x|x|_| --> |_|24|25|26|_| */
  /* |_|x|x|x|_| --> |_|23|22|21|_| */ 
  /* |_|x|x|x|_| --> |_|18|19|20|_| */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */

  else if ( relnl.z < - 1.9 * delta && relnl.z > - 2.1 * delta ) 
  {
    fill_weight_numbers( 1, &weight_numbers[0] );
    NCX_centred( relnl, delta, weight_numbers, weight_id, goflag );
  }  
}




//----------------------------------------------------------------------------
void NCX20_CX21( const coord relnl, const double delta, int* weight_id, 
	size_t* goflag ) 
//----------------------------------------------------------------------------
{
  int weight_numbers[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};

  /* Stencil config, O is the cell containing the Lagrange multiplier */
  
  /* backward z plane: z = zlocal - Delta : inf1deltaz*/
  /* |x|x|x|_|_| --> |6|7|8|_|_| */
  /* |x|x|x|_|_| --> |5|4|3|_|_| */
  /* |x|x|x|_|_| --> |0|1|2|_|_| */ 
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  
  if ( relnl.z < 1.1 * delta && relnl.z > 0.9 * delta ) 
  {
    fill_weight_numbers( -1, &weight_numbers[0] );
    NCX2_CX2( relnl, delta, weight_numbers, weight_id, goflag );
  }

  /* middle z plane: z = zlocal : sup0deltaz */
  /* |x|x|x|_|_| --> |15|16|17|_|_| */
  /* |x|x|x|_|_| --> |14|13|12|_|_| */
  /* |x|x|O|_|_| --> |_9|10|11|_|_| */ 
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  
  else if ( fabs(relnl.z) < 0.1 * delta ) 
  {
    fill_weight_numbers( 0, &weight_numbers[0] );
    NCX2_CX2( relnl, delta, weight_numbers, weight_id, goflag );
  }
  
  /* forward z plane: z = zlocal + Delta : sup1deltaz */
  /* |x|x|x|_|_| --> |24|25|26|_|_| */
  /* |x|x|x|_|_| --> |23|22|21|_|_| */
  /* |x|x|x|_|_| --> |18|19|20|_|_| */ 
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  
  else if ( relnl.z < - 0.9 * delta && relnl.z > - 1.1 * delta ) 
  {
    fill_weight_numbers( 1, &weight_numbers[0] );
    NCX2_CX2( relnl, delta, weight_numbers, weight_id, goflag );
  }
}




//----------------------------------------------------------------------------
void NCX20_CX31( const coord relnl, const double delta, int* weight_id, 
	size_t* goflag ) 
//----------------------------------------------------------------------------
{
  int weight_numbers[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};

  /* Stencil config, O is the cell containing the Lagrange multiplier */
  
  /* backward z plane: z = zlocal - Delta : inf1deltaz*/
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |x|x|x|_|_| --> |6|7|8|_|_| */
  /* |x|x|x|_|_| --> |5|4|3|_|_| */ 
  /* |x|x|x|_|_| --> |0|1|2|_|_| */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  
  if ( relnl.z < 1.1 * delta && relnl.z > 0.9 * delta ) 
  {
    fill_weight_numbers( -1, &weight_numbers[0] );
    NCX2_CX3( relnl, delta, weight_numbers, weight_id, goflag );
  }
  
  /* middle z plane: z = zlocal : sup0deltaz */
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  /* |x|x|x|_|_| --> |15|16|17|_|_| */
  /* |x|x|x|_|_| --> |14|13|12|_|_| */ 
  /* |x|x|x|_|_| --> |_9|10|11|_|_| */
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  
  else if ( fabs(relnl.z) < 0.1 * delta ) 
  {
    fill_weight_numbers( 0, &weight_numbers[0] );
    NCX2_CX3( relnl, delta, weight_numbers, weight_id, goflag );
  }
  
  /* forward z plane: z = zlocal + Delta : sup1deltaz */
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  /* |x|x|x|_|_| --> |24|25|26|_|_| */
  /* |x|x|x|_|_| --> |23|22|21|_|_| */ 
  /* |x|x|x|_|_| --> |18|19|20|_|_| */
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  
  else if ( relnl.z < - 0.9 * delta && relnl.z > - 1.1 * delta ) 
  {
    fill_weight_numbers( 1, &weight_numbers[0] );
    NCX2_CX3( relnl, delta, weight_numbers, weight_id, goflag );
  }
}




//----------------------------------------------------------------------------
void NCX20_CX41( const coord relnl, const double delta, int* weight_id, 
	size_t* goflag ) 
//----------------------------------------------------------------------------
{
  int weight_numbers[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};

  /* Stencil config, O is the cell containing the Lagrange multiplier */
  
  /* backward z plane: z = zlocal - Delta : inf1deltaz*/
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|x|x|x|_| --> |_|6|7|8|_| */
  /* |_|x|x|x|_| --> |_|5|4|3|_| */ 
  /* |_|x|x|x|_| --> |_|0|1|2|_| */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  
  if ( relnl.z < 1.1 * delta && relnl.z > 0.9 * delta ) 
  {
    fill_weight_numbers( -1, &weight_numbers[0] );
    NCX_centred( relnl, delta, weight_numbers, weight_id, goflag );
  }
  
  /* middle z plane: z = zlocal : sup0deltaz */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  /* |_|x|x|x|_| --> |_|15|16|17|_| */
  /* |_|x|x|x|_| --> |_|14|13|12|_| */ 
  /* |_|x|x|x|_| --> |_|_9|10|11|_| */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  
  else if ( fabs(relnl.z) < 0.1 * delta ) 
  {
    fill_weight_numbers( 0, &weight_numbers[0] );
    NCX_centred( relnl, delta, weight_numbers, weight_id, goflag );
  }
  
  /* forward z plane: z = zlocal + Delta : sup1deltaz */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  /* |_|x|x|x|_| --> |_|24|25|26|_| */
  /* |_|x|x|x|_| --> |_|23|22|21|_| */ 
  /* |_|x|x|x|_| --> |_|18|19|20|_| */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  
  else if ( relnl.z < - 0.9 * delta && relnl.z > - 1.1 * delta ) 
  {
    fill_weight_numbers( 1, &weight_numbers[0] );
    NCX_centred( relnl, delta, weight_numbers, weight_id, goflag );
  }
}




// NCX21_CX Bloc 3 functions

/* NCX21_CX30 = NCX20_CX31 */

/* NCX21_CX40 = NCX20_CX41 */

//----------------------------------------------------------------------------
void NCX21_CX21( const coord relnl, const double delta, int* weight_id, 
	size_t* goflag ) 
//----------------------------------------------------------------------------
{
  int weight_numbers[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};
  
  /* Stencil config, O is the cell containing the Lagrange multiplier */
  
  /* backward z plane: z = zlocal - 2Delta : inf2deltaz */
  /* |x|x|x|_|_| --> |6|7|8|_|_| */
  /* |x|x|x|_|_| --> |5|4|3|_|_| */
  /* |x|x|x|_|_| --> |0|1|2|_|_| */ 
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */

  if ( relnl.z < 2.1 * delta && relnl.z > 1.9 * delta ) 
  {
    fill_weight_numbers( -1, &weight_numbers[0] );
    NCX2_CX2( relnl, delta, weight_numbers, weight_id, goflag );
  }

  /* middle z plane: z = zlocal - Delta : inf1deltaz */
  /* |x|x|x|_|_| --> |15|16|17|_|_| */
  /* |x|x|x|_|_| --> |14|13|12|_|_| */
  /* |x|x|x|_|_| --> |_9|10|11|_|_| */ 
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  
  else if ( relnl.z < 1.1 * delta && relnl.z > 0.9 * delta ) 
  {
    fill_weight_numbers( 0, &weight_numbers[0] );
    NCX2_CX2( relnl, delta, weight_numbers, weight_id, goflag );
  }
  
  /* forward z plane: z = zlocal : sup0deltaz */
  /* |x|x|x|_|_| --> |24|25|26|_|_| */
  /* |x|x|x|_|_| --> |23|22|21|_|_| */
  /* |x|x|O|_|_| --> |18|19|20|_|_| */ 
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  
  else if ( fabs(relnl.z) < 0.1 * delta ) 
  {
    fill_weight_numbers( 1, &weight_numbers[0] );
    NCX2_CX2( relnl, delta, weight_numbers, weight_id, goflag );
  }
}




//----------------------------------------------------------------------------
void NCX21_CX31( const coord relnl, const double delta, int* weight_id, 
	size_t* goflag ) 
//----------------------------------------------------------------------------
{
  int weight_numbers[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};
  
  /* Stencil config, O is the cell containing the Lagrange multiplier */
  
  /* backward z plane: z = zlocal - 2Delta : inf2deltaz */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |x|x|x|_|_| --> |6|7|8|_|_| */
  /* |x|x|x|_|_| --> |5|4|3|_|_| */ 
  /* |x|x|x|_|_| --> |0|1|2|_|_| */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  
  if ( relnl.z < 2.1 * delta && relnl.z > 1.9 * delta ) 
  {
    fill_weight_numbers( -1, &weight_numbers[0] );
    NCX2_CX3(relnl, delta, weight_numbers, weight_id, goflag );
  }

  // middle z plane: z = zlocal - Delta : inf1deltaz
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  /* |x|x|x|_|_| --> |15|16|17|_|_| */
  /* |x|x|x|_|_| --> |14|13|12|_|_| */ 
  /* |x|x|x|_|_| --> |_9|10|11|_|_| */
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  
  else if ( relnl.z < 1.1 * delta && relnl.z > 0.9 * delta ) 
  {
    fill_weight_numbers( 0, &weight_numbers[0] );
    NCX2_CX3( relnl, delta, weight_numbers, weight_id, goflag );
  }
  
  /* forward z plane: z = zlocal : sup0deltaz */
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  /* |x|x|x|_|_| --> |24|25|26|_|_| */
  /* |x|x|O|_|_| --> |23|22|21|_|_| */ 
  /* |x|x|x|_|_| --> |18|19|20|_|_| */
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  
  else if ( fabs(relnl.z) < 0.1 * delta ) 
  {
    fill_weight_numbers( 1, &weight_numbers[0] );
    NCX2_CX3( relnl, delta, weight_numbers, weight_id, goflag );
  }
}




//----------------------------------------------------------------------------
void NCX21_CX41( const coord relnl, const double delta, int* weight_id, 
	size_t* goflag ) 
//----------------------------------------------------------------------------
{
  int weight_numbers[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};
  
  /* Stencil config, O is the cell containing the Lagrange multiplier */
  
  /* backward z plane: z = zlocal - 2Delta : inf2deltaz */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|x|x|x|_| --> |_|6|7|8|_| */
  /* |_|x|x|x|_| --> |_|5|4|3|_| */ 
  /* |_|x|x|x|_| --> |_|0|1|2|_| */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  
  if ( relnl.z < 2.1 * delta && relnl.z > 1.9 * delta ) 
  {
    fill_weight_numbers( -1, &weight_numbers[0] );
    NCX_centred( relnl, delta, weight_numbers, weight_id, goflag );
  }

  // middle z plane: z = zlocal - Delta : inf1deltaz
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  /* |_|x|x|x|_| --> |_|15|16|17|_| */
  /* |_|x|x|x|_| --> |_|14|13|12|_| */ 
  /* |_|x|x|x|_| --> |_|_9|10|11|_| */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  
  else if ( relnl.z < 1.1 * delta && relnl.z > 0.9 * delta ) 
  {
    fill_weight_numbers( 0, &weight_numbers[0] );
    NCX_centred( relnl, delta, weight_numbers, weight_id, goflag );
  }
  
  // forward z plane: z = zlocal : sup0deltaz
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  /* |_|x|x|x|_| --> |_|24|25|26|_| */
  /* |_|x|O|x|_| --> |_|23|22|21|_| */ 
  /* |_|x|x|x|_| --> |_|18|19|20|_| */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  
  else if ( fabs(relnl.z) < 0.1 * delta ) 
  {
    fill_weight_numbers( 1, &weight_numbers[0] );
    NCX_centred( relnl, delta, weight_numbers, weight_id, goflag );
  }
}




// NCX30_CX Bloc, 4 functions

//----------------------------------------------------------------------------
void NCX30_CX30( const coord relnl, const double delta, int* weight_id, 
	size_t* goflag ) 
//----------------------------------------------------------------------------
{
  int weight_numbers[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};
  
  /* Stencil config, O is the cell containing the Lagrange multiplier */
  
  /* z = localplane (plane containing the multiplier) */
  /* backward z plane: z = zlocal : sup0deltaz */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |x|x|O|_|_| --> |6|7|8|_|_| */ 
  /* |x|x|x|_|_| --> |5|4|3|_|_| */
  /* |x|x|x|_|_| --> |0|1|2|_|_| */
  
  if ( fabs(relnl.z) < 0.1 * delta ) 
  {
    fill_weight_numbers( -1, &weight_numbers[0] );
    NCX3_CX3( relnl, delta, weight_numbers, weight_id, goflag );
  }

  /* middle z plane: z = zlocal + Delta : sup1delta */
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  /* |x|x|x|_|_| --> |15|16|17|_|_| */ 
  /* |x|x|x|_|_| --> |14|13|12|_|_| */
  /* |x|x|x|_|_| --> |_9|10|11|_|_| */
  
  else if ( relnl.z < - 0.9 * delta && relnl.z > - 1.1 * delta ) 
  {
    fill_weight_numbers( 0, &weight_numbers[0] );
    NCX3_CX3( relnl, delta, weight_numbers, weight_id, goflag );
  }
  
  /* forward z plane: z = zlocal + 2. * Delta : sup2delta*/
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  /* |x|x|x|_|_| --> |24|25|26|_|_| */ 
  /* |x|x|x|_|_| --> |23|22|21|_|_| */
  /* |x|x|x|_|_| --> |18|19|20|_|_| */
  
  else if ( relnl.z < - 1.9 * delta && relnl.z > - 2.1 * delta ) 
  {
    fill_weight_numbers( 1, &weight_numbers[0] );
    NCX3_CX3( relnl, delta, weight_numbers, weight_id, goflag );
  }
}




//----------------------------------------------------------------------------
void NCX30_CX40( const coord relnl, const double delta, int* weight_id, 
	size_t* goflag ) 
//----------------------------------------------------------------------------
{
  int weight_numbers[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};
  
  /* Stencil config, O is the cell containing the Lagrange multiplier */
  
  /* z = localplane (plane containing the multiplier) */
  /* backward z plane: z = zlocal : sup0deltaz */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|x|O|x|_| --> |_|6|7|8|_| */ 
  /* |_|x|x|x|_| --> |_|5|4|3|_| */
  /* |_|x|x|x|_| --> |_|0|1|2|_| */

  if ( fabs(relnl.z) < 0.1 * delta ) 
  {
    fill_weight_numbers( -1, &weight_numbers[0] );
    NCX3_CX4( relnl, delta, weight_numbers, weight_id, goflag );
  }
 
  /* middle z plane: z = zlocal + Delta : sup1deltaz */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  /* |_|x|x|x|_| --> |_|15|16|17|_| */ 
  /* |_|x|x|x|_| --> |_|14|13|12|_| */
  /* |_|x|x|x|_| --> |_|_9|10|11|_| */
  
  else if ( relnl.z < - 0.9 * delta && relnl.z > - 1.1 * delta ) 
  {
    fill_weight_numbers( 0, &weight_numbers[0] );
    NCX3_CX4( relnl, delta, weight_numbers, weight_id, goflag );
  }
  
  /* forward z plane: z = zlocal + 2. * Delta : sup2deltaz */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  /* |_|x|x|x|_| --> |_|24|25|26|_| */ 
  /* |_|x|x|x|_| --> |_|23|22|21|_| */
  /* |_|x|x|x|_| --> |_|18|19|20|_| */
  
  else if ( relnl.z < - 1.9 * delta && relnl.z > - 2.1 * delta ) 
  {
    fill_weight_numbers( 1, &weight_numbers[0] );
    NCX3_CX4( relnl, delta, weight_numbers, weight_id, goflag );
  }
}




//----------------------------------------------------------------------------
void NCX30_CX31( const coord relnl, const double delta, int* weight_id, 
	size_t* goflag ) 
//----------------------------------------------------------------------------
{
  int weight_numbers[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};

  /* Stencil config, O is the cell containing the Lagrange multiplier */
  
  /* backward z plane: z = zlocal - Delta : inf1deltaz*/
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |x|x|x|_|_| --> |6|7|8|_|_| */ 
  /* |x|x|x|_|_| --> |5|4|3|_|_| */
  /* |x|x|x|_|_| --> |0|1|2|_|_| */
  
  if ( relnl.z < 1.1 * delta && relnl.z > 0.9 * delta ) 
  {
    fill_weight_numbers( -1, &weight_numbers[0] );
    NCX3_CX3( relnl, delta, weight_numbers, weight_id, goflag );
  }

  /* middle z plane: z = zlocal : sup0deltaz */
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  /* |x|x|x|_|_| --> |15|16|17|_|_| */ 
  /* |x|x|x|_|_| --> |14|13|12|_|_| */
  /* |x|x|x|_|_| --> |_9|10|11|_|_| */
  
  else if ( fabs(relnl.z) < 0.1 * delta ) 
  {
    fill_weight_numbers( 0, &weight_numbers[0] );
    NCX3_CX3( relnl, delta, weight_numbers, weight_id, goflag );
  }
  
  /* forward z plane: z = zlocal + Delta : sup1deltaz */
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  /* |x|x|x|_|_| --> |24|25|26|_|_| */ 
  /* |x|x|x|_|_| --> |23|22|21|_|_| */
  /* |x|x|x|_|_| --> |18|19|20|_|_| */
  
  else if ( relnl.z < - 0.9 * delta && relnl.z > - 1.1 * delta ) 
  {
    fill_weight_numbers( 1, &weight_numbers[0] );
    NCX3_CX3( relnl, delta, weight_numbers, weight_id, goflag );
  }
}




//----------------------------------------------------------------------------
void NCX30_CX41( const coord relnl, const double delta, int* weight_id, 
	size_t* goflag ) 
//----------------------------------------------------------------------------
{
  int weight_numbers[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};

  /* Stencil config, O is the cell containing the Lagrange multiplier */
  
  /* backward z plane: z = zlocal - Delta : inf1deltaz*/
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|x|x|x|_| --> |_|6|7|8|_| */ 
  /* |_|x|x|x|_| --> |_|5|4|3|_| */
  /* |_|x|x|x|_| --> |_|0|1|2|_| */
  
  if ( relnl.z < 1.1 * delta && relnl.z > 0.9 * delta ) 
  {
    fill_weight_numbers( -1, &weight_numbers[0] );
    NCX3_CX4( relnl, delta, weight_numbers, weight_id, goflag );
  }

  /* middle z plane: z = zlocal : sup0deltaz */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  /* |_|x|x|x|_| --> |_|15|16|17|_| */ 
  /* |_|x|O|x|_| --> |_|14|13|12|_| */
  /* |_|x|x|x|_| --> |_|_9|10|11|_| */
  
  else if ( fabs(relnl.z) < 0.1 * delta ) 
  {
    fill_weight_numbers( 0, &weight_numbers[0] );
    NCX3_CX4( relnl, delta, weight_numbers, weight_id, goflag );
  }
  
  /* forward z plane: z = zlocal + Delta : sup1deltaz */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  /* |_|x|x|x|_| --> |_|24|25|26|_| */ 
  /* |_|x|x|x|_| --> |_|23|22|21|_| */
  /* |_|x|x|x|_| --> |_|18|19|20|_| */
  
  else if ( relnl.z < - 0.9 * delta && relnl.z > - 1.1 * delta ) 
  {
    fill_weight_numbers( 1, &weight_numbers[0] );
    NCX3_CX4( relnl, delta, weight_numbers, weight_id, goflag );
  }
}




// NCX31_CX Bloc 3 functions

/* NCX31_CX40 = NCX30_CX41 */

//----------------------------------------------------------------------------
void NCX31_CX31( const coord relnl, const double delta, int* weight_id, 
	size_t* goflag ) 
//----------------------------------------------------------------------------
{
  int weight_numbers[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};
  
  /* Stencil config, O is the cell containing the Lagrange multiplier */
  
  /* backward z plane: z = zlocal - 2Delta : inf2deltaz */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |x|x|x|_|_| --> |6|7|8|_|_| */ 
  /* |x|x|x|_|_| --> |5|4|3|_|_| */
  /* |x|x|x|_|_| --> |0|1|2|_|_| */
  
  if ( relnl.z < 2.1 * delta && relnl.z > 1.9 * delta ) 
  {
    fill_weight_numbers( -1, &weight_numbers[0] );
    NCX3_CX3( relnl, delta, weight_numbers, weight_id, goflag );
  }

  /* middle z plane: z = zlocal - Delta : inf1deltaz */
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  /* |x|x|x|_|_| --> |15|16|17|_|_| */ 
  /* |x|x|x|_|_| --> |14|13|12|_|_| */
  /* |x|x|x|_|_| --> |_9|10|11|_|_| */
  
  else if ( relnl.z < 1.1 * delta && relnl.z > 0.9 * delta ) 
  {
    fill_weight_numbers( 0, &weight_numbers[0] );
    NCX3_CX3( relnl, delta, weight_numbers, weight_id, goflag );
  }
  
  // forward z plane: z = zlocal 
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  /* |_|_|_|_|_| --> |__|__|__|_|_| */
  /* |x|x|O|_|_| --> |24|25|26|_|_| */ 
  /* |x|x|x|_|_| --> |23|22|21|_|_| */
  /* |x|x|x|_|_| --> |18|19|20|_|_| */
  
  else if ( fabs(relnl.z) < 0.1 * delta ) 
  {
    fill_weight_numbers( 1, &weight_numbers[0] );
    NCX3_CX3( relnl, delta, weight_numbers, weight_id, goflag );
  }
}




//----------------------------------------------------------------------------
void NCX31_CX41( const coord relnl, const double delta, int* weight_id, 
	size_t* goflag ) 
//----------------------------------------------------------------------------
{
  int weight_numbers[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};
  
  /* Stencil config, O is the cell containing the Lagrange multiplier */
  
  /* backward z plane: z = zlocal - 2Delta : inf2deltaz */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|x|x|x|_| --> |_|6|7|8|_| */ 
  /* |_|x|x|x|_| --> |_|5|4|3|_| */
  /* |_|x|x|x|_| --> |_|0|1|2|_| */
  
  if ( relnl.z < 2.1 * delta && relnl.z > 1.9 * delta ) 
  {
    fill_weight_numbers(-1, &weight_numbers[0] );
    NCX3_CX4( relnl, delta, weight_numbers, weight_id, goflag );
  }

  /* middle z plane: z = zlocal - Delta : inf1deltaz */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  /* |_|x|x|x|_| --> |_|15|16|17|_| */ 
  /* |_|x|x|x|_| --> |_|14|13|12|_| */
  /* |_|x|x|x|_| --> |_|_9|10|11|_| */
  
  else if ( relnl.z < 1.1 * delta && relnl.z > 0.9 * delta ) 
  {
    fill_weight_numbers( 0, &weight_numbers[0] );
    NCX3_CX4( relnl, delta, weight_numbers, weight_id, goflag );
  }
  
  /* forward z plane: z = zlocal : sup0deltaz */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  /* |_|_|_|_|_| --> |_|__|__|__|_| */
  /* |_|x|x|x|_| --> |_|24|25|26|_| */ 
  /* |_|x|x|x|_| --> |_|23|22|21|_| */
  /* |_|x|x|x|_| --> |_|18|19|20|_| */
  
  else if ( fabs(relnl.z) < 0.1 * delta ) 
  {
    fill_weight_numbers( 1, &weight_numbers[0] );
    NCX3_CX4( relnl, delta, weight_numbers, weight_id, goflag );
  }
}




// NCX40_CX Bloc 2 functions

//----------------------------------------------------------------------------
void NCX40_CX40( const coord relnl, const double delta, int* weight_id, 
	size_t* goflag ) 
//----------------------------------------------------------------------------
{
  int weight_numbers[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};
  
  /* Stencil config, O is the cell containing the Lagrange multiplier */
  
  /* z = localplane (plane containing the multiplier) */
  /* backward z plane: z = zlocal : sup0deltaz */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|_|O|x|x| --> |_|_|6|7|8| */ 
  /* |_|_|x|x|x| --> |_|_|5|4|3| */
  /* |_|_|x|x|x| --> |_|_|0|1|2| */
  
  if ( fabs(relnl.z) < 0.1 * delta ) 
  {
    fill_weight_numbers( -1, &weight_numbers[0] );
    NCX4_CX4( relnl, delta, weight_numbers, weight_id, goflag );
  }

  // middle z plane: z = zlocal + Delta : sup1deltaz 
  /* |_|_|_|_|_| --> |_|_|__|__|__| */
  /* |_|_|_|_|_| --> |_|_|__|__|__| */
  /* |_|_|O|x|x| --> |_|_|15|16|17| */ 
  /* |_|_|x|x|x| --> |_|_|14|13|12| */
  /* |_|_|x|x|x| --> |_|_|_9|10|11| */
  
  else if ( relnl.z < - 0.9 * delta && relnl.z > - 1.1 * delta ) 
  {
    fill_weight_numbers( 0, &weight_numbers[0] );
    NCX4_CX4( relnl, delta, weight_numbers, weight_id, goflag );
  }
  
  /* forward z plane: z = zlocal + 2. * Delta : sup2deltaz */
  /* |_|_|_|_|_| --> |_|_|__|__|__| */
  /* |_|_|_|_|_| --> |_|_|__|__|__| */
  /* |_|_|O|x|x| --> |_|_|24|25|26| */ 
  /* |_|_|x|x|x| --> |_|_|23|22|21| */
  /* |_|_|x|x|x| --> |_|_|18|19|20| */
  
  else if ( relnl.z < - 1.9 * delta && relnl.z > - 2.1 * delta ) 
  {
    fill_weight_numbers( 1, &weight_numbers[0] );
    NCX4_CX4( relnl, delta, weight_numbers, weight_id, goflag );
  }
}




//----------------------------------------------------------------------------
void NCX40_CX41( const coord relnl, const double delta, int* weight_id, 
	size_t* goflag ) 
//----------------------------------------------------------------------------
{
  int weight_numbers[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};

  /* Stencil config, O is the cell containing the Lagrange multiplier */
  
  /* backward z plane: z = zlocal - Delta : inf1deltaz*/
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|_|x|x|x| --> |_|_|6|7|8| */ 
  /* |_|_|x|x|x| --> |_|_|5|4|3| */
  /* |_|_|x|x|x| --> |_|_|0|1|2| */
 
  if ( relnl.z < 1.1 * delta && relnl.z > 0.9 * delta ) 
  {
    fill_weight_numbers( -1, &weight_numbers[0] );
    NCX4_CX4( relnl, delta, weight_numbers, weight_id, goflag );
  }

  /* middle z plane: z = zlocal : sup0deltaz */
  /* |_|_|_|_|_| --> |_|_|__|__|__| */
  /* |_|_|_|_|_| --> |_|_|__|__|__| */
  /* |_|_|O|x|x| --> |_|_|15|16|17| */ 
  /* |_|_|x|x|x| --> |_|_|14|13|12| */
  /* |_|_|x|x|x| --> |_|_|_9|10|11| */
  
  else if ( fabs(relnl.z) < 0.1 * delta ) 
  {
    fill_weight_numbers( 0, &weight_numbers[0] );
    NCX4_CX4( relnl, delta, weight_numbers, weight_id, goflag );
  }
  
  // forward z plane: z = zlocal + Delta : sup1deltaz
  /* |_|_|_|_|_| --> |_|_|__|__|__| */
  /* |_|_|_|_|_| --> |_|_|__|__|__| */
  /* |_|_|x|x|x| --> |_|_|24|25|26| */ 
  /* |_|_|x|x|x| --> |_|_|23|22|21| */
  /* |_|_|x|x|x| --> |_|_|18|19|20| */
  
  else if ( relnl.z < - 0.9 * delta && relnl.z > - 1.1 * delta ) 
  {
    fill_weight_numbers( 1, &weight_numbers[0] );
    NCX4_CX4( relnl, delta, weight_numbers, weight_id, goflag );
  }
}




// NCX41_CX Bloc 1 function

//----------------------------------------------------------------------------
void NCX41_CX41( const coord relnl, const double delta, int* weight_id, 
	size_t* goflag ) 
//----------------------------------------------------------------------------
{
  int weight_numbers[9] = {-1, -1, -1, -1, -1, -1, -1, -1, -1};
  
  /* Stencil config, O is the cell containing the Lagrange multiplier */
  
  /* backward z plane: z = zlocal - 2Delta : inf2deltaz */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|_|_|_|_| --> |_|_|_|_|_| */
  /* |_|_|x|x|x| --> |_|_|6|7|8| */ 
  /* |_|_|x|x|x| --> |_|_|5|4|3| */
  /* |_|_|x|x|x| --> |_|_|0|1|2| */
  
  if ( relnl.z < 2.1 * delta && relnl.z > 1.9 * delta ) 
  {
    fill_weight_numbers( -1, &weight_numbers[0] );
    NCX4_CX4( relnl, delta, weight_numbers, weight_id, goflag );
  }

  /* middle z plane: z = zlocal - Delta : inf1deltaz */
  /* |_|_|_|_|_| --> |_|_|__|__|__| */
  /* |_|_|_|_|_| --> |_|_|__|__|__| */
  /* |_|_|x|x|x| --> |_|_|15|16|17| */ 
  /* |_|_|x|x|x| --> |_|_|14|13|12| */
  /* |_|_|x|x|x| --> |_|_|_9|10|11| */
  
  else if ( relnl.z < 1.1 * delta && relnl.z > 0.9 * delta ) 
  {
    fill_weight_numbers( 0, &weight_numbers[0] );
    NCX4_CX4( relnl, delta, weight_numbers, weight_id, goflag );
  }
  
  /* forward z plane: z = zlocal : sup0deltaz */
  /* |_|_|_|_|_| --> |_|_|__|__|__| */
  /* |_|_|_|_|_| --> |_|_|__|__|__| */
  /* |_|_|O|x|x| --> |_|_|24|25|26| */ 
  /* |_|_|x|x|x| --> |_|_|23|22|21| */
  /* |_|_|x|x|x| --> |_|_|18|19|20| */
  
  else if ( fabs(relnl.z) < 0.1 * delta ) 
  {
    fill_weight_numbers( 1, &weight_numbers[0] );
    NCX4_CX4( relnl, delta, weight_numbers, weight_id, goflag );
  }
}
#endif




//----------------------------------------------------------------------------
void assign_weight_id_quad_outward( const int NCX, const int CX, 
	const coord relnl, const double Delta, int* weight_id, 
	size_t* goflag ) 
//----------------------------------------------------------------------------
{
  double delta = Delta;
  
#if dimension == 2
  // Table taken from Peligriff
  /* (*Q2numb)(0,0) = 0 ;  */
  /* (*Q2numb)(1,0) = 4 ;  */
  /* (*Q2numb)(2,0) = 1 ;  */
  
  /* (*Q2numb)(0,1) = 7 ;  */
  /* (*Q2numb)(1,1) = 8 ;  */
  /* (*Q2numb)(2,1) = 5 ;  */
  
  /* (*Q2numb)(0,2) = 3 ;  */
  /* (*Q2numb)(1,2) = 6 ;  */
  /* (*Q2numb)(2,2) = 2 ;  */

  // weight_numbers[2] -> right column top
  // weight_numbers[5] -> right column middle
  // weight_numbers[1] -> right column bottom

  // weight_numbers[6] -> middle column top
  // weight_numbers[8] -> middle colum middle
  // weight_numbers[4] -> middle column bottom

  // weight_numbers[3] -> left column top
  // weight_numbers[7] -> left column middle
  // weight_numbers[0] -> left column bottom

  int weight_numbers[9];
  weight_numbers[3] = 3; weight_numbers[6] = 6; weight_numbers[2] = 2;
  weight_numbers[7] = 7; weight_numbers[8] = 8; weight_numbers[5] = 5;
  weight_numbers[0] = 0; weight_numbers[4] = 4; weight_numbers[1] = 1;

  switch( NCX )
  {
    // fictitous domain boundary is oriented (x+,y+)
    case 1:
      switch( CX )
      {
        case 1:
          // rel vector is oriented (x+,y+)
          NCX1_CX1( relnl, delta, weight_numbers, weight_id, goflag ); 	
	  break;
	  
        case 2:
          // rel vector is oriented (x-,y+)
          NCX1_CX2( relnl, delta, weight_numbers, weight_id, goflag );	
	  break;
	  
        case 3:
          // rel vector is oriented (x-,y-)
          NCX_centred( relnl, delta, weight_numbers, weight_id, goflag );	
	  break;
	  
        case 4:
          // rel vector is oriented (x-,y-)
          NCX1_CX4( relnl, delta, weight_numbers, weight_id, goflag );	
	  break;	  	  	  
      }
      break;
      
    // fictitous domain boundary is oriented (x-,y+)
    case 2:
      switch( CX )
      {
        case 1:
          // rel vector is oriented (x+,y+)
          // NCX2_CX1 = NCX1_CX2
          NCX1_CX2( relnl, delta, weight_numbers, weight_id, goflag );	
	  break;
	  
        case 2:
          // rel vector is oriented (x-,y+)
          NCX2_CX2( relnl, delta, weight_numbers, weight_id, goflag );	
	  break;
	  
        case 3:
          // rel vector is oriented (x-,y-)
          NCX2_CX3( relnl, delta, weight_numbers, weight_id, goflag );	
	  break;
	  
        case 4:
          // rel vector is oriented (x-,y-)
          NCX_centred( relnl, delta, weight_numbers, weight_id, goflag );	
	  break;	  	  	  	  	  	  
      }
      break;
      
    // fictitous domain boundary is oriented (x-,y-)
    case 3:
      switch( CX )
      {
        case 1:
          // rel vector is oriented (x+,y+)
          // NCX3_CX1 = NCX1_CX3 = NCX_centred
          NCX_centred( relnl, delta, weight_numbers, weight_id, goflag );	
	  break;
	  
        case 2:
          // rel vector is oriented (x-,y+)
          // NCX3_CX2 = NCX2_CX3
          NCX2_CX3(relnl, delta, weight_numbers, weight_id, goflag );	
	  break;
	  
        case 3:
          // rel vector is oriented (x-,y-)
          NCX3_CX3( relnl, delta, weight_numbers, weight_id, goflag );	
	  break;
	  
        case 4:
          // rel vector is oriented (x+,y-)
          NCX3_CX4( relnl, delta, weight_numbers, weight_id, goflag );	
	  break;
      }
      break;
      
    // fictitous domain boundary is oriented (x+,y-)
    case 4:
      switch( CX )
      {
        case 1:
          // rel vector is oriented (x+,y+)
          // NCX4_CX1 = NCX1_CX4
          NCX1_CX4( relnl, delta, weight_numbers, weight_id, goflag );	
	  break;
	  
        case 2:
          // rel vector is oriented (x-,y+)
          // NCX4_CX2 = NCX2_CX4 = NCX_centred
          NCX_centred( relnl, delta, weight_numbers, weight_id, goflag );	
	  break;
	  
        case 3:
          // rel vector is oriented (x-,y-)
          // NCX4_CX3 = NCX3_CX4
          NCX3_CX4( relnl, delta, weight_numbers, weight_id, goflag );	
	  break;
	  
        case 4:
          // rel vector is oriented (x+,y-)
          NCX4_CX4( relnl, delta, weight_numbers, weight_id, goflag );	
	  break;
      }
      break;                     
  }
#endif
  
#if dimension == 3
  switch( NCX )
  {
    // fictitous domain boundary is oriented (x+,y+,z+)
    case 10:
      switch( CX )
      {
        case 10:
          // rel vector is oriented (x+,y+,z+)
          NCX10_CX10( relnl, delta, weight_id, goflag );
	  break;
	  
        case 20:
          // rel vector is oriented (x-,y+,z+)
          NCX10_CX20( relnl, delta, weight_id, goflag );	
	  break;
	  
        case 30:
          // rel vector is oriented (x-,y-,z+)
          NCX10_CX30( relnl, delta, weight_id, goflag );
	  break;
	  
        case 40:
          // rel vector is oriented (x+,y-,z+)
          NCX10_CX40( relnl, delta, weight_id, goflag );
	  break;
	  
        case 11:
          // rel vector is oriented (x+,y+,z-)
          NCX10_CX11( relnl, delta, weight_id, goflag );
	  break;
	  
        case 21:
          // rel vector is oriented (x-,y+,z-)
          NCX10_CX21( relnl, delta, weight_id, goflag );
	  break;
	  
        case 31:
          // rel vector is oriented (x-,y-,z-)
          NCX10_CX31( relnl, delta, weight_id, goflag );
	  break;
	  
        case 41:
          // rel vector is oriented (x+,y-,z-)
          NCX10_CX41( relnl, delta, weight_id, goflag );
	  break;
      }
      break;
      
    // fictitous domain boundary is oriented (x+,y+,z-)
    case 11:
      switch( CX )
      {
        case 10:
          // rel vector is oriented (x+,y+,z+)
          // NCX11_CX10 = NCX10_CX11
          NCX10_CX11( relnl, delta, weight_id, goflag );
	  break;
	  
        case 20:
          // rel vector is oriented (x-,y+,z+)
          // NCX11_CX20 = NCX10_CX21
          NCX10_CX21( relnl, delta, weight_id, goflag );
	  break;

        case 30:
          // rel vector is oriented (x-,y-,z+)
          // NCX11_CX30 = NCX10_CX31
          NCX10_CX31( relnl, delta, weight_id, goflag );
	  break;	  

        case 40:
          // rel vector is oriented (x+,y-,z+)
          // NCX11_CX40 = NCX10_CX41
          NCX10_CX41( relnl, delta, weight_id, goflag );
	  break;
	  
        case 11:
          // rel vector is oriented (x+,y+,z-)
          NCX11_CX11( relnl, delta, weight_id, goflag );
	  break;
	  
        case 21:
          // rel vector is oriented (x-,y+,z-)
          NCX11_CX21( relnl, delta, weight_id, goflag );
	  break;
	  
        case 31:
          // rel vector is oriented (x-,y-,z-)
          NCX11_CX31( relnl, delta, weight_id, goflag );
	  break;
	  
        case 41:
          // rel vector is oriented (x+,y-,z-)
          NCX11_CX41( relnl, delta, weight_id, goflag );
	  break;
      }
      break;
      
    // fictitous domain boundary is oriented (x-,y+,z+)
    case 20:
      switch( CX )
      {
        case 10:
          // rel vector is oriented (x+,y+,z+)
          // NCX20_CX10 = NCX10_CX20
          NCX10_CX20( relnl, delta, weight_id, goflag );          
	  break;
	  
        case 20:
          // rel vector is oriented (x-,y+,z+)
          NCX20_CX20( relnl, delta, weight_id, goflag );
	  break;
	  
        case 30:
          // rel vector is oriented (x-,y-,z+)
          NCX20_CX30( relnl, delta, weight_id, goflag );
	  break;
	  
        case 40:
          // rel vector is oriented (x+,y-,z+)
          NCX20_CX40( relnl, delta, weight_id, goflag );
	  break;
	  
        case 11:
          // rel vector is oriented (x+,y+,z-)
          // NCX20_CX11 = NCX11_CX20 = NCX10_CX21
          NCX10_CX21( relnl, delta, weight_id, goflag );
	  break;
	  
        case 21:
          // rel vector is oriented (x-,y+,z-)
          NCX20_CX21( relnl, delta, weight_id, goflag );          
	  break;
	  
        case 31:
          // rel vector is oriented (x-,y-,z-)
          NCX20_CX31( relnl, delta, weight_id, goflag );          
	  break;
	  
        case 41:
          // rel vector is oriented (x+,y-,z-)
          NCX20_CX41( relnl, delta, weight_id, goflag );
	  break;
      }
      break;            

    // fictitous domain boundary is oriented (x-,y+,z-)
    case 21:
      switch( CX )
      {
        case 10:
          // rel vector is oriented (x+,y+,z+)
          // NCX21_CX10 = NCX10_CX21
          NCX10_CX21( relnl, delta, weight_id, goflag );
	  break;
	  
        case 20:
          // rel vector is oriented (x-,y+,z+)
          // NCX21_CX20 = NCX20_CX21
          NCX20_CX21( relnl, delta, weight_id, goflag );
	  break;
	  
        case 30:
          // rel vector is oriented (x-,y-,z+)
          // NCX21_CX30 = NCX20_CX31
          NCX20_CX31( relnl, delta, weight_id, goflag );
	  break;
	  
        case 40:
          // rel vector is oriented (x+,y-,z+)
          // NCX21_CX40 = NCX20_CX41
          NCX20_CX41( relnl, delta, weight_id, goflag );          
	  break;
	  
        case 11:
          // rel vector is oriented (x+,y+,z-)
          // NCX21_CX11 = NCX11_CX21
          NCX11_CX21( relnl, delta, weight_id, goflag );
	  break;
	  
        case 21:
          // rel vector is oriented (x-,y+,z-)
          NCX21_CX21( relnl, delta, weight_id, goflag );
	  break;
	  
        case 31:
          // rel vector is oriented (x-,y-,z-)
          NCX21_CX31( relnl, delta, weight_id, goflag );
	  break;
	  
        case 41:
          // rel vector is oriented (x+,y-,z-)
          NCX21_CX41( relnl, delta, weight_id, goflag );
	  break;	  	  	  	  	  	  	  
      }
      break;
      
    // fictitous domain boundary is oriented (x-,y-,z+)
    case 30:
      switch( CX )
      {
        case 10:
          // rel vector is oriented (x+,y+,z+)
          // NCX30_CX10 = NCX10_CX30
          NCX10_CX30( relnl, delta, weight_id, goflag );
	  break;
	  
        case 20:
          // rel vector is oriented (x-,y+,z+)
          // NCX30_CX20 = NCX20_CX30
          NCX20_CX30( relnl, delta, weight_id, goflag );          
	  break;
	  
        case 30:
          // rel vector is oriented (x-,y-,z+)
          NCX30_CX30( relnl, delta, weight_id, goflag );          
	  break;
	  
        case 40:
          // rel vector is oriented (x+,y-,z+)
          NCX30_CX40( relnl, delta, weight_id, goflag );          
	  break;
	  
        case 11:
          // rel vector is oriented (x+,y+,z-)
          // NCX30_CX11 = NCX11_CX30 = NCX10_CX31
          NCX10_CX31( relnl, delta, weight_id, goflag );          
	  break;
	  
        case 21:
          // rel vector is oriented (x-,y+,z-)
          // NCX30_CX21 = NCX21_CX30 = NCX20_CX31
          NCX20_CX31( relnl, delta, weight_id, goflag );
	  break;
	  
        case 31:
          // rel vector is oriented (x-,y-,z-)
          NCX30_CX31( relnl, delta,  weight_id, goflag );
	  break;
	  
        case 41:
          // rel vector is oriented (x+,y-,z-)
          NCX30_CX41( relnl, delta, weight_id, goflag );
	  break;
      }
      break;
      
    // fictitous domain boundary is oriented (x-,y-,z-)
    case 31:
      switch( CX )
      {
        case 10:
          // rel vector is oriented (x+,y+,z+)
          // NCX31_CX10 = NCX10_CX31
          NCX10_CX31( relnl, delta, weight_id, goflag );
	  break;
	  
        case 20:
          // rel vector is oriented (x-,y+,z+)
          // NCX31_CX20  = NCX20_CX31
          NCX20_CX31( relnl, delta, weight_id, goflag );
	  break;
	  
        case 30:
          // rel vector is oriented (x-,y-,z+)
          // NCX31_CX30 = NCX30_CX31
          NCX30_CX31( relnl, delta, weight_id, goflag );
	  break;
	  
        case 40:
          // rel vector is oriented (x+,y-,z+)
          // NCX31_CX40 = NCX30_CX41
          NCX30_CX41( relnl, delta, weight_id, goflag );
	  break;
	  
        case 11:
          // rel vector is oriented (x+,y+,z-)
          // NCX31_CX11 = NCX11_CX31
          NCX11_CX31( relnl, delta, weight_id, goflag );
	  break;
	  
        case 21:
          // rel vector is oriented (x-,y+,z-)
          // NCX31_CX21 = NCX21_CX31
          NCX21_CX31( relnl, delta, weight_id, goflag );
	  break;
	  
        case 31:
          // rel vector is oriented (x-,y-,z-)
          NCX31_CX31( relnl, delta, weight_id, goflag );
	  break;
	  
        case 41:
          // rel vector is oriented (x+,y-,z-)
          NCX31_CX41( relnl, delta, weight_id, goflag );
	  break;
      }
      break;
      
    // fictitous domain boundary is oriented (x+,y-,z+)
    case 40:
      switch( CX )
      {
        case 10:
          // rel vector is oriented (x+,y+,z+)
          // NCX40_CX10 = NCX10_CX40
          NCX10_CX40( relnl, delta, weight_id, goflag );
	  break;
	  
        case 20:
          // rel vector is oriented (x-,y+,z+)
          // NCX40_CX20 = NCX20_CX40
          NCX20_CX40( relnl, delta, weight_id, goflag );
	  break;
	  
        case 30:
          // rel vector is oriented (x-,y-,z+)
          // NCX40_CX30 = NCX30_CX40
          NCX30_CX40( relnl, delta, weight_id, goflag );
	  break;
	  
        case 40:
          // rel vector is oriented (x+,y-,z+)
          NCX40_CX40( relnl, delta, weight_id, goflag );
	  break;
	  
        case 11:
          // rel vector is oriented (x+,y+,z-)
          // NCX40_CX11 = NCX11_CX40 = NCX10_CX41
          NCX10_CX41( relnl, delta, weight_id, goflag );
	  break;
	  
        case 21:
          // rel vector is oriented (x-,y+,z-)
          // NCX40_CX21 = NCX21_CX40 = NCX20_CX41
          NCX20_CX41( relnl, delta, weight_id, goflag );
	  break;
	  
        case 31:
          // rel vector is oriented (x-,y-,z-)
          // NCX40_CX31 = NCX31_CX40 = NCX30_CX41
          NCX30_CX41( relnl, delta, weight_id, goflag );
	  break;
	  
        case 41:
          // rel vector is oriented (x+,y-,z-)
          NCX40_CX41( relnl, delta, weight_id, goflag );
	  break;
      }
      break;
      
    // fictitous domain boundary is oriented (x+,y-,z-)
    case 41:
      switch( CX )
      {
        case 10:
          // rel vector is oriented (x+,y+,z+)
          // NCX41_CX10 = NCX10_CX41
          NCX10_CX41( relnl, delta, weight_id, goflag );
	  break;
	  
        case 20:
          // rel vector is oriented (x-,y+,z+)
          // NCX41_CX20 = NCX20_CX41
          NCX20_CX41( relnl, delta, weight_id, goflag );
	  break;
	  
        case 30:
          // rel vector is oriented (x-,y-,z+)
          // NCX41_CX30 = NCX30_CX41
          NCX30_CX41( relnl, delta, weight_id, goflag );
	  break;
	  
        case 40:
          // rel vector is oriented (x+,y-,z+)
          // NCX41_CX40 = NCX40_CX41
          NCX40_CX41( relnl, delta, weight_id, goflag );
	  break;
	  
        case 11:
          // rel vector is oriented (x+,y+,z-)
          // NCX41_CX11 = NCX11_CX41
          NCX11_CX41( relnl, delta, weight_id, goflag );
	  break;
	  
        case 21:
          // rel vector is oriented (x-,y+,z-)
          // NCX41_CX21 = NCX21_CX41
          NCX21_CX41( relnl, delta, weight_id, goflag );
	  break;
	  
        case 31:
          // rel vector is oriented (x-,y-,z-)
          // NCX41_CX31 = NCX31_CX41
          NCX31_CX41(relnl, delta, weight_id, goflag );
	  break;
	  
        case 41:
          // rel vector is oriented (x+,y-,z-)
          NCX41_CX41( relnl, delta, weight_id, goflag );
	  break;
      }
      break;
  }
#endif
}




//----------------------------------------------------------------------------
double Q2weighting( size_t  i, double  x, double y, double z ) 
//----------------------------------------------------------------------------
{ 
  double result = 0. ;

#if dimension == 2
  
  switch( i )
  {
      case 0 :
         result = (1.0-x) * (1.0-2.0*x) * (1.0-y) * (1.0-2.0*y) ;
         break ;
      case 1 :
         result = x * (2.0*x-1.0) * (1.0-y) * (1.0-2.0*y) ;
         break ;
      case 2 :
         result = x * (2.0*x-1.0) * y * (2.0*y-1.0) ;
         break ;
      case 3 :
         result = (1.0-x) * (1.0-2.0*x) * y * (2.0*y-1.0) ;
         break ;
      case 4 :
         result = 4.0 * x * (1.0-x) * (1.0-y) * (1.0-2.0*y) ;
         break ;
      case 5 :
         result = x * (2.0*x-1.0) * 4.0 * y * (1.0-y) ;
         break ;
      case 6 :
         result = 4.0 * x * (1.0-x) * y * (2.0*y-1.0) ;
         break ;
      case 7 :
         result = (1.0-x) * (1.0-2.0*x) * 4.0 * y * (1.0-y) ;
         break ;
      case 8 :
         result = 4.0 * x * (1.0-x) * 4.0 * y * (1.0-y) ;
         break ;
  }
#endif

#if dimension == 3
  
  switch( i )
  {
      case 0 :
         result = 8.0 * (1.0-x)*(0.5-x) * (1.0-y)*(0.5-y) * (1.0-z)*(0.5-z);
         break ;
      case 1 :
         result = 16.0 * x*(1.0-x) * (1.0-y)*(0.5-y) * (1.0-z)*(0.5-z);
         break ;
      case 2 :
         result = -8.0 * x*(0.5-x) * (1.0-y)*(0.5-y) * (1.0-z)*(0.5-z) ;
         break ;
      case 3 :
         result = -16.0 * x*(0.5-x) * y*(1.0-y) * (1.0-z)*(0.5-z);
         break ;
      case 4 :
         result = 32.0 * x*(1.0-x) * y*(1.0-y) * (1.0-z)*(0.5-z) ;
         break ;
      case 5 :
         result = 16.0 * (1.0-x)*(0.5-x) * y*(1.0-y) * (1.0-z)*(0.5-z) ;
         break ;
      case 6 :
         result =-8.0 * (1.0-x)*(0.5-x) * y*(0.5-y) * (1.0-z)*(0.5-z) ;
         break ;
      case 7 :
         result = -16.0 * x*(1.0-x) * y*(0.5-y) * (1.0-z)*(0.5-z);
         break ;
      case 8 :
         result = 8.0 * x*(0.5-x) * y*(0.5-y) * (1.0-z)*(0.5-z);
         break ;
      case 9 :
         result = 16.0 * (1.0-x)*(0.5-x) * (1.0-y)*(0.5-y) * z*(1.0-z) ;
         break ;
      case 10 :
         result = 32.0 * x*(1.0-x) * (1.0-y)*(0.5-y) *  z*(1.0-z);
         break ;
      case 11 :
         result = -16.0 * x*(0.5-x) * (1.0-y)*(0.5-y) *  z*(1.0-z) ;
         break ;
      case 12 :
         result = -32.0 * x*(0.5-x) * y*(1.0-y) * z*(1.0-z) ;
         break ;
      case 13 :
         result = 64.0 * x*(1.0-x) * y*(1.0-y) *  z*(1.0-z) ;
         break ;
      case 14 :
         result = 32.0 * (1.0-x)*(0.5-x) * y*(1.0-y) *  z*(1.0-z) ;
         break ;
      case 15 :
         result =-16.0 * (1.0-x)*(0.5-x) * y*(0.5-y) *  z*(1.0-z) ;
         break ;
      case 16 :
         result = -32.0 * x*(1.0-x) * y*(0.5-y) *  z*(1.0-z) ;
         break ;
      case 17 :
         result = 16.0 * x*(0.5-x) * y*(0.5-y) * z*(1.0-z);
         break ;
      case 18 :
         result = -8.0 * (1.0-x)*(0.5-x) * (1.0-y)*(0.5-y) * z*(0.5-z) ;
         break ;
      case 19 :
         result = -16.0 * x*(1.0-x) * (1.0-y)*(0.5-y) *  z*(0.5-z);
         break ;
      case 20 :
         result = 8.0 * x*(0.5-x) * (1.0-y)*(0.5-y) *  z*(0.5-z) ;
         break ;
      case 21 :
         result = 16.0 * x*(0.5-x) * y*(1.0-y) * z*(0.5-z) ;
         break ;
      case 22 :
         result = -32.0 * x*(1.0-x) * y*(1.0-y) *  z*(0.5-z) ;
         break ;
      case 23 :
         result = -16.0 * (1.0-x)*(0.5-x) * y*(1.0-y) *  z*(0.5-z) ;
         break ;
      case 24 :
         result = 8.0 * (1.0-x)*(0.5-x) * y*(0.5-y) *  z*(0.5-z) ;
         break ;
      case 25 :
         result = 16.0 * x*(1.0-x) * y*(0.5-y) *  z*(0.5-z) ;
         break ;
      case 26 :
         result = -8.0 * x*(0.5-x) * y*(0.5-y) * z*(0.5-z);
         break ; 
  }
   
#endif
  return(result);  
}




//----------------------------------------------------------------------------
void weight_NCX1_CX1( const coord poslocal, const double Delta, double* x1, 
	double* y1) 
//----------------------------------------------------------------------------
{
  // local cell is the bottom left cell of the stencil
  *x1 = poslocal.x;
  *y1 = poslocal.y;
}




//----------------------------------------------------------------------------
void weight_NCX1_CX2( const coord poslocal, const double Delta, double* x1, 
	double* y1) 
//----------------------------------------------------------------------------
{
  // local cell is the bottom center cell of the stencil
  *x1 = poslocal.x - Delta;
  *y1 = poslocal.y;
}




//----------------------------------------------------------------------------
void weight_centred( const coord poslocal, const double Delta, double* x1, 
	double* y1) 
//----------------------------------------------------------------------------
{
  // local cell is the center cell of the stencil
  *x1 = poslocal.x - Delta;
  *y1 = poslocal.y - Delta;
}




//----------------------------------------------------------------------------
void weight_NCX1_CX4( const coord poslocal, const double Delta, double* x1, 
	double* y1) 
//----------------------------------------------------------------------------
{
  // local cell is the center left cell of the stencil
  *x1 = poslocal.x;
  *y1 = poslocal.y - Delta;
}




//----------------------------------------------------------------------------
void weight_NCX2_CX2( const coord poslocal, const double Delta, double* x1, 
	double* y1) 
//----------------------------------------------------------------------------
{
  // local cell is the bottom right cell of the stencil
  *x1 = poslocal.x - 2. * Delta;
  *y1 = poslocal.y;
}




//----------------------------------------------------------------------------
void weight_NCX2_CX3( const coord poslocal, const double Delta, double* x1, 
	double* y1) 
//----------------------------------------------------------------------------
{
  // local cell is the center right cell of the stencil
  *x1 = poslocal.x - 2. * Delta;
  *y1 = poslocal.y - Delta;
}




//----------------------------------------------------------------------------
void weight_NCX3_CX3( const coord poslocal, const double Delta, double* x1, 
	double* y1) 
//----------------------------------------------------------------------------
{
  // local cell is the top right cell of the stencil
  *x1 = poslocal.x - 2. * Delta; 
  *y1 = poslocal.y - 2. * Delta; 
}




//----------------------------------------------------------------------------
void weight_NCX3_CX4( const coord poslocal, const double Delta, double* x1, 
	double* y1) 
//----------------------------------------------------------------------------
{
  // local cell is the top center cell of the stencil
  *x1 = poslocal.x - Delta; 
  *y1 = poslocal.y - 2. * Delta; 
}




//----------------------------------------------------------------------------
void weight_NCX4_CX4( const coord poslocal, const double Delta, double* x1, 
	double* y1) 
//----------------------------------------------------------------------------
{
  // local cell is at the top left cell of the stencil 
  *x1 = poslocal.x; 
  *y1 = poslocal.y - 2. * Delta;
}




//----------------------------------------------------------------------------
double compute_weight_Quad( const int weight_id, const coord posb ,
	const coord poslocal, const int NCX, const int CX, const double Delta ) 
//----------------------------------------------------------------------------
{
  double x1 = 0., y1 = 0., weight = 0.;
  coord posref;
  
#if dimension == 2
  switch( NCX )
  {
    case 1:
      switch( CX )
      {
        case 1;
	  weight_NCX1_CX1( poslocal, Delta, &x1, &y1 ); 
	  break;
	  
        case 2;
	  weight_NCX1_CX2( poslocal, Delta, &x1, &y1 );
	  break;	  
	  
        case 3;
	  weight_centred( poslocal, Delta, &x1, &y1 ); 
	  break;

        case 4;
	  weight_NCX1_CX4( poslocal, Delta, &x1, &y1 );
	  break;	  	  
      }
      break;
      
    case 2:
      switch( CX )
      {
        case 1;
	  weight_NCX1_CX2( poslocal, Delta, &x1, &y1 ); 
	  break;
	  
        case 2;
	  weight_NCX2_CX2( poslocal, Delta, &x1, &y1 );
	  break;	  
	  
        case 3;
	  weight_NCX2_CX3( poslocal, Delta, &x1, &y1 );
	  break;

        case 4;
	  weight_centred( poslocal, Delta, &x1, &y1 );
	  break;	  	  
      }        
      break;      
  
    case 3:
      switch( CX )
      {
        case 1;
	  weight_centred( poslocal, Delta, &x1, &y1 ); 
	  break;
	  
        case 2;
	  weight_NCX2_CX3( poslocal, Delta, &x1, &y1 );
	  break;	  
	  
        case 3;
	  weight_NCX3_CX3( poslocal, Delta, &x1, &y1 );
	  break;

        case 4;
	  weight_NCX3_CX4( poslocal, Delta, &x1, &y1 );
	  break;	  	  
      }           
      break;  
  
    case 4:
      switch( CX )
      {
        case 1;
	  weight_NCX1_CX4( poslocal, Delta, &x1, &y1 );
	  break;
	  
        case 2;
	  weight_centred( poslocal, Delta, &x1, &y1 );
	  break;	  
	  
        case 3;
	  weight_NCX3_CX4( poslocal, Delta, &x1, &y1 );
	  break;

        case 4;
	  weight_NCX4_CX4( poslocal, Delta, &x1, &y1 );
	  break;	  	  
      }           
      break;    
  }
  
  posref.x = GENERAL_1DIST( posb.x, x1 ) / ( 2. * Delta );
  posref.y = GENERAL_1DIST( posb.y, y1 ) / ( 2. * Delta );
    
  return weight = Q2weighting( weight_id, posref.x, posref.y, 0. );
  
#elif dimension == 3
  double z1 = 0.;
  
  switch( NCX )
  {
    case 10:
      switch( CX )
      {
        case 10: 
	  // NCX10_CX10
          z1 = poslocal.z;
          weight_NCX1_CX1( poslocal, Delta, &x1, &y1 );
	  break;
	  
        case 20: 
          // NCX10_CX20
          z1 = poslocal.z;
          weight_NCX1_CX2( poslocal, Delta, &x1, &y1 );
	  break;	  

        case 30: 
          // NCX10_CX30
          z1 = poslocal.z;
          weight_centred( poslocal, Delta, &x1, &y1 );
	  break;

        case 40: 
          // NCX10_CX40
          z1 = poslocal.z;
          weight_NCX1_CX4( poslocal, Delta, &x1, &y1 );
	  break;
	  
        case 11: 
          // NCX10_CX11
          z1 = poslocal.z - Delta;
          weight_NCX1_CX1( poslocal, Delta, &x1, &y1 );
	  break;
	  
        case 21: 
          // NCX10_CX21
          z1 = poslocal.z - Delta;
          weight_NCX1_CX2( poslocal, Delta, &x1, &y1 );
	  break;
	  
        case 31: 
          // NCX10_CX31
          z1 = poslocal.z - Delta;
          weight_centred( poslocal, Delta, &x1, &y1 );
	  break;	  	  	

        case 41: 
          // NCX10_CX41
          z1 = poslocal.z - Delta;
          weight_NCX1_CX4( poslocal, Delta, &x1, &y1 );
	  break;      
      }
      break;
        
    case 11:
      switch( CX )
      {
        case 10: 
          // NCX11_CX10 = NCX10_CX11 
          z1 = poslocal.z - Delta;
          weight_NCX1_CX1( poslocal, Delta, &x1, &y1 );
	  break; 
	  
        case 20: 
          // NCX11_CX20
          z1 = poslocal.z - Delta;
          weight_NCX1_CX2( poslocal, Delta, &x1, &y1 );
	  break;
	  
        case 30: 
          // NCX11_CX30
          z1 = poslocal.z - Delta;
          weight_centred( poslocal, Delta, &x1, &y1 );
	  break; 	   	  

        case 40: 
          // NCX11_CX40
          z1 = poslocal.z - Delta;
          weight_NCX1_CX4( poslocal, Delta, &x1, &y1 );
	  break; 

        case 11: 
          // NCX11_CX11
          z1 = poslocal.z - 2. * Delta; 
          weight_NCX1_CX1( poslocal, Delta, &x1, &y1 );
	  break; 	  

        case 21: 
          // NCX11_CX21
          z1 = poslocal.z - 2. * Delta;
          weight_NCX1_CX2( poslocal, Delta, &x1, &y1 );
	  break; 
	  
        case 31: 
          // NCX11_CX31
          z1 = poslocal.z - 2. * Delta;
          weight_centred( poslocal, Delta, &x1, &y1 );
	  break;
	  
        case 41: 
          // NCX11_CX41
          z1 = poslocal.z - 2. * Delta;
          weight_NCX1_CX4( poslocal, Delta, &x1, &y1 );
	  break; 	   	  	       
      }
      break;  

    case 20:
      switch( CX )
      {
        case 10: 
          // NCX20_CX10 = NCX10_CX20
          z1 = poslocal.z;
          weight_NCX1_CX2( poslocal, Delta, &x1, &y1 );
	  break; 
	  
        case 20: 
          // NCX20_CX20
          z1 = poslocal.z;
          weight_NCX2_CX2( poslocal, Delta, &x1, &y1 );
	  break; 
	  
        case 30: 
          // NCX20_CX30
          z1 = poslocal.z;
          weight_NCX2_CX3( poslocal, Delta, &x1, &y1 );
	  break; 
	  
        case 40: 
          // NCX20_CX40
          z1 = poslocal.z;
          weight_centred( poslocal, Delta, &x1, &y1 );
	  break; 
	  
        case 11: 
          // NCX20_CX11 = NCX11_CX20
          z1 = poslocal.z - Delta;
          weight_NCX1_CX2( poslocal, Delta, &x1, &y1 );
	  break;
	  
        case 21: 
          // NCX20_CX21
          z1 = poslocal.z - Delta;
          weight_NCX2_CX2( poslocal, Delta, &x1, &y1 );
	  break; 
	  
        case 31: 
          // NCX20_CX31
          z1 = poslocal.z - Delta;
          weight_NCX2_CX3( poslocal, Delta, &x1, &y1 );
	  break;
	  
        case 41: 
          // NCX20_CX41
          z1 = poslocal.z - Delta;
          weight_centred( poslocal, Delta, &x1, &y1 );
	  break; 	   	  	   	  	  	  	        
      }
      break;

    case 21:
      switch( CX )
      {
        case 10: 
          // NCX21_CX10 = NCX10_CX21
          z1 = poslocal.z - Delta;
          weight_NCX1_CX2( poslocal, Delta, &x1, &y1 );
	  break; 
	  
        case 20: 
          // NCX21_CX20 = NCX20_CX21
          z1 = poslocal.z - Delta;
          weight_NCX2_CX2( poslocal, Delta, &x1, &y1 );
	  break;
	  
        case 30: 
          // NCX21_CX30
          z1 = poslocal.z - Delta;
          weight_NCX2_CX3( poslocal, Delta, &x1, &y1 );
	  break; 
	  
        case 40: 
          // NCX21_CX40
          z1 = poslocal.z - Delta;
          weight_centred( poslocal, Delta, &x1, &y1 );
	  break; 
	  
        case 11: 
          // NCX21_CX11
          z1 = poslocal.z - 2. * Delta;
          weight_NCX1_CX2( poslocal, Delta, &x1, &y1 );
	  break;
	  
        case 21: 
          // NCX21_CX21
          z1 = poslocal.z - 2. * Delta;
          weight_NCX2_CX2( poslocal, Delta, &x1, &y1 );
	  break; 
	  
        case 31: 
          // NCX21_CX31
          z1 = poslocal.z - 2. * Delta;
          weight_NCX2_CX3( poslocal, Delta, &x1, &y1 );
	  break;
	  
        case 41: 
          // NCX21_CX41
          z1 = poslocal.z - 2. * Delta;
          weight_centred( poslocal, Delta, &x1, &y1 );
	  break;
      }
      break;

    case 30:
      switch( CX )
      {	  
        case 10: 
          // NCX30_CX10
          z1 = poslocal.z;
          weight_centred( poslocal, Delta, &x1, &y1 );
	  break;
	  
        case 20: 
          // NCX30_CX20
          z1 = poslocal.z;
          weight_NCX2_CX3( poslocal, Delta, &x1, &y1 );
	  break;
	  
        case 30: 
          // NCX30_CX30
          z1 = poslocal.z;
          weight_NCX3_CX3( poslocal, Delta, &x1, &y1 );
	  break; 
	  
        case 40: 
          // NCX30_CX40
          z1 = poslocal.z;
          weight_NCX3_CX4( poslocal, Delta, &x1, &y1 );
	  break; 
	  
        case 11: 
          // NCX30_CX11
          z1 = poslocal.z - Delta;
          weight_centred( poslocal, Delta, &x1, &y1 );
	  break;
	  
        case 21: 
          // NCX30_CX21
          z1 = poslocal.z - Delta;
          weight_NCX2_CX3( poslocal, Delta, &x1, &y1 );
	  break;
	  
        case 31: 
          // NCX30_CX31
          z1 = poslocal.z - Delta;
          weight_NCX3_CX3( poslocal, Delta, &x1, &y1 );
	  break; 
	  
        case 41: 
          // NCX30_CX41
          z1 = poslocal.z - Delta;
          weight_NCX3_CX4( poslocal, Delta, &x1, &y1 );
	  break;
      }
      break;   

    case 31:
      switch( CX )
      {	  
        case 10: 
          // NCX31_CX10
          z1 = poslocal.z - Delta;
          weight_centred( poslocal, Delta, &x1, &y1 );
	  break;
	  
        case 20: 
          // NCX31_CX20
          z1 = poslocal.z - Delta;
          weight_NCX2_CX3( poslocal, Delta, &x1, &y1 );
	  break;
	  
        case 30: 
          // NCX31_CX30
          z1 = poslocal.z - Delta;
          weight_NCX3_CX3( poslocal, Delta, &x1, &y1 );
	  break; 
	  
        case 40: 
          // NCX31_CX40
          z1 = poslocal.z - Delta;
          weight_NCX3_CX4( poslocal, Delta, &x1, &y1 );
	  break; 
	  
        case 11: 
          // NCX31_CX11
          z1 = poslocal.z - 2. * Delta;
          weight_centred( poslocal, Delta, &x1, &y1 );
	  break;
	  
        case 21: 
          // NCX31_CX21
          z1 = poslocal.z - 2. * Delta;
          weight_NCX2_CX3( poslocal, Delta, &x1, &y1 );
	  break; 
	  
        case 31: 
          // NCX31_CX31
          z1 = poslocal.z - 2. * Delta;
          weight_NCX3_CX3( poslocal, Delta, &x1, &y1 );
	  break;
	  
        case 41: 
          // NCX31_CX41
          z1 = poslocal.z - 2. * Delta;
          weight_NCX3_CX4( poslocal, Delta, &x1, &y1 );
	  break;
      }
      break;

    case 40:
      switch( CX )
      {	  
        case 10: 
          // NCX40_CX10
          z1 = poslocal.z;
          weight_NCX1_CX4( poslocal, Delta, &x1, &y1 );
	  break; 
	  
        case 20: 
          // NCX40_CX20
          z1 = poslocal.z;
          weight_centred( poslocal, Delta, &x1, &y1 );
	  break;
	  
        case 30: 
          // NCX40_CX30
          z1 = poslocal.z;
          weight_NCX3_CX4( poslocal, Delta, &x1, &y1 );
	  break; 
	  
        case 40: 
          // NCX40_CX40
          z1 = poslocal.z;
          weight_NCX4_CX4( poslocal, Delta, &x1, &y1 );
	  break;
	  
        case 11: 
          // NCX40_CX11
          z1 = poslocal.z - Delta;
          weight_NCX1_CX4( poslocal, Delta, &x1, &y1 );
	  break;
	  
        case 21: 
          // NCX40_CX21
          z1 = poslocal.z - Delta;
          weight_centred( poslocal, Delta, &x1, &y1 );
	  break; 
	  
        case 31: 
          // NCX40_CX31
          z1 = poslocal.z - Delta;
          weight_NCX3_CX4( poslocal, Delta, &x1, &y1 );
	  break; 
	  
        case 41: 
          // NCX40_CX41
          z1 = poslocal.z - Delta;
          weight_NCX4_CX4( poslocal, Delta, &x1, &y1 );
	  break; 	  	  	   	   	  	   	  
      }
      break;

    case 41:
      switch( CX )
      {	  
        case 10: 
          // NCX41_CX10
          z1 = poslocal.z - Delta;
          weight_NCX1_CX4( poslocal, Delta, &x1, &y1 );
	  break;
	  
        case 20: 
          // NCX41_CX20
          z1 = poslocal.z - Delta;
          weight_centred( poslocal, Delta, &x1, &y1 );
	  break;
	  
        case 30: 
          // NCX41_CX30
          z1 = poslocal.z - Delta;
          weight_NCX3_CX4( poslocal, Delta, &x1, &y1 );
	  break;
	  
        case 40: 
          // NCX41_CX40
          z1 = poslocal.z - Delta;
          weight_NCX4_CX4( poslocal, Delta, &x1, &y1 );
	  break; 
	  
        case 11: 
          // NCX41_CX11
          z1 = poslocal.z - 2. * Delta;
          weight_NCX1_CX4( poslocal, Delta, &x1, &y1 );
	  break;
	  
        case 21: 
          // NCX41_CX21
          z1 = poslocal.z - 2. * Delta;
          weight_centred( poslocal, Delta, &x1, &y1 );
	  break;
	  
        case 31: 
          // NCX41_CX31
          z1 = poslocal.z - 2. * Delta;
          weight_NCX3_CX4( poslocal, Delta, &x1, &y1 );
	  break; 
	  
        case 41: 
          // NCX41_CX41
          z1 = poslocal.z - 2. * Delta;
          weight_NCX4_CX4( poslocal, Delta, &x1, &y1 );
	  break;
      }
      break;
  }
  
  posref.x = GENERAL_1DIST( posb.x, x1 ) / ( 2. * Delta );
  posref.y = GENERAL_1DIST( posb.y, y1 ) / ( 2. * Delta );
  posref.z = GENERAL_1DIST( posb.z, z1 ) / ( 2. * Delta );
  
  return weight = Q2weighting( weight_id, posref.x, posref.y, posref.z ); 
#endif
}




//----------------------------------------------------------------------------
size_t is_in_cubic_boundingbox( const double x, const double y, 
	const double z, const GeomParameter gp ) 
//----------------------------------------------------------------------------
{
  size_t isin = 0;

  double xmax = gp.center.x + 0.5*gp.radius;
  double xmin = gp.center.x - 0.5*gp.radius;
  double ymax = gp.center.y + 0.5*gp.radius;
  double ymin = gp.center.y - 0.5*gp.radius;

#if dimension > 2  
  double zmax = gp.center.z + 0.5*gp.radius;
  double zmin = gp.center.z - 0.5*gp.radius;
#endif

  if ((x < xmax) && (x > xmin) && (y < ymax) && (y > ymin)
#if dimension > 2
      && (z< zmax) && (z > zmin)
#endif
      )
    isin = 1;
      
  return isin;
}




/** Finds the normal to the fictitious domain and 
assign the dial accordingly */
//----------------------------------------------------------------------------
void assign_dial_fd_boundary( RigidBody* p, const coord posb, 
	const GeomParameter gp, const double Delta, int* NCX ) 
//----------------------------------------------------------------------------
{
  bool isin_xp = false, isin_yp = false, isin_zp = false;  
  double RDelta = Delta / 100.;

  switch ( p->shape )
  {
    case SPHERE:
      isin_xp = is_in_Sphere( posb.x + RDelta, posb.y, posb.z, gp );
      isin_yp = is_in_Sphere( posb.x, posb.y + RDelta, posb.z, gp );
      isin_zp = is_in_Sphere( posb.x, posb.y, posb.z + RDelta, gp );
      break;
	  
    case CIRCULARCYLINDER2D:
      isin_xp = is_in_CircularCylinder2D( posb.x + RDelta, posb.y, gp );
      isin_yp = is_in_CircularCylinder2D( posb.x, posb.y + RDelta, gp );
      break;
	  
    case CUBE:
      isin_xp = is_in_Polyhedron( posb.x + RDelta, posb.y, posb.z, gp );
      isin_yp = is_in_Polyhedron( posb.x, posb.y + RDelta, posb.z, gp );
      isin_zp = is_in_Polyhedron( posb.x, posb.y, posb.z + RDelta, gp );
      break;

    case TETRAHEDRON:
      isin_xp = is_in_Polyhedron( posb.x + RDelta, posb.y, posb.z, gp );
      isin_yp = is_in_Polyhedron( posb.x, posb.y + RDelta, posb.z, gp );
      isin_zp = is_in_Polyhedron( posb.x, posb.y, posb.z + RDelta, gp );
      break;
      
    case OCTAHEDRON:
      isin_xp = is_in_Polyhedron( posb.x + RDelta, posb.y, posb.z, gp );
      isin_yp = is_in_Polyhedron( posb.x, posb.y + RDelta, posb.z, gp );
      isin_zp = is_in_Polyhedron( posb.x, posb.y, posb.z + RDelta, gp );
      break;

    case DODECAHEDRON:
      isin_xp = is_in_Polyhedron( posb.x + RDelta, posb.y, posb.z, gp );
      isin_yp = is_in_Polyhedron( posb.x, posb.y + RDelta, posb.z, gp );
      isin_zp = is_in_Polyhedron( posb.x, posb.y, posb.z + RDelta, gp );
      break;  
            
    case ICOSAHEDRON:
      isin_xp = is_in_Polyhedron( posb.x + RDelta, posb.y, posb.z, gp );
      isin_yp = is_in_Polyhedron( posb.x, posb.y + RDelta, posb.z, gp );
      isin_zp = is_in_Polyhedron( posb.x, posb.y, posb.z + RDelta, gp );
      break;         
	  
    default:
      fprintf( stderr,"Unknown Rigid Body shape !!\n" );
  }
  
# if dimension == 2
    if ( isin_xp ) 
    {
      if ( isin_yp )
        // fictitious-boundary's normal is oriented x- and y-. 
        *NCX = 3;    
      else
	// fictitious-boundary's normal is oriented x- and y+.
	*NCX = 2;
    }
    else
    {
      if ( isin_yp )
        // fictitious-boundary's normal is oriented x+ and y-. 
        *NCX = 4;
      else
        // fictitious-boundary's normal is oriented x+ and y+.
        *NCX = 1; 
    }
# else
    if ( isin_zp ) 
    {
      // fictitious-boundary's normal is oriented z-.  
      if ( isin_xp  ) 
      {
        if ( isin_yp )
	  // fictitious-boundary's normal is oriented x-, y-, z- . 
	  *NCX = 31;
        else
	  // fictitious-boundary's normal is oriented x-, y+, z-.
	  *NCX = 21;
      }
      else
      {
        if ( isin_yp )
	  // fictitious-boundary's normal is oriented x+, y-, z-. 
	  *NCX = 41;
        else
	  // fictitious-boundary's normal is oriented x+, y+, z-.
	  *NCX = 11;
      }
    }
    else 
    {
      // fictitious-boundary's normal is oriented z+.  
      if ( isin_xp ) 
      {
        if ( isin_yp )
	  // fictitious-boundary's normal is oriented x-, y-, z+ . 
	  *NCX = 30;
        else
	  // fictitious-boundary's normal is oriented x-, y+, z+.
	  *NCX = 20;
      }
      else 
      {
        if ( isin_yp )
	  // fictitious-boundary's normal is oriented x+, y-, z+. 
	  *NCX = 40;
        else
	  // fictitious-boundary's normal is oriented x+, y+, z+.
	  *NCX = 10;
      }
    }
  
  if ( *NCX == 0 )
    fprintf( stderr, "NCX = 0, isin_xp = %u, isin_yp = %u, isin_zp = %u \n",
	isin_xp, isin_yp, isin_zp );

#endif
}
