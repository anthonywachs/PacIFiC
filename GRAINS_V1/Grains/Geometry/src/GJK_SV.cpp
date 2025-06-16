#include "GJK_SV.hh"
#include "Matrix.hh"
#include "GrainsExec.hh"


/* ========================================================================== */
/*                             Low-Level Methods                              */
/* ========================================================================== */
static inline unsigned int compareSigns( double a, double b )
{
	// Maybe there's a faster way to deal with this set of operations?
	return static_cast<unsigned int>( !( ( a > 0 ) ^ ( b > 0 ) ) );
}
    



// -----------------------------------------------------------------------------
static inline void s1d( Vector3 const y[4], 
                        unsigned int& bits, 
                        double (&lambdas)[4] )
{
	// Identify the appropriate indices
	bool s1_set = false;
	unsigned int i1 = 0xffffffff, i2 = 0xffffffff;
	for ( unsigned int i = 0; i < 4; ++i )
	{
		if ( bits & ( 1 << i ) )
		{
			if ( s1_set )
			{
				i2 = i;
				break;
			}
			else
			{
				i1 = i;
				s1_set = true;
			}
		}
	}

	// Calculate the signed volume of the simplex.
	Vector3 t = y[i2] - y[i1];
	unsigned int I = 0;
	double neg_tI = -t[0];

	if ( fabs( t[1] ) > fabs( neg_tI ) )
	{
		I = 1;
		neg_tI = -t[1];
	}

	if ( fabs( t[2] ) > fabs( neg_tI ) )
	{
		I = 2;
		neg_tI = -t[2];
	}

	double pI = ( y[i2] * t ) / Norm2( t ) * neg_tI + y[i2][I];

	// Identify the signed volume resulting from replacing each point by the 
	// origin.
	double C[2] = { -y[i2][I] + pI, 
					 y[i1][I] - pI };
	unsigned int sign_comparisons[2] = { compareSigns( neg_tI, C[0] ), 
                                         compareSigns( neg_tI, C[1] ) };

	// If all signed volumes are identical, the origin lies inside the simplex.
	if ( sign_comparisons[0] + sign_comparisons[1] == 2 )
	{
		lambdas[i1] = C[0] / neg_tI;
		lambdas[i2] = C[1] / neg_tI;
	}
	else
	{
		// The point to retain is the one whose sign matches. In the
		// first case, the origin lies past the first point.
		if ( sign_comparisons[0] )
		{
			bits &= ~(1 << i2);
			lambdas[i1] = 1.;
		}
		else
		{
			bits &= ~(1 << i1);
			lambdas[i2] = 1.;
		}
	}
}




// -----------------------------------------------------------------------------
static inline void s2d( Vector3 const y[4], 
                        unsigned int& bits, 
                        double (&lambdas)[4] )
{
	unsigned int counter = 0, point0_idx = 0, point1_idx = 0, point2_idx = 0;
	for ( unsigned int i = 0; i < 4; ++i )
	{
		if ( bits & (1 << i) )
		{
			if ( counter == 0 )
				point0_idx = i;
			else if ( counter == 1)
				point1_idx = i;
			else
				point2_idx = i;
			counter += 1;
		}
	}

	Vector3 n = static_cast<Vector3>( y[point1_idx] - y[point0_idx] ) ^ 
				static_cast<Vector3>( y[point2_idx] - y[point0_idx] );
	Vector3 p0 = ( y[point0_idx] * n / Norm2( n ) ) * n;

	// Choose maximum area plane to project onto.
	// Make sure to store the *signed* area of the plane.
	// This loop is unrolled to save a few extra ops (assigning
	// an initial area of zero, an extra abs, etc)
	unsigned int idx_x = 1;
	unsigned int idx_y = 2;
	double mu_max = ( y[point1_idx][1] * y[point2_idx][2] + 
					  y[point0_idx][1] * y[point1_idx][2] +
					  y[point2_idx][1] * y[point0_idx][2] - 
					  y[point1_idx][1] * y[point0_idx][2] -
					  y[point2_idx][1] * y[point1_idx][2] - 
					  y[point0_idx][1] * y[point2_idx][2] );

	// This term is multiplied by -1.
	double mu = ( y[point1_idx][2] * y[point0_idx][0] + 
				  y[point2_idx][2] * y[point1_idx][0] +
				  y[point0_idx][2] * y[point2_idx][0] - 
				  y[point1_idx][2] * y[point2_idx][0] -
				  y[point0_idx][2] * y[point1_idx][0] - 
				  y[point2_idx][2] * y[point0_idx][0] );
	if ( fabs( mu ) > fabs( mu_max ) )
	{
		mu_max = mu;
		idx_x = 0;
	}

	mu = ( y[point1_idx][0] * y[point2_idx][1] + 
           y[point0_idx][0] * y[point1_idx][1] +
           y[point2_idx][0] * y[point0_idx][1] - 
           y[point1_idx][0] * y[point0_idx][1] - 
           y[point2_idx][0] * y[point1_idx][1] - 
           y[point0_idx][0] * y[point2_idx][1] );
	if ( fabs( mu ) > fabs( mu_max ) )
	{
		mu_max = mu;
		idx_x = 0;
		idx_y = 1;
	}

	// Compute the signed areas of each of the simplices formed by replacing an
	// index with a projection of the origin onto the area in this plane
	double C[3] = { 0. };
	bool sign_comparisons[3] = { false };

	C[0] = ( p0[idx_x] * y[point1_idx][idx_y] + 
             p0[idx_y] * y[point2_idx][idx_x] + 
             y[point1_idx][idx_x] * y[point2_idx][idx_y] - 
             p0[idx_x] * y[point2_idx][idx_y] - 
             p0[idx_y] * y[point1_idx][idx_x] - 
             y[point2_idx][idx_x] * y[point1_idx][idx_y] );
	sign_comparisons[0] = compareSigns( mu_max, C[0] );

	C[1] = ( p0[idx_x] * y[point2_idx][idx_y] + 
             p0[idx_y] * y[point0_idx][idx_x] + 
             y[point2_idx][idx_x] * y[point0_idx][idx_y] - 
             p0[idx_x] * y[point0_idx][idx_y] -
             p0[idx_y] * y[point2_idx][idx_x] - 
             y[point0_idx][idx_x] * y[point2_idx][idx_y] );
	sign_comparisons[1] = compareSigns( mu_max, C[1] );

	C[2] = ( p0[idx_x] * y[point0_idx][idx_y] + 
             p0[idx_y] * y[point1_idx][idx_x] +
             y[point0_idx][idx_x] * y[point1_idx][idx_y] - 
             p0[idx_x] * y[point1_idx][idx_y] -
             p0[idx_y] * y[point0_idx][idx_x] - 
             y[point1_idx][idx_x] * y[point0_idx][idx_y] );
	sign_comparisons[2] = compareSigns( mu_max, C[2] );

	if ( sign_comparisons[0] + sign_comparisons[1] + sign_comparisons[2] == 3 )
	{
		lambdas[point0_idx] = C[0] / mu_max;
		lambdas[point1_idx] = C[1] / mu_max;
		lambdas[point2_idx] = C[2] / mu_max;
	}
	else
	{
		double d = 1.e9;
		Vector3 new_point;
		unsigned int new_bits = 0;
		for ( unsigned int j = 0; j < 3; ++j )
		{
			if ( !sign_comparisons[j] )
			{
				unsigned int new_used = bits;
				// Test removal of the current point.
				if ( j == 0 )
					new_used &= ~(1 << point0_idx);
				else if ( j == 1 )
					new_used &= ~(1 << point1_idx);
				else
					new_used &= ~(1 << point2_idx);

				double new_lambdas[4] = { 0. };

				s1d( y, new_used, new_lambdas );
				// Consider resetting in place if possible.
				new_point[0] = 0;
				new_point[1] = 0;
				new_point[2] = 0;
				for ( unsigned int i = 0; i < 4; ++i )
				{
					if ( new_used & (1 << i) )
						new_point += new_lambdas[i] * y[i];
				}
				double d_star = new_point * new_point;
				if ( d_star < d )
				{
					new_bits = new_used;
					d = d_star;
					for ( unsigned int i = 0; i < 4; ++i )
						lambdas[i] = new_lambdas[i];
				}
			}
		}
		bits = new_bits;
	}
}




// -----------------------------------------------------------------------------
static inline void s3d( Vector3 const y[4], 
                        unsigned int& bits, 
                        double (&lambdas)[4] )
{
	double C[4] = { 0. };

	// Compute all minors and the total determinant of the matrix M,
	// which is the transpose of the y matrix with an extra row of
	// ones at the bottom. Since the indexing is nontrivial and the
	// array is small (and we can save on some negation), all the
	// computations are done directly rather than with a loop.
	// C[0] and C[2] are negated due to the (-1)^(i+j+1) prefactor,
	// where i is always 4 because we're expanding about the 4th row.
	C[0] =  y[3][0] * y[2][1] * y[1][2] + 
			y[2][0] * y[1][1] * y[3][2] + 
			y[1][0] * y[3][1] * y[2][2] -
			y[1][0] * y[2][1] * y[3][2] - 
			y[2][0] * y[3][1] * y[1][2] -
			y[3][0] * y[1][1] * y[2][2];
	C[1] =  y[0][0] * y[2][1] * y[3][2] + 
			y[2][0] * y[3][1] * y[0][2] + 
			y[3][0] * y[0][1] * y[2][2] -
			y[3][0] * y[2][1] * y[0][2] - 
			y[2][0] * y[0][1] * y[3][2] -
			y[0][0] * y[3][1] * y[2][2];
	C[2] =  y[3][0] * y[1][1] * y[0][2] + 
			y[1][0] * y[0][1] * y[3][2] + 
			y[0][0] * y[3][1] * y[1][2] -
			y[0][0] * y[1][1] * y[3][2] - 
			y[1][0] * y[3][1] * y[0][2] -
			y[3][0] * y[0][1] * y[1][2];
	C[3] =  y[0][0] * y[1][1] * y[2][2] + 
			y[1][0] * y[2][1] * y[0][2] + 
			y[2][0] * y[0][1] * y[1][2] -
			y[2][0] * y[1][1] * y[0][2] - 
			y[1][0] * y[0][1] * y[2][2] -
			y[0][0] * y[2][1] * y[1][2];
	double dM = C[0] + C[1] + C[2] + C[3];

	unsigned int sign_comparisons[4] = {0};
	sign_comparisons[0] = compareSigns( dM, C[0] );
	sign_comparisons[1] = compareSigns( dM, C[1] );
	sign_comparisons[2] = compareSigns( dM, C[2] );
	sign_comparisons[3] = compareSigns( dM, C[3] );

	if ( ( sign_comparisons[0] + sign_comparisons[1] + 
		   sign_comparisons[2] + sign_comparisons[3] ) == 4)
	{
		for ( unsigned int i = 0; i < 4; ++i )
			lambdas[i] = C[i] / dM;
	}
	else
	{
		double d = 1.e9, d_star = 0.;
		Vector3 new_point;
		unsigned int new_bits = 0;
		for ( unsigned int j = 0; j < 4; ++j )
		{
			if ( !sign_comparisons[j] )
			{
				// Test removal of the current point.
				unsigned int new_used = bits;
				new_used &= ~(1 << j);
				double new_lambdas[4] = {0.};

				s2d( y, new_used, new_lambdas );
	
				new_point = Vector3();
				for ( unsigned int i = 0; i < 4; ++i )
				{
					if ( new_used & (1 << i) )
						new_point += new_lambdas[i] * y[i];
				}
				d_star = new_point * new_point;
				if ( d_star < d )
				{
					new_bits = new_used;
					d = d_star;
					for ( unsigned int i = 0; i < 4; ++i )
						lambdas[i] = new_lambdas[i];
				}
			}
		}
		bits = new_bits;
	}
}




// -----------------------------------------------------------------------------
static inline void computeVector( unsigned int const bits,
                                  Vector3 const y[4],
								  double const lambdas[4],
                                  Vector3& v )
{
    v.setValue( 0., 0., 0. );
	for ( unsigned int i = 0; i < 4; ++i )
	{
		if ( bits & ( 1 << i ) )
			v += lambdas[i] * y[i];
	}
}




// -----------------------------------------------------------------------------
static inline void computePoints( unsigned int const bits,
                                  Point3 const p[4],
                                  Point3 const q[4],
                                  double const lambdas[4],
                                  Point3& p1,
                                  Point3& p2 )
{
    p1.setValue( 0., 0., 0. );
    p2.setValue( 0., 0., 0. );
    for ( unsigned int i = 0; i < 4; ++i )
    {
        if ( bits & ( 1 << i ) )
        {
            p1 += lambdas[i] * p[i];
            p2 += lambdas[i] * q[i];
        } 
    }
}




// -----------------------------------------------------------------------------
static inline void sv_subalgorithm( Vector3 const y[4], 
                                    unsigned int& bits, 
                                    double (&lambdas)[4],
									Vector3& v )
{
	// The y array is never modified by this function.  The bits may be
	// modified if necessary, and the lambdas will be updated.  All the other
	// functions (if they need to make deeper calls e.g. s3d->s2d) will have to
	// make copies of bits to avoid overwriting that data incorrectly.
	unsigned int num_used = 0;
	for ( unsigned int i = 0; i < 4; ++i )
		num_used += (bits >> i) & 1;

	// Start with the most common cases.
	if ( num_used == 1 )
	{
		for ( unsigned int i = 0; i < 4; ++i )
		{
			if ( bits & (1 << i) )
				lambdas[i] = 1.;
		}
	}
	else if ( num_used == 2 )
		s1d( y, bits, lambdas );
	else if ( num_used == 3 )
		s2d( y, bits, lambdas );
	else
		s3d( y, bits, lambdas );

	computeVector( bits, y, lambdas, v );
}
    
    


// -----------------------------------------------------------------------------
// The next function is used for detecting degenerate cases that cause
// termination problems due to rounding errors.
static inline bool degenerate( unsigned int const bits,
                               Vector3 const y[4],
                               Vector3 const& w )
{
	for ( unsigned int i = 0, bit = 1; i < 4; ++i, bit <<= 1 )
	{
		if ( ( bits & bit ) && y[i] == w )
			return ( true );
	}
	return ( false );
}




/* ========================================================================== */
/*                            High-Level Methods                              */
/* ========================================================================== */
double closest_points_GJK_SV( Convex const& a, 
                              Convex const& b, 
                              Transform const& a2w,
                              Transform const& b2w,
                              Point3& pa,
                              Point3& pb,
                              int& nbIter )
{
    // GJK variables
    unsigned int bits = 0;           // identifies current simplex
    unsigned int last = 0;           // identifies last found support point
    Point3 p[4];                     // support points of A in local
    Point3 q[4];                     // support points of B in local
    Vector3 y[4];				     // support points of A-B in world
    double mu = 0.;                  // optimality gap
    int numIterations = 0;           // No. iterations
    double lambdas[4] = { 0. };      // Weights

    // Misc variables, e.g. tolerance, ...
    double relError = GrainsExec::m_colDetTolerance;    // rel error for opt gap
    double absError = 1.e-4 * relError;          // abs error for optimality gap
    bool acceleration = GrainsExec::m_colDetAcceleration;     // isAcceleration?
    double momentum = 0., oneMinusMomentum = 1.;  // in case we use acceleration

    // Initializing vectors
    Vector3 v( a2w( a.support( Vector3Null ) ) - 
               b2w( b.support( Vector3Null ) ) );
    Vector3 w( v );
    Vector3 d( v );
    double dist = Norm( v );

    while ( bits < 15 && dist > EPSILON2 && numIterations < 1000 )
    {
        ++numIterations;
        // Updating the bits, ...
        for ( unsigned int new_index = 0; new_index < 4; ++new_index )
        {
            // At least one of these must be empty, otherwise we have an overlap.
            if ( !( bits & ( 1 << new_index ) ) )
			{
				last = new_index;
				break;
			}
        }

        // Finding the suitable direction using either Nesterov or original
        // The number 8 is hard-coded. Emprically, it shows the best convergence
        // for superquadrics. For the rest of shapes, we really do not need to 
        // use Nesterov as the improvemenet is marginal.
        if ( acceleration && numIterations % 8 != 0 )
        {
            momentum = numIterations / ( numIterations + 2. );
            oneMinusMomentum = 1. - momentum;
            d = momentum * d + 
                momentum * oneMinusMomentum * v +
                oneMinusMomentum * oneMinusMomentum * w;
        }
        else
            d = v;

        p[last] = a.support( ( -d ) * a2w.getBasis() );
        q[last] = b.support( (  d ) * b2w.getBasis() );
        w = a2w( p[last] ) - b2w( q[last] );
        
        // termination criteria -- optimiality gap
        mu = dist - v * w / dist;
        if ( mu < dist * relError || mu < absError )
        {
            if ( acceleration )
            {
                if ( Norm( d - v ) < EPSILON )
                    break;
                acceleration = false;
                p[last] = a.support( ( -v ) * a2w.getBasis() );
                q[last] = b.support( (  v ) * b2w.getBasis() );
                w = a2w( p[last] ) - b2w( q[last] );
            }
            else
                break;
        }
        // termination criteria -- degenerate case
        if ( degenerate( bits, y, w ) )
            break;
        
        // if not terminated, get ready for the next iteration
        y[last] = w;
        bits |= ( 1 << last );
		sv_subalgorithm( y, bits, lambdas, v );
        dist = Norm( v );
    }
    // compute witness points
    computePoints( bits, p, q, lambdas, pa, pb );
    nbIter = numIterations;
    return ( dist );
}






// -----------------------------------------------------------------------------
double closest_points_GJK_SV( Convex const& a, 
                              Convex const& b, 
                              Transform const& a2w,
                              Transform const& b2w,
                              Vector3& v,
                              Point3& pa,
                              Point3& pb,
                              int& nbIter )
{
  	// GJK variables
    unsigned int bits = 0;           // identifies current simplex
    unsigned int last = 0;           // identifies last found support point
    Point3 p[4];                     // support points of A in local
    Point3 q[4];                     // support points of B in local
    Vector3 y[4];				     // support points of A-B in world
    double mu = 0.;                  // optimality gap
    int numIterations = 0;           // No. iterations
    double lambdas[4];               // Weights

    // Misc variables, e.g. tolerance, ...
    double relError = GrainsExec::m_colDetTolerance;    // rel error for opt gap
    double absError = 1.e-4 * relError;          // abs error for optimality gap
    bool acceleration = GrainsExec::m_colDetAcceleration;     // isAcceleration?
    double momentum = 0., oneMinusMomentum = 1.;  // in case we use acceleration

    // Initializing vectors
    Vector3 w;
    if ( v != Vector3Null )
        w = a2w( a.support( ( -v ) * a2w.getBasis() ) ) - 
            b2w( b.support( (  v ) * b2w.getBasis() ) );
    else // if we don't have a better guess
    {
        v = a2w( a.support( Vector3Null ) ) - b2w( b.support( Vector3Null ) );
        w = v;
    }
    Vector3 d( v );
    double dist = 1.;

    while ( bits < 15 && dist > EPSILON2 && numIterations < 1000 )
    {
        ++numIterations;
        // Updating the bits, ...
        for ( unsigned int new_index = 0; new_index < 4; ++new_index )
        {
            // At least one of these must be empty, otherwise we have an overlap.
            if ( !( bits & ( 1 << new_index ) ) )
			{
				last = new_index;
				break;
			}
        }

        // Finding the suitable direction using either Nesterov or original
        // The number 8 is hard-coded. Emprically, it shows the best convergence
        // for superquadrics. For the rest of shapes, we really do not need to 
        // use Nesterov as the improvemenet is marginal.
        if ( acceleration && numIterations % 8 != 0 )
        {
            momentum = numIterations / ( numIterations + 2. );
            oneMinusMomentum = 1. - momentum;
            d = momentum * d + 
                momentum * oneMinusMomentum * v +
                oneMinusMomentum * oneMinusMomentum * w;
        }
        else
            d = v;

        p[last] = a.support( ( -d ) * a2w.getBasis() );
        q[last] = b.support( (  d ) * b2w.getBasis() );
        w = a2w( p[last] ) - b2w( q[last] );
        
        // termination criteria -- optimiality gap
        mu = dist - v * w / dist;
        if ( mu < dist * relError || mu < absError )
        {
            if ( acceleration )
            {
                if ( Norm( d - v ) < EPSILON )
                    break;
                acceleration = false;
                p[last] = a.support( ( -v ) * a2w.getBasis() );
                q[last] = b.support( (  v ) * b2w.getBasis() );
                w = a2w( p[last] ) - b2w( q[last] );
            }
            else
                break;
        }
        // termination criteria -- degenerate case
        if ( degenerate( bits, y, w ) )
            break;
        
        // if not terminated, get ready for the next iteration
        y[last] = w;
        bits |= ( 1 << last );
		sv_subalgorithm( y, bits, lambdas, v );
        dist = Norm( v );
    }
    // compute witness points
    computePoints( bits, p, q, lambdas, pa, pb );
    nbIter = numIterations;
    return ( dist );
}