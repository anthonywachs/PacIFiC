#include <geomVector.hh>
//#include <GE_Point.hh>
#include <MAC_System.hh>
#include <math.h>


//----------------------------------------------------------------------
geomVector:: geomVector( )
//----------------------------------------------------------------------
{
   MAC_LABEL( "geomVector:: geomVector()" ) ;
       
   vecSize = 0;
   xx = NULL;
   
} 




//----------------------------------------------------------------------
geomVector:: geomVector( const size_t &vecSize_ )
//----------------------------------------------------------------------
{
   MAC_LABEL( "geomVector:: geomVector(const size_t &vecSize_)" ) ;

   MAC_CHECK( vecSize_ >= 1) ;
   MAC_CHECK( vecSize_ <= 3) ;
   
   vecSize = vecSize_;
   xx = new double[vecSize]; 
   for (size_t i = 0; i < vecSize; ++i) xx[i] = 0. ;
	
} 




//----------------------------------------------------------------------
geomVector:: geomVector( const double &x, const double &y )
//----------------------------------------------------------------------
{
   MAC_LABEL( "geomVector:: geomVector(x,y)" ) ;
   
   vecSize = 2;
   xx = new double[vecSize]; 
   xx[0] = x;
   xx[1] = y;
	
} 




//----------------------------------------------------------------------
geomVector:: geomVector( const double &x, const double &y, const double &z )
//----------------------------------------------------------------------
{
   MAC_LABEL( "geomVector:: geomVector(x,y,z)" ) ;
   
   vecSize = 3;
   xx = new double[vecSize]; 
   xx[0] = x;
   xx[1] = y;
   xx[2] = z;
	
} 




//----------------------------------------------------------------------
void geomVector::resize( const size_t &vecSize_ )
//----------------------------------------------------------------------
{
   MAC_LABEL( "geomVector:: resize" ) ;
    
   MAC_CHECK( vecSize_ >= 1) ;
   MAC_CHECK( vecSize_ <= 3) ;
   
   vecSize = vecSize_;
   if (xx) delete [] xx;
   xx = new double[vecSize]; 
   for (size_t i = 0; i < vecSize; ++i) xx[i] = 0.;
	
} 




//----------------------------------------------------------------------
void geomVector::set( const double &val )
//----------------------------------------------------------------------
{
   MAC_LABEL( "geomVector:: set(x,y)" ) ;
   MAC_CHECK_INV( check_invariant() );    
   
   for (size_t i = 0; i < vecSize; ++i) xx[i] = val;   
	
} 




//----------------------------------------------------------------------
void geomVector::set( const double &x, const double &y )
//----------------------------------------------------------------------
{
   MAC_LABEL( "geomVector:: set(x,y)" ) ;
   MAC_CHECK_INV( check_invariant() );    
   MAC_CHECK( vecSize == 2 ) ;
   
   xx[0] = x;
   xx[1] = y;
	
} 




//----------------------------------------------------------------------
void geomVector::set( const double &x, const double &y, const double &z )
//----------------------------------------------------------------------
{
   MAC_LABEL( "geomVector:: set(x,y,z)" ) ;
   MAC_CHECK_INV( check_invariant() );    
   MAC_CHECK( vecSize == 3 ) ;
   
   xx[0] = x;
   xx[1] = y;
   xx[2] = z;
	
} 




// //----------------------------------------------------------------------
// void geomVector::set( const GE_Point *gep )
// //----------------------------------------------------------------------
// {
//    MAC_LABEL( "geomVector:: set(gep)" ) ;
//       
//    resize( gep->nb_coordinates() );
//    for (size_t i=0;i<vecSize;++i) xx[i] = gep->coordinate(i);
// 	
// } 




//----------------------------------------------------------------------
geomVector:: geomVector( const geomVector &other )
//----------------------------------------------------------------------
{
   MAC_LABEL( "geomVector:: geomVector(other)" ) ;

   MAC_CHECK_INV( other.check_invariant() );

   vecSize = other.vecSize;  
   xx = new double[vecSize]; 
   for (size_t i = 0; i < vecSize; ++i) xx[i] = other.xx[i];
	   
} 




//----------------------------------------------------------------------
geomVector::~geomVector()
//----------------------------------------------------------------------
{
   MAC_LABEL( "geomVector::~geomVector()" ) ;

   if (xx) delete [] xx;
      
} 




//----------------------------------------------------------------------
bool 
geomVector:: isZeroVal() const
//----------------------------------------------------------------------
{
   MAC_LABEL( "geomVector:: isZero" ) ;
   MAC_CHECK_INV( check_invariant() );

   for (size_t i = 0; i < vecSize; ++i)
     if (xx[i] != 0.) return false;

   return true;
   
} 




//----------------------------------------------------------------------
double
geomVector:: calcNorm() const
//----------------------------------------------------------------------
{
   MAC_LABEL( "geomVector:: calcNorm" ) ;
   MAC_CHECK_INV( check_invariant() );

   double norm = 0.; 
   for (size_t i = 0; i < vecSize; ++i) norm += xx[i] * xx[i];

   return sqrt(norm);
   
} 




//----------------------------------------------------------------------
double
geomVector:: calcNormSquare() const
//----------------------------------------------------------------------
{
   MAC_LABEL( "geomVector:: calcNormSquare" ) ;
   MAC_CHECK_INV( check_invariant() );

   double norm = 0.; 
   for (size_t i = 0; i < vecSize; ++i) norm += xx[i] * xx[i];

   return norm;
   
} 




//----------------------------------------------------------------------
geomVector& 
geomVector::operator=( geomVector const& other )
//----------------------------------------------------------------------
{
   MAC_LABEL( "geomVector:: operator=" ) ;
   MAC_CHECK_INV( other.check_invariant()) ;

   if ( &other != this )
   {
     vecSize = other.vecSize;
     if (xx) delete [] xx;
     xx = new double[vecSize]; 
     for (size_t i = 0; i < vecSize; ++i) xx[i] = other.xx[i];	
   }

   return (*this);

}




//----------------------------------------------------------------------
bool 
geomVector::operator==( geomVector const &other ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "geomVector:: operator==" ) ;
   MAC_CHECK_INV( check_invariant() );
   MAC_CHECK_INV( other.check_invariant()) ;
   MAC_CHECK( vecSize == other.vecSize );

   for (size_t i = 0; i < vecSize; ++i)
     if (xx[i] != other.xx[i]) return false;   
	  
   return( true ) ;
   
}




//----------------------------------------------------------------------
std::ostream& operator <<( std::ostream& f, geomVector const &P )
//----------------------------------------------------------------------
{
   MAC_LABEL( "geomVector:: operator <<" ) ;
   MAC_CHECK_INV( P.check_invariant() );
   	
   f.setf(std::ios::scientific,std::ios::floatfield);

   for (size_t i = 0; i < P.vecSize; ++i) {     
     if (P.xx[i] < 0.) f.precision(5);
     else f.precision(6);
     f << P.xx[i]; 
     if ( i < P.vecSize - 1 ) f << " ";
   }    

   f.setf(std::ios::fixed);
   f.precision(6);
	
   return f;
   
}




//----------------------------------------------------------------------
std::istream& operator >>( std::istream& f, geomVector &P )
//----------------------------------------------------------------------
{
   MAC_LABEL( "geomVector:: operator >>" ) ;
   MAC_CHECK_INV( P.check_invariant() );
		
   for (size_t i = 0; i < P.vecSize; ++i)
     f >> P.xx[i];
	
   return f;
   
}




//----------------------------------------------------------------------
geomVector 
geomVector::operator*( const double &c ) const 
//----------------------------------------------------------------------
{
   MAC_LABEL( "geomVector:: operator *double" ) ;
   MAC_CHECK_INV( check_invariant() );
	
   geomVector tmpVec(vecSize);
   for (size_t i = 0; i < vecSize; ++i) 
     tmpVec.xx[i] = c * xx[i];

   return tmpVec;
   
} 




//----------------------------------------------------------------------
geomVector operator*( const double &c, const geomVector &secondVec )
//----------------------------------------------------------------------
{
   MAC_LABEL( "geomVector:: operator *" ) ;
   MAC_CHECK_INV( secondVec.check_invariant() );
	
   return secondVec*c;
   
} 




//----------------------------------------------------------------------
geomVector 
geomVector::operator+( const geomVector &secondVec ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "geomVector:: operator +" ) ;
   MAC_CHECK_INV( check_invariant() );
   MAC_CHECK_INV( secondVec.check_invariant() );
   MAC_CHECK( vecSize == secondVec.vecSize );
	
   geomVector tmpVec(vecSize);
   for (size_t i = 0; i < vecSize; ++i) 
     tmpVec.xx[i] = xx[i] + secondVec.xx[i];

   return tmpVec;
   
} 




//----------------------------------------------------------------------
geomVector 
geomVector::operator-( const geomVector &secondVec ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "geomVector:: operator -" ) ;
   MAC_CHECK_INV( check_invariant() );
   MAC_CHECK_INV( secondVec.check_invariant() );
   MAC_CHECK( vecSize == secondVec.vecSize );
	
   geomVector tmpVec(vecSize);
   for (size_t i = 0; i < vecSize; ++i) 
     tmpVec.xx[i] = xx[i] - secondVec.xx[i];

   return tmpVec;
   
} 




//----------------------------------------------------------------------
geomVector 
geomVector::operator-=( const geomVector &secondVec )
//----------------------------------------------------------------------
{
   MAC_LABEL( "geomVector:: operator -=" ) ;
   MAC_CHECK_INV( check_invariant() );
   MAC_CHECK_INV( secondVec.check_invariant() );
   MAC_CHECK( vecSize == secondVec.vecSize );
	
   for (size_t i = 0; i < vecSize; ++i) 
     xx[i] -= secondVec.xx[i];

   return *this;
   
} 




//----------------------------------------------------------------------
geomVector 
geomVector::operator+=( const geomVector &secondVec )
//----------------------------------------------------------------------
{
   MAC_LABEL( "geomVector:: operator +=" ) ;
   MAC_CHECK_INV( check_invariant() );
   MAC_CHECK_INV( secondVec.check_invariant() );
   MAC_CHECK( vecSize == secondVec.vecSize );
	
   for (size_t i = 0; i < vecSize; ++i) 
     xx[i] += secondVec.xx[i];

   return *this;
   
} 




//----------------------------------------------------------------------
geomVector 
geomVector::operator*=( const double &c )
//----------------------------------------------------------------------
{
   MAC_LABEL( "geomVector:: operator *=" ) ;
   MAC_CHECK_INV( check_invariant() );
	
   for (size_t i = 0; i < vecSize; ++i) 
     xx[i] *= c;

   return *this;
   
} 




//----------------------------------------------------------------------
geomVector 
geomVector::operator/=( const double &c )
//----------------------------------------------------------------------
{
   MAC_LABEL( "geomVector:: operator /=" ) ;
   MAC_CHECK_INV( check_invariant() );
   MAC_CHECK( c != 0. );   
	
   for (size_t i = 0; i < vecSize; ++i) 
     xx[i] /= c;

   return *this;
   
} 




//----------------------------------------------------------------------
double 
geomVector::operator,( const geomVector &secondVec ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "geomVector:: operator ," ) ;
   MAC_CHECK_INV( check_invariant() );
   MAC_CHECK_INV( secondVec.check_invariant() );
   MAC_CHECK( vecSize == secondVec.vecSize );
	
   double scalar_product = 0.;
   for (size_t i = 0; i < vecSize; ++i) 
     scalar_product += xx[i] * secondVec.xx[i];

   return scalar_product;
   
} 




//----------------------------------------------------------------------
geomVector 
geomVector::operator^( const geomVector &secondVec ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "geomVector:: operator ^" ) ;
   MAC_CHECK_INV( check_invariant() );
   MAC_CHECK_INV( secondVec.check_invariant() );
   MAC_CHECK( vecSize == secondVec.vecSize );
	
   geomVector cross_product(3);
   if (vecSize >= 2 )
   {
     cross_product.xx[2] = xx[0] * secondVec.xx[1]
     	- xx[1] * secondVec.xx[0];
     if (vecSize == 3 )
     {
       cross_product.xx[0] = xx[1] * secondVec.xx[2]
       	- xx[2] * secondVec.xx[1];
       cross_product.xx[1] = xx[2] * secondVec.xx[0]
       	- xx[0] * secondVec.xx[2];     
     }
   }

   return cross_product;
   
}   




//----------------------------------------------------------------------
void 
geomVector::addOneComp( const size_t index, const double &value )
//----------------------------------------------------------------------
{
   MAC_LABEL( "geomVector:: addOneComp" ) ;
   MAC_CHECK_INV( check_invariant() );
	
   xx[index] += value;
   
}




//----------------------------------------------------------------------
void 
geomVector::setVec( const std::vector<double> &newVec )
//----------------------------------------------------------------------
{
   MAC_LABEL( "geomVector:: setVec" ) ;
   MAC_CHECK( newVec.size() >= 1) ;
   MAC_CHECK( newVec.size() <= 3) ;
	
   if (xx) delete [] xx;
   vecSize = newVec.size();
   xx = new double[vecSize]; 
   for (size_t i = 0; i < vecSize; ++i) xx[i] = newVec[i];
      	
}  




//----------------------------------------------------------------------
double 
geomVector::calcDist( const geomVector &secondVec ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "geomVector:: calcDist(const geomVector secondVec)" ) ;
   MAC_CHECK_INV( check_invariant() );
   MAC_CHECK_INV( secondVec.check_invariant() );
   MAC_CHECK( vecSize == secondVec.vecSize );
	
   double dist = 0.;
   for (size_t i = 0; i < vecSize; ++i) 
     dist += ( xx[i] - secondVec.xx[i] ) * ( xx[i] - secondVec.xx[i] );
	
   return sqrt(dist);
   
}     




//----------------------------------------------------------------------
double 
geomVector::calcDist( geomVector const* secondVec ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "geomVector:: calcDist(const geomVector secondVec)" ) ;
   MAC_CHECK_INV( check_invariant() );
   MAC_CHECK_INV( secondVec->check_invariant() );
   MAC_CHECK( vecSize == secondVec->vecSize );
	
   double dist = 0.;
   for (size_t i = 0; i < vecSize; ++i) 
     dist += ( xx[i] - secondVec->xx[i] ) * ( xx[i] - secondVec->xx[i] );
	
   return sqrt(dist);
   
}     




// //----------------------------------------------------------------------
// double 
// geomVector::calcDist( const GE_Point* secondPoint ) const
// //----------------------------------------------------------------------
// {
//    MAC_LABEL( "geomVector:: calcDist(const GE_Point* secondPoint)" ) ;
//    MAC_CHECK_INV( check_invariant() );
//    MAC_CHECK( vecSize == secondPoint->nb_coordinates() );
// 	
//    double dist = 0.;
//    for (size_t i = 0; i < vecSize; ++i) 
//      dist += ( xx[i] - secondPoint->coordinate(i) )
//      	*( xx[i] - secondPoint->coordinate(i) );
// 	
//    return sqrt(dist);
//    
// }     




//----------------------------------------------------------------------
double 
geomVector::calcDist( double x ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "geomVector:: calcDist(double x)" ) ;
   MAC_CHECK_INV( check_invariant() );
   MAC_CHECK( vecSize == 1 );
	
   double dist = fabs( xx[0] - x ) ;
	
   return dist;
   
}    




//----------------------------------------------------------------------
double 
geomVector::calcDist( double x, double y ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "geomVector:: calcDist(double x, double y)" ) ;
   MAC_CHECK_INV( check_invariant() );
   MAC_CHECK( vecSize == 2 );   
	
   double dist = pow( xx[0] - x, 2. ) + pow( xx[1] - y, 2. ) ;   
	
   return sqrt(dist);
   
}



    
//----------------------------------------------------------------------
double 
geomVector::calcDist( double x, double y, double z ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "geomVector:: calcDist(double x, double y, double z)" ) ;
   MAC_CHECK_INV( check_invariant() );
   MAC_CHECK( vecSize == 3 );	

   double dist = pow( xx[0] - x, 2. ) + pow( xx[1] - y, 2. ) 
   	+ pow( xx[2] - z, 2. );
	
   return sqrt(dist);
   
}    




//----------------------------------------------------------------------
bool 
geomVector::check_invariant( ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "geomVector:: check_invariant" ) ;
   MAC_ASSERT( vecSize >= 1) ;
   MAC_ASSERT( vecSize <= 3) ;
   MAC_ASSERT( xx != NULL ) ;

   return( true ) ;
   
}




//----------------------------------------------------------------------
void 
geomVector::setVecZero( )
//----------------------------------------------------------------------
{
   MAC_LABEL( "geomVector:: setVecZero" ) ;
   MAC_CHECK_INV( check_invariant() );
	
   for (size_t i = 0; i < vecSize; ++i) xx[i] = 0.;

}




//----------------------------------------------------------------------
void 
geomVector::translate( const geomVector &translation_vector )
//----------------------------------------------------------------------
{
   MAC_LABEL( "geomVector:: translate" ) ;
   MAC_CHECK_INV( check_invariant() );
   MAC_CHECK_INV( translation_vector.check_invariant() );
   MAC_CHECK( vecSize == translation_vector.vecSize );
   	
   for (size_t i = 0; i < vecSize; ++i) 
     xx[i] += translation_vector.xx[i];

}



// //----------------------------------------------------------------------
// void 
// geomVector::copy_coordinates( GE_Point *gep ) const
// //----------------------------------------------------------------------
// {
//    MAC_LABEL( "geomVector:: copy_coordinates" ) ;
//    MAC_CHECK_INV( check_invariant() );
//    MAC_CHECK( vecSize == gep->nb_coordinates() );
//    	
//    for (size_t i = 0; i < vecSize; ++i) 
//      gep->set_coordinate( i, xx[i] );   
// 
// }




//----------------------------------------------------------------------
void 
geomVector::retract( const geomVector &gravity_center,const double &amount )
//----------------------------------------------------------------------
{
   MAC_LABEL( "geomVector:: retract" ) ;
   MAC_CHECK_INV( check_invariant() );
   MAC_CHECK_INV( gravity_center.check_invariant() );
   MAC_CHECK( vecSize == gravity_center.vecSize );
   	
   geomVector translation_vector = *this - gravity_center;
   double norm_translation_vector = translation_vector.calcNorm();
   for (size_t i = 0; i < vecSize; ++i) 
      xx[i] =  gravity_center.xx[i] 
      	+ ( 1. - amount / norm_translation_vector ) * translation_vector.xx[i] ;

}
