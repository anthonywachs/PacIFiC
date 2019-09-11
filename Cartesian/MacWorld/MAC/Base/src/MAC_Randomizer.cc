#include <MAC_Randomizer.hh>

#include <MAC_assertions.hh>

#include <size_t_vector.hh>
#include <iostream>

//----------------------------------------------------------------------
static const int IM1=2147483563 ;
static const int IM2=2147483399 ;
static const double AM=1./IM1 ;
static const int IMM1=IM1-1 ;
static const int IA1=40014 ;
static const int IA2=40692 ;
static const int IQ1=53668 ;
static const int IQ2=52774 ;
static const int IR1=12211 ;
static const int IR2=3791 ;
static const int NTAB=32 ;
static const int NDIV=1+IMM1/NTAB ;
static const double EPS=1.2E-7 ;
static const double RNMX=1.-EPS ;
//----------------------------------------------------------------------


//----------------------------------------------------------------------
MAC_Randomizer*
MAC_Randomizer:: create_clone( MAC_Object* a_owner ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Randomizer:: create_clone" ) ;
   MAC_Randomizer* result =  new MAC_Randomizer( a_owner, my_series ) ;

   MAC_CHECK_POST( create_clone_POST( result,  a_owner ) ) ;
   return( result ) ;
}



//----------------------------------------------------------------------
MAC_Randomizer*
MAC_Randomizer:: create( MAC_Object* a_owner, int series )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Randomizer:: create" ) ;
   MAC_Randomizer* result =  new MAC_Randomizer( a_owner, series ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}



//----------------------------------------------------------------------
MAC_Randomizer:: MAC_Randomizer( MAC_Object* a_owner,
                               int series )
//----------------------------------------------------------------------
   : MAC_Object( a_owner ),
     my_series( series ),
     IV(NTAB)
{
   MAC_LABEL( "MAC_Randomizer:: MAC_Randomizer" ) ;
   start() ;
}



//----------------------------------------------------------------------
MAC_Randomizer:: ~MAC_Randomizer( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Randomizer:: ~MAC_Randomizer" ) ;
}



//----------------------------------------------------------------------
void
MAC_Randomizer:: start( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Randomizer:: start" ) ;
   IDUM = my_series ;
   IDUM2 = my_series ;
   IY = 0 ;
   
   for( size_t j=NTAB+8; j>=1; j-- )
   {
      int K=IDUM/IQ1 ;
      IDUM=IA1*(IDUM-K*IQ1)-K*IR1 ;
      if( IDUM<0 ) IDUM=IDUM+IM1 ;
      if( j<=NTAB ) IV(j-1)=IDUM ;
   }
   IY=IV(0) ;
   go_next() ;
}



//----------------------------------------------------------------------
void
MAC_Randomizer:: go_next( void ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Randomizer:: go_next" ) ;
   
   int K=IDUM/IQ1 ;
   IDUM=IA1*(IDUM-K*IQ1)-K*IR1 ;
   if (IDUM<0) IDUM=IDUM+IM1 ;
   K=IDUM2/IQ2 ;
   IDUM2=IA2*(IDUM2-K*IQ2)-K*IR2 ;
   if(IDUM2<0) IDUM2=IDUM2+IM2 ;
   int J=1+IY/NDIV ;
   IY=IV(J-1)-IDUM2 ;
   IV(J-1)=IDUM ;
   if (IY<1) IY = IY+IMM1 ;
   value=AM*IY ;
   if (value>RNMX) value = RNMX ;
   
}

//----------------------------------------------------------------------
double
MAC_Randomizer:: item( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Randomizer:: item" ) ;

   double result = value ;
   
   MAC_CHECK_POST( 0.0 <= result ) ;
   MAC_CHECK_POST( result <= 1.0 ) ;
   
   return( result ) ;
}


//----------------------------------------------------------------------
void
MAC_Randomizer:: build_permutation( size_t_vector& vec ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Randomizer:: build_permutation" ) ;
   size_t size = vec.size() ;
   for( size_t j=0 ; j<size ; j++ )
   {
      vec(j)=j ;
   }
   
   for( size_t j=0 ; j<size-1 ; j++ )
   {
      size_t nj = size-j ;
      double r = item()*nj ;
      go_next() ;
      int ij = (int) r ;
      if( ij>0 )
      {
         size_t ja = j + ij ;
         MAC_CHECK( ja<size ) ;
         size_t tmp = vec(j) ;
         vec(j) = vec(ja) ;
         vec(ja) = tmp ;
      }
   }
   MAC_CHECK_POST( FORALL( (size_t i=0;i<vec.size();i++),
                           vec(i)<vec.size() ) ) ;
   MAC_CHECK_POST( FORALL( (size_t i=0;i<vec.size();i++),
                           FORALL( (size_t j=0;j<vec.size();j++),
                                   IMPLIES(i!=j,vec(i)!=vec(j) ) ) ) ) ;
   
}


