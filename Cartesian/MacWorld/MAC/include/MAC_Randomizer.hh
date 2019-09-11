#ifndef MAC_Randomizer_HH
#define MAC_Randomizer_HH

#include <MAC_Object.hh>
#include <doubleVector.hh>

class size_t_vector ;

// Random number generator.
//     LONG PERIOD (>2.D18) RANDOM NUMBER GENERATOR OF L'ECUYER
//     WITH BAYS-DURHAM SHUFFLE AND ADDED SAFEGUARDS. RETURNS A
//     UNIFORM RANDOM DEVIATE BETWEEN 0.0 AND 1.0. CALL WITH IDUM
//     A NEGATIVE INTEGER TO INITIALIZE. THERAFTER, DO NOT ALTER
//     IDUM BETWEEN SUCCESSIVE DEVIATES IN A SEQUENCE. RNMX SHOULD
//     APPROXIMATE THE LARGEST FLOATING VALUE THAT IS LESS THAN 1.

class MAC_Randomizer : public MAC_Object
{
   public: //---------------------------------------------------------------

   //-- Instance delivery and initialization

      static MAC_Randomizer* create( MAC_Object* a_owner,
                                    int series ) ;

      virtual MAC_Randomizer* create_clone( MAC_Object* a_owner ) const ;

  //-- Cursor movement

      // random value
      virtual double item( void ) const ;

      // Start random iterator.
      virtual void start( void ) ;
  
      // Go to next value.
      virtual void go_next( void ) ;

      void build_permutation( size_t_vector& vec ) ;
      
   protected: //------------------------------------------------------------

      virtual ~MAC_Randomizer( void ) ;

      MAC_Randomizer( MAC_Object* a_owner, int series ) ;
      
   private: //--------------------------------------------------------------

      MAC_Randomizer( void ) ;
      MAC_Randomizer( MAC_Randomizer const& other ) ;
      MAC_Randomizer const& operator=( MAC_Randomizer const& other ) ;

      int IDUM ;
      int IDUM2 ;
      int IY ;
      int my_series ;
      doubleVector IV ;
      double value ;
      
      
} ;

#endif
