#ifndef MAC_ASSERTIONS_HH
#define MAC_ASSERTIONS_HH

#include <string>

class MAC_Assertion
{
   public: //------------------------------------------------------------

      enum CheckType
      {          None          = 0,
         Precondition  = 1,
         Postcondition = 2,
         Invariant     = 4,
         Check         = 8,
         Objects       = 16
      } ;

      static MAC_Assertion& object( void ) ;

      ~MAC_Assertion( void ) ;

      static void add_handled_check( CheckType a_chec ) ;
      static bool is_handling_check( CheckType some_check ) ;

      static bool is_checking( void ) ;

      static bool test_implement_handling( const char* file,
                                           int line,
                                           std::string const& reason ) ;
      static bool action( const char* file, int line, const char* text ) ;
      static bool do_eval( bool& shortCut ) ;

      friend MAC_Assertion const& operator!( MAC_Assertion const& a ) ;
      friend bool operator&&( bool left, MAC_Assertion const& a ) ;
      friend bool operator||( bool left, MAC_Assertion const& a ) ;
      operator bool( void ) const ;

      static bool& push_bool( void ) ;
      static bool  pop_bool( void ) ;

      static CheckType current_check ;
      static bool result ;
      static bool negation ;
      static bool short_cut ;

   protected: //---------------------------------------------------------

   private: //-----------------------------------------------------------

      MAC_Assertion( void ) ;
      MAC_Assertion( MAC_Assertion const& other ) ;
      MAC_Assertion& operator=( MAC_Assertion const& other ) ;

      static MAC_Assertion unique_instance ;
      static bool eval ;
      static size_t const MAXBOOL = 256 ;
      static bool bool_table[ MAXBOOL ] ;
      static size_t nb_bool ;
      static int checking_level ;
} ;


class MAC_Marker
{
   public: //-----------------------------------------------------------

      MAC_Marker( char const* name ) ;
     ~MAC_Marker( void ) ;

      static bool is_collective( int line ) ;
      
      static size_t nb_labels( void ) ;
      static char const* label( size_t i ) ;

   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      MAC_Marker( void ) ;
      MAC_Marker( MAC_Marker const& other ) ;
      MAC_Marker const& operator=( MAC_Marker const& other ) ;

      static const char* ring[] ;
      static size_t ring_pos ;
} ;


//----------------------------------------------------------------------

#define asstLoopTest(asstFor,asstAll,predicate)  \
        MAC_Assertion::object() ; \
        { \
           if( MAC_Assertion::do_eval( MAC_Assertion::push_bool() ) ) \
           { \
              MAC_Assertion::push_bool() = MAC_Assertion::negation ; \
              MAC_Assertion::result = asstAll ; \
              for asstFor \
              { \
                 MAC_Assertion::result = predicate ; \
                 if( MAC_Assertion::result != asstAll ) break ; \
              } \
              if( MAC_Assertion::pop_bool() ) \
                 MAC_Assertion::result = !MAC_Assertion::result ;  \
          } \
          MAC_Assertion::short_cut = MAC_Assertion::pop_bool() ; \
        } \
        MAC_Assertion::result = MAC_Assertion::short_cut ? \
        MAC_Assertion::result : MAC_Assertion::result

#define MAC_ASSERTCOND(predicate,text) \
   MAC_Assertion::result = predicate, \
   MAC_Assertion::result || (MAC_Assertion::action(__FILE__,__LINE__,text))

#define MAC_CHECKUNCOND(predicate,check_type,text) \
   if( MAC_Assertion::current_check == MAC_Assertion::None ) \
   { \
      MAC_Assertion::current_check = check_type ; \
      MAC_ASSERTCOND(predicate,text) ; \
      MAC_Assertion::current_check = MAC_Assertion::None ; \
   }

#define MAC_CHECKCOND(predicate,check_type,text) \
   if( MAC_Assertion::is_handling_check( check_type ) && \
       ( MAC_Assertion::current_check == MAC_Assertion::None ) ) \
   { \
      MAC_Assertion::current_check = check_type ; \
      MAC_ASSERTCOND(predicate,text) ; \
      MAC_Assertion::current_check = MAC_Assertion::None ; \
   }

#if !defined(LEVEL) || LEVEL<0 || LEVEL>2
error
Macro LEVEL must be set when compiling Pelicans (opt: -DLEVEL=<value>).
LEVEL value meaning is :
0 : no assertion checking
1 : only preconditions will be checked (recommended).
2 : preconditions, postconditions, invariant and simple checks will be tested
    depending on dynamic level.
#endif

//----------------------------------------------------------------------
// PUBLIC MACROS
//----------------------------------------------------------------------

#define MAC_TEST_IMPLEMENTATION(X) \
   if(X) MAC_Assertion::test_implement_handling( __FILE__, __LINE__, #X )

#define MAC_ASSERT(X) MAC_ASSERTCOND(X,#X)

//--------
#if LEVEL>=2
#define MAC_CHECK_POST(X)  MAC_CHECKCOND(X,MAC_Assertion::Postcondition,#X)
#define OLD(asstname) old_##asstname
#define MAC_SAVEOLD(assttype,asstname,args) assttype old_##asstname = args
#define MAC_CHECK_INV(X) MAC_CHECKCOND(X,MAC_Assertion::Invariant,#X)
#define MAC_CHECK(X) MAC_CHECKCOND(X,MAC_Assertion::Check,#X)
#else
#define MAC_CHECK_POST(X) {}
#define OLD(asstname) {}
#define MAC_SAVEOLD(assttype,asstname,args) {}
#define MAC_CHECK_INV(X) {}
#define MAC_CHECK(X) {}
#endif
//--------

//--------
#if LEVEL>=1
#define MAC_CHECK_PRE(X) MAC_CHECKUNCOND(X,MAC_Assertion::Precondition,#X)
#define MAC_LABEL(X) MAC_Marker aSpy(X)
#define MAC_CHECK_COLLECTIVE(X) if(X){MAC_ASSERTCOND(MAC_Marker::is_collective(__LINE__),"unsynchronized collective operation") ; }
#else
#define MAC_CHECK_PRE(X) {}
#define MAC_LABEL(X) {}
#define MAC_CHECK_COLLECTIVE(X) {}
#endif
//--------

#define FORMAL(X) true
#define FORALL(asstFor,predicate)  asstLoopTest(asstFor,true,predicate)
#define EXISTS(asstFor,predicate)  asstLoopTest(asstFor,false,predicate)
#define IMPLIES(pred,res)          ( !(pred) || (res) )
#define EQUIVALENT(pred,res)       ( (pred) && (res) ) || ( !(res) && !(pred) )

#undef LEVEL

#ifndef OUTLINE
   #include <MAC_assertions.icc>
#endif

#endif
