#ifndef MAC_DOUBLE_COMPARATOR_HH
#define MAC_DOUBLE_COMPARATOR_HH

#include <MAC_Object.hh>

class MAC_ModuleExplorer ;
class MAC_ObjectRegister ;

/*
Server comparing double values.

FRAMEWORK INSTANTIATION

   1. Derive a concrete subclass, say MyComparator.
   2. Choose a name for MyComparator, say "my_comparator".
   3. Implement a private destructor.
   4. Declare all constructors private.
   5. Define the prototype to be registered :
      5.1 Implement a default constructor that initializes the
          `MAC_DoubleComparator::' subobject by calling
               `MAC_DoubleComparator( std::string const& )'
          with "my_comparator" as argument.
          Example of pseudo-code :
          | MyComparator:: MyComparator( void )
          |    : MAC_DoubleComparator( "my_comparator" ) {}
      5.2 Define and initialize a static instance by calling the default
          constructor.
             declaration (in the header file, eg MyComparator.hh) :
             | static MyComparator const* PROTOTYPE ;
             definition (in the implementation file, eg MyComparator.cc) :
             | MyComparator const* MyComparator::PROTOTYPE = new MyComparator() ;'
   6. Implement a private constructor that initializes the 
      `MAC_DoubleComparator::' subobject by calling
                 `MAC_DoubleComparator( MAC_Object* )'
      Example of pseudo-code :
      | MyComparator:: MyComparator( MAC_Object* a_owner,
      |                              MAC_ModuleExplorer const* exp )
      |    : MAC_DoubleComparator( a_owner ), ...
      | { ... }
   7. Implement the `::create_replica' method that allocates an object
      of type `MyComparator' initialized using the private constructor described
      above, and subsequently return a pointer to that object.
      Example of pseudo-code :
      | MAC_DoubleComparator const* MyComparator::create_replica(
      |                                   MAC_Object* a_owner,
      |                                   MAC_ModuleExplorer const* exp ) const
      | {
      |    MAC_LABEL( "MyComparator::create_replica" ) ;
      |    MAC_CHECK_PRE( create_replica_PRE( a_owner, exp ) ) ;
      |    MyComparator const* result = new MyComparator( a_owner, exp ) ;
      |    MAC_CHECK_POST( create_replica_POST( result, a_owner, exp ) ;
      |    return result ;
      | }
   8. Implement the `::three_way_comparison' method

PUBLISHED
*/

class MAC_DoubleComparator : public MAC_Object
{
   public: //-----------------------------------------------------------------

   //-- Instance delivery and initialization

      // Create and return an instance of `MAC_DoubleComparator::' according
      // to the data attainable by `exp'.
      static MAC_DoubleComparator const* make(
                                           MAC_Object* a_owner,
                                           MAC_ModuleExplorer const* exp ) ;

   //-- Comparison

      // if `x' equal to `y', 0 ; if smaller, <0 ; if greater, >0  
      virtual int three_way_comparison( double x, double y ) const = 0 ;
      
   protected: //--------------------------------------------------------------
      
      MAC_DoubleComparator( MAC_Object* a_owner ) ;

   //-- Plug in

      virtual ~MAC_DoubleComparator( void ) ;

      MAC_DoubleComparator( std::string const& name ) ;
    
      MAC_DoubleComparator( MAC_Object* a_owner,
                            MAC_ModuleExplorer const* exp ) ;

      virtual MAC_DoubleComparator const* create_replica(
                            MAC_Object* a_owner,
                            MAC_ModuleExplorer const* exp ) const = 0 ;

      bool is_a_prototype( void ) const ;

   //-- Preconditions, Postconditions, Invariant
      
      bool create_replica_PRE( MAC_Object* a_owner,
                               MAC_ModuleExplorer const* exp ) const ;
      bool create_replica_POST( MAC_DoubleComparator const* result,
				MAC_Object* a_owner,
				MAC_ModuleExplorer const* exp ) const ;
      
      virtual bool invariant( void ) const ;
      
   private: //----------------------------------------------------------------

      MAC_DoubleComparator( void ) ;
      MAC_DoubleComparator( MAC_DoubleComparator const& other ) ;
      MAC_DoubleComparator& operator=( MAC_DoubleComparator const& other ) ;
      
      static MAC_ObjectRegister* plugins_map( void ) ;

   //-- Attributes
      
      bool IS_PROTO ;
} ;

#endif
