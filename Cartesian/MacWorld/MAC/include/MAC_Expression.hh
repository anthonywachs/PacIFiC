#ifndef MAC_EXPRESSION_HH
#define MAC_EXPRESSION_HH

#include <MAC_Data.hh>

class MAC_ContextSimple ;
class MAC_Iterator ;
class MAC_ObjectRegister ;
class MAC_Sequence ;

/*
Context dependant data.

An expression is identified by :
   - its name, and
   - its list of arguments, where each argument is itself a `MAC_Data::'
     object.

IMPLEMENTATION :
   A register is maintained that maps the name of each concrete subclass
   with a prototype of this subclass. That prototype is an instance
   with a NULL argument list (and is thus in an invalid state) that serves
   as a model to create and initialize instances with the `create_replica'
   method.

FRAMEWORK INSTANTIATION :
   1. Derive a concrete subclass, say MyExp (derivation of 
      abstract subclasses is improbable although possible). 
   2. Choose a name for `MyExp', say "my_exp". 
   3. Implement a private destructor.
   4. Declare all constructors private.
   5. Define the prototype to be registered :
      5.1 Implement a contructor (eg the default constructor) that initializes
          the `MAC_Expression::' subobject by calling
                   `MAC_Expression( std::string const& )'
          with "my_exp" as argument.
          Example of pseudo-code :
          | MyExp:: MyExp( void ) : MAC_Expression( "my_exp" ) {}
      5.2 Define and initialize a static instance by calling the default 
          constructor.
             declaration :
             |    static MyExp const* PROTOTYPE ;
             definition :
             |    MyExp const* MyExp::PROTOTYPE = new MyExp() ;
   6. Implement the `create_replica' methods, that calls a constructor.
   7. Implement the constructor called by `create_replica'. The 
      `MAC_Expression' subobject is initialized by calling     
      `MAC_Expression( MAC_Object*, std::string const&, MAC_Sequence const* )'
      with "my_exp" as a second argument.
   8. Implement all necessary method evaluating the type and 
      value of MyExp (eg `to_double', `to_int', `data_type' for an expression
      whose value can be of type int and double).
*/

class MAC_Expression : public MAC_Data
{
   public: //-------------------------------------------------------

   //-- Instance delivery and initialization

      // Create, initialize and return an instance.
      static MAC_Expression* create( MAC_Object* a_owner,
                                     std::string const& a_name,
                                     MAC_Sequence const* argument_list,
                                     std::string const& a_comment = "" ) ;
      
      virtual MAC_Expression* create_clone( MAC_Object* a_owner ) const ;

   //-- Identification

      // name
      std::string const& name( void ) const ;
      
   //-- Context
      
      virtual void declare( MAC_List* lst ) const ;

      // IMPLEMENTATION :
      //    `MAC_Data::context_has_required_variables'( `ct' ) for all
      //    the data of `self' ?
      virtual bool context_has_required_variables( 
                                           MAC_Context const* ct ) const ;

   //-- Value

      // IMPLEMENTATION :
      //    `MAC_Data::'value_can_be_evaluated( `ct' ) for all
      //    the data of `self' ?
      virtual bool value_can_be_evaluated( MAC_Context const* ct ) const ;

      // IMPLEMENTATION :
      //    vector of all `MAC_Data::'undefined_variables( `ct' ) for all
      //    the data of `self'
      virtual stringVector const& undefined_variables(
                                           MAC_Context const* ct ) const ;

   //-- Formal calculus

      virtual MAC_Data* create_derivative( MAC_Object* a_owner,
                                           MAC_Variable const* var,
                                           MAC_Context const* ct  ) const ;

   //-- Input - Output

      virtual void print( std::ostream& os, size_t indent_width ) const ;

      static void print_prototypes( std::ostream& os, size_t indent_width ) ;

      bool external_brackets_are_set( void ) const ;
      void set_external_brackets( void ) ;
      void unset_external_brackets( void ) ;
      
   //-- Registered instances (1010.0)
      
      static stringVector const& registered_expressions( void ) ;
      
      static std::string const& usage_of( std::string const& a_name ) ;
      
      static bool valid_arguments_of( std::string const& a_name,
                                      MAC_Sequence const* argument_list ) ;
      
   protected: //-------------------------------------------------------
      
   //-- Plug in

      virtual ~MAC_Expression( void ) ;

      // for prototype registration only
      MAC_Expression( std::string const& a_name ) ;

      // Construction of an instance called `a_name', whose
      // owner is `a_owner' with the items of `argument_list' as arguments.
      MAC_Expression( MAC_Object* a_owner,
                      std::string const& a_name,
                      MAC_Sequence const* argument_list ) ;

      // Is `self' a prototype ?
      bool is_a_prototype( void ) const ;
      
      virtual MAC_Expression* create_replica( 
                              MAC_Object* a_owner,
                              MAC_Sequence const* argument_list ) const  = 0 ;

      
   //-- Characteristics

      // user documentation of `self'
      virtual std::string const& usage( void ) const = 0 ;
      
      // Is the list arguments stored in `some_arguments' valid ? 
      virtual bool valid_arguments(
                              MAC_Sequence const* some_arguments ) const = 0 ;

      // comment defined for `self'
      std::string const& comment( void ) const ;

   //-- Formal calculus

      virtual MAC_Data* create_non_const_simplification(
                                           MAC_Object* a_owner ) const ;
      virtual bool is_raw_data( void ) const ;
     
      // Create operator-specific omptimization if any.
      virtual MAC_Data* create_operator_simplification( MAC_Object* a_owner ) ;

   //-- Arguments

      // number of arguments
      size_t nb_arguments( void ) const ;
      
      // `idx'-th argument
      MAC_Data const* arg( size_t idx ) const ;

      // `idx'-th item of `some_arguments'
      static MAC_Data const* extract_arg( MAC_Sequence const* some_arguments,
                                          size_t idx ) ;

   //-- Error

      void raise_error( std::string const& message ) const ;
      
   //-- Preconditions, Postconditions, Invariant      

      virtual bool create_replica_PRE( 
                              MAC_Object const* a_owner,
                              MAC_Sequence const* argument_list ) const ;

      virtual bool create_replica_POST(
	                      MAC_Expression const* result,
                              MAC_Object const* a_owner,
                              MAC_Sequence const* argument_list ) const ;

      virtual bool valid_arguments_PRE(
                              MAC_Sequence const* some_arguments ) const ;

      virtual bool create_operator_simplification_POST(
                              MAC_Object const* a_owner,
                              MAC_Data const* result ) const ;
      
      virtual bool invariant( void ) const ;
      
   private: //-------------------------------------------------------

      MAC_Expression( void ) ;
      MAC_Expression( MAC_Expression const& other ) ;
      MAC_Expression& operator=( MAC_Expression const& other ) ;

      static stringVector& plugins_name( void ) ;
      static MAC_ObjectRegister* plugins_map( void ) ;
      

   //-- Attributes

      std::string const NAME ;
      bool HAS_BRACKETS ;
      std::string COMMENT ;
      MAC_Sequence const* const ARGUMENTS ; // sequence of MAC_Data*
      mutable MAC_Iterator* IT ;            // iterator on ARGUMENTS
} ;

#endif
