#ifndef MAC_OBJECT_HH
#define MAC_OBJECT_HH

#include <iosfwd>
#include <string>

class MAC_ListIdentity ;
class MAC_Module ;
class MAC_ObjectReader ;
class MAC_ObjectWriter ;

/*
Objects of dynamic storage duration,
   - that are referred to, accessed and manipulated exclusively 
     through pointers ;
   - that may be compared according to a total order relation ;
   - that may be hashed into a integer index, for use as keys in hash tables ;
   - whose lifetime is managed with the Ownership Method.
Any developper-written class publicly inherit from MAC_Object 

Each sub-object of type MAC_Object has one owner (possibly NULL).
The owner is set at creation and cannot be modified unless it is NULL. 
If the owner is NULL, the complete object of self is terminated, by
calling explicitely the destroy() method.
Otherwise, the complete object of self is terminated when the owner 
itself is terminated. 
When the complete object of self is terminated, all objects for which
self is the owner are also terminated. 

FRAMEWORK INSTANTIATION
   1. Derive a subclass.
   2. If concrete, implement a private destructor. 
      If abstract, implement a virtual protected destructor.
   3. If concrete, declare all constructors private.
      If abstract, declare the implemented constructors protected and
      declare all other constructors private.
   4. All implemented constructors should initialize the `MAC_Object'
      subobject by calling
         `MAC_Object( MAC_Object*)'
   5. If concrete, implement one or more static methods returning instances.
*/


class MAC_Object
{
   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization(0.1)

      // Create and a clone of `self'.
      // IMPLEMENTATION : raise a fatal error.
      virtual MAC_Object* create_clone( MAC_Object* a_owner ) const ;

      // Perform updating to maintain consistency between `self' and a 
      // set of related objects (observer pattern).
      // IMPLEMENTATION : raise a fatal error.
      virtual void update( void ) ;
      
   //-- Termination(0.2)

      // Terminate `self'. (Pseudo-destructor that has access to the real 
      // destructor : a call to destroy leads to the destruction of `self').
      void destroy( void ) const ;

      // Terminate `a_possession'.
      void destroy_possession( MAC_Object const* a_possession ) ;

   //-- Identification(0.3)

      // address
      size_t address( void ) const ;

   //-- Characteristics(0.5)

      // name of the class of which `self' is an instance
      std::string const& type_name( void ) const ;
            
      // owner (possibly NULL)
      MAC_Object const* owner( void ) const ;

      // Is `other' in the tree of the owners of `self' ?
      bool is_under_ownership_of( MAC_Object const* other ) const ;

      // Is type of `self' identical to type of `other' ?
      bool same_type( MAC_Object const* other ) const ;

   //-- Characteristics setting(0.6)
      
      // Make `a_owner' the owner of self.
      void set_owner( MAC_Object* a_owner ) ;

      // Make `a_owner' the owner of some possession `a_possession'.
      void change_owner( MAC_Object* a_owner, MAC_Object* a_possession ) ;

   //-- Comparison(0.7)

      // Is `other' comparable to `self' ? 
      // IMPLEMENTATION : `::same_type(other)'
      virtual bool comparable( MAC_Object const* other ) const ;

      // Is `other' equal to `self' ? 
      // IMPLEMENTATION : `::has_same_address(other)'
      virtual bool is_equal( MAC_Object const* other ) const ;

      // if `self' equal to `other', O ; if smaller, <0 ; if greater, >0 
      // IMPLEMENTATION : address_comparison(other)
      virtual int three_way_comparison( MAC_Object const* other ) const ;

      // hash code value
      // IMPLEMENTATION : address()
      virtual size_t hash_code( void ) const ;

      // Is `other' identical to `self' (same address) ?
      bool has_same_address( MAC_Object const* other ) const ;

      // if `self''s address equal to the address of `other', O ; 
      // if smaller, -1 ; if greater, +1 */
      int address_comparison( MAC_Object const* other ) const ;

   //-- Persistence(900.0)

      // Use `writer' to store `self' so that ist can be retrieved
      // with `::restore_state'.
      // IMPLEMENTATION : raise a fatal error.
      virtual void save_state( MAC_ObjectWriter* writer ) const ;

      // Retrieve `self' from `reader'.
      // IMPLEMENTATION : raise a fatal error.
      virtual void restore_state( MAC_ObjectReader* reader ) ;

      // Perform updating to restore the consistency between `self' and a 
      // set of related objects (observer pattern). To be used instead
      // of `::update()' in `::restore_state()' because, in the course
      // of a restart, the operations performed in 
      // `::update_for_restore_state()' are closely related to those performed
      // in `::restore_state()'.
      // IMPLEMENTATION : call `::update()'.
      virtual void update_for_restore_state( void ) ;
      
   //-- Input - Output(1001.0)

      // Write text for debugging to `os'. 
      // IMPLEMENTATION : write the type name.
      virtual void display_info( std::ostream& os, size_t indent_width ) const ;

      // Write text to `os' with `indent_width' indentation.
      // IMPLEMENTATION : `os' is left unchanged.
      virtual void print( std::ostream& os, size_t indent_width ) const ;

   //-- Hidden

      // Total number of MAC_Object actually allocated */
      static int GetNumberOf_MAC_objects( void ) ;

      // for debugging purposes only
      static std::ostream& TraceRemainingObjects( std::ostream& out ) ;

      // catch creation and deletion for particular object.
      static void catch_object( MAC_Object const* obj ) ;

      // catch creation for particular object by its creation rank.
      static void catch_object_by_rank( size_t rank ) ;

      static void start_trace_allocating( void ) ;
      static bool trace_allocating( void ) ;
      static void stop_trace_allocating( void ) ;
      static void trace_not_destroyed_object( std::ostream& out ) ;

   protected: //--------------------------------------------------------

   //-- Plug in(1002.0)

      virtual ~MAC_Object( void ) = 0 ;

      // Construction of an instance whose owner is `a_owner'
      MAC_Object( MAC_Object* a_owner ) ;

   //-- Preconditions, Postconditions, Invariant(9999.0)

      virtual bool invariant( void ) const ;

      bool create_clone_POST( MAC_Object const* result, 
                              MAC_Object const* a_owner ) const ;

      virtual bool save_state_PRE( MAC_ObjectWriter const* writer ) const ;
      virtual bool save_state_POST( MAC_ObjectWriter const* writer ) const ;

      virtual bool restore_state_PRE( MAC_ObjectReader const* reader ) const ;
      virtual bool restore_state_POST( MAC_ObjectReader const* reader ) const ;
      

      virtual bool is_equal_PRE( MAC_Object const* other ) const ;
      virtual bool is_equal_POST( bool result, 
                                  MAC_Object const* other ) const ;

      virtual bool three_way_comparison_PRE( MAC_Object const* other ) const ;
      virtual bool three_way_comparison_POST( int result, 
                                              MAC_Object const* other ) const ;

   private: //----------------------------------------------------------

      MAC_Object( void ) ;
      MAC_Object( MAC_Object const& other ) ;
      MAC_Object& operator=( MAC_Object const& other ) ;

      // Inserts a aComponent to the list of owned components.
      void insert_possession( MAC_Object* obj ) ;

   //-- Class attributes
      
      // Number of allocated MAC_Object instances remaining on the heap. 
      static size_t ALLOCATED ;
      
   //-- Attributes
      
      MAC_Object* MY_OWNER ;
      MAC_ListIdentity* POSSESSIONS ;

} ;


#endif

