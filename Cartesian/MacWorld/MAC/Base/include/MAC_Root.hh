#ifndef MAC_ROOT_HH
#define MAC_ROOT_HH

#include <MAC_Object.hh>

#include <string>

class MAC_Vector ;

/*
Object devoted to be the upper-most owner in the ownership method
of lifetime management.

Implemented as a singleton.

The programm execution consists of five stages :
   1. Initial stage (Big-Bang time) : all statics are initialized and
      the only instance of MAC_Root is created.
   2. The data deck is read and stored in memory.
   3. An instance of a concrete subclass of MAC_Application is created.
   4. Program core execution : the program execution proceeds by performing
      its specific tasks. In particular, objects are created and organized
      into ownership trees whose root node is either the unique instance
      of MAC_Root or the NULL object.
   5. Final stage : termination of the only instance of MAC_Root, leading to 
      the termination of all objects belonging to a ownership tree whose
      root node is not the NULL object.
   The only instance of MAC_Root is mainly involved in the initial and the
   final stages, and at any time an object is created to exist until the end
   of the program execution
*/

class MAC_Root : public MAC_Object
{
   public: //--------------------------------------------------------------

   //-- Instance delivery and initialization
      
      // the only instance 
      static MAC_Object* object( size_t i=0 ) ;

   //-- Termination
      
      // Terminate self: if i<j, object(i) is destroyed before object(j).
      static void cleanup( void ) ;

   protected: //-----------------------------------------------------------

   private: //-------------------------------------------------------------
  
      MAC_Root( void ) ;
     ~MAC_Root( void ) ;
      MAC_Root( MAC_Root const& other ) ;
      MAC_Root& operator=( MAC_Root const& other ) ;

      MAC_Root( MAC_Object* a_owner ) ;

      static MAC_Vector* objects( bool destroy = false ) ;

} ;

#endif
