#ifndef MAC_TEST_SUITE_HH
#define MAC_TEST_SUITE_HH

#include <MAC_Application.hh>

#include <iosfwd>

class MAC_ObjectTest ;
class MAC_List ;

/*
Unit tests performed by a suite of `MAC_ObjectTest::' concrete subclasses.

PUBLISHED 
*/

class MAC_TestSuite : public MAC_Application
{
   public: //--------------------------------------------------------------

   //-- Program core execution

      virtual void run( void ) ;
      
   protected: //-----------------------------------------------------------
      
   private: //-------------------------------------------------------------

     ~MAC_TestSuite( void ) ;
      MAC_TestSuite( MAC_TestSuite const& other ) ;
      MAC_TestSuite& operator=( MAC_TestSuite const& other ) ;

      MAC_TestSuite( MAC_Object* a_owner, MAC_ModuleExplorer const* exp ) ;

   //-- Plug in

      MAC_TestSuite( void ) ;

      virtual MAC_TestSuite* create_replica( 
                                       MAC_Object* a_owner,
				       MAC_ModuleExplorer const* exp ) const ;

   //-- Class attributes

      static MAC_TestSuite const* PROTOTYPE ;

   //-- Attributes

      MAC_List* TESTS ;
} ;

#endif 
