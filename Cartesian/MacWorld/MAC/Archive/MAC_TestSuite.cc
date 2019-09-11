/*
 *  Copyright 1995-2010 by IRSN
 *
 *  This software is an application framework, with a set of integrated  
 *  reusable components, whose purpose is to simplify the task of developing 
 *  softwares of numerical mathematics and scientific computing.
 * 
 *  This software is governed by the CeCILL-C license under French law and 
 *  abiding by the rules of distribution of free software. You can use, modify 
 *  and/or redistribute the software under the terms of the CeCILL-C license  
 *  as circulated by CEA, CNRS and INRIA at the following URL 
 *  "http://www.cecill.info". 
 *
 *  As a counterpart to the access to the source code and rights to copy,  
 *  modify and redistribute granted by the license, users are provided only 
 *  with a limited warranty and the software's author, the holder of the  
 *  economic rights, and the successive licensors have only limited liability. 
 *
 *  In this respect, the user's attention is drawn to the risks associated  
 *  with loading, using, modifying and/or developing or reproducing the  
 *  software by the user in light of its specific status of free software,
 *  that may mean that it is complicated to manipulate, and that also  
 *  therefore means that it is reserved for developers and experienced 
 *  professionals having in-depth computer knowledge. Users are therefore 
 *  encouraged to load and test the software's suitability as regards their 
 *  requirements in conditions enabling the security of their systems and/or 
 *  data to be ensured and, more generally, to use and operate it in the same 
 *  conditions as regards security. 
 *
 *  The fact that you are presently reading this means that you have had 
 *  knowledge of the CeCILL-C license and that you accept its terms.
 */

#include <MAC_TestSuite.hh>

#include <MAC_Error.hh>
#include <MAC_Exec.hh>
#include <MAC_Iterator.hh>
#include <MAC_List.hh>
#include <MAC_ObjectTest.hh>
#include <MAC_ModuleExplorer.hh>
#include <MAC_assertions.hh>

#include <iostream>
#include <sstream>

using std::endl ;

struct MAC_TestSuite_ERROR
{
   static void n0( std::string const& t_name ) ;
} ;

MAC_TestSuite const* MAC_TestSuite::PROTOTYPE = new MAC_TestSuite() ;

//-------------------------------------------------------------------------
MAC_TestSuite:: MAC_TestSuite( void )
//-------------------------------------------------------------------------
   : MAC_Application( "MAC_TestSuite" )
   , TESTS( 0 )
{
}

//-------------------------------------------------------------------------
MAC_TestSuite*
MAC_TestSuite:: create_replica( MAC_Object* a_owner,
                                MAC_ModuleExplorer const* exp ) const
//-------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_TestSuite:: create_replica" ) ;
   MAC_CHECK( create_replica_PRE( a_owner, exp ) ) ;

   MAC_TestSuite* result = new MAC_TestSuite( a_owner, exp ) ;

   MAC_CHECK( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//-------------------------------------------------------------------------
MAC_TestSuite:: MAC_TestSuite( MAC_Object* a_owner,
                               MAC_ModuleExplorer const* exp )
//-------------------------------------------------------------------------
   : MAC_Application( a_owner, exp )
   , TESTS( 0 )
{
   TESTS = MAC_List::create( this ) ;

   stringVector names( 0 ) ;

   if( exp->has_entry( "without_data_deck" ) )
   {
      names = exp->stringVector_data( "without_data_deck" ) ;
      for( size_t i=0 ; i<names.size() ; ++i )
      {
         TESTS->append( MAC_ObjectTest::object( names( i ) ) ) ;
      }
   }
   if( exp->has_module( "with_data_deck" ) )
   {
      MAC_ModuleExplorer* ee = exp->create_subexplorer( 0, "with_data_deck" ) ;
   
      ee->start_module_iterator() ;
      for( ; ee->is_valid_module() ; ee->go_next_module() )
      {
         MAC_ModuleExplorer* sexp = ee->create_subexplorer( 0 ) ;

         std::string const& t_name = sexp->string_data( "concrete_name" ) ;
         if( names.size()!=0 && names.has( t_name ) )
            MAC_TestSuite_ERROR::n0( t_name ) ;

         MAC_ObjectTest* ot = MAC_ObjectTest::object( t_name ) ;
         ot->set_data_deck_explorer( sexp ) ;
         TESTS->append( ot ) ;
         sexp->destroy() ;
      }
      ee->destroy() ; ee = 0 ;
   }
}

//-------------------------------------------------------------------------
MAC_TestSuite:: ~MAC_TestSuite( void )
//-------------------------------------------------------------------------
{
}

//-------------------------------------------------------------------------
void
MAC_TestSuite:: run( void )
//-------------------------------------------------------------------------
{   
   MAC_Iterator* it = TESTS->create_iterator( this ) ;
   
   for( it->start(); it->is_valid() ; it->go_next() )
   {
      MAC_ObjectTest* ut = static_cast<MAC_ObjectTest*>( it->item() ) ;
      MAC_CHECK( dynamic_cast<MAC_ObjectTest*>( it->item() ) != 0 ) ;

      ut->run_all_tests() ;
   }

   int exit_code = 0 ;
   if( !MAC_ObjectTest::tests_of_all_instances_are_successful() )
   {
      MAC_Exec::out() << endl << "!!!! Failure of some unit tests !!!!" << endl << endl ;
      exit_code = 5 ;
   }
   else
   {
      MAC_Exec::out() << endl << "---- Success of all unit tests ---- " << endl << endl ;
   }
   MAC_Exec::set_exit_code( exit_code ) ;
}

//internal--------------------------------------------------------------
void 
MAC_TestSuite_ERROR:: n0( std::string const& t_name )
//internal--------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << endl << "*** MAC_TestSuite error:" << endl << endl ;
   mesg << "    \"" << t_name << "\" cannot appear both:" << endl ;
   mesg << "       - in a submodule of MODULE with_data_deck" << endl ;
   mesg << "       - in the entry of keyword without_data_deck";
   MAC_Error::object()->raise_plain( mesg.str() ) ;
}
