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

#include <MAC_SimplifyPattern.hh>

#include <MAC_Exec.hh>
#include <MAC_Error.hh>
#include <MAC_Module.hh>
#include <MAC_ModuleExplorer.hh>
#include <MAC_ModulePattern.hh>
#include <MAC_assertions.hh>

#include <iostream>

MAC_SimplifyPattern const*
MAC_SimplifyPattern::PROTOTYPE = new MAC_SimplifyPattern() ;

//----------------------------------------------------------------------
MAC_SimplifyPattern:: MAC_SimplifyPattern( void )
//----------------------------------------------------------------------
   : MAC_Application( "simplify_pattern" )
   , PATTERN()
{
}

//----------------------------------------------------------------------
MAC_SimplifyPattern* 
MAC_SimplifyPattern:: create_replica(
              MAC_Object* a_owner, MAC_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_SimplifyPattern:: create_replica" ) ;
   MAC_CHECK( create_replica_PRE( a_owner, exp ) ) ;

   MAC_SimplifyPattern* result =
      new MAC_SimplifyPattern( a_owner,
                                exp->string_data( "pattern_filename" ) ) ;

   MAC_CHECK_POST( create_replica_POST( result, a_owner, exp ) ) ;
   return result ;
}

//----------------------------------------------------------------------
MAC_SimplifyPattern* 
MAC_SimplifyPattern:: create_replica_from_args(
                         MAC_Object* a_owner, stringVector& args ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_SimplifyPattern:: create_replica_from_args" ) ;
   
   if( args.size()!=1 )
      MAC_Error::object()->raise_plain( "usage : <pattern filename>" ) ;
   
   MAC_SimplifyPattern* result =
      new MAC_SimplifyPattern( a_owner, args(0) ) ;
   args.remove_at(0) ;

   MAC_CHECK_POST( create_replica_from_args_POST( result, a_owner, args ) ) ;
   return result ;
}

//----------------------------------------------------------------------
void
MAC_SimplifyPattern:: run( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_SimplifyPattern:: run" ) ;
   
   MAC_ModulePattern::open_pattern_base( PATTERN ) ;
   MAC_ModulePattern::expand_then_simplify() ;
   MAC_ModulePattern::save_pattern() ;
}

//----------------------------------------------------------------------
MAC_SimplifyPattern:: MAC_SimplifyPattern( MAC_Object* a_owner,
                                             std::string const& pattern_file )
//----------------------------------------------------------------------
   : MAC_Application( a_owner, 0 )
   , PATTERN( pattern_file )
{
}

//----------------------------------------------------------------------
MAC_SimplifyPattern:: ~MAC_SimplifyPattern( void )
//----------------------------------------------------------------------
{
   if( is_a_prototype() )
   {
      PROTOTYPE = 0 ;
   }
}
