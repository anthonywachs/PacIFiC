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

#ifndef MAC_CHECK_HH
#define MAC_CHECK_HH

#include <MAC_Application.hh>
#include <string>
#include <stringVector.hh>

class MAC_ModuleExplorer ;

/*
FRAMEWORK INSTANTIATION
*/

class MAC_Check : public MAC_Application
{

   public: //-----------------------------------------------------------

   //-- Program core execution

      virtual void run( void ) ;

   protected: //--------------------------------------------------------

   //-- Plug in
      
      virtual ~MAC_Check( void ) ; 
      
      MAC_Check( std::string const& a_name ) ;
      
      MAC_Check( MAC_Object* a_owner,
                 std::string const& a_name,
                 MAC_ModuleExplorer const* exp ) ;
      
      MAC_Check( MAC_Object* a_owner,
                 std::string const& a_name,
                 stringVector& args ) ;
      
      virtual MAC_Check* create_replica( 
                                     MAC_Object* a_owner,
				     MAC_ModuleExplorer const* exp ) const ;
      
      virtual MAC_Check* create_replica_from_args( 
                                     MAC_Object* a_owner,
                                     stringVector& args ) const ;
      
   //-- Checks

      virtual MAC_ModuleExplorer* do_check( MAC_Object* a_owner,
                                            MAC_ModuleExplorer const* exp ) const ;

   //-- Preconditions, Postconditions, Invariant    

      virtual bool do_check_PRE( MAC_Object* a_owner,
                                 MAC_ModuleExplorer const* exp ) const ;

      virtual bool do_check_POST( MAC_Object* a_owner,
                                  MAC_ModuleExplorer const* result ) const ;

   private: //----------------------------------------------------------

      MAC_Check( void ) ;
      MAC_Check( MAC_Check const& other ) ;
      MAC_Check& operator=( MAC_Check const& other ) ;

      void process( std::string const& file_to_parse ) ;
      
   //-- Class attribute
      
      static MAC_Check const* PROTOTYPE ;

   //-- Attribute
      
      std::string PATTERN ;
      stringVector MY_ARGS ;
      std::string FILE ;
      bool SILENT ;
      bool INTERACTIVE ;
      std::string NAME ;
      
      
} ;

#endif



