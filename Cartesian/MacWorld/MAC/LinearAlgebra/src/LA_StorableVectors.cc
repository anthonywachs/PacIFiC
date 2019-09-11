#include <LA_StorableVectors.hh>

#include <MAC.hh>
#include <MAC_Communicator.hh>
#include <MAC_Exec.hh>
#include <MAC_Int.hh>
#include <MAC_assertions.hh>

#include <LA_Vector.hh>

#include <iostream>


//-------------------------------------------------------------------------
LA_StorableVectors* LA_StorableVectors:: create( MAC_Object* a_owner )
//-------------------------------------------------------------------------
{
   MAC_LABEL( "LA_StorableVectors:: create" ) ;

   LA_StorableVectors* result = new LA_StorableVectors( a_owner ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}




//-------------------------------------------------------------------------
LA_StorableVectors:: LA_StorableVectors( MAC_Object* a_owner )
//-------------------------------------------------------------------------
   : MAC_Object( a_owner )
{
   MAC_CHECK_INV( invariant() ) ;
}




//-------------------------------------------------------------------------
LA_StorableVectors:: ~LA_StorableVectors( void )
//-------------------------------------------------------------------------
{
   STORED_VECTORS.clear();
}



   
//----------------------------------------------------------------------
void
LA_StorableVectors:: save_state( MAC_ObjectWriter* writer ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_StorableVectors:: save_state" ) ;
   MAC_CHECK_INV( invariant() ) ;

   writer->start_new_object( "LA_StorableVectors" ) ;
   
   writer->add_entry( "number_of_vectors", MAC_Int::create( 0,
   	STORED_VECTORS.size() ) ) ;

   // Saving vectors
   for ( list<LA_Vector*>::const_iterator
   		il=STORED_VECTORS.begin(); il!=STORED_VECTORS.end(); il++)
     (*il)->save_state( writer ) ;   
   
   writer->finalize_object() ;
}




//----------------------------------------------------------------------
void
LA_StorableVectors:: restore_state( MAC_ObjectReader* reader )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_StorableVectors:: restore_state" ) ;
   MAC_CHECK_INV( invariant() ) ;
   
   reader->start_object_retrieval( "LA_StorableVectors" ) ;

   // Retrieving stored data
   size_t const nvec = 
   	size_t(reader->data_of_entry( "number_of_vectors" )->to_int()) ;

   // Retrieving vector values
   std::string stored_name; 
   bool found = false; 
     
   // If a vector is needed but was not stored, it is initialized
   // to the values defined in the application (most likely 0)
   while( reader->next_object_class_in_current_module() == "LA_Vector" )
   {
     stored_name = reader->next_object_name_in_current_module();
     found = false;
     for ( list<LA_Vector*>::const_iterator
	il=STORED_VECTORS.begin(); il!=STORED_VECTORS.end() && !found; il++)
     {
       if ( (*il)->name() == stored_name )
       { 
	 found = true;	 
	 (*il)->restore_state( reader ) ;	 
       }
     }
     
     // If the vector was stored but is not used in the new application
     // We still need to do a "fake" read for the reader to proceed properly
     // i.e. essentially to be able to move to the next object
     if ( !found )
     {
       LA_Vector:: read_state_nonrestored( reader );
       if ( MAC_Exec::communicator()->rank() == 0 )
         MAC::out() << "      Vector \"" << stored_name << "\" stored but not"
	 	<< " reloaded" << std::endl;
     }       
   }      

   reader->end_object_retrieval() ;
}




//----------------------------------------------------------------------
void
LA_StorableVectors:: add_vector_to_store( LA_Vector* vec )
//----------------------------------------------------------------------
{
   MAC_LABEL( "LA_StorableVectors:: add_vector_to_store" ) ;
   MAC_CHECK_INV( invariant() ) ;
   
   STORED_VECTORS.push_back(vec);
}




//-------------------------------------------------------------------------
bool
LA_StorableVectors:: invariant( void ) const
//-------------------------------------------------------------------------
{
   MAC_ASSERT( MAC_Object::invariant() ) ;
   return( true ) ;
}
