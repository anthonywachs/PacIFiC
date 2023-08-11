#ifndef EXT_PETScScatter_HH
#define EXT_PETScScatter_HH

#include <LA_Scatter.hh>

#include <EXT_PETScAPI.hh>
#include <intVector.hh>
#include <size_t_vector.hh>

class EXT_PETScVector ;
class LA_SeqVector ;

/*
  `LA_Scatter::' server for `EXT_PETScImplementation::' objects.
  
  PUBLISHED
*/

class EXT_PETScScatter : public LA_Scatter
{
   public: //-----------------------------------------------------------

      static EXT_PETScScatter* create(
                     MAC_Object* a_owner,
                     EXT_PETScVector const* global_vector,
                     size_t_vector const& a_repatriated_items_table,
                     size_t_vector const& a_local_indices_table ) ;
      
   //-- Characteristics
      
      virtual LA_Implementation const* implementation( void ) const ;
      
      virtual size_t size( void ) const ;

      virtual size_t_vector const& repatriated_items( void ) const ;
      
      virtual size_t_vector const& local_indices( void ) const ;
      
      virtual MAC_DistributedPartition const* distribution( void ) const ;
      
   //-- Repatriate items
      
      virtual void get( LA_Vector const* source,
                        LA_SeqVector* dest ) const ;
      
      virtual void set( LA_SeqVector const* source,
                        LA_Vector* dest ) const ;
      
   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      EXT_PETScScatter( void ) ;
     ~EXT_PETScScatter( void ) ;
      EXT_PETScScatter( EXT_PETScScatter const& other ) ;
      EXT_PETScScatter& operator=( EXT_PETScScatter const& other ) ;

      EXT_PETScScatter( MAC_Object* a_owner,
                        EXT_PETScVector const* global_vector,
                        size_t_vector const& a_repatriated_items_table,
                        size_t_vector const& a_local_indices_table ) ;

   //-- Preconditions, Postconditions, Invariant
      
      virtual bool implementation_POST(
                            LA_Implementation const* result ) const ;
      
   //-- Attributes
      
      size_t const SIZE ;
      MAC_DistributedPartition* const DIST ;
      size_t_vector const NEEDED ;
      size_t_vector const LOCAL ;
      
      VecScatter SCATTER ;
      Vec SEQ ;
} ;

#endif
