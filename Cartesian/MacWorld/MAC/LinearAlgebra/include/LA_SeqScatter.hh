#ifndef LA_SeqScatter_HH
#define LA_SeqScatter_HH

#include <LA_Scatter.hh>

#include <size_t_vector.hh>

class LA_SeqVector ;

/*
  `LA_Scatter::' server for `LA_SeqImplementation::' objects.

  PUBLISHED
*/

class LA_SeqScatter : public LA_Scatter
{
   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization
      
      static LA_SeqScatter* create(
                    MAC_Object* a_owner,
                    size_t a_nb_rows,
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

      LA_SeqScatter( void ) ;
     ~LA_SeqScatter( void ) ;
      LA_SeqScatter( LA_SeqScatter const& other ) ;
      LA_SeqScatter& operator=( LA_SeqScatter const& other ) ;

      LA_SeqScatter( MAC_Object* a_owner,
                     size_t a_nb_rows,
                     size_t_vector const& a_repatriated_items_table,
                     size_t_vector const& a_local_indices_table ) ;

   //-- Preconditions, Postconditions, Invariant
      
      virtual bool implementation_POST(
                            LA_Implementation const* result ) const ;
      
   //-- Attributes
      
      size_t NB_ROWS ;
      size_t_vector const NEEDED ;
      size_t_vector const LOCAL ;
      mutable MAC_DistributedPartition* DIST ;
} ;

#endif
