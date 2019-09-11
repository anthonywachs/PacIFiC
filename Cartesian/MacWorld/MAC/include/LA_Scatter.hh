#ifndef LA_SCATTER_HH
#define LA_SCATTER_HH

#include <MAC_Object.hh>

class MAC_DistributedPartition ;
class size_t_vector ;
class LA_Implementation ;
class LA_SeqVector ;
class LA_Vector ;

/*
  Objects used to repatriate distributed items in local one.

  PUBLISHED
*/

class LA_Scatter : public MAC_Object
{
   public: //-----------------------------------------------------------

   //-- Characteristics

      // implementation indicator
      virtual LA_Implementation const* implementation( void ) const = 0 ;
      
      virtual size_t size( void ) const = 0 ;

      // table of items which needed to be repatriated
      virtual size_t_vector const& repatriated_items( void ) const = 0 ;
      
      // indices of `::repatriated_items' elements in the sequential vector
      virtual size_t_vector const& local_indices( void ) const = 0 ;
      
      // distribution of vector rows over the processes
      virtual MAC_DistributedPartition const* distribution( void ) const = 0 ;

   //-- Repatriate items
      
      // Repatriate distributed items in distributed `source' to local `dest'.
      virtual void get( LA_Vector const* source,
                        LA_SeqVector* dest ) const = 0 ;
      
      // Repatriate distributed items in local `source' to distributed `dest'.
      virtual void set( LA_SeqVector const* source,
                        LA_Vector* dest ) const = 0 ;
      
   protected: //--------------------------------------------------------

      virtual ~LA_Scatter( void ) ;
      LA_Scatter( MAC_Object* a_owner ) ;

   //-- Preconditions, Postconditions, Invariant
      
      virtual bool implementation_POST(
                            LA_Implementation const* result ) const ;

      virtual bool repatriated_items_POST(
                                size_t_vector const& result ) const ;

      virtual bool local_indices_POST(
                                size_t_vector const& result ) const ;
      
      virtual bool distribution_POST(
                     MAC_DistributedPartition const* result ) const ;
      
      virtual bool get_PRE( LA_Vector const* source,
                            LA_SeqVector const* dest ) const ;
      virtual bool get_POST( LA_Vector const* source,
                            LA_SeqVector const* dest ) const ;
      
      virtual bool set_PRE( LA_SeqVector const* source,
                            LA_Vector const* dest ) const ;
      virtual bool set_POST( LA_SeqVector const* source,
                            LA_Vector const* dest ) const ;
      
   private: //----------------------------------------------------------

      LA_Scatter( void ) ;
      LA_Scatter( LA_Scatter const& other ) ;
      LA_Scatter& operator=( LA_Scatter const& other ) ;

} ;

#endif
