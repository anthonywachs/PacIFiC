#ifndef PDE_SYSTEM_NUMBERING_HH
#define PDE_SYSTEM_NUMBERING_HH

#include <MAC_Object.hh>
#include <size_t_vector.hh>
#include <size_t_array2D.hh>

#include <vector>

class size_t_vector ;
class LA_Scatter ;
class LA_Vector ;
class FV_DiscreteField ;


class FV_SystemNumbering : public MAC_Object
{
   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      static FV_SystemNumbering* create( MAC_Object* a_owner,
	FV_DiscreteField const* the_field,
	size_t a_verbose_level = 0 ) ;
					        
      
   //-- Access to data

      FV_DiscreteField const* field( void ) const ;
      

   //-- Items of the global algebraic system

      size_t nb_global_unknowns( void ) const ;

      size_t nb_unknowns_on_current_process( void ) const ;


   //-- Scatters

      void define_scatter( LA_Vector const* vec ) ;

      bool scatter_is_defined( void ) const ;

      LA_Scatter const* scatter( void ) const ;


   //-- Input - Output

      void print( std::ostream& os, size_t indent_width ) const ;


   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      FV_SystemNumbering( void ) ;
     ~FV_SystemNumbering( void ) ;
      FV_SystemNumbering( FV_SystemNumbering const& other ) ;
      FV_SystemNumbering& operator=( FV_SystemNumbering const& other ) ;

      FV_SystemNumbering( MAC_Object* a_owner,
	FV_DiscreteField const* the_field,
	size_t a_verbose_level ) ;

   //-- Attributes

      FV_DiscreteField const* FIELD ;

      size_t_vector* IDX_LOCS ;
      size_t_vector* IDX_GLOBS ;

      bool OK_SCATTER ;
      LA_Scatter* SCATTER ;

      size_t NB_HANDLED_UNK ;
      size_t MAT_SIZE ;

      size_t VERB ;
} ;

#endif
