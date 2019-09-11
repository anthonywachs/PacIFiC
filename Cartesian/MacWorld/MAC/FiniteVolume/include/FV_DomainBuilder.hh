#ifndef FV_DOMAIN_BUILDER_HH
#define FV_DOMAIN_BUILDER_HH

#include <MAC_Object.hh>
#include <size_t_array2D.hh>
#include <size_t_vector.hh>
#include <stringVector.hh>
#include <map>
#include <string>
#include <vector>
#include <list>
using std::list ;
using std::pair ;
using std::string ; 

class MAC_List ;
class MAC_ListIterator ;
class MAC_ModuleExplorer ;
class MAC_Vector ;
class MAC_VectorIterator ;
class FV_Mesh ;
class FV_DiscreteField ;
struct FV_SHIFT_TRIPLET ;


enum FV_AllowedBoundaryColors { 
	FV_BC_INTERIOR,
	FV_BC_PERIODIC,
	FV_BC_LEFT,
	FV_BC_RIGHT,
	FV_BC_BOTTOM,	
	FV_BC_TOP,	
	FV_BC_BOTTOM_LEFT,
	FV_BC_BOTTOM_RIGHT,
	FV_BC_TOP_LEFT,
	FV_BC_TOP_RIGHT,
	FV_BC_BEHIND,
	FV_BC_FRONT,
	FV_BC_BEHIND_LEFT,
	FV_BC_BEHIND_RIGHT,
	FV_BC_FRONT_LEFT,
	FV_BC_FRONT_RIGHT,
	FV_BC_BEHIND_BOTTOM,
	FV_BC_BEHIND_TOP,	
	FV_BC_FRONT_BOTTOM,
	FV_BC_FRONT_TOP,	
	FV_BC_BEHIND_BOTTOM_LEFT,	
	FV_BC_BEHIND_BOTTOM_RIGHT,
	FV_BC_BEHIND_TOP_LEFT,
	FV_BC_BEHIND_TOP_RIGHT,		
	FV_BC_FRONT_BOTTOM_LEFT,		
	FV_BC_FRONT_BOTTOM_RIGHT,		
	FV_BC_FRONT_TOP_LEFT,		
	FV_BC_FRONT_TOP_RIGHT } ;
	
enum FV_AllowedDOFStatus { 
	FV_DOF_ONPROC,
	FV_DOF_BUFFERZONE,	 
	FV_DOF_HALOZONE,	
	FV_DOF_PERIODIC_BUFFERZONE,	
	FV_DOF_PERIODIC_HALOZONE,
	FV_DOF_BUFFERZONE_PERIODIC_BUFFERZONE } ;
	
	
/** @brief The Class FV_DomainBuilder.

Build the FV discretized problem.

@author A. Wachs - Particulate flow project 2010-2012 */

class FV_DomainBuilder : public MAC_Object
{
   public: //-----------------------------------------------------------

   //-- Instance delivery and initialization

      static FV_DomainBuilder* object( string const& a_name ) ;

      static FV_DomainBuilder* create( MAC_Object* a_owner, 
                                        MAC_ModuleExplorer const* exp,
                                        string const& a_name  ) ;


   //-- Characteristics

      std::string const& name( void ) const ;


   //-- Associated Module Hierarchy

      bool has_explorer( string const& path_and_name ) const ;

      MAC_ModuleExplorer* create_explorer( 
                                  MAC_Object* a_owner,
                                  string const& path_and_name ) const ;


   //-- Persistence   

      virtual void save_state( MAC_ObjectWriter* writer ) const ;

      virtual void restore_state( MAC_ObjectReader* reader ) ;

      
   //-- Access

      FV_Mesh const* primary_grid( void ) const ;
      
      static size_t nb_space_dimensions( void ) ; 
      
      list< FV_DiscreteField* > const* list_primary_fields( void ) const;
      
      static size_t get_color_number( string const& color_name ) ;
      
      static string get_color_name( size_t const& color_id ) ;       

      static string get_status_name( size_t const& status_id ) ; 
      
      static bool does_primary_color_exist( string const& color_name ) ;
      
      static bool is_main_color( size_t const& colorID ) ;
      
      static list<size_t> const* get_sub_main_color_ids( 
      	size_t const& color_id ) ; 

      static FV_SHIFT_TRIPLET const* get_shift_MacTriplet( 
      	size_t const& colorID ) ; 
      
      stringVector const* get_colors_from_macro_color( 
      	string const& color_name ) const ;

      static pair<size_t,string> normal_direction_to_main_color( 
      	size_t const& colorID ) ;


   protected: //--------------------------------------------------------

   private: //----------------------------------------------------------

      FV_DomainBuilder( void ) ;
     ~FV_DomainBuilder( void ) ;
      FV_DomainBuilder( FV_DomainBuilder const& other ) ;
      FV_DomainBuilder& operator=( FV_DomainBuilder const& other ) ;
     
      FV_DomainBuilder( MAC_Object* a_owner,
                         MAC_ModuleExplorer const* exp,
                         string const& a_name ) ;

      void build_all( void ) ;

      void build_grid( void ) ;
      
      void build_discrete_fields( void ) ;
      
      void build_special_colors( void ) ;


   //-- Class Attributes & Methods

      static std::map< string, FV_DomainBuilder* >& INSTANCES( void ) ;

      void build_allowed_DOFcolors( void );

      void build_allowed_DOFstatus( void );
      
      
   //-- Attributes

      MAC_ModuleExplorer const* EXP ;
      std::string NAME ;
      int VERB ;
      static size_t DIM ;
      FV_Mesh* GRID ;
      list< FV_DiscreteField* > PRIMARY_FVFIELDS ;
      list< pair< string, stringVector > > FVRO_COLORS ;
      static list< pair< size_t, list<size_t> > >* subMainColorIds ;
      static stringVector* allowedDOFcolors;
      static stringVector* allowedDOFstatus;        
      static list< pair< size_t, FV_SHIFT_TRIPLET > >* 
      	mainColorShiftMacTriplets ;  
} ;

#endif
