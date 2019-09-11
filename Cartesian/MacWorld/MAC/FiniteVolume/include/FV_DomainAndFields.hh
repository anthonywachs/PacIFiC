#ifndef FV_DOMAIN_AND_FIELDS_HH
#define FV_DOMAIN_AND_FIELDS_HH

//#include <FV_PostProcessingWriter_BuilderFactory.hh>
#include <FV_PostProcessingWriter.hh>
#include <MAC_Object.hh>
#include <list>
#include <stringVector.hh>
using std::list;

class MAC_Communicator ;
class MAC_Module ;
class MAC_ModuleExplorer ;
class MAC_Vector ;
class FV_DomainBuilder ;
class FV_Mesh ;
class FV_DiscreteField ;


/** @brief The Class FV_DomainAndFields.

Server to discretize the computational domain with a FV Finite Volume scheme.

@author A. Wachs - Particulate flow project 2010-2012 */

class FV_DomainAndFields : public MAC_Object
{
   public: //-----------------------------------------------------------------

   //-- Instance delivery and initialization

      static FV_DomainAndFields* create( MAC_Object* a_owner,
                                          MAC_ModuleExplorer const* exp ) ;

      // Create and return a new instance of the class and perform all
      // necessary tasks to reinterpret `self' as a subdomain of
      // a logical global domain in a processing distributed over
      // all processes in `com'.
      static FV_DomainAndFields* create( MAC_Object* a_owner,
                                          MAC_ModuleExplorer const* exp,
                                          MAC_Communicator const* com ) ;


   //-- Characteristics

      // name, uniquely determining `self'
      std::string const& name( void ) const ;


   //-- Modifiers(5)

      // Create a new discrete field called `name_of_new_field', identical
      // to the already existing discrete field called `model_name',
      // and add it to the list of discrete fields
      void duplicate_field( std::string const& model_name,
                            std::string const& name_of_new_field ) ;

      // Add a  new discrete field outside the Module Hierarchy
      // attached to `exp', and add it to the list of discrete fields 
      void append_field( FV_DiscreteField* newField ) ;


   //-- Distributed processing(10)

      // Does `self' represent a subdomain of a logical global domain
      // involved in a processing distributed over several processes ?
      bool is_distributed( void ) const ;


   //-- Problem Definition(100)

      // number of space dimensions
      size_t nb_space_dimensions( void ) const ;


   //-- Save for postprocessing(140)
      
      // Return (and create if not created yet) the post-processing writer
      FV_PostProcessingWriter* post_processing_writer( void ) const ;
 
      // Add a discrete field to the list of Paraview post-processed discrete
      // fields outside the Module Hierarchy attached to `exp'
      void append_post_processing_field( FV_DiscreteField* newField,
      	 std::string const& location, std::string const& save_name ) ;


   //-- Persistence

      virtual void save_state( MAC_ObjectWriter* writer ) const ;

      virtual void restore_state( MAC_ObjectReader* reader ) ;


   //-- Adaptation


   //-- Input - Output

      void print_grid( std::ostream& os, size_t indent_width ) const ;


   //-- Access

      FV_DiscreteField* discrete_field( std::string const& field_name ) 
      	const ;

      bool has_discrete_field( std::string const& field_name ) const ;
      
      FV_Mesh const* primary_grid( void ) const ;
      

   protected: //-------------------------------------------------------------

   private: //---------------------------------------------------------------

      FV_DomainAndFields( void ) ;
     ~FV_DomainAndFields( void ) ;
      FV_DomainAndFields( FV_DomainAndFields const& other ) ;
      FV_DomainAndFields& operator=( FV_DomainAndFields const& other ) ;

      FV_DomainAndFields( MAC_Object* a_owner,
                           MAC_ModuleExplorer const* exp,
                           MAC_Communicator const* com ) ;

   //-- Preconditions, Postconditions, Invariant

      virtual bool invariant( void ) const ;

   //-- Attributes

      MAC_ModuleExplorer const* EXP ;
      FV_DomainBuilder* BUILDER_FV ;
      MAC_Communicator const* COM ;
      list< FV_DiscreteField* > FVFIELDS ;      
      mutable list< FV_DiscreteField* > BUILDER_FVFIELDS ;
      mutable FV_PostProcessingWriter* BUILDER_WRITER ;
} ;

#endif
