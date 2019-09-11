#ifndef _FV_POSTPROCESSINGWRITER_HH_
#define _FV_POSTPROCESSINGWRITER_HH_

#include <MAC_Object.hh>
#include <FV_DiscreteField.hh>
#include <string>
#include <list>
using std::list ;
using std::string ;

class MAC_ObjectRegister ;
class MAC_Communicator ;
class MAC_Module ;
class MAC_ModuleExplorer ;
class FV_Mesh ;
class FV_DiscreteField ;
class FV_TimeIterator ;


/** @brief The Class PostProcessingWriter.

Write mesh and fields data (Paraview vtk format or matlab format)

@author A. Hammouti - SIMULTI project  2013-2014 */
//=============================================================================
class FV_PostProcessingWriter : public MAC_Object
{
   public: //-----------------------------------------------------------
    
    static FV_PostProcessingWriter* make( MAC_Object* a_owner,
		 std::string const& a_name ,			   
 		 MAC_ModuleExplorer const* exp,
	         MAC_Communicator const* com,
		 list< FV_DiscreteField const* > a_fields,
		 FV_Mesh const* a_primary_mesh,
		 bool a_binary );
    
    virtual void write_cycle( FV_TimeIterator const* t_it,
      			size_t cycle_number ) = 0 ;
    
   //-- Data clearing

     virtual void clearResultFiles() = 0 ;
      
     virtual size_t getPreviousCycleNumber() = 0 ;
      
     virtual void readTimeFile( FV_TimeIterator const* t_it,
      			size_t& cycle_number ) = 0 ; 

    protected: //----------------------------------------------------------
    
       virtual ~FV_PostProcessingWriter( void );
  
      FV_PostProcessingWriter( MAC_Object* a_owner,
		MAC_ModuleExplorer const* exp ,
                MAC_Communicator const* com,
       		list< FV_DiscreteField const* > a_fields,
       		FV_Mesh const* a_primary_mesh,
        	bool a_binary );

       FV_PostProcessingWriter( std::string const& a_name ) ;

       FV_PostProcessingWriter( MAC_Object* a_owner ) ;

       virtual FV_PostProcessingWriter* create_replica( MAC_Object* a_owner,
	        MAC_ModuleExplorer const* exp ,
                MAC_Communicator const* com,
       		list< FV_DiscreteField const* > a_fields,
       		FV_Mesh const* a_primary_mesh,
        	bool a_binary ) const = 0 ;

      bool is_a_prototype( void ) const ;
      
   //-- Preconditions, Postconditions, Invariant

      /** @name Preconditions, Postconditions, Invariant */
      //@{      
      /** @brief Test invariance */
      virtual bool invariant( void ) const ;
          
      /** @brief Postcondition of create
      @param result the created object
      @param a_owner the MAC-based object */      
      virtual bool create_replica_POST( 
      	 FV_PostProcessingWriter const* result,
      	 MAC_Object* a_owner ) const ;
      //@}          	             

      /** @brief Return a pointer to the object register */
      static MAC_ObjectRegister* plugins_map( void ) ;      


   //-- Attributes
      
      bool const IS_PROTO ;
      
 private://----------------------------------------------------------
   
       //-- Constructors & Destructor

       /** @name Constructors & Destructor */
       //@{
       /** @brief Constructor without argument */      
       FV_PostProcessingWriter( void ) ;

       /** @brief Copy constructor */       
       FV_PostProcessingWriter( FV_PostProcessingWriter const& other ) ;

       /** @brief Operator == 
       @param other the right hand side */   
       FV_PostProcessingWriter& operator=( 
	 FV_PostProcessingWriter const& other ) ;
       //@}
	
} ;

#endif
