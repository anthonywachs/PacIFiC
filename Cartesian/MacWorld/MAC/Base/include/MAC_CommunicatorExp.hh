#ifndef MAC_COMMUNICATOREXP_HH
#define MAC_COMMUNICATOREXP_HH

#include <MAC_Expression.hh>

class MAC_Communicator ;

/*
MAC_Communicator related expressions : nb_ranks() and rank().

PUBLISHED
*/

class MAC_CommunicatorExp : public MAC_Expression
{
   public: //-------------------------------------------------------

   //-- Type
      
      virtual MAC_Data::Type data_type( void ) const ;
      
   //-- Value
      
      virtual int to_int( MAC_Context const* ct = 0 ) const ;
      
   protected: //-----------------------------------------------------
      
   private: //-------------------------------------------------------

      MAC_CommunicatorExp( void ) ;
     ~MAC_CommunicatorExp( void ) ;
      MAC_CommunicatorExp( MAC_CommunicatorExp const& other ) ;
      MAC_CommunicatorExp& operator=( MAC_CommunicatorExp const& other ) ;

      MAC_CommunicatorExp( MAC_Object* a_owner,
                           std::string const& a_name,
                           MAC_Sequence const* argument_list ) ;

   //-- Plug in

      MAC_CommunicatorExp( std::string const& a_name ) ;

      virtual MAC_CommunicatorExp* create_replica( 
                                   MAC_Object * a_owner,
                                   MAC_Sequence const* argument_list ) const ;
      
   //-- Characteristics

      virtual std::string const& usage( void ) const ;

      virtual bool valid_arguments(
                                 MAC_Sequence const* some_arguments ) const ;
      

   //-- Class attributes
            
      static MAC_CommunicatorExp const* PROTOTYPE_RANK ;
      static MAC_CommunicatorExp const* PROTOTYPE_NB_RANKS ;

   //-- Attributes
      
      MAC_Communicator const* const COM ;    
} ;

#endif
