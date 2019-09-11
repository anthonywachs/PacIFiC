#ifndef MAC_MEMBERSHIP_EXP_HH
#define MAC_MEMBERSHIP_EXP_HH

#include <MAC_Expression.hh>

/* Operator in_range tests on value belonging to range.

PUBLISHED
 */

class MAC_MembershipExp : public MAC_Expression
{
   public: //---------------------------------------------------------

   //-- Type
      
      virtual MAC_Data::Type data_type( void ) const ;
      
   //-- Value
      
      virtual bool to_bool( MAC_Context const* ct ) const ;

   protected: //-------------------------------------------------------
      
   private: //---------------------------------------------------------

      MAC_MembershipExp( void ) ; 
     ~MAC_MembershipExp( void ) ;
      MAC_MembershipExp( MAC_MembershipExp const& other ) ;
      MAC_MembershipExp& operator=( MAC_MembershipExp const& other ) ;

      enum MembExp{ in_range, in_box } ;

      MAC_MembershipExp( MAC_Object* a_owner,
                         MembExp exp_id,
                         std::string const& a_name,
                         MAC_Sequence const* argument_list ) ;
      
   //-- Plug in

      MAC_MembershipExp( MembExp exp_id, std::string const& a_name ) ;
      
      virtual MAC_MembershipExp* create_replica( 
                                   MAC_Object * a_owner,
                                   MAC_Sequence const* argument_list ) const ;
      
   //-- Characteristics

      virtual std::string const& usage( void ) const ;

      virtual bool valid_arguments( MAC_Sequence const* some_arguments ) const ;
      
   //-- Class attributes

      static MAC_MembershipExp const* PROTOTYPE_in_range ;
      static MAC_MembershipExp const* PROTOTYPE_in_box ;

   //-- Attributes

      MembExp const OP ;
      MAC_Data const* const ARG0 ;
      MAC_Data const* const ARG1 ;
      MAC_Data const* const ARG2 ;
} ;

#endif
