#include <MAC_MembershipExp.hh>

#include <MAC_assertions.hh>
#include <MAC_Error.hh>
#include <MAC_Sequence.hh>
#include <doubleVector.hh>
#include <intVector.hh>

MAC_MembershipExp const*  
MAC_MembershipExp::PROTOTYPE_in_range = 
                             new MAC_MembershipExp( in_range, "in_range" ) ;

MAC_MembershipExp const*  
MAC_MembershipExp::PROTOTYPE_in_box = 
                             new MAC_MembershipExp( in_box, "in_box" ) ;

//----------------------------------------------------------------------
MAC_MembershipExp:: MAC_MembershipExp( MembExp exp_id, 
                                       std::string const& a_name  ) 
//----------------------------------------------------------------------
   : MAC_Expression( a_name )
   , OP( exp_id )
   , ARG0( 0 )
   , ARG1( 0 )
   , ARG2( 0 )
{
   MAC_LABEL( "MAC_MembershipExp:: MAC_MembershipExp" ) ;
}

//----------------------------------------------------------------------
MAC_MembershipExp*
MAC_MembershipExp:: create_replica( MAC_Object* a_owner,
                                    MAC_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_MembershipExp:: create_replica" ) ;
   MAC_CHECK( create_replica_PRE( a_owner, argument_list ) ) ;

   MAC_MembershipExp* result = new MAC_MembershipExp( a_owner, 
                                                      OP,
                                                      name(),
                                                      argument_list ) ;

   MAC_CHECK( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
MAC_MembershipExp:: MAC_MembershipExp( MAC_Object* a_owner,
                                       MembExp exp_id,
                                       std::string const& a_name,
                                       MAC_Sequence const* argument_list ) 
//----------------------------------------------------------------------
   : MAC_Expression( a_owner, a_name, argument_list )
   , OP( exp_id )
   , ARG0( arg(0) )
   , ARG1( arg(1) )
   , ARG2( exp_id == in_box ? arg(2) : 0 )
{
   MAC_LABEL( "MAC_MembershipExp:: MAC_MembershipExp" ) ;
   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
MAC_MembershipExp:: ~MAC_MembershipExp( void ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_MembershipExp:: ~MAC_MembershipExp" ) ;
   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
std::string const& 
MAC_MembershipExp:: usage( void ) const
//----------------------------------------------------------------------
{
   static std::string result ;
   switch( OP )
   {
      case in_range :
         result =
           "in_range(DS|IS,DV|IV)" ;
	 break ;
      case in_box :
         result = 
           "in_box(DV,DV,DV)" ;
	 break ;
   }
   return result ;
}

//----------------------------------------------------------------------
bool
MAC_MembershipExp:: valid_arguments( MAC_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_MembershipExp:: valid_arguments" ) ;
   MAC_CHECK( valid_arguments_PRE( some_arguments ) ) ;
   
   bool result = false ;
   switch( OP )
   {
      case in_range :
         result = some_arguments->count()==2 ;
	 if( result )
	 {
	    Type t0 = extract_arg( some_arguments, 0 )->data_type() ;
	    Type t1 = extract_arg( some_arguments, 1 )->data_type() ;
	    result = result && ( ( t0==Double && t1==DoubleVector ) || 
                                 ( t0==Int    && t1==IntVector    ) ) ;
	 }
	 break ;
      case in_box :
         result = ( some_arguments->count() == 3 ) ;
         if( result )
         {
            Type t0 = extract_arg( some_arguments, 0 )->data_type() ;
            Type t1 = extract_arg( some_arguments, 1 )->data_type() ;
            Type t2 = extract_arg( some_arguments, 2 )->data_type() ;
            result = result && ( t0 == DoubleVector && 
                                 t1 == DoubleVector && 
                                 t2 == DoubleVector ) ;
         }
	 break ;
      default : result = false ;
   }
   return result ;
}

//----------------------------------------------------------------------
MAC_Data::Type
MAC_MembershipExp:: data_type( void ) const
//----------------------------------------------------------------------
{
   return MAC_Data::Bool ;
}

//----------------------------------------------------------------------
bool
MAC_MembershipExp:: to_bool( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_MembershipExp:: to_bool" ) ;
   MAC_CHECK_PRE( to_bool_PRE( ct ) ) ;

   bool result = false ;
   switch( OP )
   {
      case in_range :
         double x ;
         double x_min ;
         double x_max ;
         if( ARG0->data_type()==Double )
         {
            x = ARG0->to_double( ct ) ;
            if( ARG1->to_double_vector( ct ).size()!=2 )
               raise_error( "second argument: two components expected" ) ;
            x_min = ARG1->to_double_vector( ct )( 0 ) ;
            x_max = ARG1->to_double_vector( ct )( 1 ) ;
         }
         else
         {
            x = ARG0->to_int( ct ) ;
            if( ARG1->to_int_vector( ct ).size()!=2 )
               raise_error( "second argument: two components expected" ) ;
            x_min = ARG1->to_int_vector( ct )(0) ;
            x_max = ARG1->to_int_vector( ct )(1) ;
         }
         if( x_min >= x_max )
            raise_error( "second argument: increasing values expected" ) ;
         result = ( x_min<=x  && x<=x_max ) ;
         break ;
      case in_box :
         doubleVector const& xx     = ARG0->to_double_vector( ct ) ;
         doubleVector const& xx_min = ARG1->to_double_vector( ct ) ;
         doubleVector const& xx_max = ARG2->to_double_vector( ct ) ;
         if( xx_min.size() != xx.size() )
            raise_error( "first and second arguments: incompatible sizes" ) ;
         if( xx_max.size() != xx.size() )
            raise_error( "first and third arguments: incompatible sizes" ) ;
         result = true ;
         for( size_t i=0 ; i<xx.size() ; ++i )
         {
            if( xx_min(i) >= xx_max(i) )
               raise_error(
                  "components of the third argument should be greater"
                  "\n    than those of the second one" ) ;
            result = result && ( xx_min(i) <= xx(i)     ) ;
            result = result && ( xx(i)     <= xx_max(i) ) ;
         }
         break ;
   }
   return result ;
}
