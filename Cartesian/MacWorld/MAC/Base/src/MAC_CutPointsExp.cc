#include <MAC_CutPointsExp.hh>

#include <MAC_assertions.hh>
#include <MAC_Error.hh>
#include <MAC_Sequence.hh>
#include <doubleVector.hh>
#include <doubleArray2D.hh>

MAC_CutPointsExp const*
MAC_CutPointsExp::PROTOTYPE_MP = new MAC_CutPointsExp(
                         "middle_point", MAC_CutPointsExp::mid_point ) ;
MAC_CutPointsExp const*
MAC_CutPointsExp::PROTOTYPE_MPS = new MAC_CutPointsExp(
                       "middle_points", MAC_CutPointsExp::mid_points ) ;
MAC_CutPointsExp const*
MAC_CutPointsExp::PROTOTYPE_X = new MAC_CutPointsExp(
                              "x_cut_points", MAC_CutPointsExp::xcut ) ;
MAC_CutPointsExp const*
MAC_CutPointsExp::PROTOTYPE_Y = new MAC_CutPointsExp(
                              "y_cut_points", MAC_CutPointsExp::ycut ) ;
MAC_CutPointsExp const*
MAC_CutPointsExp::PROTOTYPE_Z = new MAC_CutPointsExp(
                              "z_cut_points", MAC_CutPointsExp::zcut ) ;

struct MAC_CutPointsExp_ERROR
{
   static void n0( std::string const& f_name, std::string const& op_name ) ;
} ;

//----------------------------------------------------------------------
MAC_CutPointsExp:: MAC_CutPointsExp( std::string const& a_name,
                                     IS_CutOp a_op ) 
//----------------------------------------------------------------------
   : MAC_Expression( a_name )
   , OP( a_op )
{
   MAC_LABEL( "MAC_CutPointsExp:: MAC_CutPointsExp" ) ;
}

//----------------------------------------------------------------------
MAC_CutPointsExp:: MAC_CutPointsExp( MAC_Object* a_owner,
                                     std::string const& a_name,
                                     MAC_Sequence const* argument_list,
                                     IS_CutOp a_op ) 
//----------------------------------------------------------------------
   : MAC_Expression( a_owner, a_name, argument_list )
   , OP( a_op )
{
   MAC_LABEL( "MAC_CutPointsExp:: MAC_CutPointsExp" ) ;
   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
MAC_CutPointsExp:: ~MAC_CutPointsExp( void ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_CutPointsExp:: ~MAC_CutPointsExp" ) ;
   MAC_CHECK_INV( invariant() ) ;

   if( this == PROTOTYPE_X ) PROTOTYPE_X = 0 ;
   if( this == PROTOTYPE_Y ) PROTOTYPE_Y = 0 ;
   if( this == PROTOTYPE_Z ) PROTOTYPE_Z = 0 ;
   if( this == PROTOTYPE_MP ) PROTOTYPE_MP = 0 ;
   if( this == PROTOTYPE_MPS ) PROTOTYPE_MPS = 0 ;
   
}

//----------------------------------------------------------------------
MAC_CutPointsExp*
MAC_CutPointsExp:: create_replica( MAC_Object* a_owner,
                                   MAC_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_CutPointsExp:: create_replica" ) ;
   MAC_CHECK( create_replica_PRE( a_owner, argument_list ) ) ;

   MAC_CutPointsExp* result =
           new MAC_CutPointsExp( a_owner, name(), argument_list, OP ) ;

   MAC_CHECK_POST( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
std::string const& 
MAC_CutPointsExp:: usage( void ) const
//----------------------------------------------------------------------
{
   static std::string result = "unspecified" ;
   if( OP == mid_point )
   {
      result = name()+"(DS,DV)" ;
   }
   else if( OP == mid_points )
   {
      result = name()+"([DV],DV)" ;
   }
   else if( OP == xcut )
   {
      result = name()+"(DV,DS[,DS])" ;
   }
   else if( OP == ycut )
   {
      result = name()+"(DS,DV[,DS] )" ;
   }
   else if( OP == zcut )
   {
      result = name()+"(DS,DS,DV)" ;
   }
   else
   {
      MAC_CutPointsExp_ERROR::n0( "usage", name() ) ;
   }
   
   return( result ) ;
}

//----------------------------------------------------------------------
bool
MAC_CutPointsExp:: valid_arguments( MAC_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_CutPointsExp:: valid_arguments" ) ;
   MAC_CHECK( valid_arguments_PRE( some_arguments ) ) ;

   bool result = false ;
   if( OP == mid_point )
   {
      result = ( some_arguments->count() == 2 ) ;
      if( result )
      {
         MAC_Data::Type k0 = extract_arg( some_arguments, 0 )->data_type() ;
         result &= ( k0 == Double ) ;
         MAC_Data::Type k1 = extract_arg( some_arguments, 1 )->data_type() ;
         result &= ( k1 == DoubleVector ) ;
      }
   }
   else if( OP == mid_points )
   {
      result = ( some_arguments->count() == 1 || some_arguments->count() == 2 ) ;
      if( result )
      {
         MAC_Data::Type k0 = extract_arg( some_arguments, 0 )->data_type() ;
         result &= ( k0 == DoubleVector ) ;
         if( some_arguments->count() == 2 )
         {
            MAC_Data::Type k1 = extract_arg( some_arguments, 1 )->data_type() ;
            result &= ( k1 == DoubleVector ) ;
         }
      }
   }
   else if( OP == xcut )
   {
      result = ( some_arguments->count() == 2 || some_arguments->count() == 3 ) ;
      if( result )
      {
         MAC_Data::Type k0 = extract_arg( some_arguments, 0 )->data_type() ;
         result &= ( k0 == DoubleVector ) ;
         MAC_Data::Type k1 = extract_arg( some_arguments, 1 )->data_type() ;
         result &= ( k1 == Double ) ;
         if( some_arguments->count() == 3 )
         {
            MAC_Data::Type k2 = extract_arg( some_arguments, 2 )->data_type() ;
            result &= ( k2 == Double ) ;
         }
      }
   }
   else if( OP == ycut )
   {
      result = ( some_arguments->count() == 2 || some_arguments->count() == 3 ) ;
      if( result )
      {
         MAC_Data::Type k0 = extract_arg( some_arguments, 0 )->data_type() ;
         result &= ( k0 == Double ) ;
         MAC_Data::Type k1 = extract_arg( some_arguments, 1 )->data_type() ;
         result &= ( k1 == DoubleVector ) ;
         if( some_arguments->count() == 3 )
         {
            MAC_Data::Type k2 = extract_arg( some_arguments, 2 )->data_type() ;
            result &= ( k2 == Double ) ;
         }
      }
   }
   else if( OP == zcut )
   {
      result = ( some_arguments->count() == 3 ) ;
      if( result )
      {
         MAC_Data::Type k0 = extract_arg( some_arguments, 0 )->data_type() ;
         result &= ( k0 == Double ) ;
         MAC_Data::Type k1 = extract_arg( some_arguments, 1 )->data_type() ;
         result &= ( k1 == Double ) ;
         MAC_Data::Type k2 = extract_arg( some_arguments, 2 )->data_type() ;
         result &= ( k2 == DoubleVector ) ;
      }
   }
   else
   {
      MAC_CutPointsExp_ERROR::n0( "valid_arguments", name() ) ;
   }

   return( result ) ;
}

//----------------------------------------------------------------------
MAC_Data::Type
MAC_CutPointsExp:: data_type( void ) const
//----------------------------------------------------------------------
{
   MAC_Data::Type result = Undefined ;
   if( OP == mid_point )
   {
      result = Double ;
   }
   else if( OP == mid_points )
   {
      result = DoubleVector ;
   }
   else if( OP == xcut ||  OP == ycut || OP == zcut )
   {
      result = DoubleArray2D ;
   }
   else
   {
      MAC_CutPointsExp_ERROR::n0( "data_type", name() ) ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
double
MAC_CutPointsExp:: to_double( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_CutPointsExp:: to_double" ) ;
   MAC_CHECK_PRE( to_double_PRE( ct ) ) ;

   double result = MAC::bad_double() ;

   if( OP == mid_point )
   {
      double const x = arg(0)->to_double( ct ) ;
      doubleVector const& verts_table = arg(1)->to_double_vector( ct ) ;
      result = m_pt( x, verts_table ) ;
   }
   else
   {
      MAC_CutPointsExp_ERROR::n0( "to_double", name() ) ;
   }
   
   return( result ) ;
}

//----------------------------------------------------------------------
doubleVector const&
MAC_CutPointsExp:: to_double_vector( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_CutPointsExp:: to_double_vector" ) ;
   MAC_CHECK_PRE( to_double_vector_PRE( ct ) ) ;

   static doubleVector result(0) ;

   if( OP == mid_points )
   {
      if( nb_arguments()==1 )
      {
         doubleVector const& verts_table = arg(0)->to_double_vector( ct ) ;
         build_coords( verts_table, result ) ;
      }
      else if( nb_arguments()==2 )
      {
         doubleVector const& pt_table = arg(0)->to_double_vector( ct ) ;
         doubleVector const& verts_table = arg(1)->to_double_vector( ct ) ;
         result.re_initialize( pt_table.size() ) ;
         for( size_t i=0 ; i<pt_table.size() ; ++i )
         {
            result(i) = m_pt( pt_table(i), verts_table ) ;
         }
      }
      else
      {
         MAC_CutPointsExp_ERROR::n0( "to_double_vector", name() ) ;
      }
   }
   else
   {
      MAC_CutPointsExp_ERROR::n0( "to_double_vector", name() ) ;
   }
   
   return( result ) ;
}

//----------------------------------------------------------------------
doubleArray2D const&
MAC_CutPointsExp:: to_double_array2D( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_CutPointsExp:: to_double_array2D" ) ;
   MAC_CHECK_PRE( to_double_array2D_PRE( ct ) ) ;

   static doubleArray2D result(0,0) ;
   
   size_t const nb_dims = nb_arguments() ;
   if( OP == xcut )
   {
      doubleVector const& verts_table = arg(0)->to_double_vector( ct ) ;
      doubleVector coords_table(0) ;
      build_coords( verts_table, coords_table ) ;
      result.re_initialize( coords_table.size(), nb_dims ) ;
      for( size_t i=0 ; i<coords_table.size() ; ++i )
      {
         result(i,0) = coords_table(i) ;
         result(i,1) = arg(1)->to_double( ct ) ;
         if( nb_dims==3 ) result(i,2) = arg(2)->to_double( ct ) ;
      }
   }
   else if( OP == ycut )
   {
      doubleVector const& verts_table = arg(1)->to_double_vector( ct ) ;
      doubleVector coords_table(0) ;
      build_coords( verts_table, coords_table ) ;
      result.re_initialize( coords_table.size(), nb_dims ) ;
      for( size_t i=0 ; i<coords_table.size() ; ++i )
      {
         result(i,0) = arg(0)->to_double( ct ) ;
         result(i,1) = coords_table(i) ;
         if( nb_dims==3 ) result(i,2) = arg(2)->to_double( ct ) ;
      }
      
   }
   else if( OP == zcut )
   {
      doubleVector const& verts_table = arg(2)->to_double_vector( ct ) ;
      doubleVector coords_table(0) ;
      build_coords( verts_table, coords_table ) ;
      result.re_initialize( coords_table.size(), nb_dims ) ;
      for( size_t i=0 ; i<coords_table.size() ; ++i )
      {
         result(i,0) = arg(0)->to_double( ct ) ;
         result(i,1) = arg(1)->to_double( ct ) ;
         result(i,2) = coords_table(i) ;
      }
   }
   else
   {
      MAC_CutPointsExp_ERROR::n0( "to_double_array2D", name() ) ;
   }
   
   return( result ) ;
}

//----------------------------------------------------------------------
void
MAC_CutPointsExp:: check_table( doubleVector const& verts_table ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_CutPointsExp:: check_table" ) ;

   if( verts_table.size() <= 1 )
   {
      raise_error( "cutline vector should have at least two elements" ) ;
   }

   for( size_t i=0 ; i<verts_table.size()-1 ; ++i )
   {
      if( verts_table(i)>=verts_table(i+1) )
      {
         raise_error( "cutline vector should be increasing" ) ;
      }
   }
}

//----------------------------------------------------------------------
double
MAC_CutPointsExp:: m_pt( double x,
                         doubleVector const& verts_table ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_CutPointsExp:: m_pt" ) ;

   check_table( verts_table ) ;

   double result = MAC::bad_double() ;
   for( size_t i=0 ; i<verts_table.size()-1 ; ++i )
   {
      MAC_CHECK( verts_table(i)<verts_table(i+1) ) ;
      result = 0.5*( verts_table(i)+verts_table(i+1) ) ;
      if( x<verts_table(i+1) ) break ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
void
MAC_CutPointsExp:: build_coords( doubleVector const& verts_table,
                                 doubleVector& coords_table ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_CutPointsExp:: build_coords" ) ;

   check_table( verts_table ) ;
   
   coords_table.re_initialize( verts_table.size()-1 ) ;
   for( size_t i=0 ; i<verts_table.size()-1 ; ++i )
   {
      coords_table(i) = 0.5*( verts_table(i)+verts_table(i+1) ) ;
   }
}

//internal--------------------------------------------------------------
void 
MAC_CutPointsExp_ERROR:: n0( std::string const& f_name,
                             std::string const& op_name )
//internal--------------------------------------------------------------
{
   std::string mesg ;
   mesg += "*** MAC_CutPointsExp::" + f_name +"\n" ;
   mesg += "    operation " + op_name + " not implemented." ;
   MAC_Error::object()->raise_internal( mesg ) ;
}
