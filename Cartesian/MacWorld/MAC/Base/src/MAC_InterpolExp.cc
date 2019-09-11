#include <MAC_InterpolExp.hh>

#include <MAC_assertions.hh>
#include <MAC_Error.hh>
#include <MAC_Sequence.hh>
#include <MAC.hh>

#include <iostream>
#include <fstream>
#include <sstream>

MAC_InterpolExp const* 
MAC_InterpolExp::PROTOTYPE_LIN_1D = new MAC_InterpolExp(
                                                "interpol", lin_inter_1D ) ;

struct MAC_InterpolExp_ERROR
{
   static void n0( std::string const& f_name, std::string const& op_name ) ;
   static void n1( std::string const& f_name,
                   std::string const& filename,
                   size_t line_nb,
                   std::string const& line ) ;
   static void n2( std::string const& f_name, std::string const& filename ) ;
   static void n3( std::string const& f_name,
                   std::string const& filename,
                   size_t line_nb,
                   double xn,
                   double xn1 ) ;
   static void n4( std::string const& f_name,
                   std::string const& filename ) ;
} ;

//----------------------------------------------------------------------
MAC_InterpolExp:: MAC_InterpolExp( std::string const& a_name,
                                 MAC_InterpolOp a_op ) 
//----------------------------------------------------------------------
   : MAC_Expression( a_name )
   , OP( a_op )
   , FROM_FILE( false )
   , X1_IS_SET( false )
   , X1(0)
   , F_IS_SET( false )
   , FX1(0)
   , CHECK( false )
{
   MAC_LABEL( "MAC_InterpolExp:: MAC_InterpolExp" ) ;
}

//----------------------------------------------------------------------
MAC_Expression*
MAC_InterpolExp:: create_replica( MAC_Object* a_owner,
                                 MAC_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_InterpolExp:: create_replica" ) ;
   MAC_CHECK( create_replica_PRE( a_owner, argument_list ) ) ;

   MAC_InterpolExp* result = new MAC_InterpolExp( a_owner,
                                                name(), 
                                                argument_list,
                                                OP ) ;

   MAC_CHECK( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
MAC_InterpolExp:: MAC_InterpolExp( MAC_Object* a_owner,
                                   std::string const& a_name,
                                   MAC_Sequence const* argument_list,
                                   MAC_InterpolOp a_op ) 
//----------------------------------------------------------------------
   : MAC_Expression( a_owner, a_name, argument_list )
   , OP( a_op )
   , FROM_FILE( false )
   , X1_IS_SET( false )
   , X1(0)
   , F_IS_SET( false )
   , FX1(0)
   , CHECK( false )
{
   MAC_LABEL( "MAC_InterpolExp:: MAC_InterpolExp" ) ;

   switch( OP )
   {
      case lin_inter_1D :
         FROM_FILE = ( arg(0)->data_type() == String ) ;
         break ;
      default:
         MAC_InterpolExp_ERROR::n0( "MAC_InterpolExp::", name() ) ;
         break ;
   }
   
   MAC_CHECK_INV( invariant() ) ;
}

//----------------------------------------------------------------------
MAC_InterpolExp:: ~MAC_InterpolExp( void ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_InterpolExp:: ~MAC_InterpolExp" ) ;
   MAC_CHECK_INV( invariant() ) ;

   if( PROTOTYPE_LIN_1D == this )
   {
      PROTOTYPE_LIN_1D = 0 ;
   }
}

//----------------------------------------------------------------------
MAC_Data::Type
MAC_InterpolExp:: data_type( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_InterpolExp:: data_type" ) ;

   MAC_Data::Type result = MAC_Data::Undefined ;
   switch( OP )
   {
      case lin_inter_1D :
         result = MAC_Data::Double ;
         break ;
      default:
         MAC_InterpolExp_ERROR::n0( "data_type", name() ) ;
         break ;
   }     
   return( result ) ;
}

//----------------------------------------------------------------------
std::string const& 
MAC_InterpolExp:: usage( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_InterpolExp:: usage" ) ;

   static std::string result = "undefined" ;
   switch( OP )
   {
      case lin_inter_1D :
         result = name() + "(DV,DV,DV) or "+ name() + "(IS,DV) " ;
         break ;
      default:
         MAC_InterpolExp_ERROR::n0( "usage", name() ) ;
         break ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
bool
MAC_InterpolExp:: valid_arguments( MAC_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_InterpolExp:: valid_arguments" ) ;
   MAC_CHECK( valid_arguments_PRE( some_arguments ) ) ;

   bool result = false ;
   switch( OP )
   {
      case lin_inter_1D :
         result = ( some_arguments->count() > 0 ) ;
         if( result )
         {
            MAC_Data::Type k0 =
                             extract_arg( some_arguments, 0 )->data_type() ;
            if( k0 == String )
            {
               result = ( some_arguments->count() == 2 ) ;
               if( result )
               {
                  MAC_Data::Type k1 =
                             extract_arg( some_arguments, 1 )->data_type() ;
                  result = ( k1 == Double ) ;
               }
            }
            else if( k0 == DoubleVector )
            {
               result = ( some_arguments->count() == 3 ) ;
               if( result )
               {
                  MAC_Data::Type k1 =
                             extract_arg( some_arguments, 1 )->data_type() ;
                  result = ( k1 == DoubleVector ) ;
                  MAC_Data::Type k2 =
                             extract_arg( some_arguments, 2 )->data_type() ;
                  result &= ( k2 == Double ) ;
               }
            }
            else
            {
               result = false ;
            }
         }
         break ;
      default:
         MAC_InterpolExp_ERROR::n0( "valid_arguments", name() ) ;
         break ;
   }     
   return( result ) ;
}

//----------------------------------------------------------------------
double
MAC_InterpolExp:: to_double( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_InterpolExp:: to_double" ) ;
   MAC_CHECK_PRE( to_double_PRE(ct) ) ;

   double result = MAC::bad_double() ;

   switch( OP )
   {
      case lin_inter_1D :
         {
            double x = MAC::bad_double() ;
            if( FROM_FILE )
            {
               if( !X1_IS_SET )
               {
                  read_tables_1D( arg(0)->to_string( ct ), X1, FX1 ) ;
                  X1_IS_SET = true ;
                  F_IS_SET = true ;
                  CHECK = true ;
               }
               x = arg(1)->to_double( ct ) ;
            }
            else
            {
               if( !X1_IS_SET )
               {
                  X1 = arg(0)->to_double_vector(ct) ;
               }
               if( !F_IS_SET )
               {
                  FX1 = arg(1)->to_double_vector(ct) ;
               }
               if( !X1_IS_SET || !F_IS_SET )
               {
                  check_tables_1D( X1, FX1 ) ;
               }
               X1_IS_SET = true ;
               F_IS_SET = true ;
               CHECK = true ;
               x = arg(2)->to_double( ct ) ;
            }
            result = linear_interpol_1D( X1, FX1, x ) ;
         }
         break ;
      default:
         MAC_InterpolExp_ERROR::n0( "to_double", name() ) ;
         break ;
   }
   return( result ) ;
}

//----------------------------------------------------------------------
void
MAC_InterpolExp:: read_tables_1D( std::string const& filename,
                                  doubleVector& X_table,
                                  doubleVector& FX_table ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_InterpolExp:: read_tables_1D" ) ;
   
   std::ifstream in( filename.c_str() ) ;
   if( !in ) MAC_InterpolExp_ERROR::n2( name(), filename ) ;

   size_t n = 0 ;
   double x, y ;
   std::string line ;
   size_t nb = 0 ;
   
   while( getline( in, line ) )
   {
      std::istringstream sin( line ) ;
      
      while( sin >> x )
      {
         if( ! ( sin >> y ) )
         {
            MAC_InterpolExp_ERROR::n1( name(), filename, nb+1, line ) ;
         }
         X_table.append( x ) ;
         FX_table.append( y ) ;
         if( n!=0 &&  X_table(n)<=X_table(n-1) )
         {
            MAC_InterpolExp_ERROR::n3(
                  name(), filename, nb+1, X_table(n), X_table(n-1) ) ;
         }
         n++ ;
      }
      if( !sin.eof() )
      {
         MAC_InterpolExp_ERROR::n1( name(), filename, nb+1, line ) ;
      }
      nb++ ;
   }
   if( n == 0 ) MAC_InterpolExp_ERROR::n4( name(), filename ) ;
}

//----------------------------------------------------------------------
double
MAC_InterpolExp:: linear_interpol_1D( doubleVector const& X_table,
                                      doubleVector const& FX_table,
                                      double x ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_InterpolExp:: linear_interpol_1D" ) ;
   MAC_CHECK( X_table.size() == FX_table.size() ) ;

   double result = FX_table(0) ;
   if( x>X_table(0) )
   {
      size_t i=1 ;
      for( ; i<X_table.size() ; i++ )
      {
         if( x<X_table(i) )
         {
            break ;
         }
      }
      if( i<X_table.size() )
      {
         double const dy = FX_table(i)-FX_table(i-1) ;
         double const dx = X_table(i)-X_table(i-1) ;
         result = FX_table(i-1) + ( x-X_table(i-1) )*dy/dx ;
      }
      else
      {
         result = FX_table(i-1) ;
      }
   }
   return( result ) ;
}

//----------------------------------------------------------------------
void
MAC_InterpolExp:: check_tables_1D( doubleVector const& X_table,
                                   doubleVector const& FX_table ) const
//----------------------------------------------------------------------
{
   if( X_table.size() != FX_table.size() )
   {
      raise_error( "Same size expected for abscissea and ordinates tables" ) ;
   }
   for( size_t i=1 ; i<X_table.size() ; ++i )
   {
      if( X_table(i-1)>=X_table(i) )
      {
         raise_error( "Increasing first table expected" ) ;
      }
   }
}

//internal--------------------------------------------------------------
void 
MAC_InterpolExp_ERROR:: n0( std::string const& f_name,
                            std::string const& op_name )
//internal--------------------------------------------------------------
{
   std::string mesg ;
   mesg += "*** MAC_InterpolExp::"+f_name+"\n" ;
   mesg += "    operation "+op_name+" not implemented." ;
   MAC_Error::object()->raise_internal( mesg ) ;
}

//internal--------------------------------------------------------------
void 
MAC_InterpolExp_ERROR:: n1( std::string const& f_name,
                            std::string const& filename,
                            size_t line_nb,
                            std::string const& line )
//internal--------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << "*** MAC_InterpolExp::"+f_name+" error\n" ;
   mesg << "    Syntax error reading file : \""+filename+"\"\n" ;
   mesg << "    Error at line " ;
   mesg << line_nb << " : \"" ;
   mesg << line << "\"" << std::endl ;
   
   
   MAC_Error::object()->raise_plain( mesg.str() ) ;
}

//internal--------------------------------------------------------------
void 
MAC_InterpolExp_ERROR:: n2( std::string const& f_name,
                            std::string const& filename )
//internal--------------------------------------------------------------
{
   std::string mesg ;
   mesg += "*** MAC_InterpolExp::"+f_name+" error\n" ;
   mesg += "    Unable to open \""+filename+"\"." ;
   MAC_Error::object()->raise_plain( mesg ) ;
}

//internal--------------------------------------------------------------
void 
MAC_InterpolExp_ERROR:: n3( std::string const& f_name,
                            std::string const& filename,
                            size_t line_nb,
                            double xn,
                            double xn1 )
//internal--------------------------------------------------------------
{
   std::ostringstream mesg ;
   
   mesg << "*** MAC_InterpolExp::"+f_name+" error\n" ;
   mesg << "    Bad file \""+filename+"\"\n" ;
   mesg << "    At line ";
   mesg << line_nb ;
   mesg << " increasing first table expected but: "  ;
   MAC::print_double( mesg, xn ) ;
   mesg << " < " ;
   MAC::print_double( mesg, xn1 ) ;
   MAC_Error::object()->raise_plain( mesg.str() ) ;
}

//internal--------------------------------------------------------------
void 
MAC_InterpolExp_ERROR:: n4( std::string const& f_name,
                            std::string const& filename )
//internal--------------------------------------------------------------
{
   std::string mesg ;
   mesg += "*** MAC_InterpolExp::"+f_name+" error\n" ;
   mesg += "    syntax error reading file \""+filename+"\"" ;
   mesg += " : empty file." ;
    
   MAC_Error::object()->raise_plain( mesg ) ;
}
