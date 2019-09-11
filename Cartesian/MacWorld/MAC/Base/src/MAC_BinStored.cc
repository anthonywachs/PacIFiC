#include <MAC_BinStored.hh>

#include <MAC_assertions.hh>
#include <MAC_Double.hh>
#include <MAC_BoolVector.hh>
#include <MAC_DoubleVector.hh>
#include <MAC_DoubleArray2D.hh>
#include <MAC_DoubleArray3D.hh>
#include <MAC_IntVector.hh>
#include <MAC_IntArray2D.hh>
#include <MAC_IntArray3D.hh>
#include <MAC_Int.hh>
#include <MAC_Error.hh>
#include <MAC_List.hh>
#include <MAC_Module.hh>
#include <MAC_Vector.hh>
#include <MAC_Root.hh>
#include <MAC_Sequence.hh>
#include <MAC_String.hh>
#include <MAC_System.hh>

#include <boolVector.hh>
#include <doubleVector.hh>
#include <doubleArray2D.hh>
#include <doubleArray3D.hh>
#include <intArray2D.hh>
#include <intArray3D.hh>
#include <intVector.hh>

#include <iostream>
#include <fstream>
#include <sys/types.h>
#include <sys/stat.h>

//internal---------------------------------------------------------------
class MAC_ThisFileDirExp : public MAC_Expression
{
   public:
      MAC_ThisFileDirExp( MAC_Object* a_owner ) ;
     ~MAC_ThisFileDirExp( void ) ; 
      virtual Type data_type( void ) const ;
      virtual std::string const& to_string( MAC_Context const* ct = 0 ) const ;
      virtual MAC_ThisFileDirExp* create_replica( 
                                   MAC_Object * a_owner,
                                   MAC_Sequence const* argument_list ) const ;
      
      virtual std::string const& usage( void ) const ;

      virtual bool valid_arguments( MAC_Sequence const* some_arguments ) const ;
} ;
//internal---------------------------------------------------------------


//----------------------------------------------------------------------
MAC_BinStored const* MAC_BinStored::static_OpComponent = new MAC_BinStored() ;
const size_t MAC_BinStored::bad_record = (size_t)~0 ;
size_t MAC_BinStored::last_record = bad_record ;
std::string MAC_BinStored::last_file = "" ;
struct magic 
{     char magic_char[256]  ;
      int magic_int ;
      size_t magic_size_t ;
      double magic_double ;
      bool magic_bool ;
} ;
static struct magic magic_stamp = { "MAC binary output format rev 1.1",
                                    -2,
                                    (size_t)~0,
                                    1./3.0,
                                    true } ;

size_t MAC_BinStored::preamble_length = sizeof( magic ) ;
//----------------------------------------------------------------------


//----------------------------------------------------------------------
MAC_BinStored:: MAC_BinStored( void ) 
//----------------------------------------------------------------------
      : MAC_TransferExp( "binary" ),
        my_data( 0 )
{
   MAC_LABEL( "MAC_BinStored:: MAC_BinStored" ) ;
}




//----------------------------------------------------------------------
MAC_BinStored:: MAC_BinStored( MAC_Object* a_owner,
                             MAC_Sequence const* argument_list ) 
//----------------------------------------------------------------------
      : MAC_TransferExp( a_owner, "binary", argument_list ),
        my_data( 0 )
{
   MAC_LABEL( "MAC_BinStored:: MAC_BinStored" ) ;
   MAC_CHECK_PRE( argument_list!=0 ) ;
   MAC_CHECK_INV( invariant() ) ;
}




//----------------------------------------------------------------------
MAC_BinStored:: ~MAC_BinStored( void ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_BinStored:: ~MAC_BinStored" ) ;
   MAC_CHECK_INV( invariant() ) ;
   if( my_data!=0 )  my_data->destroy() ;
}




//----------------------------------------------------------------------
MAC_BinStored*
MAC_BinStored:: create_replica( MAC_Object* a_owner,
                             MAC_Sequence const* argument_list ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_BinStored:: create_replica" ) ;
   MAC_CHECK_PRE( create_replica_PRE( a_owner, argument_list ) ) ;

   MAC_BinStored* result = new MAC_BinStored( a_owner, argument_list ) ;

   MAC_CHECK_POST( create_replica_POST( result, a_owner, argument_list ) ) ;
   return( result ) ;
}




//----------------------------------------------------------------------
std::string const& 
MAC_BinStored:: usage( void ) const
//----------------------------------------------------------------------
{
   static std::string result = "binary(SS,SS,IS)" ;
   return result ;
}




//----------------------------------------------------------------------
bool
MAC_BinStored:: valid_arguments( MAC_Sequence const* some_arguments ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_BinStored:: valid_arguments" ) ;
   MAC_CHECK( valid_arguments_PRE( some_arguments ) ) ;
   
   bool result = some_arguments->count()==3 &&
      extract_arg(some_arguments,0)->data_type() == String &&
      extract_arg(some_arguments,1)->data_type() == String &&
      extract_arg(some_arguments,2)->data_type() == Int ;
   return result ;
}




//----------------------------------------------------------------------
bool
MAC_BinStored:: is_type_supported( MAC_Data::Type a_type )
//----------------------------------------------------------------------
{
   bool result =
      a_type==Double ||
      a_type==DoubleVector ||
      a_type==DoubleArray2D ||
      a_type==DoubleArray3D ||
      a_type==IntArray2D ||
      a_type==IntArray3D ||
      a_type==IntVector ||
      a_type==BoolVector ;
   MAC_CHECK_POST( EQUIVALENT( result, ( a_type==Double ||
                                         a_type==DoubleVector ||
                                         a_type==DoubleArray2D ||
                                         a_type==DoubleArray3D ||
                                         a_type==IntArray2D ||
                                         a_type==IntArray3D ||
                                         a_type==IntVector ||
                                         a_type==BoolVector ) ) ) ;
   return result ;
}




//----------------------------------------------------------------------
MAC_Data::Type
MAC_BinStored:: data_type( void ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_BinStored:: kind" ) ;
   std::string tmp = arg(0)->to_string() ;
   Type result = Undefined ;
   
   if( tmp=="Double" )
   {
      result = Double ;
   }
   else if( tmp=="DoubleVector" )
   {
      result = DoubleVector ;
   }
   else if( tmp=="DoubleArray2D" )
   {
      result = DoubleArray2D ;
   }
   else if( tmp=="DoubleArray3D" )
   {
      result = DoubleArray3D ;
   }
   else if( tmp=="IntArray2D" )
   {
      result = IntArray2D ;
   }
   else if( tmp=="IntArray3D" )
   {
      result = IntArray3D ;
   }
   else if( tmp=="IntVector" )
   {
      result = IntVector ;
   }
   else if( tmp=="BoolVector" )
   {
      result = BoolVector ;
   }
   else
   {
      MAC_Error::object()->raise_plain( "Unsupported data type " + tmp ) ;
   }
   
   return result ;
}




//----------------------------------------------------------------------
MAC_Data const*
MAC_BinStored:: data( MAC_Context const* ct ) const
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_BinStored:: data" ) ;
   MAC_CHECK( data_PRE( ct ) ) ;
   
   if( my_data==0 )
   {
      std::string const& in_file = arg(1)->to_string( ct ) ;
      if( !is_valid_binary_file( in_file ) )
      {
         MAC_Error::object()->raise_plain(
            " File " + in_file
            + " is not a valid binary file on this system " ) ;
      }
      
      size_t record_number = arg(2)->to_int( ct ) ;
      my_data = restore_from_binary( 0, in_file, record_number ) ;
      
      MAC_ASSERT( my_data->data_type() == data_type() ) ;
   }

   MAC_CHECK_POST( data_POST( my_data, ct ) ) ;
   return my_data ;
}




//----------------------------------------------------------------------
bool 
MAC_BinStored:: is_valid_binary_file( std::string const& file_name ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_BinStored:: is_valid_binary_file" ) ;
   
   std::ifstream stream( file_name.c_str(), std::ios::binary | std::ios::in ) ;
   // Check for no-empty file
   bool result = stream && !stream.eof() ;
   // Check for magic stamp
   if( result )
   {
      struct magic tmp ;
      
      stream.read( (char*)&tmp, sizeof(magic) ) ;
      result = result &&
         ( stream.gcount()==sizeof(magic) ) ;
      std::string str = tmp.magic_char ;
      result = result &&
         str == magic_stamp.magic_char &&
         tmp.magic_int == magic_stamp.magic_int &&
         tmp.magic_size_t == magic_stamp.magic_size_t &&
         tmp.magic_double == magic_stamp.magic_double &&
         tmp.magic_bool == magic_stamp.magic_bool ;
      stream.close() ;
   }
   return result ;
}




//----------------------------------------------------------------------
size_t 
MAC_BinStored:: last_record_number( std::string const& file_name ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_BinStored:: last_record_number" ) ;
   MAC_CHECK( is_valid_binary_file( file_name ) ) ;
   
   size_t result = bad_record ;
   // CACHING
   if( file_name==last_file )
   {
      result = last_record ;
   }
   else
   {
      std::ifstream stream( file_name.c_str(), 
      		std::ios::binary | std::ios::in ) ;
      // Ignore preamble
      stream.ignore( preamble_length ) ;

      BinaryRecord record ;
      size_t record_size = sizeof( BinaryRecord ) ;
      while( !stream.eof() )
      {
         stream.read( (char*)&record, record_size ) ;
      
         if( (size_t)stream.gcount()!=record_size )
         {
            break ;
         }
         result = record.number ;
         stream.ignore( record.length ) ;
      }
      stream.close() ;
   }
   return result ;
}




//----------------------------------------------------------------------
MAC_Data const*
MAC_BinStored:: restore_from_binary( MAC_Object * a_owner,
                                     std::string const& file_name,
                                     size_t record_number ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_BinStored:: restore_from_binary" ) ;
   MAC_CHECK_PRE( record_number<=last_record_number( file_name ) ) ;

   static std::string oldname ;
   static size_t old_record = 0-1 ;
   static std::ifstream stream ;
   // Garder une reference sur un fichier externe est certe plus rapide mais 
   // peut generer des erreurs innattendues si le fichier disparait en cours 
   // de route.
   static bool keep_reference = false ;
   
   static bool first = true ;
   static off_t last_modif ;
   static struct stat status ;
   if( stat(file_name.c_str(),&status)!=0 )
   {
      MAC_Error::object()->raise_plain(
         "Unable to inquire for binary file "+file_name ) ;
   }
   off_t last = status.st_size ;
   MAC_Data const* result = 0 ;
   if( !keep_reference ||
       ( first || record_number<=old_record || oldname!=file_name 
       		|| last!=last_modif ) )
   {
      if( !first && keep_reference ) stream.close() ;
      stream.open( file_name.c_str(), std::ios::binary | std::ios::in ) ;
      if( !stream )
      {
         MAC_Error::object()->raise_plain(
            "Unable to open binary file "+file_name ) ;
      }
      MAC_CHECK( stream.good() ) ;
      stream.ignore( preamble_length ) ;
      MAC_CHECK( stream.good() ) ;
      first = false ;
   }
   oldname = file_name ;
   old_record = record_number ;
   last_modif = last ;
   
   BinaryRecord record ;
   size_t record_size = sizeof( BinaryRecord ) ;
   bool found = false ;
   while( !stream.eof() )
   {
      stream.read( (char*)&record, record_size ) ;
      if( (size_t)stream.gcount()!=record_size )
      {
	 MAC_Error::object()->raise_plain( 
	    "Failed to restore corrupted binary file "+file_name ) ;
      }
      if( record.number==record_number )
      {
         found = true ;
         break ;
      }
      MAC_CHECK( record.length>0 ) ;
      stream.ignore( record.length ) ;
      MAC_CHECK( stream.good() ) ;
   }
   if( !found )
   {
      MAC_Error::object()->raise_plain( "Unable to find record" ) ;
   }
   if( record.type==Double )
   {
      MAC_ASSERT( record.length==sizeof( double ) ) ;
      double val ;
      stream.read( (char*)&val, sizeof( double ) ) ;
      result = MAC_Double::create( a_owner, val ) ;
   }
   else if( record.type==DoubleVector )
   {
      size_t nb = record.length/sizeof( double ) ;
      double *dval = new double[ nb ] ;
      stream.read( (char*)dval, record.length ) ;
      doubleVector dv( nb ) ;
      for( size_t i=0 ; i<nb ; i++ )
      {
         dv(i) = dval[i] ;
      }
      result = MAC_DoubleVector::create( a_owner, dv ) ;
      delete [] dval ;
   }
   else if( record.type==DoubleArray2D )
   {
      size_t nb = record.length/sizeof( double ) ;
      double *dval = new double[ nb ] ;
      stream.read( (char*)dval, record.length ) ;
      size_t dim0 = record.foo[0] ;
      size_t dim1 = record.foo[1] ; 
      MAC_ASSERT( nb==dim0*dim1 ) ;
      doubleArray2D da( dim0, dim1 ) ;
      for( size_t i=0 ; i<dim0 ; i++ )
      {
         for( size_t j=0 ; j<dim1 ; j++ )
         {
            da(i,j) = dval[i*dim1+j] ;
         }
      }
      result = MAC_DoubleArray2D::create( a_owner, da ) ;
      delete [] dval ;
   }
   else if( record.type==DoubleArray3D )
   {

      size_t dim0 = record.foo[0] ;
      size_t dim1 = record.foo[1] ;
      size_t dim2 = record.foo[2] ;
      size_t nb = record.length/sizeof( double ) ;
      MAC_ASSERT( nb==dim0*dim1*dim2 ) ;
      
      double *dval = new double[ nb ] ;
      stream.read( (char*)dval, record.length ) ;
      doubleArray3D da( dim0,dim1,dim2 ) ;
      size_t cpt = 0 ;
      for( size_t i=0 ; i<dim0 ; i++ )
      {
         for( size_t j=0 ; j<dim1 ; j++ )
         {
            for( size_t k=0 ; k<dim2 ; k++ )
            {
               da(i,j,k) = dval[cpt++] ;
            }
         }
      }
      result = MAC_DoubleArray3D::create( a_owner, da ) ;
      delete [] dval ;
   }
   else if( record.type==IntArray2D )
   {
      size_t nb = record.length/sizeof( int ) ;
      int *ival = new int[ nb ] ;
      stream.read( (char*)ival, record.length ) ;
      size_t dim0 = record.foo[0] ;
      size_t dim1 = record.foo[1] ;
      MAC_ASSERT( dim0*dim1==nb ) ;
      intArray2D ia( dim0, dim1 ) ;
      for( size_t i=0 ; i<dim0 ; i++ )
      {
         for( size_t j=0 ; j<dim1 ; j++ )
         {
            ia(i,j) = ival[i*dim1+j] ;
         }
      }
      result = MAC_IntArray2D::create( a_owner, ia ) ;
      delete [] ival ;
   }
   else if( record.type==IntArray3D )
   {

      size_t dim0 = record.foo[0] ;
      size_t dim1 = record.foo[1] ;
      size_t dim2 = record.foo[2] ;
      size_t nb = record.length/sizeof( int ) ;
      MAC_ASSERT( nb==dim0*dim1*dim2 ) ;
      
      int *ival = new int[ nb ] ;
      stream.read( (char*)ival, record.length ) ;
      intArray3D ia( dim0,dim1,dim2 ) ;
      size_t cpt = 0 ;
      for( size_t i=0 ; i<dim0 ; i++ )
      {
         for( size_t j=0 ; j<dim1 ; j++ )
         {
            for( size_t k=0 ; k<dim2 ; k++ )
            {
               ia(i,j,k) = ival[cpt++] ;
            }
         }
      }
      result = MAC_IntArray3D::create( a_owner, ia ) ;
      delete [] ival ;
   }
   else if( record.type==IntVector )
   {
      size_t nb = record.length/sizeof( int ) ;
      int *ival = new int[ nb ] ;
      stream.read( (char*)ival, record.length ) ;
      intVector iv( nb ) ;
      for( size_t i=0 ; i<nb ; i++ )
      {
         iv(i) = ival[i] ;
      }
      result = MAC_IntVector::create( a_owner, iv ) ;
      delete [] ival ;
   }
   else if( record.type==BoolVector )
   {
      size_t nb = record.length/sizeof( bool ) ;
      bool *bval = new bool[ nb ] ;
      stream.read( (char*)bval, record.length ) ;
      boolVector bv( nb ) ;
      for( size_t i=0 ; i<nb ; i++ )
      {
         bv(i) = bval[i] ;
      }
      result = MAC_BoolVector::create( a_owner, bv ) ;
      delete [] bval ;
   }
   if( !keep_reference ) stream.close() ;
   MAC_CHECK_POST( result!=0 ) ;
   return result ;
}




//----------------------------------------------------------------------
void
MAC_BinStored:: init_binary_file( std::string const& file_name )
//----------------------------------------------------------------------
{
   std::ofstream stream( file_name.c_str(), std::ios::binary | std::ios::out ) ;
   if( !stream )
   {
      MAC_Error::object()->raise_plain( "Unable to open "+file_name ) ;
   }
   stream.write( (char*)&magic_stamp, sizeof(magic) ) ;
   stream.close() ;
   MAC_CHECK_POST( is_valid_binary_file( file_name ) ) ;
}




//----------------------------------------------------------------------
MAC_BinStored const*
MAC_BinStored:: create_reference( MAC_Object* a_owner,
                                  MAC_Data const* a_data,
                                  std::string const& file_name,
                                  bool local_reference_file_name ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_BinStored:: create_reference" ) ;
   MAC_CHECK_PRE( a_data!=0 ) ;
   MAC_CHECK_PRE( a_data->value_can_be_evaluated(0) ) ;
   MAC_CHECK_PRE( is_type_supported( a_data->data_type() ) ) ;
   MAC_CHECK_PRE( !file_name.empty() ) ;
   MAC_CHECK_PRE( is_valid_binary_file( file_name ) ) ;

   size_t record_number = last_record_number( file_name ) ;
   if( record_number==bad_record )
   {
      record_number=0 ;
   }
   else
   {
      record_number++ ;
   }
   
   std::ofstream stream( file_name.c_str(), std::ios::binary
                         | std::ios::out | std::ios::app ) ;
   
   if( !stream )
   {
      MAC_Error::object()->raise_plain( "Unable to open "+file_name ) ;
   }
   BinaryRecord record ;
   record.type = a_data->data_type() ;
   record.number = record_number ;
   size_t nb = 0 ;
   if( record.type==Double )
   {
      nb = 1 ;
      record.length=sizeof( double ) ;
   }
   else if( record.type==DoubleVector )
   {
      nb = a_data->to_double_vector().size() ;
      record.length=nb*sizeof( double ) ;
   }
   else if( record.type==DoubleArray2D )
   {
      doubleArray2D const& ref = a_data->to_double_array2D() ;
      nb = ref.index_bound(0)*ref.index_bound(1) ;
      record.length=nb*sizeof( double ) ;
      record.foo[0] = ref.index_bound(0) ;
      record.foo[1] = ref.index_bound(1) ;
   }
   else if( record.type==DoubleArray3D )
   {
      doubleArray3D const& refd3D = a_data->to_double_array3D() ;
      nb = refd3D.index_bound(0)*
         refd3D.index_bound(1)*refd3D.index_bound(2) ;
      record.length=nb*sizeof( double ) ;
      record.foo[0] = refd3D.index_bound(0) ;
      record.foo[1] = refd3D.index_bound(1) ;
      record.foo[2] = refd3D.index_bound(2) ;
   }
   else if( record.type==IntArray2D )
   {
      intArray2D const& refi2D = a_data->to_int_array2D() ;
      nb = refi2D.index_bound(0)*
         refi2D.index_bound(1) ;
      record.length=nb*sizeof( int ) ;
      record.foo[0] = refi2D.index_bound(0) ;
      record.foo[1] = refi2D.index_bound(1) ;
   }
   else if( record.type==IntArray3D )
   {
      intArray3D const& refi3D = a_data->to_int_array3D() ;
      nb = refi3D.index_bound(0)*
         refi3D.index_bound(1)*refi3D.index_bound(2) ;
      record.length=nb*sizeof( int ) ;
      record.foo[0] = refi3D.index_bound(0) ;
      record.foo[1] = refi3D.index_bound(1) ;
      record.foo[2] = refi3D.index_bound(2) ;
   }
   else if( record.type==IntVector )
   {
      nb = a_data->to_int_vector().size() ;
      record.length=nb*sizeof( int ) ;
   }
   else if( record.type==BoolVector )
   {
      nb = a_data->to_bool_vector().size() ;
      record.length=nb*sizeof( bool ) ;
   }
   else
   {
      MAC_Error::object()->raise_plain( "Internal error" ) ;
   }
   
   stream.write( (char*)&record, sizeof( BinaryRecord ) ) ;
   if( record.type==Double )
   {
      double val = a_data->to_double() ;
      stream.write( (char*)&val, record.length ) ;
   }
   else if( record.type==DoubleVector )
   {
      double *dval = new double[ nb ] ;
      doubleVector const& dv = a_data->to_double_vector() ;
      for( size_t i=0 ; i<nb ; i++ )
      {
         dval[i] = dv(i) ;
      }
      stream.write( (char*)dval, record.length ) ;
      delete [] dval ;
   }
   else if( record.type==DoubleArray2D )
   {
      double *dval = new double[ nb ] ;
      doubleArray2D const& da = a_data->to_double_array2D() ;
      size_t nb0 = da.index_bound(0) ;
      size_t nb1 = da.index_bound(1) ;
      
      for( size_t i=0 ; i<nb0 ; i++ )
      {
         for( size_t j=0 ; j<nb1 ; j++ )
         {
            dval[i*nb1+j] = da(i,j) ;
         }
      }
      stream.write( (char*)dval, record.length ) ;
      delete [] dval ;
   }
   else if( record.type==DoubleArray3D )
   {
      double *dval = new double[ nb ] ;
      doubleArray3D const& da = a_data->to_double_array3D() ;
      size_t nb0 = da.index_bound(0) ;
      size_t nb1 = da.index_bound(1) ;
      size_t nb2 = da.index_bound(2) ;
      size_t cpt = 0 ;
      for( size_t i=0 ; i<nb0 ; i++ )
      {
         for( size_t j=0 ; j<nb1 ; j++ )
         {
            for( size_t k=0 ; k<nb2 ; k++ )
            {
               dval[cpt++] = da(i,j,k) ;
            }
         }
      }
      stream.write( (char*)dval, record.length ) ;
      delete [] dval ;
   }
   else if( record.type==IntArray3D )
   {
      int *ival = new int[ nb ] ;
      intArray3D const& ia = a_data->to_int_array3D() ;
      size_t nb0 = ia.index_bound(0) ;
      size_t nb1 = ia.index_bound(1) ;
      size_t nb2 = ia.index_bound(2) ;
      size_t cpt = 0 ;
      for( size_t i=0 ; i<nb0 ; i++ )
      {
         for( size_t j=0 ; j<nb1 ; j++ )
         {
            for( size_t k=0 ; k<nb2 ; k++ )
            {
               ival[cpt++] = ia(i,j,k) ;
            }
         }
      }
      stream.write( (char*)ival, record.length ) ;
      delete [] ival ;
   }
   else if( record.type==IntArray2D )
   {
      int *ival = new int[ nb ] ;
      intArray2D const& ia = a_data->to_int_array2D() ;
      size_t nb0 = ia.index_bound(0) ;
      size_t nb1 = ia.index_bound(1) ;
      
      for( size_t i=0 ; i<nb0 ; i++ )
      {
         for( size_t j=0 ; j<nb1 ; j++ )
         {
            ival[i*nb1+j] = ia(i,j) ;
         }
      }
      stream.write( (char*)ival, record.length ) ;
      delete [] ival ;
   }
   else if( record.type==IntVector )
   {
      int *ival = new int[ nb ] ;
      intVector const& iv = a_data->to_int_vector() ;
      for( size_t i=0 ; i<nb ; i++ )
      {
         ival[i] = iv(i) ;
      }
      stream.write( (char*)ival, record.length ) ;
      delete [] ival ;
   }
   else if( record.type==BoolVector )
   {
      bool *bval = new bool[ nb ] ;
      boolVector const& bv = a_data->to_bool_vector() ;
      for( size_t i=0 ; i<nb ; i++ )
      {
         bval[i] = bv(i) ;
      }
      stream.write( (char*)bval, record.length ) ;
      delete [] bval ;
   }
   
   MAC_List* argument_list = MAC_List::create( 0 ) ;
   argument_list->append(
      MAC_String::create( argument_list, type_name( a_data->data_type() ) ) ) ;
   if( local_reference_file_name )
   {
      MAC_List* args = MAC_List::create( 0 ) ;
      args->append( new MAC_ThisFileDirExp( args ) ) ;
      args->append( MAC_String::create( args, 
      		MAC_System::basename( file_name ) ) ) ;
      MAC_Data* join = MAC_Expression::create( argument_list, "join", args ) ;
      args->set_owner( join ) ;
      argument_list->append( join ) ;
   }
   else
   {
      argument_list->append(
         MAC_String::create( argument_list, file_name  ) ) ;
   }
   argument_list->append(
      MAC_Int::create( argument_list, record_number  ) ) ;
   
   MAC_BinStored* result = new MAC_BinStored( a_owner, argument_list ) ;
   argument_list->set_owner( result ) ;
   if( !stream.good() )
   {
      MAC_Error::object()->raise_plain( "Unable to write "+file_name ) ;
   }   
   stream.close() ;
   
   // CACHING
   last_file = file_name ;
   last_record = record_number ;
   
   MAC_CHECK_POST( result!=0 ) ;
   MAC_CHECK_POST( result->owner()==a_owner ) ;
   
   return result ;
}




//internal---------------------------------------------------------------
MAC_ThisFileDirExp:: MAC_ThisFileDirExp( MAC_Object* a_owner ) 
//internal---------------------------------------------------------------
      : MAC_Expression( a_owner, "this_file_dir", MAC_List::create( a_owner ) )
{
}




//internal---------------------------------------------------------------
MAC_ThisFileDirExp:: ~MAC_ThisFileDirExp( void ) 
//internal---------------------------------------------------------------
{
}




//internal---------------------------------------------------------------
MAC_ThisFileDirExp*
MAC_ThisFileDirExp:: create_replica(
           MAC_Object* a_owner, MAC_Sequence const* argument_list ) const
//internal---------------------------------------------------------------
{
   MAC_Error::object()->raise_internal( "Can not be called" ) ;
   return( 0 ) ;
}




//internal---------------------------------------------------------------
MAC_Data::Type
MAC_ThisFileDirExp:: data_type( void ) const
//internal---------------------------------------------------------------
{
   return( String ) ;
}




//internal---------------------------------------------------------------
std::string const&
MAC_ThisFileDirExp:: to_string( MAC_Context const* ct ) const
//internal---------------------------------------------------------------
{
   MAC_Error::object()->raise_plain(
      "*** MAC_BinStored error:\n"
      "    Expression \"binary\" can not be evaluated.\n"
      "    Use MAC_BinStored::create_reference with\n"
      "    local_reference_file_name = false" ) ;
   static std::string result ;
   return result ;
}




//internal---------------------------------------------------------------
bool
MAC_ThisFileDirExp:: valid_arguments(
                              MAC_Sequence const* some_arguments ) const
//internal---------------------------------------------------------------
{
   return true ;
}




//internal---------------------------------------------------------------
std::string const&
MAC_ThisFileDirExp:: usage( void ) const
//internal---------------------------------------------------------------
{
   MAC_Error::object()->raise_internal( "Can not be called" ) ;
   static std::string result ;
   return result ;
}
