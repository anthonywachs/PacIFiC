#include <MAC_ModuleComparator.hh>

#include <MAC.hh>
#include <MAC_Module.hh>
#include <MAC_KeywordDataPair.hh>
#include <MAC_KeywordDataIterator.hh>
#include <MAC_ModuleExplorer.hh>
#include <MAC_ModuleIterator.hh>
#include <MAC_Data.hh>
#include <MAC_DoubleVector.hh>
#include <MAC_DoubleArray2D.hh>
#include <MAC_DoubleArray3D.hh>
#include <MAC_IntVector.hh>
#include <MAC_IntArray2D.hh>
#include <MAC_String.hh>
#include <MAC_StringVector.hh>
#include <MAC_BoolVector.hh>
#include <MAC_Error.hh>
#include <MAC_assertions.hh>

#include <sstream>
#include <iostream>

using std::endl ;

//----------------------------------------------------------------------
MAC_ModuleComparator*
MAC_ModuleComparator:: create( MAC_Object* a_owner, 
                               MAC_ModuleExplorer const* exp ) 
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ModuleComparator::create(a_owner,envir)" ) ;

   MAC_ModuleComparator* result = new MAC_ModuleComparator( a_owner, exp ) ;

   MAC_CHECK_POST( result != 0 ) ;
   MAC_CHECK_POST( result->owner() == a_owner ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
MAC_ModuleComparator:: MAC_ModuleComparator( MAC_Object* a_owner, 
                                             MAC_ModuleExplorer const* exp )
//----------------------------------------------------------------------
   : MAC_Object( a_owner )
   , verbose( false )
   , valid_module( 0 )
   , valid_data( 0 )
   , ignore_data( 0 )
   , MY_DBL_EPS( 0.0 )
   , MY_DBL_MIN( 0.0 )
{
   if( exp != 0 )
   {
      if( exp->has_entry( "verbose" ) ) verbose = exp->bool_data( "verbose" ) ;
      if( exp->has_entry( "valid_module" ) )
	   valid_module = new stringVector( 
                              exp->stringVector_data( "valid_module" ) ) ;
      if( exp->has_entry( "valid_data" ) )
         valid_data = new stringVector( 
                          exp->stringVector_data( "valid_data" ) ) ;     
      if( exp->has_entry( "ignore_data" ) )
         ignore_data = new stringVector( 
                           exp->stringVector_data( "ignore_data" ) ) ;
      if( exp->has_entry( "dbl_eps" ) )
         MY_DBL_EPS = exp->double_data( "dbl_eps" ) ;
      if( exp->has_entry( "dbl_min" ) )
         MY_DBL_MIN = exp->double_data( "dbl_min" ) ;
   }
}

//----------------------------------------------------------------------
MAC_ModuleComparator:: ~MAC_ModuleComparator( void )
//----------------------------------------------------------------------
{
   if (valid_module) delete valid_module;
   if (valid_data) delete valid_data;
   if (ignore_data) delete ignore_data;
}

//----------------------------------------------------------------------
bool
MAC_ModuleComparator:: is_verbose( void)
//----------------------------------------------------------------------
{
   return verbose;
}

//----------------------------------------------------------------------
bool
MAC_ModuleComparator:: is_valid_module(std::string const& name)
//----------------------------------------------------------------------
{
   if (valid_module == 0) return true;
   for (size_t i=0;i<valid_module->size();i++) {
      if (is_verbose()) {
	 MAC::out() << "checking if " << (*valid_module)(i)
	      << " is in " << name
	      << " : " << name.find((*valid_module)(i))
	      << endl;
      }
      if (name.find((*valid_module)(i)) == 0) return true;
   }
   return false;
}

//----------------------------------------------------------------------
bool
MAC_ModuleComparator:: is_valid_data( std::string const& name,
                                      std::string const& abs_path_name )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_ModuleComparator:: is_valid_data" ) ;
   bool result = false ;
   if( valid_data == 0 )
   {
      result = ignore_data==0 || ( !ignore_data->has(name) ) &&
         ( !ignore_data->has(abs_path_name) ) ;
   }
   else
   {
      result = valid_data->has(name) || valid_data->has(abs_path_name) ;
   }
   
   return result ;
}

//----------------------------------------------------------------------
int
MAC_ModuleComparator:: compare( MAC_Module const* m1, 
                                MAC_Module const* m2, 
                                MAC_Module* result )
//----------------------------------------------------------------------
{
   MAC_ASSERT(m1 != 0);
   MAC_ASSERT(m2 != 0);
   MAC_ASSERT(result != 0);

   left=m1->name();
   right=m2->name();

   return internalCompare(m1, m2, result, "/");
}

//----------------------------------------------------------------------
int
MAC_ModuleComparator:: internalCompare( MAC_Module const* m1, 
                                        MAC_Module const* m2, 
                                        MAC_Module* result, 
                                        std::string const& path ) 
//----------------------------------------------------------------------
{
    if (is_verbose())
       MAC::out() << "checking " << path << endl;

    int nb_err = 0;
    MAC_KeywordDataIterator* kwit;
    MAC_ModuleIterator* modit;

    if (is_verbose())
	MAC::out() << "Looking for data only in " << right << ":" << m2->name() << endl;

    kwit = m2->create_entry_iterator(0);
    for( kwit->start() ; kwit->is_valid() ; kwit->go_next() ) {
	std::string const &dataname = kwit->item()->keyword();
        std::string abs_name = path+"/"+dataname ;
        
	if (!is_valid_data(dataname,abs_name)) continue;
	if (! m1->has_entry("/"+dataname)) {
	    nb_err++;
	    std::string mess("missing data in ");
	    mess += left;
	    result->add_entry(dataname, MAC_String::create(result, mess));
	    if (is_verbose()) MAC::out() << dataname << " :" << mess << endl;
	}
    }
    kwit->destroy();

    if (is_verbose())
	MAC::out() << "Looking for data only in " << left << ":" << m1->name() << endl;

    kwit = m1->create_entry_iterator(0);
    for( kwit->start() ; kwit->is_valid() ; kwit->go_next() ) {
	std::string const &dataname = kwit->item()->keyword();
        std::string abs_name = path+"/"+dataname ;

	if (!is_valid_data(dataname,abs_name)) continue;

	if (! m2->has_entry("/"+dataname)) {
	    nb_err++;
	    std::string mess("missing data in ");
	    mess += right;
	    result->add_entry(dataname, MAC_String::create(result, mess));
	    if (is_verbose()) MAC::out() << dataname << " :" << mess << endl;
	} else {
	    if (is_verbose())
		MAC::out() << "checking data " <<  m1->name() << "/" << dataname << endl;
	    nb_err += compare(dataname, m1, m2, result);
	}
    }
    kwit->destroy();

    // Look for modules only in m2
    if (is_verbose())
	MAC::out() << "Looking for modules only in " << right << ":" << m2->name() << endl;

    modit = m2->create_module_iterator(0);
    for( modit->start() ; modit->is_valid() ; modit->go_next() ) {
	std::string const &modname = modit->item()->name();

	if (!is_valid_module(path+modname+"/")) continue;

	if (! m1->has_module("/"+modname)) {
	    nb_err++;
	    std::string mess("missing module in ");
	    mess += left;
	    result->add_entry(modname, MAC_String::create(result, mess));
	    if (is_verbose()) MAC::out() << modname << " :" << mess << endl;
	}
    }
    modit->destroy();

    // Look for modules only in m1 and remember modules to check
    if (is_verbose())
	MAC::out() << "Looking for modules only in " << left << ":" << m1->name() << endl;

    modit = m1->create_module_iterator(0);
    for( modit->start() ; modit->is_valid() ; modit->go_next() ) {
	std::string const &modname = modit->item()->name();

	if (!is_valid_module(path+modname+"/")) continue;

	if (! m2->has_module("/"+modname)) {
	    nb_err++;
	    std::string mess("missing module in ");
	    mess += right;
	    result->add_entry(modname, MAC_String::create(result, mess));
	    if (is_verbose()) MAC::out() << modname << " :" << mess << endl;
	} else {
	    MAC_Module* subresult = MAC_Module::create(result, modname);
	    // recursive call ...
	    nb_err+=internalCompare(m1->module(modname), m2->module(modname), subresult, path+modname+"/");
	    if (subresult->is_empty())
		result->destroy_possession(subresult);
	    else
		result->add_module(subresult);
	}
    }
    modit->destroy();

    return nb_err;
}

//----------------------------------------------------------------------
int
MAC_ModuleComparator::compare(std::string const& dataname,
			      MAC_Module const* m1,
			      MAC_Module const* m2,
			      MAC_Module* result ) 
//----------------------------------------------------------------------
{

   const MAC_Data *d1 = m1->data_of_entry(dataname);
   const MAC_Data *d2 = m2->data_of_entry(dataname);
   MAC_ASSERT(d1 != 0);
   MAC_ASSERT(d2 != 0);

   int nb_err = 0;

   if (d1->data_type() != d2->data_type()) {
      std::string mess("different type : ");
      mess += d1->type_name(d1->data_type());
      mess += " <!!> ";
      mess += d2->type_name(d2->data_type());
      result->add_entry(dataname, MAC_String::create(result, mess));
      if (is_verbose()) MAC::out() << dataname << " :" << mess << endl;
      nb_err = 1;
   } else {
      MAC_Context const *ct1=m1->context();
      MAC_Context const *ct2=m2->context();
      int err;
      // Case of symbolic expression not evaluable in any case
      bool d1_is_evaluable = d1->value_can_be_evaluated( ct1 ) ;
      bool d2_is_evaluable = d2->value_can_be_evaluated( ct2 ) ;

      if( d1_is_evaluable != d2_is_evaluable ) 
      {
         std::string mess("not evaluable");
         result->add_entry(dataname, MAC_String::create(result, mess));
         if (is_verbose()) MAC::out() << dataname << " :" << mess << endl;
         nb_err = 1;
      }
      else if( !d1_is_evaluable )
      {
         std::ostringstream d1str ;
         std::ostringstream d2str ;
         d1->print( d1str, 0 ) ;
         d2->print( d2str, 0 ) ;
         if( d1str.str() != d2str.str() ) 
         {
            std::string mess("different expression : ");
            mess += d1str.str() ;
            mess += " <!!> ";
            mess += d2str.str() ;
            result->add_entry(dataname, MAC_String::create(result, mess));
            if (is_verbose()) MAC::out() << dataname << " :" << mess << endl;
         }
      }
      else
      {
         
         std::ostringstream mess;

         switch(d1->data_type()) {
            case MAC_Data::Double:
            {
               double v1 = d1->to_double(ct1);
               double v2 = d2->to_double(ct2);
               double adiff = 0.0;
               nb_err += (err = compare(v1, v2, adiff));
               if (err) {
                  mess << "different value -> " 
                       << v1 << " <!!> " << v2 << " <-" ;
                  result->add_entry(dataname, MAC_String::create(result, mess.str()));
               }
            }
            break;
            case MAC_Data::Int:
            {
               int v1 = d1->to_int(ct1);
               int v2 = d2->to_int(ct2);
               nb_err += (err = compare(v1, v2));
               if (err) {
                  mess << "different value -> "
                       << v1 << " <!!> " << v2 << " <-" ;
                  result->add_entry(dataname, MAC_String::create(result, mess.str()));
               }
            }
            break;
            case MAC_Data::Bool:
            {
               bool v1 = d1->to_bool(ct1);
               bool v2 = d2->to_bool(ct2);
               nb_err += (err = compare(v1, v2));
               if (err) {
                  mess << "different value -> "
                       << v1 << " <!!> " << v2 << " <-" ;
                  result->add_entry(dataname, MAC_String::create(result, mess.str()));
               }
            }
            break;
            case MAC_Data::String:
            {
               nb_err += (err=compare(d1->to_string(ct1),d2->to_string(ct2)));
               if (err) {
                  mess << "different value -> "
                       << d1->to_string(ct1) << " <!!> " << d2->to_string(ct2)
                       << "  <-" ;
                  result->add_entry(dataname, MAC_String::create(result, mess.str()));
               }
            }
            break;
            case MAC_Data::DoubleVector:
            {
               doubleVector const& v1 = d1->to_double_vector(ct1);
               doubleVector const& v2 = d2->to_double_vector(ct2);
               if( v1.size() != v2.size() ) 
               {
                  nb_err++;
                  mess << "different size -> " 
                       << v1.size() << " <!!> " << v2.size() << " <-" ;
                  result->add_entry( dataname, 
                                     MAC_String::create(result, mess.str()) ) ;
               }
               else 
               {
                  err = 0 ;
                  doubleArray2D vstat( v1.size(), 3 ) ;
                  for( size_t i=0 ; i<v1.size() ; i++ ) 
                  {
                     double status = MAC::bad_double() ;
                     err += compare( v1( i ), v2( i ), status ) ;
                     vstat( i, 0 ) = status ;
                     vstat( i, 1 ) = v1( i ) ;
                     vstat( i, 2 ) = v2( i ) ;
                  }
                  if( err ) 
                  {
                     nb_err += err;
                     mess << "different value at " << err << " position" ;
                     if( err > 1 ) mess << "s" ;
                     result->add_entry( dataname, 
                             MAC_String::create( result, mess.str() ) ) ;
                     result->add_entry( dataname+"_STAT", 
                             MAC_DoubleArray2D::create( result, vstat ) ) ;
                  }
               }
            }
            break;
            case MAC_Data::IntVector:
            {
               intVector v1 = d1->to_int_vector(ct1);
               intVector v2 = d2->to_int_vector(ct2);
               if (v1.size() != v2.size()) {
                  nb_err++;
                  mess << "different size -> " 
                       << v1.size() << " <!!> " << v2.size() << " <-";
                  result->add_entry(dataname, MAC_String::create(result, mess.str()));
               } else {
                  err=0;
                  intVector vdiff(v1.size());
                  for( size_t i=0 ; i<v1.size() ; i++ ) {
                     err += compare(v1(i),v2(i));
                     vdiff(i) = v1(i)-v2(i);
                  }
                  if (err) {
                     nb_err += err;
                     mess << "different value at " << err << " position" ;
                     if( err > 1 ) mess << "s" ;
                     result->add_entry(dataname, MAC_String::create(result, mess.str()));
                     result->add_entry(dataname+"_DIFF", MAC_IntVector::create(result, vdiff));

                  }
               }
            }
            break;
            case MAC_Data::BoolVector:
            {
               boolVector v1 = d1->to_bool_vector(ct1);
               boolVector v2 = d2->to_bool_vector(ct2);
               if (v1.size() != v2.size()) {
                  nb_err++;
                  mess << "different size -> " 
                       << v1.size() << " <!!> " << v2.size() << " <-" ;
                  result->add_entry(dataname, MAC_String::create(result, mess.str()));
               } else {
                  err=0;
                  boolVector vdiff(v1.size());
                  for( size_t i=0 ; i<v1.size() ; i++ ) {
                     err += compare(v1(i),v2(i));
                     vdiff(i) = v1(i)==v2(i);
                  }
                  if (err) {
                     nb_err += err;
                     mess << "different value at " << err << " position" ;
                     if( err > 1 ) mess << "s" ;
                     result->add_entry(dataname, MAC_String::create(result, mess.str()));
                     result->add_entry(dataname+"_DIFF", MAC_BoolVector::create(result, vdiff));
                  }
               }
            }
            break;
            case MAC_Data::StringVector:
            {
               stringVector v1 = d1->to_string_vector(ct1);
               stringVector v2 = d2->to_string_vector(ct2);
               if (v1.size() != v2.size()) {
                  nb_err++;
                  mess << "different size -> " 
                       << v1.size() << " <!!> " << v2.size() << " <-" ;
                  result->add_entry(dataname, MAC_String::create(result, mess.str()));
               } else {
                  err=0;
                  stringVector vdiff(v1.size());
                  for( size_t i=0 ; i<v1.size() ; i++ ) {
                     err += compare(v1(i),v2(i));
                  }
                  if (err) {
                     nb_err += err;
                     mess << "different value at " << err << " position" ;
                     if( err > 1 ) mess << "s" ;
                     result->add_entry(dataname, MAC_String::create(result, mess.str()));
                     result->add_entry(dataname+"_DIFF", MAC_StringVector::create(result, vdiff));
                  }
               }
            }
            break;
            case MAC_Data::DoubleArray2D:
            {
               doubleArray2D const& v1 = d1->to_double_array2D( ct1 ) ;
               doubleArray2D const& v2 = d2->to_double_array2D( ct2 ) ;
               if( ( v1.index_bound( 0 ) != v2.index_bound( 0 ) ) ||
                   ( v1.index_bound( 1 ) != v2.index_bound( 1 ) ) ) 
               {
                  nb_err++;
                  mess << "different dim -> "
                       << v1.index_bound(0) << "x" << v1.index_bound(1)
                       << " <!!> "
                       << v2.index_bound(0) << "x" << v2.index_bound(1) 
                       << " <-" ;
                  result->add_entry( dataname, 
                          MAC_String::create(result, mess.str() ) ) ;
               } 
               else 
               {
                  err = 0 ;
                  doubleArray3D vstat( v1.index_bound( 0 ), 
                                       v1.index_bound( 1 ), 3 );
                  for( size_t i=0 ; i<v1.index_bound( 0 ) ; i++ ) 
                  {
                     for( size_t j=0 ; j<v1.index_bound( 1 ) ; j++ ) 
                     {
                        double status = MAC::bad_double() ;
                        err += compare( v1( i, j ), v2( i, j ), status ) ;
                        vstat( i, j, 0 ) = status ;
                        vstat( i, j, 1 ) = v1( i, j ) ;
                        vstat( i, j, 2 ) = v2( i, j ) ;
                     }
                  }
                  if( err ) 
                  {
                     nb_err += err;
                     mess << "different value at " << err << " position" ;
                     if( err > 1 ) mess << "s" ;
                     result->add_entry( dataname, 
                             MAC_String::create( result, mess.str() ) ) ;
                     result->add_entry(dataname+"_STAT", 
                             MAC_DoubleArray3D::create( result, vstat ) ) ;
                  }
               }
            }
            break;
            case MAC_Data::DoubleArray3D:
            {
               doubleArray3D const& v1 = d1->to_double_array3D( ct1 ) ;
               doubleArray3D const& v2 = d2->to_double_array3D( ct2 ) ;
               if ( ( v1.index_bound( 0 ) != v2.index_bound( 0 ) ) ||
                    ( v1.index_bound( 1 ) != v2.index_bound( 1 ) ) ||
                    ( v1.index_bound( 2 ) != v2.index_bound( 2 ) ) )
               {
                  nb_err++;
                  mess << "different dim -> "
                       << v1.index_bound(0) << "x" << v1.index_bound(1) 
                       << "x" << v1.index_bound(2) << " <!!> "
                       << v2.index_bound(0) << "x" << v2.index_bound(1) 
                       << "x" << v2.index_bound(2) << " <-";
                  result->add_entry( dataname, 
                          MAC_String::create(result, mess.str() ) ) ;
               } 
               else 
               {
                  err = 0 ;
                  for( size_t i=0 ; i<v1.index_bound(0) ; i++ ) 
                  {
                     for( size_t j=0 ; j<v1.index_bound(1) ; j++ ) 
                     {
                        for( size_t k=0 ; k<v1.index_bound(2) ; k++ ) 
                        {
                           double status = MAC::bad_double() ;
                           err += compare( v1(i,j,k), v2(i,j,k), status ) ;
                        }
                     }
                  }
                  if( err ) 
                  {
                     nb_err += err;
                     mess << "different value at " << err << " position" ;
                     if( err > 1 ) mess << "s" ;
                     result->add_entry( dataname, 
                             MAC_String::create( result, mess.str() ) ) ;
                  }
               }
            }
            break;
            case MAC_Data::IntArray2D:
            {
               intArray2D v1 = d1->to_int_array2D(ct1);
               intArray2D v2 = d2->to_int_array2D(ct2);
               if ((v1.index_bound(0) != v2.index_bound(0))
                   || (v1.index_bound(1) != v2.index_bound(1))) {
                  nb_err++;
                  mess << "different dim -> "
                       << v1.index_bound(0) << "x" << v1.index_bound(1)
                       << " <!!> "
                       << v2.index_bound(0) << "x" << v2.index_bound(1) << " <-";
                  result->add_entry(dataname, MAC_String::create(result, mess.str()));
               } else {
                  err=0;
                  intArray2D vdiff(v1.index_bound(0), v1.index_bound(1));
                  for( size_t i=0 ; i<v1.index_bound(0) ; i++ ) {
                     for( size_t j=0 ; j<v1.index_bound(1) ; j++ ) {
                        err += compare(v1(i,j),v2(i,j));
                        vdiff(i,j) = v1(i,j)-v2(i,j);
                     }
                  }
                  if (err) {
                     nb_err += err;
                     mess << "different value at " << err << " position" ;
                     if( err > 1 ) mess << "s" ;
                     result->add_entry(dataname, MAC_String::create(result, mess.str()));
                     result->add_entry(dataname+"_DIFF", MAC_IntArray2D::create(result, vdiff));
                  }
               }
            }
            break;
            default:
            {
               mess << "unknown data type : can not check";
               result->add_entry(dataname, MAC_String::create(result, mess.str()));
               nb_err=1;
            }
         }
      }
   }
   return nb_err;
}

//----------------------------------------------------------------------
int
MAC_ModuleComparator:: compare( double v1, double v2, double& status ) 
//----------------------------------------------------------------------
{
   bool eq = MAC::double_equality( v1, v2, MY_DBL_EPS, MY_DBL_MIN, status ) ;
   return( eq ? 0 : 1 ) ;
}

//----------------------------------------------------------------------
int
MAC_ModuleComparator:: compare( bool v1, bool v2 ) 
//----------------------------------------------------------------------
{
   return (v1 == v2) ? 0 : 1;
}

//----------------------------------------------------------------------
int
MAC_ModuleComparator:: compare( int v1, int v2 ) 
//----------------------------------------------------------------------
{
   return (v1 == v2) ? 0 : 1;
}

//----------------------------------------------------------------------
int
MAC_ModuleComparator:: compare( std::string const& v1, std::string const& v2 ) 
//----------------------------------------------------------------------
{
   return ( (v1 == v2) ? 0 : 1 ) ;
}
