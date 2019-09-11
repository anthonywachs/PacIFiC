#include <MAC_Comparator.hh>

#include <MAC.hh>
#include <MAC_assertions.hh>
#include <MAC_Data.hh>
#include <MAC_Double.hh>
#include <MAC_DoubleVector.hh>
#include <MAC_Error.hh>
#include <MAC_Exec.hh>
#include <MAC_FileToModule.hh>
#include <MAC_Int.hh>
#include <MAC_KeywordDataPair.hh>
#include <MAC_KeywordDataIterator.hh>
#include <MAC_ModuleExplorer.hh>
#include <MAC_ModuleIterator.hh>
#include <MAC_StringVector.hh>
#include <MAC_System.hh>
//#include <MAC_TICio.hh>
#include <MAC_String.hh>
#include <doubleArray2D.hh>
#include <doubleArray3D.hh>
#include <doubleVector.hh>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>

using std::cout ;
using std::endl ;
using std::string ;

MAC_Comparator const* MAC_Comparator::PROTOTYPE = new MAC_Comparator() ;

struct MAC_Comparator_ERROR
{
   static void n0( std::string const& filename ) ;
   static void n1( std::string const& file1, std::string const& format1,
                   std::string const& file2, std::string const& format2 ) ;
   static void n2( MAC_ModuleExplorer const* exp, 
                   std::string const& choices ) ;
} ;

//----------------------------------------------------------------------------
MAC_Comparator:: MAC_Comparator( void )
//----------------------------------------------------------------------------
   : MAC_Application( "maccmp" )
   , FILE1()
   , FILE2()
   , FILE_OUT()
   , EXP( 0 )
   , IGNORE( 0 )
{
}

//----------------------------------------------------------------------------
MAC_Comparator*
MAC_Comparator:: create_replica( MAC_Object* a_owner,
                                 MAC_ModuleExplorer const* exp ) const
//----------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Comparator:: create_replica" ) ;
   MAC_CHECK( create_replica_PRE( a_owner, exp ) ) ;
   
   MAC_Comparator* result = new MAC_Comparator( a_owner, exp ) ;

   MAC_CHECK( create_replica_POST( result, a_owner, exp ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
MAC_Comparator:: MAC_Comparator( MAC_Object* a_owner, 
                                 MAC_ModuleExplorer const* exp ) 
//----------------------------------------------------------------------------
   : MAC_Application( a_owner, exp )
   , FILE1( exp->string_data( "left_file" ) )
   , FILE2( exp->string_data( "right_file" ) )
   , FILE_OUT()
   , EXP( 0 )
   , VERBOSE( false )
   , MY_DBL_EPS( MAC::bad_double() )
   , MY_DBL_MIN( MAC::bad_double() )
   , IGNORE( 0 )
   , FORMAT( "UNKNOWN" )
{
   if( exp->has_entry( "output_file" ) )
      FILE_OUT = exp->string_data( "output_file" ) ;
   if( exp->has_entry( "verbose" ) )
      VERBOSE = exp->bool_data( "verbose" ) ;
   if( exp->has_module( "MAC_ModuleComparator" ) )
      EXP = exp->create_subexplorer( this, "MAC_ModuleComparator" ) ;
   if( exp->has_module( "tolerance" ) )
   {
      MAC_ModuleExplorer const* ee = 
                         exp->create_subexplorer( 0 ,"tolerance" ) ;
      MY_DBL_EPS = ee->double_data( "dbl_eps" ) ;
      MY_DBL_MIN = ee->double_data( "dbl_min" ) ;
      ee->destroy() ;
   }
   if( exp->has_entry( "ignore_data" ) )
      IGNORE = exp->stringVector_data( "ignore_data" ) ;
   if( exp->has_entry( "format" ) )
   {
      FORMAT = exp->string_data( "format" ) ;
   }
}

//----------------------------------------------------------------------------
MAC_Comparator*
MAC_Comparator:: create_replica_from_args( MAC_Object* a_owner,
                                           stringVector& args ) const
//----------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Comparator:: create_replica_from_args" ) ;
   
   MAC_Comparator* result = new MAC_Comparator( a_owner, args ) ;

   MAC_CHECK( create_replica_from_args_POST( result, a_owner, args ) ) ;
   return( result ) ;
}

//----------------------------------------------------------------------------
MAC_Comparator:: MAC_Comparator( MAC_Object* a_owner, stringVector& args )
//----------------------------------------------------------------------------
   : MAC_Application( a_owner, 0 )
   , FILE1()
   , FILE2()
   , FILE_OUT()
   , EXP( 0 )
   , VERBOSE( false )
   , MY_DBL_EPS( MAC::bad_double() )
   , MY_DBL_MIN( MAC::bad_double() )
   , IGNORE( 0 )
   , FORMAT( "UNKNOWN" )
{
   size_t i_non_opt = 0 ;
   while( args.size() != 0 )
   {
      string arg = args( 0 ) ;
      args.remove_at( 0 ) ;

      // no options
      if( arg[0] == '-' )
      {         
         if( arg == "-verbose" )
         {
            VERBOSE = true ;
         }
         else if( arg=="-dbl_eps" && args.size()>0 )
         {
            std::istringstream is( args(0) ) ;
            is >> MY_DBL_EPS ;
            args.remove_at( 0 ) ;
         }
         else if( arg=="-dbl_min" && args.size()>0 )
         {
            std::istringstream is( args(0) ) ;
            is >> MY_DBL_MIN ;
            args.remove_at( 0 ) ;
         }
         else if( arg=="-ignore_data" && args.size()>0 )
         {
            IGNORE.append( args( 0 ) ) ;
            args.remove_at( 0 ) ;
         }
         else if( arg=="-format" && args.size()>0 )
         {
            FORMAT = args( 0 ) ;
            args.remove_at( 0 ) ;
         }
         else if( arg=="-output_file" && args.size()>0 )
         {
            FILE_OUT = args( 0 ) ;
            args.remove_at( 0 ) ;
         }
         else
         {         
            notify_error_in_arguments() ;
         }
      }
      else
      {
         i_non_opt++ ;
         if( i_non_opt==1 ) 
            FILE1 = arg ;
         else if( i_non_opt==2 )
            FILE2 = arg ;
         else if( i_non_opt==3 )
            FILE_OUT = arg ;
         else
            notify_error_in_arguments() ;
      }
   }
   if( ! ( ( MY_DBL_EPS==MAC::bad_double() && MY_DBL_MIN==MAC::bad_double() ) ||
           ( MY_DBL_EPS!=MAC::bad_double() && MY_DBL_MIN!=MAC::bad_double() ) ) )
       notify_error_in_arguments() ;
   if( FILE1.empty() || FILE2.empty() ) notify_error_in_arguments() ;
}

//----------------------------------------------------------------------------
MAC_Comparator:: ~MAC_Comparator( void )
//----------------------------------------------------------------------------
{
}

//----------------------------------------------------------------------------
void
MAC_Comparator:: set_preferred_motifs_formats( 
                                          MAC_ModuleExplorer const* exp ) 
//----------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Comparator:: set_preferred_motifs_formats" ) ;
   MAC_CHECK_PRE( exp != 0 ) ;
   
   pref_motifs()  = exp->stringVector_data( "motifs" ) ;
   pref_formats() = exp->stringVector_data( "formats" ) ;
   exp->test_data( "formats", "size(formats)=size(motifs)" ) ;
   
   for( size_t i=0 ; i<pref_formats().size() ; i++ )
   {
      std::string const& frmt = pref_formats()( i ) ;
      if( ! MAC_FileToModule::has( frmt ) ) 
         MAC_Comparator_ERROR::n2( exp, MAC_FileToModule::list_of_formats() ) ;
   }
}

//----------------------------------------------------------------------------
void
MAC_Comparator:: detect_file_format( MAC_ModuleExplorer const* exp,
                                     std::string const& file_name,
                                     std::string& format )
//----------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Comparator:: detect_file_format" ) ;

   if( exp!=0 && exp->has_entry( "format" ) )
   {
      format = exp->string_data( "format" ) ;
      exp->test_data_in( "format", MAC_FileToModule::list_of_formats() ) ;
   }
   else
   {
      bool found = false ;
      for( size_t n=0 ; n<pref_motifs().size() ; n++ )
      {
         if( file_name.find( pref_motifs()(n) ) < file_name.length() )
         {
            found = true ;
            format = pref_formats()(n) ;
            break ;
         }
      }
      if( !found ) MAC_FileToModule::find_file_format( file_name, format ) ;
   }
}

//----------------------------------------------------------------------------
void
MAC_Comparator:: run( void )
//----------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Comparator::run" ) ;

   std::ofstream ofs ;
   if( !FILE_OUT.empty() )
   {
      ofs.open( FILE_OUT.c_str(), std::ios::out|std::ios::trunc ) ;
      if( ! ofs )
      {
         string mess = "MAC_Comparator: unable to open file \"" ;
         mess += FILE_OUT ;
         mess += "\"" ;
         MAC_Error::object()->raise_plain( mess ) ;
      }
   }
   
   if( FORMAT == "UNKNOWN" )
   {
      std::string fmt = FORMAT ;
      detect_file_format( 0, FILE1, fmt ) ;
      detect_file_format( 0, FILE2, FORMAT ) ;
      if( FORMAT != fmt ) MAC_Comparator_ERROR::n1( FILE1, fmt, 
                                                    FILE2, FORMAT ) ;
   }

   if( VERBOSE ) MAC::out() << "|    format: " << FORMAT << endl ;
   if( FORMAT != "UNKNOWN" )
   {
      if( FILE_OUT.empty() )
         builtin_compare( MAC::out() ) ;
      else
         builtin_compare( ofs ) ;
   }
   else
   {
      if( FILE_OUT.empty() )
         do_diff( MAC::out() ) ;
      else
         builtin_compare( ofs ) ;
   }
   
   if( !FILE_OUT.empty() ) ofs.close() ;
}

//-----------------------------------------------------------------------------
void
MAC_Comparator:: builtin_compare( std::ostream& os )
//-----------------------------------------------------------------------------
{  
   MAC_LABEL( "MAC_Comparator:: builtin_compare" ) ;
   
   MAC_FileToModule const* oo = MAC_FileToModule::object( FORMAT ) ;
   MAC_Module* mod1 = oo->create_from_file( 0, "left_file",  FILE1 ) ;
   MAC_Module* mod2 = oo->create_from_file( 0, "right_file", FILE2 ) ;

   if( EXP==0 ) 
   {
      MAC_Module* comp = MAC_Module::create( this, "MAC_ModuleComparator" ) ;
      
      if( IGNORE.size()>0  )
      {
         comp->add_entry( "ignore_data", 
                          MAC_StringVector::create( comp, IGNORE ) ) ;
      }
      if( MY_DBL_EPS != MAC::bad_double()  )
      {
         comp->add_entry( "dbl_eps", MAC_Double::create( comp, MY_DBL_EPS ) ) ;
      }
      if( MY_DBL_MIN != MAC::bad_double()  )
      {
         comp->add_entry( "dbl_min", MAC_Double::create( comp, MY_DBL_MIN ) ) ;
      }
      if( comp!=0 ) EXP = MAC_ModuleExplorer::create( this, comp ) ;     
   }
   
   if( EXP!=0 && VERBOSE )
   {
      MAC::out() << "MAC_ModuleComparator configuration : " << std::endl ;
      EXP->print( MAC::out(), 2 ) ;
      MAC::out() << std::endl ;
   }
   
   MAC_Module* mdif = MAC_Module::create_as_difference( 
                                            this, "result", mod1, mod2, EXP );

   mdif->add_entry( "left_file",  MAC_String::create( mdif, FILE1 ) ) ;
   mdif->add_entry( "right_file", MAC_String::create( mdif, FILE2 ) ) ;

   int nb_errors = mdif->data_of_entry( "nb_differences" )->to_int() ;   
   if( nb_errors>0 )
   {
      if( FILE_OUT.empty() ) // for consistance with the old versions
      {
         os << "---> " << nb_errors << " difference" ;
         if( nb_errors > 1 ) os << "s" ;
         os << endl;
      }
      display_glob_info_diff( os, mdif ) ;
      display_submodule_diff( os, mdif ) ;
      MAC_Exec::set_exit_code( nb_errors ) ;
   }
   
   mod1->destroy() ; mod1 = 0 ;
   mod2->destroy() ; mod2 = 0 ;
}

//----------------------------------------------------------------------------
void
MAC_Comparator:: do_diff( std::ostream& os ) 
//----------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Comparator:: do_diff" ) ;
   if( VERBOSE ) MAC::out() << "|    do_diff " << FILE1 
                       << " and " << FILE2 << std::endl ;
   
   std::ifstream ifs1( FILE1.c_str() ) ;
   std::ifstream ifs2( FILE2.c_str() ) ;

   bool ok = ifs1.good() && ifs2.good() ;
   std::string fline ;
   std::string line1, line2 ;
   size_t iline = 0 ;
   
   while( ifs1.good() && ifs2.good() )
   {
      iline++ ;
      getline( ifs1, line1 ) ;
      getline( ifs2, line2 ) ;
      if( line1 != line2 ) 
      {
         // For linefeed in files written under windows
         MAC::replace( line1, "\r", "" ) ;
         MAC::replace( line2, "\r", "" ) ;
         
         // because windows notation is not the same than linux one...
         MAC::replace( line1, "E+00", "E+0" ) ;
         MAC::replace( line2, "E+00", "E+0" ) ;
         MAC::replace( line1, "E-00", "E-0" ) ;
         MAC::replace( line2, "E-00", "E-0" ) ;
         MAC::replace( line1, "e+00", "e+0" ) ;
         MAC::replace( line2, "e+00", "e+0" ) ;
         MAC::replace( line1, "e-00", "e-0" ) ;
         MAC::replace( line2, "e-00", "e-0" ) ;
         MAC::replace( line1, "  ", " " ) ;
         MAC::replace( line2, "  ", " " ) ;

         if( line1 != line2 )
         {
            os << "   first difference at line " << iline << endl ;
            if( line1.length() < 100 && line2.length() < 100 )
            {
               os << "   \"" << line1 << "\"" << endl ;
               os << "   \"" << line2 << "\"" << endl ;
            }
            MAC_Exec::set_exit_code( 12 ) ;
            ok = false ;
            break ;
         }
         
      }
   }
   if( ok && ( ifs1.good() || ifs2.good() ) ) 
   {
      MAC_Exec::set_exit_code( 12 ) ;
   } 
}

//---------------------------------------------------------------------------
void
MAC_Comparator:: print_usage( void ) const
//---------------------------------------------------------------------------
{
   MAC::out() << usage_title( "maccmp" )  ;
   MAC::out() << "[options...] file1 file2 [file3]" << endl << endl ;
   MAC::out() << "     Compare and display differences between pairs of files"
              << endl ;
}

//---------------------------------------------------------------------------
void
MAC_Comparator:: print_operands( void ) const
//---------------------------------------------------------------------------
{
   MAC::out() << "     -dbl_eps <value>" << endl
	<< "          argument a_dbl_eps in calls to MAC::double_equality"
        << endl 
        << "          when comparing values of type double" 
        << endl << endl ;
   MAC::out() << "     -dbl_min <value>" << endl
	<< "          argument a_dbl_min in calls to MAC::double_equality"
        << endl 
        << "          when comparing values of type double" 
        << endl << endl ;
   MAC::out() << "     -ignore_data <keyword>" << endl
	<< "          Ignore entries of keyword <keyword>" 
        << endl << endl ;
   MAC::out() << "     -format [" << MAC_FileToModule::list_of_formats() 
                                  << "]" << endl
	<< "          Specify the format of file1 and file2" 
        << endl << endl ;
   MAC::out() << operands_title() ;
   MAC::out() << "     file1, file2" << endl
	<< "          Path names of the two files to be compared" 
        << endl << endl ;
   MAC::out() << "     file3" << endl
        << "          Path names of the file containing the difference on exit"
	<< endl << endl ;
}

//---------------------------------------------------------------------------
void
MAC_Comparator:: print_exit_status( void ) const
//---------------------------------------------------------------------------
{
   MAC::out() << exit_status_title() ;
   MAC::out() << "     0    The files are identical" << endl ;
   MAC::out() << "    >0    Number of differences found" << endl ;
}

//----------------------------------------------------------------------------
void
MAC_Comparator:: display_glob_info_diff( std::ostream& os,
                                         MAC_Module const* mod ) const
//----------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Comparator:: display_glob_info_diff" ) ;

   MAC_KeywordDataIterator* it = mod->create_entry_iterator( 0 ) ;
   for( it->start() ; it->is_valid() ; it->go_next() )
   {
      string const& keywd  = it->item()->keyword() ;
      MAC_Data const* data = it->item()->data() ;
      if( keywd!="left_file" && keywd!="right_file" &&
          keywd!="nb_differences" )
      {
         MAC_ASSERT( data->data_type() == MAC_Data::String ) ;
         os << keywd << ": " << data->to_string() << endl ;
      }
   }
   it->destroy() ; it = 0 ;
}

//----------------------------------------------------------------------------
void
MAC_Comparator:: display_submodule_diff( std::ostream& os,
                                         MAC_Module const* mod ) const
//----------------------------------------------------------------------------
{
   MAC_LABEL( "MAC_Comparator:: display_submodule_diff" ) ;

   MAC_ModuleIterator* it = mod->create_module_iterator( 0 ) ;
   for( it->start() ; it->is_valid() ; it->go_next() )
   {
      MAC_Module const* mm = it->item() ;
      string modname =  mm->absolute_path_name() ;
      modname.erase( modname.begin(), modname.begin()+8 ) ; 
      MAC_KeywordDataIterator* kwit = mm->create_entry_iterator( 0 ) ;
      bool first_diff = true ;
      for( kwit->start() ; kwit->is_valid() ; kwit->go_next() )
      {
         string const& str = kwit->item()->keyword() ;
         MAC_Data const* data = kwit->item()->data() ;
         if( first_diff && data->data_type()!=MAC_Data::DoubleVector )
         {
            os << modname << endl ;
            first_diff=false ;
         }
         std::ios_base::fmtflags original_flags = os.flags() ;
         os.setf( std::ios_base::uppercase | std::ios_base::scientific ) ;
         std::streamsize p = os.precision() ;
         os << std::setprecision( 5 ) ;
         switch( data->data_type() )
         {
            case MAC_Data::DoubleArray2D :
            {
               // the compared data is of type MAC_Data::DoubleVector
               double max = -MAC::max_double() ;
               doubleArray2D const& stat = data->to_double_array2D() ;
               size_t i_max = MAC::bad_index() ;
               size_t i_other = MAC::bad_index() ;
               for( size_t i=0 ; i<stat.index_bound( 0 ) ; ++i )
               {
                  if( ( stat( i, 0 ) < 0.0 ) && ( i_other == MAC::bad_index() ) )
                  {
                     i_other = i ;
                  }
                  if( stat( i,0 ) > max )
                  {
                     i_max = i ;
                     max = stat( i, 0 ) ;
                  }
               }
               if( i_max != MAC::bad_index() &&  max > 0.0 )
               {  
                  os << "    relative error max : " << max << endl ;
                  os << "    obtained for -> " << stat( i_max, 1 ) << " <!!> ";
                  os << stat( i_max, 2 ) << " <-" << endl ;
               }
               else
               {
                  os << "    first difference at index " << i_other << endl ;
                  os << "    -> " << stat( i_other, 1 ) << " <!!> ";
                  os << stat( i_other, 2 ) << " <-" << endl ;
               }
            }
            break ;
            case MAC_Data::DoubleArray3D :
            {
               // the compared data is of type MAC_Data::DoubleArray2D
               double max = -MAC::max_double() ;
               doubleArray3D const& stat = data->to_double_array3D() ;
               size_t i_max = MAC::bad_index() ;
               size_t j_max = MAC::bad_index() ;
               size_t i_other = MAC::bad_index() ;
               size_t j_other = MAC::bad_index() ;
               for( size_t i=0 ; i<stat.index_bound( 0 ) ; ++i )
               {
                  for( size_t j=0 ; j<stat.index_bound( 1 ) ; ++j )
                  {
                     if( ( stat( i, j, 0 ) < 0.0 ) && 
                         ( i_other == MAC::bad_index() ) )
                     {
                        i_other = i ;
                        j_other = j ;
                     }
                     if( stat( i, j , 0 ) > max ) 
                     {
                        i_max = i ;
                        j_max = j ;
                        max = stat( i, j, 0 ) ;
                     }
                  }
               }
               if( i_max != MAC::bad_index() &&  max > 0.0 )
               {  
                  os << "    relative error max : " << max << endl ;
                  os << "    obtained for -> " << stat( i_max, j_max, 1 ) 
                     << " <!!> " << stat( i_max, j_max, 2 ) << " <-" << endl ;
               }
               else
               {
                  os << "    first difference at indices (" 
                     << i_other << "," << j_other << ")" << endl ;
                  os << "    -> " << stat( i_other, j_other, 1 ) << " <!!> ";
                  os << stat( i_other, j_other, 2 ) << " <-" << endl ;
               }
            }
            break ;
            case MAC_Data::String :
            {
               os << " " << str << " : " << data->to_string() << endl ;
            }
            default:
            {
            }
         }
         os << std::setprecision(p) ;
         os.flags( original_flags ) ;
      }
      kwit->destroy() ;
      display_submodule_diff( os, mm ) ; 
   }
   it->destroy() ;
}

//----------------------------------------------------------------------
stringVector&
MAC_Comparator:: pref_motifs( void )
//----------------------------------------------------------------------
{
   static stringVector result( 0 ) ;
   return( result ) ;
}

//----------------------------------------------------------------------
stringVector&
MAC_Comparator:: pref_formats( void )
//----------------------------------------------------------------------
{
   static stringVector result( 0 ) ;
   return( result ) ;
}

//internal--------------------------------------------------------------
void 
MAC_Comparator_ERROR:: n0( std::string const& filename )
//internal--------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << "*** MAC_Comparator error:" << endl << endl ;
   mesg << "    unable to open file \"" <<  filename << "\"" ;
   MAC_Error::object()->raise_plain( mesg.str() ) ;
}

//internal--------------------------------------------------------------
void 
MAC_Comparator_ERROR:: n1( std::string const& file1,
                           std::string const& format1,
                           std::string const& file2,
                           std::string const& format2 )
//internal--------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << "*** MAC_Comparator error:" << endl << endl ;
   mesg << "the two files should have the same format" << endl << endl ;
   mesg << file1 << endl << "   detected format: " << format1 << endl << endl ;
   mesg << file2 << endl << "   detected format: " << format2 << endl ;
   MAC_Error::object()->raise_plain( mesg.str() ) ;
}

//internal--------------------------------------------------------------
void 
MAC_Comparator_ERROR:: n2( MAC_ModuleExplorer const* exp,
                           std::string const& choices )
//internal--------------------------------------------------------------
{
   std::ostringstream mesg ;
   mesg << "In module: " << endl << "   " << exp->absolute_path_name() << endl ;
   mesg << "the data of keyword \"formats\"" << endl ;
   mesg << "should contains strings with one of the following values: " ;
   stringVector liste( 0 ) ;
   liste.build_from_string( choices, ',' ) ;
   for( size_t i=0 ; i<liste.size() ; i++ )
   {
      mesg << endl << "   - \"" << liste(i) << "\"" ;
   }
   MAC_Error::object()->raise_plain( mesg.str() ) ;
}



