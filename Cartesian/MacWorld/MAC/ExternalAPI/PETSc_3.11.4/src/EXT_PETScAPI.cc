#include <EXT_PETScAPI.hh>

#include <MAC_assertions.hh>
#include <MAC_Bool.hh>
#include <MAC_Root.hh>
#include <MAC_DataWithContext.hh>
#include <MAC_Error.hh>
#include <MAC_Exec.hh>
#include <MAC_Variable.hh>
#include <MAC_String.hh>
#include <MAC_System.hh>

#include <sstream>
#include <iostream>

#if( ! ( PETSC_VERSION_MAJOR    == 3 && \
         PETSC_VERSION_MINOR    == 11 && \
         PETSC_VERSION_SUBMINOR == 4 ) )
 "Bad version of PETSC ( Version 3.11.4 should be used )" ;
#endif

EXT_PETScAPI* EXT_PETScAPI:: SINGLETON = new EXT_PETScAPI() ;
MAC_Timer* EXT_PETScAPI:: timer = 0 ;


//----------------------------------------------------------------------
EXT_PETScAPI:: EXT_PETScAPI( void )
//----------------------------------------------------------------------
   : MAC_ExternalAPI( "EXT_PETScAPI", 0 )
{
   MAC_Bool* val = MAC_Bool::create( 0, true ) ;
   MAC_Exec::add_variable_to_execution_context(
                         MAC_Variable::object( "BS_with_PETSc" ), val ) ;

   std::ostringstream ver ;
   ver << PETSC_VERSION_MAJOR << "." 
       << PETSC_VERSION_MINOR << "." 
       << PETSC_VERSION_SUBMINOR ;   
   MAC_String* rev = MAC_String::create( 0, ver.str() ) ;
   MAC_Exec::add_variable_to_execution_context(
                         MAC_Variable::object( "SS_PETSc_REV" ), rev ) ;
}




//----------------------------------------------------------------------
EXT_PETScAPI:: ~EXT_PETScAPI( void )
//----------------------------------------------------------------------
{
   PetscFinalize() ;
}




//----------------------------------------------------------------------
void
EXT_PETScAPI:: initialize( int& argc, char **& argv )
//----------------------------------------------------------------------
{
   MAC_LABEL( "EXT_PETScAPI:: initialize" ) ;
   const char* help = "EXT_PETScAPI:: initialize" ;

   char** my_argv = new char* [ argc ] ;
   int my_argc=0 ;
   char** new_argv = new char* [ argc ] ;
   int new_argc = 0 ;
   
   if( argc>0 )
   {
      my_argv[my_argc++] = argv[0] ;
   }
   
   for( int i=0 ; i<argc ; ++i )
   {
      std::string str = argv[i] ;
      if(  str == "-Xpetsc" )
      {
         i++ ;
         MAC_ASSERT( i < argc ) ;
         std::string str1 = argv[i] ;
         if( str1 == "-trace" )
         {
            timer = MAC_Timer::create( MAC_Root::object() ) ;
         }
         else
         {
            my_argv[my_argc++] = argv[i] ;
         }  
      }
      else
      {
	 new_argv[new_argc++] = argv[i] ;
      }
   }
   if( my_argc > 1 || timer!=0 )
   {
      argc = new_argc ;
      argv = new_argv ;
   }
   else
   {
      delete [] new_argv ;
   }
   if( my_argc > 1 )
   {
      std::cout << "PETSc init : " << std::endl ;
      for( int i=1 ; i<my_argc ; ++i ) 
      {
         std::cout << "    " << my_argv[i] << std::endl ;
      }
   }
   PetscInitialize( &my_argc, &my_argv, PETSC_NULL, help ) ;
   delete [] my_argv ;
}




//----------------------------------------------------------------------
bool
EXT_PETScAPI:: parse_options( MAC_ModuleExplorer const* exp,
                              bool verbose )
//----------------------------------------------------------------------
{
   MAC_LABEL( "EXT_PETScAPI:: parse_options" ) ;
   bool result =  exp->has_module( "options" ) ;
   
   if( result )
   {
      MAC_ModuleExplorer * sexp = exp->create_subexplorer( 0, "options" )  ;
      for( sexp->start_entry_iterator() ;
           sexp->is_valid_entry()   ;
           sexp->go_next_entry() )
      {
         std::string const& name = "-"+sexp->keyword() ;
         MAC_Data * data = sexp->data( 0 ) ;
         std::string const& val = data->to_string() ;
         if( verbose )
            std::cout << "EXT_PETScAPI - setting option: " 
                      << name << " " << val << std::endl ;
         
         if( val.empty() )
         {
            PetscOptionsSetValue( PETSC_NULL, name.c_str(), PETSC_NULL ) ;
         }
         else
         {
            PetscOptionsSetValue( PETSC_NULL, name.c_str(), val.c_str() ) ;
         }
         
         data->destroy() ;
      }
      sexp->destroy() ;
   }
   return result ;
}




//----------------------------------------------------------------------
void
EXT_PETScAPI:: going_to_do(  char const* action )
//----------------------------------------------------------------------
{
   MAC_LABEL( "EXT_PETScAPI:: going_to_do" ) ;
   if( timer!=0 )
   {
      MAC::out()<<"["<<MAC_System::used_memory()/1024/1024<<"Mo]"
                <<"Start -> "<<action<<std::endl ;
      if( timer->is_running() )
      {
         timer->stop() ;
         MAC::out() << "*** Warning : imbricated times" << std::endl
                    << "    PETSc timer stopped at cumulative time : (s) " << timer->time() << std::endl
                    << "    Restart PETSc timer..." << std::endl ;
      }
      timer->start() ;
   }
}




//----------------------------------------------------------------------
void
EXT_PETScAPI:: verify( char const * action, int result )
//----------------------------------------------------------------------
{
   MAC_LABEL( "EXT_PETScAPI:: verify" ) ;
   if( timer!=0 )
   {
      MAC::out() << "[" << MAC_System::used_memory()/1024/1024<<"Mo]"
                 << action << " <- End." ;
      if( timer->is_running() )
      {
         timer->stop() ;
         MAC::out() << " PETSc cumulative time : (s) " << timer->time() ;
      }
      MAC::out() << std::endl << std::endl ;
   }
   if( result!=0 )
   {
      std::string mess = "Internal Petsc error encountered in " ;
      mess += action ;
      MAC_Error::object()->raise_internal( mess ) ;
   }
}

