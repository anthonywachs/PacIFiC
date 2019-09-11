%{
#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdio>
using std::istream ;
   
#include <MAC_assertions.hh>
#include <MAC_Context.hh>
#include <MAC_Double.hh>
#include <MAC_Map.hh>
#include <MAC_Lexical.hh>
#include <MAC_Expression.hh>
#include <MAC_Module.hh>
#include <MAC_KeywordDataPair.hh>
#include <MAC_Error.hh>
#include <MAC_Int.hh>
#include <MAC_IntArray2D.hh>
#include <MAC_List.hh>
#include <MAC_ListIdentity.hh>
#include <MAC_Root.hh>
#include <MAC_Data.hh>
#include <MAC_DoubleArray2D.hh>
#include <MAC_DoubleVector.hh>
#include <MAC_IntVector.hh>
#include <MAC_BoolArray2D.hh>
#include <MAC_BoolVector.hh>
#include <MAC_StringArray2D.hh>
#include <MAC_StringVector.hh>
#include <MAC_String.hh>
#include <MAC_System.hh>
#include <MAC_Variable.hh>
#include <MAC_Bool.hh>
#include <MAC.hh>

#include <stringVector.hh>

/*-----------------------------------------------.
| Public functions.                              |
`-----------------------------------------------*/
bool MAC_readFile( MAC_Module * top,
                   istream* input_stream,
                   std::string const& name,
                   bool debug = false ) ;
std::string MAC_current_parsed_module_path_name( void ) ;
int MAC_current_parsed_line( void ) ;
void MAC_re_init_lexer( void ) ;


#define YYSTYPE MAC_LexicalPtr
int yylex( void ) ;
extern int MAC_flex_debug ;

#include <stack>
using std::stack ;

static stack<bool> main_stack ;
static stack<int> nb_line_stack ;
static stack<std::string> path_stack ;
static stack<std::string> name_stack ;
static stack<istream*> file_stack ;

bool endFile( void ) ;
std::string comment( void ) ;
void MACerror( const char * ) ;
bool readFile( MAC_Module * top,
               istream* input_stream,
               std::string const& name,
               bool main_read ) ;
void switch_to_buffer( istream * file ) ;
void un_switch_to_buffer( void ) ;
void substitute_for_assignment( std::string const& key,
                                MAC_Data* data ) ;

static int MAC__NbLines = 1 ;
static std::string buff ;
static std::string relativePath ;
static std::string currentFile = "" ;
static bool parsing = false ;

static MAC_Module * YY_top_module = 0 ;
static MAC_Module * dummy_module = 0 ;
static MAC_List * modules_LILO = 0 ;

static MAC_Context const* CTX_INI = 0 ;

#define CREATE_OP(op,name,lst)\
             stringVector const& exps = \
                MAC_Expression::registered_expressions() ; \
             if( !exps.has( name ) ) \
             { std::string mess = "Unknown operator \"" ; \
               mess += name ; \
               mess += "\"\n" ; \
               mess += "valid ones are :\n" ; \
               stringVector e = exps ; \
               e.sort() ; \
               for( size_t i=0 ; i<e.size() ; ++i ) \
                  mess += "   - \""+e(i)+"\"\n" ; \
               MACerror( mess.c_str()  ) ; } \
             if( !MAC_Expression::valid_arguments_of( name, lst ) ) \
             { std::string mess =  "Valid syntax for \"" ; \
               mess += name ; \
               mess += "\" operator is : \n  " ; \
               mess += MAC_Expression::usage_of( name ) ; \
               MACerror( mess.c_str() ) ; } \
             MAC_Expression * op = MAC_Expression::create( 0, name, lst, comment() ) ; \
             op->unset_external_brackets() ;

#define CREATE_SINGLE_OP(name,arg1,resu) \
             MAC_List * lst = MAC_List::create( 0 ) ; \
             lst->append( arg1->to_data() ) ; \
             arg1->change_owner(lst,arg1->to_data() ) ; \
             CREATE_OP( op, name, lst ) ; \
             lst->set_owner( op ) ; \
             resu=MAC_Lexical::create( op )

#define CREATE_BIN_OP(name,arg1,arg2,resu) \
             MAC_List * lst = MAC_List::create( 0 ) ; \
             lst->append( arg1->to_data() ) ; \
             lst->append( arg2->to_data() ) ; \
             arg1->change_owner(lst,arg1->to_data() ) ; \
             arg2->change_owner(lst,arg2->to_data() ) ; \
             CREATE_OP( op, name, lst ) ; \
             lst->set_owner( op ) ; \
             resu=MAC_Lexical::create( op )

%}
/* Grammar description */
%token MAC__IDENTIF MAC__STRING MAC__REAL MAC__INTEGER MAC__EOF
/* Grammar key-words */
%token MAC__ZERO
%token MAC__MODULE 
%token MAC__END 
%token MAC__TRUE MAC__FALSE
%token MAC_INCLUDE
%token MAC__CONCAT
%token MAC__OR
%token MAC__AND
%token MAC__IF
%token MAC__LE
%token MAC__GE
%token MAC__NEQ
%token MAC__LAST
/* Special grammar key-words */
%start data_file
%left MAC__CONCAT
%left MAC__OR MAC__AND
%left '=' MAC__NEQ
%left '<' '>' MAC__LE MAC__GE
%left '+' '-'
%left '*' '/'
%nonassoc UMINUS
%nonassoc UNOT
%%
   /**********************************************************************
    * Top components
    */
      
data_file : /* empty */
           | data_file item_data_file 

item_data_file : module_def
               | directive

directive: '#' MAC_INCLUDE path
               {
                  if( $3->to_data()->data_type() != MAC_Data::String )
                       MACerror( "Include must refer to string expression" ) ;
                  if( YY_top_module!=dummy_module )
                  {
                    if( ! $3->to_data()->value_can_be_evaluated(
                                               YY_top_module->context() ) )
                    {
                       std::string msg = "Include path can not be evaluated\n"  ;
                       msg += "Undefined variable(s):\n" ;
                       stringVector const& undef =
                          $3->to_data()->undefined_variables(
                                                  YY_top_module->context() ) ;
                       for( size_t i=0 ; i<undef.size() ; ++i )
                       {
                          msg += "   - \""+undef(i)+"\"\n" ;
                       }
                       MACerror( msg.c_str() ) ;
                    }
                    std::string file_data =
                       $3->to_data()->to_string( YY_top_module->context() ) ;
                    readFile( YY_top_module, 0, file_data, false ) ;
                  }
               }
          | MAC__EOF { if( endFile() ) YYACCEPT ; }

path : '(' something ')' { $$ = $2 ; }
     | MAC__STRING
     
   /**********************************************************************
    * Utilities specifications
    */
    
   /*
    * Free list specifications
    */
    
free_list : /* Empty list */ 
        | free_list subsitute
        | free_list assignment
        | free_list variable_def
        | free_list module_def
        | free_list directive

assignment : MAC__IDENTIF '=' something
             {
                if( YY_top_module!=dummy_module )
                {
                   if( $1->to_data()->data_type()!=MAC_Data::String )
                      MACerror( "Key must be a string" ) ;
                   std::string const& str = $1->to_data()->to_string() ;
                   if( YY_top_module->has_entry( str ) )
                   {
                     std::string mess = "Module " ;
                     mess+=YY_top_module->name() + " already contains "+
                        str + " assignment (use ==)" ;
                     MACerror( mess.c_str() ) ;
                   }
                   MAC_Data* data = $3->to_data() ;
                   $3->change_owner(YY_top_module, data) ;
                   if( YY_top_module->has_module( str ) )
                   {
                      std::string mess = "Can't create entry and module with the same name "+str+" in module " + YY_top_module->name() ;
                      MACerror( mess.c_str() ) ;
                   }
                   YY_top_module->add_entry( str, data ) ;
                }
                $$ = 0 ;
             }

variable_def : variable '=' something
             {
                if( YY_top_module!=dummy_module )
                {
                   MAC_Variable* var = dynamic_cast<MAC_Variable*>(
                      $1->to_data() ) ;
                   $1->change_owner(YY_top_module,var) ;
                   MAC_ASSERT( var!=0 ) ;
                   if( var->data_type()!=$3->to_data()->data_type() )
                   {
                      std::string mess = "Bad type for expression assigned to " ;
                      mess+=var->name() ;
                      MACerror( mess.c_str() ) ;
                   }
                   if( !YY_top_module->context()->has_variable(var) )
                   {
                      MAC_Data* data = $3->to_data() ;
                      $3->change_owner(0,data) ;
                      YY_top_module->add_variable( var, data ) ;
                   }
                   else if( !CTX_INI->has_variable(var) )
                   {
                      std::string mess = "Module " ;
                      mess+=YY_top_module->name() + " already contains "+
                            var->name() + " variable (use ==)" ;
                      MACerror( mess.c_str() ) ;
                   }
// Module initial context is not modified...
//                   else
//                   {
//                      MAC_Data* data = $3->to_data() ;
//                      $3->change_owner(0,data) ;
//                      YY_top_module->modify_variable( var, data ) ;
//                   }
                }
                $$ = 0 ;
             }
variable_def : variable '=''=' something
             {
                if( YY_top_module!=dummy_module )
                {
                   MAC_Variable* var = dynamic_cast<MAC_Variable*>(
                      $1->to_data() ) ;
                   MAC_ASSERT( var!=0 ) ;
                   if( var->data_type()!=$4->to_data()->data_type() )
                   {
                      std::string mess = "Bad type for expression assigned to " ;
                      mess+=var->name() ;
                      MACerror( mess.c_str() ) ;
                   }
                   if( !YY_top_module->context()->has_variable(var) )
                   {
                      std::string mess = var->name() + " doesn't already exist " ;
                      MACerror( mess.c_str() ) ;
                   }
                   MAC_Data* data = $4->to_data() ;
                   $4->change_owner(0,data) ;
                   YY_top_module->modify_variable( var, data ) ;
                }
                $$ = 0 ;
             }

subsitute : MAC__IDENTIF '=''=' something
             {
                if( YY_top_module!=dummy_module )
                {
                   if( $1->to_data()->data_type()!=MAC_Data::String )
                      MACerror( "Key must be a string" ) ;
                   MAC_String * str = static_cast<MAC_String*>( $1->to_data() ) ;
                   if( !YY_top_module->has_entry( str->to_string() ) )
                   {
                      std::string mess = str->to_string() + " doesn't already exist " ;
                      MACerror( mess.c_str() ) ;
                   }
                   substitute_for_assignment( str->to_string(),
                                              $4->to_data() ) ;
                }
                $$ = 0 ;
             }

something : SimpleType 
     | vector
     | array
     | function
     | variable
     | test
     | binary_operator
     | unary_operator
     |'(' something ')'
          {
             $$ = $2 ;
             if( $$->is_data() )
             {
                MAC_Expression* op = dynamic_cast<MAC_Expression*>( $$->to_data() ) ;
                if( op != 0 ) op->set_external_brackets() ;
             }
          }
   

function: MAC__IDENTIF '(' liste_args ')' 
          {
             std::string const& op_name = $1->to_data()->to_string() ;
             if( op_name=="this_file_dir" )
             {
                $$=MAC_Lexical::create(
                   MAC_String::create( 0, relativePath ) ) ;
             }
             else
             {
                CREATE_OP( op, op_name, $3->to_list() ) ; 
                $$=MAC_Lexical::create( op ) ;
                $3->change_owner(op,$3->to_list()) ;
             }
          }

variable: '$' MAC__IDENTIF 
          {
             std::string const& str = $2->to_data()->to_string() ;
             if( str.length()<2 ||
                 ( str[0]!='I' && str[0]!='D' && str[0]!='B' && str[0]!='S' )
                 || ( str[1]!='S' && str[1]!='V' && str[1]!='A' ) )
             {
                std::string msg = "\""+str+"\" is not a valid variable name\n" ;
                msg += "A valid name is \"XY_name\"\n" ;
                msg += "   where \"X\" is the scalar type of the variable :\n" ;
                msg += "       - \"I\" : integer\n" ;
                msg += "       - \"D\" : double\n" ;
                msg += "       - \"B\" : boolean\n" ;
                msg += "       - \"S\" : string\n" ;
                msg += "   and \"Y\" defined its dimension :\n" ;
                msg += "       - \"S\" : simple (only one element)\n" ;
                msg += "       - \"V\" : vector\n" ;
                msg += "       - \"A\" : array2D\n" ;
                msg += "Examples : \"DV_coordinates\", \"SS_name\", \"IA_connectivity\"\n" ;
                MACerror( msg.c_str() ) ;
             }
             MAC_Variable const* var = MAC_Variable::object(
                $2->to_data()->to_string() ) ;
             $$=MAC_Lexical::create( var->create_clone(0) ) ;
          }

test: '(' switches_list something ')' 
          {
             MAC_List * lst = $2->to_list() ;
             $3->change_owner(lst,$3->to_data() ) ;
             lst->append( $3->to_data() ) ;
             CREATE_OP( op, "(?:)",lst ) ;
             op->set_external_brackets() ;
             $2->change_owner(op,lst) ;
             $$=MAC_Lexical::create( op ) ;
          }

switches_list :  something '?' something ':'
              {
                 MAC_LABEL( "Gram.y :: switches_list1" ) ;
                 MAC_List * lst = MAC_List::create( 0 ) ;
                 $1->change_owner(lst,$1->to_data() ) ;
                 lst->append( $1->to_data() ) ;
                 $3->change_owner(lst,$3->to_data() ) ;
                 lst->append( $3->to_data() ) ;
                 $$ = MAC_Lexical::create( lst ) ;
              }
  | switches_list something '?' something ':'
              {
                 MAC_LABEL( "Gram.y :: switches_list2" ) ;
                 MAC_List * lst = $1->to_list() ;
                 $2->change_owner(lst,$2->to_data() ) ;
                 lst->append( $2->to_data() ) ;
                 $4->change_owner(lst,$4->to_data() ) ;
                 lst->append( $4->to_data() ) ;
                 $$ = $1 ;
              }

liste_args:   {  $$=MAC_Lexical::create( MAC_List::create( 0 ) ) ; }
  | something
              {
                 MAC_List * lst = MAC_List::create( 0 ) ;
                 $1->change_owner(lst,$1->to_data() ) ;
                 lst->append( $1->to_data() ) ;
                 $$ = MAC_Lexical::create( lst ) ;
              }
  | liste_args ',' something
              {
                 MAC_List * lst = $1->to_list() ;
                 $3->change_owner(lst,$3->to_data() ) ;
                 lst->append( $3->to_data() ) ;
                 $$ = $1 ;
              }

unary_operator: '-' something %prec UMINUS{
                    if( $2->to_data()->is_raw_data() )
                    {
                       MAC_Data::Type dt = $2->to_data()->data_type() ;
                       
                       if( dt == MAC_Data::Double )
                       {
                          $$=MAC_Lexical::create(
                             MAC_Double::create( 0, - $2->to_data()->to_double() ) ) ;
                       }
                       else if( dt == MAC_Data::Int )
                       {
                          $$=MAC_Lexical::create(
                             MAC_Int::create( 0, - $2->to_data()->to_int() ) ) ;
                       }
                       else
                       {
                          MACerror( "Undefined unary minus operator" ) ;
                       }
                    }
                    else
                    {
                       CREATE_SINGLE_OP("unary_minus",$2,$$) ; 
                    } }
              | '!' something %prec UNOT { CREATE_SINGLE_OP("!",$2,$$) ; }

binary_operator: something '<' something { CREATE_BIN_OP("<",$1,$3,$$) ; }
               | something '>' something { CREATE_BIN_OP(">",$1,$3,$$) ; }
               | something '+' something { CREATE_BIN_OP("+",$1,$3,$$) ; }
               | something '-' something { CREATE_BIN_OP("-",$1,$3,$$) ; }
               | something '*' something { CREATE_BIN_OP("*",$1,$3,$$) ; }
               | something '=' something { CREATE_BIN_OP("=",$1,$3,$$) ; }
               | something '/' something { CREATE_BIN_OP("/",$1,$3,$$) ; }
               | something MAC__NEQ something { CREATE_BIN_OP("!=",$1,$3,$$) ; }
               | something MAC__CONCAT something { CREATE_BIN_OP("<<",$1,$3,$$) ; }
               | something MAC__OR  something { CREATE_BIN_OP("||",$1,$3,$$) ; }
               | something MAC__AND something { CREATE_BIN_OP("&&",$1,$3,$$) ; }
               | something MAC__LE something { CREATE_BIN_OP("<=",$1,$3,$$) ; }
               | something MAC__GE something { CREATE_BIN_OP(">=",$1,$3,$$) ; }

module_def : module
         |  error

module : module_deb free_list module_fin
         {
            MAC_ASSERT( $3==0 || YY_top_module!=$3->to_module() ) ;
            $$ = $3 ;
         }

if_module : { $$ = MAC_Lexical::create( MAC_Bool::create( 0, true ) ) ; }
          |  MAC__IF '(' something ')'
         {
            if( $3->to_data()->data_type() != MAC_Data::Bool )
               MACerror( "Conditional module if must refer to boolean expression" ) ;
            bool cond = false ;
            if( YY_top_module!=dummy_module )
            {
               if( ! $3->to_data()->value_can_be_evaluated(
                                             YY_top_module->context() ) )
               {
                  std::string msg = "Conditional module if can not be evaluated\n" ;
                  msg += "Undefined variable(s):\n" ;
                  stringVector const& undef =
                     $3->to_data()->undefined_variables(
                                                  YY_top_module->context() ) ;
                  for( size_t i=0 ; i<undef.size() ; ++i )
                  {
                     msg += "   - \""+undef(i)+"\"\n" ;
                  }
                  MACerror( msg.c_str() ) ;
               }
               cond = $3->to_data()->to_bool( YY_top_module->context() ) ;
            }
            $$ = MAC_Lexical::create( MAC_Bool::create( 0, cond ) ) ;
         }

module_deb : if_module MAC__MODULE MAC__IDENTIF
       {
          modules_LILO->append( YY_top_module ) ;
          if( $1->to_data()->to_bool() && YY_top_module!=dummy_module )
          {
             MAC_Module* old = YY_top_module ;
             std::string a_name = $3->to_data()->to_string() ;
             std::string path = "/" + a_name ;
             if( old->has_module( path ) )
             {
                YY_top_module = old->module( path ) ;
             }
             else
             {
                YY_top_module = MAC_Module::create( old, a_name ) ;
                if( old->has_entry( a_name ) )
                {
                   std::string mess = "Can't create entry and module with the same name "+a_name+" in module " + old->name() ;
                   MACerror( mess.c_str() ) ;
                }
                old->add_module( YY_top_module ) ;
             }
          }
          else
          {
             YY_top_module = dummy_module ;
          }
       }


module_fin : MAC__END MAC__MODULE  MAC__IDENTIF
       {
          if( YY_top_module!=dummy_module )
          {
             if( $3->to_data()->to_string()!=YY_top_module->name() )
             {
                std::string mess = "Module " + YY_top_module->name() +
                   " can't be closed with " + $3->to_data()->to_string() ;
                MACerror( mess.c_str() ) ;
             }
             $$=MAC_Lexical::create( YY_top_module ) ;
          }
          else
          {
             $$=0 ;
          }
          size_t last = modules_LILO->count()-1 ;
          YY_top_module = static_cast<MAC_Module*>( modules_LILO->at( last ) ) ;
          modules_LILO->remove_at( last ) ;
          if( modules_LILO->has( dummy_module ) )
          {
             YY_top_module = dummy_module ;
          }
       }  
    
   /*
    * Vector specifications
    */
vector : '<' simple_item_list '>'
         {
            MAC_List * lst = $2->to_list() ;
            MAC_Data* res=0 ;
            if( lst->count() > 0  )
            {
               MAC_Data const* item =
                  dynamic_cast<MAC_Data const*>(lst->at( 0 )) ;
               if( item!=0  )
               {
                  switch( item->data_type() )
                  {
                     case MAC_Data::Double : 
                        res = MAC_DoubleVector::create( 0, lst ) ;
                        break ;
                     case MAC_Data::Int : 
                        res = MAC_IntVector::create( 0, lst ) ;
                        break ;
                     case MAC_Data::Bool : 
                        res = MAC_BoolVector::create( 0, lst ) ;
                        break ;
                     case MAC_Data::String : 
                        res = MAC_StringVector::create( 0, lst ) ;
                        break ;
                     default :
                        break ;
                  }
               }
            }
            if( res==0 )
            {
		MACerror( "invalid list of values enclosed in < .. >" ) ;
            }
            $$ = MAC_Lexical::create( res ) ;
         }

simple_item_list : { $$ = MAC_Lexical::create( MAC_List::create( 0 ) ) ; }
	| simple_item_list item_vector
          {
             $1->to_list()->append( $2->to_data() ) ;
             $$=$1;
          }

array : '[' simple_vector_list ']'
{
            MAC_List * lst = $2->to_list() ;
            MAC_Data* res=0 ;
            if( lst->count() > 0  )
            {
               MAC_Data const* item =
                  dynamic_cast<MAC_Data const*>(lst->at( 0 )) ;
               if( item!=0  )
               {
                  switch( item->data_type() )
                  {
                     case MAC_Data::DoubleVector : 
                        res = MAC_DoubleArray2D::create( 0, lst ) ;
                        break ;
                     case MAC_Data::IntVector : 
                        res = MAC_IntArray2D::create( 0, lst ) ;
                        break ;
                     case MAC_Data::BoolVector : 
                        res = MAC_BoolArray2D::create( 0, lst ) ;
                        break ;
                     case MAC_Data::StringVector : 
                        res = MAC_StringArray2D::create( 0, lst ) ;
                        break ;
                     default :
                        break ;
                  }
               }
            }
            if( res==0 )
            {
		MACerror( "invalid list of vectors enclosed in [ .. ]" ) ;
            }
            $$ = MAC_Lexical::create( res ) ;
}

simple_vector_list : vector
          { $$ = MAC_Lexical::create( MAC_List::create( 0 ) ) ;
            $$->to_list()->append($1->to_data() ) ;
          }
	| simple_vector_list ',' vector
          {
             $1->to_list()->append( $3->to_data() ) ;
             $$=$1;
          }

item_vector: unary_operator | SimpleType
/*
 * Scalar types
 */

SimpleType: MAC__REAL
          | MAC__INTEGER
          | MAC__STRING
          | MAC__TRUE      { $$ = MAC_Lexical::create( MAC_Bool::create( 0, true ) ) ; }
	  | MAC__FALSE     { $$ = MAC_Lexical::create( MAC_Bool::create( 0, false ) ) ; }      

%%
// Read a file        
// read recursivly data file throught yyparse function
//----------------------------------------------------------------------
bool MAC_readFile( MAC_Module * top,
                   istream* input_stream,
                   std::string const& name,
                   bool debug )
//----------------------------------------------------------------------
{
   MAC_LABEL( "Gram.y::MAC_readFile" ) ;
   MAC_ASSERT( EQUIVALENT( input_stream!=0, name.empty() ) ) ;
   MAC_ASSERT( top!=0 ) ;
   MAC_ASSERT( IMPLIES( input_stream!=0, input_stream->good() ) ) ;

   if( parsing ) // Direct recursive called is forbidden
   {
      MACerror( "Try to read a new file, while a data file is already being parsed." ) ;
   }
   MAC_flex_debug = ( debug ? 1 : 0 ) ;
   CTX_INI = top->context() ;
   bool result = readFile( top, input_stream, name, true ) ;
   CTX_INI = 0 ;

   return( result ) ;
}

//----------------------------------------------------------------------
bool readFile( MAC_Module * top,
               istream* input_stream,
               std::string const& name,
               bool main_read )
//----------------------------------------------------------------------
{
   MAC_LABEL( "Gram.y::readFile" ) ;
   MAC_ASSERT( EQUIVALENT( input_stream!=0, name.empty() ) ) ;
   MAC_ASSERT( top!=0 ) ;
   MAC_ASSERT( IMPLIES( input_stream!=0, input_stream->good() ) ) ;
   
   YY_top_module = top ;
   if( modules_LILO==0 ) modules_LILO = MAC_ListIdentity::create( 0 ) ;
   if( dummy_module==0 ) dummy_module = MAC_Module::create( 0, "dummy_module" ) ;
   if( main_read )
   {
      parsing = true ;
      relativePath = "." ;
   }
   
   main_stack.push( main_read ) ;
   path_stack.push( relativePath ) ;
   nb_line_stack.push( MAC__NbLines ) ;
   name_stack.push( currentFile ) ;
   MAC__NbLines = 1 ;
   
   istream* newFile = 0 ;
   std::ifstream* createdStream = 0 ;
   if( input_stream!=0 )
   {
      newFile = input_stream ;
      currentFile = "stream reading" ;
   }
   else if( name=="-stdin" )
   {
      newFile = &std::cin ;
      currentFile = "standard input stream" ;
   }
   else
   {
      char separator = MAC_System::path_name_separator() ;
      
      bool absolute_dir = name.find_first_of( separator )==0 ||
         // Special case for windows-like path name
         ( name.length() > 2 && name[1]==':' && name[2]=='\\' ) ;
      size_t id = name.find_last_of( separator ) ;
      if( absolute_dir )
      {
         currentFile = name ;
         relativePath = name.substr( 0, id ) ;
      }
      else
      {
         currentFile = relativePath + separator + name ;
         if( id < name.length() )
         {
            relativePath = relativePath + separator + name.substr( 0, id ) ;
         }
      }
      newFile = createdStream = new std::ifstream( currentFile.c_str() ) ;
      if( newFile->fail() )
         MAC_Error::object()->raise_plain( "Unable to open file "+currentFile ) ;
      
   }
   
   file_stack.push(newFile) ;
   switch_to_buffer( newFile ) ;
   
   if( main_read )
   {
      MACparse() ;
      parsing = false ;
   }

   if( createdStream!=0 && main_read ) // Other streams are destroyed in endFile
   {
      createdStream->close() ;
      delete createdStream ; createdStream=0 ;
   }
   
   return true ;
   
}

//----------------------------------------------------------------------
bool endFile( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "Gram.y::endFile" ) ;
   MAC_ASSERT( !main_stack.empty() ) ;
   MAC_ASSERT( !path_stack.empty() ) ;
   MAC_ASSERT( !nb_line_stack.empty() ) ;
   MAC_ASSERT( !name_stack.empty() ) ;
   MAC_ASSERT( !file_stack.empty() ) ;
   
   un_switch_to_buffer() ;
   bool main_read = main_stack.top() ; main_stack.pop() ;
   relativePath = path_stack.top() ; path_stack.pop() ;
   MAC__NbLines = nb_line_stack.top() ; nb_line_stack.pop() ;
   currentFile = name_stack.top() ; name_stack.pop() ;
   istream* file = file_stack.top() ; file_stack.pop() ;
   
   bool res = false ;
   
   if( !main_read )
   { 
      delete file ; file = 0 ;
      res = false ;
   }
   else
   {
      MAC_ASSERT( main_stack.empty() ) ;
      MAC_ASSERT( path_stack.empty() ) ;
      MAC_ASSERT( nb_line_stack.empty() ) ;
      MAC_ASSERT( name_stack.empty() ) ;
      MAC_ASSERT( file_stack.empty() ) ;
      if( modules_LILO->count()!=0 )
      {
         std::string mess = "When reading MAC data structure, following modules are not correclty closed\n" ;
         for( size_t i=0 ; i<modules_LILO->count() ; i++ )
            mess += "\"" + static_cast<MAC_Module*>( modules_LILO->at( i ) )->name() + "\"\n" ;
         MACerror( mess.c_str() ) ;
      }
      
      MAC_Lexical::remove_all_lexical() ;
      modules_LILO->destroy() ; modules_LILO=0 ;
      dummy_module->destroy() ; dummy_module = 0 ;
      res = true ;
   }
   
   return res ;
   
}

//----------------------------------------------------------------------
void MAC_re_init_parser( void )
//----------------------------------------------------------------------
{
   MAC_LABEL( "MAC_re_init_parser" ) ;
   
   while(!main_stack.empty()) main_stack.pop() ;
   while(!path_stack.empty()) path_stack.pop() ;
   while(!nb_line_stack.empty()) nb_line_stack.pop() ;
   while(!name_stack.empty()) name_stack.pop() ;
   while(!file_stack.empty()) {
      //if( file_stack.top()!=0 ) delete file_stack.top() ;
      file_stack.pop() ; 
   }
   
   MAC_re_init_lexer() ;
   MAC_Lexical::remove_all_lexical() ;
   if(modules_LILO!=0 ) {
      modules_LILO->destroy() ;
      modules_LILO=0 ;
   }
   if( dummy_module!=0 ) 
   {      
      dummy_module->destroy() ;
      dummy_module = 0 ;
   }
   parsing = false ;
   
}

//
// Buffering used to help in error case 
//
//----------------------------------------------------------------------
void MAC__Buffer( std::string const& chain ) 
//----------------------------------------------------------------------
{
   size_t cr=chain.find('\n') ;
   size_t cr_last = 0 ;
   for( ;
        cr<chain.length() ;
        cr=chain.find('\n', cr+1 ) )
   {
      MAC__NbLines++ ;
      cr_last = cr+1 ;
      buff="" ;
   }
   buff += chain.substr( cr_last, chain.length()-cr_last ) ;
}

//----------------------------------------------------------------------
void MACerror( const char * s)
//----------------------------------------------------------------------
{
   MAC_re_init_parser() ;
   
   MAC_Error::object()->raise_read_syntax_error(
      currentFile,
      MAC__NbLines,
      buff,
      s ) ;
}

//----------------------------------------------------------------------
std::string comment( void )
//----------------------------------------------------------------------
{
   std::ostringstream msg ;
   msg << std::endl
       << "Last line number " << MAC__NbLines
       << " read in file " << currentFile << " was:" << std::endl
       << ">> " ;
   if( buff.length()<40 )
   {
      msg << buff << " <<" << std::endl ;
   }
   else
   {
      msg << buff.substr( 0, 40 ) << " ... (TO BE CONTINUED) <<" << std::endl ;
   }
   
   return msg.str() ;
}

//----------------------------------------------------------------------
std::string MAC_current_parsed_module_path_name( void )
//----------------------------------------------------------------------
{
   std::string result = "" ;
   if( YY_top_module!=dummy_module && YY_top_module!=0 )
   {      
      result = YY_top_module->absolute_path_name() ;
   }
   return result ;
}

//----------------------------------------------------------------------
int MAC_current_parsed_line( void )
//----------------------------------------------------------------------
{
   return MAC__NbLines ;
}

//----------------------------------------------------------------------
void substitute_for_assignment( std::string const& key,
                                MAC_Data* data )
//----------------------------------------------------------------------
{
   MAC_CHECK( YY_top_module!=dummy_module ) ;
   
   if( modules_LILO->index_limit()<1 )
   {
      std::string mess = "No path to apply substitution of " ;
      mess += key ;
      MACerror( mess.c_str() ) ;
   }
   
   MAC_Module * mod = static_cast<MAC_Module *>( modules_LILO->at( 0 ) ) ;
   MAC_Iterator * it = modules_LILO->create_iterator( 0 ) ;
   it->go_next() ;
   
   for( ; it->is_valid() ; it->go_next() )
   {
      MAC_Module * child = static_cast<MAC_Module *>( it->item() ) ;
      if( !mod->has_module( child->name() ) )
      {
         std::string mess = "No path to apply substitution of " ;
         mess += key + " for module " + child->name() ;
         mod->print( MAC::out(), 0 ) ;
         MACerror( mess.c_str() ) ;
      }
      mod = mod->module( child->name() ) ;
   }
   if( !mod->has_module( YY_top_module->name() ) )
   {
      std::string mess = "No path to apply substitution of " ;
      mess += key + " for module " + YY_top_module->name() ;
      mod->print( MAC::out(), 0 ) ;
      
      MACerror( mess.c_str() ) ;
   }
   mod = mod->module( YY_top_module->name() ) ;
   std::string root_key = "/" +key ;
   if( !mod->has_entry( root_key ) )
   {
      std::string mess = key + " doesn't already exist " ;
      MACerror( mess.c_str() ) ;
   }
   if( data->owner()!=0 ) data = data->create_clone(0) ;
   
   data->set_owner(mod) ;
   mod->replace_data_of_entry( root_key, data ) ;
   
   it->destroy() ;
      
}
