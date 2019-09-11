%option noyywrap
%{
 // File : Gram.l
 // LCH 21.07.1999
 // Object : Lexical analysis description
 //  The valid items are :
 //     MAC__INTEGER , MAC__REAL ,  MAC__STRING
#include <iostream>
using std::istream ;
#include <cstdio>
#include <cstring>

#include <MAC_Lexical.hh>
#define YYSTYPE MAC_LexicalPtr

#include <MAC_Lexical.hh>
#include <MAC_Map.hh>
#include <MAC_Module.hh>
#include <MAC_assertions.hh>
#include <MAC_Root.hh>
#include <MAC_String.hh>
#include <MAC_Double.hh>
#include <MAC_Int.hh>
#include <MAC_Yacc.hh>

#include <stack>
using std::stack ;
#ifdef _WIN32

/* Windows doesn't define the isatty() function.
   As we always read from a file, we can easily
   provide the following simple implementation: */

int isatty (int fd)
{
    return 0;
}
#endif

extern YYSTYPE MAClval ;

/* Forward declarations */
int Identifiers(char * chain) ;
extern void MAC__Buffer(std::string const& chain) ;
extern void MACerror( const char * s) ;
extern MAC_Module* YY_top_module ;

/* Change lex input method to read from C++ streams */
std::istream*my_stream = 0 ;

#define YY_INPUT(buf,result,max_size)\
MAC_ASSERT( my_stream!=0 ) ;\
if ( *my_stream )  { my_stream->read(buf,max_size);result=my_stream->gcount();} else {result=0;}

%}
%%
<<EOF>>                 { return(MAC__EOF) ;}
\/\/.*                  { // Comment description 
                          MAC__Buffer(&yytext[0]) ; 
                        }

\"[^\"]*\"              { // String description 
                          MAC__Buffer(&yytext[0]) ;
                                  // We delete extra ""
                          yytext[ strlen(yytext)-1 ]='\0' ;
                          std::string str((char*) &yytext[1]) ;
                          size_t idx ;
                          while( (idx=str.find("'"))<=str.length() ) str.replace(idx,1,"\"") ;
                          MAClval = MAC_Lexical::create( MAC_String::create( 0, str ) );
                          return(MAC__STRING) ;
                        }

[0-9]+\.[0-9]* |
[0-9]+\.[0-9]*[eE][+-]?[0-9]+      { 
                          MAC__Buffer(&yytext[0]) ;
                          MAClval = MAC_Lexical::create( MAC_Double::create( 0,
									    atof(&yytext[0]) ) );
                          return( MAC__REAL ) ;
                        }

[0-9]+                  { // Integer description
                          MAC__Buffer(&yytext[0]) ;
                          MAClval = MAC_Lexical::create( MAC_Int::create( 0, 
                                                                         atoi(&yytext[0]) ) ) ;
                          return(MAC__INTEGER) ;
                        }

[a-zA-Z_][a-zA-Z_0-9\#]* |
"<<" | 
"||" | 
"&&" | 
"<=" | 
">=" | 
"!="                     { // Identifier description 
                          MAC__Buffer(&yytext[0]) ;
                          return( Identifiers(&yytext[0]) );
                        }

[\t\n \r]               { // Ignore tabular, carriage return white space, lf
                          MAC__Buffer(&yytext[0]) ; 
                        }

[\<\>\[\]{}#+(),?\:=\-\*/$!]   { // Ponctuation
                          MAC__Buffer(&yytext[0]) ; 
                          return( yytext[0] ) ;
                        }

.                       { // Ignore extra-character after them been signified
                          MAC__Buffer(&yytext[0]) ; 
                          std::cerr << "unknown character : " << yytext[0] << std::endl ;
                          MACerror("parse error") ;
                        }
%%

  //
  // Lexical equivalences of key-words 
  //
  
#define MAC__NB_IDENTIF (MAC__LAST-MAC__ZERO-1)

  static char MAC__EQUIVALENCES[MAC__NB_IDENTIF][20] = {
  	"MODULE", "END" , "true" , "false", "include",
         "<<", "||", "&&", "if", "<=", ">=", "!=" } ;

  int Identifiers(char * chain) 
  {
     int ret = MAC__IDENTIF ;
     for( int i=0 ; i<MAC__NB_IDENTIF ; i++ )
     {
       if(strcmp(MAC__EQUIVALENCES[i],chain)==0)
       {
          ret = MAC__ZERO+1+i ;
	  MAClval = 0 ;
	  break ;
       }
     }
     if( ret == MAC__IDENTIF )
     {
        MAC_String* chainStr = MAC_String::create( 0, chain ) ;
	MAClval = MAC_Lexical::create( chainStr ) ;
     }
     return( ret ) ;
  }
  	
static stack<istream*> file_stack ;
static stack<YY_BUFFER_STATE> buffer_stack ;

void switch_to_buffer( istream * file )
{ 
   buffer_stack.push( YY_CURRENT_BUFFER ) ;
   file_stack.push( my_stream ) ;
   my_stream = file ;
   yy_switch_to_buffer( yy_create_buffer( 0, YY_BUF_SIZE ) ) ; 
}

void un_switch_to_buffer( void )
{ 	
   YY_BUFFER_STATE old_buffer = buffer_stack.top() ; buffer_stack.pop() ;
   yy_delete_buffer( YY_CURRENT_BUFFER ) ;
   yy_switch_to_buffer( old_buffer ) ;
   my_stream = file_stack.top() ; file_stack.pop() ;
}


void MAC_re_init_lexer( void )
{ 
  while(!file_stack.empty()) file_stack.pop() ;
  while(!buffer_stack.empty()) buffer_stack.pop() ;
  my_stream = 0 ;
}


