/* A Bison parser, made by GNU Bison 2.3.  */

/* Skeleton implementation for Bison's Yacc-like parsers in C

   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004, 2005, 2006
   Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor,
   Boston, MA 02110-1301, USA.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "2.3"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Using locations.  */
#define YYLSP_NEEDED 0

/* Substitute the variable and function names.  */
#define yyparse MACparse
#define yylex   MAClex
#define yyerror MACerror
#define yylval  MAClval
#define yychar  MACchar
#define yydebug MACdebug
#define yynerrs MACnerrs


/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     MAC__IDENTIF = 258,
     MAC__STRING = 259,
     MAC__REAL = 260,
     MAC__INTEGER = 261,
     MAC__EOF = 262,
     MAC__ZERO = 263,
     MAC__MODULE = 264,
     MAC__END = 265,
     MAC__TRUE = 266,
     MAC__FALSE = 267,
     MAC_INCLUDE = 268,
     MAC__CONCAT = 269,
     MAC__OR = 270,
     MAC__AND = 271,
     MAC__IF = 272,
     MAC__LE = 273,
     MAC__GE = 274,
     MAC__NEQ = 275,
     MAC__LAST = 276,
     UMINUS = 277,
     UNOT = 278
   };
#endif
/* Tokens.  */
#define MAC__IDENTIF 258
#define MAC__STRING 259
#define MAC__REAL 260
#define MAC__INTEGER 261
#define MAC__EOF 262
#define MAC__ZERO 263
#define MAC__MODULE 264
#define MAC__END 265
#define MAC__TRUE 266
#define MAC__FALSE 267
#define MAC_INCLUDE 268
#define MAC__CONCAT 269
#define MAC__OR 270
#define MAC__AND 271
#define MAC__IF 272
#define MAC__LE 273
#define MAC__GE 274
#define MAC__NEQ 275
#define MAC__LAST 276
#define UMINUS 277
#define UNOT 278




/* Copy the first part of user declarations.  */
#line 1 "Gram.y"

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



/* Enabling traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif

/* Enabling verbose error messages.  */
#ifdef YYERROR_VERBOSE
# undef YYERROR_VERBOSE
# define YYERROR_VERBOSE 1
#else
# define YYERROR_VERBOSE 0
#endif

/* Enabling the token table.  */
#ifndef YYTOKEN_TABLE
# define YYTOKEN_TABLE 0
#endif

#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef int YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif



/* Copy the second part of user declarations.  */


/* Line 216 of yacc.c.  */
#line 289 "y.tab.c"

#ifdef short
# undef short
#endif

#ifdef YYTYPE_UINT8
typedef YYTYPE_UINT8 yytype_uint8;
#else
typedef unsigned char yytype_uint8;
#endif

#ifdef YYTYPE_INT8
typedef YYTYPE_INT8 yytype_int8;
#elif (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
typedef signed char yytype_int8;
#else
typedef short int yytype_int8;
#endif

#ifdef YYTYPE_UINT16
typedef YYTYPE_UINT16 yytype_uint16;
#else
typedef unsigned short int yytype_uint16;
#endif

#ifdef YYTYPE_INT16
typedef YYTYPE_INT16 yytype_int16;
#else
typedef short int yytype_int16;
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif ! defined YYSIZE_T && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned int
# endif
#endif

#define YYSIZE_MAXIMUM ((YYSIZE_T) -1)

#ifndef YY_
# if YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(msgid) dgettext ("bison-runtime", msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(msgid) msgid
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YYUSE(e) ((void) (e))
#else
# define YYUSE(e) /* empty */
#endif

/* Identity function, used to suppress warnings about constant conditions.  */
#ifndef lint
# define YYID(n) (n)
#else
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static int
YYID (int i)
#else
static int
YYID (i)
    int i;
#endif
{
  return i;
}
#endif

#if ! defined yyoverflow || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#     ifndef _STDLIB_H
#      define _STDLIB_H 1
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's `empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (YYID (0))
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined _STDLIB_H \
       && ! ((defined YYMALLOC || defined malloc) \
	     && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef _STDLIB_H
#    define _STDLIB_H 1
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined _STDLIB_H && (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* ! defined yyoverflow || YYERROR_VERBOSE */


#if (! defined yyoverflow \
     && (! defined __cplusplus \
	 || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yytype_int16 yyss;
  YYSTYPE yyvs;
  };

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (yytype_int16) + sizeof (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

/* Copy COUNT objects from FROM to TO.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(To, From, Count) \
      __builtin_memcpy (To, From, (Count) * sizeof (*(From)))
#  else
#   define YYCOPY(To, From, Count)		\
      do					\
	{					\
	  YYSIZE_T yyi;				\
	  for (yyi = 0; yyi < (Count); yyi++)	\
	    (To)[yyi] = (From)[yyi];		\
	}					\
      while (YYID (0))
#  endif
# endif

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack)					\
    do									\
      {									\
	YYSIZE_T yynewbytes;						\
	YYCOPY (&yyptr->Stack, Stack, yysize);				\
	Stack = &yyptr->Stack;						\
	yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAXIMUM; \
	yyptr += yynewbytes / sizeof (*yyptr);				\
      }									\
    while (YYID (0))

#endif

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  2
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   304

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  41
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  28
/* YYNRULES -- Number of rules.  */
#define YYNRULES  71
/* YYNRULES -- Number of states.  */
#define YYNSTATES  120

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   278

#define YYTRANSLATE(YYX)						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,    38,     2,    31,    34,     2,     2,     2,
      32,    33,    27,    25,    37,    26,     2,    28,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,    36,     2,
      23,    22,    24,    35,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,    39,     2,    40,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    29,    30
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const yytype_uint8 yyprhs[] =
{
       0,     0,     3,     4,     7,     9,    11,    15,    17,    21,
      23,    24,    27,    30,    33,    36,    39,    43,    47,    52,
      57,    59,    61,    63,    65,    67,    69,    71,    73,    77,
      82,    85,    90,    95,   101,   102,   104,   108,   111,   114,
     118,   122,   126,   130,   134,   138,   142,   146,   150,   154,
     158,   162,   166,   168,   170,   174,   175,   180,   184,   188,
     192,   193,   196,   200,   202,   206,   208,   210,   212,   214,
     216,   218
};

/* YYRHS -- A `-1'-separated list of the rules' RHS.  */
static const yytype_int8 yyrhs[] =
{
      42,     0,    -1,    -1,    42,    43,    -1,    58,    -1,    44,
      -1,    31,    13,    45,    -1,     7,    -1,    32,    50,    33,
      -1,     4,    -1,    -1,    46,    49,    -1,    46,    47,    -1,
      46,    48,    -1,    46,    58,    -1,    46,    44,    -1,     3,
      22,    50,    -1,    52,    22,    50,    -1,    52,    22,    22,
      50,    -1,     3,    22,    22,    50,    -1,    68,    -1,    63,
      -1,    65,    -1,    51,    -1,    52,    -1,    53,    -1,    57,
      -1,    56,    -1,    32,    50,    33,    -1,     3,    32,    55,
      33,    -1,    34,     3,    -1,    32,    54,    50,    33,    -1,
      50,    35,    50,    36,    -1,    54,    50,    35,    50,    36,
      -1,    -1,    50,    -1,    55,    37,    50,    -1,    26,    50,
      -1,    38,    50,    -1,    50,    23,    50,    -1,    50,    24,
      50,    -1,    50,    25,    50,    -1,    50,    26,    50,    -1,
      50,    27,    50,    -1,    50,    22,    50,    -1,    50,    28,
      50,    -1,    50,    20,    50,    -1,    50,    14,    50,    -1,
      50,    15,    50,    -1,    50,    16,    50,    -1,    50,    18,
      50,    -1,    50,    19,    50,    -1,    59,    -1,     1,    -1,
      61,    46,    62,    -1,    -1,    17,    32,    50,    33,    -1,
      60,     9,     3,    -1,    10,     9,     3,    -1,    23,    64,
      24,    -1,    -1,    64,    67,    -1,    39,    66,    40,    -1,
      63,    -1,    66,    37,    63,    -1,    56,    -1,    68,    -1,
       5,    -1,     6,    -1,     4,    -1,    11,    -1,    12,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const yytype_uint16 yyrline[] =
{
       0,   160,   160,   161,   163,   164,   166,   191,   193,   194,
     204,   205,   206,   207,   208,   209,   211,   237,   274,   299,
     317,   318,   319,   320,   321,   322,   323,   324,   325,   336,
     352,   378,   389,   399,   410,   411,   418,   426,   450,   452,
     453,   454,   455,   456,   457,   458,   459,   460,   461,   462,
     463,   464,   466,   467,   469,   475,   476,   502,   532,   560,
     596,   597,   603,   639,   643,   649,   649,   654,   655,   656,
     657,   658
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "MAC__IDENTIF", "MAC__STRING",
  "MAC__REAL", "MAC__INTEGER", "MAC__EOF", "MAC__ZERO", "MAC__MODULE",
  "MAC__END", "MAC__TRUE", "MAC__FALSE", "MAC_INCLUDE", "MAC__CONCAT",
  "MAC__OR", "MAC__AND", "MAC__IF", "MAC__LE", "MAC__GE", "MAC__NEQ",
  "MAC__LAST", "'='", "'<'", "'>'", "'+'", "'-'", "'*'", "'/'", "UMINUS",
  "UNOT", "'#'", "'('", "')'", "'$'", "'?'", "':'", "','", "'!'", "'['",
  "']'", "$accept", "data_file", "item_data_file", "directive", "path",
  "free_list", "assignment", "variable_def", "subsitute", "something",
  "function", "variable", "test", "switches_list", "liste_args",
  "unary_operator", "binary_operator", "module_def", "module", "if_module",
  "module_deb", "module_fin", "vector", "simple_item_list", "array",
  "simple_vector_list", "item_vector", "SimpleType", 0
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const yytype_uint16 yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,    61,    60,    62,    43,    45,    42,    47,   277,
     278,    35,    40,    41,    36,    63,    58,    44,    33,    91,
      93
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const yytype_uint8 yyr1[] =
{
       0,    41,    42,    42,    43,    43,    44,    44,    45,    45,
      46,    46,    46,    46,    46,    46,    47,    48,    48,    49,
      50,    50,    50,    50,    50,    50,    50,    50,    50,    51,
      52,    53,    54,    54,    55,    55,    55,    56,    56,    57,
      57,    57,    57,    57,    57,    57,    57,    57,    57,    57,
      57,    57,    58,    58,    59,    60,    60,    61,    62,    63,
      64,    64,    65,    66,    66,    67,    67,    68,    68,    68,
      68,    68
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const yytype_uint8 yyr2[] =
{
       0,     2,     0,     2,     1,     1,     3,     1,     3,     1,
       0,     2,     2,     2,     2,     2,     3,     3,     4,     4,
       1,     1,     1,     1,     1,     1,     1,     1,     3,     4,
       2,     4,     4,     5,     0,     1,     3,     2,     2,     3,
       3,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     1,     1,     3,     0,     4,     3,     3,     3,
       0,     2,     3,     1,     3,     1,     1,     1,     1,     1,
       1,     1
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const yytype_uint8 yydefact[] =
{
       2,     0,     1,    53,     7,     0,     0,     3,     5,     4,
      52,     0,    10,     0,     0,     0,     0,     0,    69,    67,
      68,    70,    71,    60,     0,     0,     0,     0,     0,     0,
      23,    24,    25,    27,    26,    21,    22,    20,     9,     0,
       6,    57,     0,     0,    15,    12,    13,    11,     0,    14,
      54,    34,     0,    37,     0,     0,    30,    38,    63,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,    56,     0,     0,     0,     0,    35,     0,
      59,    65,    61,    66,    28,     0,     0,     0,    62,    47,
      48,    49,    50,    51,    46,    44,    39,    40,    41,    42,
      43,    45,     8,     0,    16,    58,     0,    17,    29,     0,
       0,    31,     0,    64,    19,    18,    36,    32,     0,    33
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int8 yydefgoto[] =
{
      -1,     1,     7,     8,    40,    16,    45,    46,    47,    29,
      30,    31,    32,    55,    79,    33,    34,     9,    10,    11,
      12,    50,    35,    52,    36,    59,    82,    37
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -29
static const yytype_int16 yypact[] =
{
     -29,    49,   -29,   -29,   -29,   -25,     6,   -29,   -29,   -29,
     -29,    11,   -29,   112,     2,    25,     1,    20,   -29,   -29,
     -29,   -29,   -29,   -29,   112,   112,    27,   112,    31,   219,
     -29,   -29,   -29,   -29,   -29,   -29,   -29,   -29,   -29,   112,
     -29,   -29,    33,    48,   -29,   -29,   -29,   -29,    37,   -29,
     -29,   112,   128,   -29,   175,   112,   -29,   -29,   -29,   -28,
     112,   112,   112,   112,   112,   112,   112,   112,   112,   112,
     112,   112,   112,   -29,   235,    64,    57,    88,   251,    -4,
     -29,   -29,   -29,   -29,   -29,   112,   197,    31,   -29,   265,
     276,   276,    46,    46,    -2,    -2,    46,    46,   -14,   -14,
     -29,   -29,   -29,   112,   251,   -29,   112,   251,   -29,   112,
     145,   -29,   112,   -29,   251,   251,   251,   -29,   160,   -29
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int8 yypgoto[] =
{
     -29,   -29,   -29,    47,   -29,   -29,   -29,   -29,   -29,   -24,
     -29,    61,   -29,   -29,   -29,    10,   -29,    62,   -29,   -29,
     -29,   -29,   -23,   -29,   -29,   -29,   -29,    13
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -56
static const yytype_int8 yytable[] =
{
      53,    54,     3,    57,    42,    58,    38,    13,     4,    87,
     -55,    43,    88,    71,    72,    74,    63,    64,     5,    14,
      15,    67,    68,    69,    70,    71,    72,    78,    41,   108,
      56,    86,     6,   109,    39,    26,    89,    90,    91,    92,
      93,    94,    95,    96,    97,    98,    99,   100,   101,     2,
       3,   104,    51,   107,    23,    75,     4,    76,   -55,    77,
     105,   110,    81,    44,   113,    83,     5,    17,    18,    19,
      20,    69,    70,    71,    72,    21,    22,    48,    49,   114,
       6,     0,   115,     0,     0,   116,   103,    23,   118,     0,
      24,    17,    18,    19,    20,     0,    25,     0,    26,    21,
      22,     0,    27,    28,     0,     0,     0,     0,     0,     0,
     106,    23,     0,     0,    24,    17,    18,    19,    20,     0,
      25,     0,    26,    21,    22,     0,    27,    28,     0,     0,
       0,     0,    18,    19,    20,    23,     0,     0,    24,    21,
      22,     0,     0,     0,    25,     0,    26,     0,     0,     0,
      27,    28,    80,     0,    24,     0,     0,     0,     0,    60,
      61,    62,     0,    63,    64,    65,    27,    66,    67,    68,
      69,    70,    71,    72,    60,    61,    62,     0,    63,    64,
      65,   117,    66,    67,    68,    69,    70,    71,    72,    60,
      61,    62,     0,    63,    64,    65,   119,    66,    67,    68,
      69,    70,    71,    72,     0,     0,     0,     0,    84,     0,
      85,    60,    61,    62,     0,    63,    64,    65,     0,    66,
      67,    68,    69,    70,    71,    72,     0,     0,     0,     0,
     111,     0,   112,    60,    61,    62,     0,    63,    64,    65,
       0,    66,    67,    68,    69,    70,    71,    72,     0,    60,
      61,    62,    73,    63,    64,    65,     0,    66,    67,    68,
      69,    70,    71,    72,     0,    60,    61,    62,   102,    63,
      64,    65,     0,    66,    67,    68,    69,    70,    71,    72,
      61,    62,     0,    63,    64,    65,     0,    66,    67,    68,
      69,    70,    71,    72,    63,    64,    65,     0,    66,    67,
      68,    69,    70,    71,    72
};

static const yytype_int8 yycheck[] =
{
      24,    25,     1,    27,     3,    28,     4,    32,     7,    37,
       9,    10,    40,    27,    28,    39,    18,    19,    17,    13,
       9,    23,    24,    25,    26,    27,    28,    51,     3,    33,
       3,    55,    31,    37,    32,    34,    60,    61,    62,    63,
      64,    65,    66,    67,    68,    69,    70,    71,    72,     0,
       1,    75,    32,    77,    23,    22,     7,     9,     9,    22,
       3,    85,    52,    16,    87,    52,    17,     3,     4,     5,
       6,    25,    26,    27,    28,    11,    12,    16,    16,   103,
      31,    -1,   106,    -1,    -1,   109,    22,    23,   112,    -1,
      26,     3,     4,     5,     6,    -1,    32,    -1,    34,    11,
      12,    -1,    38,    39,    -1,    -1,    -1,    -1,    -1,    -1,
      22,    23,    -1,    -1,    26,     3,     4,     5,     6,    -1,
      32,    -1,    34,    11,    12,    -1,    38,    39,    -1,    -1,
      -1,    -1,     4,     5,     6,    23,    -1,    -1,    26,    11,
      12,    -1,    -1,    -1,    32,    -1,    34,    -1,    -1,    -1,
      38,    39,    24,    -1,    26,    -1,    -1,    -1,    -1,    14,
      15,    16,    -1,    18,    19,    20,    38,    22,    23,    24,
      25,    26,    27,    28,    14,    15,    16,    -1,    18,    19,
      20,    36,    22,    23,    24,    25,    26,    27,    28,    14,
      15,    16,    -1,    18,    19,    20,    36,    22,    23,    24,
      25,    26,    27,    28,    -1,    -1,    -1,    -1,    33,    -1,
      35,    14,    15,    16,    -1,    18,    19,    20,    -1,    22,
      23,    24,    25,    26,    27,    28,    -1,    -1,    -1,    -1,
      33,    -1,    35,    14,    15,    16,    -1,    18,    19,    20,
      -1,    22,    23,    24,    25,    26,    27,    28,    -1,    14,
      15,    16,    33,    18,    19,    20,    -1,    22,    23,    24,
      25,    26,    27,    28,    -1,    14,    15,    16,    33,    18,
      19,    20,    -1,    22,    23,    24,    25,    26,    27,    28,
      15,    16,    -1,    18,    19,    20,    -1,    22,    23,    24,
      25,    26,    27,    28,    18,    19,    20,    -1,    22,    23,
      24,    25,    26,    27,    28
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const yytype_uint8 yystos[] =
{
       0,    42,     0,     1,     7,    17,    31,    43,    44,    58,
      59,    60,    61,    32,    13,     9,    46,     3,     4,     5,
       6,    11,    12,    23,    26,    32,    34,    38,    39,    50,
      51,    52,    53,    56,    57,    63,    65,    68,     4,    32,
      45,     3,     3,    10,    44,    47,    48,    49,    52,    58,
      62,    32,    64,    50,    50,    54,     3,    50,    63,    66,
      14,    15,    16,    18,    19,    20,    22,    23,    24,    25,
      26,    27,    28,    33,    50,    22,     9,    22,    50,    55,
      24,    56,    67,    68,    33,    35,    50,    37,    40,    50,
      50,    50,    50,    50,    50,    50,    50,    50,    50,    50,
      50,    50,    33,    22,    50,     3,    22,    50,    33,    37,
      50,    33,    35,    63,    50,    50,    50,    36,    50,    36
};

#define yyerrok		(yyerrstatus = 0)
#define yyclearin	(yychar = YYEMPTY)
#define YYEMPTY		(-2)
#define YYEOF		0

#define YYACCEPT	goto yyacceptlab
#define YYABORT		goto yyabortlab
#define YYERROR		goto yyerrorlab


/* Like YYERROR except do call yyerror.  This remains here temporarily
   to ease the transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  */

#define YYFAIL		goto yyerrlab

#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)					\
do								\
  if (yychar == YYEMPTY && yylen == 1)				\
    {								\
      yychar = (Token);						\
      yylval = (Value);						\
      yytoken = YYTRANSLATE (yychar);				\
      YYPOPSTACK (1);						\
      goto yybackup;						\
    }								\
  else								\
    {								\
      yyerror (YY_("syntax error: cannot back up")); \
      YYERROR;							\
    }								\
while (YYID (0))


#define YYTERROR	1
#define YYERRCODE	256


/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

#define YYRHSLOC(Rhs, K) ((Rhs)[K])
#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)				\
    do									\
      if (YYID (N))                                                    \
	{								\
	  (Current).first_line   = YYRHSLOC (Rhs, 1).first_line;	\
	  (Current).first_column = YYRHSLOC (Rhs, 1).first_column;	\
	  (Current).last_line    = YYRHSLOC (Rhs, N).last_line;		\
	  (Current).last_column  = YYRHSLOC (Rhs, N).last_column;	\
	}								\
      else								\
	{								\
	  (Current).first_line   = (Current).last_line   =		\
	    YYRHSLOC (Rhs, 0).last_line;				\
	  (Current).first_column = (Current).last_column =		\
	    YYRHSLOC (Rhs, 0).last_column;				\
	}								\
    while (YYID (0))
#endif


/* YY_LOCATION_PRINT -- Print the location on the stream.
   This macro was not mandated originally: define only if we know
   we won't break user code: when these are the locations we know.  */

#ifndef YY_LOCATION_PRINT
# if YYLTYPE_IS_TRIVIAL
#  define YY_LOCATION_PRINT(File, Loc)			\
     fprintf (File, "%d.%d-%d.%d",			\
	      (Loc).first_line, (Loc).first_column,	\
	      (Loc).last_line,  (Loc).last_column)
# else
#  define YY_LOCATION_PRINT(File, Loc) ((void) 0)
# endif
#endif


/* YYLEX -- calling `yylex' with the right arguments.  */

#ifdef YYLEX_PARAM
# define YYLEX yylex (YYLEX_PARAM)
#else
# define YYLEX yylex ()
#endif

/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)			\
do {						\
  if (yydebug)					\
    YYFPRINTF Args;				\
} while (YYID (0))

# define YY_SYMBOL_PRINT(Title, Type, Value, Location)			  \
do {									  \
  if (yydebug)								  \
    {									  \
      YYFPRINTF (stderr, "%s ", Title);					  \
      yy_symbol_print (stderr,						  \
		  Type, Value); \
      YYFPRINTF (stderr, "\n");						  \
    }									  \
} while (YYID (0))


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_value_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_value_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  if (!yyvaluep)
    return;
# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# else
  YYUSE (yyoutput);
# endif
  switch (yytype)
    {
      default:
	break;
    }
}


/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_symbol_print (FILE *yyoutput, int yytype, YYSTYPE const * const yyvaluep)
#else
static void
yy_symbol_print (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE const * const yyvaluep;
#endif
{
  if (yytype < YYNTOKENS)
    YYFPRINTF (yyoutput, "token %s (", yytname[yytype]);
  else
    YYFPRINTF (yyoutput, "nterm %s (", yytname[yytype]);

  yy_symbol_value_print (yyoutput, yytype, yyvaluep);
  YYFPRINTF (yyoutput, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_stack_print (yytype_int16 *bottom, yytype_int16 *top)
#else
static void
yy_stack_print (bottom, top)
    yytype_int16 *bottom;
    yytype_int16 *top;
#endif
{
  YYFPRINTF (stderr, "Stack now");
  for (; bottom <= top; ++bottom)
    YYFPRINTF (stderr, " %d", *bottom);
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)				\
do {								\
  if (yydebug)							\
    yy_stack_print ((Bottom), (Top));				\
} while (YYID (0))


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yy_reduce_print (YYSTYPE *yyvsp, int yyrule)
#else
static void
yy_reduce_print (yyvsp, yyrule)
    YYSTYPE *yyvsp;
    int yyrule;
#endif
{
  int yynrhs = yyr2[yyrule];
  int yyi;
  unsigned long int yylno = yyrline[yyrule];
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu):\n",
	     yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      fprintf (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr, yyrhs[yyprhs[yyrule] + yyi],
		       &(yyvsp[(yyi + 1) - (yynrhs)])
		       		       );
      fprintf (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)		\
do {					\
  if (yydebug)				\
    yy_reduce_print (yyvsp, Rule); \
} while (YYID (0))

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
# define YY_SYMBOL_PRINT(Title, Type, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef	YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif



#if YYERROR_VERBOSE

# ifndef yystrlen
#  if defined __GLIBC__ && defined _STRING_H
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static YYSIZE_T
yystrlen (const char *yystr)
#else
static YYSIZE_T
yystrlen (yystr)
    const char *yystr;
#endif
{
  YYSIZE_T yylen;
  for (yylen = 0; yystr[yylen]; yylen++)
    continue;
  return yylen;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined __GLIBC__ && defined _STRING_H && defined _GNU_SOURCE
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static char *
yystpcpy (char *yydest, const char *yysrc)
#else
static char *
yystpcpy (yydest, yysrc)
    char *yydest;
    const char *yysrc;
#endif
{
  char *yyd = yydest;
  const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif

# ifndef yytnamerr
/* Copy to YYRES the contents of YYSTR after stripping away unnecessary
   quotes and backslashes, so that it's suitable for yyerror.  The
   heuristic is that double-quoting is unnecessary unless the string
   contains an apostrophe, a comma, or backslash (other than
   backslash-backslash).  YYSTR is taken from yytname.  If YYRES is
   null, do not copy; instead, return the length of what the result
   would have been.  */
static YYSIZE_T
yytnamerr (char *yyres, const char *yystr)
{
  if (*yystr == '"')
    {
      YYSIZE_T yyn = 0;
      char const *yyp = yystr;

      for (;;)
	switch (*++yyp)
	  {
	  case '\'':
	  case ',':
	    goto do_not_strip_quotes;

	  case '\\':
	    if (*++yyp != '\\')
	      goto do_not_strip_quotes;
	    /* Fall through.  */
	  default:
	    if (yyres)
	      yyres[yyn] = *yyp;
	    yyn++;
	    break;

	  case '"':
	    if (yyres)
	      yyres[yyn] = '\0';
	    return yyn;
	  }
    do_not_strip_quotes: ;
    }

  if (! yyres)
    return yystrlen (yystr);

  return yystpcpy (yyres, yystr) - yyres;
}
# endif

/* Copy into YYRESULT an error message about the unexpected token
   YYCHAR while in state YYSTATE.  Return the number of bytes copied,
   including the terminating null byte.  If YYRESULT is null, do not
   copy anything; just return the number of bytes that would be
   copied.  As a special case, return 0 if an ordinary "syntax error"
   message will do.  Return YYSIZE_MAXIMUM if overflow occurs during
   size calculation.  */
static YYSIZE_T
yysyntax_error (char *yyresult, int yystate, int yychar)
{
  int yyn = yypact[yystate];

  if (! (YYPACT_NINF < yyn && yyn <= YYLAST))
    return 0;
  else
    {
      int yytype = YYTRANSLATE (yychar);
      YYSIZE_T yysize0 = yytnamerr (0, yytname[yytype]);
      YYSIZE_T yysize = yysize0;
      YYSIZE_T yysize1;
      int yysize_overflow = 0;
      enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
      char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
      int yyx;

# if 0
      /* This is so xgettext sees the translatable formats that are
	 constructed on the fly.  */
      YY_("syntax error, unexpected %s");
      YY_("syntax error, unexpected %s, expecting %s");
      YY_("syntax error, unexpected %s, expecting %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s");
      YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s");
# endif
      char *yyfmt;
      char const *yyf;
      static char const yyunexpected[] = "syntax error, unexpected %s";
      static char const yyexpecting[] = ", expecting %s";
      static char const yyor[] = " or %s";
      char yyformat[sizeof yyunexpected
		    + sizeof yyexpecting - 1
		    + ((YYERROR_VERBOSE_ARGS_MAXIMUM - 2)
		       * (sizeof yyor - 1))];
      char const *yyprefix = yyexpecting;

      /* Start YYX at -YYN if negative to avoid negative indexes in
	 YYCHECK.  */
      int yyxbegin = yyn < 0 ? -yyn : 0;

      /* Stay within bounds of both yycheck and yytname.  */
      int yychecklim = YYLAST - yyn + 1;
      int yyxend = yychecklim < YYNTOKENS ? yychecklim : YYNTOKENS;
      int yycount = 1;

      yyarg[0] = yytname[yytype];
      yyfmt = yystpcpy (yyformat, yyunexpected);

      for (yyx = yyxbegin; yyx < yyxend; ++yyx)
	if (yycheck[yyx + yyn] == yyx && yyx != YYTERROR)
	  {
	    if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
	      {
		yycount = 1;
		yysize = yysize0;
		yyformat[sizeof yyunexpected - 1] = '\0';
		break;
	      }
	    yyarg[yycount++] = yytname[yyx];
	    yysize1 = yysize + yytnamerr (0, yytname[yyx]);
	    yysize_overflow |= (yysize1 < yysize);
	    yysize = yysize1;
	    yyfmt = yystpcpy (yyfmt, yyprefix);
	    yyprefix = yyor;
	  }

      yyf = YY_(yyformat);
      yysize1 = yysize + yystrlen (yyf);
      yysize_overflow |= (yysize1 < yysize);
      yysize = yysize1;

      if (yysize_overflow)
	return YYSIZE_MAXIMUM;

      if (yyresult)
	{
	  /* Avoid sprintf, as that infringes on the user's name space.
	     Don't have undefined behavior even if the translation
	     produced a string with the wrong number of "%s"s.  */
	  char *yyp = yyresult;
	  int yyi = 0;
	  while ((*yyp = *yyf) != '\0')
	    {
	      if (*yyp == '%' && yyf[1] == 's' && yyi < yycount)
		{
		  yyp += yytnamerr (yyp, yyarg[yyi++]);
		  yyf += 2;
		}
	      else
		{
		  yyp++;
		  yyf++;
		}
	    }
	}
      return yysize;
    }
}
#endif /* YYERROR_VERBOSE */


/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

/*ARGSUSED*/
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
static void
yydestruct (const char *yymsg, int yytype, YYSTYPE *yyvaluep)
#else
static void
yydestruct (yymsg, yytype, yyvaluep)
    const char *yymsg;
    int yytype;
    YYSTYPE *yyvaluep;
#endif
{
  YYUSE (yyvaluep);

  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yytype, yyvaluep, yylocationp);

  switch (yytype)
    {

      default:
	break;
    }
}


/* Prevent warnings from -Wmissing-prototypes.  */

#ifdef YYPARSE_PARAM
#if defined __STDC__ || defined __cplusplus
int yyparse (void *YYPARSE_PARAM);
#else
int yyparse ();
#endif
#else /* ! YYPARSE_PARAM */
#if defined __STDC__ || defined __cplusplus
int yyparse (void);
#else
int yyparse ();
#endif
#endif /* ! YYPARSE_PARAM */



/* The look-ahead symbol.  */
int yychar;

/* The semantic value of the look-ahead symbol.  */
YYSTYPE yylval;

/* Number of syntax errors so far.  */
int yynerrs;



/*----------.
| yyparse.  |
`----------*/

#ifdef YYPARSE_PARAM
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void *YYPARSE_PARAM)
#else
int
yyparse (YYPARSE_PARAM)
    void *YYPARSE_PARAM;
#endif
#else /* ! YYPARSE_PARAM */
#if (defined __STDC__ || defined __C99__FUNC__ \
     || defined __cplusplus || defined _MSC_VER)
int
yyparse (void)
#else
int
yyparse ()

#endif
#endif
{
  
  int yystate;
  int yyn;
  int yyresult;
  /* Number of tokens to shift before error messages enabled.  */
  int yyerrstatus;
  /* Look-ahead token as an internal (translated) token number.  */
  int yytoken = 0;
#if YYERROR_VERBOSE
  /* Buffer for error messages, and its allocated size.  */
  char yymsgbuf[128];
  char *yymsg = yymsgbuf;
  YYSIZE_T yymsg_alloc = sizeof yymsgbuf;
#endif

  /* Three stacks and their tools:
     `yyss': related to states,
     `yyvs': related to semantic values,
     `yyls': related to locations.

     Refer to the stacks thru separate pointers, to allow yyoverflow
     to reallocate them elsewhere.  */

  /* The state stack.  */
  yytype_int16 yyssa[YYINITDEPTH];
  yytype_int16 *yyss = yyssa;
  yytype_int16 *yyssp;

  /* The semantic value stack.  */
  YYSTYPE yyvsa[YYINITDEPTH];
  YYSTYPE *yyvs = yyvsa;
  YYSTYPE *yyvsp;



#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  YYSIZE_T yystacksize = YYINITDEPTH;

  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;


  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yystate = 0;
  yyerrstatus = 0;
  yynerrs = 0;
  yychar = YYEMPTY;		/* Cause a token to be read.  */

  /* Initialize stack pointers.
     Waste one element of value and location stack
     so that they stay on the same level as the state stack.
     The wasted elements are never initialized.  */

  yyssp = yyss;
  yyvsp = yyvs;

  goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
 yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
	/* Give user a chance to reallocate the stack.  Use copies of
	   these so that the &'s don't force the real ones into
	   memory.  */
	YYSTYPE *yyvs1 = yyvs;
	yytype_int16 *yyss1 = yyss;


	/* Each stack pointer address is followed by the size of the
	   data in use in that stack, in bytes.  This used to be a
	   conditional around just the two extra args, but that might
	   be undefined if yyoverflow is a macro.  */
	yyoverflow (YY_("memory exhausted"),
		    &yyss1, yysize * sizeof (*yyssp),
		    &yyvs1, yysize * sizeof (*yyvsp),

		    &yystacksize);

	yyss = yyss1;
	yyvs = yyvs1;
      }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
      goto yyexhaustedlab;
# else
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
	goto yyexhaustedlab;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
	yystacksize = YYMAXDEPTH;

      {
	yytype_int16 *yyss1 = yyss;
	union yyalloc *yyptr =
	  (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
	if (! yyptr)
	  goto yyexhaustedlab;
	YYSTACK_RELOCATE (yyss);
	YYSTACK_RELOCATE (yyvs);

#  undef YYSTACK_RELOCATE
	if (yyss1 != yyssa)
	  YYSTACK_FREE (yyss1);
      }
# endif
#endif /* no yyoverflow */

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;


      YYDPRINTF ((stderr, "Stack size increased to %lu\n",
		  (unsigned long int) yystacksize));

      if (yyss + yystacksize - 1 <= yyssp)
	YYABORT;
    }

  YYDPRINTF ((stderr, "Entering state %d\n", yystate));

  goto yybackup;

/*-----------.
| yybackup.  |
`-----------*/
yybackup:

  /* Do appropriate processing given the current state.  Read a
     look-ahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to look-ahead token.  */
  yyn = yypact[yystate];
  if (yyn == YYPACT_NINF)
    goto yydefault;

  /* Not known => get a look-ahead token if don't already have one.  */

  /* YYCHAR is either YYEMPTY or YYEOF or a valid look-ahead symbol.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token: "));
      yychar = YYLEX;
    }

  if (yychar <= YYEOF)
    {
      yychar = yytoken = YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yyn == 0 || yyn == YYTABLE_NINF)
	goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  if (yyn == YYFINAL)
    YYACCEPT;

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the look-ahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the shifted token unless it is eof.  */
  if (yychar != YYEOF)
    yychar = YYEMPTY;

  yystate = yyn;
  *++yyvsp = yylval;

  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     `$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
        case 6:
#line 167 "Gram.y"
    {
                  if( (yyvsp[(3) - (3)])->to_data()->data_type() != MAC_Data::String )
                       MACerror( "Include must refer to string expression" ) ;
                  if( YY_top_module!=dummy_module )
                  {
                    if( ! (yyvsp[(3) - (3)])->to_data()->value_can_be_evaluated(
                                               YY_top_module->context() ) )
                    {
                       std::string msg = "Include path can not be evaluated\n"  ;
                       msg += "Undefined variable(s):\n" ;
                       stringVector const& undef =
                          (yyvsp[(3) - (3)])->to_data()->undefined_variables(
                                                  YY_top_module->context() ) ;
                       for( size_t i=0 ; i<undef.size() ; ++i )
                       {
                          msg += "   - \""+undef(i)+"\"\n" ;
                       }
                       MACerror( msg.c_str() ) ;
                    }
                    std::string file_data =
                       (yyvsp[(3) - (3)])->to_data()->to_string( YY_top_module->context() ) ;
                    readFile( YY_top_module, 0, file_data, false ) ;
                  }
               ;}
    break;

  case 7:
#line 191 "Gram.y"
    { if( endFile() ) YYACCEPT ; ;}
    break;

  case 8:
#line 193 "Gram.y"
    { (yyval) = (yyvsp[(2) - (3)]) ; ;}
    break;

  case 16:
#line 212 "Gram.y"
    {
                if( YY_top_module!=dummy_module )
                {
                   if( (yyvsp[(1) - (3)])->to_data()->data_type()!=MAC_Data::String )
                      MACerror( "Key must be a string" ) ;
                   std::string const& str = (yyvsp[(1) - (3)])->to_data()->to_string() ;
                   if( YY_top_module->has_entry( str ) )
                   {
                     std::string mess = "Module " ;
                     mess+=YY_top_module->name() + " already contains "+
                        str + " assignment (use ==)" ;
                     MACerror( mess.c_str() ) ;
                   }
                   MAC_Data* data = (yyvsp[(3) - (3)])->to_data() ;
                   (yyvsp[(3) - (3)])->change_owner(YY_top_module, data) ;
                   if( YY_top_module->has_module( str ) )
                   {
                      std::string mess = "Can't create entry and module with the same name "+str+" in module " + YY_top_module->name() ;
                      MACerror( mess.c_str() ) ;
                   }
                   YY_top_module->add_entry( str, data ) ;
                }
                (yyval) = 0 ;
             ;}
    break;

  case 17:
#line 238 "Gram.y"
    {
                if( YY_top_module!=dummy_module )
                {
                   MAC_Variable* var = dynamic_cast<MAC_Variable*>(
                      (yyvsp[(1) - (3)])->to_data() ) ;
                   (yyvsp[(1) - (3)])->change_owner(YY_top_module,var) ;
                   MAC_ASSERT( var!=0 ) ;
                   if( var->data_type()!=(yyvsp[(3) - (3)])->to_data()->data_type() )
                   {
                      std::string mess = "Bad type for expression assigned to " ;
                      mess+=var->name() ;
                      MACerror( mess.c_str() ) ;
                   }
                   if( !YY_top_module->context()->has_variable(var) )
                   {
                      MAC_Data* data = (yyvsp[(3) - (3)])->to_data() ;
                      (yyvsp[(3) - (3)])->change_owner(0,data) ;
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
                (yyval) = 0 ;
             ;}
    break;

  case 18:
#line 275 "Gram.y"
    {
                if( YY_top_module!=dummy_module )
                {
                   MAC_Variable* var = dynamic_cast<MAC_Variable*>(
                      (yyvsp[(1) - (4)])->to_data() ) ;
                   MAC_ASSERT( var!=0 ) ;
                   if( var->data_type()!=(yyvsp[(4) - (4)])->to_data()->data_type() )
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
                   MAC_Data* data = (yyvsp[(4) - (4)])->to_data() ;
                   (yyvsp[(4) - (4)])->change_owner(0,data) ;
                   YY_top_module->modify_variable( var, data ) ;
                }
                (yyval) = 0 ;
             ;}
    break;

  case 19:
#line 300 "Gram.y"
    {
                if( YY_top_module!=dummy_module )
                {
                   if( (yyvsp[(1) - (4)])->to_data()->data_type()!=MAC_Data::String )
                      MACerror( "Key must be a string" ) ;
                   MAC_String * str = static_cast<MAC_String*>( (yyvsp[(1) - (4)])->to_data() ) ;
                   if( !YY_top_module->has_entry( str->to_string() ) )
                   {
                      std::string mess = str->to_string() + " doesn't already exist " ;
                      MACerror( mess.c_str() ) ;
                   }
                   substitute_for_assignment( str->to_string(),
                                              (yyvsp[(4) - (4)])->to_data() ) ;
                }
                (yyval) = 0 ;
             ;}
    break;

  case 28:
#line 326 "Gram.y"
    {
             (yyval) = (yyvsp[(2) - (3)]) ;
             if( (yyval)->is_data() )
             {
                MAC_Expression* op = dynamic_cast<MAC_Expression*>( (yyval)->to_data() ) ;
                if( op != 0 ) op->set_external_brackets() ;
             }
          ;}
    break;

  case 29:
#line 337 "Gram.y"
    {
             std::string const& op_name = (yyvsp[(1) - (4)])->to_data()->to_string() ;
             if( op_name=="this_file_dir" )
             {
                (yyval)=MAC_Lexical::create(
                   MAC_String::create( 0, relativePath ) ) ;
             }
             else
             {
                CREATE_OP( op, op_name, (yyvsp[(3) - (4)])->to_list() ) ; 
                (yyval)=MAC_Lexical::create( op ) ;
                (yyvsp[(3) - (4)])->change_owner(op,(yyvsp[(3) - (4)])->to_list()) ;
             }
          ;}
    break;

  case 30:
#line 353 "Gram.y"
    {
             std::string const& str = (yyvsp[(2) - (2)])->to_data()->to_string() ;
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
                (yyvsp[(2) - (2)])->to_data()->to_string() ) ;
             (yyval)=MAC_Lexical::create( var->create_clone(0) ) ;
          ;}
    break;

  case 31:
#line 379 "Gram.y"
    {
             MAC_List * lst = (yyvsp[(2) - (4)])->to_list() ;
             (yyvsp[(3) - (4)])->change_owner(lst,(yyvsp[(3) - (4)])->to_data() ) ;
             lst->append( (yyvsp[(3) - (4)])->to_data() ) ;
             CREATE_OP( op, "(?:)",lst ) ;
             op->set_external_brackets() ;
             (yyvsp[(2) - (4)])->change_owner(op,lst) ;
             (yyval)=MAC_Lexical::create( op ) ;
          ;}
    break;

  case 32:
#line 390 "Gram.y"
    {
                 MAC_LABEL( "Gram.y :: switches_list1" ) ;
                 MAC_List * lst = MAC_List::create( 0 ) ;
                 (yyvsp[(1) - (4)])->change_owner(lst,(yyvsp[(1) - (4)])->to_data() ) ;
                 lst->append( (yyvsp[(1) - (4)])->to_data() ) ;
                 (yyvsp[(3) - (4)])->change_owner(lst,(yyvsp[(3) - (4)])->to_data() ) ;
                 lst->append( (yyvsp[(3) - (4)])->to_data() ) ;
                 (yyval) = MAC_Lexical::create( lst ) ;
              ;}
    break;

  case 33:
#line 400 "Gram.y"
    {
                 MAC_LABEL( "Gram.y :: switches_list2" ) ;
                 MAC_List * lst = (yyvsp[(1) - (5)])->to_list() ;
                 (yyvsp[(2) - (5)])->change_owner(lst,(yyvsp[(2) - (5)])->to_data() ) ;
                 lst->append( (yyvsp[(2) - (5)])->to_data() ) ;
                 (yyvsp[(4) - (5)])->change_owner(lst,(yyvsp[(4) - (5)])->to_data() ) ;
                 lst->append( (yyvsp[(4) - (5)])->to_data() ) ;
                 (yyval) = (yyvsp[(1) - (5)]) ;
              ;}
    break;

  case 34:
#line 410 "Gram.y"
    {  (yyval)=MAC_Lexical::create( MAC_List::create( 0 ) ) ; ;}
    break;

  case 35:
#line 412 "Gram.y"
    {
                 MAC_List * lst = MAC_List::create( 0 ) ;
                 (yyvsp[(1) - (1)])->change_owner(lst,(yyvsp[(1) - (1)])->to_data() ) ;
                 lst->append( (yyvsp[(1) - (1)])->to_data() ) ;
                 (yyval) = MAC_Lexical::create( lst ) ;
              ;}
    break;

  case 36:
#line 419 "Gram.y"
    {
                 MAC_List * lst = (yyvsp[(1) - (3)])->to_list() ;
                 (yyvsp[(3) - (3)])->change_owner(lst,(yyvsp[(3) - (3)])->to_data() ) ;
                 lst->append( (yyvsp[(3) - (3)])->to_data() ) ;
                 (yyval) = (yyvsp[(1) - (3)]) ;
              ;}
    break;

  case 37:
#line 426 "Gram.y"
    {
                    if( (yyvsp[(2) - (2)])->to_data()->is_raw_data() )
                    {
                       MAC_Data::Type dt = (yyvsp[(2) - (2)])->to_data()->data_type() ;
                       
                       if( dt == MAC_Data::Double )
                       {
                          (yyval)=MAC_Lexical::create(
                             MAC_Double::create( 0, - (yyvsp[(2) - (2)])->to_data()->to_double() ) ) ;
                       }
                       else if( dt == MAC_Data::Int )
                       {
                          (yyval)=MAC_Lexical::create(
                             MAC_Int::create( 0, - (yyvsp[(2) - (2)])->to_data()->to_int() ) ) ;
                       }
                       else
                       {
                          MACerror( "Undefined unary minus operator" ) ;
                       }
                    }
                    else
                    {
                       CREATE_SINGLE_OP("unary_minus",(yyvsp[(2) - (2)]),(yyval)) ; 
                    } ;}
    break;

  case 38:
#line 450 "Gram.y"
    { CREATE_SINGLE_OP("!",(yyvsp[(2) - (2)]),(yyval)) ; ;}
    break;

  case 39:
#line 452 "Gram.y"
    { CREATE_BIN_OP("<",(yyvsp[(1) - (3)]),(yyvsp[(3) - (3)]),(yyval)) ; ;}
    break;

  case 40:
#line 453 "Gram.y"
    { CREATE_BIN_OP(">",(yyvsp[(1) - (3)]),(yyvsp[(3) - (3)]),(yyval)) ; ;}
    break;

  case 41:
#line 454 "Gram.y"
    { CREATE_BIN_OP("+",(yyvsp[(1) - (3)]),(yyvsp[(3) - (3)]),(yyval)) ; ;}
    break;

  case 42:
#line 455 "Gram.y"
    { CREATE_BIN_OP("-",(yyvsp[(1) - (3)]),(yyvsp[(3) - (3)]),(yyval)) ; ;}
    break;

  case 43:
#line 456 "Gram.y"
    { CREATE_BIN_OP("*",(yyvsp[(1) - (3)]),(yyvsp[(3) - (3)]),(yyval)) ; ;}
    break;

  case 44:
#line 457 "Gram.y"
    { CREATE_BIN_OP("=",(yyvsp[(1) - (3)]),(yyvsp[(3) - (3)]),(yyval)) ; ;}
    break;

  case 45:
#line 458 "Gram.y"
    { CREATE_BIN_OP("/",(yyvsp[(1) - (3)]),(yyvsp[(3) - (3)]),(yyval)) ; ;}
    break;

  case 46:
#line 459 "Gram.y"
    { CREATE_BIN_OP("!=",(yyvsp[(1) - (3)]),(yyvsp[(3) - (3)]),(yyval)) ; ;}
    break;

  case 47:
#line 460 "Gram.y"
    { CREATE_BIN_OP("<<",(yyvsp[(1) - (3)]),(yyvsp[(3) - (3)]),(yyval)) ; ;}
    break;

  case 48:
#line 461 "Gram.y"
    { CREATE_BIN_OP("||",(yyvsp[(1) - (3)]),(yyvsp[(3) - (3)]),(yyval)) ; ;}
    break;

  case 49:
#line 462 "Gram.y"
    { CREATE_BIN_OP("&&",(yyvsp[(1) - (3)]),(yyvsp[(3) - (3)]),(yyval)) ; ;}
    break;

  case 50:
#line 463 "Gram.y"
    { CREATE_BIN_OP("<=",(yyvsp[(1) - (3)]),(yyvsp[(3) - (3)]),(yyval)) ; ;}
    break;

  case 51:
#line 464 "Gram.y"
    { CREATE_BIN_OP(">=",(yyvsp[(1) - (3)]),(yyvsp[(3) - (3)]),(yyval)) ; ;}
    break;

  case 54:
#line 470 "Gram.y"
    {
            MAC_ASSERT( (yyvsp[(3) - (3)])==0 || YY_top_module!=(yyvsp[(3) - (3)])->to_module() ) ;
            (yyval) = (yyvsp[(3) - (3)]) ;
         ;}
    break;

  case 55:
#line 475 "Gram.y"
    { (yyval) = MAC_Lexical::create( MAC_Bool::create( 0, true ) ) ; ;}
    break;

  case 56:
#line 477 "Gram.y"
    {
            if( (yyvsp[(3) - (4)])->to_data()->data_type() != MAC_Data::Bool )
               MACerror( "Conditional module if must refer to boolean expression" ) ;
            bool cond = false ;
            if( YY_top_module!=dummy_module )
            {
               if( ! (yyvsp[(3) - (4)])->to_data()->value_can_be_evaluated(
                                             YY_top_module->context() ) )
               {
                  std::string msg = "Conditional module if can not be evaluated\n" ;
                  msg += "Undefined variable(s):\n" ;
                  stringVector const& undef =
                     (yyvsp[(3) - (4)])->to_data()->undefined_variables(
                                                  YY_top_module->context() ) ;
                  for( size_t i=0 ; i<undef.size() ; ++i )
                  {
                     msg += "   - \""+undef(i)+"\"\n" ;
                  }
                  MACerror( msg.c_str() ) ;
               }
               cond = (yyvsp[(3) - (4)])->to_data()->to_bool( YY_top_module->context() ) ;
            }
            (yyval) = MAC_Lexical::create( MAC_Bool::create( 0, cond ) ) ;
         ;}
    break;

  case 57:
#line 503 "Gram.y"
    {
          modules_LILO->append( YY_top_module ) ;
          if( (yyvsp[(1) - (3)])->to_data()->to_bool() && YY_top_module!=dummy_module )
          {
             MAC_Module* old = YY_top_module ;
             std::string a_name = (yyvsp[(3) - (3)])->to_data()->to_string() ;
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
       ;}
    break;

  case 58:
#line 533 "Gram.y"
    {
          if( YY_top_module!=dummy_module )
          {
             if( (yyvsp[(3) - (3)])->to_data()->to_string()!=YY_top_module->name() )
             {
                std::string mess = "Module " + YY_top_module->name() +
                   " can't be closed with " + (yyvsp[(3) - (3)])->to_data()->to_string() ;
                MACerror( mess.c_str() ) ;
             }
             (yyval)=MAC_Lexical::create( YY_top_module ) ;
          }
          else
          {
             (yyval)=0 ;
          }
          size_t last = modules_LILO->count()-1 ;
          YY_top_module = static_cast<MAC_Module*>( modules_LILO->at( last ) ) ;
          modules_LILO->remove_at( last ) ;
          if( modules_LILO->has( dummy_module ) )
          {
             YY_top_module = dummy_module ;
          }
       ;}
    break;

  case 59:
#line 561 "Gram.y"
    {
            MAC_List * lst = (yyvsp[(2) - (3)])->to_list() ;
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
            (yyval) = MAC_Lexical::create( res ) ;
         ;}
    break;

  case 60:
#line 596 "Gram.y"
    { (yyval) = MAC_Lexical::create( MAC_List::create( 0 ) ) ; ;}
    break;

  case 61:
#line 598 "Gram.y"
    {
             (yyvsp[(1) - (2)])->to_list()->append( (yyvsp[(2) - (2)])->to_data() ) ;
             (yyval)=(yyvsp[(1) - (2)]);
          ;}
    break;

  case 62:
#line 604 "Gram.y"
    {
            MAC_List * lst = (yyvsp[(2) - (3)])->to_list() ;
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
            (yyval) = MAC_Lexical::create( res ) ;
;}
    break;

  case 63:
#line 640 "Gram.y"
    { (yyval) = MAC_Lexical::create( MAC_List::create( 0 ) ) ;
            (yyval)->to_list()->append((yyvsp[(1) - (1)])->to_data() ) ;
          ;}
    break;

  case 64:
#line 644 "Gram.y"
    {
             (yyvsp[(1) - (3)])->to_list()->append( (yyvsp[(3) - (3)])->to_data() ) ;
             (yyval)=(yyvsp[(1) - (3)]);
          ;}
    break;

  case 70:
#line 657 "Gram.y"
    { (yyval) = MAC_Lexical::create( MAC_Bool::create( 0, true ) ) ; ;}
    break;

  case 71:
#line 658 "Gram.y"
    { (yyval) = MAC_Lexical::create( MAC_Bool::create( 0, false ) ) ; ;}
    break;


/* Line 1267 of yacc.c.  */
#line 2216 "y.tab.c"
      default: break;
    }
  YY_SYMBOL_PRINT ("-> $$ =", yyr1[yyn], &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);

  *++yyvsp = yyval;


  /* Now `shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */

  yyn = yyr1[yyn];

  yystate = yypgoto[yyn - YYNTOKENS] + *yyssp;
  if (0 <= yystate && yystate <= YYLAST && yycheck[yystate] == *yyssp)
    yystate = yytable[yystate];
  else
    yystate = yydefgoto[yyn - YYNTOKENS];

  goto yynewstate;


/*------------------------------------.
| yyerrlab -- here on detecting error |
`------------------------------------*/
yyerrlab:
  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
#if ! YYERROR_VERBOSE
      yyerror (YY_("syntax error"));
#else
      {
	YYSIZE_T yysize = yysyntax_error (0, yystate, yychar);
	if (yymsg_alloc < yysize && yymsg_alloc < YYSTACK_ALLOC_MAXIMUM)
	  {
	    YYSIZE_T yyalloc = 2 * yysize;
	    if (! (yysize <= yyalloc && yyalloc <= YYSTACK_ALLOC_MAXIMUM))
	      yyalloc = YYSTACK_ALLOC_MAXIMUM;
	    if (yymsg != yymsgbuf)
	      YYSTACK_FREE (yymsg);
	    yymsg = (char *) YYSTACK_ALLOC (yyalloc);
	    if (yymsg)
	      yymsg_alloc = yyalloc;
	    else
	      {
		yymsg = yymsgbuf;
		yymsg_alloc = sizeof yymsgbuf;
	      }
	  }

	if (0 < yysize && yysize <= yymsg_alloc)
	  {
	    (void) yysyntax_error (yymsg, yystate, yychar);
	    yyerror (yymsg);
	  }
	else
	  {
	    yyerror (YY_("syntax error"));
	    if (yysize != 0)
	      goto yyexhaustedlab;
	  }
      }
#endif
    }



  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse look-ahead token after an
	 error, discard it.  */

      if (yychar <= YYEOF)
	{
	  /* Return failure if at end of input.  */
	  if (yychar == YYEOF)
	    YYABORT;
	}
      else
	{
	  yydestruct ("Error: discarding",
		      yytoken, &yylval);
	  yychar = YYEMPTY;
	}
    }

  /* Else will try to reuse look-ahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:

  /* Pacify compilers like GCC when the user code never invokes
     YYERROR and the label yyerrorlab therefore never appears in user
     code.  */
  if (/*CONSTCOND*/ 0)
     goto yyerrorlab;

  /* Do not reclaim the symbols of the rule which action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;	/* Each real token shifted decrements this.  */

  for (;;)
    {
      yyn = yypact[yystate];
      if (yyn != YYPACT_NINF)
	{
	  yyn += YYTERROR;
	  if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYTERROR)
	    {
	      yyn = yytable[yyn];
	      if (0 < yyn)
		break;
	    }
	}

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
	YYABORT;


      yydestruct ("Error: popping",
		  yystos[yystate], yyvsp);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  if (yyn == YYFINAL)
    YYACCEPT;

  *++yyvsp = yylval;


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", yystos[yyn], yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturn;

#ifndef yyoverflow
/*-------------------------------------------------.
| yyexhaustedlab -- memory exhaustion comes here.  |
`-------------------------------------------------*/
yyexhaustedlab:
  yyerror (YY_("memory exhausted"));
  yyresult = 2;
  /* Fall through.  */
#endif

yyreturn:
  if (yychar != YYEOF && yychar != YYEMPTY)
     yydestruct ("Cleanup: discarding lookahead",
		 yytoken, &yylval);
  /* Do not reclaim the symbols of the rule which action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
		  yystos[*yyssp], yyvsp);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
#if YYERROR_VERBOSE
  if (yymsg != yymsgbuf)
    YYSTACK_FREE (yymsg);
#endif
  /* Make sure YYID is used.  */
  return YYID (yyresult);
}


#line 660 "Gram.y"

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

