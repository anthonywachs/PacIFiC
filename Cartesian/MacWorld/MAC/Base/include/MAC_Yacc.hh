/* A Bison parser, made by GNU Bison 2.3.  */

/* Skeleton interface for Bison's Yacc-like parsers in C

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




#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef int YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif

extern YYSTYPE MAClval;

