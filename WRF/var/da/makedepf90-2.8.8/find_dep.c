/* A Bison parser, made by GNU Bison 2.1.  */

/* Skeleton parser for Yacc-like parsing with Bison,
   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004, 2005 Free Software Foundation, Inc.

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

/* As a special exception, when this file is copied by Bison into a
   Bison output file, you may use that output file without restriction.
   This special exception was added by the Free Software Foundation
   in version 1.24 of Bison.  */

/* Written by Richard Stallman by simplifying the original so called
   ``semantic'' parser.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output.  */
#define YYBISON 1

/* Bison version.  */
#define YYBISON_VERSION "2.1"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Using locations.  */
#define YYLSP_NEEDED 0



/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     EOSTMT = 258,
     ASSIGNMENT_OP = 259,
     GARBAGE = 260,
     CPP_INCLUDE = 261,
     F90PPR_INCLUDE = 262,
     COCO_INCLUDE = 263,
     F90PPR_DEFINE = 264,
     CPP_DEFINE = 265,
     F90PPR_UNDEF = 266,
     CPP_UNDEF = 267,
     CPP_IFDEF = 268,
     CPP_IFNDEF = 269,
     CPP_IF = 270,
     CPP_ELSE = 271,
     CPP_ELIF = 272,
     CPP_ENDIF = 273,
     F90PPR_IFDEF = 274,
     F90PPR_IFNDEF = 275,
     F90PPR_IF = 276,
     F90PPR_ELSE = 277,
     F90PPR_ELIF = 278,
     F90PPR_ENDIF = 279,
     CPP_TOENDL = 280,
     UNTERMINATED_STRING = 281,
     STRING = 282,
     WORD = 283
   };
#endif
/* Tokens.  */
#define EOSTMT 258
#define ASSIGNMENT_OP 259
#define GARBAGE 260
#define CPP_INCLUDE 261
#define F90PPR_INCLUDE 262
#define COCO_INCLUDE 263
#define F90PPR_DEFINE 264
#define CPP_DEFINE 265
#define F90PPR_UNDEF 266
#define CPP_UNDEF 267
#define CPP_IFDEF 268
#define CPP_IFNDEF 269
#define CPP_IF 270
#define CPP_ELSE 271
#define CPP_ELIF 272
#define CPP_ENDIF 273
#define F90PPR_IFDEF 274
#define F90PPR_IFNDEF 275
#define F90PPR_IF 276
#define F90PPR_ELSE 277
#define F90PPR_ELIF 278
#define F90PPR_ENDIF 279
#define CPP_TOENDL 280
#define UNTERMINATED_STRING 281
#define STRING 282
#define WORD 283




/* Copy the first part of user declarations.  */
#line 1 "find_dep.y"

/* 
 * Copyright (C) 2000-2006 Erik Edelmann <Erik.Edelmann@iki.fi>
 *
 *     This program is free software;  you  can  redistribute  it
 *     and/or modify it under the terms of the GNU General Public
 *     License version 2 as published  by  the  Free  Software  
 *     Foundation.
 *
 *     This program is distributed in the hope that  it  will  be
 *     useful, but WITHOUT ANY WARRANTY; without even the implied
 *     warranty of MERCHANTABILITY or FITNESS  FOR  A  PARTICULAR
 *     PURPOSE.   See  the  GNU  General  Public License for more
 *     details.
 *
 *     You should have received a copy of the GNU General  Public
 *     License along with this program; if not, write to the Free
 *     Software Foundation, Inc., 59  Temple  Place,  Suite  330,
 *     Boston, MA  02111-1307  USA
 */

#include <errno.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef HAVE_ALLOCA_H
#  include <alloca.h>
#endif

#include "finddep.h"
#include "utils.h"
#include "errormesg.h"
#include "global.h"
#include "macro.h"
#include "modfile_name.h"

static char *sourcefile;
static char *curr_file;
static Dependency *dep;  /* Dependencies of the file */
static List *modules;    /* Modules defined in the file */
static List *macrolist;
static char *filestack[INCLUDE_RECURSION_LIMIT+1];
static int filestack_i;
static int pp_ignore;    /* In 'false'-branch of a pre-processor 'if' */
static bool skip_to_end[20];
static int skip_i;
static bool in_interface = false;

int yyerror(const char *s);

static int modcmp(const void *m1, const void *m2);
SourceFmt get_format(const char *filename);
void handle_include(const char *incfile);

Macro *defmac;

/* Defined in lexer.l */
int yylex();
bool lex_include_file(const char *incfile);
void lex_set_format(SourceFmt fmt);



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

#if ! defined (YYSTYPE) && ! defined (YYSTYPE_IS_DECLARED)
#line 65 "find_dep.y"
typedef union YYSTYPE {
    char *string;
    int number;
} YYSTYPE;
/* Line 196 of yacc.c.  */
#line 210 "y.tab.c"
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif



/* Copy the second part of user declarations.  */


/* Line 219 of yacc.c.  */
#line 222 "y.tab.c"

#if ! defined (YYSIZE_T) && defined (__SIZE_TYPE__)
# define YYSIZE_T __SIZE_TYPE__
#endif
#if ! defined (YYSIZE_T) && defined (size_t)
# define YYSIZE_T size_t
#endif
#if ! defined (YYSIZE_T) && (defined (__STDC__) || defined (__cplusplus))
# include <stddef.h> /* INFRINGES ON USER NAME SPACE */
# define YYSIZE_T size_t
#endif
#if ! defined (YYSIZE_T)
# define YYSIZE_T unsigned int
#endif

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

#if ! defined (yyoverflow) || YYERROR_VERBOSE

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if defined (__STDC__) || defined (__cplusplus)
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#     define YYINCLUDED_STDLIB_H
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's `empty if-body' warning. */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (0)
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2005 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM ((YYSIZE_T) -1)
#  endif
#  ifdef __cplusplus
extern "C" {
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if (! defined (malloc) && ! defined (YYINCLUDED_STDLIB_H) \
	&& (defined (__STDC__) || defined (__cplusplus)))
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if (! defined (free) && ! defined (YYINCLUDED_STDLIB_H) \
	&& (defined (__STDC__) || defined (__cplusplus)))
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifdef __cplusplus
}
#  endif
# endif
#endif /* ! defined (yyoverflow) || YYERROR_VERBOSE */


#if (! defined (yyoverflow) \
     && (! defined (__cplusplus) \
	 || (defined (YYSTYPE_IS_TRIVIAL) && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  short int yyss;
  YYSTYPE yyvs;
  };

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (sizeof (short int) + sizeof (YYSTYPE))			\
      + YYSTACK_GAP_MAXIMUM)

/* Copy COUNT objects from FROM to TO.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined (__GNUC__) && 1 < __GNUC__
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
      while (0)
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
    while (0)

#endif

#if defined (__STDC__) || defined (__cplusplus)
   typedef signed char yysigned_char;
#else
   typedef short int yysigned_char;
#endif

/* YYFINAL -- State number of the termination state. */
#define YYFINAL  2
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   142

/* YYNTOKENS -- Number of terminals. */
#define YYNTOKENS  29
/* YYNNTS -- Number of nonterminals. */
#define YYNNTS  16
/* YYNRULES -- Number of rules. */
#define YYNRULES  49
/* YYNRULES -- Number of states. */
#define YYNSTATES  85

/* YYTRANSLATE(YYLEX) -- Bison symbol number corresponding to YYLEX.  */
#define YYUNDEFTOK  2
#define YYMAXUTOK   283

#define YYTRANSLATE(YYX)						\
  ((unsigned int) (YYX) <= YYMAXUTOK ? yytranslate[YYX] : YYUNDEFTOK)

/* YYTRANSLATE[YYLEX] -- Bison symbol number corresponding to YYLEX.  */
static const unsigned char yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
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
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28
};

#if YYDEBUG
/* YYPRHS[YYN] -- Index of the first RHS symbol of rule number YYN in
   YYRHS.  */
static const unsigned char yyprhs[] =
{
       0,     0,     3,     4,     7,     9,    11,    16,    19,    24,
      29,    34,    39,    44,    49,    54,    59,    63,    67,    71,
      75,    80,    84,    86,    88,    90,    92,    94,    96,    98,
     100,   102,   104,   106,   108,   110,   112,   114,   116,   118,
     120,   122,   124,   126,   127,   130,   132,   134,   136,   138
};

/* YYRHS -- A `-1'-separated list of the rules' RHS. */
static const yysigned_char yyrhs[] =
{
      30,     0,    -1,    -1,    30,    31,    -1,    33,    -1,    32,
      -1,    28,     4,    43,     3,    -1,    28,     3,    -1,    28,
      28,    43,     3,    -1,    28,    27,    43,     3,    -1,    34,
      27,    43,     3,    -1,     6,     5,    43,     3,    -1,    35,
      28,    43,     3,    -1,    36,    28,    43,     3,    -1,    37,
      28,    43,     3,    -1,    38,    28,    43,     3,    -1,    39,
      43,     3,    -1,    40,    43,     3,    -1,    41,    43,     3,
      -1,    42,    43,     3,    -1,    28,     5,    43,     3,    -1,
       5,    43,     3,    -1,     3,    -1,     1,    -1,     6,    -1,
       7,    -1,     8,    -1,    10,    -1,     9,    -1,    12,    -1,
      11,    -1,    13,    -1,    19,    -1,    14,    -1,    20,    -1,
      15,    -1,    21,    -1,    17,    -1,    23,    -1,    16,    -1,
      22,    -1,    18,    -1,    24,    -1,    -1,    43,    44,    -1,
      28,    -1,    27,    -1,     5,    -1,     4,    -1,    26,    -1
};

/* YYRLINE[YYN] -- source line where rule number YYN was defined.  */
static const unsigned short int yyrline[] =
{
       0,    81,    81,    82,    85,    86,    89,    91,    96,   136,
     143,   147,   148,   157,   174,   190,   206,   214,   219,   225,
     230,   231,   232,   233,   236,   237,   238,   241,   242,   245,
     246,   249,   250,   253,   254,   257,   258,   261,   262,   265,
     266,   269,   270,   273,   274,   277,   278,   279,   280,   281
};
#endif

#if YYDEBUG || YYERROR_VERBOSE || YYTOKEN_TABLE
/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals. */
static const char *const yytname[] =
{
  "$end", "error", "$undefined", "EOSTMT", "ASSIGNMENT_OP", "GARBAGE",
  "CPP_INCLUDE", "F90PPR_INCLUDE", "COCO_INCLUDE", "F90PPR_DEFINE",
  "CPP_DEFINE", "F90PPR_UNDEF", "CPP_UNDEF", "CPP_IFDEF", "CPP_IFNDEF",
  "CPP_IF", "CPP_ELSE", "CPP_ELIF", "CPP_ENDIF", "F90PPR_IFDEF",
  "F90PPR_IFNDEF", "F90PPR_IF", "F90PPR_ELSE", "F90PPR_ELIF",
  "F90PPR_ENDIF", "CPP_TOENDL", "UNTERMINATED_STRING", "STRING", "WORD",
  "$accept", "code", "stmt", "assignment_stmt", "keyword_stmt", "include",
  "define", "undef", "ifdef", "ifndef", "if", "elif", "else", "endif",
  "other", "misc_code", 0
};
#endif

# ifdef YYPRINT
/* YYTOKNUM[YYLEX-NUM] -- Internal token number corresponding to
   token YYLEX-NUM.  */
static const unsigned short int yytoknum[] =
{
       0,   256,   257,   258,   259,   260,   261,   262,   263,   264,
     265,   266,   267,   268,   269,   270,   271,   272,   273,   274,
     275,   276,   277,   278,   279,   280,   281,   282,   283
};
# endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives.  */
static const unsigned char yyr1[] =
{
       0,    29,    30,    30,    31,    31,    32,    33,    33,    33,
      33,    33,    33,    33,    33,    33,    33,    33,    33,    33,
      33,    33,    33,    33,    34,    34,    34,    35,    35,    36,
      36,    37,    37,    38,    38,    39,    39,    40,    40,    41,
      41,    42,    42,    43,    43,    44,    44,    44,    44,    44
};

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN.  */
static const unsigned char yyr2[] =
{
       0,     2,     0,     2,     1,     1,     4,     2,     4,     4,
       4,     4,     4,     4,     4,     4,     3,     3,     3,     3,
       4,     3,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     0,     2,     1,     1,     1,     1,     1
};

/* YYDEFACT[STATE-NAME] -- Default rule to reduce with in state
   STATE-NUM when YYTABLE doesn't specify something else to do.  Zero
   means the default is an error.  */
static const unsigned char yydefact[] =
{
       2,     0,     1,    23,    22,    43,    24,    25,    26,    28,
      27,    30,    29,    31,    33,    35,    39,    37,    41,    32,
      34,    36,    40,    38,    42,     0,     3,     5,     4,     0,
       0,     0,     0,     0,    43,    43,    43,    43,     0,    43,
       7,    43,    43,    43,    43,    43,    43,    43,    43,    43,
       0,     0,     0,     0,    21,    48,    47,    49,    46,    45,
      44,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,    16,    17,    18,    19,    11,     6,    20,     9,     8,
      10,    12,    13,    14,    15
};

/* YYDEFGOTO[NTERM-NUM]. */
static const yysigned_char yydefgoto[] =
{
      -1,     1,    26,    27,    28,    29,    30,    31,    32,    33,
      34,    35,    36,    37,    38,    60
};

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
#define YYPACT_NINF -24
static const yysigned_char yypact[] =
{
     -24,     0,   -24,   -24,   -24,   -24,    -3,   -24,   -24,   -24,
     -24,   -24,   -24,   -24,   -24,   -24,   -24,   -24,   -24,   -24,
     -24,   -24,   -24,   -24,   -24,   114,   -24,   -24,   -24,   -23,
      19,    23,    63,    64,   -24,   -24,   -24,   -24,    22,   -24,
     -24,   -24,   -24,   -24,   -24,   -24,   -24,   -24,   -24,   -24,
      26,    29,    32,    35,   -24,   -24,   -24,   -24,   -24,   -24,
     -24,    38,    41,    67,    70,    73,    76,    79,    82,    85,
     111,   -24,   -24,   -24,   -24,   -24,   -24,   -24,   -24,   -24,
     -24,   -24,   -24,   -24,   -24
};

/* YYPGOTO[NTERM-NUM].  */
static const yysigned_char yypgoto[] =
{
     -24,   -24,   -24,   -24,   -24,   -24,   -24,   -24,   -24,   -24,
     -24,   -24,   -24,   -24,    86,   -24
};

/* YYTABLE[YYPACT[STATE-NUM]].  What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule which
   number is the opposite.  If zero, do what YYDEFACT says.
   If YYTABLE_NINF, syntax error.  */
#define YYTABLE_NINF -1
static const unsigned char yytable[] =
{
       2,     3,    39,     4,    45,     5,     6,     7,     8,     9,
      10,    11,    12,    13,    14,    15,    16,    17,    18,    19,
      20,    21,    22,    23,    24,    54,    55,    56,    25,    71,
      55,    56,    72,    55,    56,    73,    55,    56,    74,    55,
      56,    75,    55,    56,    76,    55,    56,    46,    57,    58,
      59,    47,    57,    58,    59,    57,    58,    59,    57,    58,
      59,    57,    58,    59,    57,    58,    59,    57,    58,    59,
      77,    55,    56,    78,    55,    56,    79,    55,    56,    80,
      55,    56,    81,    55,    56,    82,    55,    56,    83,    55,
      56,    48,    49,    57,    58,    59,    57,    58,    59,    57,
      58,    59,    57,    58,    59,    57,    58,    59,    57,    58,
      59,    57,    58,    59,    84,    55,    56,    40,    41,    42,
      50,    51,    52,    53,     0,    61,     0,    62,    63,    64,
      65,    66,    67,    68,    69,    70,     0,    57,    58,    59,
       0,    43,    44
};

static const yysigned_char yycheck[] =
{
       0,     1,     5,     3,    27,     5,     6,     7,     8,     9,
      10,    11,    12,    13,    14,    15,    16,    17,    18,    19,
      20,    21,    22,    23,    24,     3,     4,     5,    28,     3,
       4,     5,     3,     4,     5,     3,     4,     5,     3,     4,
       5,     3,     4,     5,     3,     4,     5,    28,    26,    27,
      28,    28,    26,    27,    28,    26,    27,    28,    26,    27,
      28,    26,    27,    28,    26,    27,    28,    26,    27,    28,
       3,     4,     5,     3,     4,     5,     3,     4,     5,     3,
       4,     5,     3,     4,     5,     3,     4,     5,     3,     4,
       5,    28,    28,    26,    27,    28,    26,    27,    28,    26,
      27,    28,    26,    27,    28,    26,    27,    28,    26,    27,
      28,    26,    27,    28,     3,     4,     5,     3,     4,     5,
      34,    35,    36,    37,    -1,    39,    -1,    41,    42,    43,
      44,    45,    46,    47,    48,    49,    -1,    26,    27,    28,
      -1,    27,    28
};

/* YYSTOS[STATE-NUM] -- The (internal number of the) accessing
   symbol of state STATE-NUM.  */
static const unsigned char yystos[] =
{
       0,    30,     0,     1,     3,     5,     6,     7,     8,     9,
      10,    11,    12,    13,    14,    15,    16,    17,    18,    19,
      20,    21,    22,    23,    24,    28,    31,    32,    33,    34,
      35,    36,    37,    38,    39,    40,    41,    42,    43,     5,
       3,     4,     5,    27,    28,    27,    28,    28,    28,    28,
      43,    43,    43,    43,     3,     4,     5,    26,    27,    28,
      44,    43,    43,    43,    43,    43,    43,    43,    43,    43,
      43,     3,     3,     3,     3,     3,     3,     3,     3,     3,
       3,     3,     3,     3,     3
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
      YYPOPSTACK;						\
      goto yybackup;						\
    }								\
  else								\
    {								\
      yyerror (YY_("syntax error: cannot back up")); \
      YYERROR;							\
    }								\
while (0)


#define YYTERROR	1
#define YYERRCODE	256


/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

#define YYRHSLOC(Rhs, K) ((Rhs)[K])
#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)				\
    do									\
      if (N)								\
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
    while (0)
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
} while (0)

# define YY_SYMBOL_PRINT(Title, Type, Value, Location)		\
do {								\
  if (yydebug)							\
    {								\
      YYFPRINTF (stderr, "%s ", Title);				\
      yysymprint (stderr,					\
                  Type, Value);	\
      YYFPRINTF (stderr, "\n");					\
    }								\
} while (0)

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

#if defined (__STDC__) || defined (__cplusplus)
static void
yy_stack_print (short int *bottom, short int *top)
#else
static void
yy_stack_print (bottom, top)
    short int *bottom;
    short int *top;
#endif
{
  YYFPRINTF (stderr, "Stack now");
  for (/* Nothing. */; bottom <= top; ++bottom)
    YYFPRINTF (stderr, " %d", *bottom);
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)				\
do {								\
  if (yydebug)							\
    yy_stack_print ((Bottom), (Top));				\
} while (0)


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

#if defined (__STDC__) || defined (__cplusplus)
static void
yy_reduce_print (int yyrule)
#else
static void
yy_reduce_print (yyrule)
    int yyrule;
#endif
{
  int yyi;
  unsigned long int yylno = yyrline[yyrule];
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %lu), ",
             yyrule - 1, yylno);
  /* Print the symbols being reduced, and their result.  */
  for (yyi = yyprhs[yyrule]; 0 <= yyrhs[yyi]; yyi++)
    YYFPRINTF (stderr, "%s ", yytname[yyrhs[yyi]]);
  YYFPRINTF (stderr, "-> %s\n", yytname[yyr1[yyrule]]);
}

# define YY_REDUCE_PRINT(Rule)		\
do {					\
  if (yydebug)				\
    yy_reduce_print (Rule);		\
} while (0)

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
#  if defined (__GLIBC__) && defined (_STRING_H)
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
static YYSIZE_T
#   if defined (__STDC__) || defined (__cplusplus)
yystrlen (const char *yystr)
#   else
yystrlen (yystr)
     const char *yystr;
#   endif
{
  const char *yys = yystr;

  while (*yys++ != '\0')
    continue;

  return yys - yystr - 1;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined (__GLIBC__) && defined (_STRING_H) && defined (_GNU_SOURCE)
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
static char *
#   if defined (__STDC__) || defined (__cplusplus)
yystpcpy (char *yydest, const char *yysrc)
#   else
yystpcpy (yydest, yysrc)
     char *yydest;
     const char *yysrc;
#   endif
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
      size_t yyn = 0;
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

#endif /* YYERROR_VERBOSE */



#if YYDEBUG
/*--------------------------------.
| Print this symbol on YYOUTPUT.  |
`--------------------------------*/

#if defined (__STDC__) || defined (__cplusplus)
static void
yysymprint (FILE *yyoutput, int yytype, YYSTYPE *yyvaluep)
#else
static void
yysymprint (yyoutput, yytype, yyvaluep)
    FILE *yyoutput;
    int yytype;
    YYSTYPE *yyvaluep;
#endif
{
  /* Pacify ``unused variable'' warnings.  */
  (void) yyvaluep;

  if (yytype < YYNTOKENS)
    YYFPRINTF (yyoutput, "token %s (", yytname[yytype]);
  else
    YYFPRINTF (yyoutput, "nterm %s (", yytname[yytype]);


# ifdef YYPRINT
  if (yytype < YYNTOKENS)
    YYPRINT (yyoutput, yytoknum[yytype], *yyvaluep);
# endif
  switch (yytype)
    {
      default:
        break;
    }
  YYFPRINTF (yyoutput, ")");
}

#endif /* ! YYDEBUG */
/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

#if defined (__STDC__) || defined (__cplusplus)
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
  /* Pacify ``unused variable'' warnings.  */
  (void) yyvaluep;

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
# if defined (__STDC__) || defined (__cplusplus)
int yyparse (void *YYPARSE_PARAM);
# else
int yyparse ();
# endif
#else /* ! YYPARSE_PARAM */
#if defined (__STDC__) || defined (__cplusplus)
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
# if defined (__STDC__) || defined (__cplusplus)
int yyparse (void *YYPARSE_PARAM)
# else
int yyparse (YYPARSE_PARAM)
  void *YYPARSE_PARAM;
# endif
#else /* ! YYPARSE_PARAM */
#if defined (__STDC__) || defined (__cplusplus)
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

  /* Three stacks and their tools:
     `yyss': related to states,
     `yyvs': related to semantic values,
     `yyls': related to locations.

     Refer to the stacks thru separate pointers, to allow yyoverflow
     to reallocate them elsewhere.  */

  /* The state stack.  */
  short int yyssa[YYINITDEPTH];
  short int *yyss = yyssa;
  short int *yyssp;

  /* The semantic value stack.  */
  YYSTYPE yyvsa[YYINITDEPTH];
  YYSTYPE *yyvs = yyvsa;
  YYSTYPE *yyvsp;



#define YYPOPSTACK   (yyvsp--, yyssp--)

  YYSIZE_T yystacksize = YYINITDEPTH;

  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;


  /* When reducing, the number of symbols on the RHS of the reduced
     rule.  */
  int yylen;

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
     have just been pushed. so pushing a state here evens the stacks.
     */
  yyssp++;

 yysetstate:
  *yyssp = yystate;

  if (yyss + yystacksize - 1 <= yyssp)
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
      {
	/* Give user a chance to reallocate the stack. Use copies of
	   these so that the &'s don't force the real ones into
	   memory.  */
	YYSTYPE *yyvs1 = yyvs;
	short int *yyss1 = yyss;


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
	short int *yyss1 = yyss;
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

/* Do appropriate processing given the current state.  */
/* Read a look-ahead token if we need one and don't already have one.  */
/* yyresume: */

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

  /* Shift the look-ahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);

  /* Discard the token being shifted unless it is eof.  */
  if (yychar != YYEOF)
    yychar = YYEMPTY;

  *++yyvsp = yylval;


  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  yystate = yyn;
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
        case 7:
#line 91 "find_dep.y"
    {
        if (strcasecmp((yyvsp[-1].string), "interface") == 0)
            in_interface = true;
        free((yyvsp[-1].string));
    }
    break;

  case 8:
#line 96 "find_dep.y"
    { 
        if (strcasecmp((yyvsp[-3].string), "use") == 0) {
            if (!pp_ignore) {
                DEBUG_PRINT("Use '%s'\n", (yyvsp[-2].string));
                if (!list_find(options.ignore_mods, (yyvsp[-2].string),COMP_FUN(&strcasecmp))) {
                    if (!list_find(dep->modules, (yyvsp[-2].string), COMP_FUN(&strcasecmp)))
                        dep->modules = list_prepend(dep->modules, (yyvsp[-2].string));
                }
            }
        } else if (strcasecmp((yyvsp[-3].string), "module") == 0) {
            if (!pp_ignore && !in_interface) {
                if (!list_find(options.ignore_mods, (yyvsp[-2].string),COMP_FUN(&strcasecmp))) {
                    Module *mod;

                    mod = module_new();
                    mod->sourcefile = xstrdup(sourcefile);
                    mod->modulename = (yyvsp[-2].string);
                    mod->modfile_name = modfile_name((yyvsp[-2].string), mod->sourcefile);

                    if (list_find(modules, mod, &modcmp))
                        warning("Several modules named '%s'", (yyvsp[-2].string));
                    else
                        modules = list_prepend(modules, mod);

                    if (!list_find(dep->targets, mod->modfile_name,
                                   COMP_FUN(&strcasecmp)))
                    dep->targets=list_prepend(dep->targets, mod->modfile_name);
                }
            }
        } else if (strcasecmp((yyvsp[-3].string), "interface") == 0) {
            in_interface = true;
            free((yyvsp[-2].string));
        } else if (strcasecmp((yyvsp[-2].string), "interface") == 0 
                   && strcasecmp((yyvsp[-3].string), "end") == 0) {
            in_interface = false;
            free((yyvsp[-2].string));
        }

        free((yyvsp[-3].string));
    }
    break;

  case 9:
#line 136 "find_dep.y"
    { 
        if (strcasecmp((yyvsp[-3].string), "include") == 0) {
            handle_include((yyvsp[-2].string));
        }
        free((yyvsp[-3].string));
        free((yyvsp[-2].string));
    }
    break;

  case 10:
#line 143 "find_dep.y"
    { 
        handle_include((yyvsp[-2].string));
        free((yyvsp[-2].string));
    }
    break;

  case 12:
#line 148 "find_dep.y"
    { 
        if (!pp_ignore) {
            DEBUG_PRINT("%s defined\n", (yyvsp[-2].string));
            defmac = macro_new();
            macro_setname(defmac, (yyvsp[-2].string));
            if (!list_find(macrolist, defmac, &macrocmp))
                macrolist = list_prepend(macrolist, defmac);
        }
    }
    break;

  case 13:
#line 157 "find_dep.y"
    {
        if (!pp_ignore) {
            Macro *mac;
            List *l;
            
            mac = macro_new();
            macro_setname(mac, (yyvsp[-2].string));

            l = list_find(macrolist, mac, &macrocmp);
            if (l) {
                macrolist = list_remove(macrolist, l);
                macro_free(l->data);
                list_free(l);
            }
            macro_free(mac);
        }
    }
    break;

  case 14:
#line 174 "find_dep.y"
    {
        Macro *mac;

        mac = macro_new();
        macro_setname(mac, (yyvsp[-2].string));

        skip_i++;
        if (pp_ignore)
            pp_ignore++;
        else if (!list_find(macrolist, mac, &macrocmp))
            pp_ignore = 1;
        else
            skip_to_end[skip_i] = true;

        macro_free(mac);
    }
    break;

  case 15:
#line 190 "find_dep.y"
    {
        Macro *mac;

        mac = macro_new();
        macro_setname(mac, (yyvsp[-2].string));

        skip_i++;
        if (pp_ignore)
            pp_ignore++;
        if (list_find(macrolist, mac, &macrocmp))
            pp_ignore = 1;
        else
            skip_to_end[skip_i] = true;

        macro_free(mac);
    }
    break;

  case 16:
#line 206 "find_dep.y"
    { 
        /* #if:s can't be completely ignored, since #else:s, #elif:s and
         * #endif:s aren't.  An #if branch is allways taken, and so are any
         * following #else or #elif:s (ie. no 'skip_to_end'). */
        skip_i++;
        if (pp_ignore)  pp_ignore++;
        skip_to_end[skip_i] = false;
    }
    break;

  case 17:
#line 214 "find_dep.y"
    {
        /* Allways taken unless an #ifdef or #ifndef-branch has been taken
         * allready. */
        if (skip_to_end[skip_i] && pp_ignore == 0)  pp_ignore = 1;
    }
    break;

  case 18:
#line 219 "find_dep.y"
    {
        if (pp_ignore == 1 && skip_to_end[skip_i] == false)  
            pp_ignore = 0;
        else if (pp_ignore == 0)
            pp_ignore = 1;
    }
    break;

  case 19:
#line 225 "find_dep.y"
    {
        skip_to_end[skip_i] = false;
        if (skip_i > 0)  skip_i--;
        if (pp_ignore) pp_ignore--;
    }
    break;

  case 23:
#line 233 "find_dep.y"
    { yyerrok; }
    break;

  case 45:
#line 277 "find_dep.y"
    { free ((yyvsp[0].string)); }
    break;

  case 46:
#line 278 "find_dep.y"
    { free ((yyvsp[0].string)); }
    break;

  case 49:
#line 281 "find_dep.y"
    { 
        if (options.warn_confused) 
            warning ("Unterminated string in file '%s' on line %i", 
                     curr_file, (yyvsp[0].number));
    }
    break;


      default: break;
    }

/* Line 1126 of yacc.c.  */
#line 1500 "y.tab.c"

  yyvsp -= yylen;
  yyssp -= yylen;


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
#if YYERROR_VERBOSE
      yyn = yypact[yystate];

      if (YYPACT_NINF < yyn && yyn < YYLAST)
	{
	  int yytype = YYTRANSLATE (yychar);
	  YYSIZE_T yysize0 = yytnamerr (0, yytname[yytype]);
	  YYSIZE_T yysize = yysize0;
	  YYSIZE_T yysize1;
	  int yysize_overflow = 0;
	  char *yymsg = 0;
#	  define YYERROR_VERBOSE_ARGS_MAXIMUM 5
	  char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];
	  int yyx;

#if 0
	  /* This is so xgettext sees the translatable formats that are
	     constructed on the fly.  */
	  YY_("syntax error, unexpected %s");
	  YY_("syntax error, unexpected %s, expecting %s");
	  YY_("syntax error, unexpected %s, expecting %s or %s");
	  YY_("syntax error, unexpected %s, expecting %s or %s or %s");
	  YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s");
#endif
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
	  int yychecklim = YYLAST - yyn;
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
		yysize_overflow |= yysize1 < yysize;
		yysize = yysize1;
		yyfmt = yystpcpy (yyfmt, yyprefix);
		yyprefix = yyor;
	      }

	  yyf = YY_(yyformat);
	  yysize1 = yysize + yystrlen (yyf);
	  yysize_overflow |= yysize1 < yysize;
	  yysize = yysize1;

	  if (!yysize_overflow && yysize <= YYSTACK_ALLOC_MAXIMUM)
	    yymsg = (char *) YYSTACK_ALLOC (yysize);
	  if (yymsg)
	    {
	      /* Avoid sprintf, as that infringes on the user's name space.
		 Don't have undefined behavior even if the translation
		 produced a string with the wrong number of "%s"s.  */
	      char *yyp = yymsg;
	      int yyi = 0;
	      while ((*yyp = *yyf))
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
	      yyerror (yymsg);
	      YYSTACK_FREE (yymsg);
	    }
	  else
	    {
	      yyerror (YY_("syntax error"));
	      goto yyexhaustedlab;
	    }
	}
      else
#endif /* YYERROR_VERBOSE */
	yyerror (YY_("syntax error"));
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
	  yydestruct ("Error: discarding", yytoken, &yylval);
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
  if (0)
     goto yyerrorlab;

yyvsp -= yylen;
  yyssp -= yylen;
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


      yydestruct ("Error: popping", yystos[yystate], yyvsp);
      YYPOPSTACK;
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  if (yyn == YYFINAL)
    YYACCEPT;

  *++yyvsp = yylval;


  /* Shift the error token. */
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
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
		  yystos[*yyssp], yyvsp);
      YYPOPSTACK;
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif
  return yyresult;
}


#line 288 "find_dep.y"
 

int yyerror(const char *s)
{
    extern int yylineno;

    if (options.warn_confused)
        warning("Line %i in file '%s' confuses me\n", yylineno, curr_file);

    return 0;
}


/* Return false for failure reading file, else return true.  */

bool find_dep(char *file, Dependency *d, List **mods, const List *predef_macro)
{
    extern FILE *yyin;
    extern int yylineno;
    Macro *mac;
    const List *h;

    /* Initialize */
    sourcefile = file;
    curr_file = file;
    dep = d;
    modules = *mods;
    yylineno = 1;
    filestack_i = 0;
    pp_ignore = 0;

    yyin = open_src_file(file, options.src_path);
    if (yyin == NULL) {
        warning("Skipping file '%s': %s", file, strerror(errno));
        return false;
    }

    /* Check source format */
    if (options.src_fmt == SUFFIX)
        lex_set_format(get_format(file));
    else
        lex_set_format(options.src_fmt);

    /* Initialize a list of macros, and fill it with macrodefinitions from -D
     * flags at the command line */
    macrolist = NULL;
    for (h = predef_macro; h; h = h->next) {
        mac = macro_new();
        macro_copy(mac, (Macro *)h->data);
        if (!list_find(macrolist, mac, &macrocmp))
            macrolist = list_prepend(macrolist, mac);
    }

    yyparse();
    *mods = modules;

    /* Delete macrolist */
    for (h = macrolist; h; h = h->next)
        macro_free((Macro *)h->data);
    list_free(macrolist);

    fclose(yyin);

    return true;
}
            

int modstrcmp(const void *s, const void *m)
{
    return strcasecmp(((Module *)m)->modulename, (char *)s);
}



static int modcmp (const void *m1, const void *m2)
{
    return strcasecmp(((Module *)m1)->modulename, ((Module *)m2)->modulename);
}



static char *fixed_suffixes[] = {".f", ".F", ".for", ".FOR", ".ftn", ".FTN"};

static char *free_suffixes[]  = {".f90", ".F90", ".f95", ".F95"};
                                  
static const int fixsuffn = sizeof(fixed_suffixes)/sizeof(void *);
static const int freesuffn = sizeof(free_suffixes)/sizeof(void *);


SourceFmt get_format(const char *filename)
{
    const char *p;
    int i;

    /* Search for the end */
    for (p = filename; *p; p++);

    /* Search backwards for the last '.' */
    for (; *p != '.' && p != filename; p--);

    /* Check for any of the free suffixes */
    for (i = 0; i < freesuffn; i++) 
        if (strcmp(p, free_suffixes[i]) == 0)  return FREE;

    /* Check for any of the fixed suffixes */
    for (i = 0; i < fixsuffn; i++)
        if (strcmp(p, fixed_suffixes[i]) == 0)  return FIXED;

    return UNKNOWN;
}


void pop_filestack()
{
    curr_file = filestack[--filestack_i];
}


void handle_include(const char *incfile)
{
    if (!pp_ignore) {
        filestack[filestack_i++] = curr_file;
        curr_file = remove_citation(incfile);
        if (lex_include_file(curr_file))  {
            DEBUG_PRINT("including file '%s'\n", curr_file);
            if (!list_find(dep->includes, curr_file, COMP_FUN(&strcasecmp)))
                dep->includes = list_prepend(dep->includes, curr_file);
        } else {
            pop_filestack();
        }
    }
}


Module *module_new()
{
    Module *m;

    m = (Module *) xmalloc(sizeof(Module));
    m->modulename = m->modfile_name = m->sourcefile = NULL;
    return m;
}


Dependency *dependency_new()
{
    Dependency *d;

    d = (Dependency *) xmalloc(sizeof(Dependency));
    d->sourcefile = NULL;
    d->targets = d->modules = d->includes = NULL;
    return d;
}

