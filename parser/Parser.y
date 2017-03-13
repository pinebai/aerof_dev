%{
#include <math.h>
#include <stdio.h>
#include <string.h>
extern int  yyCmdflex(void);
extern void yyCmdferror(const char  *);
#include "parser/Assigner.h"
#include "parser/Dictionary.h"
#include "ResizeArray.h"

void reduceList(int, int *(list[]), ... );

ResizeArray<Assigner *> underVar(0);
int underLevel = 0;

%}

%union
{
  int       ival;
  double    dval;
  int	    token;
  char     *sval; 
  int      **list;
  Assigner *asgn;
}

%token UNDER EoF IntConstant DblConstant Symbol String END BOGUS

%left '+' '-'
%left '*' '/'
%left '^'

%type <ival>	IntConstant IntExpr
%type <dval>    DblConstant DblExpr
%type <token>   Symbol
%type <asgn>	Assignable
%type <sval>	String
%type <list>	List1
%type <list>	List2
%type <list>	List3
%type <list>	List4
%type <list>	List5
%type <list>	List6
%type <list>	List7

%%

File: ValidInput EoF
          { return 0; }
	| ValidInput END
	{ return 0; }

ValidInput: 
	| ValidInput Assignment
	| ValidInput GroupInput

Assignment: Assignable '=' Symbol ';'
	{ $1->assignToken($3); }
	| Assignable '=' IntExpr ';'
	{ $1->assignInt($3); }
	| Assignable '=' DblExpr ';'
	{ $1->assignDouble($3); }
	| Assignable '=' String ';'
	{ $1->assignString($3); }
	| Assignable '=' Symbol IntExpr ';'
	{ $1->assignTokenIntPair($3,$4); }
        | Assignable '=' List1 ';'
	{ $1->assignList(1,$3); }
        | Assignable '=' List2 ';'
	{ $1->assignList(2,$3); }
        | Assignable '=' List3 ';'
	{ $1->assignList(3,$3); }
        | Assignable '=' List4 ';'
	{ $1->assignList(4,$3); }
        | Assignable '=' List5 ';'
	{ $1->assignList(5,$3); }
        | Assignable '=' List6 ';'
	{ $1->assignList(6,$3); }
        | Assignable '=' List7 ';'
	{ $1->assignList(7,$3); }

GroupInput: UNDER Assignable '{'
	{
          underVar[underLevel] = $2; underLevel++;
       }
	 ValidInput '}'
	{ underLevel--; }


Assignable:
	Symbol
	{ 
           $$ = 0;
           for(int i = underLevel; i-- ; ) {
             $$ = underVar[i]->findSubToken($1);
             if($$ != 0)  
                break; 
           } 
           if($$ == 0)
             $$ = findSysObj($1); 
        }
	| Assignable '.' Symbol
	{ 
          $$ = $1->findSubToken($3); 
          if($$ == 0) {
	     fprintf(stderr, "ERROR: Structure element not found: %s\n",
                  dictionary->word($3).c_str());
             exit(-1);
          }
        }
        | Assignable '[' IntExpr ']'
        { $$ = $1->findIndexObject($3); 
          if($$ == 0) {
            fprintf(stderr, "ERROR: Object is not an array\n");
            exit(-1);
          }
        }
        


IntExpr:
	IntConstant
	| '(' IntExpr ')'
	{ $$ = $2; }
	| IntExpr '+' IntExpr
	{ $$ = $1 + $3; }
	| IntExpr '-' IntExpr
	{ $$ = $1 - $3; }
	| IntExpr '*' IntExpr
	{ $$ = $1 * $3; }
	| IntExpr '/' IntExpr
	{ $$ = $1 / $3; }

DblExpr:
	DblConstant
	| '(' DblExpr ')'
        { $$ = $2; }
        | DblExpr '+' DblExpr
        { $$ = $1 + $3; }
        | DblExpr '+' IntExpr
	{ $$ = $1 + $3; }
        | IntExpr '+' DblExpr
	{ $$ = $1 + $3; }
        | DblExpr '-' DblExpr
        { $$ = $1 - $3; }
        | DblExpr '-' IntExpr
        { $$ = $1 - $3; }
        | IntExpr '-' DblExpr
        { $$ = $1 - $3; }
        | DblExpr '*' DblExpr
        { $$ = $1 * $3; }
        | IntExpr '*' DblExpr
        { $$ = $1 * $3; }
        | DblExpr '*' IntExpr
        { $$ = $1 * $3; }
        | DblExpr '/' DblExpr
        { $$ = $1 / $3; }
        | IntExpr '/' DblExpr
        { $$ = $1 / $3; }
        | DblExpr '/' IntExpr
        { $$ = $1 / $3; }

List1: '{' Symbol '}'
	{ 
          int size = 1;
          int *list = NULL;
          reduceList(size, &list,$2);
          $$ = &list;
	}

List2: '{' Symbol ',' Symbol '}'
	{ 
          int size = 2;
          int *list = NULL;
          reduceList(size, &list,$2,$4);
          $$ = &list;
	}

List3: '{' Symbol ',' Symbol ',' Symbol '}'
	{ 
          int size = 3;
          int *list = NULL;
          reduceList(size, &list,$2,$4,$6);
          $$ = &list;
	}

List4: '{' Symbol ',' Symbol ',' Symbol ',' Symbol '}'
	{ 
          int size = 4;
          int *list = NULL;
          reduceList(size, &list,$2,$4,$6,$8);
          $$ = &list;
	}

List5: '{' Symbol ',' Symbol ',' Symbol ',' Symbol ',' Symbol '}'
	{ 
          int size = 5;
          int *list = NULL;
          reduceList(size, &list,$2,$4,$6,$8,$10);
          $$ = &list;
	}

List6: '{' Symbol ',' Symbol ',' Symbol ',' Symbol ',' Symbol ',' Symbol '}'
	{ 
	  int size = 6;
	  int *list = NULL;
	  reduceList(size, &list,$2,$4,$6,$8,$10,$12); 
	  $$ = &list;
	}
List7: '{'  Symbol ',' Symbol ',' Symbol ',' Symbol ',' Symbol ',' Symbol ',' Symbol '}'
	{ 
	  int size = 7;
	  int *list = NULL;
	  reduceList(size, &list, $2,$4,$6,$8,$10,$12,$14); 
	  $$ = &list;
	}

%%

void reduceList(int symNum, int *(list[]), ...) {
  int tk = 0;
  va_list syms;
  va_start(syms,list);
  *list = (int *) malloc(symNum*sizeof(int));
  if (*list == NULL){
    fprintf(stderr, "Error allocating memory for variable set\n");
    exit(-1);
  }
  memset(*list,0,symNum*sizeof(int));
  for(int i = 0; i < symNum; i++){
    tk = va_arg(syms,int);
    (*list)[i] = tk;
  }
  va_end(syms);
}

