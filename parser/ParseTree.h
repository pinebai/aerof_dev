/*
 * ParseTree.h
 *
 *  Created on: Feb 10, 2009
 *      Author: Michel Lesoinne
 */

#ifndef PARSETREE_H_
#define PARSETREE_H_

#include <map>

#include "Assigner.h"
typedef int Token;
class ParseNode;

class ParseTree : public Assigner {
    std::map<Token,int> intMap;
    std::map<Token,double> doubleMap;
    std::map<Token,Token> tokenMap;
    std::map<Token,const char *> stringMap;
    std::map<Token,ParseTree *> subTokenMap;
    std::map<int,ParseTree *> subIndexedMap;

    std::map<Token, ParseNode *> nodeMap;
  public:
    ParseTree(const char *name);
    virtual ~ParseTree();


    Assigner *findSubToken(int);

    Assigner *findIndexObject(int);

    void addInt(Token token, int val);
    void addDouble(Token token, double val);
    void addToken(Token token, int val);
    void addString(Token token, const char *val);
    Assigner *getSubToken(Token token, Token subToken);

    void implement(Assigner *);
};

class ParseNode : public Assigner {
    int token;
    ParseTree &tree;
  public:
    ParseNode(ParseTree &t, int tk);
    void assignInt(int val) {
      tree.addInt(token, val);
    }
    void assignDouble(double val) {
      tree.addDouble(token, val);
    }
    void assignToken(int val) {
      tree.addToken(token, val);
    }
    void assignString(const char *val) {
      tree.addString(token, val);
    }

    Assigner *findSubToken(Token tk) {
      return tree.getSubToken(token, tk);
    }

    Assigner *findIndexObject(int tk);
};

template <class T>
class ClassParseTree  : public Assigner {
  T *ptr;
  ParseTree T::*tree;
public:
    ClassParseTree(ClassAssigner *, const char *n, T *ptr, ParseTree T::*sp);
    Assigner *findSubToken(int);
};

template <class T>
ClassParseTree<T>::ClassParseTree(ClassAssigner *ca, const char *n,
                                  T *_ptr, ParseTree T::*_tree)
: Assigner(n)
{
 ptr = _ptr;
 tree = _tree;
 ca->addSmb(n, this);
}

template <class T>
Assigner *
ClassParseTree<T>::findSubToken(int token) {
  ParseTree &t = ptr->*tree;
  return t.findSubToken(token);
}

#endif /* PARSETREE_H_ */
