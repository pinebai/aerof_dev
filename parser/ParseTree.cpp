/*
 * ParseTree.cpp
 *
 *  Created on: Feb 10, 2009
 *      Author: michel
 */
#include <iostream>

#include "ParseTree.h"
#include "Dictionary.h"

ParseTree::ParseTree(const char *name) : Assigner(name) {

}

ParseTree::~ParseTree() {

}

void ParseTree::addInt(Token token, int val) {
 intMap[token] = val;
}
void ParseTree::addDouble(Token token, double val) {
  doubleMap[token] = val;
}
void ParseTree::addToken(Token token, int val) {
  tokenMap[token] = val;
}
void ParseTree::addString(Token token, const char *val) {
  stringMap[token] = val;
}

Assigner *ParseTree::findSubToken(Token token) {
  /*map<Token,ParseTree*>::iterator it = subTokenMap.find(token);
  ParseTree *result;
  if(it == subTokenMap.end()) {
    result = new ParseTree();
    subTokenMap[token] = result;
  }
  else
    result = it->second;
  return result;*/

  map<Token,ParseNode*>::iterator it = nodeMap.find(token);
  ParseNode *result;
  if(it == nodeMap.end()) {
    result = new ParseNode(*this, token);
    nodeMap[token] = result;
  }
  else
    result = it->second;
  return result;

}

Assigner *ParseTree::getSubToken(Token token, Token subToken)
{
  map<Token,ParseTree*>::iterator it = subTokenMap.find(token);
  ParseTree *result;
  if(it == subTokenMap.end()) {
    result = new ParseTree(dictionary->word(token).c_str());
    subTokenMap[token] = result;
  }
  else
    result = it->second;
  return result->findSubToken(subToken);
}

Assigner *ParseTree::findIndexObject(int index) {
 /* map<Token,ParseTree*>::iterator it = subTokenMap.find(index);
  ParseTree *result;
  if(it == subTokenMap.end()) {
    result = new ParseTree();
    subTokenMap[index] = result;
  }
  else
    result = it->second;
  return result;*/
  return 0;
}

ParseNode::ParseNode(ParseTree &t, int tk) :
  Assigner(dictionary->word(tk).c_str()),
  tree(t), token(tk) {}

Assigner *ParseNode::findIndexObject(int tk) { return 0; }

/** propagate the parse tree to an actual object */
void ParseTree::implement(Assigner *assigner) {
  for(map<Token,ParseTree*>::iterator it = subTokenMap.begin();
      it != subTokenMap.end() ; ++it) {
    Assigner *subAssigner = assigner->findSubToken(it->first);
    if(subAssigner)
      it->second->implement(subAssigner);
    else
      std::cerr << "Did not find a token: " << dictionary->word(it->first).c_str()<< std::endl;
  }

  for(map<Token,int>::iterator it = intMap.begin();
        it != intMap.end() ; ++it) {
    Assigner *subAssigner = assigner->findSubToken(it->first);
    if(subAssigner)
      subAssigner->assignInt(it->second);
    else
          std::cerr << "Did not find a token!" << dictionary->word(it->first).c_str()<< std::endl;
  }
  for(map<Token,double>::iterator it = doubleMap.begin();
        it != doubleMap.end() ; ++it) {
    Assigner *subAssigner = assigner->findSubToken(it->first);
    if(subAssigner)
      subAssigner->assignDouble(it->second);
    else
          std::cerr << "Did not find a token!" << dictionary->word(it->first).c_str() << std::endl;
  }
  for(map<Token,Token>::iterator it = tokenMap.begin();
        it != tokenMap.end() ; ++it) {
    Assigner *subAssigner = assigner->findSubToken(it->first);
    if(subAssigner)
      subAssigner->assignInt(it->second);
    else
          std::cerr << "Did not find a token!" << dictionary->word(it->first).c_str() << std::endl;
  }
  for(map<Token,const char *>::iterator it = stringMap.begin();
        it != stringMap.end() ; ++it) {
    Assigner *subAssigner = assigner->findSubToken(it->first);
    if(subAssigner)
      subAssigner->assignString(it->second);
    else
          std::cerr << "Did not find a token!" << dictionary->word(it->first).c_str() << std::endl;
  }
}
