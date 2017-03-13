#include <cstdio>
#include "Dictionary.h"
#include "Assigner.h"

SysSmbMap *sysSmb = 0;
Dictionary *dictionary = 0;

SysSmbMap::~SysSmbMap()
{
  for(map<int, Assigner *>::iterator it=forest.begin();it!=forest.end();++it)
    {
      delete it->second;
    }
}

int addSysSymbol(const char *name, Assigner *a)
{
  static SysSmbMap theSysSmb;
  sysSmb =  &theSysSmb;
  map<int, Assigner *>::iterator it = sysSmb->forest.find(findSysToken(name));
  if( it != sysSmb->forest.end()) delete it->second;
  sysSmb->forest[findSysToken(name)] = a;
  return 1;
}

int findSysToken(const char *str)
{
  static Dictionary theDictionary;
  dictionary = &theDictionary;
  return dictionary->token(str);
}

Assigner *
findSysObj(int tk)
{
  map<int, Assigner *>::iterator it = sysSmb->forest.find(tk);
  if(it == sysSmb->forest.end()) {
    fprintf(stderr, "Error: Symbol not found: %s\n", dictionary->word(tk).c_str());
    // Exit so that AERO-f does not seg-fault.
    exit(-1);
    //return 0;
  }
  return it->second;
}
