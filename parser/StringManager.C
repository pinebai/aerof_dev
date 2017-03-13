#include "StringManager.h"

StringManager theStringManager;

StringManager::~StringManager()
{
  for(std::list<char*>::iterator it=strings.begin(); it != strings.end();++it)
    {
      delete[] *it;
    }
}

char * StringManager::stringNew(int bufferLength)
{
  char * result = new char[bufferLength];
  strings.push_back(result);
  return result;
}
