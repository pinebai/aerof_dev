#ifndef _STRING_MANAGER_H_
#define _STRING_MANAGER_H_
#include <list>

class StringManager
{
 public:
  char *stringNew(int bufferLength);
  ~StringManager();
  StringManager(){}
 private:
  std::list<char*> strings;
  StringManager(const StringManager &);
  StringManager& operator=(const StringManager &);
};

extern StringManager theStringManager;

#endif /*_STRING_MANAGER_H_*/
