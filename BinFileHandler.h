#ifndef _BIN_FILE_HANDLER_H_
#define _BIN_FILE_HANDLER_H_

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <limits.h>   // PATH_MAX
#include <string.h>   // strdup
#include <libgen.h>   // dirname()

#define MAXLINE 500

//------------------------------------------------------------------------------

class BinFileHandler {
public:
#if defined(__sgi) 
  typedef long long OffType; 
  inline OffType fseek(FILE *fp, OffType offset, int whence)  { return fseek64(fp, offset, whence); }
  inline OffType ftell(FILE *fp)  { return ftell64(fp); }
#elif defined(linux)
  typedef __off64_t  OffType;
  inline OffType fseek(FILE *fp, OffType offset, int whence)  { return fseeko64(fp, offset, whence); }
  inline OffType ftell(FILE *fp)  { return ftello64(fp); }
#elif defined(__sun)
  typedef long OffType; 
#elif defined(__APPLE__) && defined(__MACH__)
  typedef off_t OffType;
  inline OffType fseek(FILE *fp, OffType offset, int whence)  { return fseeko(fp, offset, whence); }
  inline OffType ftell(FILE *fp)  { return ftello(fp); }
#else
#error Update the definition of OffType for your machine
#endif

private:

  OffType cpos;

  int headersize;

  double version;

  bool swapBytes;

  int fileid;

  FILE *file;
    void makepath(const char * filename);

public:

  BinFileHandler(const char *, const char *, double = 0.0);
  ~BinFileHandler();

  template<class Scalar>
  void swapVector(Scalar *, int);

  template<class Scalar>
  void read(Scalar *, int);

  template<class Scalar>
  void write(Scalar *, int);

  void seek(OffType);

  OffType tell();

  double getVersion() const { return version; }

};

//------------------------------------------------------------------------------

template <class Scalar>
void BinFileHandler::read(Scalar *p, int nobjs)
{
  int len;
  if (file)
    len = fread(p, sizeof(Scalar), nobjs, file);
  else len = ::read(fileid, p, nobjs*sizeof(Scalar));

  cpos += nobjs*sizeof(Scalar);

  if (swapBytes) swapVector(p, nobjs);

}

//------------------------------------------------------------------------------

template <class Scalar>
void BinFileHandler::write(Scalar *p, int nobjs)
{

  if (swapBytes) swapVector(p, nobjs);
  int len;
  if (file) len = fwrite(p, sizeof(Scalar), nobjs, file);
  else len = ::write(fileid, p, nobjs*sizeof(Scalar));

  cpos += nobjs*sizeof(Scalar);

  if (swapBytes) swapVector(p, nobjs);

}

//------------------------------------------------------------------------------

template <class Scalar>
void BinFileHandler::swapVector(Scalar *p, int nobjs)
{

  for (int obj = 0; obj < nobjs; ++obj) {

    Scalar x = p[obj];

    char *px = (char *) &x;
    char *pp = (char *) (p+obj);

    for (int c = 0; c < sizeof(Scalar); ++c)
      pp[sizeof(Scalar)-1-c] = px[c];

  }

}


//------------------------------------------------------------------------------

inline
void BinFileHandler::seek(BinFileHandler::OffType size) 
{ 

  size += headersize;

  if (file) fseek(file, size, SEEK_SET);
  else lseek(fileid, size, SEEK_SET);

  cpos = size;

}


//------------------------------------------------------------------------------

inline
BinFileHandler::OffType BinFileHandler::tell() 
{ 

  BinFileHandler::OffType pos;

  if (file) 
    pos = ftell(file);
#ifdef sgi
  else 
    pos = ::tell(fileid);
#else
  else
    pos = cpos;
#endif

  pos -= headersize;

  return pos;

}

//------------------------------------------------------------------------------

inline
BinFileHandler::BinFileHandler(const char *name, const char *flag, double ver) :
  file(0), fileid(0), swapBytes(0), version(ver) 
{

  // lei lei, 28 March 2016: calls mkdir() recursively if dirname(name) does not exist
  makepath(name);

  int ierr = 0;

  if (strcmp(flag, "r") == 0) {
    fileid = open(name, O_RDONLY, 0644);
    if (fileid == -1) ierr = 1;
  }
  else if (strcmp(flag, "w") == 0) {
    fileid = open(name, O_WRONLY | O_CREAT | O_TRUNC, 0644);
    if (fileid == -1) ierr = 1;
  }
  else if (strcmp(flag, "ws") == 0) {
    fileid = open(name, O_WRONLY | O_CREAT | O_TRUNC | O_SYNC, 0644);    
    if (fileid == -1) ierr = 1;
  }
  else if (strcmp(flag, "w+") == 0) {
    fileid = open(name, O_WRONLY | O_CREAT, 0644);
    if (fileid == -1) ierr = 1;
  }
  else if (strcmp(flag, "ws+") == 0) {
    fileid = open(name, O_WRONLY | O_CREAT | O_SYNC, 0644);
    if (fileid == -1) ierr = 1;
  }
  else if (strcmp(flag, "rb") == 0) {
    file = fopen(name, "rb");
    if (!file) ierr = 1;
  }
  else if (strcmp(flag, "wb") == 0) {
    file = fopen(name, "wb");
    if (!file) ierr = 1;
  }
  else {
    fprintf(stderr, "*** Error: wrong flag (%s) for \'%s\'\n", flag, name);
    exit(1);
  }
  
 if (ierr) {
   fprintf(stderr, "*** Error: unable to open \'%s\' with flag (%s)\n", name, flag);
    exit(1);
  }

  headersize = sizeof(int) + sizeof(double);
    
  int one = 1;
  if (strcmp(flag, "r") == 0 || strcmp(flag, "rb") == 0) {
    read(&one, 1);
    if (one != 1) swapBytes = 1;
    read(&version, 1);
  } 
  else if (strcmp(flag, "w") == 0 || strcmp(flag, "ws") == 0 || strcmp(flag, "wb") == 0) {
    write(&one, 1);
    write(&version, 1);
  }
  else if (strcmp(flag, "w+") == 0 || strcmp(flag, "ws+") == 0)
    seek(0);

  cpos = headersize;

}

//------------------------------------------------------------------------------

inline
BinFileHandler::~BinFileHandler() 
{ 

  if (file) fclose(file); 
  else close(fileid);

}

//------------------------------------------------------------------------------

inline
int computeNumberOfDigits(int num)
{

  int digits = 1;

  while (num >= 10) {
    num /= 10;
    ++digits;
  }

  return digits;

}

//------------------------------------------------------------------------------

inline
char *computeClusterSuffix(int num, int maxNum)
{

  int numZeros = computeNumberOfDigits(maxNum) - computeNumberOfDigits(num);

  char zeros[100];
  char *suffix = new char[100];

  strcpy(zeros, "");
  for (int k=0; k<numZeros; ++k)
    strcat(zeros, "0");

  sprintf(suffix, "%s%d", zeros, num);

  return suffix;

}

//------------------------------------------------------------------------------

//TODO: add error handling code later
inline
void BinFileHandler::makepath(const char *filename) {
  char *s = strdup(filename);
  char *path = dirname(s);

  char tmp[PATH_MAX];
  char *p = NULL;
  size_t len;

  snprintf(tmp, sizeof(tmp), "%s", path);
  len = strlen(tmp);
  if (tmp[len - 1] == '/')
    tmp[len - 1] = '\0';
  for (p = tmp + 1; *p; p++) {
    if (*p == '/') {
      *p = '\0';
      mkdir(tmp, S_IRWXU);
      *p = '/';
    }
  }
  mkdir(tmp, S_IRWXU);
  free(s);
}

#endif
