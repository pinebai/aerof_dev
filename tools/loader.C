#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <dlfcn.h>

#ifdef USE_MPI
#include <mpi.h>
#endif

#define MAXLINE 500

int processInputFile(char *, int *&, int *&, char ***&);
int processInputLine(int, char **, int *&, int *&, char ***&);
void decodeArgumentLine(char *, int &, int &, char **&);
int startCode(int, char **, const char * = "entrypoint");

//------------------------------------------------------------------------------

#ifdef USE_MPI
#ifndef MPI_INTEGER
#define MPI_INTEGER MPI_INT
#endif
void dummy()
{

  MPI_Comm comm1;
  MPI_Abort(MPI_COMM_WORLD, 3);
  MPI_Comm_split(MPI_COMM_WORLD, 1, 1, &comm1);
  MPI_Allreduce(0, 0, 0, MPI_INTEGER, MPI_MAX, comm1);
  MPI_Intercomm_create(comm1, 0, comm1, 1, 1, &comm1);
  MPI_Comm_remote_size(comm1, 0);
  MPI_Barrier(comm1);
  MPI_Waitall(0, 0, 0);
  MPI_Isend(0, 0, MPI_INTEGER, 0, 0, comm1, 0);
  MPI_Recv(0, 0, MPI_INTEGER, MPI_ANY_SOURCE, 0, comm1, 0);
  MPI_Get_count(0, MPI_INTEGER, 0);
  MPI_Irecv(0, 0, MPI_INTEGER, 0, 0, comm1, 0);
  MPI_Bcast(0, 0, MPI_INTEGER, 0, comm1);

}
#endif

//------------------------------------------------------------------------------

int main(int argc, char **argv)
{

  int size = 1, rank = 0;

#ifdef USE_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);
#endif

  int i, numCodes;
  int *size1, *argc1;
  char ***argv1;
  int ierr;

  if (argv[1] == 0) {
    if (rank == 0) {
      fprintf(stderr, "Incorrect usage\n");
      fprintf(stderr, "Usage: %s entry [, entry ...]\n", argv[0]);
      fprintf(stderr, "       loader -f <filename>\n");
    }
    exit(-1);
  }
  else if (strncmp(argv[1], "-f", 2) == 0)
    numCodes = processInputFile(argv[2], size1, argc1, argv1);
  else
    numCodes = processInputLine(argc-1, argv+1, size1, argc1, argv1);

  /*
  if (rank == 0) {
    for (i=0; i<numCodes; ++i) {
      fprintf(stderr, "code %d has %d argument(s) and uses %d cpu(s):\n", 
	      i, argc1[i], size1[i]);
      for (int j=0; j<argc1[i]; ++j)
	fprintf(stderr, "  \'%s\'\n",argv1[i][j]);
    }
  }
  */

  int numCpus = 0;
  for (i=0; i<numCodes; ++i)
    numCpus += size1[i];

  if (numCpus != size) {
    if (rank == 0)
      fprintf(stderr, "*** Error: incorrect number of cpus (%d vs %d)\n", numCpus, size);
    exit(-1);
  }

  int *index = new int[numCodes + 1];
  index[0] = 0;
  for (i=1; i<numCodes+1; ++i)
    index[i] = index[i-1] + size1[i-1];

//  for (i=0; i<numCodes; ++i)
//    if (rank >= index[i] && rank < index[i+1])
//      startCode(argc1[i], argv1[i]);
  for (i=0; i<numCodes; ++i)
    if (rank >= index[i] && rank < index[i+1]){
      ierr = startCode(argc1[i], argv1[i]);
      if(ierr) ierr = startCode(argc1[i], argv1[i], "entrypoint_");
      if(ierr) exit(-1);
    }

#ifdef USE_MPI
  MPI_Finalize();
#endif
  
  return 0;

}

//------------------------------------------------------------------------------

int processInputFile(char *name, int *&size, int *&argc1, char ***&argv1)
{

  char line[MAXLINE];

  FILE *fp = fopen(name, "r");
  if (!fp) {
    fprintf(stderr, "*** Error: could not open \'%s\'\n", name);
    exit(-1);
  }

  int numCodes = 0;

  while (fgets(line, MAXLINE, fp) != 0)
    ++numCodes;

  fclose(fp);

  size = new int[numCodes];
  argc1 = new int[numCodes];
  argv1 = new char **[numCodes];

  fp = fopen(name, "r");

  numCodes = 0;

  while (fgets(line, MAXLINE, fp) != 0) {
    decodeArgumentLine(line, size[numCodes], argc1[numCodes], argv1[numCodes]);
    ++numCodes;
  }

  fclose(fp);

  return numCodes;

}

//------------------------------------------------------------------------------

int processInputLine(int argc, char **argv, int *&size, int *&argc1, char ***&argv1)
{

  int numCodes = 1;

  int i;
  for (i=0; i<argc; ++i)
    if (strcmp(argv[i], ",") == 0)
      ++numCodes;

  int *index = new int[numCodes];
  size = new int[numCodes];
  argc1 = new int[numCodes];
  argv1 = new char **[numCodes];

  int codeNum = 0;
  argc1[0] = 0;
  index[0] = 0;

  for (i=0; i<argc; ++i) {
    if (strcmp(argv[i], ",") == 0) {
      ++codeNum;
      argc1[codeNum] = 0;
      index[codeNum] = i+1;
    }
    else
      ++argc1[codeNum];
  }

  for (i=0; i<numCodes; ++i) {
    argv1[i] = new char *[ argc1[i] ];
    argc1[i] -= 1;
    argv1[i][ argc1[i] ] = 0;
  }

  for (i=0; i<numCodes; ++i) {
    size[i] = atoi(argv[ index[i] ]);
    for (int j=0; j<argc1[i]; ++j) {
      argv1[i][j] = new char[strlen(argv[ index[i]+j+1 ]) + 1];
      sprintf(argv1[i][j], "%s", argv[ index[i]+j+1 ]);
    }
  }
    
  if (index) delete [] index;

  return numCodes;

}

//------------------------------------------------------------------------------

void decodeArgumentLine(char *line, int &size, int &argc, char **&argv)
{

  if (line == 0) {
    fprintf(stderr, "*** Error: input line incorrect\n");
    exit(-1);
  }

  int i, len = strlen(line);

  char s[MAXLINE];

  i = 0;
  argc = 0;

  do {
    if (line[i] != ' ' ) {
      sscanf(line + i, "%s", s);
      ++argc;
      i += strlen(s);
    }
    else
      ++i;
  } while (i != len - 1);

  argv = new char *[argc];
  argc -= 1;
  argv[argc] = 0;

  i = 0;
  argc = -1;

  do {
    if (line[i] != ' ' ) {
      sscanf(line + i, "%s", s);

      if (argc == -1) {
	size = atoi(s);
	argc = 0;
      } else {
	argv[argc] = new char[strlen(s) + 1];
	sprintf(argv[argc], "%s", s);
	++argc;
      }

      i += strlen(s);
    }
    else
      ++i;
  } while (i != len - 1);

}

//------------------------------------------------------------------------------

int startCode(int argc, char **argv, const char *routine)
{

  char *name = argv[0];

  void *handle = dlopen(name, RTLD_NOW);

  char *msg = dlerror();

  if (msg) {
    fprintf(stderr,"*** Error: dynamic loading of \'%s\': %s\n", name, msg);
    dlclose(handle);
    return 1;
    //exit(-1);
  }

  int (*fct)(int, char **) = (int (*)(int, char **)) dlsym(handle, routine);

  if (!fct) {
    fprintf(stderr,"*** Error: could not find \'%s\' in \'%s\'\n", routine, name);
    dlclose(handle);
    return 1;
    //exit(-1);
  }

  int ierr = (*fct)(argc, argv);

  return ierr;

}

//------------------------------------------------------------------------------
/*
#ifdef USE_MPI
void exit(int status)
{

  MPI_Abort(MPI_COMM_WORLD, status);

}
#endif
*/
//------------------------------------------------------------------------------
