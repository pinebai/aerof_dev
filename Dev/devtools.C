/*
 * devtools.C
 *
 *  Created on: Nov 16, 2016
 *      Author: Lukas Scheucher
 *      Mail:   lscheuch@stanford.edu
 */
#include "./devtools.h"
#include <stdlib.h>

//#include <stdio.h>
#include <iostream>
#include <fstream> //file io

#include <execinfo.h>

//declare a NULL pointer, such that no extra header has to be included just for that putpose
#ifndef NULL
#define NULL   ((void *) 0)
#endif

#define BUFFSIZE 100

//trigger colored output to terminal
//You might want to deactivate this if you plan to redirect terminal output to files
#define COLORS

#define RED "\033[91m"     //escape sequence for red    font on black background
#define YELLOW "\033[93m"  //escape sequence for yellow font on black background
#define WHITE "\033[00m"   //escape sequence for white  font on black background

namespace Dev{

/*****************************************************************************
 *  The purpose of this function is to print more useful warning messages    *
 *  to the screen. It utilized the printf function of the Communicator class,*
 *  so it is ensured that only one processor prints the output               *
 *                                                     Lukas Scheucher 11/16 *
 *****************************************************************************/
void Warning(Communicator *com,const char* message,bool printBacktrace)
{
  //Print the warning message
#ifdef COLORS
  com->printf(0,YELLOW);
#endif
  com->printf(0,"WARNING!: ");
  com->printf(0,message);
#ifdef COLORS
  com->printf(0,WHITE);
#endif
  com->printf(0,"\n");

  if (printBacktrace==true)
    printStack(com);

}

//------------------------------------------------------------------------------

//TODO
void Warning(MPI_Comm comm, const char* message,bool printBacktrace)
{
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_rank(comm, &size);

//Print the warning message
if (rank==0)
{
#ifdef COLORS
  std::cout<<YELLOW;
#endif
  std::cout<<"WARNING!: ";
  std::cout<<message;
#ifdef COLORS
  std::cout<<WHITE;
#endif
  std::cout<<"\n";

  if (printBacktrace==true)
    printStack(comm);//printStack function with MPI communicator
}

}


/*****************************************************************************
 *  The purpose of this function is to print more useful error messages    *
 *  to the screen. It utilized the printf function of the Communicator class,*
 *  so it is ensured that only one processor created the output              *
 *                                                     Lukas Scheucher 11/16 *
 *****************************************************************************/
void Error(Communicator *com,const char* message,bool printBacktrace)
{
#ifdef COLORS
  com->printf(0,RED);
#endif
  com->printf(0,"ERROR!: ");
  com->printf(0,message);
#ifdef COLORS
  com->printf(0,WHITE);
#endif
  com->printf(0,"\n");

  if (printBacktrace==true)
    printStack(com);

  //Exit program execution
  exit(-1);

}

//------------------------------------------------------------------------------

void Error(MPI_Comm comm, const char* message,bool printBacktrace)
{
  int rank, size;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_rank(comm, &size);

#ifdef COLORS
  if (rank==0) fprintf(stdout,RED);
#endif
  if (rank==0) fprintf(stdout,"ERROR!: ");
  if (rank==0) fprintf(stdout,message);
#ifdef COLORS
  if (rank==0) fprintf(stdout,WHITE);
#endif
  if (rank==0) fprintf(stdout,"\n");

  if (printBacktrace==true)
    printStack(comm);

  //Exit program execution
  sleep(2000);
  exit(-1);

}

//------------------------------------------------------------------------------

//printing the stacktrace
void printStack(Communicator *com)
{

  //get the stacktrace size
  void *buffer[BUFFSIZE];
  int nptrs = backtrace(buffer, BUFFSIZE);

  //retrieve the stacktrace symbols and check for validity
  char **strings = backtrace_symbols(buffer, nptrs);
  if (strings == NULL) {
    perror("backtrace_symbols");
    exit(EXIT_FAILURE);
  }

  //print the lines of stacktrace one by one
  //skip the first two entries, since they are part of error.C
  for (int j = 2; j < nptrs; j++)
  {

#ifdef NDEBUG
    com->printf(0,"%s\n", strings[j]);
#else
    char syscom[256];
    sprintf(syscom,"addr2line %p -e ",buffer[j]);
    sprintf(syscom+ strlen(syscom),getenv("_"));//append current executables name
    com->system(syscom);
#endif

  }

  free(strings);
}

//------------------------------------------------------------------------------

//printing the stacktrace
void printStack(MPI_Comm comm)
{

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_rank(MPI_COMM_WORLD, &size);

  //get the stacktrace size
  void *buffer[BUFFSIZE];
  int nptrs = backtrace(buffer, BUFFSIZE);

  //retrieve the stacktrace symbols and check for validity
  char **strings = backtrace_symbols(buffer, nptrs);
  if (strings == NULL) {
    perror("backtrace_symbols");
    exit(EXIT_FAILURE);
  }

  //print the lines of stacktrace one by one
  //skip the first two entries, since they are part of error.C
  for (int j = 2; j < nptrs; j++)
  {

#ifdef NDEBUG
    if (rank==0) fprintf(stdout,"%s\n", strings[j] );
#else
    char syscom[256];
    sprintf(syscom,"addr2line %p -e ",buffer[j]);
    sprintf(syscom+ strlen(syscom),getenv("_"));//append current executables name
    if (rank == 0) std::system(syscom);
#endif

  }

  free(strings);
}

//------------------------------------------------------------------------------

/*****************************************************************************
 * This functions prints the distributed solution vector 'vec' to the        *
 * standard output. It is assumed that every node holds exacly one vector,   *
 * otherwise this function will not work yet.                                *
 *****************************************************************************/
template<int dim>
void printDistSVec(Communicator* com,DistSVec<double,dim> &vec)
{
  com->barrier();

  int rank=com->cpuNum();
  int size=com->size();

  int go=0;
  if(rank==0)
  {
    char message[50];
    sprintf(message,"This is processor number: %i\n",rank);
    std::cout<<"This Vector has a length of "<<vec(0).len<<std::endl<<std::flush;
    //vec(0).print(message);
    fflush(stdout);
    //sleep(1);

    go=1;
    com->sendTo(rank+1,rank,&go,1);
  }
  else
  {
    com->recFrom(rank-1, rank-1, &go, 1);

    char message[50];
    sprintf(message,"This is processor number: %i\n",rank);
    std::cout<<"This Vector has a length of "<<vec(0).len<<std::endl<<std::flush;//TODO delete line
    //vec(0).print(message);
    fflush(stdout);
    //sleep(1);

    if (rank<size-1)
      com->sendTo(rank+1,rank,&go,1);
  }
  com->barrier();
}

//------------------------------------------------------------------------------

template<int dim>
void writeDistSVec(Communicator* com, char* filename,DistSVec<double,dim> &vec)
{
  int rank=com->cpuNum();

  char fullfile[50];
  sprintf(fullfile,"%s_%i",filename,rank);

  std::ofstream myfile;
  myfile.open(fullfile);

  myfile.flush();
  for (int i=0; i<vec(0).len; ++i)
  {
    myfile << i << ": ";
    for (int j=0; j<dim; ++j)
      myfile << vec[i][j] << " ";
    myfile <<std::endl;
  }

  myfile.flush();

  myfile.close();
}

//------------------------------------------------------------------------------

template<class Scalar>
void writeDenseMatrix(Communicator* com, char* filename,GenFullM<Scalar> &mat)
{
  int rank=com->cpuNum();

  char fullfile[50];
  sprintf(fullfile,"%s_%i",filename,rank);

  std::ofstream myfile;
  myfile.open(fullfile);

  myfile.flush();


  int i,j;
  for(i = 0 ; i < mat.numRow() ; ++i) {
    for(j=0; j < mat.numCol() ; ++j)
      myfile<<mat[i][j]<<" ";
    myfile<<"\n";
  }

  myfile.flush();

  myfile.close();

}//end writeDenseMatrix(~)

//------------------------------------------------------------------------------

template<class Scalar>
void writeDenseMatrix(char* filename,GenFullM<Scalar> &mat)
{

  char fullfile[50];
  sprintf(fullfile,"%s_%0",filename);

  std::ofstream myfile;
  myfile.open(fullfile);

  myfile.flush();


  int i,j;
  for(i = 0 ; i < mat.numRow() ; ++i) {
    for(j=0; j < mat.numCol() ; ++j)
      myfile<<mat[i][j]<<" ";
    myfile<<"\n";
  }

  myfile.flush();

  myfile.close();

}//end writeDenseMatrix(~)

//------------------------------------------------------------------------------

template void printDistSVec(Communicator* com,DistSVec<double,5> &vec);
template void printDistSVec(Communicator* com,DistSVec<double,6> &vec);
template void printDistSVec(Communicator* com,DistSVec<double,7> &vec);

template void writeDistSVec(Communicator* com,char* filename,DistSVec<double,5> &vec);
template void writeDistSVec(Communicator* com,char* filename,DistSVec<double,6> &vec);
template void writeDistSVec(Communicator* com,char* filename,DistSVec<double,7> &vec);

template void writeDenseMatrix(Communicator* com, char* filename,GenFullM<double> &mat);
template void writeDenseMatrix(char* filename,GenFullM<double> &mat);
template void writeDenseMatrix(Communicator* com, char* filename,GenFullM<std::complex<double> > &mat);
template void writeDenseMatrix(char* filename,GenFullM<std::complex<double> > &mat);

}//end namespace Dev

