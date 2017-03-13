#ifndef _ERROR_HANDLER_H_
#define _ERROR_HANDLER_H_

#include <Communicator.h>
#include <vector>

struct ErrorHandler{

  enum Error {UNPHYSICAL = 0, SATURATED_LS = 1, BAD_RIEMANN = 2, REDUCE_TIMESTEP = 3, PRESSURE_CLIPPING = 4, DENSITY_CLIPPING = 5, REDO_TIMESTEP = 6, LARGE_VELOCITY = 7, RAPIDLY_CHANGING_PRESSURE = 8, RAPIDLY_CHANGING_DENSITY = 9, REDUCE_TIMESTEP_TIME = 10, SIZE = 11};
  enum Type {ALL = 0, SOLVER = 1};
  int localErrors[SIZE];
  int globalErrors[SIZE];
  std::vector<int> *solverErrors;
  Communicator *com;

  ErrorHandler(Communicator *comIn){
    com = comIn;
    int solverErrorsArray[] = {UNPHYSICAL,SATURATED_LS,BAD_RIEMANN,PRESSURE_CLIPPING,DENSITY_CLIPPING,LARGE_VELOCITY,RAPIDLY_CHANGING_PRESSURE,RAPIDLY_CHANGING_DENSITY};
    solverErrors = new std::vector<int>(solverErrorsArray,solverErrorsArray + sizeof(solverErrorsArray)/sizeof(int));
    for (int i=0; i<SIZE; i++) localErrors[i]=0;
    for (int i=0; i<SIZE; i++) globalErrors[i]=0;
  }

  void reduceError(){ 
    int i, toto[solverErrors->size()];
    i = 0;
    for(std::vector<int>::iterator it = solverErrors->begin(); it != solverErrors->end(); ++it) toto[i++] = localErrors[*it];
    com->globalSum(solverErrors->size(),toto);
    i = 0;
    for(std::vector<int>::iterator it = solverErrors->begin(); it != solverErrors->end(); ++it) globalErrors[*it] = toto[i++];
    return;
  }

  void clearError(int type=ALL){
    if(type==ALL) for(int i=0; i<SIZE; i++) {globalErrors[i]=0; localErrors[i]=0; }
    if(type==SOLVER) for(int i=0; i<solverErrors->size(); i++) {globalErrors[solverErrors->at(i)]=0; localErrors[solverErrors->at(i)]=0;}
    return;
  }

  void printError(int type=ALL){
    char str[200];
    //sprintf(str,"");
    if(type==ALL) for(int i=0; i<SIZE; i++) {sprintf(str,"%s%i, ",str,globalErrors[i]);}
    if(type==SOLVER) for(int i=0; i<solverErrors->size(); i++) {sprintf(str,"%s%i, ",str,globalErrors[solverErrors->at(i)]);}
    std::printf("%s\n",str);
    fflush(stdout);
    com->barrier();
    return;
  }

};


#endif
