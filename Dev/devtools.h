/*
 * devtools.h
 *
 *  Created on: Nov 16, 2016
 *      Author: Lukas Scheucher
 *      Mail:   lscheuch@stanford.edu
 */

#ifndef DEVTOOLS_H_
#define DEVTOOLS_H_


#include "Communicator.h"
#include "../DistVector.h"
#include "../DenseMatrix.h"
#include <mpi.h>


/************************************************************************
 * This namepace contains some auxiliary functions that I found useful  *
 * for my implementation. Once mature, a lot of this stuff should       *
 * probably better be moved to the corresponding classes, for a         *
 * better code design. Keeping it here, decreased the compile times     *
 * during development.                                                  *
 ************************************************************************/
namespace Dev{
//print a warning to the terminal
//optional oputput of stacktrace
void Warning(Communicator *com,const char* message,bool printBacktrace=false);

void Warning(MPI_Comm comm, const char* message,bool printBacktrace=false);

//print an error to the terminal
//optional oputput of stacktrace
//exit program execution
void Error(Communicator *com,const char* message,bool printBacktrace=true);

void Error(MPI_Comm comm, const char* message,bool printBacktrace=true);

//print the stacktrace
void printStack(Communicator *com);

//print the stacktrace
void printStack(MPI_Comm comm);


template<int dim>
void printDistSVec(Communicator* com,DistSVec<double,dim> &vec);

template<int dim>
void writeDistSVec(Communicator* com, char* filename,DistSVec<double,dim> &vec);

template<class Scalar>
void writeDenseMatrix(Communicator* com, char* filename,GenFullM<Scalar> &mat);

template<class Scalar>
void writeDenseMatrix(char* filename,GenFullM<Scalar> &mat);

}


#endif /* ERROR_H_ */
