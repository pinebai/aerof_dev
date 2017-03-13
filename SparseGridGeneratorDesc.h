#ifndef _SPARSEGRIDGENERATOR_DESC_H_
#define _SPARSEGRIDGENERATOR_DESC_H_

#include <IoData.h>
#include <SparseGrid.h>
#include <SparseGridCluster.h>
#include <cmath>
#include <ctime>
#include <Communicator.h>

class RefVal;
class VarFcn;
class LocalRiemann;

//------------------------------------------------------------------------------
  void functionTest(double *in, double *res, double *parameters){ 
    res[0] = in[0]*in[0]*in[1]*in[1];
  }

//------------------------------------------------------------------------------
class SparseGridGeneratorDesc {

private:

  VarFcn *varFcn;
  LocalRiemannGfmparGasJWL *lriemannGasJwl;
  Communicator *com;

  void memberFunctionTest(double *in, double *res, double *parameters){ 
    res[0] = in[0]*in[0]*in[1]*in[1];
  }

public:

  SparseGridGeneratorDesc(IoData &ioData, Communicator *comm){
    com = comm;
    varFcn = createVarFcn(ioData);
    lriemannGasJwl = new LocalRiemannGfmparGasJWL(varFcn, 0, 1, NULL, MultiFluidData::RK2,1.0,1.0,1.0);
  }

  ~SparseGridGeneratorDesc(){
    delete varFcn;
    delete lriemannGasJwl;
    com = 0;
  }

//------------------------------------------------------------------------------

  void tabulate(IoData & ioData){
    bool debugInfo = false;
    srand(time(NULL));

    if(debugInfo)
      com->fprintf(stdout, "### SparseGridGeneratorDesc::tabulate -- started\n");

    if(ioData.mf.riemannComputation == MultiFluidData::TABULATION2){

      double *refIn = new double[2]; double *refOut = new double[1];
      refIn[0] = ioData.ref.rv.density;
      refIn[1] = ioData.ref.rv.entropy;
      refOut[0] = ioData.ref.rv.velocity;
      //com->fprintf(stdout, "refIn/refOut are %e %e / %e\n", refIn[0], refIn[1], refOut[0]);

      double *parameters = new double[2];
      parameters[0] = 1.0;
      parameters[1] = ioData.eqs.fluidModel.jwlModel.rhoref;

      if(false){ // for testing SparseGrid structure with Riemann Invariants of JWL EOS
        SparseGrid sparseGrid(ioData.mf.sparseGrid, parameters, refIn, refOut);
        com->fprintf(stdout, "### SparseGridGeneratorDesc::tabulate -- 2\n");
        //sparseGrid.tabulate(functionTest);
        //sparseGrid.tabulate(&SparseGridGeneratorDesc::memberFunctionTest,*this);
        sparseGrid.tabulate(&LocalRiemannGfmparGasJWL::riemannInvariantGeneral2ndOrder_wrapper,*lriemannGasJwl);
        sparseGrid.printToFile(refIn, refOut, ioData.output.transient.sparseGrid);
        com->fprintf(stdout, "### SparseGridGeneratorDesc::tabulate -- 3\n");
        SparseGrid sparseGridCopy;
        sparseGridCopy.readFromFile(refIn, refOut, ioData.output.transient.sparseGrid);
        int number = 5;
        sparseGridCopy.test(&LocalRiemannGfmparGasJWL::riemannInvariantGeneral2ndOrder_wrapper,*lriemannGasJwl, 2, &number, parameters);
      }

      // actual SparseGrids tabulation
      SparseGridCluster sgCluster;
      int sp = strlen(ioData.output.transient.prefix) + 1;
      char *sparseGridOutFileName = new char[sp + strlen(ioData.output.transient.sparseGrid)];
      sprintf(sparseGridOutFileName, "%s%s", 
	    ioData.output.transient.prefix, ioData.output.transient.sparseGrid);
      sgCluster.generate(ioData.mf.sparseGrid, parameters, &LocalRiemannGfmparGasJWL::riemannInvariantGeneral2ndOrder_wrapper,*lriemannGasJwl, sparseGridOutFileName, refIn, refOut, com);

      // test the first tabulation of the cluster
      com->barrier();
      SparseGrid sparseGridCopy;
      char *sparseGridReadFileName = new char[sp + strlen(ioData.output.transient.sparseGrid)+1];
      sprintf(sparseGridReadFileName, "%s%s%d", 
	    ioData.output.transient.prefix, ioData.output.transient.sparseGrid, com->cpuNum()+1);
      sparseGridCopy.readFromFile(refIn, refOut, sparseGridReadFileName);
      int number = 5;
      sparseGridCopy.test(&LocalRiemannGfmparGasJWL::riemannInvariantGeneral2ndOrder_wrapper,*lriemannGasJwl, 1, &number, parameters);

      delete [] refIn; delete [] refOut; delete [] parameters;

    }else if(ioData.mf.riemannComputation == MultiFluidData::TABULATION5){
      double *parameters = NULL;
      double *refIn = new double[5]; double *refOut = new double[2]; // 2outputs
      refIn[0] = ioData.ref.rv.density;
      refIn[1] = ioData.ref.rv.pressure;
      refIn[2] = ioData.ref.rv.density;
      refIn[3] = ioData.ref.rv.pressure;
      refIn[4] = ioData.ref.rv.velocity;
      refOut[0] = ioData.ref.rv.density;
      refOut[1] = ioData.ref.rv.density; // 2outputs
      //com->fprintf(stdout, "refIn are %e %e %e\n", refIn[0],refIn[4],refIn[1]);

/*
      SparseGrid sparseGrid(ioData.mf.sparseGrid, parameters, refIn, refOut);
      if(true){
        sparseGrid.tabulate(&LocalRiemannGfmparGasJWL::eriemanngj_wrapper,*lriemannGasJwl);
        sparseGrid.printToFile(refIn, refOut, ioData.output.transient.sparseGrid);
      }

      SparseGrid sparseGridCopy;
      sparseGridCopy.readFromFile(refIn, refOut, ioData.output.transient.sparseGrid);
      int numTest = 5;
      if(true) 
        sparseGridCopy.test(&LocalRiemannGfmparGasJWL::eriemanngj_wrapper,*lriemannGasJwl,1,&numTest, parameters);

      if(true)
        sparseGridCopy.test(&LocalRiemannGfmparGasJWL::eriemanngj_wrapper,*lriemannGasJwl,2,&numTest, parameters);
*/

      SparseGridCluster sgCluster;
      int sp = strlen(ioData.output.transient.prefix) + 1;
      char *sparseGridOutFileName = new char[sp + strlen(ioData.output.transient.sparseGrid)];
      sprintf(sparseGridOutFileName, "%s%s", 
	    ioData.output.transient.prefix, ioData.output.transient.sparseGrid);
      sgCluster.generate(ioData.mf.sparseGrid, parameters, &LocalRiemannGfmparGasJWL::eriemanngj_wrapper,*lriemannGasJwl, sparseGridOutFileName, refIn, refOut, com);
      delete [] refIn; delete [] refOut; delete parameters;
    }else{
      com->fprintf(stdout, "### SparseGridGeneratorDesc::nothing done!\n");
    }
    if(debugInfo)
      com->fprintf(stdout, "### SparseGridGeneratorDesc::tabulate -- finished\n");
  }

//------------------------------------------------------------------------------

  VarFcn *createVarFcn(IoData &ioData){
    VarFcn *vf = 0;
    if(ioData.mf.riemannComputation == MultiFluidData::TABULATION2){
      //vf = new VarFcnJWLEuler3D(ioData);
      vf = new VarFcn(ioData);
      if(vf->getType(0) != VarFcnBase::JWL){
        fprintf(stdout, "*** Error: a JWL EOS is needed for this tabulation\n"); 
        exit(1);
      }
    }else if(ioData.mf.riemannComputation == MultiFluidData::TABULATION5){
      //vf = new VarFcnJWLInGasEuler3D(ioData);
      vf = new VarFcn(ioData);
      if((vf->getType(0) != VarFcnBase::PERFECTGAS && vf->getType(0) != VarFcnBase::STIFFENEDGAS) || vf->getType(1) != VarFcnBase::JWL){
        fprintf(stdout, "*** Error: a SG EOS AND  a JWL EOS are needed for this tabulation\n"); 
        exit(1);
      }
    }

    if(!vf){
      com->fprintf(stdout, "*** Error: no valid choice for the VarFcn\n");
      exit(1);
    }
    return vf;
  }
    

};

//------------------------------------------------------------------------------


#endif
