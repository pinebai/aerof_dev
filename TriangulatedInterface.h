#pragma once

#include <IoData.h>

#include <Domain.h>

class TriangulatedInterface {

 public:
 
  TriangulatedInterface();

  ~TriangulatedInterface();

  int NumNodes();

  int NumElems();
  
  int (*getIndices())[3];

  double* getVertexLocations();

  void initializeAsSquare(int);

  void setExactSquare(double x, double v);

  void initializeIntersector(IoData& iod,  Communicator*, Domain*,
                              DistSVec<double,3>& X);


  void update(double dt);

  class DistLevelSetStructure* getIntersector();

  class LevelSetStructure* getSubLSS(int iSub);

  template <int dim>
  void integraterk2(Domain* domain, DistSVec<double,3>& X,
                    DistSVec<double,dim>& V,
                    DistNodalGrad<dim, double>& grad,
                    DistVec<int>& fid,
                    double dt) {

    int rnk;
    MPI_Comm_rank(MPI_COMM_WORLD,&rnk);
    double* xyzn = new double[numNodes*3];
    memcpy(xyzn, xyz,sizeof(double)*3*numNodes);
    step(domain, X, V, grad,fid);

    //if (rnk == 0)
    //  std::cout << "dt = " << dt << std::endl;   
 
    for (int i = 0; i < numNodes; ++i) {
      for (int k = 0; k < 3; ++k)
        xyz[i*3+k] += xyzdot[i*3+k]*dt*0.5;
    }
    step(domain, X, V, grad,fid);

    for (int i = 0; i < numNodes; ++i) {
      for (int k = 0; k < 3; ++k) {
        xyz[i*3+k] = xyzn[i*3+k] + dt*xyzdot[i*3+k];
      }
      //if (rnk == 0)
      //  std::cout << i << " " << xyzdot[i*3] << " " << xyzdot[i*3+1] << " " << xyzdot[i*3+2] <<
      //          " " << xyzn[i*3] << " " << xyzn[i*3+1] << " " << xyzn[i*3+2] << std::endl;
    }

    delete [] xyzn;    
    
  }

  template <int dim>
  void step(Domain* domain, DistSVec<double,3>& X,
            DistSVec<double,dim>& V,
            DistNodalGrad<dim, double>& grad,
            DistVec<int>& fid) {
 
    int* cnt1 = new int[numNodes];
    int* cnt2 = new int[numNodes];
    memset(cnt1, 0, sizeof(int)*numNodes);
    memset(cnt2, 0, sizeof(int)*numNodes);

    double* xyzdot1 = new double[3*numNodes];
    double* xyzdot2 = new double[3*numNodes];

  //  std::cout << "numNodes = " << numNodes << std::endl;
    memset(xyzdot1, 0, sizeof(double)*3*numNodes);
    memset(xyzdot2, 0, sizeof(double)*3*numNodes);

    for (int iSub = 0; iSub < V.numLocSub(); ++iSub) {

      SubDomain* S = domain->getSubDomain()[iSub];
      SVec<double,3>& subX = X(iSub);
      for (int i = 0; i < numNodes; ++i) {

        Vec3D p(xyz[i*3], xyz[i*3+1], xyz[i*3+2]);
        Elem* E = S->searchPoint(p, X(iSub));
        SVec<double,dim>& dx = grad.getX()(iSub);
        SVec<double,dim>& dy = grad.getY()(iSub);
        SVec<double,dim>& dz = grad.getZ()(iSub);

        if (!E)
          continue;

        int nn = E->numNodes();
        for (int m = 0; m < nn; ++m) {
    
          int nd = E->nodeNum()[m];

          double v[3];
          v[0] = dx[nd][1]*(p[0]-subX[nd][0])+
                 dy[nd][1]*(p[1]-subX[nd][1])+
                 dz[nd][1]*(p[2]-subX[nd][2]);

          v[1] = dx[nd][2]*(p[0]-subX[nd][0])+
                 dy[nd][2]*(p[1]-subX[nd][1])+
                 dz[nd][2]*(p[2]-subX[nd][2]);
          
          v[2] = dx[nd][3]*(p[0]-subX[nd][0])+
                 dy[nd][3]*(p[1]-subX[nd][1])+
                 dz[nd][3]*(p[2]-subX[nd][2]);
      
          v[0] += V(iSub)[nd][1];
          v[1] += V(iSub)[nd][2];
          v[2] += V(iSub)[nd][3];
          //std::cout << p[0] << " " << p[1] << " " << p[2] << " " << v[0] << std::endl;
          //std::cout << V(iSub)[nd][1] <<  " " << dx[nd][1] << " " << subX[nd][0] << std::endl;
          if (fid(iSub)[nd] == 0) {
            for (int k = 0; k < 3; ++ k) {
              xyzdot1[i*3+k] += v[k];
            }
            ++cnt1[i];
          } else {
            for (int k = 0; k < 3; ++ k) {
              xyzdot2[i*3+k] += v[k];
            }
            ++cnt2[i];

          }
        }
      }

    }

    Communicator* com = domain->getCommunicator();
    com->globalSum(numNodes*3,xyzdot1);
    com->globalSum(numNodes , cnt1);
    com->globalSum(numNodes*3,xyzdot2);
    com->globalSum(numNodes , cnt2);

    for (int i = 0; i < numNodes; ++i) {

      for (int k = 0; k < 3; ++k)
        xyzdot1[i*3+k] /= std::max(1,cnt1[i]);
      for (int k = 0; k < 3; ++k)
        xyzdot2[i*3+k] /= std::max(1,cnt2[i]);


      for (int k = 0; k < 3; ++k)
        xyzdot[i*3+k] = 0.5*(xyzdot1[i*3+k]+xyzdot2[i*3+k]);

//      for (int k = 0; k < 3; ++k)
//        xyzdot[i*3+k] = (xyzdot1[i*3+k]+xyzdot2[i*3+k])/(cnt1[i]+cnt2[i]);


      if (constrained[i] & 1) {
        xyzdot[i*3+1] = 0.0;
      }
      if (constrained[i] & 2) {
        xyzdot[i*3+2] = 0.0;
      }
        
    }

    delete [] cnt1;
    delete [] cnt2;    
    delete [] xyzdot1;
    delete [] xyzdot2;
    
  }

 
 private:

  int (*abc)[3];

  int *constrained;

  double* xyz;

  double* xyzdot;

  int numNodes;

  int numElems;

  class DistIntersectorPhysBAM* pIntersector;
};
