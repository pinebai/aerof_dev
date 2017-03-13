#ifndef _MULTIGRID_LEVEL_H_
#define _MULTIGRID_LEVEL_H_

#include <DistInfo.h>
#include <DistVector.h>
#include <MultiGridSmoothingMatrix.h>
#include <DistMvpMatrix.h>
#include <DistTimeState.h>
#include <SparseMatrix.h>
#include <MultigridCommon.h>
#include <AgglomeratedFace.h>
#include <Aerof_unordered_set.h>
#include <list>
#include <set>


class Domain;
class Connectivity;
class EdgeSet;
class EdgeDef;
template<class Scalar, int dim> class DistSVec;

typedef Aerof_unordered_set<int>::type PriorityNodes;

struct NeighborDomain {

  int id;
  // Vector pair of nodes (global, local)
  std::vector<std::pair<int,int> > sharedNodes;
};

struct MultigridSubdomain {

  enum Topology { TopoVertex = 0, TopoLine = 1, TopoFace = 2, TopoInterior = 3, TopoUnknown = 4 };
  
  std::vector<int> ownedNodes;
  std::vector<int> sharedNodes;

  int numNodes;

  std::map<int, NeighborDomain*> neighbors;

  int* localNodeMapping;

  int* ctrlVolCount;

  int* locToGlobMap;

  double* nodeVolumes;
  Vec3D* X;

  Topology* nodeTopology;
};

template<class Scalar>
class MultiGridLevel {

 public:

  enum AgglomerationType { AgglomerationLocal,
			   AgglomerationGlobal };

 private:

  AgglomerationType agglomType;


    int my_dim, neq1,neq2;

    double turbRelaxCutoff;

    Domain& domain;
    CommPattern<int> * nodeIdPattern;
    CommPattern<double> * nodeVolPattern;
    CommPattern<double> * nodeVecPattern;
    CommPattern<double> * nodeVecPatternEq1;
    CommPattern<double> * nodeVecPatternEq2;
    CommPattern<double> * nodeNormalsPattern;
    CommPattern<double> * nodePosnPattern;
    CommPattern<double> * matPattern;
    CommPattern<double> * matPattern2;
    CommPattern<double> * offDiagMatPattern;
    CommPattern<double> * edgeAreaPattern;
    CommPattern<double> * edgeVecPattern;

    Connectivity * subToSub;
    SubDTopo* levelSubDTopo;

    Connectivity ** sharedNodes;
    Connectivity** nodeToNodeMaskILU;
    int maxNodesPerSubD;
    mutable std::map<int,std::map<int,int> > locToGlobMap;

    MultiGridMethod mgMethod;

    DistGeoState* myGeoState;

    IoData* myIoData;
    
    std::map<int,std::set<int> >* toTransfer;

    std::list<Vec3D>** nodeNormals;

    bool useVolumeWeightedAverage;
    
    MultigridSubdomain* mgSubdomains;

    inline int rcvChannel(int glSub,int neighb) const {
      return levelSubDTopo->getChannelID(neighb,glSub);
    }

    inline int sndChannel(int glSub,int neighb) const {
      return levelSubDTopo->getChannelID(glSub,neighb);
    }

  protected:
    DistInfo * nodeDistInfo;
    DistInfo * edgeDistInfo;
    DistInfo* faceDistInfo;
    DistInfo* agglomFaceDistInfo;
    DistInfo* faceNormDistInfo;
    DistInfo* inletNodeDistInfo;

    int numLocSub;
    bool ownsData;
    Connectivity ** connectivity;
    EdgeSet ** edges;
    FaceSet** faces;
    AgglomeratedFaceSet** agglomeratedFaces;
    EdgeDef*** sharedEdges;
    int** numSharedEdges;

    DistVec<int>* pNodeMapping;
    DistVec<int>* pEdgeMapping;
    DistVec<int>* faceMapping;
    DistVec<Vec3D>* edgeNormals;
    DistSVec<double,3>* globalFaceNormals;

    DistSVec<int, 2>* lineMap;
    DistSVec<double, 3>* Xn;

    DistVec<int>* finestNodeMapping;

    DistVec<int> *lineIDMap;
    
    DistVec<int> *lineLocIDMap;

    DistVec<double>* ctrlVol; 
    DistVec<double>* ctrlVolCount;
    
    DistVec<int>* mgLevelStatus;

    int* numLines;

    int** lineLengths;

    int* numNodes;

    std::vector<int>* lineids;

    void identifyEdges(CommPattern<int> &edgeNumPat, int mySub);

    void sndEdgeInfo(CommPattern<int> &edgeNumPat, int mySub);

    void rcvEdgeInfo(CommPattern<int> &edgeNumPat, int mySub,int dim);

    MultiGridLevel* parent;
  
    DistSVec<int,2>* fv_comp_tag;

    void computeNodeNormalClasses();
 
    DistVec<int>* nodeNormalCount;
  
    double mesh_topology_threshold;

    double total_mesh_volume;

    DistVec<double>* edgeArea;

    template<class T,int dim>
      void assembleInternal(DistMat<T,dim> & A) const;

    template<class T,int dim>
      void assembleInternal(DistSVec<T,dim> & A) const;

    template<class T>
      void assembleInternal(DistVec<T> & A) const;

    template<class T>
      void assembleInternalMax(DistVec<T> & A) const;

  public:
    MultiGridLevel(MultiGridMethod,MultiGridLevel*,Domain& domain, DistInfo& refinedNodeDistInfo, DistInfo& refinedEdgeDistInfo);

    MultiGridLevel(MultiGridLevel* parent,int level_num,
                   double ref_length,
		   Domain& domain,const char* fn_base,
		   int dim, int neq1, int neq2,
		   GeoSource& geoSource);

    // level1 is finer than level2
    //MultiGridLevel(MultiGridLevel& level1, MultiGridLevel& level2);

    ~MultiGridLevel();

    MultigridSubdomain* getMgSubDomains() { return mgSubdomains; }

    MultiGridLevel* getParent() { return parent; }

    int getMyLevel() const {
      if (!parent) return 0;
      else return parent->getMyLevel()+1;
    }

    enum Topology { TopoVertex = 0, TopoLine = 1, TopoFace = 2, TopoInterior = 3, TopoUnknown = 4 };

    DistVec<Topology>* getFinestTopology() {

      if (!parent) return NULL;
      if (parent->parent) return parent->getFinestTopology();
      return nodeTopology;
    }

    int mapFineToCoarse(int iSub,int i) {
      if (parent)
        return (*pNodeMapping)(iSub)[parent->mapFineToCoarse(iSub,i)];
      else
        return i;
    }

    int mapFaceFineToCoarse(int iSub,int i) {
      if (parent)
        return (*faceMapping)(iSub)[parent->mapFaceFineToCoarse(iSub,i)];
      else
        return i;
    }

    void mergeFinerInto(MultiGridLevel& finer);
  
    double getTotalMeshVolume() { return parent ? parent->getTotalMeshVolume() : total_mesh_volume;  }

    MultiGridLevel* getFinestLevel() {

      if (parent)
        return parent->getFinestLevel();
      else
        return this;
    }   

    template <int dim>
    void setupBcs(DistBcData<dim>&, DistBcData<dim>&,DistSVec<Scalar,dim>&);

    DistGeoState& getGeoState() const { return *myGeoState; }

    DistSVec<double,3>& getXn() const { return myGeoState->getXn(); }
    
    void mapNodeList(int iSub,Aerof_unordered_set<int>::type&);

    enum SeedNodeChoice { Random, Mavripilis };

    void setUseVolumeWeightedAverage(bool b) { useVolumeWeightedAverage = b; }

    DistInfo& getNodeDistInfo()       { return *nodeDistInfo; }
    DistInfo& getEdgeDistInfo()       { return *edgeDistInfo; }
    DistInfo& getFaceDistInfo()       { return *faceDistInfo; }
    DistInfo& getAgglomFaceDistInfo()       { return *agglomFaceDistInfo; }
    DistInfo& getInletNodeDistInfo()       { return *inletNodeDistInfo; }
    Connectivity ** getConnectivity() { return connectivity; }
    EdgeSet ** getEdges()             { return edges; }
    FaceSet ** getFaces()             { return faces; }
    AgglomeratedFaceSet ** getAgglomeratedFaces()             { return agglomeratedFaces; }
    CommPattern<int>& getIdPat()      { return *nodeIdPattern; }
    Connectivity ** getSharedNodes()  { return sharedNodes; }

    DistVec<int>& getNodeMapping() { return *pNodeMapping; }

    DistVec<Vec3D>& getEdgeNormals() { return *edgeNormals; }
    DistVec<double>& getCtrlVol() const { return *ctrlVol; }

    EdgeDef*** getSharedEdges() { return sharedEdges; }
    int** getNumSharedEdges() { return numSharedEdges; }

    int* getNodeType(int iSub) { return nodeType[iSub]; }
 
    DistSVec<int,2>& getFVCompTag() const { return *fv_comp_tag; }

    int getNodeWithMostAgglomeratedNeighbors(std::vector<std::set<int> >& C,//Connectivity* C, 
                                             PriorityNodes& P,
                                             Vec<int>& nodeMapping,int iSub);
    void copyRefinedState(const DistInfo& refinedNodeDistInfo, const DistInfo& refinedEdgeDistInfo, const DistInfo& refinedFaceDistInfo,const DistInfo& refinedInletNodeDistInfo,DistGeoState& refinedGeoState, Domain& domain, int dim, int lneq1,int lneq2);

    int getNumNeighborSubdomains(int mySub);

    int* getNeighboringSubdomains(int mySub);

    void agglomerate(const DistInfo& refinedNodeDistInfo,
                     const DistInfo& refinedEdgeDistInfo,
                     CommPattern<int>& refinedNodeIdPattern,
                     Connectivity** refinedSharedNodes,
                     Connectivity ** nToN, EdgeSet ** edges,
                     EdgeDef***, int**,
                     Domain& domain,int dim,int neq1,int neq2,
                     DistVec<Vec3D>& refinedEdgeNormals,
                     DistVec<double>* refinedEdgeAreas,
                     DistVec<double>& refinedVol,
                     DistVec<int>*,double beta,int maxNumNodesPerAgglom);

    void computeRestrictedQuantities(const DistGeoState& refinedGeoState);

    template<class Scalar2, int dim> void Restrict(const MultiGridLevel<Scalar>& fineGrid,
                                                   const DistSVec<Scalar2, dim>& fineData,
                                                   DistSVec<Scalar2, dim>& coarseData,
						   bool average = true,
						   bool apply_relaxation = false) const;

    template<class Scalar2, int dim> void RestrictFaceVector(const MultiGridLevel<Scalar>& fineGrid,
                                                   const DistSVec<Scalar2, dim>& fineData,
                                                   DistSVec<Scalar2, dim>& coarseData) const;
    template<class Scalar2>void Restrict(const MultiGridLevel<Scalar>& fineGrid,
                                         const DistVec<Scalar2>& fineData,
                                         DistVec<Scalar2>& coarseData) const;
    template<class Scalar2, int dim> void Prolong(MultiGridLevel<Scalar>& coarseGrid, const DistSVec<Scalar2,dim>& coarseInitialData,
                                                  const DistSVec<Scalar2,dim>& coarseData, DistSVec<Scalar2,dim>& fineData,DistSVec<Scalar2,dim>& fine_ref,double relax_factor, VarFcn*, class DistLevelSetStructure* coarse = NULL,class DistLevelSetStructure* fine = NULL) const;

    template<class Scalar2, int dim>
      void ExtrapolateProlongation(MultiGridLevel<Scalar>& fineGrid, 
				   const DistSVec<Scalar2,dim>& coarseInitialData,
				   DistSVec<Scalar2,dim>& coarseData,
				   DistLevelSetStructure* coarseLSS, 
				   DistLevelSetStructure* fineLSS);

    template<class Scalar2, int dim>
    void ProjectResidual(DistSVec<Scalar2,dim>& r) const;

    template<class Scalar2, int dim>
    void RestrictOperator(const MultiGridLevel<Scalar>& fineGrid,
                          DistMat<Scalar2,dim>& fineOperator,
                          DistMat<Scalar2,dim>& coarseOperator);

    bool isLine(int iSub,int edgei,int edgej, int* lineid, int* loci, int* locj);

    int NumLines(int iSub) const { return numLines[iSub]; }

    int* getLineData(int iSub, int lineid) { return &lineids[iSub][lineid*8]; } 
    int lineLength(int iSub, int lineid) { return lineLengths[iSub][lineid]; } 

    template <class Scalar2>
    void assemble(DistVec<Scalar2>& V);
    template <class Scalar2,int dim>
    void assemble(DistSVec<Scalar2,dim>& V);
    template <class Scalar2,int dim>
    void assemble(DistMat<Scalar2,dim>& V);
    template <class Scalar2,int dim>
    void assembleMax(DistSVec<Scalar2,dim>& V);
    template <class Scalar2>
    void assembleMax(DistVec<Scalar2>& V);

    template <class Scalar2,int dim>
    void computeMatVecProd(DistMat<Scalar2,dim>& mat,
                           DistSVec<Scalar2,dim>& p,
                           DistSVec<Scalar2,dim>& prod);

    template <class Scalar2,int dim>
    void computeGreenGaussGradient(DistSVec<Scalar2,dim>& V,
                                   DistSVec<Scalar2,dim>& dX,
                                   DistSVec<Scalar2,dim>& dY,
                                   DistSVec<Scalar2,dim>& dZ, 
                                   DistLevelSetStructure* = NULL);
    
    void WriteTopFile(const std::string& fileName);

    template <int dim>
    void writeXpostFile(const std::string& fileName,
                        DistSVec<Scalar,dim>& val,int id);
    
    void writeXpostFile(const std::string& fileName,
                        DistVec<Scalar>& val);

    Connectivity* createEdgeBasedConnectivity(int iSub);
    compStruct* createRenumbering(int iSub,Connectivity *nodeToNode,
 		 		  int typeRenum, int print);

    template <int dim>
    SparseMat<Scalar,dim>* createMaskILU(int iSub,int fill, 
                                         int renum, int *ndType);

    void writePVTUFile(const char* filename);
    void writePVTUAgglomerationFile(const char* filename);

    template <int dim>
    void writePVTUSolutionFile(const char* filename, DistSVec<double,dim>&);

    void setNodeType();
 
    DistVec<double>& getEdgeArea() { return *edgeArea; }
    DistVec<double>* getEdgeAreaPointer() { return edgeArea; }

    void setTurbRelaxCutoff(double t) { turbRelaxCutoff = t; }

  private:

    DistVec<Topology>* nodeTopology;
    DistVec<Topology>* coarseNodeTopology;
    DistVec<Vec3D>* topologyNormal;
    DistVec<Vec3D>* coarseTopologyNormal;

    int** nodeType;

};

#endif
