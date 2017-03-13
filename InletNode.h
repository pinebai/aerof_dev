#ifndef _INLET_NODE_H_
#define _INLET_NODE_H_

//libraries
#include <Node.h>
//#include <Elem.h>
#include <Vector.h>
#include <Vector3D.h>


//friends (class, structure, template class)

template<int dim> class Extrapolation;
template<int dim> class BcData;

class VarFcn;
class Connectivity;
class ElemSet;
class GeoState;
class IoData;

//------------------------------------------------------------------------------

class InletNode {
	
	int node;				//node number in the subdomain
	int numTets;			//# of tets connected to the node
	int numFaces;			//# of faces connected to the node
	int *tets;				//list of the tets connected to the node
        int *tets2;				//list of the secondary tetrahedra used for extrapolation
	int *faces;				//list of the faces connected to the node
        bool master;                    //tells if a shared node is the master one for extrapolation (this is
					//	different from the master flag for nodes and edges in subdomain)	
	
  public:
  	
  	InletNode();
  	~InletNode();
  	
  	int &operator[] (int i) const {return tets[i]; }
  	void setInletNode(int, int, int, int*, int*);
	void addSecondTets(int*);
        void chooseExtrapolation(int, double, int*, int*, double*, double*, double*);
  	int getNodeNum()   { return node; }
  	int getNumTets()   { return numTets; }
  	int getNumFaces()  { return numFaces; }
  	int *getTets()     { return tets; }
  	int *getTets2()     { return tets2; }
  	int *getFaces()    { return faces; }
        bool getMaster()   {return master; }
  	
  	void checkInletNodes(int, int*);
  	
	template<int dim>
	void computeZeroExtrapolation(VarFcn*, bool, Vec3D&, double*, double*, double*,
					 double*, int*, int*,  SVec<double,dim>&,
					 SVec<double,3>&, int, int* =0);
        template<int dim>
        void computeZeroExtrapolation(VarFcn*, bool, Vec3D&, double*, double*, int, double*,
                                         double*, int*, int*,  SVec<double,dim>&,
                                         SVec<double,3>&, int, int* =0);
	
	template<int dim>
	void assignFreeStreamValues(int, double *, double *, double *);

        template<int dim>
        void computeDifference(VarFcn*, SVec<double,dim>&, SVec<double,dim>&);


};

//------------------------------------------------------------------------------

class InletNodeSet {
	
	
        int locSubNum;
	int numInletNodes;
	InletNode *inletNodes;
	
  public:
  
  	InletNodeSet();
  	~InletNodeSet();
  	
  	InletNode &operator[] (int i) const { return inletNodes[i]; }
  	int size() const { return numInletNodes; }
  	void setup(int, int, IoData &);
  	void checkInletNodes(int*);
  	

        template<int dim>
  	void recomputeRHS(VarFcn*, Extrapolation<dim>*, ElemSet &, SVec<double,dim>&,
			  BcData<dim>& , GeoState&, SVec<double,dim>&, SVec<double,3>&, int *);
        template<int dim>
        void recomputeRHS(VarFcn*, Extrapolation<dim>*, ElemSet &, SVec<double,dim>&,
                          Vec<int> &, BcData<dim>& , GeoState&, SVec<double,dim>&,
                          SVec<double,3>&, int *);
        template<int dim>
        void recomputeResidual(SVec<double,dim> &F, SVec<double,dim> &Finlet);
  	
  	template<int dim>
  	void getExtrapolationValue(Extrapolation<dim>*,SVec<double,dim> &, SVec<double,dim> &,
  			   VarFcn* , BcData<dim>& , GeoState&, ElemSet &, int*, SVec<double,3>&);
  	
  	template<int dim>
        void applyExtrapolationToSolutionVector(Extrapolation<dim>*, SVec<double,dim> &, SVec<double,dim> &, int*);
  	
  	template<int dim>
	void printVariable(SVec<double,dim>&, Connectivity*, VarFcn*);	

        template<int dim>
        void printInletVariable(SVec<double,dim>&, Connectivity*);
                                                                                                                                                                                                     
        template<int dim>
        void checkExtrapolationValue(SVec<double,dim>&, Connectivity*, int*, VarFcn* , BcData<dim>& , GeoState&);
                                                                                                                                                                                                     

};

//------------------------------------------------------------------------------

#ifdef TEMPLATE_FIX
#include <InletNode.C>
#endif


#endif
