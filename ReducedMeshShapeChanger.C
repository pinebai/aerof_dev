#include <ReducedMeshShapeChanger.h>
#include <GappyPreprocessing.h>

template<int dim>
ReducedMeshShapeChanger<dim>::ReducedMeshShapeChanger(Communicator *_com, IoData &_ioData, Domain &dom, DistGeoState *_geoState) : 
GappyPreprocessing<dim>(_com, _ioData, dom, _geoState) {

	this->globalNodes = new std::vector <int> [1];
}

template<int dim>
ReducedMeshShapeChanger<dim>::~ReducedMeshShapeChanger() { 

	for (int i = 0; i < this->nReducedNodes; ++i) 
		delete [] xyz[i];
	delete [] xyz;

}

template<int dim>
void ReducedMeshShapeChanger<dim>::buildReducedModel() {

	// read in full nodes
//	readReducedNodes(this->input->reducedfullnodemap );

//	fillXYZ();	// compute xyz coordinates

//	if (this->thisCPU == 0) {
//		readWriteTopFile(this->input->mesh);
//	}
	// read in TOP file while writing out new top file
	// write out dwall file

//	this->outputWallDistanceReduced();
}

template<int dim>
void ReducedMeshShapeChanger<dim>::readReducedNodes(const char *reducedNodeFileName)  {

	// INPUT: reduced node file name
	// OUTPUT: this->nReducedNodes,this->globalNodes[0]

//	FILE *reducedNodeFile = fopen(reducedNodeFileName, "r");

//	int currentReducedNode;
//	this->nReducedNodes = 0;

//	while (fscanf(reducedNodeFile, "%d",&currentReducedNode) != EOF) {
//		(this->globalNodes)[0].push_back(currentReducedNode-1);	// reads in the reduced node plus one
//		++(this->nReducedNodes);
//	}

}

template<int dim>
void ReducedMeshShapeChanger<dim>::fillXYZ()  {
/*
	// initialize
	xyz = new double * [this->nReducedNodes];
	int * cpus = new int [this->nReducedNodes];
	int * subdomains = new int [this->nReducedNodes];
	int * localnodes = new int [this->nReducedNodes];

	for (int i = 0; i < this->nReducedNodes; ++i) {
		cpus[i] = 0;
		subdomains[i] = 0;
		localnodes[i] = 0;
		xyz[i] = new double [3];
		for (int j = 0; j < 3; ++j) {
			xyz[i][j] = 0.0;
		}
	}

	int globalNodeNum;

	for (int iSub = 0 ; iSub < this->numLocSub ; ++iSub) {	// all subdomains
		int *locToGlobNodeMap = this->subD[iSub]->getNodeMap();
		bool *locMasterFlag = this->nodeDistInfo.getMasterFlag(iSub); // master nodes on subdomain
		for (int iLocNode = 0; iLocNode < this->subD[iSub]->numNodes(); ++iLocNode) {	// all local this->globalNodes in subdomain
			if (locMasterFlag[iLocNode]) { 	// only do for master nodes
				globalNodeNum = locToGlobNodeMap[iLocNode];	// global node number
				for (int iReducedNode = 0; iReducedNode < this->nReducedNodes; ++iReducedNode) {
					if (globalNodeNum == this->globalNodes[0][iReducedNode]) {
						this->computeXYZ(iSub,iLocNode,xyz[iReducedNode]);
						cpus[iReducedNode] = this->thisCPU;
						subdomains[iReducedNode] = iSub;
						localnodes[iReducedNode] = iLocNode;
						break;
					}
				}
			}
		}
	}

	// make consistent across all processors
	for (int iReducedNode = 0; iReducedNode < this->nReducedNodes; ++iReducedNode) {
		this->com->globalSum(3, xyz[iReducedNode]);
		this->com->globalSum(1, cpus+iReducedNode);
		this->com->globalSum(1, subdomains+iReducedNode);
		this->com->globalSum(1, localnodes+iReducedNode);
		this->globalNodeToCpuMap.insert(pair<int, int > (this->globalNodes[0][iReducedNode], cpus[iReducedNode]));
		this->globalNodeToLocSubDomainsMap.insert(pair<int, int > (this->globalNodes[0][iReducedNode], subdomains[iReducedNode]));
		this->globalNodeToLocalNodesMap.insert(pair<int, int > (this->globalNodes[0][iReducedNode], localnodes[iReducedNode]));
	}
*/
}

template<int dim>
void ReducedMeshShapeChanger<dim>::readWriteTopFile(const char *inputTopFileName)  {
/*
	ifstream inputTopFile (inputTopFileName);	// old reduced mesh top file

	// new reduced mesh top file
	char *outputTopFileName = new char[strlen(this->ioData->output.transient.prefix) +
		strlen(this->ioData->output.rom.mesh)+1];
	strcpy(outputTopFileName,this->ioData->output.transient.prefix);
	strcat(outputTopFileName,this->ioData->output.rom.mesh);

	ofstream outputTopFile(outputTopFileName);
	outputTopFile.setf(ios::scientific);
	outputTopFile.precision(6);

	string line;
	getline(inputTopFile, line);
	outputTopFile << line; 

	for (int iReducedNode = 0; iReducedNode < this->nReducedNodes; ++iReducedNode) {
		getline(inputTopFile, line);	// read over node
		outputTopFile << endl << iReducedNode + 1 << " " << xyz[iReducedNode][0] << " "
			<< xyz[iReducedNode][1] << " " <<xyz[iReducedNode][2];
	}

	while (inputTopFile.good() ) {
		getline(inputTopFile, line);
		outputTopFile << endl << line;
	}
*/
}
