#ifndef _POLYGONRECONSTRUCTIONDATA_H_
#define _POLYGONRECONSTRUCTIONDATA_H_

class Elem;
class LevelSetStructure;
template<class Scalar, int dim> class SVec;

struct Vec3D;

//------------------------------------------------------------------------------

struct PolygonReconstructionData { //for force computation under the embedded framework
    PolygonReconstructionData() : numberOfEdges(0) {}
    int numberOfEdges,nodeToLookFrom;
    int edge[4];
    int edgeWithVertex[4][2];

    void AssignTriangleSingle(const int n1, const int n2, const int n3, const int n4, const int l1, const int l2, const int l3){
        numberOfEdges=3; nodeToLookFrom=n1;
        edge[0]=l1; edge[1]=l2; edge[2]=l3;
        edgeWithVertex[0][0]=n1; edgeWithVertex[0][1]=n2;
        edgeWithVertex[1][0]=n1; edgeWithVertex[1][1]=n3;
        edgeWithVertex[2][0]=n1; edgeWithVertex[2][1]=n4;
    }

    void AssignTriangleMulti(const int n1, const int n2, const int n3, const int n4, const int l1, const int l2, const int l3){
        numberOfEdges=3; nodeToLookFrom=n2;
        edge[0]=l1; edge[1]=l2; edge[2]=l3;
        edgeWithVertex[0][0]=n2; edgeWithVertex[0][1]=n1;
        edgeWithVertex[1][0]=n3; edgeWithVertex[1][1]=n1;
        edgeWithVertex[2][0]=n4; edgeWithVertex[2][1]=n1;
    }

    void AssignQuadrilateral(const int n1, const int n2, const int n3, const int n4, const int l1, const int l2, const int l3, const int l4){
        numberOfEdges=4; nodeToLookFrom=n1;
        edge[0]=l1; edge[1]=l2; edge[2]=l3; edge[3]=l4;
        edgeWithVertex[0][0]=n1; edgeWithVertex[0][1]=n3;
        edgeWithVertex[1][0]=n1; edgeWithVertex[1][1]=n4;
        edgeWithVertex[2][0]=n2; edgeWithVertex[2][1]=n4;
        edgeWithVertex[3][0]=n2; edgeWithVertex[3][1]=n3;
    }

    void AssignQuadTriangle(const int n1, const int n2, const int n3, const int n4, const int l1, const int l2, const int l3){ //for PhysBAM only
        numberOfEdges=3; nodeToLookFrom=n1;
        edge[0]=l1; edge[1]=l2; edge[2]=l3;
        edgeWithVertex[0][0]=n1; edgeWithVertex[0][1]=n3;
        edgeWithVertex[1][0]=n1; edgeWithVertex[1][1]=n4;
        edgeWithVertex[2][0]=n2; edgeWithVertex[2][1]=n3;
    }

    void AssignConnectingQuadTriangle(const int n1, const int n2, const int n3, const int n4, const int l1, const int l2){ //for PhysBAM only
        numberOfEdges=3; nodeToLookFrom=n3;
        edge[0]=l1; edge[1]=l1; edge[2]=l2;
        edgeWithVertex[0][0]=n1; edgeWithVertex[0][1]=n2;
        edgeWithVertex[1][0]=n2; edgeWithVertex[1][1]=n1;
        edgeWithVertex[2][0]=n3; edgeWithVertex[2][1]=n4;
    }

    void AssignTwoEdges(const int n1, const int n2, const int n3, const int n4, const int l1, const int l2){ //for PhysBAM only
        numberOfEdges=4; nodeToLookFrom=n4;
        edge[0]=l1; edge[1]=l2; edge[2]=l2; edge[3]=l1;
        edgeWithVertex[0][0]=n1; edgeWithVertex[0][1]=n2;
        edgeWithVertex[1][0]=n1; edgeWithVertex[1][1]=n3;
        edgeWithVertex[2][0]=n3; edgeWithVertex[2][1]=n1;
        edgeWithVertex[3][0]=n2; edgeWithVertex[3][1]=n1;
    }

    void AssignTwoEdgesPartTwo(const int n1, const int n2, const int n3, const int n4, const int l1, const int l2, const int l3){ //for PhysBAM only
        numberOfEdges=4; nodeToLookFrom=n4;
        edge[0]=l1; edge[1]=l3; edge[2]=l3; edge[3]=l2;
        edgeWithVertex[0][0]=n2; edgeWithVertex[0][1]=n1;
        edgeWithVertex[1][0]=n2; edgeWithVertex[1][1]=n3;
        edgeWithVertex[2][0]=n3; edgeWithVertex[2][1]=n2;
        edgeWithVertex[3][0]=n3; edgeWithVertex[3][1]=n1;
    }
};

//------------------------------------------------------------------------------

int getPolygons(Elem &elem, LevelSetStructure &LSS, PolygonReconstructionData* polygons);
void getPolygonNormal(SVec<double,3>& X, Vec3D &normal, LevelSetStructure &LSS, PolygonReconstructionData &polygon);

#endif
