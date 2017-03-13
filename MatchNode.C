#include <MatchNode.h>

#ifdef OLD_STL
#include <defalloc.h>
#include <algo.h>
#else
#include <algorithm>
using std::stable_sort;
#endif

//------------------------------------------------------------------------------

struct MatchNode {

  int index[3];

  double gap[3];

  bool operator<(const MatchNode &node) const { return index[0] < node.index[0]; }

};

//------------------------------------------------------------------------------

template<class NodeMap>
void MatchNodeSet::renumberNodes(NodeMap &nodemap)
{

  // remap into local subdomain numbering

  int i, j;
  for (i=0; i<numNodes; ++i)
    index[i][0] = nodemap[ index[i][0] ];

  // reorder the matched nodes according to their local subdomain numbering

  MatchNode *orderNodes = new MatchNode[numNodes];
  
  for (i=0; i<numNodes; ++i) {
    for (j=0; j<3; ++j) {
      orderNodes[i].index[j] = index[i][j];
      orderNodes[i].gap[j] = gap[i][j];
    }
  }

#ifdef OLD_STL
  sort(orderNodes, orderNodes + numNodes);
#else
  stable_sort(orderNodes, orderNodes + numNodes);
#endif

  for (i=0; i<numNodes; ++i) {
    for (j=0; j<3; ++j) {
      index[i][j] = orderNodes[i].index[j];
      gap[i][j] = orderNodes[i].gap[j];
    }
  }

  if (orderNodes) delete [] orderNodes;

}

//------------------------------------------------------------------------------

template<int dim>
void MatchNodeSet::send(double scale, double (*f)[dim], double (*buffer)[dim])
{

  for (int i=0; i<numNodes; ++i)
    for (int k=0; k<dim; ++k)
      buffer[ index[i][2] ][k] = scale * f[ index[i][0] ][k];

}

//------------------------------------------------------------------------------

template<int dim>
void MatchNodeSet::sendWithMasterFlag(double scale, double (*f)[dim], double (*buffer)[dim], bool *masterFlag, int locOffset)
{

  for (int i=0; i<numNodes; ++i)
    for (int k=0; k<dim; ++k) {
      if(masterFlag[ locOffset + index[i][0] ])
        buffer[ index[i][2] ][k] = scale * f[ index[i][0] ][k];
      else
        buffer[ index[i][2] ][k] = 0.0;
    }

}
