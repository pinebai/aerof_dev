/* KDTree.h

 */

#pragma once

#include <list>
#include <cstring>

#include <cfloat>

class KDTree {

 public:

  KDTree();

  ~KDTree();

  struct ScalarGlob {
    double x[3];
    int id;
  };

  struct KDTreeNode {

    int axis;
    ScalarGlob data;
    KDTreeNode* left,*right;
  };
  // Test comment
  template <int d>
  static int globComp(const void* a, const void* b) {
    
    return ( ((ScalarGlob*)a)->x[d] - ((ScalarGlob*)b)->x[d] ) > 0.0 ? 1 : -1;
  }

  void construct(SVec<double,3>& X) {

    Vec<ScalarGlob>* Xlocal = new Vec<ScalarGlob>(X.size());
    for (int i = 0; i < X.size(); ++i) {
      memcpy((*Xlocal)[i].x,X[i],sizeof((*Xlocal)[i].x));
      (*Xlocal)[i].id = i;
    }
    
    root = construct(*Xlocal, 0,0,X.size());
    delete Xlocal;
  }    

  void collectPointsInRadius(double x0[3],double radius,  
			     std::list<int>& ids) {

    //ids.clear();

    for (int k = 0; k < 3; ++k) {
      if (x0[k] + radius < boundingBox[0][k] ||
	  x0[k] - radius > boundingBox[1][k])
	return;
    }

    collectPointsInRadius(x0,radius, ids, root);
  }

 private:
     
  KDTreeNode* root;

  double boundingBox[2][3];

  void collectPointsInRadius(double x0[3],double radius,
			     std::list<int>& ids,KDTreeNode* node) {

    double dist = 0.0;
    for (int k = 0; k < 3; ++k)
      dist += pow(node->data.x[k] - x0[k],2);
    
    dist = sqrt(dist);
    if (dist <= radius)
      ids.push_back(node->data.id);
    
    if (node->left) {

      dist = node->data.x[node->axis] - x0[node->axis];
      if (dist >= 0.0) {
	collectPointsInRadius(x0,radius,ids,node->left);
	if (dist <= radius)
	  collectPointsInRadius(x0,radius,ids,node->right);
      } else {
	collectPointsInRadius(x0,radius,ids,node->right);
	if (-dist <= radius)
	  collectPointsInRadius(x0,radius,ids,node->left);
      }
    }     
    
  }

  KDTreeNode* construct(Vec<ScalarGlob>& Xlocal,int depth,int min,int max) {

    if (max-min < 1)
      return NULL;

    int axis;
    double bbox[2][3] = {{FLT_MAX,FLT_MAX,FLT_MAX},
			 {-FLT_MAX,-FLT_MAX,-FLT_MAX}};

    for (int i = min; i < max; ++i) {

      for (int k = 0; k < 3; ++k) {
	if (Xlocal[i].x[k] < bbox[0][k])
	  bbox[0][k] = Xlocal[i].x[k];
	
	if (Xlocal[i].x[k] > bbox[1][k])
	  bbox[1][k] = Xlocal[i].x[k];
      }
    }

    if (depth == 0)
      memcpy(boundingBox,bbox,sizeof(bbox));
    
    double size[3];
    for (int k = 0; k < 3; ++k)
      size[k] = bbox[1][k]-bbox[0][k];
    
    if (size[0] > size[1] &&
	size[0] > size[2])
      axis = 0;
    else if (size[1] > size[2])
      axis = 1;
    else
      axis = 2;
      
    if (axis == 0)
      Xlocal.sort(KDTree::globComp<0>, min,max);
    else if (axis == 1)
      Xlocal.sort(KDTree::globComp<1>, min,max);
    else // axis == 2
      Xlocal.sort(KDTree::globComp<2>, min,max);
    
    int median = (min+max)/2;
    
    KDTreeNode* node = new KDTreeNode;
    node->axis = axis;
    memcpy(&node->data,&(Xlocal[median]),sizeof(node->data));
    
    node->left = construct(Xlocal, depth+1,min,median);
    
    node->right = construct(Xlocal, depth+1,median+1,max);
    

    return node;
    
  }


};
