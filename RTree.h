/* RTree.h
 *
 */

#pragma once

// The first parameter is the object
template <class T>
class RTree 
{

	struct Node 
	{
    Node* child1,*child2;
    double bbox[6];
    T* obj;
  };

 public:

  RTree() { root = NULL; }

  ~RTree() { }
  
  Node* root;

  template <void (T::*bbox)(SVec<double,3>& X, double*)>
  void construct(SVec<double,3>& X, T** objects, int num) 
  {
	  root = constructInternal<bbox>(X, objects, num);
  }

  /* ---------- */

  void destruct() 
  {
    if (root) destructInternal(root);
  }
  
  /* ---------- */

  template <void (T::*bbox)(SVec<double,3>& X, double*)>
	  void reconstruct(SVec<double,3>& X,T** objects, int num) 
  {
	  if(root) destructInternal(root);

    root = constructInternal<bbox>(X,objects,num);
  }

  /* ---------- */

  template <bool (T::*inside)(SVec<double,3>&, const Vec3D&),
            class T2, bool (T2::*validity)(T*)>
	  T* search(T2* obj, SVec<double,3>& X, Vec3D& loc) 
  {
    return search<inside, T2, validity>(obj,root, X, loc);
  }

  /* ---------- */

  template <bool (T::*inside)(SVec<double,3>&, const Vec3D&),
            class T2, bool (T2::*validity)(T*)>
	  T* search(T2* obj, Node* n, SVec<double,3>& X, Vec3D& loc) 
  {

    if (loc[0] < n->bbox[0] || loc[0] > n->bbox[1] ||
        loc[1] < n->bbox[2] || loc[1] > n->bbox[3] || 
		 loc[2] < n->bbox[4] || loc[2] > n->bbox[5]) 
	 {
      return NULL;
    }

    if(n->obj) 
	 {
      if (((n->obj)->*inside)(X, loc) && (obj->*validity)(n->obj))
        return n->obj;
      else
        return NULL;
    }

    T* p = search<inside, T2, validity>(obj,n->child1,X,loc);
    if (p) return p;
    p = search<inside, T2, validity>(obj,n->child2,X,loc);
    return p;
  }

  /* ---------- */

private:

  template <void (T::*bbox)(SVec<double,3>& X, double*)>
  Node* constructInternal(SVec<double,3>& X, T** objects, int num) 
  {

    Node* nn = new Node;

    memset(nn, 0, sizeof(Node));

    double boundingbox[6] = {FLT_MAX,-FLT_MAX,FLT_MAX,-FLT_MAX, FLT_MAX,-FLT_MAX};
    double ebb[6];

	  for(int i = 0; i < num; ++i) 
	  {
      (objects[i]->*bbox)(X,ebb);

		  for(int j = 0; j < 3; ++j) 
		  {
        boundingbox[j*2] = std::min(boundingbox[j*2], ebb[j*2]);
        boundingbox[j*2+1] = std::max(boundingbox[j*2+1], ebb[j*2+1]);
      }
    }
    memcpy(nn->bbox, boundingbox, sizeof(boundingbox));

	  if(num == 1) 
	  {
      nn->obj = objects[0];
      return nn;
    }

    double dx[3] = {boundingbox[1]-boundingbox[0],
                    boundingbox[3]-boundingbox[2],
                    boundingbox[5]-boundingbox[4] };
    int dim = 0; 

	  if(dx[1] > dx[0] && dx[1] > dx[2]) 
		  dim = 1;
	  else if (dx[2] > dx[0]) 
		  dim = 2;

    double threshold = boundingbox[dim*2] + 0.5*dx[dim];
     
    std::vector<T*> list1,list2;
	  for(int i = 0; i < num; ++i) 
	  {
		  (objects[i]->*bbox)(X,ebb);

      if ((ebb[dim*2]+ebb[dim*2+1])*0.5 < threshold)
        list1.push_back(objects[i]);
      else
        list2.push_back(objects[i]); 
    }

	  if(list1.size() == 0) 
	  {
      list1.push_back(list2[list2.size()-1]);
      list2.pop_back();
	  } 
	  else if(list2.size() == 0) 
	  {
      list2.push_back(list1[list1.size()-1]);
      list1.pop_back();

    }

    nn->child1 = constructInternal<bbox>(X,&list1[0],list1.size());
    nn->child2 = constructInternal<bbox>(X,&list2[0],list2.size());
   
    return nn; 
  }

  void destructInternal(Node* nn) 
  {
	  if(!nn->obj) 
	  {
      destructInternal(nn->child1);
      destructInternal(nn->child2);
    }
    delete nn;
  }
  
};

