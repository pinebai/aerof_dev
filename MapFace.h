#ifndef _MAP_FACE_H_
#define _MAP_FACE_H_

#include <algorithm>
#include <map>

struct MaxFace {
    static const int MaxNumNd = 4;
    int nnd; // Number of nodes in this face
    // rotDir indicates if we change the rotation direction of the face
    int rotDir; 
    int nd[MaxNumNd];

    MaxFace(int n, int *d) { 
      nnd = n;
      for(int i=0; i < nnd; ++i)
        nd[i] = d[i];
    }
    
    bool operator<(const MaxFace &f) const {
      if(nnd != f.nnd) return nnd < f.nnd;
      for(int i = 0; i < nnd; ++i)
        if(nd[i] < f.nd[i]) return true;
	else if(nd[i] > f.nd[i]) return false;
      return false;
    }
    // reorder the nodes so that the first one is the smallest one
    // and the rotation order is such that the next one is smaller than the
    // previous one.
    // This routine would not work if the face is degenerate
    // as the solution would not be unique
    void reorder() {
      int iMin =0;
      int i;
      int nd2[MaxNumNd];
      for(i = 0; i < nnd; ++i)
        nd2[i] = nd[i];
      for(i = 1; i < nnd; ++i)
        if(nd[i] < nd[iMin]) iMin = i;
      int im1 = iMin-1;
      if(im1 < 0) im1 += nnd;
      int ip1 = iMin+1;
      if(ip1 >= nnd) ip1 -=nnd;
      if(nd[im1] < nd[ip1])
        rotDir = -1;
      else
        rotDir = 1; 
      for(i = 0; i < nnd; ++i) {
        int j = iMin+i*rotDir;
	if(j >= nnd) j -= nnd;
	if(j < 0) j += nnd;
	nd[i] = nd2[j];
      }
    }
};

typedef std::map<MaxFace, int> MapFaces;

#endif
