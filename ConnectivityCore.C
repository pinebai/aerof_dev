#include <cstdio>
#include <alloca.h>

#include <Connectivity.h>

#include <ResizeArray.h>
#include <BinFileHandler.h>
#include <Face.h>

#ifdef _FEM_CODE_
#include <Utils.d/dofset.h>
#include <Elem.h>
#endif

#define MAX_ALLOCA_SIZE 65536

Connectivity::Connectivity(int _size, int *_pointer, int *_target) 
{
 size      = _size;
 pointer   = _pointer;
 target    = _target;
 numtarget = pointer[size];
 weight    = 0;
}

Connectivity::Connectivity(ElemSet *els)
{
 int i;
 weight = (float *)0;

 size = els->size();

 // Find out the number of targets we will have
 pointer = new int[size+1] ;
 int pp = 0;
 for(i=0; i < size; ++i) {
   pointer[i] = pp;
   pp += (*els)[i] ? (*els)[i].numNodes() : 0;
  }
 pointer[size] = pp;
 numtarget = pp;

 // Create the target array
 target = new int[pp];

 // Fill it in
 for(i=0; i < size; ++i)
  if((*els)[i]) (*els)[i].nodes(target+pointer[i]);
}

//HB
Connectivity::Connectivity(FaceSet* fels)
{
 weight = (float *)0;
                                                                                                    
 size = fels->size();
                                                                                                    
 // Find out the number of targets we will have
 pointer = new int[size+1] ;
 int pp = 0;
 for(int i=0; i < size; ++i) {
   pointer[i] = pp;
   pp += (*fels)[i].numNodes();
  }
 pointer[size] = pp;
 numtarget = pp;
                                                                                                    
 // Create the target array
 target = new int[pp];
                                                                                                    
 // Fill it in
 for(int i=0; i < size; ++i)
  (*fels)[i].nodes(target+pointer[i]);
}

Connectivity::Connectivity(int _size, int *_count)
{
 size    = _size;
 pointer = new int[size+1];
 pointer[0] = 0;
 int i;
 for(i=0; i < _size; ++i)
    pointer[i+1] = pointer[i] + _count[i];
 numtarget = pointer[size];
 target = new int[numtarget];
 weight = 0;
}

Connectivity::~Connectivity()
{
  if (pointer) delete [] pointer;
  if (target) delete [] target;
  if (weight) delete [] weight;
}

// reverse() return a new connectivity that is the reverse of the present one
Connectivity *
Connectivity::reverse(float *w)
{
 Connectivity *res = new Connectivity();

 // The reverse connectivity has the same size as the original
 res->numtarget = numtarget;
 res->target = new int[numtarget];
 res->weight = w;

 // Find the max of target
 int maxtarg = 0;
 int i;
 for(i=0; i < numtarget; ++i)
   if(target[i] > maxtarg) maxtarg = target[i];

 res->size = maxtarg+1;
 res->pointer = new int[res->size+1];
 for(i = 0; i <= res->size; ++i)
   res->pointer[i] = 0;

 // Now do a first pass to fill in res->pointer
 for(i=0; i < numtarget; ++i)
   res->pointer[target[i]]++;

 for(i = 1; i <= res->size; ++i)
   res->pointer[i] += res->pointer[i-1];

 // Second pass fills in target
 int j;
 for(i=0; i < size; ++i) {
   int jend = pointer[i+1];
   for(j = pointer[i]; j < jend; ++j) {
     int tg = target[j];
     res->target[--res->pointer[tg]] = i;
   }
 }

 return res;
}

int
Connectivity::offset(int i, int j)
{
 int ii;
 for(ii = pointer[i]; ii < pointer[i+1];++ii)
  if(target[ii] == j) return ii;

 return -1; // We didn't find a connection between i and j
}

Connectivity* 
Connectivity::transcon( Connectivity* tc)
{
 int i,j,k;

 // First find the biggest target so we can size arrays correctly
 int tgmax=-1;

 for(i =0; i < tc->numtarget; ++i)
   if(tc->target[i] > tgmax) tgmax=tc->target[i];
 tgmax++; // Important adjustment

 // Now we can size the array that flags if a target has been visited
 int *flags = reinterpret_cast<int *>(alloca(sizeof(int)*tgmax));
 
 for(i = 0; i < tgmax; ++i)
   flags[i] = -1;

 // Compute the new pointers
 int *np = new int[size+1];
 int cp = 0;
 for(i = 0; i < size; ++i) {
   np[i] = cp;
   for(j = pointer[i]; j < pointer[i+1]; ++j) {
      int intermed = target[j];
      for (k = 0; k < tc->num(intermed); ++k)
        if(flags[(*tc)[intermed][k]] != i) {
          flags[(*tc)[intermed][k]] = i;
          cp ++;
        }
   }
 }
 np[size] = cp;

 // Now allocate and fill the new target
 for(i = 0; i < tgmax; ++i)
   flags[i] = -1;
 int *ntg = new int[cp];

 cp = 0;

 for(i = 0; i < size; ++i) {
   for(j = pointer[i]; j < pointer[i+1]; ++j) {
     int intermed = target[j];
     for (k = 0; k < tc->num(intermed); ++k)
       if(flags[(*tc)[intermed][k]] != i) {
         flags[(*tc)[intermed][k]] = i;
         ntg[cp] = (*tc)[intermed][k];
         cp ++;
     }
   }
 }

 Connectivity *res = new Connectivity();
 res->size      = size;
 res->pointer   = np;
 res->numtarget = np[size];
 res->target    = ntg;
 res->weight    = 0;
 return res;
}

Connectivity*
Connectivity::transconOne( Connectivity* tc)
// D. Rixen :for every pointer
// of tc only one entry in target of tc is considered
// (used to associate a mpc term for a interface node
//  to only one sub, with a preference for already included sub)
{
 int i,j,k;

 // First find the biggest target so we can size arrays correctly
 int tgmax=-1;

 for(i =0; i < tc->numtarget; ++i)
   if(tc->target[i] > tgmax) tgmax=tc->target[i];
 tgmax++; // Important adjustment

 // Now we can size the array that flags if a target has been visited
 int *flags = reinterpret_cast<int *>(alloca(sizeof(int)*tgmax));

 // For every pointer, build the number of occurrence of targets
 // and choose one target per tc->pointer (max occurrence)
 // At the same time build new pointers np

 int cp = 0;
 int *np = new int[size+1];
 int ii;
 for(i = 0; i < size; ++i) {
   np[i] = cp;
   for(ii = 0; ii < tgmax; ++ii)
    flags[ii] = 0;
   //-- store number of occurence in flag
   for(j = pointer[i]; j < pointer[i+1]; ++j){
      int intermed = target[j];
      for (k = 0; k < tc->num(intermed); ++k)
          flags[(*tc)[intermed][k]]++;
   }
   //-- set pointer and size target
   for(j = pointer[i]; j < pointer[i+1]; ++j){
      int intermed = target[j];
      //-- target with max occ.
      int targMaxOcc;
      int maxOcc = 0;
      for (k = 0; k < tc->num(intermed); ++k){
          if(flags[(*tc)[intermed][k]]==-1){
            maxOcc=0;
            break;
          }
          else if(flags[(*tc)[intermed][k]]> maxOcc){
            maxOcc = flags[(*tc)[intermed][k]];
            targMaxOcc = (*tc)[intermed][k];
          }
      }
      if(maxOcc!=0){cp++;flags[targMaxOcc]=-1;}
    }
 }
 np[size] = cp;
 // Now allocate and fill new target
 int *ntg = new int[cp];
 cp = 0;
 for(i = 0; i < size; ++i) {
   for(ii = 0; ii < tgmax; ++ii)
    flags[ii] = 0;
   //-- store number of occurence in flag
   for(j = pointer[i]; j < pointer[i+1]; ++j){
      int intermed = target[j];
      for (k = 0; k < tc->num(intermed); ++k)
          flags[(*tc)[intermed][k]]++;
   }
   //-- fill target
   for(j = pointer[i]; j < pointer[i+1]; ++j){
      int intermed = target[j];
      //-- target with max occ.
      int targMaxOcc;
      int maxOcc = 0;
      for (k = 0; k < tc->num(intermed); ++k){
          if(flags[(*tc)[intermed][k]]==-1){
            maxOcc=0;
            break;
          }
          else if(flags[(*tc)[intermed][k]]> maxOcc){
            maxOcc = flags[(*tc)[intermed][k]];
            targMaxOcc = (*tc)[intermed][k];
          }
      }
      if(maxOcc!=0){
         ntg[cp]=targMaxOcc;cp++;flags[targMaxOcc]=-1;}
    }
 }

 Connectivity *res = new Connectivity();
 res->size      = size;
 res->pointer   = np;
 res->numtarget = np[size];
 res->target    = ntg;
 res->weight    = 0;
 return res;
}



int
Connectivity::num(int nd, int *mask)
{
 int res=0;
 int jstrt = pointer[nd];
 int jstop = pointer[nd+1];
 int j;
 for(j = jstrt; j < jstop; ++j)
    if(mask[target[j]]) res++;
 return res;
}


void
Connectivity::findPseudoDiam(int *s, int *e, int *mask)
{
 int i,k,nw;
 // Select the node with the lowest connectivity
 int cmin = numtarget+1;
 int cmax = 0;
 for(i = 0; i < size; ++i) {
  if(mask[i]) {
     if((nw = num(i,mask)) < cmin) {
       cmin=nw;
       *s = i;
     }
    if(nw > cmax) cmax = nw;
  }
 }

  //int *ls  = reinterpret_cast<int *>(alloca(numtarget*sizeof(int)));
  //int *xls = reinterpret_cast<int *>(alloca((numtarget+1)*sizeof(int)));
 int*  ls = (numtarget>MAX_ALLOCA_SIZE) ? new int[numtarget  ] : reinterpret_cast<int *>(alloca( numtarget   *sizeof(int)));
 int* xls = (numtarget>MAX_ALLOCA_SIZE) ? new int[numtarget+1] : reinterpret_cast<int *>(alloca((numtarget+1)*sizeof(int)));


 // Created rooted level structure
 int w;
 int h = rootLS(*s, xls, ls, w, mask);
 int subconsize = xls[h];
 int *sorted = reinterpret_cast<int *>(alloca((cmax+1)*sizeof(int)));
 *e = ls[xls[h-1]]; // Give a default end point in case h == subconsize.
 while(h < subconsize) {
   for(k=0; k <= cmax; ++k) sorted[k] = -1;
   // Find in the last level the node with the minimum connectivity
   //int maxweight = subconsize;
   int kstrt = xls[h-1];
   int kstop = xls[h];
   for(k=kstrt; k < kstop; ++k)
      {
        sorted[num(ls[k],mask)] = ls[k];
      }
   int w_e = subconsize;
   for(k = 0; k <= cmax; ++k)
    if(sorted[k] >= 0) {
      int nh = rootLS(sorted[k], xls, ls, w, mask);
      if(w < w_e) {
         if(nh > h) {
           *s = sorted[k];
           h = nh;
           break;
         }
         *e = sorted[k];
         w_e = w;
      }
    }
   if(k > cmax) break;
 }

 if(numtarget > MAX_ALLOCA_SIZE &  ls!=0) delete []  ls;
 if(numtarget > MAX_ALLOCA_SIZE & xls!=0) delete [] xls;

 return;
}

int
Connectivity::rootLS(int root, int *xls, int *ls, int &w, int *mask)
{
 int i, j;
 w = 0;

 int *locMask = new int[size];

 if(mask)
   for(i = 0; i < size; ++i) locMask[i] = mask[i];
 else
   for(i = 0; i < size; ++i) locMask[i] = 1;

 locMask[root] = 0;
 ls[0] = root;
 xls[0] = 0;
 int nlvl = 1;
 int nf = 1;
 while(nf > xls[nlvl-1]) {
   xls[nlvl] = nf;
   int lbegin = xls[nlvl-1];
   int lend   = xls[nlvl];
   for (i = lbegin; i <lend; ++i) {
     int n1 = ls[i];
     int jstart = pointer[n1];
     int jstop  = pointer[n1+1];
     for(j=jstart; j < jstop; ++j) {
       int n2 = target[j];
       if(locMask[n2]) {
          locMask[n2] = 0;
          ls[nf++] = n2;
       }
     }
   }
   if(nf-xls[nlvl] > w) w = nf-xls[nlvl];
   nlvl++;
 }

 delete [] locMask;

 return nlvl-1;
}

#ifdef _FEM_CODE_
compStruct
#else
compStruct *
#endif
Connectivity::renumByComponent(int renumAlg)
{

// size = total number of nodes

  int *globalRenum = new int[size];

  /*

#include <LinkFC.h>
extern "C" { void METIS_NodeND(int *, int *, int *, int *, int *, int *, int *); };

  if (renumAlg == 3) {

    int nnz = pointer[size] - size;

    int *xadj = reinterpret_cast<int *>(alloca((size+1) * sizeof(int)));
    int *adjncy = reinterpret_cast<int *>(alloca(nnz * sizeof(int)));
    int *perm = reinterpret_cast<int *>(alloca(size * sizeof(int)));

    xadj[0] = 0;
    for (int _i=0; _i<size; ++_i) {
      xadj[_i+1] = xadj[_i];
      for (int _j=pointer[_i]; _j<pointer[_i+1]; ++_j) {
	if (target[_j] != _i) {
	  ++xadj[_i+1];
	  adjncy[xadj[_i+1]-1] = target[_j];
	}
      }
    }

    if (xadj[size] != nnz) printf("error for xadj\n");

    int numflag = 0;
    int options[8] = {1, 3, 2, 2, 0, 0, 0, 1};
    METIS_NodeND(&size, xadj, adjncy, &numflag, options, perm, globalRenum);


    compStruct ret;
    ret.numComp = 1;
    ret.xcomp   = 0; 
    ret.renum   = globalRenum;
    
    return ret;

  }
  */

  int *mark        = reinterpret_cast<int *>(alloca(size*sizeof(int)));
  int *ls          = reinterpret_cast<int *>(alloca(size*sizeof(int)));
  int *xls         = reinterpret_cast<int *>(alloca((size+1)*sizeof(int)));

  ResizeArray<int> xcomp(0,2);
// Initialize mark to zero, accounting for missing node #s
// Initialize globalMask

  int inode;
  for(inode = 0; inode < size; ++inode) {
    mark[inode] = (num(inode) != 0) ? 1 : 0;
    globalRenum[inode] = -1;
  }

// Loop over nodes checking which ones are marked
// and belong to the same component.

  int j, k, nextNum = 0, count = 0;
  int *locMask = reinterpret_cast<int *>(alloca(size*sizeof(int)));
  for(inode = 0; inode < size; ++inode)
     locMask[inode] = 0;

  int currentNumber = 0;
  xcomp[0] = currentNumber;

  for(inode = 0; inode < size; ++inode)
  {
    if(mark[inode] == 1)
    {
      // Find all neighbors of inode
 
      int w;
      int h = rootLS(inode, xls, ls, w);

      // Declare and set local mask

      for(j=0; j<xls[h]; ++j) {
          locMask[ls[j]] = 1;
          mark[ls[j]] = 0;
      }

      // call renumbering for local mask
      int *lrenum = reinterpret_cast<int *>(alloca(size*sizeof(int)));
      switch(renumAlg) {
           case 1:
              renumSloan(locMask, nextNum, lrenum);
              break;
           case 2:
              renumRCM(locMask, nextNum, lrenum);
              break;
           default:
              break;
      }

      // reset locMask to zero
      for(j=0; j<xls[h]; ++j)
        locMask[ls[j]] = 0;

      // Assemble local mask into global mask
      if(renumAlg) { 
        for(j=0; j<xls[h]; ++j) {
            k = ls[j];
            globalRenum[k] = lrenum[k];
        }
      }
      else {
        for(j=0; j<xls[h]; ++j) {
           k = ls[j];
           globalRenum[k] = currentNumber + j;
        }
      }
      currentNumber += xls[h];
      count += 1;
      xcomp[count] = currentNumber;
    }
  }

#ifdef _FEM_CODE_
  compStruct ret;
  ret.numComp = count;
  ret.xcomp   = xcomp.yield();
  ret.renum   = globalRenum;
  int ii;
  if(renumAlg == 0) {
     for(ii = 0; ii < size; ++ii) ret.renum[ii] = ii;
  }
#else
  compStruct *ret = new compStruct;
  ret->numComp = count;
  ret->xcomp   = xcomp.yield(); 
  ret->renum   = globalRenum;
  int ii;
  if(renumAlg == 0) {
     for(ii = 0; ii < size; ++ii) ret->renum[ii] = ii;
  }
#endif

  return ret;
}

int *
Connectivity::renumRCM(int *mask, int &nextNum, int *renum)
{
 int i,j,k;
 if(mask == 0)
  {
    mask = reinterpret_cast<int *>(alloca(sizeof(int)*size));
    for(i=0; i < size; ++i)
       mask[i] = (num(i)) ? 1 : 0;
  }

 int s_node, e_node;

 findPseudoDiam( &s_node, &e_node, mask);

 if(renum == NULL) renum = new int[size];
 int *order  = reinterpret_cast<int *>(alloca(sizeof(int)*size));
 int *degree = reinterpret_cast<int *>(alloca(sizeof(int)*size));

 // get the degree of all the nodes
 for(i =0; i < size; ++i)
   degree[i] = num(i,mask);

 order[0] =e_node;
 mask[e_node] = -mask[e_node]; // mark we have seen this node
 int lastNode=1; // number of nodes which have been assigned a number
 for(i = 0; i < lastNode; ++i) {
   int curNode = order[i];
   int firstNeighb = lastNode;
   // Look at the neighbors of this node and add that to the list
   for(j = pointer[curNode]; j < pointer[curNode+1]; ++j) {
     int neighbNode = target[j];
     if(mask[neighbNode] > 0) {
       order[lastNode] = neighbNode;
       mask[neighbNode] = -mask[neighbNode];
       lastNode += 1;
     }
   }
   // now sort the added nodes by degree
   if(firstNeighb < lastNode-1) {
     for(j = firstNeighb; j < lastNode-1; ++j)
      for(k = j+1; k < lastNode; ++k) 
       if(degree[order[k]] < degree[order[j]]) {
         int tmp = order[k];
         order[k] = order[j];
         order[k] = tmp;
       }
   }
 }
 // now reverse the order
 for(i = lastNode; i--; ) {
   renum[order[i]] = (nextNum++);
 }
 return renum;
}


int *
Connectivity::renumSloan(int *mask, int &nextNum, int *renum)
{
 int i,j,k;
 int s_node, e_node;
 float w1=1.0, w2=2.0;

 if(mask == 0)
  {
    mask = reinterpret_cast<int *>(alloca(sizeof(int)*size));
    for(i=0; i < size; ++i)
       mask[i] = (num(i)) ? 1 : 0;
  }

 findPseudoDiam( &s_node, &e_node, mask);

 int *ls  = reinterpret_cast<int *>(alloca(sizeof(int)*numtarget));
 int *xls = reinterpret_cast<int *>(alloca(sizeof(int)*(numtarget+1)));

 int w;
 int h = rootLS(e_node, xls,ls,w,mask);
 // now give a distance to each point

 if(renum == NULL) renum = new int[size];
 int *distance = renum;
 int *status = distance;

 for(i = 0; i < size; ++i)
   distance[i] = -1;

 for(i = 0; i < h; ++i) {
    for(j = xls[i]; j < xls[i+1]; ++j)
      distance[ls[j]] = i;
  }

 float *priority = reinterpret_cast<float *>(alloca(sizeof(float)*size));
 // initialize the priority values
 for(i = 0; i < size; ++i) {
   if(mask[i]) {
     priority[i] = w1*distance[i] - w2*num(i);
     distance[i] = -2; // status and distance are stored in the same array
   }
 }

 // initalize the queue with the starting point
 int *q = reinterpret_cast<int *>(alloca(sizeof(int)*size));
       // maximum size is all the points in the q
 q[0] = s_node;
 int nn=1;
 status[s_node] = -1;

 // While the queue is not empty
 while(nn > 0) {
   // Find in the queue the point with maximum priority
   int next = 0;
   float maxpri = priority[q[next]];
   for(i = 1; i < nn; ++i) {
      if(priority[q[i]] > maxpri) {
         maxpri = priority[q[i]];
         next = i;
      }
   }
   // Remove the next node numbered from the queue
   int nextnode = q[next];
   q[next] = q[nn-1];
   nn = nn-1;
   int istart = pointer[nextnode], istop = pointer[nextnode+1];
   if(status[nextnode] == -1) {
     status[nextnode] = 0; // So we don't step on our feet.
     // Preactive node. Examine its neighbors
     for(i = istart; i < istop; ++i) {
        int neighbor = target[i];
        priority[neighbor] += w2; // update the priority
        if(status[neighbor] == -2) { // if neighbor was inactive, it becomes pre-active
           // NB the next loop will set the same nodes as active
           q[nn] = neighbor;
           nn++;
           status[neighbor] = -1;
        }
     }
   }
   status[nextnode] = nextNum;
   nextNum++;
   // Now scan the preactive neighbors, make them active and their neighbors become
   // preactive if they were inactive
   for(i = istart; i < istop; ++i) {
      int neighbor = target[i];
      if(status[neighbor] == -1) {
        status[neighbor] = 0;
        priority[neighbor] += w2;
        int kstart = pointer[neighbor], kstop = pointer[neighbor+1];
        for(k = kstart; k < kstop ; ++k) {
            int kn = target[k];
            priority[kn] += w2;
            if(status[kn] == -2) { // This node is inactive. must become preactive
               status[kn] = -1;
               q[nn] = kn;
               nn++;
            }
        }
      }
    }
 }
 return status;
}

void
Connectivity::print(FILE *f)
{
 int i, fortran = 0;
 int maxdist = 0;
 fprintf(f, "Connectivity size %d\n",size);
 fflush(f);
 for(i = 0; i < size; ++i) {
   fprintf(f, "%d ->", i+fortran);
   int j;
   for(j = pointer[i]; j < pointer[i+1]; ++j) {
     fprintf(f, " %d", target[j]+fortran);
     if(i-target[j] > maxdist ) maxdist = i-target[j];
   }
   fprintf(f,"\n");
 }
 fprintf(f, "Max dist %d\n", maxdist);
 fflush(f);
}


int
Connectivity::findMaxDist(int *renum)
{
 int i;
 int maxdist = 0;
 int localDist;
 double avgDist=0;
 for(i = 0; i < size; ++i) {
   localDist = 0;
   int j;
   for(j = pointer[i]; j < pointer[i+1]; ++j) {
     if(renum[i]-renum[target[j]] > maxdist )
	 maxdist = renum[i]-renum[target[j]];
     if(renum[i]-renum[target[j]] > localDist)
         localDist = renum[i]-renum[target[j]];
   }
   avgDist += localDist;
 }
 fprintf(stderr, "Total   distance: %e\n", avgDist);
 fprintf(stderr, "Average distance: %e\n", avgDist/size);
 fprintf(stderr, "Maximum distance: %d\n", maxdist);
 return maxdist;
}

#ifdef _FEM_CODE_
int
Connectivity::findProfileSize(EqNumberer *eqn, int unroll)
{
 int i;
 int profileSize = 0;
 int localSize   = 0;
 for(i = 0; i < size; ++i) {
    localSize  = 0;
    int fdof       = eqn->firstdof(i);
    int localWidth = eqn->weight(i);

    int j;
    for(j = pointer[i]; j < pointer[i+1]; ++j) {
      int fjdof = eqn->firstdof(target[j]);
      if(fjdof < 0) continue;

      // This is for skyline unrolling of level unroll
      fjdof = fjdof - (fjdof % unroll);

      if(fdof - fjdof > localSize)
        localSize = fdof - fjdof;
    }
    localSize   *= localWidth;
    localSize   += ((localWidth*(localWidth+1))/2);
    profileSize += localSize;
 }
 fprintf(stderr,"---------------------\n");
 fprintf(stderr, "Profile Size: %d\n", profileSize);

 return profileSize;

}
#endif

Connectivity*
Connectivity::merge(Connectivity *con2)
{
 int size1 = csize();
 int size2 = con2->csize();

 int i;
 int *cp = new int[size1 + size2+1];
 int fp = 0;
 for(i = 0; i < size1; ++i) {
   cp[i] = fp;
   fp += num(i);
 }
 for(i = 0; i < size2; ++i) {
   cp[i+size1] = fp;
   fp += con2->num(i);
 }
 cp[size1 + size2] = fp;

 int *ct = new int[fp];
 fp = 0;
 for(i=0; i<size1; ++i) {
   int j;
   for(j =0; j < num(i); ++j)
     ct[fp++] = (*this)[i][j];
 }
 for(i = 0; i < size2; ++i) {
   int j;
   for(j =0; j < con2->num(i); ++j)
     ct[fp++] = (*con2)[i][j];
 }

 return new Connectivity(size1+size2, cp, ct);
}

void
Connectivity::write(BinFileHandler &file)
{
 file.write(&size, 1);
 file.write(&numtarget, 1);
 file.write(pointer, size+1);
 if (numtarget > 0)
   file.write(target, numtarget);
}

Connectivity::Connectivity(BinFileHandler &file)
{

  file.read(&size, 1);
  file.read(&numtarget, 1);

  pointer = new int[size+1];
  target = new int[numtarget];
  weight = 0;

  file.read(pointer, size+1);
  file.read(target, numtarget);

}

int
Connectivity::numNonZeroP()
{
 int count = 0;
 int i;
 for(i=0; i < size; ++i)
   if(pointer[i+1] != pointer[i]) ++count;
 return count;
}

