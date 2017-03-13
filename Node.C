#include <Node.h>

#include <BinFileHandler.h>

//------------------------------------------------------------------------------

int NodeSet::read(BinFileHandler &file, int numRanges, int (*ranges)[2], 
		   int *locToGlobMap, int *locToClusMap) 
{

  // read in number of nodes in cluster
  int numClusNodes;
  file.read(&numClusNodes, 1);

  // read in the offset for the first node
  BinFileHandler::OffType start;
  file.read(&start, 1);

  int count = 0;

  // read in ranges
  for (int iRange = 0; iRange < numRanges; ++iRange)  {

    // compute number of nodes in range
    int nNodes = ranges[iRange][1] - ranges[iRange][0] + 1;

    // seek to correct position in file
    file.seek(start + ranges[iRange][0] * (sizeof(int) + 3 * sizeof(double)));
    
    for (int i = 0; i < nNodes; i++) {
      // read in global node number
      file.read(locToGlobMap + count, 1);   
      // compute cluster node number
      locToClusMap[count] = ranges[iRange][0] + i;
      // read in node position  
      file.read(&v[count][0], 3);
      count++;
    }

  }

  if (count != len) {
    fprintf(stderr, "*** Error: wrong number of nodes read (%d instead of %d)\n",
	    count, len);
    exit(1);
  }

  return numClusNodes;

}

//------------------------------------------------------------------------------
