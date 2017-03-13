#include <Connectivity.h>

//------------------------------------------------------------------------------

template<class Map>
void Connectivity::renumberTargets(Map &theMap)
{

  for (int i = 0; i < numtarget; ++i)
    target[i] = theMap[ target[i] ];

}
  
//------------------------------------------------------------------------------
