#include <BCond.h>
#include <iostream>

//------------------------------------------------------------------------------

void BCondSet::read(BinFileHandler &file)
{
  // read in number of boundary conditions in cluster 
  file.read(&numBC, 1);
  bcs = new BCond[numBC];

  for (int i = 0; i < numBC; i++) {
    // read in cluster node number
    file.read(&bcs[i].nnum, 1);
    // read dof number
    file.read(&bcs[i].dofnum, 1);
    // read value
    file.read(&bcs[i].val, 1);
  }
}

void BCondSet::print()
{
  std::cerr << "numBC = " << numBC << std::endl;
  for(int i=0; i<numBC; ++i) {
    std::cerr << " BC " << i << ": nnum = " << bcs[i].nnum << ", dofnum = " 
              << bcs[i].dofnum << ", val = " << bcs[i].val << std::endl;
  }
}

//------------------------------------------------------------------------------
