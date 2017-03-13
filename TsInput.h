#ifndef TS_INPUT_H
#define TS_INPUT_H

#include <string>

class IoData;

//------------------------------------------------------------------------------

struct TsInput {
  char *solutions;
  char *multiSolutions;
  char *positions;
  char *levelsets;
  char *podFile;
  char *displacements;
//  char *snapFile;
//  char *snapRefSolutionFile;

// Gappy offline
//  char *podFileState;
//  char *podFileRes;
//  char *podFileJac;
//  char *podFileResHat;
//  char *podFileJacHat;

// Gappy online
//  char *sampleNodes;
//  char *jacMatrix;
//  char *resMatrix;
//  char *staterom;
//  char *reducedfullnodemap;
//  char *mesh;

  char *wallsurfacedisplac;
  char *shapederivatives;

  TsInput(IoData &);
  ~TsInput();

  static char * absolutePath(const std::string & rawPath, const std::string & prefix);

private:

  // Disallow copy and assignment
  TsInput(const TsInput &);
  TsInput & operator=(const TsInput &);
};

//------------------------------------------------------------------------------

#endif
