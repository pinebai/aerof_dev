//
// Created by lei on 7/17/16.
//
#include <parser/Assigner.h>
#include "als_io.h"

EmbeddedAlternatingLeastSquareData::EmbeddedAlternatingLeastSquareData() {
    maxBasisSize = 500;
    relativeMinimumEnergy = 0.95;
    maxIteration = 20;
    leastSquareSolver = SVD;
}

EmbeddedAlternatingLeastSquareData::setup(const char *name, ClassAssigner *father){
    ClassAssigner *ca = new ClassAssigner(name, 5, father);

    new ClassInt<EmbeddedAlternatingLeastSquareData>(ca, "MaxBasisSize", this, &EmbeddedAlternatingLeastSquareData::maxBasisSize);
    new ClassDouble<EmbeddedAlternatingLeastSquareData>(ca, "RelativeMinimumEnergy", this, &EmbeddedAlternatingLeastSquareData::relativeMinimumEnergy);
    new ClassInt<EmbeddedAlternatingLeastSquareData>(ca, "MaxIteration", this, EmbeddedAlternatingLeastSquareData::maxIteration);
    new ClassToken<EmbeddedAlternatingLeastSquareData> (ca, "LeastSquareSolver", this, reinterpret_cast<int EmbeddedAlternatingLeastSquareData::*>(&EmbeddedAlternatingLeastSquareData::leastSquareSolver), 2,
                                                        "QR", 0, "SVD", 1);
    snapshots.setup("Snapshots", ca);
}