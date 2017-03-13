//
// Created by lei on 7/17/16.
// input file option handling related to embedded alternating least square.
// separate from IoDataCore.C and IoData.h for better readability.

#ifndef PROJECT_ALS_IO_H
#define PROJECT_ALS_IO_H

// forward declaration: telling compiler to
// not give warning about an undefined class
class ClassAssigner;

struct EmbeddedAlternatingLeastSquareData {
    int maxBasisSize;
    double relativeMinimumEnergy;
    int maxIteration;
    enum LeastSquareSolver {QR = 0, SVD = 1} leastSquareSolver;

    EmbeddedAlternatingLeastSquareData();
    ~EmbeddedAlternatingLeastSquareData() {}

    void setup(const char *, ClassAssigner * = 0);
};

#endif //PROJECT_ALS_IO_H
