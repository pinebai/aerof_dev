//
// Created by lei on 2/3/16.
//

#ifndef _EMBEDDED_ALTERNATING_LEAST_SQUARE_H
#define _EMBEDDED_ALTERNATING_LEAST_SQUARE_H

#include <NonlinearRom.h>
#include <vector>
#include <DistInfo.h>

template<int dim>
class EmbeddedAlternatingLeastSquare : public NonlinearRom<dim> {

public: // TODO: change back to protected when done testing

    VecSet<DistSVec<char, dim> > *mask; //<! create by C++ new when needed
    int nSnapShotFiles;
    int numSnapshots; //<! total number of snapshots across all files, set when read state mask
    std::vector<int> numMasters; //<! numMasters[i] = #{master nodes} in subdomain i
    int reducedDimension;
    // IoData parameters, loaded in constructor
    int maxBasisSize;
    double relativeMinimumEnergy;
    int maxIteration;

    /**
     * Wrapper to NonlinearRom::readSnapshotsFiles
     * @param keyword determines what snapshots is used.
     * @param preprocess defaults to false; see NonlinearRomDatabaseConstruction::constructDatabase()
     * @return the number of files being read.
     * The result is stored in private variable NonlinearRom::snapshots
     */
    int readSnapshotsFilesHelper(char *keyword, bool preprocess = false);

    /**
     * read maskSnapshot from file if existed
     * note that maskSnapshot is of type VecSet<DistVec<bool> >.
     * The result is stored in private variable mask.
     * TODO: read/write vector<bool> seems to produce problem
     */
    int readStateMaskFile();

    /**
     * strips distSVec of slave nodes, put it into a continuous region of memory, pointed by mem
     */
    void freeSlaves(double *&mem, const VecSet<DistSVec<double, dim> > &X, const int M, const int N);
    void freeSlaves(char *&mem, const VecSet<DistSVec<char, dim> > &X, const int M, const int N);

    /**
     * Refills distSVec with mem, a continuous region of memory, put it in X
     */
    void summonSlaves(double *&mem, VecSet<DistSVec<double, dim> > &X, const int M, const int N);
    void summonZombies(double *&mem, VecSet<DistSVec<double, dim> > &X, const int M, const int N);

    /**
     * convert a matrix stored in a continuous region of memory from row-major to column-major
     */
    void transpose(double* &buff1, double* &buff2, int nrow, int ncol);

    /**
     * returns the number of master nodes in a vector
     * res[i] = # of master nodes in subdomain i
     */
    std::vector<int> countMasters(DistInfo &distinfo);

    /**
     * test if two VecSet<DistSVec<double, dim> >& X, Y are equal
     */
    bool isEqualStateMatrices(const VecSet<DistSVec<double, dim> > &X,
                         const VecSet<DistSVec<double, dim> > &Y);
    bool isEqualMaskMatrices(const VecSet<DistSVec<char, dim> > &X,
                         const VecSet<DistSVec<char, dim> > &Y);
    void outputBasis(const VecSet<DistSVec<double, dim> > &U);
    /**
     * read rob matrix, Phi, and reference state U_ref repsectively and return them in paramter
     * used in both offline and online phase:
     */
    void readBasisFiles(VecSet<DistSVec<double, dim> > &U);
    void readReferenceStateFiles(DistSVec<double, dim> &U);

    /**
     * determing the initial rank of SVD by energy criterion
     */
    int initialization(VecSet<DistSVec<double, dim> > &basisInit);

    /**
     * single cluster initialization. instead of
     * kmean(); writeClusteredSnapshots("state");
     * simply call
     * dummyClustering(); writeClusteredSnapshots("state");
     * to set all relevant variables.
     */
    void dummyClustering();

    /**
     * duplicate functionality of constructDatabase from NonlinearRomDatabaseConstruction::constructDatabase
     */
    void constructDatabase();
public:

    EmbeddedAlternatingLeastSquare(Communicator *, IoData &, Domain &/*, DistGeoState * */);

    ~EmbeddedAlternatingLeastSquare();
    using NonlinearRom<dim>::determinePath;
    using NonlinearRom<dim>::qr;
    using NonlinearRom<dim>::readSnapshotFiles;
    using NonlinearRom<dim>::readClusteredSnapshots;
    using NonlinearRom<dim>::readClusterCenters;
    using NonlinearRom<dim>::readNearestSnapsToCenters;
    /**
     * Compute the reduced order basis given snapshot and mask.
     * Input: mask, snapshots
     */
    //TODO: change int dim to double minimum_energy
    void ReducedOrderBasisConstruction();
    void ReducedOrderBasisConstruction(int _dim);
    void ReducedOrderBasisConstructionTesting(int _dim);

    /**
     * testing if snapshot I/O is correct; testing if armadillo compile correctly.
     */
    void testingSnapshotIO();
    void testingALS();
    void testingInitialization();
};


#include "EmbeddedAlternatingLeastSquare.C"

#endif // _EMBEDDED_ALTERNATING_LEAST_SQUARE_H
