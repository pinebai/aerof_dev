#ifndef RESTRICTION_MAPPING_H_
#define RESTRICTION_MAPPING_H_

#include <DistVector.h>

class Domain;

#include <vector>
#include <set>
#include <map>

//------------------------------------------------------------------------------

template <int dim>
class RestrictionMapping {
public:
  int localSubdomainCount() const { return localSubdomainCount_; }
  const DistInfo & originDistInfo() const { return originDistInfo_; }
  const DistInfo & restrictedDistInfo() const { return restrictedDistInfo_; }
	const std::vector<std::vector<int> > & getRestrictedToOriginLocNode() { return restrictedToOrigin_;}

	// from one DistInfo to another
	const DistSVec<double, dim> & restriction(const DistSVec<double, dim>
			&inFull, DistSVec<double, dim> &outRestrict) const;
	const DistSVec<double, dim> & expansion(const DistSVec<double, dim>
			&inRestrict, DistSVec<double, dim> &outFull) const;
 
  double dotProduct(const DistSVec<double, dim> & originVec, const DistSVec<double, dim> & restrictedVec) const;

  template <typename InputIterator>
  RestrictionMapping(Domain * domain, InputIterator globalIndexBegin, InputIterator globalIndexEnd);


  void recomputeConnectedTopology();

private:
  Domain *dom;

  typedef std::map<int, int> NumberingMap;

  int localSubdomainCount_;
  const DistInfo & originDistInfo_;
  DistInfo restrictedDistInfo_;
  
  std::set<int> sampleNodes_;
  std::vector<NumberingMap> originToRestricted_;
  std::vector<std::vector<int> > restrictedToOrigin_;	// local node number in original

  // Disallow copy and assignment
  RestrictionMapping(const RestrictionMapping &); // = delete;
  const RestrictionMapping & operator=(const RestrictionMapping &); // = delete;
};

#ifdef TEMPLATE_FIX
#include <RestrictionMapping.C>
#endif

#endif
