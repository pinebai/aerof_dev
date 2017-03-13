#pragma once

#include <set>

#ifdef USE_TR1
#include <tr1/unordered_map>
#elif defined(AEROF_WITH_BOOST_UNORDERED_SET)
#include <boost/unordered_set.hpp>
#else
#endif
 
template <class T>
struct Aerof_unordered_set {
  typedef
#ifdef USE_TR1
  std::tr1::unordered_set<T> 
#elif defined(AEROF_WITH_BOOST_UNORDERED_SET)
  boost::unordered_set<T>
#else
  std::set<T>
#endif
  type;
};
