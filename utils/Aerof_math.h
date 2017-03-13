#pragma once

template <class T>
inline int aerof_isnan(const T& t) {

  return (t != t);
}
