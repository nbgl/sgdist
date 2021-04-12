#include <cmath>

#include "sgdistimpl.h"


template<typename T>
T sgdist(T const* mass1,
         T const* mass2,
         T dx,
         size_t n) {
  T cumul_mass = 0.0;
  T cumul_M = 0.0;
  T cumul_S = 0.0;
  T mass_balance = 0.0;

  size_t i1 = 0;
  size_t i2 = 0;
  while (true) {
    T mass_moved, diff;
    if (mass_balance >= 0) {
      if (i1 >= n) break;
      T new_mass = mass1[i1];
      mass_moved = std::fmin(new_mass, mass_balance);
      mass_balance -= new_mass;
      diff = (long)(i1 - i2);
      ++i1;
    } else {
      if (i2 >= n) break;
      T new_mass = mass2[i2];
      T mass_moved = std::fmin(new_mass, -mass_balance);
      mass_balance += new_mass;
      diff = (long)(i1 - i2);
      ++i2;
    }
    cumul_mass += mass_moved;
    T new_cumul_M = std::fma(mass_moved / cumul_mass, diff - cumul_M, cumul_M);
    cumul_S = std::fma(
      mass_moved, (diff - cumul_M) * (diff - new_cumul_M), cumul_S);
    cumul_M = new_cumul_M;
  }

  return dx * std::sqrt(cumul_S);
}
