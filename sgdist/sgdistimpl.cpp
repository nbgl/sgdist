#include <cmath>
#include <cstddef>

template<typename T>
T otdist(T const* pos1,
         T const* mass1,
         size_t n1,
         T const* pos2,
         T const* mass2,
         size_t n2) {
  T cumul = 0.0;
  T mass_balance = 0.0;

  size_t i1 = 0;
  size_t i2 = 0;
  T currpos1 = pos1[i1];
  T currpos2 = pos1[i2];
  while (true) {
    T mass_moved;
    if (mass_balance >= 0) {
      if (i1 >= n1) break;
      T new_mass = mass1[i1];
      currpos1 = pos1[i1];
      mass_moved = std::fmin(new_mass, mass_balance);
      mass_balance -= new_mass;
      ++i1;
    } else {
      if (i2 >= n2) break;
      T new_mass = mass2[i2];
      currpos2 = pos2[i2];
      T mass_moved = std::fmin(new_mass, -mass_balance);
      mass_balance += new_mass;
      ++i2;
    }
    T diff = currpos1 - currpos2;
    T diffsq = diff * diff;
    cumul = std::fma(mass_moved, diffsq, cumul);
  }

  return std::sqrt(cumul);
}

template float otdist<float>(
  float const* pos1, float const* mass1, size_t n1,
  float const* pos2, float const* mass2, size_t n2);
template double otdist<double>(
  double const* pos1, double const* mass1, size_t n1,
  double const* pos2, double const* mass2, size_t n2);



template<typename T>
T sgdist(T const* pos1,
         T const* mass1,
         size_t n1,
         T const* pos2,
         T const* mass2,
         size_t n2) {
  T cumul_mass = 0.0;
  T cumul_M = 0.0;
  T cumul_S = 0.0;
  T mass_balance = 0.0;

  size_t i1 = 0;
  size_t i2 = 0;
  T currpos1 = pos1[i1];
  T currpos2 = pos1[i2];
  while (true) {
    T mass_moved;
    if (mass_balance >= 0) {
      if (i1 >= n1) break;
      T new_mass = mass1[i1];
      currpos1 = pos1[i1];
      mass_moved = std::fmin(new_mass, mass_balance);
      mass_balance -= new_mass;
      ++i1;
    } else {
      if (i2 >= n2) break;
      T new_mass = mass2[i2];
      currpos2 = pos2[i2];
      T mass_moved = std::fmin(new_mass, -mass_balance);
      mass_balance += new_mass;
      ++i2;
    }
    cumul_mass += mass_moved;
    T diff = currpos1 - currpos2;
    T new_cumul_M = std::fma(mass_moved / cumul_mass, diff - cumul_M, cumul_M);
    cumul_S = std::fma(
      mass_moved, (diff - cumul_M) * (diff - new_cumul_M), cumul_S);
    cumul_M = new_cumul_M;
  }

  return std::sqrt(cumul_S);
}

template float sgdist<float>(
  float const* pos1, float const* mass1, size_t n1,
  float const* pos2, float const* mass2, size_t n2);
template double sgdist<double>(
  double const* pos1, double const* mass1, size_t n1,
  double const* pos2, double const* mass2, size_t n2);
