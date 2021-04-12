# distutils: language = c++
# distutils: sources = sgdistimpl.cpp

cdef extern from "sgdistimpl.h":
    T sgdist[T](const T* mass1, const T* mass2, T dx, size_t n)
