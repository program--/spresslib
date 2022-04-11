#pragma once
#ifndef _SPRESS_OMP_H_
#define _SPRESS_OMP_H_

#ifdef _OPENMP
#include <omp.h>
#if _OPENMP >= 201307
#define OMP_VER_4
#endif
#endif

#ifdef OMP_VER_4
#define SAFE_SIMD     _Pragma("omp simd")
#define SAFE_FOR_SIMD _Pragma("omp for simd")
#define SAFE_PAR_FOR  _Pragma("omp parallel for")
#else
#define SAFE_SIMD
#define SAFE_FOR_SIMD
#define SAFE_PAR_FOR
#endif

#endif