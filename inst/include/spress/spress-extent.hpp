#pragma once
#ifndef _SPRESS_EXTENT_H_
#define _SPRESS_EXTENT_H_

#include <hilbert.hpp>
#include <vector>

#include "spress-omp.hpp"

using std::size_t;
using std::vector;

namespace spress {
namespace extent {

template<typename numeric_t>
static inline void subextent(
  size_t     nn,
  size_t     h,
  numeric_t  origin_xmax,
  numeric_t  origin_xmin,
  numeric_t  origin_ymax,
  numeric_t  origin_ymin,
  numeric_t& xmax,
  numeric_t& xmin,
  numeric_t& ymax,
  numeric_t& ymin
)
{
    numeric_t xdiff = (origin_xmax - origin_xmin) / (nn - 1),
              ydiff = (origin_ymax - origin_ymin) / (nn - 1);

    size_t x, y;
    hilbert::curve::indexToPosition(nn, h, &x, &y);

    xmax = origin_xmin + (xdiff * x) + xdiff;
    xmin = origin_xmin + (xdiff * x) - xdiff;
    ymax = origin_ymin + (ydiff * y) + ydiff;
    ymin = origin_ymin + (ydiff * y) - ydiff;

    xmax = xmax > origin_xmax ? origin_xmax : xmax;
    xmin = xmin < origin_xmin ? origin_xmin : xmin;
    ymax = ymax > origin_ymax ? origin_ymax : ymax;
    ymin = ymin < origin_ymin ? origin_ymin : ymin;
}

template<typename numeric_t>
static inline void subextents(
  size_t             nn,
  vector<size_t>&    h,
  vector<numeric_t>& origin_xmax,
  vector<numeric_t>& origin_xmin,
  vector<numeric_t>& origin_ymax,
  vector<numeric_t>& origin_ymin,
  vector<numeric_t>& xmax,
  vector<numeric_t>& xmin,
  vector<numeric_t>& ymax,
  vector<numeric_t>& ymin
)
{
    SAFE_FOR_SIMD
    for (size_t i = 0; i < h.size(); i++) {
        subextent(
          nn,
          h[i],
          origin_xmax[i],
          origin_xmin[i],
          origin_ymax[i],
          origin_ymin[i],
          xmax[i],
          xmin[i],
          ymax[i],
          ymin[i]
        );
    }
}

} // namespace extent
} // namespace spress

#endif