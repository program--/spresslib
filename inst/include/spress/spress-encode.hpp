#pragma once
#ifndef _SPRESS_ENCODE_H_
#define _SPRESS_ENCODE_H_

#include <cstdint>
#include <hilbert.hpp>
#include <vector>

#include "spress-extent.hpp"
#include "spress-omp.hpp"
#include "spress-utils.hpp"

using std::size_t;
using std::uint16_t;
using std::vector;

namespace spress {

template<typename numeric_t>
void encode(
  size_t             n,
  size_t             k,
  vector<numeric_t>& x,
  vector<numeric_t>& y,
  numeric_t          xmax,
  numeric_t          xmin,
  numeric_t          ymax,
  numeric_t          ymin,
  vector<uint16_t>&  h
)
{
    if (n > 7) {
        // TODO: Better handling
        n = 7;
    }

    uint16_t nn  = static_cast<uint16_t>(1 << n);
    size_t   len = x.size();

    if (h.size() < (k * len)) {
        h.resize(k * len);
    }

    // First Hilbert Order
    SAFE_FOR_SIMD
    for (size_t i = 0; i < len; i++) {
        hilbert::curve::positionToIndex(
          nn,
          spress::utils::coord_to_dim(nn, x[i], xmax, xmin),
          spress::utils::coord_to_dim(nn, y[i], ymax, ymin),
          &h.at(i)
        );
    }

    // Precision depth == 1
    if (k > 1) {
        // for each coordinate
        SAFE_PAR_FOR
        for (size_t i = 0; i < len; i++) {
            // for 1..k-1
            numeric_t orig_xmax = xmax, orig_xmin = xmin, orig_ymax = ymax,
                      orig_ymin = ymin;

            SAFE_SIMD
            for (size_t depth = 1; depth < k; depth++) {
                numeric_t sub_xmax, sub_xmin, sub_ymax, sub_ymin;

                // subextent for parent extent ->
                extent::subextent(
                  nn,
                  h[i + (len * (depth - 1))],
                  orig_xmax,
                  orig_xmin,
                  orig_ymax,
                  orig_ymin,
                  sub_xmax,
                  sub_xmin,
                  sub_ymax,
                  sub_ymin
                );

                hilbert::curve::positionToIndex(
                  nn,
                  spress::utils::coord_to_dim(nn, x[i], sub_xmax, sub_xmin),
                  spress::utils::coord_to_dim(nn, y[i], sub_ymax, sub_ymin),
                  &h.at(i + (len * depth))
                );

                // Prep for next iteration
                orig_xmax = sub_xmax;
                orig_xmin = sub_xmin;
                orig_ymax = sub_ymax;
                orig_ymin = sub_ymin;
            }
        }
    }
}

} // namespace spress

#endif