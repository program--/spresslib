#pragma once
#ifndef _SPRESS_DECODE_H_
#define _SPRESS_DECODE_H_

#include <cstdint>
#include <hilbert.hpp>
#include <vector>

#include "spress-extent.hpp"
#include "spress-omp.hpp"
#include "spress-utils.hpp"

namespace spress {

template<typename numeric_t>
void decode(
  size_t             n,
  size_t             k,
  vector<uint16_t>&  h,
  numeric_t          xmax,
  numeric_t          xmin,
  numeric_t          ymax,
  numeric_t          ymin,
  vector<numeric_t>& x,
  vector<numeric_t>& y
)
{
    if (n > 7) {
        // TODO: better handling
        n = 7;
    }

    uint16_t nn  = static_cast<uint16_t>(1 << n);
    size_t   len = h.size() / k;

    if (x.size() < len) {
        x.resize(len);
    }

    if (y.size() < len) {
        y.resize(len);
    }

    // we are essentially iterating through rows
    SAFE_PAR_FOR
    for (size_t i = 0; i < len; i++) {
        numeric_t orig_xmax = xmax, orig_xmin = xmin, orig_ymax = ymax,
                  orig_ymin = ymin;

        // we are getting the extent of all k - 1 indices,
        // but the last one will give the actual coordinates
        SAFE_SIMD
        for (size_t depth = 0; depth < k - 1; depth++) {
            // Get the subextents
            numeric_t sub_xmax, sub_xmin, sub_ymax, sub_ymin;

            extent::subextent(
              nn,
              h[i + (len * depth)],
              orig_xmax,
              orig_xmin,
              orig_ymax,
              orig_ymin,
              sub_xmax,
              sub_xmin,
              sub_ymax,
              sub_ymin
            );

            orig_xmax = sub_xmax;
            orig_xmin = sub_xmin;
            orig_ymax = sub_ymax;
            orig_ymin = sub_ymin;
        }

        // Get coordinates now and add to `x`/`y`
        uint16_t xpos, ypos;
        hilbert::curve::indexToPosition(
          nn, h[i + (len * (k - 1))], &xpos, &ypos
        );
        x[i] = spress::utils::dim_to_coord(nn, xpos, orig_xmax, orig_xmin);
        y[i] = spress::utils::dim_to_coord(nn, ypos, orig_ymax, orig_ymin);
    }
}

} // namespace spress

#endif