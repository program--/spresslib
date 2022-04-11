// clang-format off
#include "spress/spress-decode.hpp"
#include "spress/spress-encode.hpp"
#include <cpp11.hpp>
#include <vector>

using cpp11::doubles;
using cpp11::doubles_matrix;
using cpp11::integers;
using cpp11::integers_matrix;
using std::vector;

[[cpp11::register]]
integers_matrix<> spr_encode_(size_t n, size_t k, doubles x, doubles y, doubles extent)
{
    size_t           len = x.size();
    vector<double>   xx(x.begin(), x.end()), yy(y.begin(), y.end());
    vector<uint16_t> h(k * len);

    spress::encode(
      n,
      k,
      xx,
      yy,
      extent["xmax"],
      extent["xmin"],
      extent["ymax"],
      extent["ymin"],
      h
    );

    cpp11::writable::integers_matrix<> ret(len, k);
    SAFE_PAR_FOR
    for (size_t i = 0; i < len; i++) {
        SAFE_SIMD
        for (size_t j = 0; j < k; j++) {
            ret(i, j) = h[i + (len * j)];
        }
    }

    return ret;
}

[[cpp11::register]]
doubles_matrix<> spr_decode_(size_t n, size_t k, integers h, doubles extent)
{
    size_t           len = h.size() / k;
    vector<uint16_t> hh(h.begin(), h.end());
    vector<double>   x(len), y(len);

    spress::decode(
      n,
      k,
      hh,
      extent["xmax"],
      extent["xmin"],
      extent["ymax"],
      extent["ymin"],
      x,
      y
    );

    cpp11::writable::doubles_matrix<> ret(len, 2);
    SAFE_FOR_SIMD
    for (size_t i = 0; i < len; i++) {
        ret(i, 0) = x[i];
        ret(i, 1) = y[i];
    }

    return ret;
}

// integers_matrix<> spr_encode_sf_(size_t n, size_t k, SEXP sfobj) {
//     
// }