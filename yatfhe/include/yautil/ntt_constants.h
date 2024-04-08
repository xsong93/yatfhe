//
// Created by Xintong Song on 2024/1/15.
//

#ifndef HLS_YATFHE_NTT_CONSTANTS_H
#define HLS_YATFHE_NTT_CONSTANTS_H

#include <cstdint>
#include "yatfhe/torus.h"

// todo: adapt ntt/intt to arbitrary polynomial length

extern NttType phi_normal_2[1024];
extern NttType phi_inverse_2[1024];

extern NttType wb_normal_2[1024];
extern NttType wb_inverse_2[1024];
extern NttType scale_2;

#endif //HLS_YATFHE_NTT_CONSTANTS_H
