//
// Created by Xintong on 4/18/25.
//

#ifndef HLS_YATFHE_KEY_PATTERNS_H
#define HLS_YATFHE_KEY_PATTERNS_H


#include <array>

constexpr std::array<std::array<int, 4>, 4> KEY_PATTERNS2 = {{
      {1, 0, 0, 0},  // 0b00
      {0, 1, 0, 0},  // 0b01
      {0, 0, 1, 0},  // 0b10
      {0, 0, 0, 1}   // 0b11
}};


#endif //HLS_YATFHE_KEY_PATTERNS_H
