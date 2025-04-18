//
// Created by Xintong on 4/18/25.
//

#ifndef HLS_YATFHE_KEY_PATTERNS_H
#define HLS_YATFHE_KEY_PATTERNS_H


#include <array>

constexpr std::array<std::array<int, 4>, 4> KEY_PATTERNS2 = {{
      {0, 0, 0, 1},  // 0b00: Encrypt slot 3
      {0, 1, 0, 0},  // 0b01: Encrypt slot 1
      {0, 0, 1, 0},  // 0b10: Encrypt slot 2
      {1, 0, 0, 0}   // 0b11: Encrypt slot 0
}};


#endif //HLS_YATFHE_KEY_PATTERNS_H
