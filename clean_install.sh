#! /bin/bash

rm -rf CMAKE_BUILD

# remove headers
sudo rm -rf /usr/local/include/yatfhe
sudo rm -rf /usr/local/include/yautil

# remove libs
sudo rm -f /usr/local/lib/libyatfhe_lib.a
sudo rm -f /usr/local/lib/libyatfhe_lib.so