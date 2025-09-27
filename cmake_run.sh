#! /bin/bash

rm -rf CMAKE_BUILD
mkdir CMAKE_BUILD && cd CMAKE_BUILD
cmake .. -DINSTALL=OFF
sudo make install -j32