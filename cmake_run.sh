# /bin/bash

rm -rf CMAKE_BUILD
mkdir CMAKE_BUILD && cd CMAKE_BUILD
cmake ..
sudo make install -j8