# /bin/bash

cmake -DCMAKE_INSTALL_PREFIX=/usr/local -B CMAKE_BUILD/
cd CMAKE_BUILD/
sudo make install -j8