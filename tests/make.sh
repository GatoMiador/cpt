#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd $SCRIPT_DIR

export WINEDEBUG=-all

wine mingw/bin/g++.exe -c -o dll.o psim.cpp
rm -f psim.dll
wine mingw/bin/g++.exe -o psim.dll dll.o -s -shared -Wl,--subsystem,windows

wine mingw/bin/g++.exe -c -o dll_mono.o psim_mono.cpp
rm -f psim_mono.dll
wine mingw/bin/g++.exe -o psim_mono.dll dll_mono.o -s -shared -Wl,--subsystem,windows

exit 0
