#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd $SCRIPT_DIR

export WINEDEBUG=-all

wine mingw/bin/g++.exe -c -o dll.o psim.cpp
rm psim.dll
wine mingw/bin/g++.exe -o psim.dll dll.o -s -shared -Wl,--subsystem,windows

exit 0
