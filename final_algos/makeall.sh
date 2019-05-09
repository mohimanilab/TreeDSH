#!/bin/bash
echo "Compiling Treehash..."
g++ -std=c++11 -O2 common_funcs.cpp mass_spec_files/*.cpp treehash.cpp -o treehash
echo "Compiling Minhash..."
g++ -std=c++11 -O2 common_funcs.cpp mass_spec_files/*.cpp minhash.cpp -o minhash
echo "Compiling LSH..."
g++ -std=c++11 -O2 common_funcs.cpp mass_spec_files/*.cpp lsh.cpp -o lsh
echo "Compiling Bruteforce..."
g++ -std=c++11 -O2 common_funcs.cpp mass_spec_files/*.cpp bruteforce.cpp -o bruteforce
echo "Done"
