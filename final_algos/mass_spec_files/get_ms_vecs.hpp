#ifndef get_ms_vecs_hpp
#define get_ms_vecs_hpp

#include "ppseq.hpp"
#include "io.hpp"
#include "theoretical_util.hpp"

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <chrono>
#include <array>
#include <vector>
#include <list>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <set>
#include <utility>
#include <numeric>

#define WATER_MASS 18.0105647
#define PROTON_MASS 1.00727649
#define PI 3.14159265

using namespace std;

pair<pair<vector<vector<bool> >, vector<vector<bool> > >, vector<string> > get_mass_spec_data(string specfile, string pepfile, vector<int> ids); 

#endif
