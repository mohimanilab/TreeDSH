#ifndef common_funcs_hpp
#define common_funcs_hpp

#include <iostream>
#include <utility>
#include <vector>
#include <math.h>
#include <unordered_map>
using namespace std;

struct probs_header {
    double p00, p01, p10, p11;
    double p0X, p1X, p0Y, p1Y;
    double q00, q01, q10, q11;
};
typedef struct probs_header probs;

int my_rand(int* x);
probs set_pq_values(double a, double b, double c, double d);
pair<vector<vector<bool> >, vector<vector<bool> > > make_data(probs P, int N, int T);
vector <int> get_permut(int T, int seed);
void print_probs(probs P);
void print_results(int N, double r, int b, double expfp, double fp, double tp, int t1, int t2);

#endif
