#include "common_funcs.hpp"
using namespace std;

#define MODULO 100000

// Initial random seed
int init_seed = 420;

int my_rand(int *x) {
    long long tmp = *x;
    long long lol = 1 << 31;
    long long tmp2 = (tmp * 48271) % (lol - 1);
    *x = int(tmp2);
    return int(tmp2);
}


probs set_pq_values(double a, double b, double c, double d) {
    
    probs ret;

    ret.p00 = a;
    ret.p01 = b;
    ret.p10 = c;
    ret.p11 = d;

    ret.p0X = ret.p00 + ret.p01;
    ret.p1X = ret.p10 + ret.p11;
    ret.p0Y = ret.p00 + ret.p10;
    ret.p1Y = ret.p01 + ret.p11;

    ret.q00 = ret.p0X * ret.p0Y;
    ret.q01 = ret.p0X * ret.p1Y;
    ret.q10 = ret.p1X * ret.p0Y;
    ret.q11 = ret.p1X * ret.p1Y;

    return ret;
}


pair<vector<vector<bool> >, vector<vector<bool> > > make_data(probs P, int N, int T) {

    vector <vector<bool> > X;
    vector <vector<bool> > Y;

    for (int i = 0; i < N; i++) {
    
        vector <bool> x;
        vector <bool> y;
        for (int j = 0; j < T; j++) {
            
            int num = P.p0X * MODULO;
            int rand_num = my_rand(&init_seed) % MODULO;

            if (rand_num <= num) {
                num = (P.p00/P.p0X) * MODULO;
                rand_num = my_rand(&init_seed) % MODULO;
                if (rand_num <= num) y.push_back(0);
                else y.push_back(1);   
                x.push_back(0);
            }
            else {
                num = (P.p10/P.p1X) * MODULO;
                rand_num = my_rand(&init_seed) % MODULO;
                if (rand_num <= num) y.push_back(0);
                else y.push_back(1);   
                x.push_back(1);
            }
        }

        X.push_back(x);
        Y.push_back(y);
    } 

    return make_pair(X, Y);
}


vector <int> get_permut(int T, int seed) {

    vector <int> permut;
    for (int i = 0; i < T; i++) permut.push_back(i);
    
    for (int i = T-1; i > 0; i--) {    
        int rand_num = my_rand(&seed) % i;
        swap(permut[i], permut[rand_num]);
    }
    return permut;
}


void print_probs(probs P) {

    cout<<"P0X = "<<P.p0X<<endl;
    cout<<"P1X = "<<P.p1X<<endl;
    cout<<"P0Y = "<<P.p0Y<<endl;
    cout<<"P1Y = "<<P.p1Y<<endl;
    cout<<"log(q00/p00) = "<<log(P.q00/P.p00)<<endl;
    cout<<"log(q11/p11) = "<<log(P.q11/P.p11)<<endl;
    cout<<endl;
}


void print_results(int N, double r, int b, double expfp, double fp, double tp, int t1, int t2) {

    cout<<N<<"\t";
    cout<<r<<"\t";
    cout<<b<<"\t";
    cout<<expfp<<"\t";
    cout<<fp<<"\t";
    cout<<tp<<"\t";
    cout<<t1<<"\t";
    cout<<t2<<"\t";
    cout<<t1+t2;
    if (tp >= 0.99) cout<<"\tGOOD TP\n";
    else cout<<"\n";
}
