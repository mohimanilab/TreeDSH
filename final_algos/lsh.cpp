#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <set>
#include <map>
#include <unordered_map>
#include <chrono>
#include <string>
#include <memory>
#include "common_funcs.hpp"
#include "mass_spec_files/get_ms_vecs.hpp"

using namespace std;

typedef chrono::high_resolution_clock Clock;

// mass spectrometry related variables
string specfile = "PRM_plus.mgf";
string pepfile = "uniprot-proteome_homo_sapiens_reference_proteome_5_23_2016.fasta";
vector <int> ids{4910}; 
vector <string> bestpeps{"SGETEDTFIADLVVGLCTGQIK"};

// peptide names
vector <string> pep_names;

// spectra length
int max_spec_len = 0;

// The 2 sets of vectors that need to be compared
vector <vector <bool> > X;
vector <vector <bool> > Y;

// Length of a vector
int T = 10000;

// Number of vectors in each set
int N;

// The structure of all probabilities
probs P;

// Hashes and LSH parameters
int r, b;
vector <vector <int> > hashx;
vector <vector <int> > hashy;
vector <map <vector<bool>, pair<vector <int>, vector <int> > > > buckets;
vector <set<int> > match_hashes;

// fp and tp rates;
double fp, tp, expfp;


probs get_ms_pq_values() {
    
    double n00 = 0, n01 = 0, n10 = 0, n11 = 0;
    
    for (int i = 0; i < X.size(); i++) {
        double *t_prm;
        double *t_srm;
        double t_mass;
        theoretical_vector(0, bestpeps[i], t_prm, t_srm, t_mass);

        unordered_set <int> pep_ones;
        for (int j = 0; j < bestpeps[i].length() - 1; j++) pep_ones.insert(t_prm[j]);
        for (int j = 0; j < bestpeps[i].length() - 1; j++) pep_ones.insert(t_srm[j] + max_spec_len);

        for (int j = 0; j < 2 * max_spec_len; j++) {
            if (X[i][j] == 0) {
                if (pep_ones.find(j) == pep_ones.end()) n00 += 1;
                else n01 += 1;
            }
            else {
                if (pep_ones.find(j) == pep_ones.end()) n10 += 1;
                else n11 += 1; 
            }
        }
    }

    n00 /= (X.size() * 2 * max_spec_len);
    n01 /= (X.size() * 2 * max_spec_len);
    n10 /= (X.size() * 2 * max_spec_len);
    n11 /= (X.size() * 2 * max_spec_len);

    return set_pq_values(n00, n01, n10, n11);
}


vector <vector<int> > lsh_hasher(vector<vector<bool> > vecs, int num_hashes) {
    
    vector <vector <int> > hashes;
    int num_vecs = vecs.size();
    int length_vec = vecs[0].size();
    int seed = 1234;

    for (int i = 0; i < num_hashes; i++) {
        int randnum = my_rand(&seed) % length_vec;
        
        vector <int> cur_hashes;
        for (int j = 0; j < num_vecs; j++) {
            cur_hashes.push_back(vecs[j][randnum]);
        }
        hashes.push_back(cur_hashes);
    }

    return hashes;
}


void get_expfp() {

    double s = P.p1X * P.p1Y + P.p0X * P.p0Y;
    expfp = 1 - pow(1 - pow(s, r), b);
}


void hashing_and_bucketing() {

    for (int i = 0; i < X.size(); i++) {
        int sum = 0;
        for (int j = 0; j < X[i].size(); j++) {
            if (X[i][j] == 1) sum++;
        }
    }

    buckets.resize(b);
    hashx = lsh_hasher(X, r * b);
    hashy = lsh_hasher(Y, r * b); 
    
    for (int i = 0; i < X.size(); i++) {
        vector <int> hashes;
        for (int j = 0; j < r * b; j++) hashes.push_back(hashx[j][i]);

        vector <bool> cur_vec;
        int bucket_num = 0;
        for (int j = 0; j < r * b; j++) {
            cur_vec.push_back(hashes[j]);

            if (j % r == r - 1) {
                
                if (buckets[bucket_num].find(cur_vec) == buckets[bucket_num].end()) {
                    vector <int> a{i};
                    vector <int> b;
                    buckets[bucket_num].insert(make_pair(cur_vec, make_pair(a, b)));
                }
                else buckets[bucket_num][cur_vec].first.push_back(i);

                cur_vec.clear();
                bucket_num++;
            }
        }
    }
    
    for (int i = 0; i < Y.size(); i++) {
        vector <int> hashes;
        for (int j = 0; j < r * b; j++) hashes.push_back(hashy[j][i]);
        
        vector <bool> cur_vec;
        int bucket_num = 0;
        for (int j = 0; j < r * b; j++) {
            cur_vec.push_back(hashes[j]);

            if (j % r == r - 1) {
                
                if (buckets[bucket_num].find(cur_vec) != buckets[bucket_num].end())
                    buckets[bucket_num][cur_vec].second.push_back(i);
                
                cur_vec.clear();
                bucket_num++;
            }
        }
    } 
}


void obtain_all_matches() {
     
    match_hashes.resize(N);
    long long total_matches = 0;
    long long real_matches = 0;
    
    for (int i = 0; i < b; i++) {
        for (auto it = buckets[i].begin(); it != buckets[i].end(); ++it) {            
            vector <int> vecsX = it->second.first;
            vector <int> vecsY = it->second.second;
            for (int j = 0; j < vecsX.size(); j++) {
                for (int j2 = 0; j2 < vecsY.size(); j2++) {
                    match_hashes[vecsX[j]].insert(vecsY[j2]);
                }
            }
        }
    }

    // find real and total matches
    for (int i = 0; i < N; i++) {
        total_matches += match_hashes[i].size();
        for (auto it = match_hashes[i].begin(); it != match_hashes[i].end(); ++it) {
            int dp = 0;
            int j = *it;
            for (int t = 0; t < T; t++) {
                if (X[i][t] && Y[j][t]) dp++;
            }
            if (dp > double(0.7 * P.p11 * T) && i == j) real_matches++;
        }
    }
    
    //cout<<total_matches<<" matches found\n";
    //cout<<real_matches<<" real matches found\n";
    
    fp = double(total_matches-real_matches) / (double (N) * double (N-1));
    tp = double(real_matches) / double(N);
}


void obtain_msgf_matches() {
    
    cout<<"Obtaining matches...\n";
    
    match_hashes.resize(X.size());
    long long total_matches = 0;
    long long real_matches = 0;
   
    for (int i = 0; i < b; i++) {
        for (auto it = buckets[i].begin(); it != buckets[i].end(); ++it) {            
            
            vector <bool> thehash = it->first;
            int sum = 0;
            for (int j = 0; j < thehash.size(); j++) sum += thehash[j];
            
            vector <int> vecsX = it->second.first;
            vector <int> vecsY = it->second.second;
            for (int j = 0; j < vecsX.size(); j++) {
                //cout<<vecsX[j]<<"\t";
                for (int j2 = 0; j2 < vecsY.size(); j2++) {
                    match_hashes[vecsX[j]].insert(vecsY[j2]);
                }
            }
            //cout<<vecsY.size()<<"\t";
            //cout<<sum<<endl;
        }
    }

    int max_spec_len = X[0].size()/2; 

    // find real and total matches
    for (int i = 0; i < X.size(); i++) {
        total_matches += match_hashes[i].size();
        cout<<total_matches<<endl;
        for (auto it = match_hashes[i].begin(); it != match_hashes[i].end(); ++it) {
            int j = *it;
            string pep = pep_names[j];
            double* t_prm = nullptr;
            double* t_srm = nullptr;
            double t_mass;
            theoretical_vector(0, pep, t_prm, t_srm, t_mass);
            int dp = 0;
            for (int k = 0; k < pep.length() - 1; k++) dp += X[i][t_prm[k]];
            for (int k = 0; k < pep.length() - 1; k++) dp += X[i][t_srm[k] + max_spec_len];
        }
    }
}


void simulated_data_setup(int argc, char** argv) {

    // Use default probability matrix if values not provided
    if (argc < 8) P = set_pq_values(0.215, 0.0025, 0.255, 0.5275);
    else P = set_pq_values(atof(argv[4]), atof(argv[5]), atof(argv[6]), atof(argv[7]));
    pair <vector<vector<bool> >, vector<vector<bool> > > vecs;
    vecs = make_data(P, N, T);
    X = vecs.first;
    Y = vecs.second;
    return;
}


void msgf_data_setup() {
	
    pair <pair <vector<vector<bool> >, vector<vector<bool> > >, vector <string> > vecs;
    vecs = get_mass_spec_data(specfile, pepfile, ids);
    X = vecs.first.first;
    Y = vecs.first.second;
    pep_names = vecs.second;
    max_spec_len = X[0].size() / 2;
    P = get_ms_pq_values();
    return;
}


int main(int argc, char** argv) {

    if (argc < 4) {
        cout<<"Too few arguments\n";
        return 0;
    }

    N = atoi(argv[1]);
    b = atoi(argv[2]);
    r = atoi(argv[3]);

    //cout<<"Creating data\n";
    auto cstart = Clock::now(); 
    
    // For simulated data 
    simulated_data_setup(argc, argv);

    // For mass spectrometry data
    // msgf_data_setup();
    
    get_expfp();
    auto cend = Clock::now();
    //cout<<"Done in "<<chrono::duration_cast<chrono::milliseconds>(cend-cstart).count()<<" ms\n\n";
    
    //cout<<"Making hashes and putting them in buckets\n";
    cstart = Clock::now();
    hashing_and_bucketing();
    cend = Clock::now();
    //cout<<"Done in "<<chrono::duration_cast<chrono::milliseconds>(cend-cstart).count()<<" ms\n\n"; 
    double t1 = chrono::duration_cast<chrono::milliseconds>(cend-cstart).count();

    //cout<<"Creating matches\n";
    cstart = Clock::now();
    
    // For simulated data
    obtain_all_matches();

    // For mass spectrometry data
    //obtain_msgf_matches();

    cend = Clock::now();
    //cout<<"Done in "<<chrono::duration_cast<chrono::milliseconds>(cend-cstart).count()<<" ms\n\n"; 
    double t2 = chrono::duration_cast<chrono::milliseconds>(cend-cstart).count();

    print_results(N, r, b, expfp, fp, tp, t1, t2);

    return 0;
}
