#include <iostream>
#include <map>
#include <set>
#include <vector>
#include <queue>
#include <math.h>
#include <algorithm>
#include <chrono>
#include "common_funcs.hpp"
#include "mass_spec_files/get_ms_vecs.hpp"
using namespace std;

typedef chrono::high_resolution_clock Clock;

// mass spectrometry related variables
string specfile = "../summer/brownian/PRM_plus.mgf";
string pepfile = "../summer/brownian/uniprot-proteome_homo_sapiens_reference_proteome_5_23_2016.fasta";
vector <int> ids{4910}; 
vector <string> bestpeps{"SGETEDTFIADLVVGLCTGQIK"};

// names of peptides
vector <string> pep_names;

// spectra length
int max_spec_len = 0;

struct hash_tree_node {
    vector <int> vecs;
    vector <bool> path;
    double prob;
    int cur;

    hash_tree_node(vector <int> v, vector <bool> pa, double pr, int c) {
        vecs = v;
        path = pa;
        prob = pr;
        cur = c;
    }
};
typedef struct hash_tree_node htn;

// The 2 sets of vectors that need to be compared
vector <vector <bool> > X;
vector <vector <bool> > Y;

// Length of a vector
int T = 10000;

// Number of vectors in each set
int N;

// The structure of all probabilities
probs P;

// stopping limit and bands
double stopping_limit;
int b;

// random permutation of [0, T-1]
vector <int> permut;

// the key is the path from the root to a particular node in the hash tree. The
// vector pair is set of vectors in X and Y which are hashed to that particular
// path
vector <map <vector <bool>, pair<vector <int>, vector <int> > > > all_paths;

// matches for every vector
vector <set<int> > match_hashes;

// fp and tp rates
double fp, tp, expfp, exptp, alpha, beta;


int totalmctime = 0;

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

void get_expfp() {

    queue <htn> q;
    
    // inserting initial node into queue
    vector <int> init_vec;
    vector <bool> init_path;
    double init_prob = 0;
    int init_cur = 0;
    htn init_node(init_vec, init_path, init_prob, init_cur);
    q.push(init_node);

    while (!q.empty()) {

        htn node = q.front();
        q.pop();

        // node just exceeded the stopping limit
        if (node.prob <= stopping_limit) {
            double curfp = 1;
            double curtp = 1;
            for (int i = 0; i < node.path.size(); i++) {
                if (node.path[i]) {
                    curfp *= (P.p1X * P.p1Y);
                    curtp *= P.p11;
                }
                else {
                    curfp *= (P.p0X * P.p0Y);
                    curtp *= P.p00;
                }
            }
            beta += curfp;
            alpha += curtp;
        }

        else {

            // update paths for the new nodes 
            vector <bool> path0 = node.path;
            path0.push_back(0);
            vector <bool> path1 = node.path;
            path1.push_back(1);

            // new probabilities for the new nodes
            double prob0 = node.prob + log(P.q00/P.p00);
            double prob1 = node.prob + log(P.q11/P.p11);

            // creating new nodes to be pushed into the hash tree
            htn node0(node.vecs, path0, prob0, node.cur + 1);
            htn node1(node.vecs, path1, prob1, node.cur + 1);
            q.push(node0);
            q.push(node1);
        }
    }

    expfp = 1 - pow(1-beta, b);
    exptp = 1 - pow(1-alpha, b);
}


void make_hash_trees(int band, int seed) {

    int veclen = X[0].size();
    permut = get_permut(veclen, seed);

    queue <htn> q;

    // inserting initial node into queue
    vector <bool> init_path;
    double init_prob = 0;
    int init_cur = 0;
    
    vector <int> init_vecX;
    for (int i = 0; i < X.size(); i++) init_vecX.push_back(i);
    htn init_nodeX(init_vecX, init_path, init_prob, init_cur);
    
    vector <int> init_vecY;
    for (int i = 0; i < Y.size(); i++) init_vecY.push_back(i);
    htn init_nodeY(init_vecY, init_path, init_prob, init_cur);
    
    q.push(init_nodeX);

    auto cstart = Clock::now();
    // Preparing the hash tree for set X
    while (!q.empty()) {
        
        htn node = q.front();
        q.pop();

        // node just exceeded the stopping limit
        if (node.prob <= stopping_limit) {
            vector <int> empty_vec;
            all_paths[band].insert(make_pair(node.path, make_pair(node.vecs, empty_vec)));
        }

        else {
            vector <int> node_vecs = node.vecs;
            int idx = node.cur;

            // prepare the vector sets for new nodes
            vector <int> vecs0;
            vector <int> vecs1;
            for (int i = 0; i < node_vecs.size(); i++) {
                if (X[node_vecs[i]][permut[idx]]) vecs1.push_back(node_vecs[i]);
                else vecs0.push_back(node_vecs[i]);
            }
            
            // update paths for the new nodes 
            vector <bool> path0 = node.path;
            path0.push_back(0);
            vector <bool> path1 = node.path;
            path1.push_back(1);

            // new probabilities for the new nodes
            double prob0 = node.prob + log(P.q00/P.p00);
            double prob1 = node.prob + log(P.q11/P.p11);

            // creating new nodes to be pushed into the hash tree
            htn node0(vecs0, path0, prob0, node.cur + 1);
            htn node1(vecs1, path1, prob1, node.cur + 1);
            if (vecs0.size() > 0) q.push(node0);
            if (vecs1.size() > 0) q.push(node1);
        }
    }
    auto cend = Clock::now();
    
    int timespec = chrono::duration_cast<chrono::microseconds>(cend-cstart).count();
    totalmctime += timespec;

    // insert initial node for the hash tree for set Y
    q.push(init_nodeY);

    // Preparing the hash tree for set Y
    while (!q.empty()) {
        
        htn node = q.front();
        q.pop();

        // node just exceeded the stopping limit
        if (node.prob <= stopping_limit) {

            // we only hash the vectors in current node if there exists 
            // vectors from set X with the same path
            if (all_paths[band].find(node.path) != all_paths[band].end())
                all_paths[band][node.path].second = node.vecs;
        }

        else {
            vector <int> node_vecs = node.vecs;
            int idx = node.cur;

            // prepare the vector sets for new nodes
            vector <int> vecs0;
            vector <int> vecs1;
            for (int i = 0; i < node_vecs.size(); i++) {
                if (Y[node_vecs[i]][permut[idx]]) vecs1.push_back(node_vecs[i]);
                else vecs0.push_back(node_vecs[i]);
            }
            
            // update paths for the new nodes 
            vector <bool> path0 = node.path;
            path0.push_back(0);
            vector <bool> path1 = node.path;
            path1.push_back(1);

            // new probabilities for the new nodes
            double prob0 = node.prob + log(P.q00/P.p00);
            double prob1 = node.prob + log(P.q11/P.p11);

            // creating new nodes to be pushed into the hash tree
            htn node0(vecs0, path0, prob0, node.cur + 1);
            htn node1(vecs1, path1, prob1, node.cur + 1);
            if (vecs0.size() > 0) q.push(node0);
            if (vecs1.size() > 0) q.push(node1);
        }
    }

}


void obtain_all_matches() {
    
    match_hashes.resize(N);
    long long total_matches = 0;
    long long real_matches = 0;
    
    // find the pairs created by every hashed path in the two hash trees
    // also calculate dot product only for time analysis
    for (int i = 0; i < b; i++) {
        for (auto it = all_paths[i].begin(); it != all_paths[i].end(); ++it) {
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

    fp = double(total_matches-real_matches) / (double (N) * double(N-1));
    tp = double(real_matches) / double (N);
}

void obtain_msgf_matches() {
    
    cout<<"Obtaining matches...\n";
    
    match_hashes.resize(X.size());
    long long total_matches = 0;
    long long real_matches = 0;
   
    for (int i = 0; i < b; i++) {
        for (auto it = all_paths[i].begin(); it != all_paths[i].end(); ++it) {            
            
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
    for (int i = 0; i < X.size(); i++) {
        total_matches += match_hashes[i].size();
        cout<<total_matches<<endl;

        int maxdp = 0;
        string bestpep;
        bool foundpep = false;

        vector <int> dps;
        for (auto it = match_hashes[i].begin(); it != match_hashes[i].end(); ++it) {
            int j = *it;
            string pep = pep_names[j];
            if (pep == bestpeps[i]) foundpep = true;

            int dp = 0;
            for (int k = 0; k < X[0].size(); k++) {
                if (X[i][k] && Y[j][k]) dp++;
            } 
            dps.push_back(dp);
        }

        
        if (foundpep) cout<<"Found\n";
        else cout<<"Not found\n";
    }
}



void print_dot_products() {

    int same_sum = 0;
    int diff_sum = 0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            int sum = 0;
            for (int k = 0; k < T; k++) {
                if (X[i][k] && Y[j][k]) sum++;
            }
            cout<<i<<" "<<j<<" "<<sum<<endl;
            if (i == j) same_sum += sum;
            else diff_sum += sum;
        }
    }

    cout<<"Same indices average: "<<same_sum/N<<endl;
    cout<<"Different indices average: "<<diff_sum/(N*(N-1))<<endl;
}


void get_only_spectra() {

    prmSrm ps(specfile);

    unordered_map <int, int> reverse_index;

    // get the indices from the spectra IDs
    for (int i = 0; i < ps.n_spectra; i++) {
        int idx = ps.spectra_index[i];
        reverse_index.insert(make_pair(idx, i));
    }

    // get the maximum spectra length
    for (int i = 0; i < ids.size(); i++) {
        int idx = reverse_index[ids[i] + 1];
        max_spec_len = max(max_spec_len, ps.spectra_length[idx]);
    }

    // Prepare set X using spectra 
    for (int i = 0; i < ids.size(); i++) {
        
        int idx = reverse_index[ids[i] + 1];
        
        vector <bool> spec;
        int sum = 0;
        
        int len = ps.spectra_length[idx];
        double *prm = ps.prm[idx];
        double *srm = ps.srm[idx];

        for (int j = 0; j < len; j++) {
            if (prm[j] > 5.3) {
                spec.push_back(1);
                sum++;
            }
            else spec.push_back(0);
        }
        for (int j = len; j < max_spec_len; j++) spec.push_back(0);

        for (int j = 0; j < len; j++) {
            if (srm[j] > 5.3) {
                spec.push_back(1);
                sum++;
            }
            else spec.push_back(0);
        }
        for (int j = len; j < max_spec_len; j++) spec.push_back(0);

        X.push_back(spec);
    }
}


void simulated_data_setup(char** argv) {

    P = set_pq_values(atof(argv[4]), atof(argv[5]), atof(argv[6]), atof(argv[7]));
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
    stopping_limit = atof(argv[3]);

    //cout<<"TREEHASH: Working with b = "<<b<<" and threshold = "<<stopping_limit<<endl<<endl;

    //cout<<"Creating data "<<endl;
    auto cstart = Clock::now();   
    
    // Simulated data
    simulated_data_setup(argv);

    // Mass spectra data
    // msgf_data_setup();    

    get_expfp();
    auto cend = Clock::now();
    //cout<<"Done in "<<chrono::duration_cast<chrono::milliseconds>(cend-cstart).count()<<" ms\n\n";
 
    all_paths.resize(b);
    //cout<<"Creating hash trees "<<endl;
    cstart = Clock::now();
    for (int i = 0; i < b; i++) {
        make_hash_trees(i, 1234 + i * 2345);
    }
    cend = Clock::now();
    //cout<<"Done in "<<chrono::duration_cast<chrono::milliseconds>(cend-cstart).count()<<" ms\n\n";
    int t1 = chrono::duration_cast<chrono::milliseconds>(cend-cstart).count();
    //cout<<"Spectra hash time: "<<totalmctime / b<<endl;


    //cout<<"Creating matches "<<endl;
    cstart = Clock::now();

    // Obtain matches in simulated data
    obtain_all_matches();

    // Obtain matches in msgf data
    // obtain_msgf_matches();
    
    cend = Clock::now();
    //cout<<"Done in "<<chrono::duration_cast<chrono::milliseconds>(cend-cstart).count()<<" ms\n\n";
    int t2 = chrono::duration_cast<chrono::milliseconds>(cend-cstart).count();
    //cout<<t2<<" milliseconds"<<endl;

    print_results(N, stopping_limit, b, expfp, fp, tp, t1, t2);

    return 0;
}
