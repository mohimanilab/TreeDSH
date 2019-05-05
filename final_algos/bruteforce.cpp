#include "common_funcs.hpp"
#include "mass_spec_files/get_ms_vecs.hpp"
using namespace std;

typedef chrono::high_resolution_clock Clock;

// mass spectrometry related variables
string specfile = "../summer/brownian/PRM_plus.mgf";
string pepfile = "../summer/brownian/uniprot-proteome_homo_sapiens_reference_proteome_5_23_2016.fasta";
vector <int> ids{4910}; 
vector <string> bestpeps{"SGETEDTFIADLVVGLCTGQIK"};

vector <vector <bool> > X;
vector <vector <bool> > Y;
vector <string> pep_names;

int max_spec_len = 0;

int main() {
     
    // Mass spectra data
    pair <pair <vector<vector<bool> >, vector<vector<bool> > >, vector <string> > vecs;
    vecs = get_mass_spec_data(specfile, pepfile, ids);
    X = vecs.first.first;
    Y = vecs.first.second;
    pep_names = vecs.second;
    max_spec_len = X[0].size() / 2;
    
    auto cstart = Clock::now(); 
    for (int i = 0; i < X.size(); i++) {
        
        vector <int> dps;
        for (int j = 0; j < Y.size(); j++) {
            int dp = 0;
            for (int k = 0; k < 2 * max_spec_len; k++) {
                if (X[i][k] && Y[j][k]) dp++;
            }
            dps.push_back(dp);
        }

        int maxdp = 0;
        int maxi = 0;
        for (int j = 0; j < Y.size(); j++) {
            if (dps[j] > maxdp) {
                maxdp = dps[j];
                maxi = j;
            }
        }

        if (pep_names[maxi] == bestpeps[i]) cout<<"Found\n";
        else cout<<"Not Found\n";        
    }
    auto cend = Clock::now();

    int timems = chrono::duration_cast<chrono::milliseconds>(cend-cstart).count();
    cout<<timems<<endl;

    return 0;
}
