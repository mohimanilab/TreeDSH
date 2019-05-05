#include "ppseq.hpp"
#include "io.hpp"
#include "theoretical_util.hpp"
#include "get_ms_vecs.hpp"
using namespace std;

typedef chrono::high_resolution_clock Clock; 

// important variables
ppseq *ppDB; // peptide database
int num_proteins, minlen = 6, maxlen = 60; // number of proteins, min and max length of peptides
map <int, int> reverse_index; 
unordered_set <int> spectra_masses;
vector <vector <double> > protein_prefix; //obtained in get_protein_prefixes()
vector <vector <int> > protein_prefix_bad; //obtained in get_protein_prefixes()
vector <unordered_set <string> > all_peptides; //obtained in get_all_peptides();


// Amino acid mass table
map <char,double> aa_mass_map {
    {'G',57.02147},
    {'A',71.03712},
    {'S',87.03203},
    {'P',97.05277},
    {'V',99.06842},
    {'T',101.04768},
    {'C',103.00919+57.021},
    {'I',113.08407},
    {'L',113.08407},
    {'N',114.04293},
    {'D',115.02695},
    {'Q',128.05858},
    {'K',128.09497},
    {'E',129.04260},
    {'M',131.04049},
    {'H',137.05891},
    {'F',147.06842},
    {'R',156.10112},
    {'Y',163.06333},
    {'W',186.07932}
};


//Invalid amino acids
bool is_bad_char(char x) {
    if (aa_mass_map.find(x) == aa_mass_map.end()) return true;
    else return false;
}


void get_protein_prefixes() {

    protein_prefix.resize(num_proteins);
    protein_prefix_bad.resize(num_proteins);

    for (int i = 0; i < num_proteins; i++) {
        string seq = ppDB[i].seq;
        int len = ppDB[i].length;

        protein_prefix[i].resize(len + 1);
        protein_prefix_bad[i].resize(len + 1);
        
        if (len == 0) continue; // if protein is empty
        
        for (int j = 1; j <= len; j++) {
            if (is_bad_char(seq[j-1])) {
                protein_prefix[i][j] = protein_prefix[i][j-1];
                protein_prefix_bad[i][j] = protein_prefix_bad[i][j-1] + 1;
            }
            else {
                protein_prefix[i][j] = protein_prefix[i][j-1] + aa_mass_map[seq[j-1]];
                protein_prefix_bad[i][j] = protein_prefix_bad[i][j-1];
            }
        }
    }
}


void get_all_peptides() {

    cout<<"\nObtaining all peptides...\n";
    auto cstart = Clock::now();

    all_peptides.resize(70); // 70 different possible peptide lengths

    for (int len = minlen; len <= maxlen; len++) {
        
        for (int i = 0; i < num_proteins; i++) {

            int protlen = ppDB[i].length;

            for (int j = len; j <= protlen; j++) {
                if (protein_prefix_bad[i][j] == protein_prefix_bad[i][j - len]) {
                    double mass = protein_prefix[i][j] - protein_prefix[i][j - len];
                    if (spectra_masses.find(round(mass)) != spectra_masses.end()) {
                        string pep;
                        for (int k = j - len; k < j; k++) pep.push_back(ppDB[i].seq[k]);
                        all_peptides[len].insert(pep);
                    }
                }
            }
        }

    }

    auto cend = Clock::now();
    cout<<"Obtained all valid peptides from length "<<minlen<<" to "<<maxlen<<" in "<<chrono::duration_cast<chrono::seconds>(cend - cstart).count()<<" seconds\n";
}


pair<pair<vector<vector<bool> >, vector<vector<bool> > >, vector<string> > get_mass_spec_data(string specfile, string pepfile, vector<int> ids) {

    prmSrm ps(specfile);
    read_ppDatabase(pepfile, ppDB, num_proteins);

    vector <vector<bool> > X;
    vector <vector<bool> > Y;
    vector <string> ret_peps;

    // get the indices from the spectra IDs
    for (int i = 0; i < ps.n_spectra; i++) {
        int idx = ps.spectra_index[i];
        reverse_index.insert(make_pair(idx, i));
    }

    // get the maximum spectra length
    int max_spec_len = 0;
    for (int i = 0; i < ids.size(); i++) {
        int idx = reverse_index[ids[i] + 1];
        max_spec_len = max(max_spec_len, ps.spectra_length[idx]);
    }

    // Prepare set X using spectra 
    for (int i = 0; i < ids.size(); i++) {
        
        int idx = reverse_index[ids[i] + 1];
        int charge = ps.charges[idx];
        double pre_mass = ps.precursor_mass[idx];
        spectra_masses.insert(round(charge * (pre_mass - PROTON_MASS) - WATER_MASS));
        
        vector <bool> spec;
        
        int len = ps.spectra_length[idx];
        double *prm = ps.prm[idx];
        double *srm = ps.srm[idx];

        /*        
        vector <pair <double, int> > bestpos;
        for (int j = 0; j < len; j++) bestpos.push_back(make_pair(prm[j], j));
        for (int j = 0; j < len; j++) bestpos.push_back(make_pair(srm[j], j + max_spec_len));
        sort(bestpos.rbegin(), bestpos.rend());

        for (int j = 0; j < 2 * max_spec_len; j++) spec.push_back(0);
        for (int j = 0; j < 60; j++) spec[bestpos[j].second] = 1;
        */

        for (int j = 0; j < len; j++) {
            if (prm[j] > 5.3) spec.push_back(1);
            else spec.push_back(0);
        }
        for (int j = len; j < max_spec_len; j++) spec.push_back(0);

        for (int j = 0; j < len; j++) {
            if (srm[j] > 5.3) spec.push_back(1);
            else spec.push_back(0);
        }
        for (int j = len; j < max_spec_len; j++) spec.push_back(0);

        X.push_back(spec);
    }

    get_protein_prefixes();
    get_all_peptides();

    cout<<"Preparing peptide vectors...\n";
    // Prepare set Y using peptides
    for (int len = minlen; len <= maxlen; len++) {
        
        for (auto it = all_peptides[len].begin(); it != all_peptides[len].end(); ++it) {
            
            string pep = *it;
            
            double* t_prm = nullptr;
            double* t_srm = nullptr;
            double t_mass = 0;
            theoretical_vector(0, pep, t_prm, t_srm, t_mass);
            
            vector <bool> pepvec(2 * max_spec_len);
            for (int i = 0; i < pep.length() - 1; i++) pepvec[t_prm[i]] = 1;
            for (int i = 0; i < pep.length() - 1; i++) pepvec[t_srm[i] + maxlen] = 1;

            Y.push_back(pepvec);
            ret_peps.push_back(pep);
        }
    }

    ofstream of("peptides.txt");
    for (int i = 0; i < ret_peps.size(); i++) of<<ret_peps[i]<<endl;
    of.close();

    return make_pair(make_pair(X, Y), ret_peps);
}
