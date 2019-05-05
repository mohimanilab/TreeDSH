//
//  theoretical_util.cpp
//  wminhash
//
//  Created by Liu Cao on 5/17/18.
//  Copyright © 2018 liucao. All rights reserved.
//

#include "theoretical_util.hpp"
#include "io.hpp"

#include <map>
#include <cmath>
#include <iostream>
#include <string>

double pep_mass(std::string str_pep){
    // length of peptide
    unsigned long k = str_pep.size();
    // mass dictionary
    std::map<char,double> dict_aadict {
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
    double total_mass = 0;
    for(int i=0;i<k;i++){
        total_mass += round(dict_aadict[str_pep[i]]);
    }
    return total_mass;
}

void test_pep_mass(){
    std::string str = "TVPSAG";
    printf("%f\n", pep_mass(str));
    str = "T";
    printf("%f\n", pep_mass(str));
}

void theoretical_vector(int charge, std::string str_pep, double* &prm, double* &srm, double &total_mass){
    // (M+nH)n+
    // length of peptide
    double mass_H = round(1.00728);
    double mass_H2O = round(18.01056);
    unsigned long k = str_pep.size();
    // mass dictionary
    std::map<char,double> dict_aadict {
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
    // compute total mass
    total_mass = pep_mass(str_pep);
    
    // main output: prm, srm
    prm = new double[k-1];
    srm = new double[k-1];
    std::fill(prm, prm+k-1, 0);
    std::fill(srm, srm+k-1, 0);
    //double b_mass = mass_H*charge;
    //double y_mass = mass_H*charge + mass_H2O + total_mass;
    double b_mass = 0;
    double y_mass = total_mass;
    double m = 0;
    for(int i=0;i<k-1;i++){
        m = round(dict_aadict[str_pep[i]]);
        b_mass += m;
        y_mass -= m;
        prm[i] = b_mass;
        srm[i] = y_mass;
    }
}

void test_theoretical_vector(){
    std::string str = "FEIINAIYEPTEEECEWKPDEEDEISEELKEK";
    auto n = str.size();
    double* prm;
    double* srm;
    int charge = 4;
    double total_mass;
    theoretical_vector(charge, str, prm, srm, total_mass);
    printf("total mass is %f\n\n", total_mass);
    printf("PRM:\n");
    for(int i=0; i<n-1; i++){
        printf("%f\n",prm[i]);
    }
    printf("\nSRM:\n");
    for(int i=0; i<n-1; i++){
        printf("%f\n",srm[i]);
    }
}

double dotprod_sparse_dense(double v_sparse[], int n1, double v_dense[], int n2){
    double score = 0;
    for(int i=0;i<n1;i++){
        score += v_dense[ int(round(v_sparse[i])) ];
    }
    return score;
}

void test_dotprod_sparse_dense(){
    double* v_dense= new double[10];
    std::fill(v_dense,v_dense+10,3);
    double* v_sparse= new double[3];
    std::fill(v_sparse,v_sparse+3,1);
    v_sparse[0] = 3;
    v_sparse[1] = 4;
    v_sparse[2] = 9;
    printf("score is : %f\n", dotprod_sparse_dense(v_sparse,3,v_dense,10));
    delete[] v_sparse;
    delete[] v_dense;
}

double score1(double prm[], double srm[], int length, int charge, std::string str_pep){
    double* theoretical_prm=nullptr;
    double* theoretical_srm=nullptr;
    double theoretical_total_mass = 0;
    unsigned long length_pep = str_pep.size();
    
    theoretical_vector(charge, str_pep, theoretical_prm, theoretical_srm, theoretical_total_mass);
    double s = 0;
    s += dotprod_sparse_dense(theoretical_prm, length_pep-1, prm, length);
    s += dotprod_sparse_dense(theoretical_srm, length_pep-1, srm, length);
    return s;
}

void test_score1(){
    std::string str_prmFilename = "/Users/apple/Desktop/Research/proteomics/data/PRM_plus.mgf";
    prmSrm experimental_spectra(str_prmFilename);
    std::string str_pep = "FEIINAIYEPTEEECEWKPDEEDEISEELKEK";

    for(int i=0;i<experimental_spectra.n_spectra;i++){

        printf("%d ", experimental_spectra.spectra_length[i]);
/*
        if(experimental_spectra.spectra_index[i] == 4404){
            double *prm = experimental_spectra.prm[i];
            int len = experimental_spectra.spectra_length[i];
            printf("\n4404 score: %f\n", score1(experimental_spectra.prm[i], experimental_spectra.srm[i], experimental_spectra.spectra_length[i], experimental_spectra.charges[i], str_pep));
            
            printf("PRM: ");
            for (int j = 0; j < len; j++) {
                printf("%f ", experimental_spectra.prm[i][j]);
            }

            printf("\nPrecursor mass: %f ", experimental_spectra.precursor_mass[i]);
        }

*/
    }
}


