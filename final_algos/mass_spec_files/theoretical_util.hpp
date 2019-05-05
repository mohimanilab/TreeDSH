//
//  theoretical_util.hpp
//  wminhash
//
//  Created by Liu Cao on 5/17/18.
//  Copyright © 2018 liucao. All rights reserved.
//

#ifndef theoretical_util_hpp
#define theoretical_util_hpp

#include <stdio.h>
#include <string>

double pep_mass(std::string str_pep);
void test_pep_mass();


void theoretical_vector(int charge, std::string str_pep, double* &prm, double* &srm, double &total_mass);
void test_theoretical_vector();

double dotprod_sparse_dense(double v_sparse[], int n1, double v_dense[], int n2);
void test_dotprod_sparse_dense();

double score1(double prm[], double srm[], int length, int charge, std::string str_pep);
void test_score1();

#endif /* theoretical_util_hpp */
