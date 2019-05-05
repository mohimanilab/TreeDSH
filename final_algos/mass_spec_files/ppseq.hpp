//
//  ppseq.hpp
//  minhash
//
//  Created by Liu Cao on 2/1/18.
//  Copyright © 2018 liucao. All rights reserved.
//



#ifndef ppseq_hpp
#define ppseq_hpp

#include <stdio.h>
#include <string>
#include "io.hpp"

struct ppseq{
    std::string seq;
    unsigned long length;
    std::string name;
};

void count_seq(std::string str_fileName);
void test_count_seq(std::string str_fileName);
void read_ppDatabase(std::string str_fileName, ppseq * &ppDB, int &size);
void test_read_ppDatabase(std::string str_fileName);
void pairwise_scoring(prmSrm &prmDB, ppseq * ppDB, int ppDB_size);
void test_pairwize_scoring(std::string filename1, std::string filename2);
#endif /* ppseq_hpp */
