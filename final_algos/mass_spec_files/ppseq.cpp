//
//  ppseq.cpp
//  minhash
//
//  Created by Liu Cao on 2/1/18.
//  Copyright © 2018 liucao. All rights reserved.
//



#include "ppseq.hpp"
#include <string>
#include <fstream>
#include <chrono>
#include <regex>
#include <iostream>
#include "io.hpp"
#include <map>
#include <algorithm>


typedef std::chrono::high_resolution_clock Clock;


void read_ppDatabase(std::string str_fileName, ppseq * &ppDB, int &size){
    std::ifstream inFile(str_fileName);
    std::string line;
    int count = 0;
    
    //std::regex rgx("^>sp\\|([0-9a-zA-Z]+)\\|.*");
    //std::smatch mat;
    
    // count number of peptide
    auto tic = Clock::now();
    while(std::getline(inFile, line)){
        //printf("%d,",count_line);
        if(line.compare(0,1,">") == 0){
            count++;
        }
    }
    auto toc = Clock::now();
    printf("\nTotal number of proteins is %d\n", count);
    printf("Total running time is %d milliseconds.\n",int(std::chrono::duration_cast<std::chrono::milliseconds>(toc - tic).count()));
    inFile.clear();
    inFile.seekg(0, std::ios::beg);
    
    ppDB = new ppseq[count+1];
    size = count;
    
    // read peptide
    tic = Clock::now();
    count = 0;
    std::string protseq("");
    while(std::getline(inFile, line)){
        //printf("%d,",count_line);
        if(line.compare(0,1,">") == 0){
            ppDB[count].length = protseq.size();
            ppDB[count].seq = protseq;
            protseq = std::string("");
            
            count++;
            ppDB[count].name = line.substr(1,line.size()-1);
        }
        else{
            protseq.append(line);
        }
    }
    ppDB[count].length = protseq.size();
    ppDB[count].seq = protseq;
    protseq = std::string("");
    
    toc = Clock::now();
    printf("Total running time of reading peptide is %d milliseconds.\n",int(std::chrono::duration_cast<std::chrono::milliseconds>(toc - tic).count()));
}

void test_read_ppDatabase(std::string str_fileName){
    int size = 0;
    ppseq *ppDB = nullptr;
    read_ppDatabase(str_fileName, ppDB, size);
    for(int i=1;i<size+1;i++){
        std::cout << ppDB[i].name << " " << ppDB[i].seq  << " " << ppDB[i].length << std::endl;
    }
    delete[] ppDB;
}

void pairwise_scoring(prmSrm &prmDB, ppseq * ppDB, int ppDB_size){
    std::map<char,int> dict_aadict {
        {'G',57},
        {'A',71},
        {'S',87},
        {'P',97},
        {'V',99},
        {'T',101},
        {'C',103+57},
        {'I',113},
        {'L',113},
        {'N',114},
        {'D',115},
        {'Q',128},
        {'K',128},
        {'E',129},
        {'M',131},
        {'H',137},
        {'F',147},
        {'R',156},
        {'Y',163},
        {'W',186}
    };
    double sum = 0;
    int pre_mass = 0;
    
    auto tic = Clock::now();
    for(int i=0; i< prmDB.n_spectra; i++){
        for(int j=0; j<ppDB_size;j++){
            for(int k=0;k< ppDB[j].length;k++){
                pre_mass += dict_aadict[ppDB[j].seq[k]];
                if(pre_mass < prmDB.spectra_length[i]){
                    sum += prmDB.prm[i][pre_mass] + prmDB.srm[i][pre_mass];
                }
            }
            pre_mass = 0;
            sum = 0;
        }
        printf("%d,",i);
    }
    auto toc = Clock::now();
    
    printf("Total pairwise scoring running time is %d milliseconds.\n",int(std::chrono::duration_cast<std::chrono::milliseconds>(toc - tic).count()));
    
}

void test_pairwize_scoring(std::string filename1, std::string filename2){
    prmSrm prmDB(filename1);
    int size = 0;
    ppseq *ppDB = nullptr;
    read_ppDatabase(filename2, ppDB, size);
    for(int i=1;i<size+1;i++){
        std::cout << ppDB[i].name << " " << ppDB[i].seq  << " " << ppDB[i].length << std::endl;
    }
    pairwise_scoring(prmDB, ppDB, size);
    delete[] ppDB;
}

void count_seq(std::string str_fileName){
    std::ifstream inFile(str_fileName);
    std::string line;
    int count = 0;
    
    std::regex rgx("^>sp\\|([0-9a-zA-Z]+)\\|.*");
    std::smatch mat;
    
    // count number of protein
    auto tic = Clock::now();
    while(std::getline(inFile, line)){
        //printf("%d,",count_line);
        if(line.compare(0,1,">") == 0){
            count++;
        }
    }
    auto toc = Clock::now();
    printf("\nTotal number of proteins is %d\n", count);
    printf("Total running time is %d milliseconds.\n",int(std::chrono::duration_cast<std::chrono::milliseconds>(toc - tic).count()));
    inFile.clear();
    inFile.seekg(0, std::ios::beg);
    
    ppseq ppDatabase[count+1];
    
    // count the length of each protein
    tic = Clock::now();
    count = 0;
    std::string protseq("");
    while(std::getline(inFile, line)){
        //printf("%d,",count_line);
        if(line.compare(0,1,">") == 0){
            ppDatabase[count].length = protseq.size();
            protseq = std::string("");
            
            count++;
            std::regex_search(line,mat,rgx);
            ppDatabase[count].name = mat[1].str();
        }
        else{
            protseq.append(line);
        }
    }
    ppDatabase[count].length = protseq.size();
    protseq = std::string("");
    
    toc = Clock::now();
    unsigned long total_length =0;
    for(int i=1;i<count+1;i++){
        total_length += ppDatabase[i].length;
    }
    printf("\nTotal length of proteins is %lu\n", total_length);
    printf("Total running time is %d milliseconds.\n",int(std::chrono::duration_cast<std::chrono::milliseconds>(toc - tic).count()));
}

void test_count_seq(std::string str_fileName){
    count_seq(str_fileName);
}
