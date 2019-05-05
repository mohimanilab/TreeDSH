//
//  io.cpp
//  minhash
//
//  Created by Liu Cao on 1/31/18.
//  Copyright © 2018 liucao. All rights reserved.
//

#include "io.hpp"
#include <string>
#include <fstream>
#include <iostream>
#include <chrono>
#include <regex>
#include <cmath>

typedef std::chrono::high_resolution_clock Clock;

prmSrm::prmSrm(std::string str_prmFileName){
    // set up
    std::ifstream inFile(str_prmFileName);
    std::string line;
    int count = 0;
    //int count_line = 0;
    
    //printf("%s", str_prmFileName.c_str());
    // count how many spectra in the file
    auto tic = Clock::now();
    while(std::getline(inFile, line)){
        //printf("%d,",count_line);
        if(!line.compare(0,1,"B")){
            count++;
        }
    }
    auto toc = Clock::now();
    printf("\nTotal number of spectra is %d\n", count/2);
    printf("Total running time is %lf seconds.\n",double(std::chrono::duration_cast<std::chrono::seconds>(toc - tic).count()));
    
    // reset the file handler to the beginning of the file
    inFile.clear();
    inFile.seekg(0, std::ios::beg);
    
    // initialization
    tic = Clock::now();
    this->n_spectra = count / 2;
    this->charges = new int[n_spectra]{};
    this->precursor_mass = new double[this->n_spectra]{};
    this->spectra_length = new int[this->n_spectra]{};
    this->spectra_index = new int[this->n_spectra]{};
    
    this->prm = new double*[this->n_spectra]{};
    this->srm = new double*[this->n_spectra]{};
    
    std::regex rgx("^TITLE=[SP]RM_SpecIndex=([0-9]+) .*");
    std::smatch mat;
    int spect = 0;
    // read PRM and SRM
    for(int i =0; i<this->n_spectra; i++){
        // PRM
        // BEGIN
        std::getline(inFile, line);
        // TITLE
        std::getline(inFile, line);
        std::regex_search(line, mat, rgx);
        spect = stoi(mat[1].str());
        this->spectra_index[i] = spect; // initialzie spect index
        // PEPMASS
        std::getline(inFile, line);
        this->precursor_mass[i] = stod(line.substr(8,line.size()-8)); // initialzie precursor mass
        // SCANS
        std::getline(inFile, line);
        // CHARGE
        std::getline(inFile, line);
        this->charges[i] = stoi(line.substr(7,1));
        this->spectra_length[i] = std::ceil(this->precursor_mass[i]*this->charges[i]); // initialize prm srm length
        this->prm[i] = new double[this->spectra_length[i]]{}; // initialize prm
        this->srm[i] = new double[this->spectra_length[i]]{}; // initialize srm
        // first pair
        std::getline(inFile, line);
        int n_line = 1;
        while(line.compare(0,1,"E") != 0){
            prm[i][n_line] = stod( line.substr(line.find("\t"),line.size()) ); // initialize prm[i]
            n_line ++;
            std::getline(inFile, line);
        }
        //printf("%i,PRM finish\n",i);
        //SRM
        // BEGIN
        std::getline(inFile, line);
        // TITLE
        std::getline(inFile, line);
        // PEPMASS
        std::getline(inFile, line);
        // SCANS
        std::getline(inFile, line);
        // CHARGE
        std::getline(inFile, line);
        // first pair
        std::getline(inFile, line);
        n_line = 1;
        /*
        double *srm_temp = new double[this->spectra_length[i]]{};
        while(line.compare(0,1,"E") != 0){
            srm_temp[n_line] = stod( line.substr(line.find("\t"),line.size()) );
            n_line ++;
            std::getline(inFile, line);
        }
        for(int j=1;j<spectra_length[i];j++){
            //std::cout << int((pm-1)*charge - 18) << ",";
            this->srm[i][ int( round( ((this->precursor_mass[i]-1)*this->charges[i]-18 )*0.9995 )) -j ] = srm_temp[j];
        }
        */
        while(line.compare(0,1,"E") != 0){
            this->srm[i][n_line] = stod( line.substr(line.find("\t"),line.size()) );
            n_line++;
            std::getline(inFile, line);
        }
        //printf("%i,SRM finish\n",i);
        //delete[] srm_temp;
    }
    toc = Clock::now();
    printf("Total initialization time is %lf seconds.\n",double(std::chrono::duration_cast<std::chrono::seconds>(toc - tic).count()));
}

prmSrm::~prmSrm(){
    //for(int i=10;i<11;i++){
    //    for(int j=0;j<spectra_length[i];j++){
    //        std::cout << prm[i][j] << ", " << srm[i][j] << std::endl;
    //    }
    //}
    for(int i =0; i<n_spectra; i++){
        //printf("spectra_index:%d, charge:%d, pm:%lf, spectra_length:%d\n",spectra_index[i],charges[i],precursor_mass[i],spectra_length[i]);
        delete[] prm[i];
        delete[] srm[i];
    }
    delete[] prm;
    delete[] srm;
    delete[] charges;
    delete[] precursor_mass;
    delete[] spectra_index;
    delete[] spectra_length;
}


