//
//  io.hpp
//  minhash
//
//  Created by Liu Cao on 1/31/18.
//  Copyright © 2018 liucao. All rights reserved.
//



#ifndef io_hpp
#define io_hpp

#include <stdio.h>
#include <string>


class prmSrm{
public:
    prmSrm(std::string str_prmFileName);
    ~prmSrm();
//private:
    int n_spectra;
    int * charges;
    double ** prm;
    double ** srm;
    double * precursor_mass;
    int * spectra_length;
    int * spectra_index;

};

#endif /* io_hpp */

