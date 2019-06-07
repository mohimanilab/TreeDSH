//
//  main.cpp
//  DSHTree
//
//  Created by Sean Chang on 4/25/19.
//  Copyright © 2019 Sean Chang. All rights reserved.
//

#include <iostream>
#include <cstdlib>
#include <math.h>
#include <cmath>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <chrono>
#include <algorithm>

typedef std::chrono::high_resolution_clock Clock;

//Keeping track of the number of accepts, rejects, alpha, and beta
int accepts = 0;
int rejects = 0;
double TFP = 0.0;
double TTP = 0.0;
double gx = 0.0;
double gy = 0.0;
int totalX = 0;
int totalY = 0;

//Debugging purposes
int numNodes = 0;
int xnodes = 0;
int ynodes = 0;
int deepest = 0;

//Toggle this boi
double C1;

double p1;
double p2;
double p3;
double p4;

std::vector<bool> found;
int unique;
std::vector<std::vector<bool>> ufp;
int funique;
std::vector<int> checkedbuckets;

//2D vector that holds doubles
typedef std::vector<std::vector<double>> arr;

//p and q are 2D vectors that represent a probability matrix
struct probs_header {
    arr p, q;
    std::vector<double> px, py, qx, qy;
};
typedef struct probs_header probs;

//Special string structure so we can tell which string is paired with which string
struct istring {
    std::string str;
    int id;
};

//Represents a bucket
//xstr and ystr are the strings at that bucket.
//xmatches and ymatches represent the strings that go into the bucket.
struct bucket {
    std::string xstr, ystr;
    std::vector<istring*> xmatches, ymatches;
};

//Tree node
//Contains a vector of its children, what string it represents, and all the buckets at that node
struct tnode {
    std::vector<tnode*> children;
    std::string str;
    std::vector<bucket*> bvec;
};

//Copied from common_funcs.cpp
//Xcode linker why don't you work
#define MODULO 100000

// Initial random seed
int init_seed = 420;

int my_rand(int *x) {
    long long tmp = *x;
    long long lol = 1 << 31;
    long long tmp2 = (tmp * 48271) % (lol - 1);
    *x = int(tmp2);
    return int(tmp2);
}

//Copied from common_funcs.cpp
std::vector <int> get_permut(int T, int seed) {
    
    std::vector <int> permut;
    for (int i = 0; i < T; i++) permut.push_back(i);
    
    for (int i = T-1; i > 0; i--) {
        int rand_num = my_rand(&seed) % i;
        std::swap(permut[i], permut[rand_num]);
    }
    return permut;
}

//Converted from make_data in common_funcs.cpp
//Note: There is some hardcoding going on in here that makes this only work with 2x2 matrices
std::pair<std::vector<istring*>, std::vector<istring*>> make_data(probs *P, int N, int T) {
    
    std::vector<istring*> X;
    std::vector<istring*> Y;
    
    for (int i = 0; i < N; i++) {
        
        istring *x = new istring();
        x -> id = i;
        x -> str = "";
        istring *y = new istring();
        y -> id = i;
        y -> str = "";
        for (int j = 0; j < T; j++) {
            double p0X = ((P -> p).at(0).at(0) + (P -> p).at(0).at(1));
            int num = p0X * MODULO;
            int rand_num = my_rand(&init_seed) % MODULO;
            
            if (rand_num <= num) {
                num = ((P -> p).at(0).at(0)/p0X) * MODULO;
                rand_num = my_rand(&init_seed) % MODULO;
                if (rand_num <= num) y -> str += "0";
                else y -> str += "1";
                x -> str += "0";
            }
            else {
                double p1X = ((P -> p).at(1).at(0) + (P -> p).at(1).at(1));
                num = ((P -> p).at(1).at(0)/p1X) * MODULO;
                rand_num = my_rand(&init_seed) % MODULO;
                if (rand_num <= num) y -> str += "0";
                else y -> str += "1";
                x -> str += "1";
            }
        }
        
        X.push_back(x);
        Y.push_back(y);
    }
    
    return make_pair(X, Y);
}

//Compute the q matrix given the p matrix and set up a probability distribution with them
probs* set_pq_values(arr pmatrix) {
    
    probs *ret = new probs();
    arr qmatrix;
    
    ret -> p = pmatrix;
    
    std::vector<double> rowsums(pmatrix.size(), 0.0);
    std::vector<double> colsums(pmatrix.size(), 0.0);
    
    for(int i = 0; i < pmatrix.size(); i++) {
        for(int j = 0; j < pmatrix.size(); j++) {
            rowsums.at(i) += pmatrix.at(i).at(j);
            colsums.at(j) += pmatrix.at(i).at(j);
        }
    }
    
    ret -> px = rowsums;
    ret -> py = colsums;
    
    for(int i = 0; i < pmatrix.size(); i++) {
        std::vector<double> newrow;
        for(int j = 0; j < pmatrix.size(); j++) {
            newrow.push_back(rowsums.at(i) * colsums.at(j));
        }
        qmatrix.push_back(newrow);
    }
    
    ret -> q = qmatrix;
    
    std::vector<double> qx(qmatrix.size(), 0.0);
    std::vector<double> qy(qmatrix.size(), 0.0);
    
    for(int i = 0; i < qmatrix.size(); i++) {
        for(int j = 0; j < qmatrix.size(); j++) {
            qx.at(i) += qmatrix.at(i).at(j);
            qy.at(j) += qmatrix.at(i).at(j);
        }
    }
    
    ret -> qx = qx;
    ret -> qy = qy;
    
    return ret;
}

//Generate a probability matrix with given size
probs* generateProbability(int size) {
    double *P = new double[size * size];
    double sum = 0;
    for(int i = 0; i < size * size; i++) {
        P[i] = double(rand() / (double)RAND_MAX);
        sum += P[i];
    }
    for(int i = 0; i < size * size; i++) {
        P[i] = P[i] / sum;
    }
    int counter = 0;
    arr pmatrix;
    for(int i = 0; i < size; i++) {
        std::vector<double> newrow;
        for(int j = 0; j < size; j++) {
            newrow.push_back(P[counter]);
            counter++;
        }
        pmatrix.push_back(newrow);
    }
    probs *prob = set_pq_values(pmatrix);
    delete [] P;
    return prob;
}

//Generate 1000 pairs of length 100 strings based on a probability distribution
std::pair<std::vector<istring*>, std::vector<istring*>> generateStrings(probs* prob, int numStrings, int strLength) {
    std::vector<istring*> xs, ys;
    for(int i = 0; i < numStrings; ++i) {
        istring *x = new istring();
        istring *y = new istring();
        x -> id = i;
        y -> id = i;
        std::string xstr = "";
        std::string ystr = "";
        for(int j = 0; j < strLength; j++) {
            double roll = double(rand() / (double)RAND_MAX);
            int newX = 0;
            int newY = 0;
            while(roll > 0) {
                roll -= (prob -> p).at(newX).at(newY);
                if(roll > 0) {
                    if(newX < (prob -> p).size() - 1) {
                        newX++;
                    }
                    else {
                        newX = 0;
                        newY++;
                    }
                }
            }
            xstr += std::to_string((long long int)newX);
            ystr += std::to_string((long long int)newY);
        }
        x -> str = xstr;
        y -> str = ystr;
        xs.push_back(x);
        ys.push_back(y);
    }
    return std::make_pair(xs, ys);
}

//Should we accept as a bucket?
bool accept(int N, double lstar, double P, double Q) {
    return P / Q >= pow(N, 1 + 1 - lstar) * p1;
}

//Should we reject as a bucket?
/*bool reject(double P, double Q, double C, double l1, double l2) {
    //return pow(P, l1) * pow(Q, l2) < 1 / (double) C;
    return P < 1/C1;
}*/

bool reject(double P, double Qx, double Qy, int N, double mu, double nu, double b, double lstar) {
    return P/Qx < pow(N, 1 - lstar) * p2 || P/Qy < pow(N, 1 - lstar) * p3 || /*P < pow(N, 1 - 2 * lstar) * p4*/ pow(P, 1 + mu + nu - b) * pow(Qx, -1 * mu) * pow(Qy, -1 * nu) < pow(N, -1 * lstar) * p4;
}

//Given R and tne corresponding P and Q values, calculate rij
double computerij(double R, double pij, double qij, double l1, double l2) {
    return R * pow(pij, l1) * pow(qij, l2);
}

//We do this computation a lot so it gets it's own function
double Fsum(probs* prob, double mu, double t) {
    double count = 0.0;
    for(int i = 0; i < (prob -> p).size(); i++) {
        for(int j = 0; j < (prob -> p).size(); j++) {
            count += pow((prob -> p).at(i).at(j), mu + 1) * pow((prob -> q).at(i).at(j), t);
        }
    }
    return count;
}

//Approximate t with Newton's method
double approxT(probs* prob, double mu) {
    double t = 0.0;
    
    //Newton's method
    //Shouldn't be anything greater than 30
    for(int i = 0; i < 30; i++) {
        double slope = 0.0;
        for(int k = 0; k < (prob -> p).size(); k++) {
            for(int j = 0; j < (prob -> p).size(); j++) {
                //if((prob -> q).at(i).at(j) != 0) {
                    slope += pow((prob -> p).at(k).at(j), mu + 1) * pow((prob -> p).at(k).at(j), t) * log((prob -> q).at(k).at(j));
                //}
            }
        }
        /*pow(prob.p00, mu + 1) * pow(prob.p00, t) * log(prob.q00) +
        pow(prob.p01, mu + 1) * pow(prob.q01, t) * log(prob.q01) +
        pow(prob.p10, mu + 1) * pow(prob.q10, t) * log(prob.q10) +
        pow(prob.p11, mu + 1) * pow(prob.q11, t) * log(prob.q11);*/
        double current = Fsum(prob, mu, t);
        double deltaT = (current - 1) / slope;
        double t_old = t;
        t -= deltaT;
        //Line search if necessary
        int itr = 0;
        double updated = Fsum(prob, mu, t);
        while(abs(current - 1) < abs(updated - 1) && itr < 1000) {
            t = (t + t_old) / 2;
            updated = Fsum(prob, mu, t);
            itr++;
        }
    }
    return t;
}

//Calculate lambda1 and lambda2 from a probability distribution
std::pair<double, double> calculateLambdas(probs* prob, bool star) {
    double R = 1.0;
    double minval = 0.0;
    double lam1 = 1.0;
    double lam2 = 1.0;
    //Go from 0 to 99 since approxT takes mu and solves p^(mu + 1) * q^t = 1
    for(double i = 0.0; i <= 99.0; i += 0.1) {
        double l2 = approxT(prob, i);
        //std::cout << l2 << "\n";
        arr rijs;
        for(int k = 0; k < (prob -> p).size(); k++) {
            std::vector<double> newrow;
            for(int j = 0; j < (prob -> p).size(); j++) {
                newrow.push_back(computerij(R, (prob -> p).at(k).at(j), (prob -> q).at(k).at(j), i + 1, l2));
                //std::cout <<computerij(R, (prob -> p).at(k).at(j), (prob -> q).at(k).at(j), i + 1, l2) << "\n";
            }
            rijs.push_back(newrow);
        }
        
        //rij * log(rij)
        double val1 = 0.0;
        //rij * log(pij)
        double val2 = 0.0;
        //rij * log(qij)
        double val3 = 0.0;
        
        for(int k = 0; k < (prob -> p).size(); k++) {
            for(int j = 0; j < (prob -> p).size(); j++) {
                val1 += rijs.at(k).at(j) * log(rijs.at(k).at(j));
                val2 += rijs.at(k).at(j) * log((prob -> p).at(k).at(j));
                val3 += rijs.at(k).at(j) * log((prob -> q).at(k).at(j));
            }
        }
        
        //Calculate the max
        //double candidate = (std::max(val1 / (val3 - val1), val1 / (2 * val1 - val2)));
        double l1 = i + 1;
        double candidate = (l1 - l2) / (2 * l1 + l2 - 1);
        
        //Is it lower?
        //std::cout<<val1<<" "<<val2<<" "<<val3<<" "<<candidate<<"\n";
        if(val3 < val1 && 2 * val1 < val2 && candidate > minval) {
            lam1 = i + 1;
            lam2 = l2;
            minval = candidate;
        }
    }
    if(star) {
        //std::cout << lam1 << "," << lam2 << "\n";
        return std::make_pair(minval, minval);
    }
    return std::make_pair(lam1, lam2);
}

/*std::vector<double> newLambda(probs* prob) {
    double lstar = 0.0;
    double mustar = 0.0;
    double nustar = 0.0;
    double bstar = 0.0;
    int psize = (prob -> p).size();
    for(double mun = 0.1; mun <= 10; mun += 0.1) {
        for(double nun = 0.1; nun <= 10; nun += 0.1) {
            for(double bn = 0.1; bn <= 10; bn += 0.1) {
                double mu = mun;
                double nu = nun;
                double b = bn;
                while(mu >= 0 && nu >= 0 && b >= 0) {
                    double numerator = 0.0;
                    double d1 = 0.0;
                    double d2 = 0.0;
                    double d3 = 0.0;
                    for(int i = 0; i < psize; i++) {
                        for(int j = 0; j < psize; j++) {
                            double temp = pow((prob -> p).at(i).at(j), 1 + mu + nu - b) * pow((prob -> px).at(i), -1 * mu) * pow((prob -> px).at(j), -1 * nu);
                            numerator += temp - 1;
                            d1 += temp * log((prob -> p).at(i).at(j)/(prob -> px).at(i));
                            d2 += temp * log((prob -> p).at(i).at(j)/(prob -> px).at(j));
                            d3 += temp * log(1.0/(prob -> p).at(i).at(j));
                        }
                    }
                    double dmu = numerator / d1;
                    double dnu = numerator / d2;
                    double db = numerator / d3;
                    mu -= dmu;
                    nu -= dnu;
                    b -= db;
                    if(dmu < 0.001 && dnu < 0.001 && db < 0.001) {
                        break;
                    }
                }
                if(mu >= 0 && nu >= 0 && b >= 0 && (1 + mu + nu)/(1 + mu + nu - b) > lstar) {
                    double check = 0.0;
                    for(int i = 0; i < psize; i++) {
                        for(int j = 0; j < psize; j++) {
                            check += pow((prob -> p).at(i).at(j), 1 + mu + nu - b) * pow((prob -> px).at(i), -1 * mu) * pow((prob -> px).at(j), -1 * nu);
                        }
                    }
                    check -= 1;
                    std::cout << "Check: " << check << "\n";
                    if(abs(check) <= 0.01) {
                        mustar = mu;
                        nustar = nu;
                        bstar = b;
                        lstar = (1 + mu + nu)/(1 + mu + nu - b);
                    }
                }
            }
        }
    }
    return {mustar, nustar, bstar, lstar};
}*/

double calculateEta(probs *prob, double mu, double nu) {
    double eta = -1;
    double closest = 0.01;
    int psize = (prob -> p).size();
    for(double i = 0.01; i <= std::min(mu, nu); i += 0.01) {
        double sum = 0.0;
        for(int r = 0; r < psize; r++) {
            for(int c = 0; c < psize; c++) {
                sum += pow((prob -> p).at(r).at(c), 1 + mu + nu - i) * pow((prob -> qx).at(r), -1 * mu) * pow((prob -> qx).at(c), -1 * nu);
            }
        }
        if(abs(1 - sum) < closest) {
            closest = abs(1 - sum);
            eta = i;
        }
    }
    return eta;
}

std::vector<double> newLambda(probs *prob) {
    double lstar = 0.0;
    double mustar = 0.0;
    double nustar = 0.0;
    double etastar = 0.0;
    for(double m = 0.01; m <= 10; m += 0.01) {
        for(double n = 0.01; n <= 10; n += 0.01) {
            double e = calculateEta(prob, m, n);
            if(e != -1) {
                double candidate = (1 + m + n)/(1 + m + n - e);
                if(candidate > lstar) {
                    lstar = candidate;
                    mustar = m;
                    nustar = n;
                    etastar = e;
                }
            }
        }
    }
    return {mustar, nustar, etastar, lstar};
}

//The probability that x and y were generated given p
double matchProbability(probs *prob, std::string x, std::string y) {
    /*if(x.length() == 0) {
        return 1.0;
    }
    int xind = std::stoi(x.substr(0, 1));
    int yind = std::stoi(y.substr(0, 1));
    return (prob -> p).at(xind).at(yind) * matchProbability(prob, x.substr(1, x.length() - 1), y.substr(1, y.length() - 1));*/
    double product = 1.0;
    for(int i = 0; i < x.length(); i++) {
        int xind = std::stoi(x.substr(i, 1));
        int yind = std::stoi(y.substr(i, 1));
        product *= (prob -> p).at(xind).at(yind);
    }
    return product;
}

double matchProbabilityQ(probs *prob, std::string x, std::string y) {
    /*if(x.length() == 0) {
     return 1.0;
     }
     int xind = std::stoi(x.substr(0, 1));
     int yind = std::stoi(y.substr(0, 1));
     return (prob -> p).at(xind).at(yind) * matchProbability(prob, x.substr(1, x.length() - 1), y.substr(1, y.length() - 1));*/
    double product = 1.0;
    for(int i = 0; i < x.length(); i++) {
        int xind = std::stoi(x.substr(i, 1));
        int yind = std::stoi(y.substr(i, 1));
        product *= (prob -> q).at(xind).at(yind);
    }
    return product;
}

//Recursive call for generating a tree
void generateTreeHelper(probs* prob, double oldP, double oldQ, double qx, double qy, double gamma_x, double gamma_y, int xind, int yind, int N, std::vector<double> info, tnode *xtree, tnode *ytree) {
    
    //Debugging purposes
    numNodes++;
    
    //Update the probabilities
    double newP = oldP * (prob -> p).at(xind).at(yind);
    double newQ = oldQ * (prob -> q).at(xind).at(yind);
    double newqx = qx * (prob -> qx).at(xind);
    double newqy = qy * (prob -> qy).at(yind);
    double newgamma_x = gamma_x * (prob -> px).at(xind);
    double newgamma_y = gamma_y * (prob -> py).at(yind);
    
    /*NEW CONTROL FLOW
      First create the xtree and ytree nodes if they haven't already been created.
      Check for accept/reject
      Recursive call if not
     */
    //Create new nodes if they don't already exist
    if((xtree -> children).at(xind) == NULL) {
        tnode *newX = new tnode();
        newX -> str = (xtree -> str) + std::to_string((long long int)xind);
        for(int i = 0; i < (prob -> p).size(); i++) {
            (newX -> children).push_back(NULL);
        }
        (xtree -> children).at(xind) = newX;
        xnodes++;
    }
    if((ytree -> children).at(yind) == NULL) {
        tnode *newY = new tnode();
        newY -> str = (ytree -> str) + std::to_string((long long int)yind);
        for(int i = 0; i < (prob -> p).size(); i++) {
            (newY -> children).push_back(NULL);
        }
        (ytree -> children).at(yind) = newY;
        ynodes++;
    }
    
    //These are the nodes that we are looking at in this call
    tnode* currX = (xtree -> children).at(xind);
    tnode* currY = (ytree -> children).at(yind);
    
    //Check for acceptance
    if(accept(N, info.at(3), newP, newQ)) {
        //Create a newbucket
        bucket *nb = new bucket();
        nb -> xstr = currX -> str;
        nb -> ystr = currY -> str;
        //std::cout<<currX -> str<<" "<<currY -> str<<"\n";
        currX -> bvec.push_back(nb);
        currY -> bvec.push_back(nb);
        
        //Update number of accepts, alpha, and beta
        accepts++;
        TFP += newQ;
        TTP += newP;
        
        gx += newgamma_x;
        gy += newgamma_y;
        
        if((currX -> str).length() > deepest) {
            deepest = (currX -> str).length();
        }
        return;
    }
    //Check for rejection
    else if(reject(newP, newqx, newqy, N, info.at(0), info.at(1), info.at(2), info.at(3))) {
        //Update number of rejections
        rejects++;
        return;
    }
    else {
        //Recursive calls.
        for(int i = 0; i < (prob -> p).size(); i++) {
            for(int j = 0; j < (prob -> p).size(); j++) {
                generateTreeHelper(prob, newP, newQ, newqx, newqy, newgamma_x, newgamma_y, i, j, N, info, currX, currY);
            }
        }
        return;
    }
}

//Create the root and recursively generate the tree
std::pair<tnode*, tnode*> generateTree(probs* prob, int N, std::vector<double> info) {
    //Setting up the root stuff
    tnode *xtree = new tnode();
    xtree -> str = "";
    tnode *ytree = new tnode();
    ytree -> str = "";
    for(int i = 0; i < (prob -> p).size(); i++) {
        (xtree -> children).push_back(NULL);
        (ytree -> children).push_back(NULL);
    }
    
    //Recursively generate the tree.
    for(int i = 0; i < (prob -> p).size(); i++) {
        for(int j = 0; j < (prob -> p).size(); j++) {
            generateTreeHelper(prob, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, i, j, N, info, xtree, ytree);
        }
    }
    return std::make_pair(xtree, ytree);
}

//Feed a string into a tree and put it in every bucket along the way to the leaves.
void feedString(istring *istr, std::string str, tnode *tree, bool isX, std::vector<int> *perm) {
    /*if(str.length() == 0) {
        return;
    }
    if(tree == NULL) {
        return;
    }
    std::vector<bucket*> buckets = tree -> bvec;
    for(auto i = buckets.begin(); i != buckets.end(); ++i) {
        if(isX) {
            ((*i) -> xmatches).push_back(istr);
        }
        else {
            ((*i) -> ymatches).push_back(istr);
        }
    }

    int nextInd = std::stoi(str.substr(0, 1));
    feedString(istr, str.substr(1, str.length() - 1), (tree -> children).at(nextInd), isX);*/
    tnode *curr = tree;
    int i = 0;
    while(curr != NULL) {
        std::vector<bucket*> buckets = curr -> bvec;
        for(auto j = buckets.begin(); j != buckets.end(); ++j) {
            if(isX) {
                ((*j) -> xmatches).push_back(istr);
                totalX++;
            }
            else {
                ((*j) -> ymatches).push_back(istr);
                totalY++;
            }
        }
        int nextInd = std::stoi(str.substr((*perm).at(i), 1));
        curr = (curr -> children).at(nextInd);
        i++;
    }
}

//Return true if we should prune the tree
bool pruneTree(tnode* tree) {
    if(tree == NULL) {
        return true;
    }
    bool isLeaf = true;
    for(int i = 0; i < (tree -> children).size(); i++) {
        if((tree -> children).at(i) != NULL) {
            isLeaf = false;
        }
    }
    if(isLeaf) {
        if((tree -> bvec).size() == 0) {
            delete tree;
            xnodes--;
	    return true;
        }
        return false;
    }
    else {
        bool prune = true;
        for(int i = 0; i < (tree -> children).size(); i++) {
            bool pruneChild = pruneTree((tree -> children).at(i));
            if(pruneChild) {
                (tree -> children).at(i) = NULL;
            }
            prune = prune && pruneChild;
        }
        prune = prune && ((tree -> bvec).size() == 0);
        if(prune) {
	    xnodes--;
            delete tree;
        }
        return prune;
    }
}

double gammaprime(tnode *tree, std::vector<double> ps, double gamma) {
    if(tree == NULL) {
        return 0;
    }
    double subgamma = gamma;
    for(int i = 0; i < (tree -> children).size(); i++) {
        subgamma += gammaprime((tree -> children).at(i), ps, gamma * ps.at(i));
    }
    return subgamma;
}

//Feed every string into a tree
std::pair<long long, long long> feedStrings(std::vector<istring*> xstrings, std::vector<istring*> ystrings, tnode *xtree, tnode *ytree, std::vector<int> *perm) {
    auto cstart = Clock::now();
    for(auto i = xstrings.begin(); i != xstrings.end(); ++i) {
        feedString(*i, (*i) -> str, xtree, true, perm);
    }
    auto cend = Clock::now();
    long long xtime = std::chrono::duration_cast<std::chrono::milliseconds>(cend-cstart).count();
    cstart = Clock::now();
    for(auto i = ystrings.begin(); i != ystrings.end(); ++i) {
        feedString(*i, (*i) -> str, ytree, false, perm);
    }
    cend = Clock::now();
    long long ytime = std::chrono::duration_cast<std::chrono::milliseconds>(cend-cstart).count();
    return std::make_pair(xtime, ytime);
}

//Analyze the number of true and false positives in a single bucket
std::pair<long long, long long> analyzeBucket(bucket *b) {
    std::vector<istring*> xs = b -> xmatches;
    std::vector<istring*> ys = b -> ymatches;
    long long TP = 0;
    long long FP = 0;
    
    for(auto i = xs.begin(); i != xs.end(); ++i) {
        istring *x = *i;
        for(auto j = ys.begin(); j != ys.end(); ++j) {
            istring *y = *j;
            if(x -> id == y -> id) {
                if(!found.at(x -> id)) {
                    unique++;
                    found.at(x -> id) = true;
                }
                TP++;
            }
            else {
                if(!ufp.at(x -> id).at(y -> id)) {
                    funique++;
                    ufp.at(x -> id).at(y -> id) = true;
                }
                FP++;
            }
        }
    }
    /*for(auto i = xs.begin(); i != xs.end(); ++i) {
        istring *x = *i;
        std::cout << x -> id << ", ";
    }
    std::cout << "\n";
    for(auto i = ys.begin(); i != ys.end(); ++i) {
        istring *y = *i;
        std::cout << y -> id << ", ";
    }*/
    
    return std::make_pair(TP, FP);
}

std::pair<int, int> _analyzeBucket(bucket *b, probs *prob) {
    std::vector<istring*> xs = b -> xmatches;
    std::vector<istring*> ys = b -> ymatches;
    int TP = 0;
    int FP = 0;
    
    for(auto i = xs.begin(); i != xs.end(); ++i) {
        istring *x = *i;
        for(auto j = ys.begin(); j != ys.end(); ++j) {
            istring *y = *j;
            matchProbability(prob, x -> str, y -> str);
            if(x -> id == y -> id) {
                TP++;
            }
            else {
                FP++;
            }
        }
    }
    
    return std::make_pair(TP, FP);
}

//Analyze all the buckets in a tree (which is all the buckets period).
std::pair<long long, long long> analyzeTree(tnode *tree) {
    if(tree == NULL) {
        return std::make_pair(0, 0);
    }
    long long TP = 0;
    long long FP = 0;
    std::vector<bucket*> buckets = tree -> bvec;
    for(auto i = buckets.begin(); i != buckets.end(); ++i) {
        std::pair<long long, long long> binfo = analyzeBucket(*i);
        TP += binfo.first;
        FP += binfo.second;
    }
    for(int i = 0; i < (tree -> children).size(); i++) {
        std::pair<long long, long long> info = analyzeTree((tree -> children).at(i));
        TP += info.first;
        FP += info.second;
    }
    return std::make_pair(TP, FP);
}

//Analyze all the buckets in a tree (which is all the buckets period).
std::pair<int, int> _analyzeTree(tnode *tree, probs* prob) {
    if(tree == NULL) {
        return std::make_pair(0, 0);
    }
    int TP = 0;
    int FP = 0;
    std::vector<bucket*> buckets = tree -> bvec;
    for(auto i = buckets.begin(); i != buckets.end(); ++i) {
        std::pair<int, int> binfo = _analyzeBucket(*i, prob);
        TP += binfo.first;
        FP += binfo.second;
    }
    for(int i = 0; i < (tree -> children).size(); i++) {
        std::pair<int, int> info = _analyzeTree((tree -> children).at(i), prob);
        TP += info.first;
        FP += info.second;
    }
    return std::make_pair(TP, FP);
}

//Helper function that swaps indicies i and j in string s
std::string swapInd(std::string s, int i, int j) {
    std::string temp = s.substr(i, 1);
    s.replace(i, 1, s.substr(j, 1));
    s.replace(j, 1, temp);
    return s;
}

//Fisher-Yates to shuffle x and y in the same way
std::pair<std::string, std::string> strShuffle(std::string x, std::string y) {
    for(int i = deepest; i > 0; --i) {
        int j = rand()%(i + 1);
        x = swapInd(x, i, j);
        y = swapInd(y, i, j);
    }
    return std::make_pair(x, y);
}

//Permute every pair randomly
void permuteStrings(std::vector<istring*> xs, std::vector<istring*> ys) {
    for(int i = 0; i < xs.size(); ++i) {
        istring *x = xs.at(i);
        istring *y = ys.at(i);
        std::pair<std::string, std::string> newStrs = strShuffle(x -> str, y -> str);
        x -> str = newStrs.first;
        y -> str = newStrs.second;
    }
}

//Clear the xmatches and ymatches from every bucket to reset for a new iteration
void clearBuckets(tnode *tree) {
    if(tree == NULL) {
        return;
    }
    std::vector<bucket*> buckets = tree -> bvec;
    for(auto i = buckets.begin(); i != buckets.end(); ++i) {
        ((*i) -> xmatches).clear();
        ((*i) -> ymatches).clear();
    }
    for(int i = 0; i < (tree -> children).size(); i++) {
        clearBuckets((tree -> children).at(i));
    }
}

void freeTree(tnode* tree, bool freeBuckets) {
    if(tree == NULL) {
        return;
    }
    if(freeBuckets) {
        std::vector<bucket*> bux = tree -> bvec;
        for(int i = 0; i < bux.size(); i++) {
            delete bux.at(i);
        }
    }
    std::vector<tnode*> children = tree -> children;
    for(int i = 0; i < children.size(); i++) {
        freeTree(children.at(i), freeBuckets);
    }
    delete tree;
}

double runTest(probs *prob, int N, std::vector<istring*> xs, std::vector<istring*> ys, std::vector<double> info, double best) {
    TTP = 0;
    TFP = 0;
    deepest = 0;
    accepts = 0;
    rejects = 0;
    gx = 0.0;
    gy = 0.0;
    numNodes = 0;
    xnodes = 0;
    ynodes = 0;
    totalX = 0;
    totalY = 0;
    
    auto cstart = Clock::now();
    //std::cout << "Start\n";
    std::pair<tnode*, tnode*> trees = generateTree(prob, N, info);
    auto cend = Clock::now();
    //std::cout << "End\n";
    //std::cout<<numNodes<<" "<<accepts<<" "<<xnodes<<" "<<ynodes<<"\n";
    
    long long treeGen = std::chrono::duration_cast<std::chrono::milliseconds>(cend-cstart).count();
    
    tnode *xtree = trees.first;
    tnode *ytree = trees.second;
    
    double tb = log(0.01)/log(1 - TTP);
    
    if(std::find(checkedbuckets.begin(), checkedbuckets.end(), accepts) != checkedbuckets.end()) {
        std::cout << p1 << "," << p2 << "," << p3 << "," << p4 << ",checked,checked,checked,checked,checked,checked,checked,checked,checked,checked,checked,checked,checked,checked,checked,checked,checked,checked,checked\n";
        freeTree(xtree, true);
        freeTree(ytree, false);
        return -1;
    }
    
    if(tb > 2000 || treeGen + TFP * tb * N * N * 0.015 > 1.2 * best) {
        //std::cout << "Pruned: " << p1 << " " << p2 << " " << p3 << " " << p4 << "\n";
        std::cout << p1 << "," << p2 << "," << p3 << "," << p4 << ",pruned,pruned,pruned,pruned,pruned,pruned,pruned,pruned,pruned,pruned,pruned,pruned,pruned,pruned,pruned,pruned,pruned,pruned,pruned\n";
        freeTree(xtree, true);
        freeTree(ytree, false);
        return -1;
     }
    
    long long insert = 0;
    long long xtime = 0;
    long long ytime = 0;
    long long analyzeFP = 0;
    
    long long positives = 0;
    
    if(TTP == 0) {
        //std::cout << "Failed on: " << p1 << ", " << p2 << ", " << p3 << "," << p4 << "\n";
        std::cout << p1 << "," << p2 << "," << p3 << "," << p4 << ",failed,failed,failed,failed,failed,failed,failed,failed,failed,failed,failed,failed,failed,failed,failed,failed,failed,failed,failed\n";
        freeTree(xtree, true);
        freeTree(ytree, false);
        return -1;
    }
    

    //std::cout << "Total nodes before pruning: " << xnodes + ynodes << "\n";
    pruneTree(xtree);
    pruneTree(ytree);
    //std::cout << "Total nodes after pruning: " << xnodes + ynodes << "\n";
    
    double gammaprimex = gammaprime(xtree, prob -> px, 1.0);
    double gammaprimey = gammaprime(ytree, prob -> py, 1.0);
    
    
    /*if(tb > 1000 || (TFP) * N * N * tb > 50000) {
        std::cout << "Pruning: " << p1 << ", " << p2 << ", " << p3 << "," << p4 << "\n";
        return -1;
    }*/
    
    double b = 0.0;
    
    double tprate = 1.0;
    
    int strl = (xs.at(0) -> str).length();
    double ETP = 0.0;
    double EFP = 0.0;
    
    //std::cout << "Theoretical b: " << tb << "\n";
    
    //for(double i = 0.0; i < tb; i++) {
    while(tprate > 0.01) {
    //while(unique < 2700) {
        cstart = Clock::now();
        std::vector<int> *perm = new std::vector<int>();
        *perm = get_permut(strl, 1234 + b * 2345);
        
        /*for(int j = 0; j < strl; j++) {
            (*perm).push_back(j);
        }
        
        for(int k = strl - 1; k > 0; --k) {
            int l = rand()%(k + 1);
            int temp = (*perm).at(k);
            (*perm).at(k) = (*perm).at(l);
            (*perm).at(l) = temp;
        }*/
        
        std::pair<long long, long long> times = feedStrings(xs, ys, xtree, ytree, perm);
        cend = Clock::now();
        insert += std::chrono::duration_cast<std::chrono::milliseconds>(cend-cstart).count();
        xtime += times.first;
	ytime += times.second;
	cstart = Clock::now();
        std::pair<long long, long long> info = analyzeTree(xtree);
        cend = Clock::now();
        long long TP = info.first;
        long long FP = info.second;
        positives += TP + FP;
        ETP += TP;
        EFP += FP;
        //std::cout << TP << " " << FP << "\n";
        analyzeFP += 0.015 * (TP + FP);
        tprate *= 1.0 - ((double)TP/(double)N);
        /*cstart = Clock::now();
        permuteStrings(xs, ys);
        cend = Clock::now();
        permute += std::chrono::duration_cast<std::chrono::milliseconds>(cend-cstart).count();*/
        /*std::cout << "Insert: " << insert << "\n";
        std::cout << TP << " " << FP << "\n";
        std::cout << "Positives: " << positives << "\n";
        std::cout << "Time so far: " << treeGen + insert + 0.015 * positives << "\n";*/
        //std::cout << b << ": " << unique << "\n";
        delete perm;
        //std::cout << "b: " << b << ", tprate: " << tprate << "\n";
        b++;
        clearBuckets(xtree);
    }
    
    /*std::cout << "Bands: " << b << "\n";
    std::cout << "Deepest: " << deepest << "\n";
    std::cout << "Tree generation: " << treeGen << "ms\n";
    std::cout << "Data insertion: " << insert << "ms\n";
    std::cout << "Analyze FP: " << analyzeFP << "ms\n";
    std::cout << "Permute: " << permute << "ms\n";
    std::cout << "Total: " << treeGen + insert + analyzeFP + permute << "ms\n";*/
    //std::cout << p1 << "," << p2 << "," << p3 << "," << p4 << "," << treeGen + insert + 0.015 * positives << "," << insert << "," << treeGen << "," << 0.015 * positives << "," << tb << "," << ETP << "," << ETP/(b * 100000 * 100000) << "," << TTP << "," << funique << "," << 1 - pow(1 - (funique/(3000.0 * 2999)), 1/b) << "," << TFP << "," << accepts << "," << numNodes << "," << xnodes << "," << ynodes << "," << gx << "," << gy << "," << totalX << "," << totalY << "\n";
    std::cout << p1 << "," << p2 << "," << p3 << "," << p4 << "," << treeGen + insert + 0.015 * positives << "," << insert << "," << treeGen << "," << 0.015 * positives << "," << b << "," << ETP << "," << 1 - pow(tprate, 1/b) << "," << TTP << "," << EFP << "," << EFP/(b * N * N) << "," << TFP << "," << accepts << "," << numNodes << "," << xnodes << "," << ynodes << "," << gx << "," << gy << "," << totalX << "," << totalY << "\n";
    //std::cout << "Ratio: " << insert/(b * gammaprimex * N + b * gammaprimey * N) << "\n";
    //std::cout << "X Ratio: " << xtime / (b * N * (gx + 1)) << " Y Ratio: " << ytime / (b * N * (gy + 1)) << "\n";
    //std::cout << "GammaX: " << gx << " GammaPrimeX: " << gammaprimex << " GammaY: " << gy << " GammaPrimeY: " << gammaprimey << "\n";
    freeTree(xtree, true);
    freeTree(ytree, false);
    for(int i = 0; i < found.size(); i++) {
        found.at(i) = false;
    }
    for(int i = 0; i < ufp.size(); i++) {
        for(int j = 0; j < ufp.size(); j++) {
            ufp.at(i).at(j) = false;
        }
    }
    /*std::cout << positives << "\n";
    std::cout << treeGen << "\n";
    std::cout << insert<< "\n";
    std::cout << 0.015 * positives << "\n";
    std::cout << treeGen + insert + 0.015 * positives << "\n";
    std::cout << gx << " " << gy << "\n";
    std::cout << "-------------------\n";*/
    unique = 0;
    funique = 0;
    checkedbuckets.push_back(accepts);
    return treeGen + insert + 0.015 * positives;
}

std::string ionToString(std::vector<std::pair<int, std::string>> strInfo) {
    std::string str = "";
    for(int i = 0; i < 2000; i++) {
        str += "3";
    }
    for(int j = 0; j < strInfo.size(); j++) {
        str.replace(strInfo.at(j).first, 1, strInfo.at(j).second);
    }
    return str;
}

void realData(probs *prob, int N, std::vector<double> info) {
    std::string line;
    std::ifstream data ("/Users/sean/Documents/sample_spectra_10K_4rank.mgf");
    int id = 0;
    bool isX = true;
    std::vector<istring*> xstrings, ystrings;
    if(data.is_open()) {
        while(getline (data, line) ) {
            if(line == "BEGIN IONS") {
                std::vector<std::pair<int, std::string>> strInfo;
                while(getline (data, line)) {
                    if(line == "END IONS") {
                        std::string entry = ionToString(strInfo);
                        if(isX) {
                            istring *x = new istring();
                            x -> str = entry;
                            x -> id = id;
                            xstrings.push_back(x);
                            isX = false;
                        }
                        else {
                            if(matchProbability(prob, xstrings.at(id) -> str, entry)/matchProbabilityQ(prob, xstrings.at(id) -> str, entry) < 914.329) {
                                xstrings.pop_back();
                             }
                            else {
                                istring *y = new istring();
                                y -> str = entry;
                                y -> id = id;
                                ystrings.push_back(y);
                                id++;
                            }
                            isX = true;
                        }
                        break;
                    }
                    //std::cout << line.substr(0, line.find("\t")) << "\n";
                    //std::cout << stod(line.substr(0, line.find("\t"))) + 0.5 << "\n";
                    int index = (int)(stod(line.substr(0, line.find("\t"))) + 0.5);
                    std::string value = (line.substr(line.length() - 1, 1));
                    strInfo.push_back(std::make_pair(index, value));
                }
                
            }
        }
        data.close();
        std::cout << "Done reading data\n";
        std::cout << xstrings.size() << "\n";
        /*std::vector<double> vals;
        for(int i = 0; i < xstrings.size(); i++) {
            vals.push_back(matchProbability(prob, xstrings.at(i) -> str, ystrings.at(i) -> str)/matchProbabilityQ(prob, xstrings.at(i) -> str, ystrings.at(i) -> str));
        }
        std::sort(vals.begin(), vals.end());
        std::cout << vals.at(2000) << "\n";*/
        double besttime = 50000.0;
        double bp1 = 0.0;
        double bp2 = 0.0;
        double bp3 = 0.0;
        double bp4 = 0.0;
        
        //std::vector<double> as = {1.0/4, 1.0/2, 1, 2, 4, 8, 16};
        std::vector<double> mults = {1.0/16, 1.0/8, 1.0/4, 1.0/2, 1, 2, 4, 8, 16};
        std::vector<double> as = {4, 8, 16};
        
        for(int i = 0; i < 5000; i++) {
            found.push_back(false);
        }
        for(int i = 0; i < xstrings.size(); i++) {
            std::vector<bool> row(xstrings.size(), false);
            ufp.push_back(row);
        }
        
        std::cout << "C1,C2,C3,C4,Total complexity,Complexity of hashing,Complexity of tree construction,Complexity of FP,#bands needed,Overall empirical TP,TP at one band,Theoretical TP in one band,Overall empirical FP, empirical FP in one band,Theoretical FP in one band,#buckets,#nodes in the tree,total number of x nodes,total number of y nodes,gamma_x,gamma_y,mapped_x,mapped_y\n";
        
        for(int a = 0; a < 9; a++) {
            for(int b = 0; b < 9; b++) {
                    for(int d = 0; d < 9; d++) {
                        p1 = mults.at(a);
                        //p1 = as.at(a);
                        p2 = mults.at(b);
                        p3 = mults.at(b);
                        p4 = mults.at(d);
                        
                        
                        //if(a != 0.25 && b != 0.25 && d != 0.25) {
                            double candidate = runTest(prob, N, xstrings, ystrings, info, besttime);
                            if(candidate < besttime && candidate > 0) {
                                besttime = candidate;
                                bp1 = p1;
                                bp2 = p2;
                                bp3 = p3;
                                bp4 = p4;
                            }
                        //}
                    }
            }
            /*std::cout << "Best time: " << besttime << "\n";
            std::cout << "With params: " << bp1 << " " << bp2 << " " << bp3 << " " << bp4 << "\n";*/
        }
        std::cout << "Final Best time: " << besttime << "\n";
        std::cout << "With params: " << bp1 << " " << bp2 << " " << bp3 << " " << bp4 << "\n";
        /*p1 = 0.25;
        p2 = 0.5;
        p3 = 0.125;
        p4 = 0.125;
        for(int i = 0; i < 5000; i++) {
            found.push_back(false);
        }
        for(int i = 0; i < xstrings.size(); i++) {
            std::vector<bool> row(xstrings.size(), false);
            ufp.push_back(row);
        }
        runTest(prob, N, xstrings, ystrings, info);
        std::cout << unique << "\n";*/
    }
    else {
        std::cout << "Unable to open file\n";
    }
}

int main(int argc, const char * argv[]) {
    //Seed the RNG
    srand((unsigned int)time(NULL));
    
    /*if(argc < 4) {
        std::cout << "usage: ./simulated <N> <S> <i>\n";
        return -1;
    }
    
    int N = atoi(argv[1]);
    int S = atoi(argv[2]);
    double i = atof(argv[3]);*/
    
    int N = 2000;
    int S = 10000;
    
    double i = 0;
    
    double besttime = 500000000.0;
    double bp1 = 0.0;
    double bp2 = 0.0;
    double bp3 = 0.0;
    double bp4 = 0.0;
    
    std::vector<double> mults = {1.0/16.0, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16};
    std::vector<double> as = {1};
    std::vector<double> bs = {0.5};
    std::vector<double> cs = {0.0625};
    std::vector<double> ds = {0.0625};
    
    //1.4382877544533643, 5.023511952707102, 5.023511952707102, 3.36634674522192
    
    double t = (double)i/19;
    double p00 = 0.345 * (1 - t) + 0.019625 * t;
    double p01 = 0;
    double p10 = 0.31 * (1 - t) + 0.036875 * t;
    double p11 = 0.345 * (1 - t) + 0.9435 * t;
    arr pmatrix;
    pmatrix.push_back({p00, p01});
    pmatrix.push_back({p10, p11});
    probs *prob = set_pq_values(pmatrix);
    //std::vector<double> info = {2.341076852, 6.247385441, 1.85346908, 1.239621294};
    std::vector<double> info = {5.023511952707102, 5.023511952707102, 3.36634674522192, 1.4382877544533643};
    std::cout << info.at(0) << " " << info.at(1) << " " << info.at(2) << " " << info.at(3) << "\n";
    std::pair<std::vector<istring*>, std::vector<istring*>> strs = make_data(prob, N, S);
    std::vector<istring*> xs = strs.first;
    std::vector<istring*> ys = strs.second;
    
    for(int j = 0; j < N; j++) {
        found.push_back(false);
    }
    for(int j = 0; j < N; j++) {
        std::vector<bool> row(N, false);
        ufp.push_back(row);
    }
    
    std::cout << "C1,C2,C3,C4,Total complexity,Complexity of hashing,Complexity of tree construction,Complexity of FP,#bands needed,Overall empirical TP,TP at one band,Theoretical TP in one band,Overall empirical FP, empirical FP in one band,Theoretical FP in one band,#buckets,#nodes in the tree,total number of x nodes,total number of y nodes,gamma_x,gamma_y,mapped_x,mapped_y\n";
    
    for(int a = 0; a < mults.size(); a++) {
        for(int b = 0; b < mults.size(); b++) {
            for(int c = 0; c < mults.size(); c++) {
                for(int d = 0; d < mults.size(); d++) {
                    p1 = mults.at(a);
                    p2 = mults.at(b);
                    p3 = mults.at(c);
                    p4 = mults.at(d);
                    //std::cout << info.at(0) << " " << info.at(1) << " " << info.at(2) << " " << info.at(3) << "\n";
                    double candidate = runTest(prob, N, xs, ys, info, besttime);
                    if(candidate < besttime && candidate > 0) {
                        besttime = candidate;
                        bp1 = p1;
                        bp2 = p2;
                        bp3 = p3;
                        bp4 = p4;
                    }
                }
            }
        }
        //std::cout << "Best time: " << besttime << "\n";
        //std::cout << "With params: " << bp1 << " " << bp2 << " " << bp3 << " " << bp4 << "\n";
    }
    std::cout << "Final Best time: " << besttime << "\n";
    std::cout << "With params: " << bp1 << " " << bp2 << " " << bp3 << " " << bp4 << "\n";

    
    
    
    /*int N = 5000;
    
    arr pmatrix;
    pmatrix.push_back({0.000125, 0.00005, 0.000000097, 0.000405});
    pmatrix.push_back({0.00005, 0.00021, 0.0000062, 0.002});
    pmatrix.push_back({0.000000097, 0.0000062, 0.000027, 0.000355});
    pmatrix.push_back({0.000405, 0.002, 0.000355, 0.994165416});
    
    probs *prob = set_pq_values(pmatrix);
    
    //std::vector<double> info  = newLambda(prob);
    //std::vector<double> info = {1.02, 1.02, 0.76, 1.33333};
    std::vector<double> info = {1.0390874586281873, 1.1226949489947105, 0.7779293468903407, 1.3263327592227026};
    
    std::cout << info.at(0) << " " << info.at(1) << " " << info.at(2) << " " << info.at(3) << "\n";
    
    realData(prob, N, info);*/
    
    
    return 0;
}
