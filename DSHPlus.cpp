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

//Keeping track of the number of accepts, rejects, alpha, and beta
int accepts = 0;
int rejects = 0;
double TFP = 0.0;
double TTP = 0.0;

//Debugging purposes
int numNodes = 0;
int xnodes = 0;
int ynodes = 0;

//2D vector that holds doubles
typedef std::vector<std::vector<double>> arr;

//p and q are 2D vectors that represent a probability matrix
struct probs_header {
    arr p, q;
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
    
    for(int i = 0; i < pmatrix.size(); i++) {
        std::vector<double> newrow;
        for(int j = 0; j < pmatrix.size(); j++) {
            newrow.push_back(rowsums.at(i) * colsums.at(j));
        }
        qmatrix.push_back(newrow);
    }
    
    ret -> q = qmatrix;
    
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
            xstr += std::to_string(newX);
            ystr += std::to_string(newY);
        }
        x -> str = xstr;
        y -> str = ystr;
        xs.push_back(x);
        ys.push_back(y);
    }
    return std::make_pair(xs, ys);
}

//Should we accept as a bucket?
bool accept(int N, double C, double P, double Q) {
    return P / Q > (double) (N * N) / C;
}

//Should we reject as a bucket?
bool reject(double P, double Q, double C, double l1, double l2) {
    return pow(P, l1) * pow(Q, l2) < 1 / (double) C;
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
        for(int i = 0; i < (prob -> p).size(); i++) {
            for(int j = 0; j < (prob -> p).size(); j++) {
                slope += pow((prob -> p).at(i).at(j), mu + 1) * pow((prob -> p).at(i).at(j), t) * log((prob -> q).at(i).at(j));
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
std::pair<double, double> calculateLambdas(probs* prob) {
    double R = 1.0;
    double minval = 100.0;
    double lam1 = 1.0;
    double lam2 = 1.0;
    //Go from 0 to 99 since approxT takes mu and solves p^(mu + 1) * q^t = 1
    for(double i = 0.0; i <= 99.0; i += 0.1) {
        double l2 = approxT(prob, i);
        arr rijs;
        for(int k = 0; k < (prob -> p).size(); k++) {
            std::vector<double> newrow;
            for(int j = 0; j < (prob -> p).size(); j++) {
                newrow.push_back(computerij(R, (prob -> p).at(k).at(j), (prob -> q).at(k).at(j), i + 1, l2));
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
        double candidate = (std::max(val1 / (val3 - val1), val1 / (2 * val1 - val2)));
        
        //Is it lower?
        if(val3 < val1 && 2 * val1 < val2 && candidate < minval) {
            lam1 = i + 1;
            lam2 = l2;
            minval = candidate;
        }
    }
    return std::make_pair(lam1, lam2);
}

//Recursive call for generating a tree
void generateTreeHelper(probs* prob, double oldP, double oldQ, int xind, int yind, int N, double C, double l1, double l2, tnode *xtree, tnode *ytree) {
    //Debugging purposes
    numNodes++;
    
    //Update the probabilities
    double newP = oldP * (prob -> p).at(xind).at(yind);
    double newQ = oldQ * (prob -> q).at(xind).at(yind);
    
    /*NEW CONTROL FLOW
      First create the xtree and ytree nodes if they haven't already been created.
      Check for accept/reject
      Recursive call if not
     */
    //Create new nodes if they don't already exist
    if((xtree -> children).at(xind) == NULL) {
        tnode *newX = new tnode();
        newX -> str = (xtree -> str) + std::to_string(xind);
        for(int i = 0; i < (prob -> p).size(); i++) {
            (newX -> children).push_back(NULL);
        }
        (xtree -> children).at(xind) = newX;
        xnodes++;
    }
    if((ytree -> children).at(yind) == NULL) {
        tnode *newY = new tnode();
        newY -> str = (ytree -> str) + std::to_string(yind);
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
    if(accept(N, C, newP, newQ)) {
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
        return;
    }
    //Check for rejection
    else if(reject(newP, newQ, C, l1, l2)) {
        //Update number of rejections
        rejects++;
        return;
    }
    else {
        //Recursive calls.
        for(int i = 0; i < (prob -> p).size(); i++) {
            for(int j = 0; j < (prob -> p).size(); j++) {
                generateTreeHelper(prob, newP, newQ, i, j, N, C, l1, l2, currX, currY);
            }
        }
        return;
    }
}

//Create the root and recursively generate the tree
std::pair<tnode*, tnode*> generateTree(probs* prob, int N, double C, double l1, double l2) {
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
            generateTreeHelper(prob, 1.0, 1.0, i, j, N, C, l1, l2, xtree, ytree);
        }
    }
    return std::make_pair(xtree, ytree);
}

//Feed a string into a tree and put it in every bucket along the way to the leaves.
void feedString(istring *istr, std::string str, tnode *tree, bool isX) {
    if(str.length() == 0) {
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
    feedString(istr, str.substr(1, str.length() - 1), (tree -> children).at(nextInd), isX);
}

//Feed every string into a tree
void feedStrings(std::vector<istring*> xstrings, std::vector<istring*> ystrings, tnode *xtree, tnode *ytree) {
    for(auto i = xstrings.begin(); i != xstrings.end(); ++i) {
        feedString(*i, (*i) -> str, xtree, true);
    }
    for(auto i = ystrings.begin(); i != ystrings.end(); ++i) {
        feedString(*i, (*i) -> str, ytree, false);
    }
}

//Analyze the number of true and false positives in a single bucket
std::pair<int, int> analyzeBucket(bucket *b) {
    std::vector<istring*> xs = b -> xmatches;
    std::vector<istring*> ys = b -> ymatches;
    int TP = 0;
    int FP = 0;
    
    for(auto i = xs.begin(); i != xs.end(); ++i) {
        istring *x = *i;
        for(auto j = ys.begin(); j != ys.end(); ++j) {
            istring *y = *j;
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
std::pair<int, int> analyzeTree(tnode *tree) {
    if(tree == NULL) {
        return std::make_pair(0, 0);
    }
    int TP = 0;
    int FP = 0;
    std::vector<bucket*> buckets = tree -> bvec;
    for(auto i = buckets.begin(); i != buckets.end(); ++i) {
        std::pair<int, int> binfo = analyzeBucket(*i);
        TP += binfo.first;
        FP += binfo.second;
    }
    for(int i = 0; i < (tree -> children).size(); i++) {
        std::pair<int, int> info = analyzeTree((tree -> children).at(i));
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
    for(int i = x.length() - 1; i > 0; --i) {
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

int main(int argc, const char * argv[]) {
    //Seed the RNG
    srand((unsigned int)time(NULL));
    
    //Generate the probability distribution
    //probs *prob = generateProbability(2);
    arr pmatrix;
    pmatrix.push_back({0.215, 0.0025});
    pmatrix.push_back({0.255, 0.5275});
    probs *prob = set_pq_values(pmatrix);
    
    //Calculate lambda1 and lambda2
    std::pair<double, double> lambdas = calculateLambdas(prob);
    double l1 = lambdas.first;
    double l2 = lambdas.second;
    
    //Check if lambdas are too high. May change calculateLambdas to return lambda* as well and case on that.
    if(l1 > 40) {
        std::cout << "Lambdas too high\n";
        return 0;
    }
    
    //Generate the x and y trees.
    std::pair<tnode*, tnode*> trees = generateTree(prob, 1000, 200000, l1, l2);
    
    std::cout<<numNodes<<" "<<accepts<<" "<<xnodes<<" "<<ynodes<<"\n";
    
    //Stop immediately if the distribution made a really bad tree (no accept nodes)
    if(accepts == 0) {
        std::cout << "Bad Tree\n";
        return 0;
    }
    std::cout << "Tree generated!\n";
    tnode *xtree = trees.first;
    tnode *ytree = trees.second;
    
    //Generate the random pairs
    std::pair<std::vector<istring*>, std::vector<istring*>> strs = generateStrings(prob, 1000, 100);
    std::vector<istring*> xs = strs.first;
    std::vector<istring*> ys = strs.second;
    
    //Initial values for the TP rate, the number of TP and FP, and b
    double tprate = 1.0;
    unsigned long EFP = 0;
    unsigned long ETP = 0;
    int b = 0;
    
    //Run until the TP rate is 90%
    while(tprate > 0.1) {
        feedStrings(xs, ys, xtree, ytree);
        std::pair<int, int> info = analyzeTree(xtree);
        int TP = info.first;
        int FP = info.second;
        //std::cout<<TP<<" "<<FP<<"\n";
        ETP += TP;
        EFP += FP;
        if(TP + FP == 0) {
            std::cout << "TP + FP = 0\n";
            return 0;
        }
        tprate *= (1 - ((double)TP/(TP + FP)));
        b++;
        if(tprate > 0.05) {
            permuteStrings(xs, ys);
            clearBuckets(xtree);
        }
    }
    std::cout << "Theoretical False Positive: " << TFP << "\n";
    std::cout << "Theoretical True Positive: " << TTP << "\n";
    std::cout << "Empirical False Positive: " << EFP << "\n";
    std::cout << "True Positive Rate: " << ((double)ETP/(ETP + EFP)) << "\n";
    std::cout << "Number of bands: " << b << "\n";
    return 0;
}
