Generating figure 6 using benchmark.sh: Need to provide 11 arguments:
1) N (size of sets)
2,3) Value of b and threshold for treehash
4,5) Value of b and r for minhash
6,7) Value of b and r for LSH
8,9,10,11) Values of p00, p01, p10, p11, must sum to 1.
Values of b, r, and threshold must be tuned according to the probability 
matrix such that TP rate is at least 0.99 (these hyperparameters are
independent of the value of N).

The results of this script are stored in test.out, and are printed on 3 lines,
the first one has results for treehash, the second one is for minhash, and the
third one is for LSH. Each line is of the following form:
1) N 
2) r/thresh
3) b 
4) expected fp rate 
5) observed fp rate 
6) observed tp rate 
7) time for hashing 
8) time for bucketing and finding matches
9) total time
10) GOOD TP if observed tp rate >= 0.99, blank otherwise.


Generating figure 7 using bigbenchmark.sh: Need to provide 4 arguments
1) b
2) threshold for treehash
3) r for minhash
4) r for LSH
The last 3 arguments must be tuned to obtain 0.99 accuracy for the fixed value
of b. Will run for N = 500, 1000, 2000, 5000, 10000, 20000, 50000, 1000000.


Generating figure 12 using treehash:
1) First change treehash.cpp. In the main function comment out 
simulated_data_steup() and uncomment msgf_data_setup(). Similarly comment out 
obtain_all_matches() and uncomment obtain_msgf_matches(). 
2) Remake treehash as in makeall.sh.
3) ./treehash takes in 3 arguments for msgf data: N, b, threshold. Ignore value
of N, and tune b to obtain the true positive for given value of threshold.
