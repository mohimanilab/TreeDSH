Treehash=[7115,5671,5024,4405,3787,3128,2809,2339,2054,1800,1647,1607,1297,1257,1213,1075,1192,1181,1167,1181]
List=[]
Minhash=[11376,9964,7112,4898,4368,3362,2631,2581,1869,1503,1194,1161,1040,1106,768,640,631,398,404,589]
LSH=[7279,5994,4906,4003,3489,2924,2797,2319,2425,2281,2017,2015,1814,1992,2051,2287,2352,3099,4284,9615]
for x in range(0, 20):
    List.append(x/19)
    
import matplotlib.pyplot as plt
plt.plot(List, LSH, 'r-') # plotting t, a separately 

plt.plot(List, Minhash, 'b-') # plotting t, b separately 

plt.plot(List, Treehash, 'g-') # plotting t, c separately 

plt.xlabel('P(x)')
plt.ylabel('Simulation time (ms)')
plt.legend(['LSH','Minhash','TreeDSH'])