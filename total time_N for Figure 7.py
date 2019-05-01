LogN=[2.698970004,3,3.301029996,3.698970004,4,4.301029996]
LogTreehash=[2.365487985,2.423245874,2.514547753,2.752816431,3.255272505,3.777934049]

LogMinhash=[1.763427994,1.991226076,2.328379603,3.027757205,3.70918513,4.331062701]
LogLSH=[1.447158031,1.851258349,2.307496038,3.172018809,3.855276804,4.492564843]



import matplotlib.pyplot as plt
plt.plot(LogN, LogLSH, 'r-') # plotting t, a separately 

plt.plot(LogN, LogMinhash, 'b-') # plotting t, b separately 

plt.plot(LogN, LogTreehash, 'g-') # plotting t, c separately 

plt.xlabel('$\log_{10}(N)$')
plt.ylabel('$\log_{10}$(Total Time)')
plt.legend(['LSH','Minhash','TreeDSH'])
#