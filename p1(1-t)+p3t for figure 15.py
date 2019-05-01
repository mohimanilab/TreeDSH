
def minmax(a, n): 
  
    # inbuilt function to find the position of minimum  
    minpos = a.index(min(a)) 
      
    # inbuilt function to find the position of maximum  
    maxpos = a.index(max(a))  
      
    # printing the position  
    print ("The maximum is at position", maxpos + 1)  
    print ("The minimum is at position", minpos + 1)
      

import numpy as np
List=[]
Listtreehash=[]
Listminhash=[]
ListLSH=[]
Listperc=[]
Listdist1=[]
Listdist2=[]
Listdist3=[]
Listdist4=[]
Listdist5=[]
Listp00=[]
Listp01=[]
Listp10=[]
Listp11=[]
for x in range(0, 1):
    print(x)
    
#p3 0.0025,0,0.0225,0.975
#p2 0.345,0.0.31,0.345
#p1 0.215,0.0025,0.255,0.5275
#maxx stands for the number of quadruples 
maxx= int(input("Enter an integer number 'maxx' as the resolution for each of p00,p01,p10,p11. (some number around 200).  maxx=? "))

percm=0
for x in range(0, maxx+1):
            p00=0.215+(x/maxx)*(0.019625-0.215)
            p01=0.0025+(x/maxx)*(0-0.0025)
            p10=0.255+(x/maxx)*(0.036875-0.255)
            p11=0.5275+(x/maxx)*(0.9435-0.5275)
            q00=(p00+p01)*(p00+p10)
            q01=(p00+p01)*(p01+p11)
            q10=(p10+p11)*(p00+p10)
            q11=(p10+p11)*(p01+p11)
            if abs(p00*p11-p10*p01)<0.00001:
                #it is proved that in this case treehash=LSH=minhash=2
                treehash=-1+2
                LSH=-1+2
                minhash=-1+2
                print(p00,p01,p10,p11,'(p00p11=p01p10),minhash=1') 
                dist1=(LSH-minhash)
                dist2=-(LSH-minhash)
                dist3=min(LSH,minhash)-treehash
                dist4=minhash-treehash
                dist5=LSH-treehash
                Listtreehash.append(treehash)
                Listminhash.append(minhash)
                ListLSH.append(LSH)
                List.append(x/maxx)
                Listdist1.append(dist1)
                Listdist2.append(dist2)
                Listdist3.append(dist3)
                Listdist4.append(dist4)
                Listdist5.append(dist5)
                Listp00.append(p00)
                Listp01.append(p01)
                Listp10.append(p10)
                Listp11.append(p11)
                #print("perc=",(1-treehash/min(LSH,minhash)),"treehash=",treehash,"LSH=",LSH,"minhash=",minhash,p00,p01,p10,p11,abs(p00*p11-p10*p01),x,y,z)
           #if p00*p11>p10*p01, q00*(p00/q00)**mu+q11*(p11/q11)**mu=1 should be solved
           #if p00*p11<p10*p01, q01*(p01/q01)**mu+q10*(p10/q10)**mu=1 should be solved
            elif p00*p11>p10*p01:
                mu=1
                dmu=1
                while mu>=1:
                    mu=mu-1*(q00*(p00/q00)**mu+q11*(p11/q11)**mu-1)/(np.log(p00/q00)*q00*(p00/q00)**mu+np.log(p11/q11)*q11*(p11/q11)**mu)
                    dmu=abs(1*(q00*(p00/q00)**mu+q11*(p11/q11)**mu-1)/(np.log(p00/q00)*q00*(p00/q00)**mu+np.log(p11/q11)*q11*(p11/q11)**mu))
                    if dmu<0.00001:
                        break
                #print(mu)
                treehash=-1+1+(mu-1)/mu
                LSH=-1+min(1+np.log(p01+p10)/np.log(q10+q01),1+np.log(p00+p11)/np.log(q00+q11))
                minhash=-1+min(1+np.log(p00/(1-p11))/np.log(q00/(1-q11)),1+np.log(p11/(1-p00))/np.log(q11/(1-q00)),1+np.log(p01/(1-p10))/np.log(q01/(1-q10)),1+np.log(p10/(1-p01))/np.log(q10/(1-q01)))
                if 1+np.log(p00/(1-p11))/np.log(q00/(1-q11))<=min(1+np.log(p11/(1-p00))/np.log(q11/(1-q00)),1+np.log(p01/(1-p10))/np.log(q01/(1-q10)),1+np.log(p10/(1-p01))/np.log(q10/(1-q01))):
                    print('00/11')
                print(p00,p01,p10,p11,'mihash=log(p00/(1-p11))/log(q00/(1-q11))=',np.log(p00/(1-p11))/np.log(q00/(1-q11)))
                if 1+np.log(p01/(1-p10))/np.log(q01/(1-q10))<=min(1+np.log(p00/(1-p11))/np.log(q00/(1-q11)),1+np.log(p11/(1-p00))/np.log(q11/(1-q00)),1+np.log(p10/(1-p01))/np.log(q10/(1-q01))):
                    print('01/10')
                print(p00,p01,p10,p11,'mihash=log(p01/(1-p10))/log(q01/(1-q10))=',np.log(p01/(1-p10))/np.log(q01/(1-q10)))
                if 1+np.log(p11/(1-p00))/np.log(q11/(1-q00))<=min(1+np.log(p00/(1-p11))/np.log(q00/(1-q11)),1+np.log(p01/(1-p10))/np.log(q01/(1-q10)),1+np.log(p10/(1-p01))/np.log(q10/(1-q01))):
                    print('11/00')
                print(p00,p01,p10,p11,'minhash=log(p11/(1-p00))/log(q11/(1-q00))=',np.log(p11/(1-p00))/np.log(q11/(1-q00)))
                if 1+np.log(p10/(1-p01))/np.log(q10/(1-q01))<=min(1+np.log(p00/(1-p11))/np.log(q00/(1-q11)),1+np.log(p11/(1-p00))/np.log(q11/(1-q00)),1+np.log(p01/(1-p10))/np.log(q01/(1-q10))):
                    print('10/01')
                print(p00,p01,p10,p11,'minhash=log(p10/(1-p01))/log(q10/(1-q01))=',np.log(p10/(1-p01))/np.log(q10/(1-q01))) 
                if 1+np.log(p01+p10)/np.log(q10+q01)<1+np.log(p00+p11)/np.log(q00+q11):
                    print('01+10')
                if 1+np.log(p01+p10)/np.log(q10+q01)>=1+np.log(p00+p11)/np.log(q00+q11):
                    print('00+11')
                print(p00,p01,p10,p11,'LSH=log(p01+p10)/log(q10+q01)=',np.log(p01+p10)/np.log(q10+q01))
                print(p00,p01,p10,p11,'LSH=np.log(p00+p11)/np.log(q00+q11)=',np.log(p00+p11)/np.log(q00+q11))
                dist1=(LSH-minhash)
                dist2=-(LSH-minhash)
                dist3=min(LSH,minhash)-treehash
                dist4=minhash-treehash
                dist5=LSH-treehash
                Listtreehash.append(treehash)
                Listminhash.append(minhash)
                ListLSH.append(LSH)
                List.append(x/maxx)
                Listdist1.append(dist1)
                Listdist2.append(dist2)
                Listdist3.append(dist3)
                Listdist4.append(dist4)
                Listdist5.append(dist5)
                Listp00.append(p00)
                Listp01.append(p01)
                Listp10.append(p10)
                Listp11.append(p11)
               # print("perc=",perc,"treehash=",treehash,"LSH=",LSH,"minhash=",minhash,p00,p01,p10,p11,abs(p00*p11-p10*p01),x,y,z) 
            else:
                mu=1
                dmu=1
                while mu>=1:
                    mu=mu-(q01*(p01/q01)**mu+q10*(p10/q10)**mu-1)/(np.log(p01/q01)*q01*(p01/q01)**mu+np.log(p10/q10)*q10*(p10/q10)**mu)
                    dmu=abs((q01*(p01/q01)**mu+q10*(p10/q10)**mu-1)/(np.log(p01/q01)*q01*(p01/q01)**mu+np.log(p10/q10)*q10*(p10/q10)**mu))
                    if dmu<0.00001:
                        break
                #print(mu)
                treehash=-1+1+(mu-1)/mu
                LSH=-1+min(1+np.log(p01+p10)/np.log(q10+q01),1+np.log(p00+p11)/np.log(q00+q11))
                minhash=-1+min(1+np.log(p00/(1-p11))/np.log(q00/(1-q11)),1+np.log(p11/(1-p00))/np.log(q11/(1-q00)),1+np.log(p01/(1-p10))/np.log(q01/(1-q10)),1+np.log(p10/(1-p01))/np.log(q10/(1-q01)))
                if 1+np.log(p00/(1-p11))/np.log(q00/(1-q11))<=min(1+np.log(p11/(1-p00))/np.log(q11/(1-q00)),1+np.log(p01/(1-p10))/np.log(q01/(1-q10)),1+np.log(p10/(1-p01))/np.log(q10/(1-q01))):
                    print('00/11')
                print(p00,p01,p10,p11,'mihash=log(p00/(1-p11))/log(q00/(1-q11))=',np.log(p00/(1-p11))/np.log(q00/(1-q11)))
                if 1+np.log(p01/(1-p10))/np.log(q01/(1-q10))<=min(1+np.log(p00/(1-p11))/np.log(q00/(1-q11)),1+np.log(p11/(1-p00))/np.log(q11/(1-q00)),1+np.log(p10/(1-p01))/np.log(q10/(1-q01))):
                    print('01/10')
                print(p00,p01,p10,p11,'mihash=log(p01/(1-p10))/log(q01/(1-q10))=',np.log(p01/(1-p10))/np.log(q01/(1-q10)))
                if 1+np.log(p11/(1-p00))/np.log(q11/(1-q00))<=min(1+np.log(p00/(1-p11))/np.log(q00/(1-q11)),1+np.log(p01/(1-p10))/np.log(q01/(1-q10)),1+np.log(p10/(1-p01))/np.log(q10/(1-q01))):
                    print('11/00')
                print(p00,p01,p10,p11,'minhash=log(p11/(1-p00))/log(q11/(1-q00))=',np.log(p11/(1-p00))/np.log(q11/(1-q00)))
                if 1+np.log(p10/(1-p01))/np.log(q10/(1-q01))<=min(1+np.log(p00/(1-p11))/np.log(q00/(1-q11)),1+np.log(p11/(1-p00))/np.log(q11/(1-q00)),1+np.log(p01/(1-p10))/np.log(q01/(1-q10))):
                    print('10/01')
                print(p00,p01,p10,p11,'minhash=log(p10/(1-p01))/log(q10/(1-q01))=',np.log(p10/(1-p01))/np.log(q10/(1-q01))) 
                if 1+np.log(p01+p10)/np.log(q10+q01)<1+np.log(p00+p11)/np.log(q00+q11):
                    print('01+10')
                if 1+np.log(p01+p10)/np.log(q10+q01)>=1+np.log(p00+p11)/np.log(q00+q11):
                    print('00+11')
                print(p00,p01,p10,p11,'LSH=log(p01+p10)/log(q10+q01)=',np.log(p01+p10)/np.log(q10+q01))
                print(p00,p01,p10,p11,'LSH=np.log(p00+p11)/np.log(q00+q11)=',np.log(p00+p11)/np.log(q00+q11))
                dist1=(LSH-minhash)
                dist2=-(LSH-minhash)
                dist3=min(LSH,minhash)-treehash
                dist4=minhash-treehash
                dist5=LSH-treehash
                #percc=1-(treehash)/min(LSH,minhash)
                #if percc<0:
                 #   perc=-(-percc//0.00001)*0.00001
               # else:
               # perc=(percc//0.00001)*0.00001
                Listtreehash.append(treehash)
                Listminhash.append(minhash)
                ListLSH.append(LSH)
                List.append(x/maxx)
                Listdist1.append(dist1)
                Listdist2.append(dist2)
                Listdist3.append(dist3)
                Listdist4.append(dist4)
                Listdist5.append(dist5)
                Listp00.append(p00)
                Listp01.append(p01)
                Listp10.append(p10)
                Listp11.append(p11)


import matplotlib.pyplot as plt
plt.plot(List, ListLSH, 'r-') # plotting t, a separately 

plt.plot(List, Listminhash, 'b-') # plotting t, b separately 

plt.plot(List, Listtreehash, 'g-') # plotting t, c separately 

plt.xlabel('P(x)')
plt.ylabel('Theoretical guarantees')
plt.legend(['LSH','Minhash','Treehash'])










