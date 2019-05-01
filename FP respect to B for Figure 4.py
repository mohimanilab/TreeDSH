TP=0.9  
import numpy as np
Listx=[]
Listy=[]

p00=0.25
p01=0.215
p10=0.53
p11=0.005
q00=(p00+p01)*(p00+p10)
q01=(p11+p01)*(p00+p01)
q10=(p10+p11)*(p10+p00)
q11=(p11+p01)*(p11+p10)
#For different values of threshold(thr), we find b and FP as we know alpha and beta from the code "build the tree".
thr=400 #0.23 to 0.33
b=np.log(1-TP)/np.log(1-0.010534338309392361)
FP=1-(1-0.0000209)**b
Listx.append(np.log(b))
Listy.append(np.log(FP))
print("treehash,thr=",thr,"b=",b,"FP=",FP)

thr=200 #0.23 to 0.33
b=np.log(1-TP)/np.log(1-0.017962598443482712)
FP=1-(1-0.000074195785)**b
Listx.append(np.log(b))
Listy.append(np.log(FP))
print("treehash,thr=",thr,"b=",b,"FP=",FP)


thr=80 #0.23 to 0.33
b=np.log(1-TP)/np.log(1-0.03418904486382274)
FP=1-(1-0.00034452970224510885)**b
Listx.append(np.log(b))
Listy.append(np.log(FP))
print("treehash,thr=",thr,"b=",b,"FP=",FP)

thr=40 #0.23 to 0.33
b=np.log(1-TP)/np.log(1-0.052401106600888175)
FP=1-(1-0.000942033218426659)**b
Listx.append(np.log(b))
Listy.append(np.log(FP))
print("treehash,thr=",thr,"b=",b,"FP=",FP)

thr=20 #0.23 to 0.33
b=np.log(1-TP)/np.log(1-0.08919364472887166)
FP=1-(1-0.003327224609890983)**b
Listx.append(np.log(b))
Listy.append(np.log(FP))
print("treehash,thr=",thr,"b=",b,"FP=",FP)

thr=10 #0.23 to 0.33
b=np.log(1-TP)/np.log(1-0.15223524400810445)
FP=1-(1-0.011831826017281794)**b
Listx.append(np.log(b))
Listy.append(np.log(FP))
print("treehash,thr=",thr,"b=",b,"FP=",FP)

thr=5 #0.33 to 0.48
b=np.log(1-TP)/np.log(1-0.26024361269060503)
FP=1-(1-0.042230698888406804)**b
Listx.append(np.log(b))
Listy.append(np.log(FP))
print("treehash,thr=",thr,"b=",b,"FP=",FP)

thr=3 #0.68
b=np.log(1-TP)/np.log(1-0.36719663845000006)
FP=1-(1-0.09583354506620656)**b
Listx.append(np.log(b))
Listy.append(np.log(FP))
print("treehash,thr=",thr,"b=",b,"FP=",FP)

thr=1.5 #to 1
b=np.log(1-TP)/np.log(1-0.60985)
FP=1-(1-0.31912908000000006)**b
Listx.append(np.log(b))
Listy.append(np.log(FP))
print("treehash,thr=",thr,"b=",b,"FP=",FP)



 


p00=0.25
p01=0.215
p10=0.53
p11=0.005
q00=(p00+p01)*(p00+p10)
q01=(p11+p01)*(p00+p01)
q10=(p10+p11)*(p10+p00)
q11=(p11+p01)*(p11+p10)

minhash1=1+np.log(p00/(1-p11))/np.log(q00/(1-q11))
minhash2=1+np.log(p11/(1-p00))/np.log(q11/(1-q00))
minhash3=1+np.log(p01/(1-p10))/np.log(q01/(1-q10))
minhash4=1+np.log(p10/(1-p01))/np.log(q10/(1-q01))
print('minhash=',minhash1,minhash2,minhash3,minhash4)
Listxminhash=[]
Listyminhash=[]

xp=p01/(1-p10)
xq=q01/(1-q10)

#For different values of r, we find b and FP 
r=1
b=np.log(1-TP)/np.log(1-(xp)**r)
FP=1-(1-(xq)**r)**b           
Listxminhash.append(np.log(b))
Listyminhash.append(np.log(FP))
print("minhash,r=",r,"b=",b,"FP=",FP)

r=2
b=np.log(1-TP)/np.log(1-(xp)**r)
FP=1-(1-(xq)**r)**b           
Listxminhash.append(np.log(b))
Listyminhash.append(np.log(FP))
print("minhash,r=",r,"b=",b,"FP=",FP)

r=3
b=np.log(1-TP)/np.log(1-(xp)**r)
FP=1-(1-(xq)**r)**b          
Listxminhash.append(np.log(b))
Listyminhash.append(np.log(FP))
print("minhash,r=",r,"b=",b,"FP=",FP)

r=4
b=np.log(1-TP)/np.log(1-(xp)**r)
FP=1-(1-(xq)**r)**b          
Listxminhash.append(np.log(b))
Listyminhash.append(np.log(FP))
print("minhash,r=",r,"b=",b,"FP=",FP)


r=5
b=np.log(1-TP)/np.log(1-(xp)**r)
FP=1-(1-(xq)**r)**b          
Listxminhash.append(np.log(b))
Listyminhash.append(np.log(FP))
print("minhash,r=",r,"b=",b,"FP=",FP)

r=6
b=np.log(1-TP)/np.log(1-(xp)**r)
FP=1-(1-(xq)**r)**b          
Listxminhash.append(np.log(b))
Listyminhash.append(np.log(FP))
print("minhash,r=",r,"b=",b,"FP=",FP)


ListxLSH=[]
ListyLSH=[]



#For different values of r, we find b and FP 
xp=p01+p10
xq=q01+q10
r=1
b=np.log(1-TP)/np.log(1-(xp)**r)
FP=1-(1-(xq)**r)**b           
ListxLSH.append(np.log(b))
ListyLSH.append(np.log(FP))
print("LSH,r=",r,"b=",b,"FP=",FP)

r=2
b=np.log(1-TP)/np.log(1-(xp)**r)
FP=1-(1-(xq)**r)**b           
ListxLSH.append(np.log(b))
ListyLSH.append(np.log(FP))
print("LSH,r=",r,"b=",b,"FP=",FP)

r=3
b=np.log(1-TP)/np.log(1-(xp)**r)
FP=1-(1-(xq)**r)**b           
ListxLSH.append(np.log(b))
ListyLSH.append(np.log(FP))
print("LSH,r=",r,"b=",b,"FP=",FP)

r=4
b=np.log(1-TP)/np.log(1-(xp)**r)
FP=1-(1-(xq)**r)**b           
ListxLSH.append(np.log(b))
ListyLSH.append(np.log(FP))
print("LSH,r=",r,"b=",b,"FP=",FP)


r=5
b=np.log(1-TP)/np.log(1-(xp)**r)
FP=1-(1-(xq)**r)**b           
ListxLSH.append(np.log(b))
ListyLSH.append(np.log(FP))
print("LSH,r=",r,"b=",b,"FP=",FP)


r=6
b=np.log(1-TP)/np.log(1-(xp)**r)
FP=1-(1-(xq)**r)**b           
ListxLSH.append(np.log(b))
ListyLSH.append(np.log(FP))
print("LSH,r=",r,"b=",b,"FP=",FP)

r=7
b=np.log(1-TP)/np.log(1-(xp)**r)
FP=1-(1-(xq)**r)**b           
ListxLSH.append(np.log(b))
ListyLSH.append(np.log(FP))
print("LSH,r=",r,"b=",b,"FP=",FP)

r=8
b=np.log(1-TP)/np.log(1-(xp)**r)
FP=1-(1-(xq)**r)**b           
ListxLSH.append(np.log(b))
ListyLSH.append(np.log(FP))
print("LSH,r=",r,"b=",b,"FP=",FP)

r=9
b=np.log(1-TP)/np.log(1-(xp)**r)
FP=1-(1-(xq)**r)**b           
ListxLSH.append(np.log(b))
ListyLSH.append(np.log(FP))
print("LSH,r=",r,"b=",b,"FP=",FP)


r=10
b=np.log(1-TP)/np.log(1-(xp)**r)
FP=1-(1-(xq)**r)**b           
ListxLSH.append(np.log(b))
ListyLSH.append(np.log(FP))
print("LSH,r=",r,"b=",b,"FP=",FP)


r=11
b=np.log(1-TP)/np.log(1-(xp)**r)
FP=1-(1-(xq)**r)**b           
ListxLSH.append(np.log(b))
ListyLSH.append(np.log(FP))
print("LSH,r=",r,"b=",b,"FP=",FP)


r=12
b=np.log(1-TP)/np.log(1-(xp)**r)
FP=1-(1-(xq)**r)**b           
ListxLSH.append(np.log(b))
ListyLSH.append(np.log(FP))
print("LSH,r=",r,"b=",b,"FP=",FP)

r=13
b=np.log(1-TP)/np.log(1-(xp)**r)
FP=1-(1-(xq)**r)**b           
ListxLSH.append(np.log(b))
ListyLSH.append(np.log(FP))
print("LSH,r=",r,"b=",b,"FP=",FP)


r=14
b=np.log(1-TP)/np.log(1-(xp)**r)
FP=1-(1-(xq)**r)**b           
ListxLSH.append(np.log(b))
ListyLSH.append(np.log(FP))
print("LSH,r=",r,"b=",b,"FP=",FP)



r=16
b=np.log(1-TP)/np.log(1-(xp)**r)
FP=1-(1-(xq)**r)**b           
ListxLSH.append(np.log(b))
ListyLSH.append(np.log(FP))
print("LSH,r=",r,"b=",b,"FP=",FP)









import matplotlib.pyplot as plt
plt.plot(Listx, Listy, 'ro-') # plotting t, a separately 
plt.legend(Listx,'af')
plt.plot(Listxminhash, Listyminhash, 'bo-') # plotting t, b separately 
plt.legend(Listxminhash,'f')
plt.plot(ListxLSH, ListyLSH, 'go-') # plotting t, c separately 
plt.legend(ListxLSH,'fs')
plt.xlabel('log(b)')
plt.ylabel('log(FP)')
plt.legend(['Treehash','Minhash','LSH'])
#plt.title('False positive is plotted respect to b in logarithmic scale')
plt.show()
LSH1=1+np.log(p01+p10)/np.log(q10+q01)
LSH2=1+np.log(p00+p11)/np.log(q00+q11)
print(LSH1,LSH2)
mu=1
dmu=1
while mu>=1:
    mu=mu-(q01*(p01/q01)**mu+q10*(p10/q10)**mu-1)/(np.log(p01/q01)*q01*(p01/q01)**mu+np.log(p10/q10)*q10*(p10/q10)**mu)
    dmu=abs((q01*(p01/q01)**mu+q10*(p10/q10)**mu-1)/(np.log(p01/q01)*q01*(p01/q01)**mu+np.log(p10/q10)*q10*(p10/q10)**mu))
    if dmu<0.00001:
        break
print(mu)
treehash=1+(mu-1)/mu
print('treehash=',treehash)

               

