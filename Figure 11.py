Listthresh=[]
Listnewthresh=[]
ListComplex=[]
ListConstant=[]
Listlenx=[]
Listleny=[]
Listthresh=[]

import numpy as np
p00=0.215
p01=0.0025
p10=0.255
p11=0.5275
q00=(p00+p01)*(p00+p10)
q01=(p11+p01)*(p00+p01)
q10=(p10+p11)*(p10+p00)
q11=(p11+p01)*(p11+p10)
pleftTree = p00/q00
prightTree = p11/q11
rmin=min(q00/p00,q11/p11)
print(rmin)
if p00*p11>p10*p01:
    mu=1
    dmu=1
    while mu>=1:
        mu=mu-1*(q00*(p00/q00)**mu+q11*(p11/q11)**mu-1)/(np.log(p00/q00)*q00*(p00/q00)**mu+np.log(p11/q11)*q11*(p11/q11)**mu)
        dmu=abs(1*(q00*(p00/q00)**mu+q11*(p11/q11)**mu-1)/(np.log(p00/q00)*q00*(p00/q00)**mu+np.log(p11/q11)*q11*(p11/q11)**mu))
        if dmu<0.00001:
            break
                #print(mu)
    treehashh=-1+1+(mu-1)/mu
    print('+',mu)    
else:
    mu=1
    dmu=1
    while mu>=1:
        mu=mu-(q01*(p01/q01)**mu+q10*(p10/q10)**mu-1)/(np.log(p01/q01)*q01*(p01/q01)**mu+np.log(p10/q10)*q10*(p10/q10)**mu)
        dmu=abs((q01*(p01/q01)**mu+q10*(p10/q10)**mu-1)/(np.log(p01/q01)*q01*(p01/q01)**mu+np.log(p10/q10)*q10*(p10/q10)**mu))
        if dmu<0.00001:
            break
                #print(mu)
    treehashh=-1+1+(mu-1)/mu
    print('-',mu)     
print('f',treehashh)      


#
#
#if pleftTree>1 and prightTree>1 :
m=300
for t in range(1, m):   
    thresh = 1+t/10+t**2/100
    newthresh=np.log(thresh)
#    thresh=1.5
    pleftTree = p00/q00
    prightTree = p11/q11


    pleftEval = p00+p01
    prightEval = p10+p11
        
    qleftEval = p00+p10
    qrightEval = p01+p11

    class Node:
        val = 0
        left = None
        right = None

        def __init__(self, v):
            self.val = v
            self.right = None
            self.left = None




    def buildTree(node):
        if node.val > thresh:
            return
    # if node.val*pleftTree > thresh:
        node.left = Node(node.val * pleftTree)
        buildTree(node.left)
   # if node.val*prightTree > thresh:
        node.right = Node(node.val * prightTree)
        buildTree(node.right)


    def eval(sum, node, pr, pl, val, depth):
        if node.left is None and node.right is None:
            return sum+val*(depth-1)
        sl = 0
        sr = 0
        if node.left is not None:
            sl = eval(sum, node.left, pr, pl, val * pl, depth+1)
        if node.right is not None:
            sr = eval(sum, node.right, pr, pl, val * pr, depth+1)
        return sl+sr

    if __name__ == "__main__":
        root = Node(1)
        buildTree(root)


        lenx=eval(0.0, root, prightEval, pleftEval, 1.0,1)
        leny=eval(0.0, root, qrightEval, qleftEval, 1.0,1)
        Listthresh.append(thresh)
        Listnewthresh.append(newthresh)
        Listlenx.append(lenx)
        Listleny.append(leny)
# 
    
    

#
import matplotlib.pyplot as plt
plt.plot(Listnewthresh, Listlenx, 'r-') 
plt.plot(Listnewthresh, Listleny, 'b-') 
plt.xlabel('T')
plt.ylabel('Expected depth of TreeDSH')
plt.legend(['$len_x(T)$','$len_y(T)$'])
