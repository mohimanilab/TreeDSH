import numpy as np
Listthresh=[]
Listthreshnew=[]
Listalpha=[]
Listconstant=[]
Listalphatop=[]
Listalphabot=[]
Listbetatop=[]
Listbetabot=[]
Listbeta=[]
Listcomp=[]
p00=0.215
p01=0.0025   
p10=0.255
p11=0.5275
    
q00=(p00+p01)*(p00+p10)
q01=(p11+p01)*(p00+p01)
q10=(p10+p11)*(p10+p00)
q11=(p11+p01)*(p11+p10)
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
    print('11',mu)    
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
    print('10',mu)    
rmin=min(q01/p01,q10/p10,q00/p00,q11/p11)
print(rmin)
rmin=min(q00/p00,q11/p11)
print(rmin)
print(mu)
m=300
for x in range(0,m+1):
    thresh=1+x/10+x**2/100 
    threshnew=np.log(thresh)
    pleftTree = p00/q00
    prightTree = p11/q11

    if pleftTree>1 and prightTree>1 :

        pleftEval = p00
        prightEval = p11

        qleftEval = q00
        qrightEval = q11

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


        def eval(sum, node, pr, pl, val):
            if node.left is None and node.right is None:
                return sum+val
            sl = 0
            sr = 0
            if node.left is not None:
                sl = eval(sum, node.left, pr, pl, val * pl)
            if node.right is not None:
                sr = eval(sum, node.right, pr, pl, val * pr)
            return sl+sr


        if __name__ == "__main__":
            root = Node(1)
            buildTree(root)

    Listthresh.append(thresh)
    Listthreshnew.append(threshnew)
    Listalpha.append(eval(0.0, root, prightEval, pleftEval, 1.0))
    Listbeta.append(eval(0.0, root, qrightEval, qleftEval, 1.0))
    Listcomp.append(np.log(eval(0.0, root, prightEval, pleftEval, 1.0))/np.log(eval(0.0, root, qrightEval, qleftEval, 1.0)))
    Listalphatop.append(((1/thresh))**(mu-1))
    Listalphabot.append(((1/thresh)*rmin)**(mu-1))
    Listbetatop.append(((1/thresh))**(mu))
    Listbetabot.append(((1/thresh)*rmin)**(mu))
    Listconstant.append((mu-1)/mu)


#kist of diffrent values of  threshold
print('Listthresh=',Listthresh)
#threshnew=np.log(thresh)
print('Listthreshnew=',Listthreshnew)
#alpha as a function of threshold
print('Listalpha=',Listalpha)
#beta as a function of threshold
print('Listbeta=',Listbeta)
#complexity as a function of threshold
print('Listcomp=',Listcomp)
#upper bound on alpha as a function of threshold
print('Listalphatop=',Listalphatop)
#lower bound on alpha as a function of threshold
print('Listalphabot=',Listalphabot)
#upper bound on beta as a function of threshold
print('Listbetatop=',Listbetatop)
#lower bound on beta as a function of threshold
print('Listbetabot=',Listbetabot)
#complexity of treeDSH a function of threshold
print('Listconstant=',Listconstant)
    