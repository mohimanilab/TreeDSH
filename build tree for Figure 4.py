#For different values of threshold(thresh), we build the tree and find alpha and beta
import numpy as np
p00=0.25
p01=0.215
p10=0.53
p11=0.005
q00=(p00+p01)*(p00+p10)
q01=(p11+p01)*(p00+p01)
q10=(p10+p11)*(p10+p00)
q11=(p11+p01)*(p11+p10)
pleftTree = p01/q01
prightTree = p10/q10
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


thresh = 80

  

pleftEval = p01
prightEval = p10
qleftEval = q01
qrightEval = q10

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
        node.left = Node(node.val * pleftTree)
        buildTree(node.left)
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

        complexity=(np.log(eval(0.0, root, prightEval, pleftEval, 1.0)))/(np.log(eval(0.0, root, qrightEval, qleftEval, 1.0)))
        print("alpha=", eval(0.0, root, prightEval, pleftEval, 1.0))
        print("beta=", eval(0.0, root, qrightEval, qleftEval, 1.0))




        
        
        
        