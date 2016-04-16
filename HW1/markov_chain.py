import numpy as np
from os import walk
mypath = 'proteins/'  # use path to data files
_, _, filenames = next(walk(mypath), (None, None, []))

mSeq = len(filenames)        # read in each sequence
o,x = [],[]
for i in range(mSeq):
    f = open( mypath+'/'+filenames[i] , 'r')
    o.append( f.readline()[:-1] )  # strip trailing '\n'
    x.append( f.readline()[:-1] )
    f.close()

xvals, ovals = set(),set()  # extract the symbols used in x and o
for i in range(mSeq):
    xvals |= set(x[i])
    ovals |= set(o[i])
xvals = list( np.sort( list(xvals) ) )
ovals = list( np.sort( list(ovals) ) )
dx,do = len(xvals),len(ovals)

for i in range(mSeq):       # and convert to numeric indices
    x[i] = np.array([xvals.index(s) for s in x[i]])
    o[i] = np.array([ovals.index(s) for s in o[i]])

def initialStateDistribution(x):
    counts = [0]*dx
    for xi in x:
        try:
            counts[xi[0]]+=1
        except:
            counts[xi[0]]=1
    p0 = np.array(counts)*1.0/sum(counts)
    return p0

def stateTransitionProbabilities(x):
    T = np.zeros((dx,dx))
    for xp in x:
        for j in range(0,len(xp)-1,1):
            xi = xp[j]
            xj = xp[j+1]
            T[xi,xj]+=1
    T = np.transpose(np.transpose(T)/np.sum(T,axis=1)*1.0)
    return T

def emissionProbabilities(x,o):
    O = np.zeros((dx,do))
    for xi,oi in zip(x,o):
        for xij,oij in zip(xi,oi):
            O[xij,oij]+=1
    O = O.T/np.sum(O.T,axis=0)*1.0
    O = O.T
    return O

def markovMarginals(x,o,p0,Tr,Ob):
    '''Compute p(o) and the marginal probabilities p(x_t|o) for a Markov model
       defined by P[xt=j|xt-1=i] = Tr(i,j) and P[ot=k|xt=i] = Ob(i,k) as numpy matrices'''
    dx,do = Ob.shape   # if a numpy matrix
    L = len(o)
    f = np.zeros((L,dx))
    r = np.zeros((L,dx))
    p = np.zeros((L,dx))
    Z = np.zeros(L)
    f[0,:] = (Ob[:,o[0]] * p0) # compute initial forward message
    #print f[0,:]
    Z[0] = f[0,:].sum()
    log_pO = np.log(Z[0])  # update probability of sequence so far
    
    f[0,:] /= Z[0] # normalize (to match definition of f)
    
    for t in range(1,L): # compute forward messages
        fp = f[t-1,:][:,np.newaxis]
        f[t,:] = Ob[:,o[t]] * (np.dot(Tr.T,fp).T)
        Z[t] = f[t,:].sum()
        log_pO += np.log(Z[t])
        f[t,:] /= Z[t]
    print log_pO
    r[L-1,:] = np.ones(dx) # initialize reverse messages
    #p[L,:] = ... # and marginals READ: (this term is taken care of in the loop itself)

    for t in range(L-1,-1,-1):
        if t<=L-2:
            r[t,:] = Tr.dot(Ob[:,o[t+1]]*(r[t+1,:]))
        r[t,:] /= r[t,:].sum()
        p[t,:] = f[t,:]*r[t,:]
        p[t,:] /= p[t,:].sum()*1.0
    return p,log_pO

    
p0 = initialStateDistribution(x)
Tr = stateTransitionProbabilities(x)
stationary_distribution = np.linalg.matrix_power(Tr,1000)
stationary_distribution = (stationary_distribution.T/np.sum(stationary_distribution.T,axis=0)).T
stationary_distribution = stationary_distribution[0]
Ob = emissionProbabilities(x,o)
print markovMarginals(x,o[4],p0,Tr,Ob)
