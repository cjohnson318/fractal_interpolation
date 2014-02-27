import numpy as np
import scipy, random 

def d( x ):
    # expects an enumerable, subtracts the right endpoint from the left
    return float(x[-1]-x[0])


def an( x, i ):
    # the (0,0) element in the rotation matrix of the iterated function system (IFS)
    return ( x[i] - x[i-1] )/d(x)


def dn( x, i ):
    # the (0) element in the translation vector of the IFS
    return ( x[i-1] - 0*x[i] )/d(x)


def cn( x, y, i, sn ):
    # the (1,0) element in the rotation matrix of the IFS
    return ( y[i] - y[i-1] )/d(x) - sn*( y[-1] - y[0] )/d(x)


def en( x, y, i, sn ):
    # the (1) element in the translation vector of the IFS
    return ( x[-1]*y[i-1] - x[0]*y[i])/d(x) - sn*( x[-1]+y[0] - x[0]*y[-1] )/d(x)

def Wn( X, U, i, sn ):
    '''
    the iterated function sytem
      R is the rotation matrix
      T is the translation vector
    computes
      R*X + T
    '''
    # rotation matrix
    R = np.matrix([[ an(U[:,0],i), 0 ],\
                   [ cn(U[:,0],U[:,1],i,sn), sn ]])
    # transalation vector
    T = np.matrix([[ dn(U[:,0],i) ],\
                   [ en(U[:,0],U[:,1],i,sn) ]])
    # calculate R*X + T
    tmp = R * np.matrix(X).T + T
    # return the new points
    xp, yp = np.array( tmp.T )[0]
    return xp, yp
    
def T( U, X, deterministic=True ):
    '''
    final transformation that converts points from the
    attractor to points that interpolate the data
    '''
    M = U.shape[0]
    N = X.shape[0]
    if M % 2 == 0:
        s = M - 1
    else:
        s = M - 2
    u, v = list(), list()
    # interpolates between each two points
    for i in range( 1, M ):
        if deterministic:
            x_prime = np.median( X[(i-1)*s:(i)*s,0] )
            y_prime = np.median( X[(i-1)*s:(i)*s,1] )
        else:
            x_prime = random.choice( X[(i-1)*s:(i)*s,0] )
            y_prime = random.choice( X[(i-1)*s:(i)*s,1] )
        dU = U[i,0]-U[i-1,0]
        dX = X[(i)*s,0] - X[(i-1)*s,0]
        u.append( U[i-1,0] + dU*( x_prime - X[(i-1)*s,0] ) / dX )
        v.append( y_prime )
    X = np.vstack((u,v)).T
    X = np.vstack((U,X))
    X = X[ X[:,0].argsort() ]
    return X
    
def G( U, sn, balance=False ):
    # the fractal interpolating function
    X = U.copy()
    x, y = list( X[:,0] ), list( X[:,1] )
    M = U.shape[0]
    N = X.shape[0]
    # for each data point..
    for i in range(N):
        # call an IFS for each segment
        for j in range( 1,M ):
            xp, yp = Wn( X[i], U, j, sn )
            x.append( xp )
            y.append( yp )
            if balance:
                xp, yp = Wn( X[i], U, j, -sn )
                x.append( xp )
                y.append( yp )
    x = np.array(x)
    y = np.array(y)
    # this puts the interpolated
    # data points at the bottom of X
    X = np.vstack((x,y)).T
    X = X[ X[:,0].argsort() ]
    # these two lines rearrage X so that the interpolated
    # data points are between the original data points
    null, indices = np.unique( X[:,0], return_index=True )
    X = X[ indices ]
    return X

def FIF( U, sn, balance=False, deterministic=True ):
    # create attractor
    X = G( U, sn, balance )
    # form interpolation
    X = T( U, X, deterministic )
    return X
