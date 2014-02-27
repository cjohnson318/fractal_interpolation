fractal_interpolation
=====================

This implements a technique for curve fitting by fractal interpolation found in a paper by Manousopoulos, Drakopoulos, and Theoharis, found [here](http://graphics.di.uoa.gr/Downloads/papers/journals/p30.pdf). I also used infromation about non-linear fractal interpolating functions found [here](http://www.mi.sanu.ac.rs/vismath/kobes/).


Usage
=====

Here, `U` is the original data. **Line 5** "fractalizes" example data, and **line 6** performs the interpolation.

    u = np.array([0.0,0.4,0.7,1.0])
    v = np.array([0.0,0.5,0.2,0.0])
    
    U = np.vstack((u,v)).T
    U = G( U, 0.1, balance=0 )
    X = FIF( U, 0.01, balance=1 )
    
    plot( U[:,0], U[:,1], '.-' )
    plot( X[:,0], X[:,1], '.-' ) 
