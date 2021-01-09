import numpy as np


#%%
#   This class implements a non-linear conjugate gradient descent for Total Variation
#
#
class nonlin_conjugateGradient():

    ## data
    k_data = None           # K-space image <N,NCha,NSpk,NSli>
    k_samples = None        # position of the samples in K-space (codified in complex values) <N,N>
    coilprofile = None      # coil profiles <NIx,NIy,NSli,NCha>
    dcf = None              # density compensation function <N,NSpk,NSli>

    num_spokes_session = []  # each entry indicates the number of spokes in each session
    spokes_per_vol = None   # Number of spokes per volume


    ## Conjugate gradient params
    # line search parameters
    maxline_iter = 150
    gradToll = 1e-3
    l1Smooth = 1e-15
    alpha = 0.01
    beta = 0.6
    t0 = 1
    k = 0

    # generatl CG parameters
    maxIter = 1e3
    regParam = 0

    # objective parameters
    E = None
    W = None
    y = None

    def __init__(self, x0, E, W, y, regParam):
        self.E = E
        self.W = W
        self.y = y
        self.regParam = regParam


    def run(self, x0):
        

        # set initial parameters
        x_k = x0
        g_k = self.gradient( x_k )
        dx = - g_k
        k = 0
        # ITERATE!
        while(True):
            
            # compute function evaluation
            f0 = self.objective(x_k) 
            t = self.t0
            f1 = self.objective( x_k + t*dx )

            # run line search
            line_iter = 0
            while ( f1 > f0 - self.alpha*t*np.abs( g_k.dot( dx ) ) ):
                # evaluate at next point
                t = t * self.beta
                f1 = self.objective( x_k + t*dx )
                line_iter += 1

                if ( line_iter == self.maxline_iter ):
                    print 'Line Search error. Max iterations found.'
                    break

            # do an adaptation to the line search length
            if line_iter > 2:
                self.t0 *= self.beta
            elif line_iter < 1:
                self.t0 *= 1./float(self.beta)

            # update
            x_k += t*dx

            print 'Iter: %d, f(x)=%0.6f' %(k,f1)
            if (k> self.maxIter) | (np.linalg.norm(dx) < self.gradToll):
                print 'Gato Dominguez!'
                break

            k += 1

    def objective(self, x):

        # residual norm
        residual = self.E.dot( x ) - self.y
        resObj = np.linalg.norm(residual,2)

        # TV regularization 
        if self.regParam > 0:
            w = self.W.dot(x)
            TVobj = np.sum(np.sqrt( np.conj(w)*w + self.l1Smooth ))
        else:
            TVobj = 0

        return resObj  + self.regParam*TVobj

    def gradient(self, x):

        # residual gradient
        residual = self.E.dot( x ) - self.y
        resGrad = 2 * self.E.T.dot( residual )

        # TV regularization 
        if self.regParam > 0:
            w = self.W.dot(x)
            TVGrad = self.W.T.dot( np.sqrt( np.conj(w)*w + self.l1Smooth ) )
        else:
            TVGrad = 0

        return resGrad  + self.regParam*TVGrad
