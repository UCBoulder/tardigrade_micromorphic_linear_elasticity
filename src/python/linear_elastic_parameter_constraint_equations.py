"""
A file which defines the constraint equations for the parameters required for fitting 
"""
import numpy as np

def evaluate_g1( mu, s1, **kwargs ):
    """
    Evaluate the first constraint equation and also return the jacobians

    :param float l: The value of the modulus mu
    :param float s1: The constraint value
    """

    return mu - s1**2, {'mu':1., 's1':- 2 * s1 }

def evaluate_g2( kappa, nu, sigma, s2, **kwargs ):
    """
    Evaluate the second constraint equation and also return the jacobians

    :param float kappa: The value of the modulus kappa
    :param float nu: The value of the modulus nu
    :param float sigma: The value of the modulus sigma
    :param float s2: The constraint value
    """

    return kappa + nu - 2 * sigma - s2**2, {'kappa':1., 'nu':1., 'sigma':-2., 's2':-2 * s2 }

def evaluate_g3( kappa, nu, sigma, mu, s3, **kwargs ):
    """
    Evaluate the third constraint equation and also return the jacobians

    :param float kappa: The value of the modulus kappa
    :param float nu: The value of the modulus nu
    :param float sigma: The value of the modulus sigma
    :param float mu: The value of the modulus mu
    :param float s3: The constraint value
    """

    return ( kappa + nu - 2 * sigma ) * mu - 2 * sigma**2 - s3**2,\
           {'kappa':mu, 'nu':mu, 'sigma':-2 * mu - 4 * sigma, 'mu':( kappa + nu - 2 * sigma), 's3':-2 * s3 }

def evaluate_g4( l, mu, s4, **kwargs ):
    """
    Evaluate the fourth constraint equation and also return the jacobians
    Evaluate the constraint equation for lambda and mu and
    also return the jacobians

    :param float l: The value of the modulus lambda
    :param float mu: The value of the modulus mu
    :param float s4: The value of the constraint
    """

    return 3 * l + 2 * mu - s4**2, { 'lambda':3., 'mu':2., 's4':-2 * s4 }

def evaluate_g5( kappa, nu, eta, tau, sigma, s5, **kwargs ):
    """
    Evaluate the fifth constraint equation and also return the jacobians

    :param float kappa: The value of the modulus kappa
    :param float nu: The value of the modulus nu
    :param float eta: The value of the modulus eta
    :param float tau: The value of the modulus tau
    :param float sigma: The value of the modulus sigma
    :param float s5: The value of the constraint
    """

    return kappa + nu + 3 * eta - 3 * tau - 2 * sigma - s5**2,\
           { 'kappa':1., 'nu':1., 'eta':3., 'tau':-3., 'sigma':-2., 's5':-2*s5 }

def evaluate_g6( kappa, nu, eta, tau, sigma, l, mu, s6, **kwargs ):
    """
    Evaluate the sixth constraint equation and also return the jacobians
    
    :param float kappa: The value of the modulus kappa
    :param float nu: The value of the modulus nu
    :param float eta: The value of the modulus eta
    :param float tau: The value of the modulus tau
    :param float sigma: The value of the modulus sigma
    :param float l: The value of the modulus lambda
    :param float mu: The value of the modulus mu
    :param float s6: The value of the constraint
    """

    return ( kappa + nu + 3 * eta - 3 * tau - 2 * sigma ) * ( 3 * l + 2 * mu ) - ( 3 * tau + 2 * sigma )**2 - s6**2,\
           { 'kappa':( 3 * l + 2 * mu ), 'nu':( 3 * l + 2 * mu ), 'eta':3 * ( 3 * l + 2 * mu ),\
             'tau':-3 * ( 3 * l + 2 * mu) - 6 * ( 3 * tau + 2 * sigma ),\
             'sigma':-2 * ( 3 * l + 2 * mu ) - 4 * ( 3 * tau + 2 * sigma ),\
             'lambda':3 * ( kappa + nu + 3 * eta - 3 * tau - 2 * sigma ),\
             'mu':2 * ( kappa + nu + 3 * eta - 3 * tau - 2 * sigma ),\
             's6':-2 * s6 }

def evaluate_g7( kappa, nu, s7, **kwargs ):
    """
    Evaluate the seventh constraint equation and also return the jacobian

    :param float kappa: The value of the modulus kappa
    :param float nu: The value of the modulus nu
    :param float s7: The value of the constraint
    """

    return kappa - nu - s7**2, {'kappa':1., 'nu':-1., 's7':-2*s7 }

def evaluate_g8( mu, kappa, nu, sigma, s8, **kwargs ):
    """
    Evaluate the eighth constraint equation and also return the jacobian

    :param float mu: The value of the modulus mu
    :param float kappa: The value of the modulus kappa
    :param float nu: The value of the modulus nu
    :param float sigma: The value of the modulus sigma
    :param float s8: The value of the constraint
    """

    return 4 * mu * ( kappa + nu - 2 * sigma ) - 2 * sigma - s8**2,\
           { 'mu':4 * ( kappa + nu - 2 * sigma ), 'kappa': 4 * mu, 'nu': 4 * mu, 'sigma':-8 * mu - 2, 's8':-2 * s8 }

def evaluate_g9( tau7, tau8, tau9, tau10, tau11, s9, **kwargs ):
    """
    Evaluate the ninth constraint equation and also return the Jacobian

    :param float tau7: The seventh tau parameter
    :param float tau8: The eighth tau parameter
    :param float tau9: The ninth tau parameter
    :param float tau10: The tenth tau parameter
    :param float tau11: The eleventh tau parameter
    :param float s9: The value of the constraint
    """

    return tau7 + 2 * tau8 - abs( tau9 + tau10 + tau11 ) - s9**2,\
            { 'tau7':1., 'tau8':2., 'tau9':float( -np.sign( tau9 ) ),\
              'tau10':float( -np.sign( tau10 ) ),\
              'tau11':( -np.sign( tau11 ) ), 's9':-2*s9 }

def evaluate_g10( tau7, tau8, tau9, tau10, tau11, s10, **kwargs ):
    """
    Evaluate the tenth constraint equation and also return the Jacobian

    :param float tau7: The seventh tau parameter
    :param float tau8: The eighth tau parameter
    :param float tau9: The ninth tau parameter
    :param float tau10: The tenth tau parameter
    :param float tau11: The eleventh tau parameter
    :param float s10: The value of the constraint
    """

    term1 = np.sqrt( ( tau9 - tau10 )**2 + ( tau10 - tau11 )**2 + ( tau11 - tau9 )**2 )
    dterm1dint = ( 0.5 / term1 )

    term = ( 1. / np.sqrt( 2. ) ) * term1
    dtermdint = ( 1. / np.sqrt( 2. ) ) * dterm1dint

    g10 = tau7 - tau8 - term - s10**2
    
    dg10dtau7  = 1.
    dg10dtau8  = -1.
    dg10dtau9  = -2 * dtermdint * ( (tau9 - tau10 )  - ( tau11 - tau9 ) )
    dg10dtau10 = -2 * dtermdint * ( ( tau10 - tau11 ) - ( tau9 - tau10 ) )
    dg10dtau11 = -2 * dtermdint * ( ( tau11 - tau9 ) - ( tau10 - tau11 ) )
    dg10ds10   = -2 * s10

    return g10, { 'tau7':dg10dtau7, 'tau8':dg10dtau8, 'tau9':dg10dtau9, 'tau10':dg10dtau10, 'tau11':dg10dtau11, 's10':dg10ds10 }

def evaluate_g11( tau1, tau2, tau3, tau4, tau5, tau6, tau7, tau8, tau9, tau10, tau11, s11, tol=1e-9, **kwargs ):
    """
    Construct the eleventh constraint equation

    :param np.array tau1-11: The taus
    :param float s11: The constraint value
    """

    taus = np.array( [ tau1, tau2, tau3, tau4, tau5, tau6, tau7, tau8, tau9, tau10, tau11 ] )

    T, dTdtaus = construct_T_matrix( taus )

    g11 = np.trace( T ) - s11**2

    dg11dtau = np.einsum( 'ij,ijk->k', np.eye( 3 ), dTdtaus )
    dg11ds11 = -2*s11

    return g11, dict( [( 'tau{0}'.format(i+1), dg11dtau[i] ) for i in range( taus.size ) ] + [ ('s11', dg11ds11) ] )

def evaluate_g12( tau1, tau2, tau3, tau4, tau5, tau6, tau7, tau8, tau9, tau10, tau11, s12, tol=1e-9, **kwargs ):
    """
    Construct the twelth constraint equation

    :param np.array tau1-11: The taus
    :param float s12: The constraint value
    :param float tol: The tolerance. Because this calculation involves an inverse, the
        resulting T matrix must be invertable.
    """

    taus = np.array( [ tau1, tau2, tau3, tau4, tau5, tau6, tau7, tau8, tau9, tau10, tau11 ] )

    T, dTdtaus = construct_T_matrix( taus )
    invT = np.linalg.inv( T )
    detT = np.linalg.det( T )

    g12 = np.trace( invT.T * detT ) - s12**2 - tol

    dg12dT = detT * ( np.trace( invT.T ) * invT - np.dot( invT, invT ) ).T

    dg12dtaus = np.einsum( 'kl,klm', dg12dT, dTdtaus )
    dg12ds12 = -2*s12

    return g12, dict( [( 'tau{0}'.format(i+1), dg12dtaus[i] ) for i in range( taus.size ) ] + [ ('s12', dg12ds12) ] )

def evaluate_g13( tau1, tau2, tau3, tau4, tau5, tau6, tau7, tau8, tau9, tau10, tau11, s13, tol=1e-9, **kwargs ):
    """
    Construct the twelth constraint equation

    :param np.array tau1-11: The taus
    :param float s13: The constraint value
    :param float tol: The tolerance. Because this calculation involves an inverse, the
        resulting T matrix must be invertable.
    """
 
    taus = np.array( [ tau1, tau2, tau3, tau4, tau5, tau6, tau7, tau8, tau9, tau10, tau11 ] )

    T, dTdtaus = construct_T_matrix ( taus )
    invT = np.linalg.inv( T )
    detT = np.linalg.det( T )
    g13 = detT - s13**2

    dg13dT = detT * invT.T

    dg13dtaus = np.einsum( 'ij,ijk->k', dg13dT, dTdtaus )

    dg13ds13 = -2*s13

    return g13, dict( [( 'tau{0}'.format(i+1), dg13dtaus[i] ) for i in range( taus.size ) ] + [ ('s13', dg13ds13) ] )


def construct_T_matrix( taus ):
    """
    Construct the tau matrix from the taus

    :param np.array taus: The taus ( should be of length 11 )
    """

    if ( taus.size != 11 ):
        raise ValueError( "The tau vector must have 11 elements" )

    t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11 = taus

    T = np.array ( [ [ t1 + t2 + 3 * t3 + t7 + t10, 3 * t1 + t4 + t5 + t8 + t11, 3 * t2 + t5 + t6 + t8 + t9 ],\
                     [ 3 * t1 + t2 + t3 + t8 + t11, t1 + 3 * t4 + t5 + t7 + t9, t2 + 3 * t5 + t6 + t8 + t10 ],\
                     [ t1 + 3 * t2 + t3 + t8 + t9, t1 + t4 + 3 * t5 + t8 + t10, t2 + t5 + 3 * t6 + t7 + t11 ] ] );

    jacobian = np.zeros( ( 3, 3, 11 ) )

    jacobian[ :, :, 0 ] = np.array( [ [ 1, 3, 0 ],
                                      [ 3, 1, 0 ],
                                      [ 1, 1, 0 ] ] ).astype( float )

    jacobian[ :, :,  1 ] = np.array( [ [ 1, 0, 3 ],
                                       [ 1, 0, 1 ],
                                       [ 3, 0, 1 ] ] ).astype( float )

    jacobian[ :, :,  2 ] = np.array( [ [ 3, 0, 0 ],
                                       [ 1, 0, 0 ],
                                       [ 1, 0, 0 ] ] ).astype( float )

    jacobian[ :, :,  3 ] = np.array( [ [ 0, 1, 0 ],
                                       [ 0, 3, 0 ],
                                       [ 0, 1, 0 ] ] ).astype( float )

    jacobian[ :, :,  4 ] = np.array( [ [ 0, 1, 1 ],
                                       [ 0, 1, 3 ],
                                       [ 0, 3, 1 ] ] ).astype( float )

    jacobian[ :, :,  5 ] = np.array( [ [ 0, 0, 1 ],
                                       [ 0, 0, 1 ],
                                       [ 0, 0, 3 ] ] ).astype( float )

    jacobian[ :, :,  6 ] = np.array( [ [ 1, 0, 0 ],
                                       [ 0, 1, 0 ],
                                       [ 0, 0, 1 ] ] ).astype( float )

    jacobian[ :, :,  7 ] = np.array( [ [ 0, 1, 1 ],
                                       [ 1, 0, 1 ],
                                       [ 1, 1, 0 ] ] ).astype( float )

    jacobian[ :, :,  8 ] = np.array( [ [ 0, 0, 1 ],
                                       [ 0, 1, 0 ],
                                       [ 1, 0, 0 ] ] ).astype( float )

    jacobian[ :, :,  9 ] = np.array( [ [ 1, 0, 0 ],
                                       [ 0, 0, 1 ],
                                       [ 0, 1, 0 ] ] ).astype( float )

    jacobian[ :, :, 10 ] = np.array( [ [ 0, 1, 0 ],
                                       [ 1, 0, 0 ],
                                       [ 0, 0, 1 ] ] ).astype( float )

    return T, jacobian
