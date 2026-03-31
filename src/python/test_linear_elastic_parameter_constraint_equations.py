from linear_elastic_parameter_constraint_equations import *
import numpy as np
from numpy.random import rand

def checkDictionary( answerDict, resultDict ):
    """
    Check two dictionaries to make sure they are the same.

    :param dict answerDict: The answer dictionary
    :param dict resultDict: The result dictionary
    """

    answerKeys = [ v for v in answerDict.keys() ]
    resultKeys = [ v for v in resultDict.keys() ]

    if ( len( answerKeys ) != len( resultKeys ) ):
        assert( False )

    for key in answerDict.keys():
        assert( key in resultDict.keys())
        assert( np.isclose( resultDict[ key ], answerDict[ key ] ) )

def finiteDifferenceConstraint( fxn, vals, keys, eps=1e-6 ):
    """
    Compute the finite difference of a constraint function

    :param function fxn: The function to evaluate.
    :param np.array vals: The values of the parameters in the order
        they appear in the function call
    :param dict keys: The keynames for the variables.
    :param float eps: The delta to apply to the values
    """

    jacobian = dict( [ ( key, None ) for key in keys ] )

    vals = np.array( vals )

    for i, v in enumerate( vals ):

        delta = np.zeros( vals.size )
        delta[ i ] = eps * abs( v ) + eps

        fp = fxn( *( vals + delta) )[ 0 ]
        fm = fxn( *( vals - delta) )[ 0 ]

        jacobian[ keys[ i ] ] = ( fp - fm ) / ( 2 * delta[ i ] )

    return jacobian

def finiteDifferenceTMatrix( vals, eps = 1e-6 ):
    """
    Compute the finite difference of the T matrix

    :param np.array vals: The values of the taus
    """
    
    jacobian = np.zeros( ( 3, 3, vals.size ) )

    for i, v in enumerate( vals ):

        delta = np.zeros( vals.size )
        delta[ i ] = eps * abs( v ) + eps

        fp = construct_T_matrix( vals + delta )[ 0 ]
        fm = construct_T_matrix( vals - delta )[ 0 ]

        jacobian[ :, :, i ] = ( fp - fm ) / ( 2 * delta[ i ] )

    return jacobian

def test_evaluate_g1():

    mu = rand()
    s = rand()

    answerR = mu - s**2
    answerJ = finiteDifferenceConstraint( evaluate_g1, [ mu, s ], [ 'mu', 's1' ] )

    resultR, resultJ = evaluate_g1( mu, s )

    assert( np.isclose( answerR, resultR ) )

    checkDictionary( answerJ, resultJ )

def test_evaluate_g2():

    k  = rand()
    n  = rand()
    si = rand()
    s  = rand()

    answerR = k + n - 2 * si - s**2
    answerJ = finiteDifferenceConstraint( evaluate_g2, [ k, n, si, s ], [ 'kappa', 'nu', 'sigma', 's2' ] )

    resultR, resultJ = evaluate_g2( k, n, si, s )

    assert( np.isclose( answerR, resultR ) )

    checkDictionary( answerJ, resultJ )

def test_evaluate_g3():

    k  = rand()
    n  = rand()
    si = rand()
    m  = rand()
    s  = rand()

    answerR = ( k + n - 2 * si ) * m - 2 * si**2 - s**2
    answerJ = finiteDifferenceConstraint( evaluate_g3, [ k, n, si, m, s ], [ 'kappa', 'nu', 'sigma', 'mu', 's3' ] )

    resultR, resultJ = evaluate_g3( k, n, si, m, s )

    assert( np.isclose( answerR, resultR ) )

    checkDictionary( answerJ, resultJ )

def test_evaluate_g4():

    l = rand()
    m = rand()
    s = rand()

    answerR = ( 3 * l + 2 * m ) - s**2
    answerJ = finiteDifferenceConstraint( evaluate_g4, [ l, m, s ], [ 'lambda', 'mu', 's4' ] )

    resultR, resultJ = evaluate_g4( l, m, s )

    assert( np.isclose( answerR, resultR ) )

    checkDictionary( answerJ, resultJ )

def test_evaluate_g5():

    k  = rand()
    n  = rand()
    e  = rand()
    t  = rand()
    si = rand()
    s  = rand()

    answerR = ( k + n + 3 * e - 3 * t - 2 * si ) - s**2
    answerJ = finiteDifferenceConstraint( evaluate_g5, [ k, n, e, t, si, s ], [ 'kappa', 'nu', 'eta', 'tau', 'sigma', 's5' ] )

    resultR, resultJ = evaluate_g5( k, n, e, t, si, s )

    assert( np.isclose( answerR, resultR ) )

    checkDictionary( answerJ, resultJ )

def test_evaluate_g6():

    k  = rand()
    n  = rand()
    e  = rand()
    t  = rand()
    si = rand()
    l  = rand()
    m  = rand()
    s  = rand()

    answerR = ( k + n + 3 * e - 3 * t - 2 * si ) * ( 3 * l + 2 * m ) - ( 3 * t + 2 * si )**2 - s**2
    answerJ = finiteDifferenceConstraint( evaluate_g6, [ k, n, e, t, si, l, m, s ],\
                  [ 'kappa', 'nu', 'eta', 'tau', 'sigma', 'lambda', 'mu', 's6' ] )

    resultR, resultJ = evaluate_g6( k, n, e, t, si, l, m, s )

    assert( np.isclose( answerR, resultR ) )

    checkDictionary( answerJ, resultJ )

def test_evaluate_g7():

    k = rand()
    n = rand()
    s = rand()

    answerR = k - n - s**2
    answerJ = finiteDifferenceConstraint( evaluate_g7, [ k, n, s ], [ 'kappa', 'nu', 's7' ] )

    resultR, resultJ = evaluate_g7( k, n, s )

    assert( np.isclose( answerR, resultR ) )

    checkDictionary( answerJ, resultJ )

def test_evaluate_g8():

    m  = rand()
    k  = rand()
    n  = rand()
    si = rand()
    s  = rand()

    answerR = 4 * m * ( k + n - 2 * si ) - 2 * si - s**2
    answerJ = finiteDifferenceConstraint( evaluate_g8, [ m, k, n, si, s ], [ 'mu', 'kappa', 'nu', 'sigma', 's8' ] )

    resultR, resultJ = evaluate_g8( m, k, n, si, s )

    assert( np.isclose( answerR, resultR ) )

    checkDictionary( answerJ, resultJ )

def test_evaluate_g9():

    taus = np.random.rand( 11 )

    t1  = taus[  0 ]
    t2  = taus[  1 ]
    t3  = taus[  2 ]
    t4  = taus[  3 ]
    t5  = taus[  4 ]
    t6  = taus[  5 ]
    t7  = taus[  6 ]
    t8  = taus[  7 ]
    t9  = taus[  8 ]
    t10 = taus[  9 ]
    t11 = taus[ 10 ]
    s   = rand()

    answerR = t7 + 2 * t8 - abs( t9 + t10 + t11 ) - s**2
    answerJ = finiteDifferenceConstraint( evaluate_g9, [ t7, t8, t9, t10, t11, s ], [ 'tau7', 'tau8', 'tau9', 'tau10', 'tau11', 's9' ] )

    resultR, resultJ = evaluate_g9( t7, t8, t9, t10, t11, s )

    assert( np.isclose( answerR, resultR ) )

    checkDictionary( answerJ, resultJ )

def test_evaluate_g10():

    taus = np.random.rand( 11 )

    t1  = taus[  0 ]
    t2  = taus[  1 ]
    t3  = taus[  2 ]
    t4  = taus[  3 ]
    t5  = taus[  4 ]
    t6  = taus[  5 ]
    t7  = taus[  6 ]
    t8  = taus[  7 ]
    t9  = taus[  8 ]
    t10 = taus[  9 ]
    t11 = taus[ 10 ]
    s   = rand()

    answerR = t7 - t8 - ( 1. / np.sqrt( 2 ) ) * np.sqrt( ( t9 - t10 )**2 + ( t10 - t11 )**2 + ( t11 - t9 )**2 ) - s**2
    answerJ = finiteDifferenceConstraint( evaluate_g10, [ t7, t8, t9, t10, t11, s ], [ 'tau7', 'tau8', 'tau9', 'tau10', 'tau11', 's10' ] )

    resultR, resultJ = evaluate_g10( t7, t8, t9, t10, t11, s )

    assert( np.isclose( answerR, resultR ) )

    checkDictionary( answerJ, resultJ )

def test_construct_T_matrix():

    taus = np.random.rand( 11 )

    t1  = taus[  0 ]
    t2  = taus[  1 ]
    t3  = taus[  2 ]
    t4  = taus[  3 ]
    t5  = taus[  4 ]
    t6  = taus[  5 ]
    t7  = taus[  6 ]
    t8  = taus[  7 ]
    t9  = taus[  8 ]
    t10 = taus[  9 ]
    t11 = taus[ 10 ]

    answerR = np.array( [ [ t1 + t2 + 3 * t3 + t7 + t10, 3 * t1 + t4 + t5 + t8 + t11, 3 * t2 + t5 + t6 + t8 + t9 ],
                          [ 3 * t1 + t2 + t3 + t8 + t11, t1 + 3 * t4 + t5 + t7 + t9, t2 + 3 * t5 + t6 + t8 + t10 ],
                          [ t1 + 3 * t2 + t3 + t8 + t9, t1 + t4 + 3 * t5 + t8 + t10, t2 + t5 + 3 * t6 + t7 + t11 ] ] )
    answerJ = finiteDifferenceTMatrix( taus )

    resultR, resultJ = construct_T_matrix( taus )

    assert( np.allclose( answerR, resultR ) )

    assert( np.allclose( answerJ, resultJ ) )

def test_evaluate_g11():

    taus = np.random.rand( 11 )

    t1  = taus[  0 ]
    t2  = taus[  1 ]
    t3  = taus[  2 ]
    t4  = taus[  3 ]
    t5  = taus[  4 ]
    t6  = taus[  5 ]
    t7  = taus[  6 ]
    t8  = taus[  7 ]
    t9  = taus[  8 ]
    t10 = taus[  9 ]
    t11 = taus[ 10 ]
    s   = rand()

    T, dTdtau = construct_T_matrix( taus );

    answerR = np.trace( T ) - s**2 - 1e-9
    answerJ = finiteDifferenceConstraint( evaluate_g11, [ v for v in taus] + [ s ], [ 'tau{0}'.format(i+1) for i in range( taus.size ) ] + [ 's11' ] )

    resultR, resultJ = evaluate_g11( *taus, s )

    checkDictionary( answerJ, resultJ )

def test_evaluate_g12():

    taus = np.random.rand( 11 )

    t1  = taus[  0 ]
    t2  = taus[  1 ]
    t3  = taus[  2 ]
    t4  = taus[  3 ]
    t5  = taus[  4 ]
    t6  = taus[  5 ]
    t7  = taus[  6 ]
    t8  = taus[  7 ]
    t9  = taus[  8 ]
    t10 = taus[  9 ]
    t11 = taus[ 10 ]
    s   = rand()

    T, dTdtau = construct_T_matrix( taus );

    answerR = np.trace( np.linalg.inv( T ).T * np.linalg.det( T ) ) - s**2 - 1e-9

    answerJ = finiteDifferenceConstraint( evaluate_g12, [ v for v in taus] + [ s ], [ 'tau{0}'.format(i+1) for i in range( taus.size ) ] + [ 's12' ] )

    resultR, resultJ = evaluate_g12( *taus, s )

    assert( np.allclose( answerR, resultR ) )

    checkDictionary( answerJ, resultJ )

def test_evaluate_g13():

    taus = np.random.rand( 11 )

    t1  = taus[  0 ]
    t2  = taus[  1 ]
    t3  = taus[  2 ]
    t4  = taus[  3 ]
    t5  = taus[  4 ]
    t6  = taus[  5 ]
    t7  = taus[  6 ]
    t8  = taus[  7 ]
    t9  = taus[  8 ]
    t10 = taus[  9 ]
    t11 = taus[ 10 ]
    s   = rand()

    T, dTdtau = construct_T_matrix( taus );

    answerR = np.linalg.det( T ) - s**2 - 1e-9

    answerJ = finiteDifferenceConstraint( evaluate_g13, [ v for v in taus] + [ s ], [ 'tau{0}'.format(i+1) for i in range( taus.size ) ] + [ 's13' ] )

    resultR, resultJ = evaluate_g13( *taus, s )

    assert( np.allclose( answerR, resultR ) )

    checkDictionary( answerJ, resultJ )
