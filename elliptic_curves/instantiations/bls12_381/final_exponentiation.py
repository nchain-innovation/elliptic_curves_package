from elliptic_curves.instantiations.bls12_381.parameters import u

# Final exponentiation --------------------------------------------------------------------------------------------------------

def easy_exponentiation(miller_loop_output):
        '''
        Easy exponentation for BLS12_381: f -> f^{(q^6-1)(q^2+1)}
        '''
        
        a = miller_loop_output.invert() * miller_loop_output.conjugate()
        b = a.frobenius(n=2)

        return a * b

def hard_exponentiation(miller_loop_output):
    '''
    Hard exponentation for BLS12_381
    '''

    t0 = miller_loop_output.power(2)
    t1 = t0.power(u)
    t2 = t1.power(u//2)
    t3 = miller_loop_output.conjugate()
    t1 = t1 * t3
    t1 = t1.conjugate()
    t1 = t1 * t2
    t2 = t1.power(u)
    t3 = t2.power(u)
    t1 = t1.conjugate()
    t3 = t1 * t3
    t1 = t1.invert()
    t1 = t1.frobenius(n=3)
    t2 = t2.frobenius(n=2)
    t1 = t1 * t2
    t2 = t3.power(u)
    t2 = t2 * t0
    t2 = t2 * miller_loop_output
    t1 = t1 * t2 
    t2 = t3.frobenius(n=1)

    return t1 * t2
# -----------------------------------------------------------------------------------------------------------------------------