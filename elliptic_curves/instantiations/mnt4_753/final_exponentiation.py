from .parameters import u

# Final exponentiation --------------------------------------------------------------------------------------------------------

def easy_exponentiation(miller_loop_output):
    """
    Easy exponentiation for MNT4_753 is f -> f^{q^2-1}
    """

    out = miller_loop_output.frobenius(2) * miller_loop_output.invert()
    return out

def hard_exponentiation(cyclotomic_element):
    """
    Hard exponentiation for MNT4_753 is f -> f^{q + u + 1}
    """

    return cyclotomic_element.frobenius(1) * cyclotomic_element.power(u) * cyclotomic_element
    
# -----------------------------------------------------------------------------------------------------------------------------