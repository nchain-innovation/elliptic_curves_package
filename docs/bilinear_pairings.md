# Finite fields

The library implements a class `BilinearPairing` to handle the evaluation of bilinear pairings on pairing-friendly curves. The current implementation is only able to handle bilinear pairings which entail the computation of a single Miller loop, and that do not require further multiplications at the end of the Miller loop (e.g., BLS12 curves and MNT4 curves are supported).

To implement the calculation of a bilinear pairing on a pairing friendly curve, a user only needs to implement the finite fields required by the curve (for which the finite field class are provided), and then implement the final exponentiation for the particular curve of interest.

The code allows some flexibility to customise the way multiplications are carried out in the Miller loop, so to optimise code. However, optimisation is voluntary, and the code can be used out-of-the-box once the finite fields and the final exponentiation are implemented.

Pairings over two curves have been instatiated as way of example:
- BLS12_381
- MNT4_753

The code below shows how to use the instantiations

```python
from elliptic_curves.instantiations.bls12_381 import Fr, bls12_381, easy_exponetiation, hard_exponentiation

# Fr is the scalar field of BLS12_381

# The generator of the group G1
g1 = bls12_381.g1
# The generator of the group G2
g2 = bls12_381.g1

# We can compute the pairing
pairing_g1_g2 = bls12_381.pairing(g1,g2)

# Let's test that it is bilinear
random_element = Fr.generate_random_element().to_list()[0]
assert(pairing_g1_g2.power(random_element) == bls12_381.pairing(g1.multiply(random_element),g2))
assert(pairing_g1_g2.power(random_element) == bls12_381.pairing(g1,g2.multiply(random_element)))

# We can also compute the miller loop

# We can compute it on the base curve and choose the way to eliminate denominators
miller_loop_output_base_curve = bls12_381.miller_loop_on_base_curve(P=g1,Q=g2,denominator_elimination='quadratic')
# We can compute it on the twisted curve and choose not to eliminate denominators
miller_loop_output_twisted_curve = bls12_381.miller_loop_on_base_curve(P=g1,Q=g2,denominator_elimination=None)

# The give rise to the same pairing
assert(hard_exponentiation(easy_exponentiation(miller_loop_output_twisted_curve)) == hard_exponentiation(easy_exponentiation(miller_loop_output_base_curve)))
```