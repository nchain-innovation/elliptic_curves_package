import os, sys
sys.path.append(os.path.normpath(os.path.join(os.path.dirname(__file__),'../../')))

from elliptic_curves.fields.fq import base_field_from_modulus
from elliptic_curves.fields.quadratic_extension import quadratic_extension_from_base_field_and_non_residue
from elliptic_curves.fields.cubic_extension import cubic_extension_from_base_field_and_non_residue

from elliptic_curves.models.ec import elliptic_curve_from_curve
from elliptic_curves.models.curve import Curve, BilinearPairingCurve
from elliptic_curves.models.bilinear_pairings import BilinearPairing

from elliptic_curves.instantiations.bls12_381.parameters import *
from elliptic_curves.instantiations.bls12_381.final_exponentiation import easy_exponentiation, hard_exponentiation

# Field instantiation
Fq = base_field_from_modulus(q=q)
NON_RESIDUE_FQ = Fq.from_list(NON_RESIDUE_FQ)
Fq2 = quadratic_extension_from_base_field_and_non_residue(base_field=Fq,non_residue=NON_RESIDUE_FQ)
NON_RESIDUE_FQ2 = Fq2.from_list(NON_RESIDUE_FQ2)
Fq6 = cubic_extension_from_base_field_and_non_residue(base_field=Fq2,non_residue=NON_RESIDUE_FQ2)
NON_RESIDUE_FQ6 = Fq6.from_list(NON_RESIDUE_FQ6)
Fq12 = quadratic_extension_from_base_field_and_non_residue(base_field=Fq6,non_residue=NON_RESIDUE_FQ6)

# Define mul_by_line_eval method for Fq12
# We are not customising line_evaluations, so they will be elements in the miller_ouput_type
def mul_by_line_eval(self,line_eval):
    return self * line_eval

Fq12.mul_by_line_eval = mul_by_line_eval

# Scalar field of the curve
Fr = base_field_from_modulus(q=r)

# Curves
bls12_381_curve = Curve(a = Fq(a), b = Fq(b))
bls12_381_twisted_curve = Curve(a = Fq2.zero(), b = Fq(b) * NON_RESIDUE_FQ2)

# Curve classes
BLS12_381, _ = elliptic_curve_from_curve(curve=bls12_381_curve)
BLS12_381_Twist, _ = elliptic_curve_from_curve(curve=bls12_381_twisted_curve)

# Twisting morphisms
def to_twisted_curve(self):
    '''
    Untwisting morphism Psi : E' --> E, (x',y') --> (x'/omega^2, y'/omega^3)
    Here omega in F_q^k such that the equation of the twist is y^2 = x^3 + a omega^4 x + b omega^6
    For BLS12, this means: y^2 = x^3 + b omega^6.

    We are in the M-twist case, with equation: y^2 = x^3 + b * xi, hence omega^6 = xi.
    F_q^12 = F_q^6[w] / (w^2 - v) = F_q^2[w,v] / (w^2 - v, v^3 - xi)
    Hence, omega^6 = w, and the untwisting morphism is Psi(x',y') = (x' / w^2, y' / w^3) in E(F_q^12)
    Thus, the twisting morphism is Phi : E --> E' : (x,y) --> (x * w^2, y * w^3) in E'(F_q^12)
    '''
    
    return BLS12_381_Twist(self.x * Fq12.u().power(2), self.y * Fq12.u().power(3))

def to_base_curve(self):
    '''
    Ref. to_twisted_curve => omega = w
    '''

    return BLS12_381(self.x * Fq12.u().power(-2), self.y * Fq12.u().power(-3))

BLS12_381.to_twisted_curve = to_twisted_curve
BLS12_381_Twist.to_base_curve = to_base_curve

# BilinearPairing
bls12_381 = BilinearPairingCurve(
    q = q,
    r = r,
    h1 = h1,
    h2 = h2,
    curve = bls12_381_curve,
    twisted_curve = bls12_381_twisted_curve,
    g1 = BLS12_381(x = Fq(g1_X), y = Fq(g1_Y)),
    g2 = BLS12_381_Twist(x = Fq2(Fq(g2_X0),Fq(g2_X1)), y = Fq2(Fq(g2_Y0),Fq(g2_Y1))),
    exp_t_minus_one=exp_t_minus_one
)

bls12_381_bilinear_pairing = BilinearPairing(
    bilinear_pairing_curve = bls12_381,
    miller_output_type = Fq12,
    easy_exponentiation = easy_exponentiation,
    hard_exponentation = hard_exponentiation
)