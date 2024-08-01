import os, sys
sys.path.append(os.path.normpath(os.path.join(os.path.dirname(__file__),'../../')))

from fields.fq import base_field_from_modulus
from fields.quadratic_extension import quadratic_extension_from_base_field_and_non_residue

from models.ec import elliptic_curve_from_curve
from models.curve import Curve, BilinearPairingCurve
from models.bilinear_pairings import BilinearPairing

from instantiations.mnt4_753.parameters import *
from instantiations.mnt4_753.final_exponentiation import easy_exponentiation, hard_exponentiation

# Field instantiation
Fq = base_field_from_modulus(q=q)
NON_RESIDUE_FQ = Fq.from_list(NON_RESIDUE_FQ)
Fq2 = quadratic_extension_from_base_field_and_non_residue(base_field=Fq,non_residue=NON_RESIDUE_FQ)
NON_RESIDUE_FQ2 = Fq2.from_list(NON_RESIDUE_FQ2)
Fq4 = quadratic_extension_from_base_field_and_non_residue(base_field=Fq2,non_residue=NON_RESIDUE_FQ2)

# Define mul_by_line_eval method for Fq4
# We are not customising line_evaluations, so they will be elements in the miller_ouput_type
def mul_by_line_eval(self,line_eval):
    return self * line_eval

Fq4.mul_by_line_eval = mul_by_line_eval

# Curves
mnt4_753_curve = Curve(a = Fq(a), b = Fq(b))
mnt4_753_twisted_curve = Curve(a = Fq2(NON_RESIDUE_FQ.scalar_mul(a),Fq.zero()), b = Fq2(Fq.zero(),NON_RESIDUE_FQ.scalar_mul(b)))

# Curve classes
MNT4_753, _ = elliptic_curve_from_curve(curve=mnt4_753_curve)
MNT4_753_Twist, _ = elliptic_curve_from_curve(curve=mnt4_753_twisted_curve)

# Twisting morphisms
def to_twisted_curve(self):
    '''
    The untwisting morphism Psi : E' --> E, (x',y') --> (x'/omega^2, y'/omega^3)

    The general equation of a twist is: y^2 = x^3 + a omega^4 x + b omega^6.

    The equation of the twist of MNT4_753 is: y^2 = x^3 + a * 13 * x + b * 13 * u, hence omega^4 = 13.
    
    F_q^4 = F_q[u,r] / (r^2 - u, u^2 - 13) => omega = r and the twisting morphism is Phi : E --> E' : (x,y) --> (x * r^2, y * r^3) in E'(F_q^12)
    '''

    return MNT4_753_Twist(self.x * Fq4.u().power(2), self.y * Fq4.u().power(3))

def to_base_curve(self):
    '''
    Ref. to_twist() => omega = r
    '''

    return MNT4_753(self.x * Fq4.u().power(-2), self.y * Fq4.u().power(-3))

MNT4_753.to_twisted_curve = to_twisted_curve
MNT4_753_Twist.to_base_curve = to_base_curve

# BilinearPairing
mnt4_753 = BilinearPairingCurve(
    q = q,
    r = r,
    h1 = h1,
    h2 = h2,
    curve = mnt4_753_curve,
    twisted_curve = mnt4_753_twisted_curve,
    g1 = MNT4_753(x = Fq(g1_X), y = Fq(g1_Y)),
    g2 = MNT4_753_Twist(x = Fq2(Fq(g2_X0),Fq(g2_X1)), y = Fq2(Fq(g2_Y0),Fq(g2_Y1))),
    exp_t_minus_one=exp_t_minus_one
)

mnt4_753_bilinear_pairing = BilinearPairing(
    bilinear_pairing_curve = mnt4_753,
    miller_output_type = Fq4,
    easy_exponentiation = easy_exponentiation,
    hard_exponentation = hard_exponentiation
)