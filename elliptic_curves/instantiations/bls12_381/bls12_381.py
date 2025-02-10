from elliptic_curves.fields.prime_field import PrimeField
from elliptic_curves.fields.quadratic_extension import QuadraticExtension
from elliptic_curves.fields.cubic_extension import CubicExtension

from elliptic_curves.models.ec import (
    ShortWeierstrassEllipticCurvePoint,
    ShortWeierstrassEllipticCurve,
)
from elliptic_curves.models.miller_loop_engine import MillerLoopEngine
from elliptic_curves.models.types import (
    G1Point,
    G2Point,
    MillerLoopCurve,
    DenominatorElimination,
)
from elliptic_curves.models.pairing_engine import PairingEngine
from elliptic_curves.models.bilinear_pairings import BilinearPairingCurve
from elliptic_curves.data_structures.vk import VerifyingKeyGeneric
from elliptic_curves.data_structures.proof import ProofGeneric

from elliptic_curves.instantiations.bls12_381.parameters import (
    a,
    b,
    q,
    r,
    NON_RESIDUE_FQ,
    NON_RESIDUE_FQ2,
    NON_RESIDUE_FQ6,
    g1_X,
    g1_Y,
    g2_X0,
    g2_X1,
    g2_Y0,
    g2_Y1,
    exp_miller_loop,
    u,
    h1,
    h2,
)

# Field instantiation
Fq = PrimeField(q)
NON_RESIDUE_FQ = Fq.from_list(NON_RESIDUE_FQ)
Fq2 = QuadraticExtension(base_field=Fq, non_residue=NON_RESIDUE_FQ)
NON_RESIDUE_FQ2 = Fq2.from_list(NON_RESIDUE_FQ2)
Fq6 = CubicExtension(base_field=Fq2, non_residue=NON_RESIDUE_FQ2)
NON_RESIDUE_FQ6 = Fq6.from_list(NON_RESIDUE_FQ6)
Fq12 = QuadraticExtension(base_field=Fq6, non_residue=NON_RESIDUE_FQ6)

# Scalar field of the curve
Fr = PrimeField(r)

# Curves
bls12_381 = ShortWeierstrassEllipticCurve(Fq(a), Fq(b))
g1_generator = bls12_381(Fq(g1_X), Fq(g1_Y), False)
g1_curve = bls12_381.set_generator(g1_generator, h1, Fr)

bls12_381_twisted = ShortWeierstrassEllipticCurve(Fq2.zero(), Fq(b) * NON_RESIDUE_FQ2)
g2_generator = bls12_381_twisted(
    Fq2(Fq(g2_X0), Fq(g2_X1)), Fq2(Fq(g2_Y0), Fq(g2_Y1)), False
)
g2_curve = bls12_381_twisted.set_generator(g2_generator, h2, Fr)


# Morphisms
def twisting_morphism(g1_point: G1Point) -> G2Point:
    """
    Untwisting morphism Psi : E' --> E, (x',y') --> (x'/omega^2, y'/omega^3)
    Here omega in F_q^k such that the equation of the twist is y^2 = x^3 + a omega^4 x + b omega^6
    For BLS12, this means: y^2 = x^3 + b omega^6.

    We are in the M-twist case, with equation: y^2 = x^3 + b * xi, hence omega^6 = xi.
    F_q^12 = F_q^6[w] / (w^2 - v) = F_q^2[w,v] / (w^2 - v, v^3 - xi)
    Hence, omega^6 = w, and the untwisting morphism is Psi(x',y') = (x' / w^2, y' / w^3) in E(F_q^12)
    Thus, the twisting morphism is Phi : E --> E' : (x,y) --> (x * w^2, y * w^3) in E'(F_q^12)
    """
    assert not g1_point.is_infinity()
    return ShortWeierstrassEllipticCurvePoint(
        g2_curve, g1_point.x * Fq12.u().power(2), g1_point.y * Fq12.u().power(3), False
    )


def untwisting_morphism(g2_point: G2Point) -> G1Point:
    """Ref. to_twisted_curve => omega = w"""
    assert not g2_point.is_infinity()
    return ShortWeierstrassEllipticCurvePoint(
        g1_curve,
        g2_point.x * Fq12.u().power(-2),
        g2_point.y * Fq12.u().power(-3),
        False,
    )


# Engines
bls_miller_engine = MillerLoopEngine(
    g1_curve,
    g2_curve,
    twisting_morphism,
    untwisting_morphism,
    u,
    exp_miller_loop,
    MillerLoopCurve.G2_CURVE,
    Fq12,
    DenominatorElimination.QUADRATIC,
)


class Bls12PairingEngine(PairingEngine):
    def easy_exponentiation(self, element: QuadraticExtension):
        """Easy exponentation for BLS12_381: f -> f^{(q^6-1)(q^2+1)}."""

        a = element.invert() * element.conjugate()
        b = a.frobenius(n=2)

        return a * b

    def hard_exponentiation(self, element: QuadraticExtension):
        """Hard exponentation for BLS12_381."""

        t0 = element.power(2)
        t1 = t0.power(u)
        t2 = t1.power(u // 2)
        t3 = element.conjugate()
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
        t2 = t2 * element
        t1 = t1 * t2
        t2 = t3.frobenius(n=1)

        return t1 * t2


bls_pairing_engine = Bls12PairingEngine(bls_miller_engine)

# Bilinear pairing curve
BLS12_381 = BilinearPairingCurve(g1_curve, g2_curve, bls_pairing_engine)

# VerifyingKey
VerifyingKeyBls12381 = VerifyingKeyGeneric(BLS12_381)

# Proof
ProofBls12381 = ProofGeneric(BLS12_381)
