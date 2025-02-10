from elliptic_curves.fields.prime_field import PrimeField
from elliptic_curves.fields.quadratic_extension import QuadraticExtension

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


from elliptic_curves.instantiations.mnt4_753.parameters import (
    a,
    b,
    q,
    r,
    NON_RESIDUE_FQ,
    NON_RESIDUE_FQ2,
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
Fq4 = QuadraticExtension(base_field=Fq2, non_residue=NON_RESIDUE_FQ2)

# Scalar field of the curve
Fr = PrimeField(r)

# Curves
mnt4_753 = ShortWeierstrassEllipticCurve(Fq(a), Fq(b))
g1_generator = mnt4_753(Fq(g1_X), Fq(g1_Y), False)
g1_curve = mnt4_753.set_generator(g1_generator, h1, Fr)

mnt4_753_twisted = ShortWeierstrassEllipticCurve(
    Fq2(NON_RESIDUE_FQ.scalar_mul(a), Fq.zero()),
    Fq2(Fq.zero(), NON_RESIDUE_FQ.scalar_mul(b)),
)
g2_generator = mnt4_753_twisted(
    Fq2(Fq(g2_X0), Fq(g2_X1)), Fq2(Fq(g2_Y0), Fq(g2_Y1)), False
)
g2_curve = mnt4_753_twisted.set_generator(g2_generator, h2, Fr)


# Twisting morphisms
def twisting_morphism(g1_point: G1Point) -> G2Point:
    """
    The untwisting morphism Psi : E' --> E, (x',y') --> (x'/omega^2, y'/omega^3)

    The general equation of a twist is: y^2 = x^3 + a omega^4 x + b omega^6.

    The equation of the twist of MNT4_753 is: y^2 = x^3 + a * 13 * x + b * 13 * u, hence omega^4 = 13.

    F_q^4 = F_q[u,r] / (r^2 - u, u^2 - 13) => omega = r and the twisting morphism is Phi : E --> E' : (x,y) --> (x * r^2, y * r^3) in E'(F_q^12)
    """
    assert not g1_point.is_infinity()
    return ShortWeierstrassEllipticCurvePoint(
        g2_curve, g1_point.x * Fq4.u().power(2), g1_point.y * Fq4.u().power(3), False
    )


def untwisting_morphism(g2_point: G2Point) -> G1Point:
    """Ref. to_twist() => omega = r."""
    assert not g2_point.is_infinity()
    return ShortWeierstrassEllipticCurvePoint(
        g1_curve, g2_point.x * Fq4.u().power(-2), g2_point.y * Fq4.u().power(-3), False
    )


# Engines
mnt_miller_engine = MillerLoopEngine(
    g1_curve,
    g2_curve,
    twisting_morphism,
    untwisting_morphism,
    u,
    exp_miller_loop,
    MillerLoopCurve.G2_CURVE,
    Fq4,
    DenominatorElimination.QUADRATIC,
)


class Mnt4PairingEngine(PairingEngine):
    def easy_exponentiation(self, element: QuadraticExtension):
        """Easy exponentiation for MNT4_753 is f -> f^{q^2-1}."""

        out = element.frobenius(2) * element.invert()
        return out

    def hard_exponentiation(self, element: QuadraticExtension):
        """Hard exponentiation for MNT4_753 is f -> f^{q + u + 1}."""

        return element.frobenius(1) * element.power(u) * element


mnt_pairing_engine = Mnt4PairingEngine(mnt_miller_engine)

# Bilinear pairing curve
MNT4_753 = BilinearPairingCurve(g1_curve, g2_curve, mnt_pairing_engine)

# VerifyingKey
VerifyingKeyMnt4753 = VerifyingKeyGeneric(MNT4_753)

# Proof
ProofMnt4753 = ProofGeneric(MNT4_753)
