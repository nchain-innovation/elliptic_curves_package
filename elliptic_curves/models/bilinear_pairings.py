from elliptic_curves.models.pairing_engine import PairingEngine
from elliptic_curves.models.types import G1Point, G2Point, G1, G2


class BilinearPairingCurve:
    def __init__(
        self,
        g1_curve: G1,
        g2_curve: G2,
        pairing_engine: PairingEngine,
    ):
        self.g1_curve = g1_curve
        self.g1_field = g1_curve.a.field
        self.g2_curve = g2_curve
        self.g2_field = g2_curve.a.field
        self.scalar_field = g1_curve.scalar_field
        self.pairing_engine = pairing_engine
        self.miller_loop_engine = pairing_engine.miller_loop_engine

        return

    def twisting_morphism(self, g1_point: G1Point) -> G2Point:
        return self.pairing_engine.miller_loop_engine.twisting_morphism(g1_point)

    def untwisting_morphism(self, g2_point: G2Point) -> G1Point:
        return self.pairing_engine.miller_loop_engine.untwisting_morphism(g2_point)

    def miller_loop(self, g1_points: list[G1Point], g2_points: list[G2Point]):
        return self.miller_loop_engine.miller_loop(g1_points, g2_points)

    def pairing(self, g1_points: list[G1Point], g2_points: list[G2Point]):
        return self.pairing_engine.pairing(g1_points, g2_points)

    def get_modulus(self):
        return self.g1_curve.a.get_modulus()

    def get_order_scalar_field(self):
        return self.scalar_field.get_modulus()
