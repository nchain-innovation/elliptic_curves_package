from abc import ABC, abstractmethod
from elliptic_curves.models.miller_loop_engine import MillerLoopEngine
from elliptic_curves.models.types import G1Point, G2Point


class PairingEngine(ABC):
    def __init__(self, miller_loop_engine: MillerLoopEngine):
        self.miller_loop_engine = miller_loop_engine

    @abstractmethod
    def easy_exponentiation(self, element):
        pass

    @abstractmethod
    def hard_exponentiation(self, element):
        pass

    def pairing(self, g1_points: list[G1Point], g2_points: list[G2Point]):
        r"""Compute \prod pairing(g1_points[i], g2_points[i])"""
        assert len(g1_points) == len(g2_points)

        prepared_g1_points = [point for point in g1_points if not point.is_infinity()]
        prepared_g2_points = [point for point in g2_points if not point.is_infinity()]
        miller_loop_output = self.miller_loop_engine.miller_loop(
            prepared_g1_points, prepared_g2_points
        )
        return self.hard_exponentiation(self.easy_exponentiation(miller_loop_output))
