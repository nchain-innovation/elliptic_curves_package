from typing import TypeVar
from enum import Enum
from elliptic_curves.models.ec import (
    ShortWeierstrassEllipticCurvePoint,
    ShortWeierstrassEllipticCurveWithGenerator,
)

G1 = TypeVar("G1", bound=ShortWeierstrassEllipticCurveWithGenerator)
G2 = TypeVar("G2", bound=ShortWeierstrassEllipticCurveWithGenerator)
G1Point = TypeVar("G1Point", bound=ShortWeierstrassEllipticCurvePoint)
G2Point = TypeVar("G2Point", bound=ShortWeierstrassEllipticCurvePoint)


class DenominatorElimination(Enum):
    NONE = "None"
    QUADRATIC = "Quadratic"
    CUBIC = "Cubic"


class MillerLoopCurve(Enum):
    G1_CURVE = "Base"
    G2_CURVE = "Twisted"
