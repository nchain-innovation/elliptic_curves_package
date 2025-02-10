from typing import Callable

from elliptic_curves.models.types import (
    G1,
    G2,
    G1Point,
    G2Point,
    MillerLoopCurve,
    DenominatorElimination,
)


class MillerLoopEngine:
    def __init__(
        self,
        g1_curve: G1,
        g2_curve: G2,
        twisting_morphism: Callable[[G1Point], G2Point],
        untwisting_morphism: Callable[[G2Point], G1Point],
        val_miller_loop: int,
        exp_miller_loop: list[int],
        miller_loop_curve: MillerLoopCurve,
        target_group,
        denominator_elimination: DenominatorElimination = DenominatorElimination.NONE,
    ):
        if exp_miller_loop[-1] == 0:
            raise ValueError(
                "The most significant element of exp_miller_loop must be non-zero"
            )

        self.g1_curve = g1_curve  # G1
        self.g2_curve = g2_curve  # G2
        self.twisting_morphism = twisting_morphism
        self.untwisting_morphism = untwisting_morphism
        self.val_miller_loop = val_miller_loop
        self.exp_miller_loop = (
            exp_miller_loop  # Signed binary decomposition of `val_miller_loop`
        )
        self.miller_loop_curve = miller_loop_curve
        self.target_group = target_group  # GT
        self.denominator_elimination = denominator_elimination

    def miller_loop(self, g1_points: list[G1Point], g2_points: list[G2Point]):
        r"""
        Compute the \prod miller_loop(P[i],Q[i]).

        This implementation is not optimised for the specific curve.
        """
        assert len(g1_points) == len(g2_points)
        n_points = len(g1_points)

        out = self.target_group.identity()
        prepared_g1_points = (
            g1_points
            if self.miller_loop_curve is MillerLoopCurve.G1_CURVE
            else [self.twisting_morphism(point) for point in g1_points]
        )
        prepared_g2_points = (
            g2_points
            if self.miller_loop_curve is MillerLoopCurve.G2_CURVE
            else [self.untwisting_morphism(point) for point in g2_points]
        )
        running_g2_points = [
            point.copy_with_same_curve() if self.exp_miller_loop[-1] == 1 else -point
            for point in prepared_g2_points
        ]

        for i in range(len(self.exp_miller_loop) - 2, -1, -1):
            out = out.power(2)

            line_eval_update = [
                running_g2_points[j].line_evaluation(
                    running_g2_points[j], prepared_g1_points[j]
                )
                for j in range(n_points)
            ]
            running_g2_points = [
                running_g2_points[j].multiply(2) for j in range(n_points)
            ]

            match self.denominator_elimination:
                case DenominatorElimination.QUADRATIC:
                    pass
                case DenominatorElimination.CUBIC:
                    raise ValueError("To do!")
                case DenominatorElimination.NONE:
                    line_eval_update = [
                        line_eval_update[j]
                        * running_g2_points[j]
                        .line_evaluation(-running_g2_points[j], prepared_g1_points[j])
                        .invert()
                        for j in range(n_points)
                    ]
                case _:
                    raise ValueError("Undefined denominator elimination")

            for j in range(n_points):
                out = out * line_eval_update[j]

            if self.exp_miller_loop[i]:
                line_eval_update = [
                    running_g2_points[j].line_evaluation(
                        prepared_g2_points[j].multiply(self.exp_miller_loop[i]),
                        prepared_g1_points[j],
                    )
                    for j in range(n_points)
                ]
                running_g2_points = [
                    running_g2_points[j]
                    + prepared_g2_points[j].multiply(self.exp_miller_loop[i])
                    for j in range(n_points)
                ]

                match self.denominator_elimination:
                    case DenominatorElimination.QUADRATIC:
                        pass
                    case DenominatorElimination.CUBIC:
                        raise ValueError("To do!")
                    case DenominatorElimination.NONE:
                        line_eval_update = [
                            line_eval_update[j]
                            * running_g2_points[j]
                            .line_evaluation(
                                -running_g2_points[j], prepared_g1_points[j]
                            )
                            .invert()
                            for j in range(n_points)
                        ]
                    case _:
                        raise ValueError("Undefined denominator elimination")

                for j in range(n_points):
                    out = out * line_eval_update[j]

        return out
