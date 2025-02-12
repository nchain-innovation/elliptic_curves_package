"""Utilities to facilitate interaction between this library and zkscript."""

from elliptic_curves.data_structures.zkscript import (
    MultiAdditionGradients,
    UnrolledMultiplicationGradients,
    MsmWithFixedBasesGradients,
)
from elliptic_curves.models.ec import ShortWeierstrassEllipticCurvePoint


def multi_addition_gradients(
    points: list[ShortWeierstrassEllipticCurvePoint],
) -> MultiAdditionGradients:
    r"""Compute the gradients required to compute points[0] + .. + points[-1].

    The gradients are computed in the following order: if `out = multi_addition_gradients(points)`, then
    `out[i]` is the gradient required to compute the sum `points[i+1] + \sum_(j=0)^(i) points[i]`.

    Args:
        points (list[ShortWeierstrassEllipticCurvePoint]): The list of points to compute the gradients for.
    """
    running_sum = points[0]
    gradients = []
    for point in points[1:]:
        try:
            gradients.append(running_sum.gradient(point))
        except Exception as _:
            gradients.append(None)
        running_sum += point

    return MultiAdditionGradients(gradients)


def unrolled_multiplication_gradients(
    scalar: int,
    base: ShortWeierstrassEllipticCurvePoint,
    scalar_binary_expansion: list[int] | None = None,
) -> UnrolledMultiplicationGradients:
    """Compute the gradients required to compute `scalar * base`.

    The function outputs `[]` if `scalar == 0`.

    Args:
        scalar (int): The scalar we want to multiply base by.
        base (ShortWeierstrassEllipticCurvePoint): The base of the multiplication.
        scalar_binary_expansion (list[int]): The signed binary expansion of scalar.
    """
    if scalar == 0:
        return UnrolledMultiplicationGradients(None)

    scalar_binary_expansion = (
        scalar_binary_expansion
        if scalar_binary_expansion is not None
        else [int(bit) if scalar > 0 else -int(-bit) for bit in bin(abs(scalar))[2:]][
            ::-1
        ]
    )

    return UnrolledMultiplicationGradients(base.gradients(scalar_binary_expansion))


def multi_scalar_multiplication_with_fixed_bases_gradients(
    scalars: list[int], bases: list[ShortWeierstrassEllipticCurvePoint]
) -> MsmWithFixedBasesGradients:
    r"""Compute the gradients required to compute the msm \sum_i scalars[i] * bases[i].

    Args:
        scalars (list[int]): `scalars[i]` is the scalar by which we multiply `bases[i]`.
        bases (list[ShortWeierstrassEllipticCurvePoint]): `bases[i]` is the i-th base of the msm.

    Returns:
        An instance of MsmWithFixedBasesGradients, where
            - `MsmWithFixedBasesGradients.gradients_additions[i]` is the gradient required to compute the sum
                    `bases[n-i-2] + (\sum_(j=n-i-1)^(n-1) bases[j]`
            - `MsmWithFixedBasesGradients.gradients_multiplications[i]` are the gradients required to compute the multiplication
                    `scalars[i] * bases[i]`
    """
    gradients_multiplications = [
        unrolled_multiplication_gradients(scalar, base)
        for (scalar, base) in zip(scalars, bases)
    ]

    gradients_additions = multi_addition_gradients(
        points=[
            base.multiply(scalar)
            for (scalar, base) in zip(reversed(scalars), reversed(bases))
        ]
    )

    return MsmWithFixedBasesGradients(gradients_multiplications, gradients_additions)
