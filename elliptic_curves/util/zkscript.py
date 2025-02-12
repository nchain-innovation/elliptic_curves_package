"""Utilities to facilitate interaction between this library and zkscript."""

from elliptic_curves.models.ec import ShortWeierstrassEllipticCurvePoint

def multi_addition_gradients(points: list[ShortWeierstrassEllipticCurvePoint]) -> list[list[int]]:
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
            gradients.append(running_sum.gradient(point).to_list())
        except Exception as _:
            gradients.append([])
        running_sum += point
    
    return gradients

def unrolled_multiplication_gradients(scalar: int, base: ShortWeierstrassEllipticCurvePoint) -> list[list[list[int]]] | None:
    """Compute the gradients required to compute `scalar * base`.
    
    The function outputs `None` if `scalar == 0`.

    Args:
        scalar (int): The scalar by which we multiply base.
        base (ShortWeierstrassEllipticCurvePoint): The base of the multiplication.
    """
    assert scalar > 0, f"The scalar must be positive, scalar: {scalar}"

    scalar_binary_expansion = [int(bin(scalar)[j]) for j in range(2, len(bin(scalar)))][::-1]
    gradients = [[s.to_list() for s in el] for el in base.gradients(scalar_binary_expansion)] if scalar else None

    return gradients

def multi_scalar_multiplication_with_fixed_bases_gradients(
    scalars: list[int], bases: list[ShortWeierstrassEllipticCurvePoint]
) -> tuple[list[list[int]], list[list[list[list[int]]]]]:
    r"""Compute the gradients required to compute the msm \sum_i scalars[i] * bases[i].
    
    Args:
        scalars (list[int]): `scalars[i]` is the scalar by which we multiply `bases[i]`.
        bases (list[ShortWeierstrassEllipticCurvePoint]): `bases[i]` is the i-th base of the msm.
    
    Returns:
        A tuple `(gradients_additions, gradients_multiplications)` where:
            - `gradients_additions[i]` is the gradient required to compute the sum
                    `bases[n-i-2] + (\sum_(j=n-i-1)^(n-1) bases[j]`
            - `gradients_multiplications[i]` are the gradients required to compute the multiplication
                    `scalars[i] * bases[i]`
    """
    gradients_multiplications = [
        unrolled_multiplication_gradients(
            scalar,
            base
        ) for (scalar, base) in zip(scalars, bases)
    ]

    gradients_additions = multi_addition_gradients(points = [
        base.multiply(scalar) for (scalar, base) in zip(reversed(scalars), reversed(bases))
        ]
    )

    return gradients_additions, gradients_multiplications