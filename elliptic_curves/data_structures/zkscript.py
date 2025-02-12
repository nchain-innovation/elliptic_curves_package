from dataclasses import dataclass
from elliptic_curves.fields.cubic_extension import CubicExtensionElement
from elliptic_curves.fields.prime_field import PrimeField
from elliptic_curves.fields.quadratic_extension import QuadraticExtension


class ZkScriptVerifyingKey:
    """Class encapsulating the data required to generate a Groth16 verifier Bitcoin script.

    Args:
        alpha_beta: The value of pairing(vk.alpha, vk.beta).
        minus_gamma: The value of -vk.gamma.
        minus_delta: The value of -vk.delta.
        gradients_minus_gamma: The gradients required to compute minus_gamma * curve.miller_loop_engine.exp_miller_loop
        gradients_minus_delta: The gradients required to compute minus_delta * curve.miller_loop_engine.exp_miller_loop
    """

    def __init__(
        self,
        alpha_beta: list[int],
        minus_gamma: list[int],
        minus_delta: list[int],
        gamma_abc: list[list[int]],
        gradients_minus_gamma: list[list[list[int]]],
        gradients_minus_delta: list[list[list[int]]],
    ):
        self.alpha_beta = alpha_beta
        self.minus_gamma = minus_gamma
        self.minus_delta = minus_delta
        self.gamma_abc = gamma_abc
        self.gradients_minus_gamma = gradients_minus_gamma
        self.gradients_minus_delta = gradients_minus_delta
        return


class ZkScriptProof:
    """Class encapsulating the data required to generate the unlocking script for a Groth16 verifier.

    Args:
        a (list[int]): The value `a` in the proof.
        b (list[int]): The value `b` in the proof.
        c (list[int]): The value `c` in the proof.
        public_statements (list[int]): The public statements for which the proof has been created (w/o the initial 1).
        gradients_b: The gradients required to compute proof.b * curve.miller_loop_engine.exp_miller_loop.
        gradients_minus_gamma: The gradients required to compute minus_gamma * curve.miller_loop_engine.exp_miller_loop.
        gradients_minus_delta: The gradients required to compute minus_delta * curve.miller_loop_engine.exp_miller_loop.
        inverse_miller_loop: The inverse of miller_loop([proof.a, sum(vk.gamma_abc), proof.c], [proof.b, -vk.gamma, -vk.delta]).
        gradients_multiplications: The gradients required to compute the multiplications public_statements[i] * vk.gamma_abc[i+1]
        gradients_additions: The gradients required to compute the sum sum_(i=1)^l public_statements[i-1] * vk.gamma_abc[i]
        gradient_gamma_abc_zero: The gradient to compute vk.gamma_abc[0] + sum_(i=1)^l public_statement[i-1] * vk.gamma_abc[i].
    """

    def __init__(
        self,
        a: list[int],
        b: list[int],
        c: list[int],
        public_statements: list[int],
        gradients_b: list[list[list[int]]],
        gradients_minus_gamma: list[list[list[int]]],
        gradients_minus_delta: list[list[list[int]]],
        inverse_miller_loop: list[int],
        gradients_multiplications: list[list[list[list[int]]]],
        gradients_additions: list[list[int]],
        gradient_gamma_abc_zero: list[int],
    ):
        self.a = a
        self.b = b
        self.c = c
        self.public_statements = public_statements
        self.gradients_b = gradients_b
        self.gradients_minus_gamma = gradients_minus_gamma
        self.gradients_minus_delta = gradients_minus_delta
        self.inverse_miller_loop = inverse_miller_loop
        self.gradients_multiplications = gradients_multiplications
        self.gradients_additions = gradients_additions
        self.gradient_gamma_abc_zero = gradient_gamma_abc_zero
        return


@dataclass
class MultiAdditionGradients:
    """Class encapsulating the gradients required to compute a multi addition of elliptic curve points.

    Args:
        gradients (list[PrimeField | QuadraticExtension | CubicExtensionElement | None]): If an element of the
            list is not `None`, it is the gradient required to compute a step in the sum. Else, it means that
            one of the two points being summed at that step is the point at infinity.
    """

    gradients: list[PrimeField | QuadraticExtension | CubicExtensionElement | None]

    def as_data(self) -> list[list[int]]:
        """Return the gradients in the format required by the multi_addition script in zkscript."""

        return [
            [] if gradient is None else gradient.to_list()
            for gradient in self.gradients
        ]


@dataclass
class UnrolledMultiplicationGradients:
    """Class encapsulating the gradients required to compute a scalar point multiplication.

    Args:
        gradients (list[list[PrimeField | QuadraticExtension | CubicExtensionElement]] | None):
            If `None`, then the scalar we are multiply by is zero. Else, it is a list of the required gradients
            generated with base.gradients(scalar_binary_expansion)
    """

    gradients: (
        list[list[PrimeField | QuadraticExtension | CubicExtensionElement]] | None
    )

    def as_data(self) -> list[list[list[int]]]:
        """Return the gradients in the format required by the unrolled_multiplication script in zkscript."""
        if self.gradients is None:
            return []

        return [[s.to_list() for s in gradients] for gradients in self.gradients]


@dataclass
class MsmWithFixedBasesGradients:
    """Class encapsulating the data required to compute a multi scalar multiplication.

    Args:
        gradients_multiplications (list[UnrolledMultiplicationGradients]): The gradients required by the multiplications bases[i] * scalar[i].
        gradients_additions (MultiAdditionGradients): The gradients required to compute the multi addition sum_(i=0)^n bases[i] * scalar[i].
    """

    gradients_multiplications: list[UnrolledMultiplicationGradients]
    gradients_additions: MultiAdditionGradients

    def as_data(self) -> tuple[list[list[list[list[int]]]], list[list[int]]]:
        """Return the gradients in the format required by the multi_scalar_multiplication_with_fixed_bases script in zkscript."""
        return [
            gradients.as_data() for gradients in self.gradients_multiplications
        ], self.gradients_additions.as_data()
