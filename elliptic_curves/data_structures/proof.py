from elliptic_curves.models.types import G1Point, G2Point
from elliptic_curves.models.bilinear_pairings import BilinearPairingCurve
from elliptic_curves.data_structures.vk import PreparedVerifyingKey
from elliptic_curves.data_structures.zkscript import ZkScriptProof
from elliptic_curves.util.zkscript import (
    multi_scalar_multiplication_with_fixed_bases_gradients,
    unrolled_multiplication_gradients,
)


class PreparedProof:
    """Prepared proof encapsulating the pre-computed data required to prove that a Groth16 proof is valid.

    Args:
        proof: The proof.
        curve: The curve over which Groth16 is instantiated.
        public_statements (list[int]): The public statements for which the proof has been created (w/o the initial 1).
        gradients_b: The gradients required to compute proof.b * curve.miller_loop_engine.exp_miller_loop.
        gradients_minus_gamma: The gradients required to compute minus_gamma * curve.miller_loop_engine.exp_miller_loop.
        gradients_minus_delta: The gradients required to compute minus_delta * curve.miller_loop_engine.exp_miller_loop.
        inverse_miller_loop: The inverse of miller_loop([proof.a, sum(vk.gamma_abc), proof.c], [proof.b, -vk.gamma, -vk.delta]).
        msm_key: Instance of MsmWithFixedBasesGradients to compute sum_(i=1)^l public_statement[i-1] * vk.gamma_abc[i].
        gradient_gamma_abc_zero: The gradient to compute vk.gamma_abc[0] + sum_(i=1)^l public_statement[i-1] * vk.gamma_abc[i].
    """

    def __init__(
        self,
        proof,
        curve: BilinearPairingCurve,
        public_statements: list[int],
        gradients_b,
        gradients_minus_gamma,
        gradients_minus_delta,
        inverse_miller_loop,
        msm_key,
        gradient_gamma_abc_zero,
    ):
        self.proof = proof
        self.curve = curve
        self.public_statements = public_statements
        self.gradients_b = gradients_b
        self.gradients_minus_gamma = gradients_minus_gamma
        self.gradients_minus_delta = gradients_minus_delta
        self.inverse_miller_loop = inverse_miller_loop
        self.msm_key = msm_key
        self.gradient_gamma_abc_zero = gradient_gamma_abc_zero
        return


class Proof:
    """Class encapsulating a Groth16 proof."""

    def __init__(self, curve: BilinearPairingCurve, a: G1Point, b: G2Point, c: G1Point):
        self.curve = curve
        self.a = a
        self.b = b
        self.c = c

    def prepare(
        self, prepared_vk: PreparedVerifyingKey, public_statements: list[int]
    ) -> PreparedProof:
        """Turn an instance of `Self` into an instance of `PreparedProof`."""
        public_statements_extended = [1, *public_statements]

        # Compute \sum_(i=0)^l a_i * gamma_abc[i]
        n_pub = len(public_statements)
        assert (
            len(prepared_vk.vk.gamma_abc) == n_pub + 1
        ), "Wrong number of public inputs"

        sum_gamma_abc = prepared_vk.vk.gamma_abc[1].multiply(
            public_statements_extended[1]
        )
        for i in range(2, n_pub + 1):
            sum_gamma_abc += prepared_vk.vk.gamma_abc[i].multiply(
                public_statements_extended[i]
            )

        # Gradients for the pairing
        gradients_b = unrolled_multiplication_gradients(
            self.curve.miller_loop_engine.val_miller_loop,
            self.b,
            self.curve.miller_loop_engine.exp_miller_loop,
        )

        # Inverse of the Miller loop output
        inverse_miller_loop = self.curve.miller_loop(
            [self.a, sum_gamma_abc + prepared_vk.vk.gamma_abc[0], self.c],
            [self.b, prepared_vk.minus_gamma, prepared_vk.minus_delta],
        ).invert()

        # Gradient for addition of gamma_abc[0]
        gradient_gamma_abc_zero = sum_gamma_abc.gradient(prepared_vk.vk.gamma_abc[0])

        # Gradients for msm
        msm_key = multi_scalar_multiplication_with_fixed_bases_gradients(
            public_statements,
            prepared_vk.vk.gamma_abc[1:],
        )

        return PreparedProof(
            self,
            self.curve,
            public_statements,
            gradients_b,
            prepared_vk.gradients_minus_gamma,
            prepared_vk.gradients_minus_delta,
            inverse_miller_loop,
            msm_key,
            gradient_gamma_abc_zero,
        )

    def prepare_for_zkscript(
        self,
        prepared_vk: PreparedVerifyingKey,
        public_statements: list[int],
        prepared_proof: PreparedProof | None = None,
    ) -> ZkScriptProof:
        """Turn an instance of `Self` into an instance of `ZkScriptProof`."""
        prepared_proof = (
            prepared_proof
            if prepared_proof is not None
            else self.prepare(prepared_vk, public_statements)
        )

        gradients_multiplications, gradients_additions = (
            prepared_proof.msm_key.as_data()
        )

        return ZkScriptProof(
            self.a.to_list(),
            self.b.to_list(),
            self.c.to_list(),
            public_statements,
            prepared_proof.gradients_b.as_data(),
            prepared_vk.gradients_minus_gamma.as_data(),
            prepared_vk.gradients_minus_delta.as_data(),
            prepared_proof.inverse_miller_loop.to_list(),
            gradients_multiplications,
            gradients_additions,
            prepared_proof.gradient_gamma_abc_zero.to_list(),
        )


class ProofGeneric:
    """Generic proof class encapsulating the curve over which Groth16 proof is instantiated."""

    def __init__(self, curve: BilinearPairingCurve):
        self.curve = curve

    def __call__(self, a, b, c):
        return Proof(self.curve, a, b, c)

    def deserialise(self, serialised: list[bytes]) -> Proof:
        """Function to deserialise a proof.

        This function is based on arkworks unchecked deserialisation of a proof. [https://github.com/arkworks-rs/groth16/blob/master/src/data_structures.rs#L9]

        A proof is formed by: A, B, C, and each element is serialised in turn
            A, C -> elements in G1
            B -> element in G2
        """
        length_g1 = (
            (self.curve.g1_curve.a.get_modulus().bit_length() + 8)
            // 8
            * self.curve.g1_field.get_extension_degree_over_prime_field()
        )
        length_g2 = (
            (self.curve.g2_curve.a.get_modulus().bit_length() + 8)
            // 8
            * self.curve.g2_field.get_extension_degree_over_prime_field()
        )

        index = 0
        a = self.curve.g1_curve.deserialise(
            serialised[: index + 2 * length_g1], self.curve.g1_field
        )
        index += 2 * length_g1
        b = self.curve.g2_curve.deserialise(
            serialised[index : index + 2 * length_g2], self.curve.g2_field
        )
        index += 2 * length_g2
        c = self.curve.g1_curve.deserialise(
            serialised[index : index + 2 * length_g1], self.curve.g1_field
        )

        return Proof(
            self.curve,
            a,
            b,
            c,
        )
