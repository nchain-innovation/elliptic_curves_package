from elliptic_curves.models.types import G1Point, G2Point
from elliptic_curves.models.bilinear_pairings import BilinearPairingCurve
from elliptic_curves.data_structures.vk import PreparedVerifyingKey
from elliptic_curves.data_structures.zkscript import ZkScriptProof


class PreparedProof:
    def __init__(
        self,
        proof,
        curve: BilinearPairingCurve,
        public_statements: list[int],
        gradients_b,
        gradients_minus_gamma,
        gradients_minus_delta,
        inverse_miller_loop,
        gradients_msm,
        gradients_public_statements,
    ):
        self.proof = proof
        self.curve = curve
        self.public_statements = public_statements
        self.gradients_b = gradients_b
        self.gradients_minus_gamma = gradients_minus_gamma
        self.gradients_minus_delta = gradients_minus_delta
        self.inverse_miller_loop = inverse_miller_loop
        self.gradients_msm = gradients_msm
        self.gradients_public_statements = gradients_public_statements
        return


class Proof:
    def __init__(self, curve: BilinearPairingCurve, a: G1Point, b: G2Point, c: G1Point):
        self.curve = curve
        self.a = a
        self.b = b
        self.c = c

    def prepare(
        self, prepared_vk: PreparedVerifyingKey, public_statements: list[int]
    ) -> PreparedProof:
        public_statements_extended = [1, *public_statements]

        # Compute \sum_(i=0)^l a_i * gamma_abc[i]
        n_pub = len(public_statements)
        assert (
            len(prepared_vk.vk.gamma_abc) == n_pub + 1
        ), "Wrong number of public inputs"

        sum_gamma_abc = prepared_vk.vk.gamma_abc[0]
        for i in range(1, n_pub + 1):
            sum_gamma_abc += prepared_vk.vk.gamma_abc[i].multiply(
                public_statements_extended[i]
            )

        # Gradients for the pairing
        gradients_b = self.b.gradients(self.curve.miller_loop_engine.exp_miller_loop)

        # Inverse of the Miller loop output
        inverse_miller_loop = self.curve.miller_loop(
            [self.a, sum_gamma_abc, self.c],
            [self.b, prepared_vk.minus_gamma, prepared_vk.minus_delta],
        ).invert()

        # Compute gradients for partial sums: gradients between a_i * gamma_abc[i] and \sum_(j=0)^(i-1) a_j * gamma_abc[j]
        gradients_msm = []
        for i in range(n_pub, 0, -1):
            sum_gamma_abc -= prepared_vk.vk.gamma_abc[i].multiply(
                public_statements_extended[i]
            )
            if (
                sum_gamma_abc.is_infinity()
                or prepared_vk.vk.gamma_abc[i]
                .multiply(public_statements_extended[i])
                .is_infinity()
            ):
                gradients_msm.append([])
            else:
                gradient = sum_gamma_abc.gradient(
                    prepared_vk.vk.gamma_abc[i].multiply(public_statements_extended[i])
                )
                gradients_msm.append(gradient)

        # Gradients for multiplications pub[i] * gamma_abc[i]
        gradients_public_statements = []
        for i in range(1, n_pub + 1):
            if public_statements_extended[i] == 0:
                gradients_public_statements.append([])
            else:
                # Binary expansion of pub[i]
                exp_pub_i = [
                    int(bin(public_statements_extended[i])[j])
                    for j in range(2, len(bin(public_statements_extended[i])))
                ][::-1]

                gradients_public_statements.append(
                    prepared_vk.vk.gamma_abc[i].gradients(exp_pub_i)
                )

        return PreparedProof(
            self,
            self.curve,
            public_statements,
            gradients_b,
            prepared_vk.gradients_minus_delta,
            prepared_vk.gradients_minus_delta,
            inverse_miller_loop,
            gradients_msm,
            gradients_public_statements,
        )

    def prepare_for_zkscript(
        self,
        prepared_vk: PreparedVerifyingKey,
        public_statements: list[int],
        prepared_proof: PreparedProof | None = None,
    ) -> ZkScriptProof:
        prepared_proof = (
            prepared_proof
            if prepared_proof is not None
            else self.prepare(prepared_vk, public_statements)
        )

        gradients_msm = []
        for gradient in prepared_proof.gradients_msm:
            try:
                gradients_msm.append(gradient.to_list())
            except Exception as _:
                gradients_msm.append([])

        gradients_public_statements = []
        for gradients in prepared_proof.gradients_public_statements:
            try:
                gradients_public_statements.append(
                    [
                        list(map(lambda s: s.to_list(), gradient))
                        for gradient in gradients
                    ]
                )
            except Exception as _:
                gradients_public_statements.append([])

        return ZkScriptProof(
            self.a.to_list(),
            self.b.to_list(),
            self.c.to_list(),
            public_statements,
            [
                list(map(lambda s: s.to_list(), gradient))
                for gradient in prepared_proof.gradients_b
            ],
            [
                list(map(lambda s: s.to_list(), gradient))
                for gradient in prepared_vk.gradients_minus_gamma
            ],
            [
                list(map(lambda s: s.to_list(), gradient))
                for gradient in prepared_vk.gradients_minus_delta
            ],
            prepared_proof.inverse_miller_loop.to_list(),
            gradients_msm,
            gradients_public_statements,
        )


class ProofGeneric:
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
