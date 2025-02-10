from elliptic_curves.models.types import G1Point, G2Point
from elliptic_curves.models.bilinear_pairings import BilinearPairingCurve
from elliptic_curves.data_structures.zkscript import ZkScriptVerifyingKey


class PreparedVerifyingKey:
    def __init__(
        self,
        vk,
        curve: BilinearPairingCurve,
        alpha_beta,
        minus_gamma,
        minus_delta,
        gradients_minus_gamma,
        gradients_minus_delta,
    ):
        self.vk = vk
        self.curve = curve
        self.alpha_beta = alpha_beta
        self.minus_gamma = minus_gamma
        self.minus_delta = minus_delta
        self.gradients_minus_gamma = gradients_minus_gamma
        self.gradients_minus_delta = gradients_minus_delta
        return


class VerifyingKey:
    def __init__(
        self,
        curve: BilinearPairingCurve,
        alpha: G1Point,
        beta: G2Point,
        gamma: G2Point,
        delta: G2Point,
        gamma_abc: list[G1Point],
    ):
        self.curve = curve
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        self.delta = delta
        self.gamma_abc = gamma_abc
        return

    def prepare(self):
        alpha_beta = self.curve.pairing([self.alpha], [self.beta])
        minus_gamma = -self.gamma
        minus_delta = -self.delta
        gradients_minus_gamma = minus_gamma.gradients(
            self.curve.miller_loop_engine.exp_miller_loop
        )
        gradients_minus_delta = minus_delta.gradients(
            self.curve.miller_loop_engine.exp_miller_loop
        )

        return PreparedVerifyingKey(
            self,
            self.curve,
            alpha_beta,
            minus_gamma,
            minus_delta,
            gradients_minus_gamma,
            gradients_minus_delta,
        )

    def prepare_for_zkscript(
        self, prepared_vk: PreparedVerifyingKey | None = None
    ) -> ZkScriptVerifyingKey:
        prepared_vk = prepared_vk if prepared_vk is not None else self.prepare()

        return ZkScriptVerifyingKey(
            prepared_vk.alpha_beta.to_list(),
            prepared_vk.minus_gamma.to_list(),
            prepared_vk.minus_delta.to_list(),
            [point.to_list() for point in self.gamma_abc],
            [
                list(map(lambda s: s.to_list(), gradient))
                for gradient in prepared_vk.gradients_minus_gamma
            ],
            [
                list(map(lambda s: s.to_list(), gradient))
                for gradient in prepared_vk.gradients_minus_delta
            ],
        )


class VerifyingKeyGeneric:
    def __init__(self, curve: BilinearPairingCurve):
        self.curve = curve
        return

    def __call__(
        self,
        alpha: G1Point,
        beta: G2Point,
        gamma: G2Point,
        delta: G2Point,
        gamma_abc: list[G1Point],
    ):
        return VerifyingKey(self.curve, alpha, beta, gamma, delta, gamma_abc)

    def deserialise(self, serialised: list[bytes]) -> VerifyingKey:
        """Deserialise a verifying key.

        This function is based on the deserialisation of VK in arkworks. [https://github.com/arkworks-rs/groth16/blob/master/src/data_structures.rs#L32]

        vk is a list of: alpha_g1, beta_g2, gamma_g2, delta_g2, gamma_abc_g1, and each element is serialised in turn
            alpha_g1 -> element in G1
            beta_g2, gamma_g2, delta_g2 -> elements in G2
            gamma_abc_g1 -> list of elements in G1 (as it is a vector, is prepended with the length of the list, encoded as an 8-byte little-endian number )
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
        alpha = self.curve.g1_curve.deserialise(
            serialised[: index + 2 * length_g1], self.curve.g1_field
        )
        index += 2 * length_g1
        beta = self.curve.g2_curve.deserialise(
            serialised[index : index + 2 * length_g2], self.curve.g2_field
        )
        index += 2 * length_g2
        gamma = self.curve.g2_curve.deserialise(
            serialised[index : index + 2 * length_g2], self.curve.g2_field
        )
        index += 2 * length_g2
        delta = self.curve.g2_curve.deserialise(
            serialised[index : index + 2 * length_g2], self.curve.g2_field
        )
        index += 2 * length_g2

        # Check correct length of gamma_abc
        n_abc = int.from_bytes(
            bytes=bytearray(serialised[index : index + 8]), byteorder="little"
        )
        index += 8

        gamma_abc = []
        for _ in range(n_abc):
            gamma_abc.append(
                self.curve.g1_curve.deserialise(
                    serialised[index : index + 2 * length_g1], self.curve.g1_field
                )
            )
            index += 2 * length_g1

        assert index == len(serialised)
        return VerifyingKey(self.curve, alpha, beta, gamma, delta, gamma_abc)
