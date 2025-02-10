class ZkScriptVerifyingKey:
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
        gradients_msm: list[list[int]],
        gradients_public_statements: list[list[list[list[int]]]],
    ):
        self.a = a
        self.b = b
        self.c = c
        self.public_statements = public_statements
        self.gradients_b = gradients_b
        self.gradients_minus_gamma = gradients_minus_gamma
        self.gradients_minus_delta = gradients_minus_delta
        self.inverse_miller_loop = inverse_miller_loop
        self.gradients_msm = gradients_msm
        self.gradients_public_statements = gradients_public_statements
        return
