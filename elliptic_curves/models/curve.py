from elliptic_curves.models.bilinear_pairings import BilinearPairing

class Curve:
    '''
    Curve class to keep track of curve parameters
    '''

    def __init__(self, a, b):
        assert(not(a.power(3).scalar_mul(4) - b.power(2).scalar_mul(27)).is_zero())                    # Assert curve is not singular

        self.a = a
        self.b = b
        
        return

    def evaluate_equation(self, x, y, z = None):
        Field = type(x)
        
        if z is not None:
            return y.power(2) * z - x.power(3) - self.a * x * z.power(2) - self.b * z.power(3)
        else:
            return y.power(2) - x.power(3) - self.a * x - self.b * Field.identity()
        
class BilinearPairingCurve(BilinearPairing):
    def __init__(self, q, r, val_miller_loop, exp_miller_loop, h1, h2, curve, twisted_curve, g1, g2, miller_output_type, easy_exponentiation, hard_exponentiation):
        self.q = q
        self.r = r
        # Value for which we compute the Miller loop: a s.t. e(P,Q) requires computing f_{a,Q}(P)
        self.val_miller_loop = val_miller_loop
        # Signed binary expansion of val_miller_loop
        self.exp_miller_loop = exp_miller_loop
        self.h1 = h1
        self.h2 = h2
        self.curve = curve
        self.twisted_curve = twisted_curve
        self.g1 = g1
        self.g2 = g2
        self.miller_output_type = miller_output_type
        self.easy_exponentiation = easy_exponentiation
        self.hard_exponentiation = hard_exponentiation

        return
    
    def deserialise_vk(self, serialised: list[bytes]):
        '''
        Deserialise the serialisation of a verifying key. This function is based on the deserialisation of VK in arkworks. [https://github.com/arkworks-rs/groth16/blob/master/src/data_structures.rs#L32]

        vk is a list of: alpha_g1, beta_g2, gamma_g2, delta_g2, gamma_abc_g1, and each element is serialised in turn
            alpha_g1 -> element in G1
            beta_g2, gamma_g2, delta_g2 -> elements in G2
            gamma_abc_g1 -> list of elements in G1 (as it is a vector, is prepended with the length of the list, encoded as an 8-byte little-endian number )
        '''
        G1 = type(self.g1)
        G2 = type(self.g2)

        field_G1 = type(G1.CURVE.a)
        field_G2 = type(G2.CURVE.a)

        length_G1 = (field_G1.get_modulus().bit_length() + 8) // 8 * field_G1.EXTENSION_DEGREE
        length_G2 = (field_G2.get_modulus().bit_length() + 8) // 8 * field_G2.EXTENSION_DEGREE

        index = 0
        alpha = G1.deserialise(serialised[:index+2*length_G1],field_G1)
        index += 2*length_G1
        beta = G2.deserialise(serialised[index:index+2*length_G2],field_G2)
        index += 2*length_G2
        gamma = G2.deserialise(serialised[index:index+2*length_G2],field_G2)
        index += 2*length_G2
        delta = G2.deserialise(serialised[index:index+2*length_G2],field_G2)
        index += 2*length_G2

        # Check correct length of gamma_abc
        n_abc = int.from_bytes(bytes=bytearray(serialised[index:index+8]),byteorder='little')
        index += 8

        gamma_abc = []
        for i in range(n_abc):
            gamma_abc.append(G1.deserialise(serialised[index:index+2*length_G1],field_G1))
            index += 2*length_G1

        assert(index == len(serialised))
        return {'alpha' : alpha,
                'beta': beta,
                'gamma': gamma,
                'delta': delta,
                'gamma_abc': gamma_abc}

    def deserialise_proof(self, serialised: list[bytes]):
        """
        Function to deserialise a proof. This function is based on arkworks deserialisation of a proof. [https://github.com/arkworks-rs/groth16/blob/master/src/data_structures.rs#L9]

        A proof is formed by: A, B, C, and each element is serialised in turn
            A, C -> elements in G1
            B -> element in G2
        """
        G1 = type(self.g1)
        G2 = type(self.g2)

        field_G1 = type(G1.CURVE.a)
        field_G2 = type(G2.CURVE.a)

        length_G1 = (field_G1.get_modulus().bit_length() + 8) // 8 * field_G1.EXTENSION_DEGREE
        length_G2 = (field_G2.get_modulus().bit_length() + 8) // 8 * field_G2.EXTENSION_DEGREE

        index = 0
        a = G1.deserialise(serialised[index:index+2*length_G1],field_G1)
        index += 2*length_G1
        b = G2.deserialise(serialised[index:index+2*length_G2],field_G2)
        index += 2*length_G2
        c = G1.deserialise(serialised[index:index+2*length_G1],field_G1)
        index += 2*length_G1

        assert(index == len(serialised))

        return {'a': a,
                'b': b,
                'c': c}
    
    def prepare_groth16_proof(self, pub, proof, vk, miller_loop_type, denominator_elimination):
        """
        Take a a list of public statements, a proof and a vk, returns the data needed to generate the unlocking script for the Groth16 Bitcoin Script verifier [https://github.com/nchain-innovation/zkscript_package/blob/main/zkscript/groth16/model/groth16.py#L141]
		
		Miller loop type is either 'base_curve' or 'twisted_curve'
        """
        assert(miller_loop_type in ['base_curve','twisted_curve'])
        assert(denominator_elimination in [None,'quadratic','cubic'])

        exp_miller_loop = self.exp_miller_loop

        gamma = vk['gamma']
        delta = vk['delta']
        gamma_abc = vk['gamma_abc']

        A = proof['a']
        B = proof['b']
        C = proof['c']

        pub_extended = [1] + pub

        # Compute \sum_(i=0)^l a_i * gamma_abc[i]
        n_pub = len(pub_extended) - 1
        sum_gamma_abc = gamma_abc[0]
        for i in range(1,n_pub+1):
            sum_gamma_abc += gamma_abc[i].multiply(pub_extended[i])

        # Lambdas for the pairing
        lambdas_B_exp_miller_loop = [list(map(lambda s: s.to_list(),el)) for el in B.get_lambdas(exp_miller_loop)]
        lambdas_minus_gamma_exp_miller_loop = [list(map(lambda s: s.to_list(),el)) for el in (-gamma).get_lambdas(exp_miller_loop)]
        lambdas_minus_delta_exp_miller_loop = [list(map(lambda s: s.to_list(),el)) for el in (-delta).get_lambdas(exp_miller_loop)]

		# Inverse of the Miller loop output
        match miller_loop_type:
            case 'base_curve':
                inverse_miller_loop = self.triple_miller_loop_on_base_curve(A,sum_gamma_abc,C,B,-gamma,-delta,denominator_elimination).invert().to_list()
            case 'twisted_curve':
                inverse_miller_loop = self.triple_miller_loop_on_twisted_curve(A,sum_gamma_abc,C,B,-gamma,-delta,denominator_elimination).invert().to_list()

        # Compute lamdbas for partial sums: gradients between a_i * gamma_abc[i] and \sum_(j=0)^(i-1) a_j * gamma_abc[j]
        lamdbas_partial_sums = []
        for i in range(n_pub,0,-1):
            sum_gamma_abc -= gamma_abc[i].multiply(pub_extended[i])
            if sum_gamma_abc.is_infinity() or gamma_abc[i].multiply(pub_extended[i]).is_infinity():
                lamdbas_partial_sums.append([])
            else:
                lam = sum_gamma_abc.get_lambda(gamma_abc[i].multiply(pub_extended[i]))
                lamdbas_partial_sums.append(lam.to_list())

		# Lambdas for multiplications pub[i] * gamma_abc[i]
        lambdas_multiplications = []
        for i in range(1,n_pub+1):
            if pub_extended[i] == 0:
                lambdas_multiplications.append([])
            else:
                # Binary expansion of pub[i]
                exp_pub_i = [int(bin(pub_extended[i])[j]) for j in range(2,len(bin(pub_extended[i])))][::-1]

                lambdas_multiplications.append([list(map(lambda s: s.to_list(),el)) for el in gamma_abc[i].get_lambdas(exp_pub_i)])
        

        out = {
            'pub': pub,
            'A': A.to_list(),
            'B': B.to_list(),
            'C': C.to_list(),
            'lambdas_B_exp_miller_loop': lambdas_B_exp_miller_loop,
            'lambdas_minus_gamma_exp_miller_loop': lambdas_minus_gamma_exp_miller_loop,
            'lambdas_minus_delta_exp_miller_loop': lambdas_minus_delta_exp_miller_loop,
            'inverse_miller_loop': inverse_miller_loop,
            'lamdbas_partial_sums': lamdbas_partial_sums,
            'lambdas_multiplications': lambdas_multiplications
        }

        return out
