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
    def __init__(self, q, r, t_minus_one, exp_t_minus_one, h1, h2, curve, twisted_curve, g1, g2, miller_output_type, easy_exponentiation, hard_exponentiation):
        self.q = q
        self.r = r
        self.t_minus_one = t_minus_one
        self.exp_t_minus_one = exp_t_minus_one
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
        Deserialise the serialisation of a verifying key. This function is based on the deserialisation of VK in arkworks. [ref]

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
        Function to deserialise a proof. This function is based on arkworks deserialisation of a proof. [ref]

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