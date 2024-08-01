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
        
class BilinearPairingCurve:
    def __init__(self, q, r, h1, h2, curve, twisted_curve, g1, g2, exp_t_minus_one):
        self.q = q
        self.r = r
        self.h1 = h1
        self.h2 = h2
        self.curve = curve
        self.twisted_curve = twisted_curve
        self.g1 = g1
        self.g2 = g2
        self.exp_t_minus_one = exp_t_minus_one

        return
    