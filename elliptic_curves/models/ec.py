from copy import deepcopy

# The two classes below are not meant to be directly used by the user. They should be exported using the function below.
class EllipticCurve:
    CURVE = None

    def __init__(self, x, y):
        """
        Point on elliptic curve specified by curve class
        """
        Curve = type(self)
        assert((x is None and y is None) or Curve.CURVE.evaluate_equation(x,y).is_zero())   # Model point at infinity as (None,None)

        self.x = x
        self.y = y

        return

    def __eq__(P,Q):
        return P.__dict__ == Q.__dict__

    def __repr__(self):
        return f'({self.x},{self.y})'

    def __neg__(self):
        if self.is_infinity():
            out = deepcopy(self)
        else:
            out = deepcopy(self)
            out.y = -out.y
        
        return out

    def __add__(P,Q):
        assert(type(P) == type(Q))
        Curve = type(P)

        if P.is_infinity():
            out = deepcopy(Q)
        elif Q.is_infinity():
            out = deepcopy(P)
        else:
            if P == -Q:
                out = Curve.point_at_infinity()
            else:
                out = deepcopy(P)

                lambdaCoeff = P.get_lambda(Q)
                out.x = lambdaCoeff.power(2) - P.x - Q.x
                out.y = lambdaCoeff * (P.x - out.x) - P.y

        return out

    def __sub__(P,Q):
        return P + (-Q)

    def get_lambda(self,Q):
        r"""
        Compute the gradient of the line through self and Q.
        If self == Q, return the gradient of the tangent line at P.

        Remark: self and Q must not be the point at infinity.
        """
        assert(type(self) == type(Q))
        Curve = type(self)
        Field = type(self.x)

        if self == Q:
            return (self.x.power(2).scalar_mul(3) + Curve.CURVE.a * Field.identity()) * self.y.scalar_mul(2).power(-1)
        else:
            return (Q.y - self.y) * (Q.x - self.x).power(-1)

    def get_lambdas(self, expU: list[int]):
        r"""
        Computes the lambdas of the multiplication: u * Q, where u = \sum expU[i] * 2**i
        lambdas[i] is the (list of) lambda(s) computed at the i-th step of the iteration, going down from log(u)-2 to 0:
        lambdas[0] = lambda(s) computed when i = log(u)-2.
        If expU[i] != 0, then lambdas[i] is a list where the first element is the lambda for the doubling, and the second is the one for the sum/subtraction.
        """

        lambdas = []

        if expU[-1] == 1:
            T = deepcopy(self)
        elif expU[-1] == -1:
            T = -deepcopy(self)
        else:
            raise ValueError('The most significant element of expE must be non-zero')

        for i in range(len(expU)-2,-1,-1):
            toAdd = []
            toAdd.append(T.get_lambda(T))
            T = T + T

            if expU[i] == 1:
                toAdd.append(T.get_lambda(self))
                T = T + self
            elif expU[i] == -1:
                toAdd.append(T.get_lambda(-self))
                T = T - self
            else:
                pass

            lambdas.extend([toAdd])

        return lambdas

    def point_at_infinity():
        r"""
        We model the point at infinity as (None,None)
        """
        return EllipticCurve(x=None,y=None)

    def is_infinity(self) -> bool:
        return (self.x is None) and (self.y is None)

    def multiply(self, n: int):
        Curve = type(self)

        if self.is_infinity():
            result = deepcopy(self)
        else:
            if n == 0:
                result = Curve.point_at_infinity()
            else:
                val = deepcopy(self)
                result = Curve.point_at_infinity()

                if n < 0:
                    n = -n
                    val = -val
                
                while n > 0:
                    if n % 2 == 1:
                        result = result + val
                    val = val + val
                    n = n // 2

        return result

    def to_projective(self):
        Field = type(self.x)
        
        if self.is_infinity():
            return EllipticCurveProjective.point_at_infinity(Field)
        else:
            return EllipticCurveProjective(
                x=deepcopy(self.x),
                y=deepcopy(self.y),
                z=Field.identity()
                )

    def line_evaluation(self,Q,P):
        r"""
        Evaluate the line through self and Q at P. If self == Q, the line is the tanget at self. If self == -Q, the line is the vertical
        
        The line is y - self.y = lambda * (x - self.x), where lambda = self.getLambda(Q)
        Remark: self, Q and P must not be the point at infinity.
        """
        if self.is_infinity() or Q.is_infinity() or P.is_infinity():
            raise ValueError("Self, Q and P must not be the point at infinity!")
        
        Field_self = type(self.x)
        Field_Q = type(Q.x)
        Field_P = type(P.x)

        # Handle the case in which self, Q and P live on the same curve, but with coordinates in different extension fields
        if Field_self.EXTENSION_DEGREE > max(Field_Q.EXTENSION_DEGREE, Field_P.EXTENSION_DEGREE):
            Field = Field_self
        elif Field_Q.EXTENSION_DEGREE > Field_P.EXTENSION_DEGREE:
            Field = Field_Q
        else:
            Field = Field_P

        if self == -Q:
            out = P.x * Field.identity() - Q.x * Field.identity()
        else:
            lam = self.get_lambda(Q) * Field.identity()
            out = P.y * Field.identity() - self.y * Field.identity() - lam * (P.x  * Field.identity()  - self.x  * Field.identity())

        return out
    
    def deserialise(serialised: list[bytes], field):
        """
        Function that a list of integers and inteprets it as a point on the elliptic curve self and returns its serialisation.
        This function is based on the deserialisation function for the trait SWCurveConfig of arkworks, only uncompressed mode. [ref]
        
        It works as follows: serialised is a list of ints representing the little-endian encoding of (x,y). The encoding is:
            [LE(x), LE(y)_mod]
        where both elements are of length equal to the byte length of the field over which the curve is defined, and
            LE(y)_mod[:n-1] = LE(y), LE(y)_mod[-1] = LE(y)[-1] | flags
        where flags is the OR of:
            1 << 7 if y > -y (lexicographic order)
            1 << 6 if Point at infinity
        """
        is_infinity = (serialised[-1] >> 6) & 1
        is_largest = (serialised[-1] >> 7) & 1
        
        if is_infinity:
                return EllipticCurve.point_at_infinity()
            
        else:        
            serialised_x = serialised[:len(serialised)//2]
            x = field.deserialise(serialised_x)
            serialised_y = serialised[len(serialised)//2:]
            serialised_y[-1] = serialised_y[-1]  & ~(1 << 7)
            y = field.deserialise(serialised_y)

            y_is_largest = False
            for el, minus_el in zip(y.to_list()[::-1],(-y).to_list()[::-1]):
                if el > minus_el:
                    y_is_largest = True
                    break
            
            if (y_is_largest and not is_largest) or (not y_is_largest and is_largest):
                y = -y
        
        return EllipticCurve(x=x,y=y)

    def to_list(self) -> list[int]:
        """
        Returns the list of coordinates defining self. First the x-coordinate, then the y-coordinate
        """
        out = []
        out.extend(self.x.to_list())
        out.extend(self.y.to_list())
        
        return out

class EllipticCurveProjective:
    CURVE = None

    def __init__(self, x, y, z):
        """
        Projective point on the elliptic curve specified by curve class
        """
        Curve = type(self)
        assert(Curve.CURVE.evaluate_equation(x,y,z).is_zero())

        self.x = x
        self.y = y
        self.z = z

        return

    def __eq__(P,Q):
        if P.z.is_zero():
            if Q.z.is_zero():
                return True
            else:
                return False
        else:
            if Q.z.is_zero():
                return False
            else:
                return P.to_affine() == Q.to_affine()

    def __repr__(self):
        return f'[{self.x} : {self.y} : {self.z}]'

    def __neg__(self):
        out = deepcopy(self)
        out.y = -out.y
        return out

    def __add__(P,Q):
        assert(type(P) == type(Q))
        Field = type(P.z)
        Curve = type(P)

        if P.z.is_zero():
            return deepcopy(Q)
        elif Q.z.is_zero():
            return deepcopy(P)
        else:
            sumAff = P.to_affine() + Q.to_affine()
            if P != Q:
                denominator = P.z * Q.z * (Q.x - P.x).power(3)                                      # Order chosen to be consistent with affine sum
            else:
                denominator = P.z * Q.z * (P.x.scalar_mul(2)).power(3)

            if sumAff.is_infinity():                                                                # Points are inverse of one-another
                out = Curve.point_at_infinity(field=Field)
            else:
                out = deepcopy(P)
                out.x = sumAff.x * denominator
                out.y = sumAff.y * denominator
                out.z = denominator

            return out

    def __sub__(P,Q):
        return P + (-Q)

    def point_at_infinity(field):
        return EllipticCurveProjective(field.zero(),field.identity(),field.zero())

    def is_infinity(self) -> bool:
        return (self.x.is_zero()) and (self.y.x == 1) and (self.z.is_zero())

    def multiply(self, n: int):
        Curve = type(self)
        Field = type(self.x)

        if self.is_infinity():
            result = deepcopy(self)
        else:

            if n == 0:
                result = Curve.point_at_infinity(field=Field)
            else:
                val = deepcopy(self)
                result = Curve.point_at_infinity(field=Field)

                if n < 0:
                    n = -n
                    val = -val
                
                while n > 0:
                    if n % 2 == 1:
                        result = result + val
                    val = val + val
                    n = n // 2

        return result

    def to_affine(self):
        if not self.z.is_zero():
            return EllipticCurve(x=self.x * self.z.invert(),y=self.y * self.z.invert())

    def to_list(self) -> list[int]:
        """
        Returns the list of coordinates defining self. First the x-coordinate, then the y-coordinate, then the z-coordinate
        """
        out = []
        out.extend(self.x.to_list())
        out.extend(self.y.to_list())
        out.extend(self.z.to_list())
        
        return out

def elliptic_curve_from_curve(curve):
    """
    Exports EllipticCurve and EllipticCurveProjective for a give curve
    """

    class AffineEllipticCurve(EllipticCurve):
        CURVE = curve

        def point_at_infinity():
            r"""
            We model the point at infinity as (None,None)
            """
            return AffineEllipticCurve(x=None,y=None)
        
        def to_projective(self):
            Field = type(self.x)
            
            if self.is_infinity():
                return ProjectiveEllipticCurve.point_at_infinity(Field)
            else:
                return ProjectiveEllipticCurve(
                    x=deepcopy(self.x),
                    y=deepcopy(self.y),
                    z=Field.identity()
                    )
        
        def deserialise(serialised: list[bytes], field):
            """
            See comments for function above
            """
            is_infinity = (serialised[-1] >> 6) & 1
            is_largest = (serialised[-1] >> 7) & 1
            
            if is_infinity:
                    return AffineEllipticCurve.point_at_infinity()
                
            else:        
                serialised_x = serialised[:len(serialised)//2]
                x = field.deserialise(serialised_x)
                serialised_y = serialised[len(serialised)//2:]
                serialised_y[-1] = serialised_y[-1]  & ~(1 << 7)
                y = field.deserialise(serialised_y)

                y_is_largest = False
                for el, minus_el in zip(y.to_list()[::-1],(-y).to_list()[::-1]):
                    if el > minus_el:
                        y_is_largest = True
                        break
                
                if (y_is_largest and not is_largest) or (not y_is_largest and is_largest):
                    y = -y
            
            return AffineEllipticCurve(x=x,y=y)

    class ProjectiveEllipticCurve(EllipticCurveProjective):
        CURVE = curve
        
        def point_at_infinity(field):
            return ProjectiveEllipticCurve(field.zero(),field.identity(),field.zero())

        def to_affine(self):
            if not self.z.is_zero():
                return AffineEllipticCurve(x=self.x * self.z.invert(),y=self.y * self.z.invert())
    
    return AffineEllipticCurve, ProjectiveEllipticCurve

        








