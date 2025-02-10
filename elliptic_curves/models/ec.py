from typing import Self


def lexicographic_greatest(element: list[int], other: list[int]) -> bool:
    """Return `True` is `element` is the lexicographic largest. Else, returns `False`."""
    is_largest = None
    for el, ot in zip(element.to_list()[::-1], other.to_list()[::-1]):
        if el < ot:
            is_largest = False
            break
        if el > ot:
            is_largest = True
            break
    return is_largest


class ShortWeierstrassEllipticCurve:
    """Short Weierstrass curve in affine form.."""

    def __init__(self, a, b):
        assert a.is_same_field(b)
        assert not (
            a.power(3).scalar_mul(4) - b.power(2).scalar_mul(27)
        ).is_zero()  # Assert curve is not singular

        self.a = a
        self.b = b

        return

    def __call__(self, x, y, infinity):
        return ShortWeierstrassEllipticCurvePoint(self, x, y, infinity)

    def __eq__(element, other):
        result = True
        result &= element.a == other.a
        return result & (element.b == other.b)

    def evaluate_equation(self, x, y):
        return y.power(2) - x.power(3) - self.a * x - self.b

    def infinity(self):
        """Return the point at infinity: (0, 0, True)."""
        return ShortWeierstrassEllipticCurvePoint(
            self, self.a.field.zero(), self.b.field.zero(), True
        )

    def set_generator(self, generator, cofactor: int, scalar_field):
        return ShortWeierstrassEllipticCurveWithGenerator(
            self.a.copy_with_same_field(),
            self.b.copy_with_same_field(),
            generator,
            cofactor,
            scalar_field,
        )

    def from_list(self, coordinates, field) -> list[int]:
        """Reads a the list of coordinates into a ShortWeierstrassEllipticCurvePoint. First the x-coordinate, then the y-coordinate."""
        if coordinates == [None, None]:
            return self.infinity()
        else:
            length = len(coordinates)
            return ShortWeierstrassEllipticCurvePoint(
                self,
                field.from_list(coordinates[: length // 2]),
                field.from_list(coordinates[length // 2 :]),
                False,
            )

    def deserialise(self, serialised: list[bytes], field):
        """Deserialise a list of bytes into a ShortWeierstrassEllipticCurvePoint.

        This function is based on the unchecked deserialisation function for the trait SWCurveConfig of arkworks.
        See [https://github.com/arkworks-rs/algebra/blob/master/ec/src/models/short_weierstrass/mod.rs#L115].

        It works as follows: serialised is a list of ints representing the little-endian encoding of (x,y). The encoding is:
            [LE(x), LE(y)_mod]
        where both elements are of length equal to the byte length of the field over which the curve is defined, and
            LE(y)_mod[:n-1] = LE(y), LE(y)_mod[-1] = LE(y)[-1] | flags
        where flags is the OR of:
            1 << 7 if y > -y (lexicographic order)
            1 << 6 if Point at infinity

        In uncompressed mode, the flag is disregarded if it's not the infinity flag.
        """
        is_infinity = (serialised[-1] >> 6) & 1

        if is_infinity:
            return self.infinity()
        else:
            length = len(serialised)
            x = field.deserialise(serialised[: length // 2])

            serialised_y = serialised[length // 2 :]
            serialised_y[-1] = serialised_y[-1] & ~(1 << 7)  # Remove flag
            y = field.deserialise(serialised_y)

            return ShortWeierstrassEllipticCurvePoint(self, x, y, False)


class ShortWeierstrassEllipticCurveWithGenerator(ShortWeierstrassEllipticCurve):
    def __init__(self, a, b, generator, cofactor: int, scalar_field):
        assert a.is_same_field(b)
        assert not (a.power(3).scalar_mul(4) - b.power(2).scalar_mul(27)).is_zero()
        assert generator.multiply(scalar_field.get_modulus()).is_infinity()

        self.a = a
        self.b = b
        self.generator = generator
        self.cofactor = cofactor
        self.scalar_field = scalar_field

        return

    def generate_random_point_and_multiplier(self):
        multiplier = self.scalar_field.generate_random_point().to_int()
        return self.generator.multiply(multiplier), multiplier

    def generate_random_point(self):
        rnd_point, _ = self.generate_random_point_and_multiplier()
        return rnd_point

    def get_generator(self):
        return self.generator


class ShortWeierstrassEllipticCurvePoint:
    """Point on an affine curve in Short Weierstrass form."""

    def __init__(self, curve, x, y, infinity: bool):
        assert infinity or curve.evaluate_equation(x, y).is_zero()

        self.curve = curve
        self.x = x
        self.y = y
        self.infinity = infinity

        return

    def is_infinity(self) -> bool:
        return self.infinity

    def is_same_curve(self, other):
        return self.curve == other.curve

    def copy_with_same_curve(self):
        return (
            self.curve.infinity()
            if self.is_infinity()
            else ShortWeierstrassEllipticCurvePoint(
                self.curve,
                self.x.copy_with_same_field(),
                self.y.copy_with_same_field(),
                False,
            )
        )

    def gradient(self, other: Self):
        """Compute the gradient of the line through `self` and `other`.

        If `self` == `other`, return the gradient of the tangent line at P.
        """
        assert self.is_same_curve(other)
        assert not self.is_infinity()
        assert not other.is_infinity()

        if self == other:
            return (self.x.power(2).scalar_mul(3) + self.curve.a) * self.y.scalar_mul(
                2
            ).power(-1)
        else:
            return (other.y - self.y) * (other.x - self.x).power(-1)

    def __eq__(element, other: Self):
        result = True
        result &= element.curve == other.curve
        result &= element.x == other.x
        return result & (element.y == other.y)

    def __repr__(self):
        return f"ShortWeierstrassCurve({self.x},{self.y})"

    def __neg__(self):
        negated = self.copy_with_same_curve()
        negated.y = -negated.y
        return negated

    def __add__(element, other: Self):
        assert element.is_same_curve(other)

        if element.is_infinity():
            return other.copy_with_same_curve()
        elif other.is_infinity():
            return element.copy_with_same_curve()
        else:
            if element == -other:
                return element.curve.infinity()
            else:
                gradient = element.gradient(other)
                x_coordinate = gradient.power(2) - element.x - other.x
                y_coordinate = gradient * (element.x - x_coordinate) - element.y

                out = element.curve.infinity()
                out.x = x_coordinate
                out.y = y_coordinate
                out.infinity = False
                return out

    def __sub__(element, other: Self):
        return element + (-other)

    def gradients(self, bits: list[int]):
        r"""
        Computes the gradients of the multiplication: e * Q, where e = \sum bits[i] * 2**i
        gradients[i] is the (list of) gradients(s) computed at the i-th step of the iteration, going down from log(e)-2 to 0:
        gradients[0] = gradient(s) computed when i = log(e)-2.
        If bits[i] != 0, then gradients[i] is a list where the first element is the gradient for the doubling, and the second is the one for the sum/subtraction.
        """

        gradients = []

        if bits[-1] == 1:
            T = self.copy_with_same_curve()
        elif bits[-1] == -1:
            T = -self
        else:
            raise ValueError("The most significant element of `bits` must be non-zero")

        for i in range(len(bits) - 2, -1, -1):
            to_add = []
            to_add.append(T.gradient(T))
            T = T + T

            if bits[i] == 1:
                to_add.append(T.gradient(self))
                T = T + self
            elif bits[i] == -1:
                to_add.append(T.gradient(-self))
                T = T - self
            else:
                pass

            gradients.extend([to_add])

        return gradients

    def multiply(self, n: int):
        if self.is_infinity():
            return self.copy_with_same_curve()
        else:
            if n == 0:
                return self.curve.infinity()
            else:
                val = self.copy_with_same_curve()
                result = self.curve.infinity()

                if n < 0:
                    n = -n
                    val = -val

                while n > 0:
                    if n % 2 == 1:
                        result = result + val
                    val = val + val
                    n = n // 2

        return result

    def line_evaluation(self, Q, P):
        r"""
        Evaluate the line through `self` and `Q` at `P`. If `self` == `Q`, the line is the tangent at `self`. If `self` == `-Q`, the line is the vertical.

        The line is y - self.y = gradient * (x - self.x), where gradient = self.gradient(Q)
        Remark: `self`, `Q` and `P` must not be the point at infinity.
        """
        assert self.is_same_curve(Q)
        assert self.is_same_curve(P)
        assert not self.is_infinity()
        assert not Q.is_infinity()
        assert not P.is_infinity()

        """
        # Handle the case in which self, Q and P live on the same curve, but with coordinates in different extension fields
        if self.field.get_extension_degree_over_prime_field() > max(
            Q.field.get_extension_degree_over_prime_field(),
            P.field.get_extension_degree_over_prime_field()
        ):
            identity = self.field.identity()
        elif Q.field.get_extension_degree_over_prime_field() > P.field.get_extension_degree_over_prime_field():
            identity = Q.field.identity()
        else:
            identity = P.field.identity()
        """

        if self == -Q:
            return P.x - Q.x
        else:
            gradient = self.gradient(Q)
            out = P.y - self.y - gradient * (P.x - self.x)

        return out

    def to_list(self) -> list[int]:
        """Returns the list of coordinates defining self. First the x-coordinate, then the y-coordinate."""
        if self.is_infinity():
            return [None, None]
        else:
            return [*self.x.to_list(), *self.y.to_list()]
