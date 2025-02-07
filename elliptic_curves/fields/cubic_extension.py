from math import gcd

from elliptic_curves.fields.prime_field import PrimeFieldElement


class CubicExtension:
    """
    Field F, cubic extension of base field B.
    F = B[v] / (v^3 - non_residue)
    """

    EXTENSION_DEGREE_OVER_BASE_FIELD = 3

    def __init__(self, base_field, non_residue: PrimeFieldElement):
        self.base_field = base_field
        self.non_residue = non_residue
        return

    def __call__(
        self, x0: PrimeFieldElement, x1: PrimeFieldElement, x2: PrimeFieldElement
    ):
        return CubicExtensionElement(self, x0, x1, x2)

    def __eq__(field, other):
        result = True
        result &= type(other) is type(field)
        if not result:
            return result
        else:
            result &= field.base_field == other.base_field
            return (
                result
                if not result
                else result & (field.non_residue == other.non_residue)
            )

    def get_modulus(self):
        return self.base_field.get_modulus()

    def get_extension_degree_over_prime_field(self):
        return (
            CubicExtension.EXTENSION_DEGREE_OVER_BASE_FIELD
            * self.base_field.get_extension_degree_over_prime_field()
        )

    def identity(self):
        return CubicExtensionElement(
            self,
            self.base_field.identity(),
            self.base_field.zero(),
            self.base_field.zero(),
        )

    def zero(self):
        return CubicExtensionElement(
            self, self.base_field.zero(), self.base_field.zero(), self.base_field.zero()
        )

    def v(self):
        return CubicExtensionElement(
            self,
            self.base_field.zero(),
            self.base_field.identity(),
            self.base_field.zero(),
        )

    def generate_random_point(self):
        return CubicExtensionElement(
            self,
            self.base_field.generate_random_point(),
            self.base_field.generate_random_point(),
            self.base_field.generate_random_point(),
        )

    def from_list(self, elements: list[int]):
        """Reads a list of integers into an element of CubicExtension."""
        extension_degree_over_prime_field = self.get_extension_degree_over_prime_field()

        assert len(elements) == extension_degree_over_prime_field

        return CubicExtensionElement(
            self,
            self.base_field.from_list(
                elements[: extension_degree_over_prime_field // 3]
            ),
            self.base_field.from_list(
                elements[
                    extension_degree_over_prime_field
                    // 3 : extension_degree_over_prime_field * 2 // 3
                ]
            ),
            self.base_field.from_list(
                elements[extension_degree_over_prime_field * 2 // 3 :]
            ),
        )

    def deserialise(self, serialised: list[bytes]):
        """Reads a sequence of bytes into an element of CubicExtension"""
        length = (
            (self.get_modulus().bit_length() + 8)
            // 8
            * self.get_extension_degree_over_prime_field()
        )
        assert len(serialised) == length

        return CubicExtensionElement(
            self,
            self.base_field.deserialise(serialised[: length // 3]),
            self.base_field.deserialise(serialised[length // 3 : length * 2 // 3]),
            self.base_field.deserialise(serialised[length * 2 // 3 :]),
        )


class CubicExtensionElement:
    def __init__(self, field: CubicExtension, x0, x1, x2):
        self.field = field
        self.x0 = x0
        self.x1 = x1
        self.x2 = x2
        return

    def is_same_field(self, other):
        result = True
        result &= type(other) is type(self)
        return result if not result else result & (self.field == other.field)

    def get_modulus(self):
        """Return the modulus (i.e., characteristic) of the base field."""
        return self.field.get_modulus()

    def copy_with_same_field(self):
        return CubicExtensionElement(
            self.field,
            self.x0.copy_with_same_field(),
            self.x1.copy_with_same_field(),
            self.x2.copy_with_same_field(),
        )

    def __eq__(element, other):
        result = True
        result &= element.is_same_field(other)
        result &= element.x0 == other.x0
        result &= element.x1 == other.x1
        return result & (element.x2 == other.x2)

    def __add__(element, other):
        if type(other) is PrimeFieldElement:
            # Multiplication by element in prime field
            assert element.get_modulus() == other.get_modulus()
            return CubicExtensionElement(
                element.field, element.x0 + other, element.x1, element.x2
            )
        elif element.is_same_field(other):
            # Same field
            return CubicExtensionElement(
                element.field,
                element.x0 + other.x0,
                element.x1 + other.x1,
                element.x2 + other.x2,
            )
        elif (
            gcd(
                element.field.get_extension_degree_over_prime_field(),
                other.field.get_extension_degree_over_prime_field(),
            )
            == 1
        ):
            # Neither field can be an extension of the other
            raise ValueError("Operation not implemented")
        else:
            # One of the two fields is an extension of the other
            a = element.copy_with_same_field()
            b = other.copy_with_same_field()
            # Ensure the first element is the one in the largest extension degree
            if (
                a.field.get_extension_degree_over_prime_field()
                < b.field.get_extension_degree_over_prime_field()
            ):
                a, b = b, a
            if a.field.EXTENSION_DEGREE_OVER_BASE_FIELD == 2:
                if (
                    a.field.get_extension_degree_over_prime_field() // 2
                    < b.field.get_extension_degree_over_prime_field()
                ):
                    # Neither field can be an extension of the other
                    raise ValueError("Operation not implemented")
                target_type = type(a)
                return target_type(a.field, a.x0 + b, a.x1)
            elif a.field.EXTENSION_DEGREE_OVER_BASE_FIELD == 3:
                if (
                    a.field.get_extension_degree_over_prime_field() // 2
                    < b.field.get_extension_degree_over_prime_field()
                ):
                    # Neither field can be an extension of the other
                    raise ValueError("Operation not implemented")
                target_type = type(a)
                return target_type(a.field, a.x0 + b, a.x1, a.x2)
            else:
                raise ValueError("Operation not implemented")

    def __sub__(element, other):
        return element + (-other)

    def __neg__(self):
        return CubicExtensionElement(
            self.field,
            -self.x0,
            -self.x1,
            -self.x2,
        )

    def __mul__(element, other):
        if type(other) is PrimeFieldElement:
            # Multiplication by element in prime field
            assert element.get_modulus() == other.get_modulus()
            return CubicExtensionElement(
                element.field,
                element.x0 * other,
                element.x1 * other,
                element.x2 * other,
            )
        elif element.is_same_field(other):
            # Same field
            return CubicExtensionElement(
                element.field,
                element.x0 * other.x0
                + (element.x1 * other.x2 + element.x2 * other.x1)
                * element.field.non_residue,
                element.x0 * other.x1
                + element.x1 * other.x0
                + element.x2 * other.x2 * element.field.non_residue,
                element.x0 * other.x2 + element.x1 * other.x1 + element.x2 * other.x0,
            )
        elif (
            gcd(
                element.field.get_extension_degree_over_prime_field(),
                other.field.get_extension_degree_over_prime_field(),
            )
            == 1
        ):
            # Neither field can be an extension of the other
            raise ValueError("Operation not implemented")
        else:
            # One of the two fields is an extension of the other
            a = element.copy_with_same_field()
            b = other.copy_with_same_field()
            # Ensure the first element is the one in the largest extension degree
            if (
                a.field.get_extension_degree_over_prime_field()
                < b.field.get_extension_degree_over_prime_field()
            ):
                a, b = b, a
            if a.field.EXTENSION_DEGREE_OVER_BASE_FIELD == 2:
                if (
                    a.field.get_extension_degree_over_prime_field() // 2
                    < b.field.get_extension_degree_over_prime_field()
                ):
                    # Neither field can be an extension of the other
                    raise ValueError("Operation not implemented")
                target_type = type(a)
                return target_type(a.field, a.x0 * b, a.x1 * b)
            elif a.field.EXTENSION_DEGREE_OVER_BASE_FIELD == 3:
                if (
                    a.field.get_extension_degree_over_prime_field() // 3
                    < b.field.get_extension_degree_over_prime_field()
                ):
                    # Neither field can be an extension of the other
                    raise ValueError("Operation not implemented")
                target_type = type(a)
                return target_type(a.field, a.x0 * b, a.x1 * b, a.x2 * b)
            else:
                raise ValueError("Operation not implemented")

    def __repr__(self):
        return f"Fq{self.field.get_extension_degree_over_prime_field()}(\n{self.x0},\n{self.x1},\n{self.x2}\n)"

    def invert(self):
        assert not self.is_zero()

        a = self.x0 * self.x0 - self.x1 * self.x2 * self.field.non_residue
        b = self.x2 * self.x2 * self.field.non_residue - self.x0 * self.x1
        c = self.x1 * self.x1 - self.x0 * self.x2
        d = (
            self.x1 * self.field.non_residue * c
            + self.x0 * a
            + self.x2 * self.field.non_residue * b
        )
        e = d.invert()

        return CubicExtensionElement(self.field, a * e, b * e, c * e)

    def is_zero(self):
        return self.x0.is_zero() and self.x1.is_zero() and self.x2.is_zero()

    def scalar_mul(self, n: int):
        return CubicExtensionElement(
            self.field,
            self.x0.scalar_mul(n),
            self.x1.scalar_mul(n),
            self.x2.scalar_mul(n),
        )

    def power(self, n: int):
        if self.is_zero():
            if n != 0:
                return self.copy_with_same_field()
            else:
                return ValueError("0^0 is not defined")

        if n == 0:
            return self.field.identity()

        val = self.copy_with_same_field()
        result = self.field.identity()

        if n < 0:
            n = -n
            val = val.invert()

        while n > 0:
            if n % 2 == 1:
                result = result * val
            val = val * val
            n = n // 2

        return result

    def frobenius(self, n: int):
        """Frobenius: f -> f^q^n."""
        gamma_x1 = self.field.non_residue.power(
            (
                self.get_modulus()
                ** (n % self.field.get_extension_degree_over_prime_field())
                - 1
            )
            // 3
        )
        gamma_x2 = self.field.non_residue.power(
            2
            * (
                self.get_modulus()
                ** (n % self.field.get_extension_degree_over_prime_field())
                - 1
            )
            // 3
        )

        return CubicExtensionElement(
            self.field,
            self.x0.frobenius(n),
            self.x1.frobenius(n) * gamma_x1,
            self.x2.frobenius(n) * gamma_x2,
        )

    def to_list(self):
        """Convert element to list of its coordinates. The order is: x0, x1, x2."""
        return [*self.x0.to_list(), *self.x1.to_list(), *self.x2.to_list()]

    def to_bytes(self):
        """Serialise the CubicExtensionElement element as its little-endian byte representation. The order is: x0, x1, x2."""
        return self.x0.to_bytes() + self.x1.to_bytes() + self.x2.to_bytes()

    def serialise(self):
        """Serialise the CubicExtensionElement element as the list of bytes of its little-endian representation. The order is: x0, x1, x2."""
        return list(self.to_bytes())
