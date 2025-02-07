from math import gcd

from elliptic_curves.fields.prime_field import PrimeFieldElement


# The following class is not meant to be used by the user. It should be re-exported using the function below
class QuadraticExtension:
    """
    Field F, quadratic extension of base field B.
    F = B[u] / (u^2 - non_residue)
    """

    EXTENSION_DEGREE_OVER_BASE_FIELD = 2

    def __init__(self, base_field, non_residue: PrimeFieldElement):
        self.base_field = base_field
        self.non_residue = non_residue
        return

    def __call__(self, x0: PrimeFieldElement, x1: PrimeFieldElement):
        return QuadraticExtensionElement(self, x0, x1)

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
            QuadraticExtension.EXTENSION_DEGREE_OVER_BASE_FIELD
            * self.base_field.get_extension_degree_over_prime_field()
        )

    def identity(self):
        return QuadraticExtensionElement(
            self, self.base_field.identity(), self.base_field.zero()
        )

    def zero(self):
        return QuadraticExtensionElement(
            self, self.base_field.zero(), self.base_field.zero()
        )

    def u(self):
        return QuadraticExtensionElement(
            self, self.base_field.zero(), self.base_field.identity()
        )

    def generate_random_point(self):
        return QuadraticExtensionElement(
            self,
            self.base_field.generate_random_point(),
            self.base_field.generate_random_point(),
        )

    def from_list(self, elements: list[int]):
        """Reads a list of integers into an element of QuadraticExtension."""
        extension_degree_over_prime_field = self.get_extension_degree_over_prime_field()

        assert len(elements) == extension_degree_over_prime_field

        return QuadraticExtensionElement(
            self,
            self.base_field.from_list(
                elements[: extension_degree_over_prime_field // 2]
            ),
            self.base_field.from_list(
                elements[extension_degree_over_prime_field // 2 :]
            ),
        )

    def deserialise(self, serialised: list[bytes]):
        """Reads a sequence of bytes into an element of QuadraticExtension"""
        length = (
            (self.get_modulus().bit_length() + 8)
            // 8
            * self.get_extension_degree_over_prime_field()
        )
        assert len(serialised) == length

        return QuadraticExtensionElement(
            self,
            self.base_field.deserialise(serialised[: length // 2]),
            self.base_field.deserialise(serialised[length // 2 :]),
        )


class QuadraticExtensionElement:
    def __init__(self, field: QuadraticExtension, x0, x1):
        self.field = field
        self.x0 = x0
        self.x1 = x1
        return

    def is_same_field(self, other):
        result = True
        result &= type(other) is type(self)
        return result if not result else result & (self.field == other.field)

    def get_modulus(self):
        """Return the modulus (i.e., characteristic) of the base field."""
        return self.field.get_modulus()

    def copy_with_same_field(self):
        return QuadraticExtensionElement(
            self.field,
            self.x0.copy_with_same_field(),
            self.x1.copy_with_same_field(),
        )

    def __eq__(element, other):
        result = True
        result &= element.is_same_field(other)
        result &= element.x0 == other.x0
        return result & (element.x1 == other.x1)

    def __add__(element, other):
        if type(other) is PrimeFieldElement:
            # Multiplication by element in prime field
            assert element.get_modulus() == other.get_modulus()
            return QuadraticExtensionElement(
                element.field, element.x0 + other, element.x1
            )
        elif element.is_same_field(other):
            # Same field
            return QuadraticExtensionElement(
                element.field,
                element.x0 + other.x0,
                element.x1 + other.x1,
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
                    a.field.get_extension_degree_over_prime_field() // 3
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
        return QuadraticExtensionElement(self.field, -self.x0, -self.x1)

    def __mul__(element, other):
        if type(other) is PrimeFieldElement:
            # Multiplication by element in prime field
            assert element.get_modulus() == other.get_modulus()
            return QuadraticExtensionElement(
                element.field, element.x0 * other, element.x1 * other
            )
        elif element.is_same_field(other):
            # Same field
            return QuadraticExtensionElement(
                element.field,
                element.x0 * other.x0
                + element.x1 * other.x1 * element.field.non_residue,
                element.x0 * other.x1 + element.x1 * other.x0,
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
        return f"Fq{self.field.get_extension_degree_over_prime_field()}(\n{self.x0},\n{self.x1}\n)"

    def conjugate(self):
        return QuadraticExtensionElement(self.field, self.x0, -self.x1)

    def invert(self):
        assert not self.is_zero()

        z = self.x0.power(2) - self.x1.power(2) * self.field.non_residue
        z = z.invert()

        conjugate = self.conjugate()

        return QuadraticExtensionElement(self.field, conjugate.x0 * z, conjugate.x1 * z)

    def is_zero(self):
        return self.x0.is_zero() and self.x1.is_zero()

    def scalar_mul(self, n: int):
        return QuadraticExtensionElement(
            self.field, self.x0.scalar_mul(n), self.x1.scalar_mul(n)
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
        gamma = self.field.non_residue.power(
            (
                self.get_modulus()
                ** (n % self.field.get_extension_degree_over_prime_field())
                - 1
            )
            // 2
        )

        return QuadraticExtensionElement(
            self.field,
            self.x0.frobenius(n),
            self.x1.frobenius(n) * gamma,
        )

    def to_list(self):
        """Convert element to list of its coordinates. The order is: x0, x1."""

        return [*self.x0.to_list(), *self.x1.to_list()]

    def to_bytes(self):
        """Serialise the QuadraticExtensionElement element as its little-endian byte representation. The order is: x0, x1."""
        return self.x0.to_bytes() + self.x1.to_bytes()

    def serialise(self):
        """Serialise the QuadraticExtensionElement element as the list of bytes of its little-endian representation. The order is: x0, x1."""
        return list(self.to_bytes())
