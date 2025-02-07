from secrets import randbelow


class PrimeField:
    EXTENSION_DEGREE_OVER_BASE_FIELD = 1

    def __init__(self, modulus: int):
        self.modulus = modulus
        return

    def __call__(self, x: int):
        return PrimeFieldElement(self, x)

    def __eq__(field, other):
        result = True
        result &= type(other) is type(field)
        return result if not result else result & (field.modulus == other.modulus)

    def get_modulus(self):
        return self.modulus

    def get_extension_degree_over_prime_field(self):
        return 1

    def identity(self):
        return PrimeFieldElement(self, 1)

    def zero(self):
        return PrimeFieldElement(self, 0)

    def generate_random_point(self):
        return PrimeFieldElement(self, randbelow(self.modulus))

    def from_list(self, elements: list[int]):
        """Reads a list into an element of Fq."""
        assert len(elements) == self.EXTENSION_DEGREE_OVER_BASE_FIELD

        return PrimeFieldElement(self, elements[0])

    def deserialise(self, serialised: list[bytes]):
        """Reads a sequence of bytes into an element of Fq."""
        length = (self.modulus.bit_length() + 8) // 8
        assert len(serialised) == length * self.EXTENSION_DEGREE_OVER_BASE_FIELD

        return PrimeFieldElement(
            self, int.from_bytes(bytes=serialised, byteorder="little")
        )


class PrimeFieldElement:
    def __init__(self, field: PrimeField, x: int):
        self.x = x % field.modulus
        self.field = field
        return

    def is_same_field(self, other):
        result = True
        result &= type(other) is type(self)
        return result if not result else result & (self.field == other.field)

    def get_modulus(self):
        """Return the modulus (i.e., characteristic) of the base field."""
        return self.field.get_modulus()

    def copy_with_same_field(self):
        return PrimeFieldElement(self.field, self.x)

    def __eq__(element, other):
        result = True
        result &= element.is_same_field(other)
        return result & (element.x == other.x)

    def __add__(element, other):
        try:
            if other.field.get_extension_degree_over_prime_field() == 1:
                return PrimeFieldElement(element.field, element.x + other.x)
            else:
                return other + element
        except Exception as _:
            raise ValueError("Multiplication not implemented")

    def __sub__(element, other):
        return element + (-other)

    def __neg__(self):
        return PrimeFieldElement(self.field, -self.x)

    def __mul__(element, other):
        try:
            if other.field.get_extension_degree_over_prime_field() == 1:
                return PrimeFieldElement(element.field, element.x * other.x)
            else:
                return other * element
        except Exception as _:
            raise ValueError("Multiplication not implemented")

    def __repr__(self):
        return f"Fq({self.x})"

    def invert(self):
        return PrimeFieldElement(self.field, pow(self.x, -1, self.get_modulus()))

    def power(self, n: int):
        return PrimeFieldElement(self.field, pow(self.x, n, self.get_modulus()))

    def is_zero(self):
        return self.x == 0

    def scalar_mul(self, n: int):
        return PrimeFieldElement(self.field, self.x * n)

    def frobenius(self, _: int):
        """Frobenius morphism: f --> f^q^n."""
        return self.copy_with_same_field()

    def to_int(self) -> int:
        return self.x

    def to_list(self) -> list[int]:
        """Convert element to list of its coordinates."""
        return [self.to_int()]

    def to_bytes(self) -> bytes:
        """Serialise the Fq element as its little-endian byte representation"""
        length = (self.get_modulus().bit_length() + 8) // 8
        return bytearray(self.x.to_bytes(length=length, byteorder="little"))

    def serialise(self) -> list[bytes]:
        """Serialise the Fq element as the list of bytes of its little-endian representation"""
        return list(self.to_bytes())
