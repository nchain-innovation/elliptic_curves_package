from copy import deepcopy
from secrets import randbelow

# The following class is not meant to be used by the user. It should be re-exported using the function below
class Fq:
    EXTENSION_DEGREE = 1
    MODULUS = None

    def __init__(self, x: int):
        Field = type(self)
        self.x = x % Field.MODULUS

        return

    def __eq__(x,y):
        return x.__dict__== y.__dict__

    def __add__(x,y):
        assert(type(x) == type(y))
        Field = type(x)

        return Field(x.x + y.x)

    def __sub__(x,y):
        assert(type(x) == type(y))
        Field = type(x)

        return Field(x.x - y.x)

    def __neg__(self):
        Field = type(self)

        return Field(-self.x)

    def __mul__(x,y):
        Field = type(x)

        if type(y) == Field:
            return Field(x.x * y.x)
        else:
            return y.scalar_mul(x.x)

    def __repr__(self):
        return f'{self.x}'

    def invert(self):
        Field = type(self)

        return Field(pow(self.x,-1,Field.MODULUS))

    def power(self, n:int):
        Field = type(self)

        return Field(pow(self.x,n,Field.MODULUS))

    def identity():
        return Fq(1)

    def zero():
        return Fq(0)

    def is_zero(self):
        return self.x == 0

    def scalar_mul(self, n:int):
        Field = type(self)

        return Field(self.x * n)

    def generate_random_point():
        return Fq(randbelow(Fq.MODULUS))
    
    def frobenius(self, n: int):
        """
        Frobenius morphism: f --> f^q^n
        """

        return deepcopy(self)
    
    def get_modulus():
        """
        Return the MODULUS (i.e., characteristic) of the field
        """
        return Fq.MODULUS
    
    def to_list(self):
        """
        Convert element to list of its coordinates
        """
        return [self.x]
    
    def from_list(L: list[int]):
        """
        Reads a list into an element of Fq
        """

        assert(len(L) == Fq.EXTENSION_DEGREE)

        return Fq(L[0])
    
    def serialise(self):
        """
        Serialise the Fq element as its little-endian byte representation
        """
        Field = type(self)
        length = (Field.MODULUS.bit_length()+8)//8
        little_endian = bytearray(self.x.to_bytes(length=length,byteorder='little'))
        
        return list(little_endian)

    def deserialise(L: list[bytes]):
        """
        Reads a sequence of bytes into an element of Fq2
        """
        length = (Fq.MODULUS.bit_length() + 8)//8
        assert(len(L) == length * Fq.EXTENSION_DEGREE)

        x0 = int.from_bytes(bytes=L,byteorder='little')

        return Fq(x0)
    
def base_field_from_modulus(q: int):
    """
    Function to export class Fq with MODULUS set to q
    """

    class Field(Fq):
        MODULUS = q
        EXTENSION_DEGREE = 1

        def identity():
            return Field(1)

        def zero():
            return Field(0)
        
        def generate_random_point():
            return Field(randbelow(Field.MODULUS))
        
        def get_modulus():
            """
            Return the MODULUS (i.e., characteristic) of the field
            """
            return Field.MODULUS
        
        def from_list(L: list[int]):
            """
            Reads a list into an element of Fq
            """

            assert(len(L) == Field.EXTENSION_DEGREE)

            return Field(L[0])
        
        def deserialise(L: list[bytes]):
            """
            Reads a sequence of bytes into an element of Fq2
            """
            length = (Field.MODULUS.bit_length() + 8)//8
            assert(len(L) == length * Field.EXTENSION_DEGREE)

            x0 = int.from_bytes(bytes=L,byteorder='little')

            return Field(x0)
        
    return Field