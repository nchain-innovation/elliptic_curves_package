from copy import deepcopy

# The following class is not meant to be used by the user. It should be re-exported using the function below
class QuadraticExtension:
    """
    Field F, quadratic extension of base field B.
    F = B[u] / (u^2 - non_residue)
    """
    NON_RESIDUE = None
    BASE_FIELD = None
    EXTENSION_DEGREE = None
    EXTENSION_DEGREE_OVER_BASE_FIELD = None

    def __init__(self, x0, x1):
        self.x0 = x0
        self.x1 = x1

        return

    def __eq__(x,y):
        return x.__dict__== y.__dict__

    def __add__(x,y):
        assert(type(x) == type(y))
        Field = type(x)

        return Field(x.x0 + y.x0, x.x1 + y.x1)

    def __sub__(x,y):
        assert(type(x) == type(y))
        Field = type(x)

        return Field(x.x0 - y.x0, x.x1 - y.x1)

    def __neg__(self):
        Field = type(self)

        return Field(-self.x0,-self.x1)

    def __mul__(x,y):
        Field = type(x)
        Field_Y = type(y)

        if Field_Y == Field: # Same type
            return Field(x.x0 * y.x0 + x.x1 * y.x1 * Field.NON_RESIDUE, x.x0 * y.x1 + x.x1 * y.x0)
        else:
            a = deepcopy(x)
            b = deepcopy(y)
            # Ensure first element is in the largest extension
            if Field.EXTENSION_DEGREE < Field_Y.EXTENSION_DEGREE:
                a, b = b, a
                Field, Field_Y = Field_Y, Field
            # Try the quadratic extension case
            try:
                return Field(a.x0 * b, a.x1 * b)
            # Otherwise, raise exception
            except:
                raise ValueError('Multiplication not implemented')
            
    def __repr__(self):
        return f'({self.x0},{self.x1})'

    def conjugate(self):
        Field = type(self)

        return Field(self.x0,-self.x1)

    def invert(self):
        assert(not self.is_zero())
        Field = type(self)

        z = self.x0.power(2) - self.x1.power(2) * Field.NON_RESIDUE
        z = z.invert()

        conjugate = self.conjugate()

        return Field(conjugate.x0 * z, conjugate.x1 * z)

    def identity():
        return QuadraticExtension(QuadraticExtension.BASE_FIELD.identity(),QuadraticExtension.BASE_FIELD.zero())

    def zero():
        return QuadraticExtension(QuadraticExtension.BASE_FIELD.zero(),QuadraticExtension.BASE_FIELD.zero())

    def is_zero(self):
        return self.x0.is_zero() and self.x1.is_zero()

    def scalar_mul(self, n:int):
        Field = type(self)

        return Field(self.x0.scalar_mul(n), self.x1.scalar_mul(n))

    def u():
        return QuadraticExtension(QuadraticExtension.BASE_FIELD.zero(),QuadraticExtension.BASE_FIELD.identity())

    def power(self,n: int):
        if self.is_zero():
            if n != 0:
                return deepcopy(self)
            else:
                return ValueError('0^0 is not defined')
        
        Field = type(self)
        if n == 0:
            return Field.identity()
            
        val = deepcopy(self)
        result = Field.identity()

        if n < 0:
            n = -n
            val = val.invert()

        while n > 0:
            if n % 2 == 1:
                result = result * val
            val = val * val
            n = n // 2

        return result

    def generate_random_point():
        x0 = QuadraticExtension.BASE_FIELD.generate_random_point()
        x1 = QuadraticExtension.BASE_FIELD.generate_random_point()

        return QuadraticExtension(x0,x1)
    
    def frobenius(self, n:int):
        """
        Frobenius: f -> f^q^n
        """
        Field = type(self)
        gamma = Field.NON_RESIDUE.power((Field.get_modulus()**(n % Field.EXTENSION_DEGREE)-1)//2)

        return Field(self.x0.frobenius(n), self.x1.frobenius(n) * gamma)
    
    def get_modulus():
        """
        Get MODULUS (i.e., characteristic) of the field
        """

        return QuadraticExtension.BASE_FIELD.get_modulus()
    
    def to_list(self):
        """
        Convert element to list of its coordinates. The order is: x0, x1
        """
        
        x0_list = self.x0.to_list()
        x1_list = self.x1.to_list()

        out = []
        out.extend(x0_list)
        out.extend(x1_list)

        return out
    
    def from_list(L: list[int], q: int):
        """
        Reads a list into an element of QuadraticExtension
        """
        assert(len(L) == QuadraticExtension.EXTENSION_DEGREE)

        tmp = []
        index = 0
        for i in range(2):
            tmp.append(QuadraticExtension.BASE_FIELD.from_list(L[index:index+QuadraticExtension.EXTENSION_DEGREE//2]))
            index += QuadraticExtension.EXTENSION_DEGREE//2
        
        return QuadraticExtension(tmp[0],tmp[1])

    def serialise(self):
        """
        Serialise the QuadraticExtension element as list of its little-endian byte representation. The order is: x0, x1
        """
        little_endian_x0 = self.x0.serialise()
        little_endian_x1 = self.x1.serialise()

        out = []
        out.extend(list(little_endian_x0))
        out.extend(list(little_endian_x1))

        return out
    
    def deserialise(L: list[bytes]):
        """
        Reads a sequence of bytes into an element of QuadraticExtension
        """

        length = (QuadraticExtension.get_modulus().bit_length() + 8)//8
        assert(len(L) == length * QuadraticExtension.EXTENSION_DEGREE)

        tmp = []
        index = 0
        for i in range(2):
            tmp.append(QuadraticExtension.BASE_FIELD.deserialise(L[index:index+(length*QuadraticExtension.EXTENSION_DEGREE//2)]))
            index += length * QuadraticExtension.EXTENSION_DEGREE//2

        return QuadraticExtension(tmp[0],tmp[1])
    
def quadratic_extension_from_base_field_and_non_residue(base_field, non_residue):
    """
    Function to export class QuadraticExtension with BASE_FIELD = base_field and NON_RESIDUE = non_residue
    """

    class QuadraticExtensionField(QuadraticExtension):
        NON_RESIDUE = non_residue
        BASE_FIELD = base_field
        EXTENSION_DEGREE = 2 * BASE_FIELD.EXTENSION_DEGREE
        EXTENSION_DEGREE_OVER_BASE_FIELD = 2

        def identity():
            return QuadraticExtensionField(QuadraticExtensionField.BASE_FIELD.identity(),QuadraticExtensionField.BASE_FIELD.zero())

        def zero():
            return QuadraticExtensionField(QuadraticExtensionField.BASE_FIELD.zero(),QuadraticExtensionField.BASE_FIELD.zero())
        
        def u():
            return QuadraticExtensionField(QuadraticExtensionField.BASE_FIELD.zero(),QuadraticExtensionField.BASE_FIELD.identity())
        
        def generate_random_point():
            x0 = QuadraticExtensionField.BASE_FIELD.generate_random_point()
            x1 = QuadraticExtensionField.BASE_FIELD.generate_random_point()

            return QuadraticExtensionField(x0,x1)
        
        def get_modulus():
            """
            Return the MODULUS (i.e., characteristic) of the field
            """
            return QuadraticExtensionField.BASE_FIELD.get_modulus()
        
        def from_list(L: list[int]):
            """
            Reads a list into an element of QuadraticExtensionField
            """
            assert(len(L) == QuadraticExtensionField.EXTENSION_DEGREE)

            tmp = []
            index = 0
            for i in range(2):
                tmp.append(QuadraticExtensionField.BASE_FIELD.from_list(L[index:index+QuadraticExtensionField.EXTENSION_DEGREE//2]))
                index += QuadraticExtensionField.EXTENSION_DEGREE//2
            
            return QuadraticExtensionField(tmp[0],tmp[1])
        
        def deserialise(L: list[bytes]):
            """
            Reads a sequence of bytes into an element of QuadraticExtensionField
            """

            length = (QuadraticExtensionField.get_modulus().bit_length() + 8)//8
            assert(len(L) == length * QuadraticExtensionField.EXTENSION_DEGREE)

            tmp = []
            index = 0
            for i in range(2):
                tmp.append(QuadraticExtensionField.BASE_FIELD.deserialise(L[index:index+(length*QuadraticExtensionField.EXTENSION_DEGREE//2)]))
                index += length * QuadraticExtensionField.EXTENSION_DEGREE//2

            return QuadraticExtensionField(tmp[0],tmp[1])
        
    return QuadraticExtensionField

        

    
