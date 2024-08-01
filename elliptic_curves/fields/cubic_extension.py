from copy import deepcopy

# The following class is not meant to be used by the user. It should be re-exported using the function below
class CubicExtension:
    """
    Field F, quadratic extension of base field B.
    F = B[u] / (u^2 - non_residue)
    """
    NON_RESIDUE = None
    BASE_FIELD = None
    EXTENSION_DEGREE = None

    def __init__(self, x0, x1, x2):
        self.x0 = x0
        self.x1 = x1
        self.x2 = x2

        return

    def __eq__(x,y):
        return x.__dict__== y.__dict__

    def __add__(x,y):
        assert(type(x) == type(y))
        Field = type(x)

        return Field(x.x0 + y.x0, x.x1 + y.x1, x.x2 + y.x2)
    
    def __sub__(x,y):
        assert(type(x) == type(y))
        Field = type(x)

        return Field(x.x0 - y.x0, x.x1 - y.x1, x.x2 - y.x2)

    def __neg__(self):
        Field = type(self)

        return Field(-self.x0,-self.x1,-self.x2)

    def __mul__(x,y):
        Field = type(x)
        Field_Y = type(y)

        if Field_Y == Field: # Same type
            return Field(
            x.x0 * y.x0 + (x.x1 * y.x2  + x.x2 * y.x1) * Field.NON_RESIDUE,
            x.x0 * y.x1 + x.x1 * y.x0 + x.x2 * y.x2 * Field.NON_RESIDUE, 
            x.x0 * y.x2 + x.x1 * y.x1 + x.x2 * y.x0
            )
        else:
            a = deepcopy(x)
            b = deepcopy(y)
            # Ensure first element is in the largest extension
            if Field.EXTENSION_DEGREE < Field_Y.EXTENSION_DEGREE:
                a, b = b, a
                Field, Field_Y = Field_Y, Field
            # Try the cubic extension case
            try:
                return Field(a.x0 * b, a.x1 * b, a.x2 * b)
            # Otherwise, raise exception
            except:
                raise ValueError('Multiplication not implemented')

    def __repr__(self):
        return f'({self.x0},{self.x1},{self.x2})'

    def invert(self):
        assert(not self.is_zero())
        Field = type(self)

        a = self.x0 * self.x0 - self.x1 * self.x2 * Field.NON_RESIDUE
        b = self.x2 * self.x2 * Field.NON_RESIDUE - self.x0 * self.x1
        c = self.x1 * self.x1 - self.x0 * self.x2
        d = self.x1 * Field.NON_RESIDUE * c + self.x0 * a + self.x2 * Field.NON_RESIDUE * b
        e = d.invert()

        return Field(a * e, b * e, c * e)

    def identity():
        return CubicExtension(CubicExtension.BASE_FIELD.identity(),CubicExtension.BASE_FIELD.zero(),CubicExtension.BASE_FIELD.zero())

    def zero():
        return CubicExtension(CubicExtension.BASE_FIELD.zero(),CubicExtension.BASE_FIELD.zero(),CubicExtension.BASE_FIELD.zero())

    def is_zero(self):
        return self.x0.is_zero() and self.x1.is_zero() and self.x2.is_zero()

    def scalar_mul(self, n:int):
        Field = type(self)
        return Field(self.x0.scalar_mul(n), self.x1.scalar_mul(n), self.x2.scalar_mul(n))

    def v():
        return CubicExtension(CubicExtension.BASE_FIELD.zero(),CubicExtension.BASE_FIELD.identity(),CubicExtension.BASE_FIELD.zero())

    def generate_random_point():
        return CubicExtension(
            CubicExtension.BASE_FIELD.generate_random_point(),
            CubicExtension.BASE_FIELD.generate_random_point(),
            CubicExtension.BASE_FIELD.generate_random_point()
        )

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
    
    def frobenius(self, n:int):
        """
        Frobenius: f -> f^q^n
        """
        Field = type(self)
        gamma_x1 = Field.NON_RESIDUE.power((Field.get_modulus()**(n % Field.EXTENSION_DEGREE)-1)//3)
        gamma_x2 = Field.NON_RESIDUE.power(2*(Field.get_modulus()**(n % Field.EXTENSION_DEGREE)-1)//3)

        return Field(
            self.x0.frobenius(n),
            self.x1.frobenius(n) * gamma_x1,
            self.x2.frobenius(n) * gamma_x2
        )
    
    def get_modulus():
        """
        Get MODULUS (i.e., characteristic) of the field
        """

        return CubicExtension.BASE_FIELD.get_modulus()
    
    def to_list(self):
        """
        Convert element to list of its coordinates. The order is: x0, x1, x2
        """
        x0_list = self.x0.to_list()
        x1_list = self.x1.to_list()
        x2_list = self.x2.to_list()

        out = []
        out.extend(x0_list)
        out.extend(x1_list)
        out.extend(x2_list)

        return out
    
    def from_list(L: list[int]):
        """
        Reads a list into an element of CubicExtension
        """
        assert(len(L) == CubicExtension.EXTENSION_DEGREE)

        tmp = []
        index = 0
        for i in range(3):
            tmp.append(CubicExtension.BASE_FIELD.from_list(L[index:index+CubicExtension.EXTENSION_DEGREE//3]))
            index += CubicExtension.EXTENSION_DEGREE//3
        
        return CubicExtension(tmp[0],tmp[1],tmp[2])
    
    def serialise(self):
        """
        Serialise the CubicExtension element as list of its little-endian byte representation. The order is: x0, x1, x2
        """
        little_endian_x0 = self.x0.serialise()
        little_endian_x1 = self.x1.serialise()
        little_endian_x2 = self.x2.serialise()

        out = []
        out.extend(list(little_endian_x0))
        out.extend(list(little_endian_x1))
        out.extend(list(little_endian_x2))

        return out
    
    def deserialise(L: list[bytes]):
        """
        Reads a sequence of bytes into an element of CubicExtension
        """
        length = (CubicExtension.get_modulus().bit_length() + 8)//8
        assert(len(L) == length * CubicExtension.EXTENSION_DEGREE)


        tmp = []
        index = 0
        for i in range(3):
            tmp.append(CubicExtension.BASE_FIELD.deserialise(L[index:index+length*CubicExtension.EXTENSION_DEGREE//3]))
            index += length * CubicExtension.EXTENSION_DEGREE//3
        
        return CubicExtension(tmp[0],tmp[1],tmp[2])
    
def cubic_extension_from_base_field_and_non_residue(base_field, non_residue):
    """
    Function to export class CubicExtension with BASE_FIELD = base_field and NON_RESIDUE = non_residue
    """

    class CubicExtensionField(CubicExtension):
        NON_RESIDUE = non_residue
        BASE_FIELD = base_field
        EXTENSION_DEGREE = 3 * BASE_FIELD.EXTENSION_DEGREE

        def identity():
            return CubicExtensionField(CubicExtensionField.BASE_FIELD.identity(),CubicExtensionField.BASE_FIELD.zero(),CubicExtensionField.BASE_FIELD.zero())

        def zero():
            return CubicExtensionField(CubicExtensionField.BASE_FIELD.zero(),CubicExtensionField.BASE_FIELD.zero(),CubicExtensionField.BASE_FIELD.zero())
        
        def v():
            return CubicExtensionField(
                CubicExtensionField.BASE_FIELD.zero(),
                CubicExtensionField.BASE_FIELD.identity(),
                CubicExtensionField.BASE_FIELD.zero()
                )
        
        def generate_random_point():
            x0 = CubicExtensionField.BASE_FIELD.generate_random_point()
            x1 = CubicExtensionField.BASE_FIELD.generate_random_point()
            x2 = CubicExtensionField.BASE_FIELD.generate_random_point()

            return CubicExtensionField(x0,x1,x2)
        
        def get_modulus():
            """
            Return the MODULUS (i.e., characteristic) of the field
            """
            return CubicExtensionField.BASE_FIELD.get_modulus()
        
        def from_list(L: list[int]):
            """
            Reads a list into an element of CubicExtensionField
            """
            assert(len(L) == CubicExtensionField.EXTENSION_DEGREE)

            tmp = []
            index = 0
            for i in range(3):
                tmp.append(CubicExtensionField.BASE_FIELD.from_list(L[index:index+CubicExtensionField.EXTENSION_DEGREE//3]))
                index += CubicExtensionField.EXTENSION_DEGREE//3
            
            return CubicExtensionField(tmp[0],tmp[1],tmp[2])
    
        def serialise(self):
            """
            Serialise the CubicExtensionField element as list of its little-endian byte representation. The order is: x0, x1, x2
            """
            little_endian_x0 = self.x0.serialise()
            little_endian_x1 = self.x1.serialise()
            little_endian_x2 = self.x2.serialise()

            out = []
            out.extend(list(little_endian_x0))
            out.extend(list(little_endian_x1))
            out.extend(list(little_endian_x2))

            return out
        
        def deserialise(L: list[bytes]):
            """
            Reads a sequence of bytes into an element of CubicExtensionField
            """
            length = (CubicExtensionField.get_modulus().bit_length() + 8)//8
            assert(len(L) == length * CubicExtensionField.EXTENSION_DEGREE)


            tmp = []
            index = 0
            for i in range(3):
                tmp.append(CubicExtensionField.BASE_FIELD.deserialise(L[index:index+length*CubicExtensionField.EXTENSION_DEGREE//3]))
                index += length * CubicExtensionField.EXTENSION_DEGREE//3
            
            return CubicExtensionField(tmp[0],tmp[1],tmp[2])
        
    return CubicExtensionField

        

    
