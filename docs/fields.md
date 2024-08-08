# Finite fields

The library implements three basic classes to work with finite fields:
- `Fq`, which is the class implementing finite fields of prime cardinality
- `QuadraticExtension`, which is the class implementing quadratic extension fields
- `CubicExtension`, which is the class implementing cubic extensions

These classes are not meant to be directly used by the user, rather they are meant to be exported with the functions provided. Below are some examples.

## Prime fields

```python
from elliptic_curves.fields.fq import base_field_from_modulus

# By using the above function, we can instantiate a prime field
Fq = base_field_from_modulus(q=17)

# Then, we can define elements of the field
a = Fq(5)
b = Fq(7)

# Sample them at random
c = Fq.generate_random_point()

# Perfom operations on them
assert(a + b == Fq(12))
assert(c.power(2) == c * c)
assert(c * c.invert() == c)
```

## Quadratic extensions

```python
from elliptic_curves.fields.fq import base_field_from_modulus
from elliptic_curves.fields.quadratic_extension import quadratic_extension_from_base_field_and_non_residue

Fq = base_field_from_modulus(q=17)
# Let's build a quadratic extension of Fq; (-3)^(8) = -1 mod 17 --> -3 is a non quadratic residue
Fq2 = quadratic_extension_from_base_field_and_non_residue(base_field=Fq,non_residue=Fq(-3))

# Now we can work with Fq2 as we did with Fq
a = Fq2.generate_random_point()
b = Fq2(Fq(1),Fq(1))
u = Fq2.u() # The element (0,1)

assert(a.invert() * a == a)
assert(a.frobenius(2) == a.power(Fq2.get_modulus()**2))
assert(b + u == Fq2(Fq(1),Fq(2)))

# As Fq2 is an extension of Fq, we can multiply elements of Fq with elements of Fq2
assert(Fq(2) * b == Fq2(Fq(2),Fq(2)))
```

## Cubic extensions

```python
from elliptic_curves.fields.fq import base_field_from_modulus
from elliptic_curves.fields.cubic_extension import cubic_extension_from_base_field_and_non_residue
from elliptic_curves.fields.quadratic_extension import quadratic_extension_from_base_field_and_non_residue

Fq = base_field_from_modulus(q=19)
# Let's build a quadratic extension of Fq; (-3)^(6) = 7 mod 19 --> -3 is a non cubic residue
Fq3 = cubic_extension_from_base_field_and_non_residue(base_field=Fq,non_residue=Fq(-3))

# Now we can work with Fq3 as we did with Fq and Fq2
a = Fq3.generate_random_point()
b = Fq3(Fq(1),Fq(1),Fq(1))
v = Fq3.v() # The element (0,1,0)

assert(a.invert() * a == a)
assert(a.frobenius(3) == a.power(Fq2.get_modulus()**3))
assert(v.power(3) == Fq3(Fq(-3),Fq.zero(),Fq.zero()))

# As Fq3 is an extension of Fq, we can multiply elements of Fq with elements of Fq3
assert(Fq(2) * b == Fq3(Fq(2),Fq(2),Fq(2)))

# However, if we tried to multiply an element of Fq3 with an element of Fq2 we would get an error
Fq2 = quadratic_extension_from_base_field_and_non_residue(base_field=Fq,non_residue=Fq(-1))
try:
    Fq2.identity() * Fq3.identity()
except:
    print('The multiplication returned an error because it cannot be carried out')
```