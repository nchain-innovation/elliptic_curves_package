# Finite fields

The library implements two basic classes to work with elliptic curves in Short-Weierstrass form:
- `EllipticCurve`, which is the class implementing affine elliptic curves in SW form
- `EllipticCurveProjective`, which is the class implementing projective elliptic curves in SW form

These classes are not meant to be directly used by the user, rather they are meant to be exported with the functions provided. Below is an example:

```python
from elliptic_curves.fields.fq import base_field_from_modulus

from elliptic_curves.models.ec import elliptic_curve_from_curve
from elliptic_curves.models.curve import Curve

# Prime field
Fq = base_field_from_modulus(q=19)

# We define a Curve object to keep track of the curve parameters
curve = Curve(a = Fq(0), b = Fq(8))

# We now export affine and projective version of curve
affine_curve, projective_curve = elliptic_curve_from_curve(curve=curve)

# We can now create points
P = affine_curve(x = Fq(1), y = Fq(3))
Q = affine_curve(x = Fq(1), y = Fq(-3))

# We can sum them
R = P + Q
S = P + P
assert(S == P.multiply(2))

# We can get the lambda of the line through two points
lambda_P_S = P.get_lamdba(S)
assert(lambda_P_S * (S.x - P.x)  == (S.y - P.y))

# We can convert a point to its list of coordinates
list_coordinates_P = P.to_list()
assert(list_coordinates_P == [1,3])

list_coordinates_Q = Q.to_list()
assert(list_coordinates_P == [1,16])
```