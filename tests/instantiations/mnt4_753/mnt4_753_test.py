import os, sys

from elliptic_curves.fields.fq import base_field_from_modulus
from elliptic_curves.instantiations.mnt4_753.mnt4_753 import Fq4, mnt4_753, mnt4_753_bilinear_pairing

g1 = mnt4_753.g1
g2 = mnt4_753.g2
Fr = base_field_from_modulus(q=mnt4_753.r)

pairing_g1_g2_serialised = [180, 141, 84, 186, 49, 46, 169, 53, 16, 153, 219, 27, 34, 98, 128, 24, 37, 74, 39, 167, 161, 218, 126, 199, 67, 184, 51, 137, 211, 169, 15, 253, 216, 67, 175, 233, 248, 94, 247, 223, 74, 127, 23, 207, 65, 99, 69, 70, 126, 118, 31, 139, 78, 69, 18, 28, 45, 52, 16, 217, 210, 68, 9, 33, 2, 56, 23, 107, 121, 78, 42, 56, 216, 157, 89, 151, 99, 90, 207, 213, 193, 103, 239, 241, 227, 232, 146, 107, 111, 173, 164, 48, 234, 253, 0, 77, 9, 236, 2, 6, 48, 40, 104, 172, 16, 112, 229, 104, 203, 208, 19, 77, 54, 115, 80, 46, 27, 199, 95, 246, 174, 128, 129, 95, 13, 15, 30, 215, 27, 169, 43, 64, 152, 84, 11, 154, 64, 183, 35, 169, 78, 222, 117, 118, 218, 110, 12, 95, 106, 72, 31, 51, 207, 72, 87, 205, 252, 28, 14, 220, 34, 98, 229, 39, 172, 192, 27, 154, 243, 147, 20, 66, 128, 47, 37, 130, 103, 158, 93, 205, 211, 209, 61, 137, 13, 170, 96, 49, 61, 0, 78, 120, 203, 70, 50, 122, 88, 84, 124, 83, 102, 224, 159, 111, 216, 209, 230, 225, 140, 75, 167, 224, 64, 106, 39, 233, 221, 19, 63, 78, 8, 133, 142, 122, 32, 203, 28, 215, 23, 130, 189, 186, 22, 226, 87, 122, 202, 32, 205, 66, 77, 36, 163, 203, 75, 104, 120, 82, 21, 150, 109, 20, 207, 215, 34, 24, 76, 145, 96, 198, 69, 154, 48, 100, 61, 67, 247, 134, 206, 10, 235, 204, 164, 22, 133, 226, 10, 24, 190, 83, 90, 47, 75, 98, 1, 60, 201, 126, 61, 9, 85, 24, 152, 99, 215, 188, 54, 66, 12, 42, 205, 142, 112, 119, 129, 139, 202, 85, 128, 53, 32, 136, 123, 191, 61, 153, 158, 34, 133, 221, 142, 85, 249, 16, 191, 15, 168, 84, 242, 177, 87, 213, 29, 189, 227, 232, 84, 136, 161, 132, 36, 202, 173, 247, 138, 37, 186, 19, 191, 188, 145, 253, 61, 49, 153, 96, 33, 199, 36, 244, 219, 252, 130, 228, 206, 26, 10, 233, 102, 22, 117, 226, 127, 32, 120, 31, 138, 237, 169, 0]
pairing_g1_g2 = Fq4.deserialise(pairing_g1_g2_serialised)

def test_pairing() -> bool:
    miller_output_twisted_curve = mnt4_753_bilinear_pairing.miller_loop_on_twisted_curve(g1,g2,'quadratic')

    pairing_base_curve = mnt4_753_bilinear_pairing.pairing(g1,g2)
    pairing_twisted_curve = miller_output_twisted_curve.power((mnt4_753.q**4-1)//mnt4_753.r)

    assert(pairing_base_curve == pairing_g1_g2)
    assert(pairing_base_curve == pairing_twisted_curve)

    for i in range(5):
        l = Fr.generate_random_point().x
        assert(mnt4_753_bilinear_pairing.pairing(g1.multiply(l),g2) == pairing_base_curve.power(l))
        assert(mnt4_753_bilinear_pairing.pairing(g1,g2.multiply(l)) == pairing_base_curve.power(l))
    
    return True

assert(test_pairing())

print("MNT4_753: all tests successful")

