# Elliptic_curves

Python library with implementations of:
- [Finite fields](./docs/fields.md) (prime fields: `PrimeField`, quadratic extensions: `QuadraticExtension`, and cubic extensions: `CubicExtension`)
- [Elliptic curves](./docs/elliptic_curves.md) in Short-Weierstrass form: `ShortWeierstrassEllipticCurve`
- [Bilinear pairings](./docs/bilinear_pairings.md): `BilinearPairingCurve`

The structure of the library follows in part that of the [Arkworks](https://github.com/arkworks-rs) library.

## Instantiations currently implemented:

The library currently contains instantiations of the following curves:
- BLS12_381
- MNT4_753

# Installation.

To install the elliptic_curves package in the system python directory, please run the command below.

```bash
pip install -e .
```

To install in the python virtual environment

```bash
export VIRTUAL_ENV=<SOME PATH>
python -m venv $VIRTUAL_ENV
source $VIRTUAL_ENV/bin/activate

pip3 install -e . 
```

# Execute the tests in the system environment
```bash
cd tests/instantiations/bls13_381
python bls12_381_test.py

cd tests/instantiations/mnt4_753
python mnt4_753_test.py
```

# Execute the tests in the python virtual environment

```bash
cd tests/instantiations/bls13_381
python3 bls12_381_test.py

cd tests/instantiations/mnt4_753
python3 mnt4_753_test.py
```

# Disclaimer

The code and resources within this repository are intended for research and educational purposes only.
Please note:
- No guarantees are provided regarding the security or the performance of the code.
- Users are responsible for validating the code and understanding its implications before using it in any capacity.
- There may be edge cases causing bugs or unexpected behaviours. Please contact us if you find any bug.

# LICENSE

The code in the current repository is licensed under the MIT license.