# Elliptic_curves

Python library with implementations of:
- [Finite fields](./docs/fields.md) (prime fields: 'Fq', quadratic extensions: 'QuadraticExtension', and cubic extensions: 'CubicExtension')
- [Elliptic curves](./docs/elliptic_curves.md) in Short-Weierstrass form: 'EllipticCurve' and 'ElliptiCurveProjective'
- [Bilinear pairings](./docs/bilinear_pairings.md): 'BilinearPairingCurve'

## Instantiations currenly implemented:

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