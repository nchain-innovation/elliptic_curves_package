# Python implementation of the following models:
- Finite fields (Fq, Quadratic extensions and Cubic extensions)
- Elliptic curves in Short-Weierstrass form
- Bilinear pairings

## Instantiations of:
- BLS12_381
- MNT4_753

The user wishing to implement new bilinear pairing needs to fetch the relevant defining parameters (modulus, field extensions, curve parameters) and then implement the final exponentiation of the bilinear pairing. See BLS12_381 or MNT4_753 for an example.

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

# Exectue the tests in the python virtual environment

```bash
cd tests/instantiations/bls13_381
python3 bls12_381_test.py

cd tests/instantiations/mnt4_753
python3 mnt4_753_test.py