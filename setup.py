from setuptools import setup

setup(
    name="elliptic_curves",
    version="0.1.0",
    author="Federico Barbacovi",
    packages=[
        "elliptic_curves",
        "elliptic_curves.instantiations",
        "elliptic_curves.models",
        "elliptic_curves.fields",
        "elliptic_curves.data_structures",
    ],
)
