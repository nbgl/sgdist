from setuptools import setup
from setuptools.extension import Extension

from Cython.Build import cythonize


setup(
    name="sgdist",
    version="0.0",
    packages=['sgdist'],
    ext_modules=cythonize([Extension(
        name="sgdist.sgdistglue",
        sources=["sgdist/sgdistglue.pyx"],
        extra_compile_args=["-O3"],
    )]),
    zip_safe=False,
)
