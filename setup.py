from setuptools import setup, find_packages, Extension
from Cython.Build import cythonize


setup(
    name='pasta',
    version='0.5 alpha',
    packages=find_packages(),
    ext_modules=cythonize(
        [Extension('forstr', sources=['pasta/forstr.pyx'])]),
    install_requires=['numpy>=1.11.3', 'spglib>=1.9.9.25'],
    author='Yunxing Zuo, WeiJi Hsiao, and Jianshu Jie',
    maintainer='WeiJi Hsiao and Jianshu Jie'
)
