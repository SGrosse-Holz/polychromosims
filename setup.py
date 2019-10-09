from setuptools import setup, find_packages

with open('README.md') as f:
  readme = f.read()

setup(
    name='polychromosims',
    version='0.0.0',
    description='polychrom utilities',
    long_description=readme,
    author='Simon Grosse-Holz',
    packages=find_packages(exclude=('tests', 'docs')),
    )
