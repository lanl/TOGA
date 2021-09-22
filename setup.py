

import setuptools

with open('README.md', 'r', encoding='utf-8') as fh:
    long_description = fh.read()

setuptools.setup(
    name='togamgxs',
    version='0.0.1',
    description='Monte Carlo MGXS generation and optimization',
    long_description=long_description,
    package_dir={"":"src"},
    packages=setuptools.find_packages(where='src'),
    python_requres=">=3.6"
)
