
from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='py-adhore',
    version='1.0',
    packages=['src'],
    license='GPL',
    author='Arthur Zwanepoel',
    author_email='arzwa@psb.vib-ugent.be',
    description='py-adhore',
    py_modules=['pyadhore'],
    include_package_data=True,
    install_requires=[
        'click>=7.0',
        'coloredlogs>=10.0',
        'numpy>=1.16',
        'matplotlib>=3.0.2',
        'pandas==0.24.1',
    ],
    entry_points='''
        [console_scripts]
        py-adhore=pyadhore:cli
    ''',
)
