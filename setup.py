from setuptools import setup, find_packages

setup(
    name='zfel',
    version='0.1a',
    packages=find_packages(),
    url='https://github.com/slaclab/zfel',
    description='1-D FEL code',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    install_requires=[
        'matplotlib',
        'numpy',
        'scipy',
    ],
)
