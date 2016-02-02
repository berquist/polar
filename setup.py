from setuptools import setup

setup(name='polar',
      version='0.1',
      description='Code for polarizability with finite field method',
      url='http://github.com/fhqgfss/polar',
      author='Yilin Zhao',
      author_email='zhao229@mcmaster.ca',
      license='MIT',
      packages=['polar'],
      test_suite='nose.collector',
      install_requires=[
        'numpy',
        'nose',
      ],
      zip_safe=False)
