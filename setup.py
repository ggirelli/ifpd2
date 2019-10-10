"""A setuptools based setup module.

See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages

# Get the long description from the README file
with open('README.md', "r") as f:
	long_description = f.read()

setup(name='ifpd2',
	version='0.0.1',
	description='''An iFISH probe design pipeline (v2).''',
	long_description=long_description,
	long_description_content_type='text/markdown',
	url='https://github.com/ggirelli/ifpd2',
	author='Gabriele Girelli',
	author_email='gabriele.girelli@scilifelab.se',
	license='MIT',
	classifiers=[
		'Development Status :: 3 - Alpha',
		'Intended Audience :: Science/Research',
		'Topic :: Scientific/Engineering :: Bio-Informatics',
		'License :: OSI Approved :: MIT License',
		'Programming Language :: Python :: 3 :: Only',
		'Programming Language :: Python :: 3.6'
	],
	keywords='biology cell DNA RNA FISH fluorescence hybridization bioimaging genome',
	packages=["ifpd2"],
	install_requires=[
		'numpy>=1.15.4',
		'pandas>=0.23.4',
		'tqdm>=4.19.8'
	],
	scripts=[
		"bin/ifpd2_query",
		"bin/ifpd2_mkbindb",
		"bin/ifpd2_dbchk"
	],
	test_suite="nose.collector",
	tests_require=["nose"],
)
