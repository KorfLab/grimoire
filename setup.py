import setuptools

with open('README.md', 'r') as fh:
	long_description = fh.read()

setuptools.setup(
	name='grimoire',
	version='0.0.1',
	author='Ian Korf',
	author_email='ifkorf@gmail.com',
	description='Biological sequence analysis tools',
	long_description=long_description,
	long_description_content_type='text/markdown',
	url='https://github.com/KorfLab/grimoire',
	packages=setuptools.find_packages(),
	classifiers=[
		'Programming Language :: Python :: 3',
		'License :: OSI Approved :: MIT License',
		'Operating System :: OS Independent',
	],
)

