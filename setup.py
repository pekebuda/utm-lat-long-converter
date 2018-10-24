from setuptools import setup, find_packages

setup(
    name='geocon',
    version='0.0.1',
    description='A UTM -> Lat/Long (or vice versa) converter.',
    # long_description=open('README.rst').read(),
    url='https://github.com/pekebuda/utm-lat-long-converter-python',
    author='Fernando Diaz Laclaustra',
    author_email='peke_buda@hotmail.com',
    license='MIT',
    keywords="hello world insanely useful",
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
    ],
    packages=find_packages(),
    install_requires=[],
    entry_points={
        'console_scripts': [],  # <name>=<package>:<function>
    }
)
