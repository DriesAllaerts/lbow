from setuptools import find_packages, setup

with open('requirements.txt') as f:
    REQUIREMENTS = f.read().splitlines()

setup(
    name='pyslaw',
    version='0.1.0',
    description='Python project for solving linear atmospheric gravity waves problems',
    author='Dries Allaerts',
    license='',
    packages=find_packages(include=['pyslaw']),
    install_requires=REQUIREMENTS,
    extras_require={
        'interactive': ['jupyter', 'matplotlib']
        },
    #setup_requires=['pytest-runner', 'flake8'],
    #tests_require=['pytest']
)
