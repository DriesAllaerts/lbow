from setuptools import find_packages, setup

with open('requirements.txt') as f:
    REQUIREMENTS = f.read().splitlines()

setup(
    name='lbow',
    version='1.0',
    description='Python package for solving linear buoyancy waves problems',
    author='Dries Allaerts',
    license='',
    packages=find_packages(include=['lbow']),
    install_requires=REQUIREMENTS,
    extras_require={
        'interactive': ['jupyter', 'matplotlib']
        },
    #setup_requires=['pytest-runner', 'flake8'],
    #tests_require=['pytest']
)
