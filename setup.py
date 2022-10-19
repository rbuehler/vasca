from setuptools import find_packages, setup

setup(
    packages=find_packages(),
    entry_points={"console_scripts": ["vasca_pipe=vasca.vasca_pipe:run_from_file"]},
)
