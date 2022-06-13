from setuptools import find_packages, setup

setup(
    packages=find_packages(),
    entry_points={"console_scripts": ["uvva_pipe=uvva.uvva_pipe:run_from_file"]},
)
