from setuptools import setup, find_packages
import os


def read_file(filename):
    if os.path.exists(filename):
        with open(filename, "r", encoding="utf-8") as f:
            return f.read()
    return ""


setup(
    name="placeholder_name",
    version="0.0.1",
    packages=find_packages(),
    package_data={
        "placeholder_name": [
            "datafiles/name_dicts/*.json",
        ],
    },
    include_package_data=True,
    zip_safe=False,
    install_requires=read_file("requirements.txt").splitlines(),
    author="De Novo Chem Team",
    author_email="carson.britt@denovochem.com",
    description="Name-to-SMILES conversion",
    long_description=read_file("README.md"),
    long_description_content_type="text/markdown",
    url="https://github.com/denovochem/name_to_smiles",
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Topic :: Software Development :: Libraries",
    ],
    python_requires=">=3.10",
)
