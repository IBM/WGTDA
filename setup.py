# Copyright 2023 IBM Corp. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

from setuptools import find_packages, setup

setup(
    name="wgtda",
    version="0.1.0",
    description="A python package to identify biomarkers in bulk RNAseq data",
    author="IBM Research Africa",
    author_email=["n.nyase@gmail.com", "lebohang.mashatola@ibm.com"],
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/IBM/wgtda",
    install_requires=[
        "numpy>=1.21.0",
        "pandas>=1.3.0",
        "dcor>=0.6",
        "scipy>=1.7.0",
        "matilda @ git+https://github.com/IBM/matilda.git",
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.9",
)
