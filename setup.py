import setuptools
with open("README.md", "r") as fh:
    long_description = fh.read()
setuptools.setup(
    name="barseqcount",
    version="0.0.5",
    author="Damien Marsic",
    author_email="damien.marsic@aliyun.com",
    description="Analysis of DNA barcode sequencing experiments",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/damienmarsic/barseqcount",
    package_dir={'': 'barseqcount'},
    packages=setuptools.find_packages(),
    py_modules=["barseqcount"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires='>=3.6',
)
