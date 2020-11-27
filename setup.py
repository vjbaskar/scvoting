import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="scvoting-vjbaskar", # Replace with your own username
    version="0.0.2",
    author="Vijay",
    author_email="vjbaskar@gmail.com",
    description="Projection of on scRNAseq on another and voting. Sam Watcham's thesis chapter 3",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/pypa/scvoting",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
