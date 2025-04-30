from setuptools import setup, find_packages

setup(
    name="pynnotate",  
    version="0.1.0", 
    description="A script for downloading and processing GenBank features.",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/fernandacaron/pynnotate",  
    packages=find_packages(),
    py_modules=["pynnotate"],
    install_requires=[
        "biopython",
        "pandas",
    ],
    entry_points={
        "console_scripts": [
            "pynnotate=pynnotate.main:main",
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",
)
