from setuptools import setup, find_packages

setup(
    name="annopython",  # Nome do pacote
    version="0.1",  # Versão inicial
    description="A script for downloading and processing GenBank features.",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/fernandacaron/annopython",  
    packages=find_packages(),
    py_modules=["annopython"],
    install_requires=[
        "biopython",
        "pandas",
    ],
    entry_points={
        "console_scripts": [
            "annopython=annopython.main:main",
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",
)
