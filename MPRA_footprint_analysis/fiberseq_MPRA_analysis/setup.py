from setuptools import setup, find_packages

setup(
    name="footprint_analysis",
    version="1.0.0",
    packages=find_packages(),
    description="A package for SNP footprint analysis with multiprocessing support",
    author="Your Name",
    author_email="your.email@example.com",
    install_requires=[
        "numpy>=1.20.0",
        "pandas>=1.3.0",
        "matplotlib>=3.3.0",
        "seaborn>=0.11.0",
        "scipy>=1.7.0",
        "tqdm>=4.60.0",  # Required by FiberHMM_functions
    ],
    python_requires=">=3.6",
    entry_points={
        "console_scripts": [
            "footprint-analysis=footprint_analysis.main:main",
        ],
    },
)