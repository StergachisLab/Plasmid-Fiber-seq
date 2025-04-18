from setuptools import setup, find_packages

setup(
    name="footprint_analysis",
    version="1.0.0",
    packages=find_packages(),
    description="A package for SNP footprint analysis with multiprocessing support",
    author="Your Name",
    author_email="your.email@example.com",
    install_requires=[
        "numpy",
        "pandas",
        "matplotlib",
        "seaborn",
        "scipy",
    ],
    python_requires=">=3.6",
    entry_points={
        "console_scripts": [
            "footprint-analysis=footprint_analysis.main:main",
        ],
    },
)
