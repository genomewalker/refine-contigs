from setuptools import setup
import versioneer

requirements = [
    "pandas>=1.2.0",
    "scipy>=1.5.2",
    "networkx>=2.5",
    "tqdm==4.50.0",
    "PyYAML>=5.4",
    "biopython>=1.68",
    "pyfaidx>=0.5.9",
    "pyranges>=0.0.97",
]

setup(
    name="refine-contigs",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description="A simple contig merging tool using minimus2",
    license="GNUv3",
    author="Antonio Fernandez-Guerra",
    author_email="antonio@metagenomics.eu",
    url="https://github.com/genomewalker/refine-contigs",
    packages=["refine_contigs"],
    entry_points={"console_scripts": ["refineC=refine_contigs.__main__:main"]},
    install_requires=requirements,
    keywords="refine-contigs",
    classifiers=[
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
    ],
)
