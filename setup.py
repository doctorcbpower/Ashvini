import setuptools

setuptools.setup(
    name="Ashvini",
    version="0.1",
    author="Anand Menon",
    #    author_email="",
    description="light-weight galaxy formation and evolution model",
    packages=setuptools.find_packages(),
    python_requires=">=3.6",
    install_requires=[
        "numpy",
        "scipy",
        "astropy",
        "h5py",
        "joblib",
        "tqdm",
        "pyyaml",
    ],
)
