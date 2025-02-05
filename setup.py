from setuptools import setup, find_packages

setup(
    name='fipy',
    version='0.1',
    packages=find_packages(),  # This will find the 'fipy' package
    install_requires=[
        'numpy==1.23.5',
        'pandas==1.5.3',
        'pyteomics==4.6',
        'setuptools==65.6.3',
    ],
    # Define CLI commands using entry_points.
    entry_points={
        'console_scripts': [
            'fipy = scripts.fipy:main',   
            'utils = scripts.utils:main', 
        ],
    },
    include_package_data=True,  # Ensures non-code files are included
    package_data={
        # Include any CSV files in the fipy/data directory.
        'fipy': ['data/*.csv'],
    },
)
