from setuptools import setup, find_packages

setup(
    name='your_project',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'numpy==1.23.5',
        'pandas==1.5.3',
        'pyteomics==4.6',
        'setuptools==65.6.3',
    ],
    scripts=['scripts/utils.py', 'scripts/fipy.py'],  # Add your scripts here
    # Or, if you want the scripts to be installed as CLI commands, use entry_points
    entry_points={
        'console_scripts': [
            'fipy = scripts.fipy:main',  # Assuming `main()` is in fipy.py
            'utils = scripts.utils:main',  # If you have a main() in utils.py
        ],
    },
)
