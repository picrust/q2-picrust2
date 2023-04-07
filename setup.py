from setuptools import setup, find_packages

setup(
    name="q2-picrust2",
    version="2023.2",
    packages=find_packages(),
    package_data={'q2_picrust2': ['citations.bib']},
    author="Gavin Douglas",
    author_email="gavinmdouglas@gmail.com",
    description="Metagenome inference based on marker gene data.",
    license='GPL3',
    url="https://github.com/picrust/picrust2",
    entry_points={
        'qiime2.plugins': ['q2-picrust2=q2_picrust2.plugin_setup:plugin']
    },
    zip_safe=False,
)
