#!/usr/bin/env python
"""
MultiQC plugin for feelnc, used in the TAGADA pipeline
For more information about Tagada, see https://github.com/FAANG/analysis-TAGADA
"""

from setuptools import setup, find_packages

version = '1.O'

setup(
    name = 'feelnc_mqc_plugin',
    version = version,
    author = 'Cervin Guyomar',
    author_email = 'cervin.guyomar@inrae.fr',
    description = "",
    long_description = __doc__,
    keywords = 'bioinformatics',
    url = '',
    download_url = '',
    license = '',
    packages = find_packages(),
    include_package_data = True,
    install_requires = [
        'multiqc'
    ],
    entry_points = {
        'multiqc.modules.v1': [
            'feelnc = feelnc_mqc_plugin.modules.feelnc:MultiqcModule',
        ],
        'multiqc.hooks.v1': [
            'execution_start = feelnc_mqc_plugin.custom_code:feelnc_mqc_plugin_execution_start'
        ]
    },
    classifiers = [
    ],
)
