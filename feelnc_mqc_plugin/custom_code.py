#!/usr/bin/env python
""" feelnc_mqc_plugin functions

We can add any custom Python functions here and call them
using the setuptools plugin hooks.
"""

from __future__ import print_function
from pkg_resources import get_distribution
import logging

from multiqc.utils import report, util_functions, config

# Initialise the main MultiQC logger
log = logging.getLogger('multiqc')

# Save this plugin's version number (defined in setup.py) to the MultiQC config
config.feelnc_mqc_plugin_version = get_distribution("feelnc_mqc_plugin").version

# Add default config options for the things that are used in MultiQC_NGI
def feelnc_mqc_plugin_execution_start():
    """ Code to execute after the config files and
    command line flags have been parsedself.

    This setuptools hook is the earliest that will be able
    to use custom command line flags.
    """

    log.info("Running feelnc plugin v{}".format(config.feelnc_mqc_plugin_version))

    # Add to the main MultiQC config object.
    # User config files have already been loaded at this point
    #   so we check whether the value is already set. This is to avoid
    #   clobbering values that have been customised by users.

    # config.ignore_images="false"

    if 'feelnc/rf_summary' not in config.sp:
        config.update_dict( config.sp, { 'feelnc/rf_summary': { 'fn': '*RF_summary.txt' } } )
    if 'feelnc/lnc_classes' not in config.sp:
        config.update_dict( config.sp, { 'feelnc/lnc_classes': { 'fn': '*classes.txt' } } )
    if 'feelnc/roc' not in config.sp:
        config.update_dict( config.sp, { 'feelnc/roc': { 'fn': '*_TGROC.png' } } )
    if 'feelnc/lnc_classes_log' not in config.sp:
        config.update_dict( config.sp, { 'feelnc/lnc_classes_log': { 'fn': '*feelncclassifier.log' } } )
    if 'feelnc/classification_summary' not in config.sp:
        config.update_dict( config.sp, { 'feelnc/classification_summary': { 'fn': 'feelnc_classification_summary.txt' } } )

    if 'feelnc/filter_log' not in config.sp:
        config.update_dict( config.sp, { 'feelnc/filter_log': { 'fn': '*feelncfilter.log' } } )
    # # Some additional filename cleaning
    # config.fn_clean_exts.extend([
    #     '.my_tool_extension',
    #     '.removeMetoo'
    # ])

    # # Ignore some files generated by the custom pipeline
    # config.fn_ignore_paths.extend([
    #     '*/my_awesome_pipeline/fake_news/*',
    #     '*/my_awesome_pipeline/red_herrings/*',
    #     '*/my_awesome_pipeline/noisy_data/*',
    #     '*/my_awesome_pipeline/rubbish/*'
    # ])
