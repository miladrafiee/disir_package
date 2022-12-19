#!/usr/bin/env python3d

from setuptools import setup

setup(
    name='disir_package',
    version='0.0.1',
    packages=['disir_package',
              'disir_package.common'],
    package_dir={'disir_package': 'disir_package'},
    entry_points={'console_scripts':
                      ['disir_package=disir_package.disir_package:run_disir_package']
                  }
)

