#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  4 14:19:42 2017

@author: Sebastijan Mrak <smrak@gmail.com>
"""

from setuptools import setup


setup(name='gpstec',
      description='Utility functions for madrigal GPS TEC maps',
      author='Sebastijan Mrak',
      url='https://github.com/aldebaran1/gpstec.git',
      install_requires=['madrigalWeb'],
      packages=['gpstec']

)