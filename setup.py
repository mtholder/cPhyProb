#! /usr/bin/env python
# Copyright (c) 2005-7 by Mark T. Holder,  University of Kansas
# (see bottom of file)
import ez_setup
import sys
ez_setup.use_setuptools()

from setuptools import setup, find_packages, Extension

args = sys.argv
help_msg = """
cPhyProb specific options:

Compilation flags
    --debug         compile in debug mode
    --no-inline     use this if your compiler does not deal with the "inline"
                    keyword.
Control how rescaling to avoid underflow is implemented using one of:
    --const-step-rescaling      The multiplier is a compile time constant
    --int-incr-rescaling        The multiplier is the integer close to the real
                                    number that would rescale the max. cond.
                                    likelihood to 1.0.
    (no arg)                    The default is to rescale with a real number 
                                    that makes the max. cond. likelihood 1.0.
"""
if ("--help" in args) or ("-h" in args):
    sys.stderr.write(help_msg)
if "--debug" in args:
    kDebugPrint = "1"
    extra_compile_args = ["-Wall"]
else:
    kDebugPrint = "0"
    extra_compile_args = []
if "--no-inline" in args:
    no_inline_val = 1
    args.remove("--no-inline")
else:
    no_inline_val = 0
    
if "--const-step-rescaling" in args:
    rescale_type_define = "CONSTANT_STEP_RESCALING"
    args.remove("--const-step-rescaling")
elif "--int-incr-rescaling" in args:
    rescale_type_define = "INT_INCR_RESCALING"
    args.remove("--int-incr-rescaling")
else:
    rescale_type_define = "TO_EXACTLY_ONE_RESCALING"
    
setup(name = "cPhyProb", 
      version = "0.01",
      packages = find_packages(),
      maintainer = "Mark Holder", 
      maintainer_email = "mtholder@gmail.com", 
      description = "Calculations of likelihoods for phylogenetics using C extensions", 
      test_suite = "cPhyProb.tests",
      ext_modules = [
        Extension("cPhyProb.ccore.dsct_model", 
            define_macros=[("PRINTING_LOTS", kDebugPrint), 
                           (rescale_type_define, 1),
                           ("BUILDING_FOR_PYTHON",1),
                           ("NO_INLINE", no_inline_val),
                          ],
            extra_compile_args=extra_compile_args,
            sources=[
                "cPhyProb/ccore/dsct_model.c",
                "cPhyProb/ccore/py_dsct_model.c"
                ])],
      license = "GNU General Public License (see LICENCE.txt and GPL.txt)",
      keywords = "phylogenetics bioinformatics CIPRES",
      url = "http://people.ku.edu/~mtholder/cPhyProb",
      classifiers = [
            "Development Status :: 2 - Pre-Alpha",
            "Environment :: Console",
            "Intended Audience :: Developers",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: BSD License",
            "License :: OSI Approved :: GNU Library or  General Public License (GPL)",
            "Natural Language :: English",
            "Operating System :: OS Independent",
            "Programming Language :: Python",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
            ],
    )

################################################################################
# cPhyProb is a package implementing some probability calculations used in
#   calculating likelihoods on phylogenies. 
# 
# Copyright (C) 2005-2007  Mark Holder mtholder@gmail.com
# 
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU  General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
# 
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  General Public License for more details.
# 
# You should have received a copy of the GNU  General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
################################################################################
