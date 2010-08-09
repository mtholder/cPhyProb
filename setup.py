#! /usr/bin/env python
# Copyright (c) 2005-7 by Mark T. Holder,  University of Kansas
# (see bottom of file)
import ez_setup
import sys, os
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
if 'CFLAGS' in os.environ:
    extra_compile_args.extend(os.environ['CFLAGS'].split())

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


include_dirs = []
library_dirs = []
libraries = []
preprocessor_defines =[("PRINTING_LOTS", kDebugPrint),
                       (rescale_type_define, 1),
                       ("BUILDING_FOR_PYTHON",1),
                       ("NO_INLINE", no_inline_val),
                      ]
if "--use-beagle" in args:
    args.remove("--use-beagle")
    beagle_full_version = "1.0.1"
    beagle_version = ".".join(beagle_full_version.split(".")[:2])
    beagle_lib_dir = os.environ.get('BEAGLE_LIB_DIR')
    if not beagle_lib_dir:
        beagle_lib_dir = os.path.join("usr", "local")
        if not os.path.exists(beagle_lib_dir):
            sys.exit("BEAGLE_LIB_DIR is not defined in the env, and the default path /usr/local was not found!")
    elif not os.path.exists(beagle_lib_dir):
        sys.exit("Could not find the directory %s" % beagle_lib_dir)

    beagle_lib_inc_dir = os.path.join(beagle_lib_dir, "include", "libhmsbeagle-%s" % beagle_version)
    if not os.path.exists(beagle_lib_inc_dir):
        sys.exit("Could not find the directory %s" % beagle_lib_inc_dir)

    beagle_lib_lib_dir = os.path.join(beagle_lib_dir, "lib")
    if not os.path.exists(beagle_lib_lib_dir):
        sys.exit("Could not find the directory %s" % kbeagle_lib_lib_dir)
    include_dirs.append(beagle_lib_inc_dir)
    library_dirs.append(beagle_lib_lib_dir)
    cuda_lib_dir = os.environ.get('BEAGLE_DEPENDENCY_LIB_DIR')
    if cuda_lib_dir:
        library_dirs.append(cuda_lib_dir)
        cuda_lib = os.environ.get('BEAGLE_DEPENDENCY_LIB')
        if cuda_lib:
            libraries.append(cuda_lib)
    libraries.append('hmsbeagle-%s' % beagle_full_version)
    preprocessor_defines.append(('USE_BEAGLE_LIB', 1))

using_ncl = "--use-ncl" in args
if using_ncl:
    args.remove("--use-ncl")
    ncl_install_dir = os.environ.get('NCL_INSTALL_DIR')
    if not ncl_install_dir:
        ncl_install_dir = os.path.join("usr", "local")
        if not os.path.exists(ncl_install_dir):
            sys.exit("NCL_INSTALL_DIR is not defined in the env, and the default path /usr/local was not found!")
    elif not os.path.exists(ncl_install_dir):
        sys.exit("Could not find the directory %s" % ncl_install_dir)

    ncl_inc_dir = os.path.join(ncl_install_dir, "include")
    if not os.path.exists(ncl_inc_dir):
        sys.exit("Could not find the directory %s" % ncl_inc_dir)

    ncl_lib_dir = os.path.join(ncl_install_dir, "lib", "ncl")
    if not os.path.exists(ncl_lib_dir):
        sys.exit("Could not find the directory %s" % ncl_lib_dir)
    include_dirs.append(ncl_inc_dir)
    library_dirs.append(ncl_lib_dir)
    libraries.append('ncl')

ext_sources = [ #"cPhyProb/ccore/dsct_model.c",
            "cPhyProb/ccore/asrv.c",
            "cPhyProb/ccore/beagle_wrap.c",
            "cPhyProb/ccore/cphyprob_inlines.c",
            "cPhyProb/ccore/cphyprob.c",
            "cPhyProb/ccore/non_beagle_impl.c",
            "cPhyProb/ccore/phylo_util.c",
            "cPhyProb/ccore/state_set_lookup.c",
          ]
if using_ncl:
    ext_sources.append("cPhyProb/ccore/from_ncl_factory.c")
setup(name = "cPhyProb",
      version = "0.01",
      packages = find_packages(),
      maintainer = "Mark Holder",
      maintainer_email = "mtholder@gmail.com",
      description = "Calculations of likelihoods for phylogenetics using C extensions",
      test_suite = "cPhyProb.tests",
      ext_modules = [
        Extension("cPhyProb.ccore.dsct_model",
            define_macros=preprocessor_defines,
            include_dirs=include_dirs,
            library_dirs=library_dirs,
            libraries=libraries,
            extra_compile_args=extra_compile_args,
            sources=ext_sources)],
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
