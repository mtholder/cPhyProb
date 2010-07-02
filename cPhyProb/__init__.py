#! /usr/bin/env python
# Copyright (c) 2005-7 by Mark T. Holder,  University of Kansas
# (see bottom of file)
"""Library for probability calculations used in the calculate likelihoods on
phylogenies.

Computationally intensive tasks are implemented as C-extensions.
"""
__all__ = [
    "ccore",
    "discrete_model",
    "discrete_char_type",
    "tests",
    ]

import logging, os

class ComputationalResourceFlags:
    """Enumeration of flags for describing computational resources."""
    # must be kept in sync with BeagleFlags in beagle.h!
    SINGLE  = 1L << 0    #/**< single precision computation */
    DOUBLE  = 1L << 1    #/**< double precision computation */
    SYNCH   = 1L << 2    #/**< synchronous computation */
    ASYNCH  = 1L << 3    #/**< asynchronous computation */
    REAL_EIGEN = 1L << 4    #/**< complex eigenvalue computation */
    COMPLEX = 1L << 5    #/**< complex eigenvalue computation */

    SCALING_MANUAL = 1 << 6  #< Manual scaling */
    SCALING_AUTO = 1 << 7 # Auto-scaling on */
    SCALING_ALWAYS = 1 << 8 # Scale at every updatePartials */

    RAW_SCALER = 1L << 9    #/**< save log scalers */
    LOG_SCALER = 1L << 10    #/**< save log scalers */

    SSE     = 1L << 11   # SSE */
    VECTOR_NONE = 1 << 12 # no sse

    OPENMP  = 1L << 13   # OpenMP threading */
    THREADING_NONE = 1L << 14 #

    CPU     = 1L << 15   #/**< CPU */
    GPU     = 1L << 16   #/**< GPU */
    FPGA    = 1L << 17   #/**< FPGA */
    CELL    = 1L << 18   # Cell */
    def to_str_list(x):
        s = []
        if x & ComputationalResourceFlags.DOUBLE:
            s.append("DOUBLE")
        if x & ComputationalResourceFlags.SINGLE:
            s.append("SINGLE")
        if x & ComputationalResourceFlags.ASYNCH:
            s.append("ASYNCH")
        if x & ComputationalResourceFlags.SYNCH:
            s.append("SYNCH")
        if x & ComputationalResourceFlags.COMPLEX:
            s.append("COMPLEX")
        if x & ComputationalResourceFlags.LOG_SCALER:
            s.append("LOG_SCALER")
        if x & ComputationalResourceFlags.CPU:
            s.append("CPU")
        if x & ComputationalResourceFlags.GPU:
            s.append("GPU")
        if x & ComputationalResourceFlags.FPGA:
            s.append("FPGA")
        if x & ComputationalResourceFlags.SSE:
            s.append("SSE")
        if x & ComputationalResourceFlags.CELL:
            s.append("CELL")
        if x & ComputationalResourceFlags.OPENMP:
            s.append("OPENMP")
        return s
    to_str_list = staticmethod(to_str_list)
    def resource_descrip_to_str(t):
        return "(%s, %s, [%s], [%s])" % (t[0], t[1],
                                       ComputationalResourceFlags.to_str_list(t[2]),
                                       ComputationalResourceFlags.to_str_list(t[3]))
    resource_descrip_to_str = staticmethod(resource_descrip_to_str)



###############################################################################
## LOGGING

_LOGGING_LEVEL_ENVAR="C_PHY_PROB_LOGGING_LEVEL"
_LOGGING_FORMAT_ENVAR="C_PHY_PROB_LOGGING_FORMAT"

def get_logging_level():
    if _LOGGING_LEVEL_ENVAR in os.environ:
        if os.environ[_LOGGING_LEVEL_ENVAR].upper() == "NOTSET":
            level = logging.NOTSET
        elif os.environ[_LOGGING_LEVEL_ENVAR].upper() == "DEBUG":
            level = logging.DEBUG
        elif os.environ[_LOGGING_LEVEL_ENVAR].upper() == "INFO":
            level = logging.INFO
        elif os.environ[_LOGGING_LEVEL_ENVAR].upper() == "WARNING":
            level = logging.WARNING
        elif os.environ[_LOGGING_LEVEL_ENVAR].upper() == "ERROR":
            level = logging.ERROR
        elif os.environ[_LOGGING_LEVEL_ENVAR].upper() == "CRITICAL":
            level = logging.CRITICAL
        else:
            level = logging.NOTSET
    else:
        level = logging.NOTSET
    return level

def get_logger(name="C_PHY_PROB"):
    """
    Returns a logger with name set as given, and configured
    to the level given by the environment variable _LOGGING_LEVEL_ENVAR.
    """
    logger_set = False
#     package_dir = os.path.dirname(module_path)
#     config_filepath = os.path.join(package_dir, _LOGGING_CONFIG_FILE)
#     if os.path.exists(config_filepath):
#         try:
#             logging.config.fileConfig(config_filepath)
#             logger_set = True
#         except:
#             logger_set = False
    logger = logging.getLogger(name)
    if not logger_set:
        level = get_logging_level()
        rich_formatter = logging.Formatter("[%(asctime)s] %(filename)s (%(lineno)d): %(levelname) 8s: %(message)s")
        simple_formatter = logging.Formatter("%(levelname) 8s: %(message)s")
        raw_formatter = logging.Formatter("%(message)s")
        default_formatter = None
        logging_formatter = default_formatter
        if _LOGGING_FORMAT_ENVAR in os.environ:
            if os.environ[_LOGGING_FORMAT_ENVAR].upper() == "RICH":
                logging_formatter = rich_formatter
            elif os.environ[_LOGGING_FORMAT_ENVAR].upper() == "SIMPLE":
                logging_formatter = simple_formatter
            elif os.environ[_LOGGING_FORMAT_ENVAR].upper() == "NONE":
                logging_formatter = None
            else:
                logging_formatter = default_formatter
        else:
            logging_formatter = default_formatter
        if logging_formatter is not None:
            logging_formatter.datefmt='%H:%M:%S'
        logger.setLevel(level)
        ch = logging.StreamHandler()
        ch.setLevel(level)
        ch.setFormatter(logging_formatter)
        logger.addHandler(ch)
    return logger
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
