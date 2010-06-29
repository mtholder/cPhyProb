#! /usr/bin/env python
# Copyright (c) 2005-7 by Mark T. Holder,  University of Kansas
# (see bottom of file)
"""Classes for encapsulating datatypes, and factory for encapsulating data 
structures for characters at the leaves of a tree and the conditional 
likelihood array for a tree.
"""

from cPhyProb.ccore.dsct_model import cleaf_data_ctor, cla_ctor, full_la_ctor, \
    csslookup_ctor

_DNA_TYPE = None

def DNAType():
    global _DNA_TYPE
    if _DNA_TYPE is None:
        ambig_codes = (("R", "AG"),
                       ("Y","CT"),
                       ("M","AC"),
                       ("W","AT"),
                       ("S","CG"),
                       ("K","GT"),
                       ("V","ACG"),
                       ("H","ACT"),
                       ("D","AGT"),
                       ("B","CGT"),)
        aliases = (("-", "N"),("X", "N"), ("?", "N"))
        _DNA_TYPE = DiscreteCharType("ACGT", 
                                     missing="N", 
                                     ambig_codes=ambig_codes,
                                     aliases=aliases)
    return _DNA_TYPE

class DiscreteCharType(object):
    def __init__(self, states, ambig_codes=(), aliases=None, ignore_case=True, missing=None):
        """Creates an immutable DiscreteCharType object from the list of
        `states`.

        The states list determines the ordinate of states.  `states` must be an
        iterable collection of strings.
        `ambig_codes` must be a iterable collection of key-value pairs. Where the
            keys are new state symbols or state labels, and the values or 
            lists of previously declared states.
        `aliases` is an iterable collection of pairs of strings:
            (alias, real-name).
        If `missing` is specified, it is treated as the first ambig_code
        """
        self._n_states = len(states)
        self._states = tuple(states)
        self._ignore_case = ignore_case
        if self._n_states < 2:
            raise ValueError("The number of states must be greater than 1")
        labels = []
        expansions = []
        symbols_to_state_sets = {}
        symbol_to_ind = {}
        for n, s in enumerate(states):
            if ignore_case:
                s = s.upper()
            if s in symbols_to_state_sets:
                m = "The state %s occurred twice is the state list" % s
                raise ValueError(m)
            
            symbols_to_state_sets[s] = (n,)
            symbol_to_ind[s] = len(labels)
            labels.append(s)
            expansions.append((n,))
        alias_list = aliases and list(aliases) or []
        if missing is not None:
            if missing in symbols_to_state_sets:
                m = "The missing state label (%s) was declared as a state" % missing
                raise ValueError(m)
            ac = [(missing, states)]
            ac.extend(ambig_codes)
        else:
            ac = ambig_codes
        for key, val in ac:
            if ignore_case:
                key = key.upper()
            if key in symbols_to_state_sets:
                m = "The ambiguity code (%s) was declared as a state" % key
                raise ValueError(m)
            expanded = set()
            for v_el in val:
                if ignore_case:
                    v_el = v_el.upper()
                try:
                    e_el = symbols_to_state_sets[v_el]
                    for e in e_el:
                        expanded.add(e)
                except KeyError:
                    raise ValueError("The ambiguity code expansion (for code %s) contains an unknown label (%s)" % (key, v_el))
            exp_list = list(expanded)
            exp_list.sort()
            texp = tuple(exp_list)
            if texp in expansions:
                i = expansions.index(texp)
                alias_list.append((key, labels[i]))
            else:
                symbols_to_state_sets[key] = texp
                symbol_to_ind[key] = len(labels)
                labels.append(key)
                expansions.append(texp)
        for key, val in alias_list:
            if ignore_case:
                key = key.upper()
                val = val.upper()
            if key in symbols_to_state_sets:
                m = "The alias (%s) was declared as a state" % key
                raise ValueError(m)
            if val not in symbols_to_state_sets:
                m = "The alias %s maps to an unknown state (%s)" % (key, val)
                raise ValueError(m)
            symbol_to_ind[key] = symbol_to_ind[val]
        self._symbol_to_ind = symbol_to_ind
        self._all_symbols = tuple(labels)
        self._state_code_lookup = tuple(expansions)
        expansions = []
        for i in self._state_code_lookup:
            expansions.append(list(i))
        self._csslookup = csslookup_ctor(self._n_states, expansions)

    def to_indices(self, seq):
        "Converts a sequence of character symbols to the corresponding indices."
        l = []
        sti = self._symbol_to_ind
        for s in seq:
            l.append(sti[s])
        return l

    def get_n_states(self):
        return self._n_states
    n_states = property(get_n_states)

    def get_states(self):
        return self._states
    states = property(get_states)

    def get_ignore_case(self):
        return self._ignore_case
    ignore_case = property(get_ignore_case)

    def get_symbol_to_ind(self):
        return dict(self._symbol_to_ind)
    symbol_to_ind = property(get_symbol_to_ind)

    def get_all_symbols(self):
        return self._all_symbols
    all_symbols = property(get_all_symbols)

    def get_state_code_lookup(self):
        return self._state_code_lookup
    state_code_lookup = property(get_state_code_lookup)

    def get_csslookup(self):
        return  self._csslookup
    csslookup = property(get_csslookup)

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
