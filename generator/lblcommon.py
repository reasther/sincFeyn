"""
Shared emission helper for the four-loop lepton-loop generators: each node becomes its own C
function (see build.node_function).
"""
import os, build

def use_form():
    """Node polynomials as FORM Format O3 code, optimised jointly over each node's Wick blocks
    (build.form_optimise), when the FORM executable is found; GEN_FORM=0 or no FORM: the Python
    factored trees (compact.py then shares sub-sums across the blocks).  At four loops FORM gives the
    faster kernels since the line identifications (l4_c13 1.9x).  The Dirac tables are unaffected."""
    import formdirac
    build.USE_FORM = os.environ.get('GEN_FORM', '1') != '0' and formdirac.available() is not None
    build.FORM_JOINT = os.environ.get('GEN_FORM_JOINT', '1') != '0'
    if not build.USE_FORM: print('lepton-loop kernels: Python node polynomials (no FORM or GEN_FORM=0)', flush=True)
    return build.USE_FORM

def node_function(nd, tab, npar=None):
    return build.node_function(nd, tab)
