"""Identity of a ground-truth sample set.

A ground-truth value is y = QoI(u(z)) for ONE exact parametric problem, so two
sets are interchangeable only if every field in IDENTITY agrees. Two guards,
both cheap:

  * the folder name carries an 8-character hash of IDENTITY, so sets that differ
    in any field cannot land in the same directory (and a stale directory can
    never be resumed into or merged with a different problem);
  * each folder holds a GT_meta.json spelling the fields out, together with the
    exact string that was hashed -- so nobody has to reverse-engineer the hash.

Deliberately NOT in IDENTITY:
  nb_samples  a set is allowed to grow; the draw is prefix-stable, which
              tests/test_gt_meta.py checks
  solver      mixing solvers is fine PROVIDED each sample converged: the linear
  rtol/atol   solve is orders of magnitude below the discretisation error. They
              are recorded as a quality claim, not as identity.

Note the identity field is the number of parameters d, which happens to be 2*nb_cosines on many examples provided.

WARNING: This was only partially reviewed after vibecoding. 
"""
import hashlib
import json
import os

__author__ = ["Jean-Luc Bouchot"]
__copyright__ = "Copyright 2019 - 2026, INRIA, LMU Munisch, Chair C for Mathematics (Analysis), RWTH Aachen and Seminar for Applied Mathematics, ETH Zurich and School of Mathematics and Statistics, Beijing Institute of Technology"
__credits__ = ["Jean-Luc Bouchot", "Benjamin, Bykowski", "Falk Pulsmeyer", "Holger Rauhut", "Christoph Schwab"]
__license__ = "GPL"
__version__ = "0.1.0-dev"
__maintainer__ = "Jean-Luc Bouchot"
__email__ = "jlbouchot@gmail.com"
__status__ = "Development"
__created__ = "2026/11/03"
__lastmodified__ = "2026/11/07"

## IDENTITY is the hash contract: explicit and ordered, because changing it
## changes every folder name. The two tables below say where each field comes
## from, and _check_tables() below asserts they still describe IDENTITY, so a
## new coefficient parameter cannot be added to one place and forgotten in the
## others.
IDENTITY = ("d", "power", "weight_cosine", "fluctuation_importance", "abar",
            "rhs", "qoi", "mesh", "degree", "seed")

## Read from the experiment's .ini, with the section they live in. nb_cosines
## is here but not in IDENTITY: it is what d is derived from (d = 2*nb_cosines).
CFG_FIELDS = {"nb_cosines": "main", "power": "main", "weight_cosine": "main",
              "fluctuation_importance": "main", "abar": "main", "qoi": "main",
              "degree": "pde"}

## Properties of the reference itself, not of the experiment: command line only.
GT_ONLY = ("mesh", "seed", "rhs")

META = "GT_meta.json"

## Everything in an .ini is text, so each field says how to read it: a value
## that came from a file must hash the same as the same value typed on the
## command line.
CAST = {"nb_cosines": int, "power": float, "weight_cosine": float,
        "fluctuation_importance": float, "abar": float, "qoi": str,
        "degree": int}

def _check_tables():
    """IDENTITY must be exactly: d, plus the cfg fields other than nb_cosines,
    plus the GT-only fields."""
    if set(CAST) != set(CFG_FIELDS):
        raise AssertionError("CAST and CFG_FIELDS disagree: %s"
                            % sorted(set(CAST) ^ set(CFG_FIELDS)))
    derived = {"d"} | (set(CFG_FIELDS) - {"nb_cosines"}) | set(GT_ONLY)
    if derived != set(IDENTITY):
        raise AssertionError(
            "gt_meta tables disagree with IDENTITY:\n  only in tables: %s\n"
            "  only in IDENTITY: %s"
            % (sorted(derived - set(IDENTITY)), sorted(set(IDENTITY) - derived)))
    return True

def check_parser(parser):
    """Every field a user may override must have a --flag. Catches the case
    where a parameter joins the identity but stays invisible on the CLI."""
    flags = set()
    for a in parser._actions:                      # noqa: SLF001 (argparse has no public view)
        flags.update(o.lstrip("-") for o in a.option_strings)
    want = set(CFG_FIELDS) | set(GT_ONLY)
    missing = sorted(want - flags)
    if missing:
        raise AssertionError("no command-line option for: %s" % ", ".join(missing))
    return True

def _fmt(v):
    """Stable text for a value: %.12g so 3.9605 does not hash differently from
    the same number read back from json."""
    return "%.12g" % v if isinstance(v, float) else str(v)

def canonical(identity):
    """The exact string that gets hashed."""
    missing = [k for k in IDENTITY if k not in identity]
    if missing:
        raise KeyError("ground-truth identity is missing %s" % ", ".join(missing))
    return ";".join("%s=%s" % (k, _fmt(identity[k])) for k in IDENTITY)

def gt_hash(identity):
    return hashlib.sha1(canonical(identity).encode("utf-8")).hexdigest()[:8]

def folder(identity):
    """Where this problem's samples belong, e.g. GTresults_d20H3000_1a2b3c4d."""
    return "GTresults_d%sH%s_%s" % (identity["d"], identity["mesh"], gt_hash(identity))

def write(path, identity, quality=None):
    """(Re)write path/GT_meta.json. ``quality`` is the convergence claim."""
    doc = {"identity": {k: identity[k] for k in IDENTITY},
           "canonical": canonical(identity),
           "hash": gt_hash(identity),
           "quality": quality or {}}
    with open(os.path.join(path, META), "w") as fh:
        json.dump(doc, fh, indent=2, sort_keys=True)
        fh.write("\n")
    return doc

def read(path):
    """The stored document, or None when there is no meta file."""
    f = os.path.join(path, META)
    if not os.path.exists(f):
        return None
    with open(f) as fh:
        return json.load(fh)

def check(path, identity):
    """Raise unless path holds samples for exactly this problem. Returns the
    stored document, or None if the folder is new (no meta yet)."""
    doc = read(path)
    if doc is None:
        return None
    if doc.get("hash") != gt_hash(identity):
        raise SystemExit(
            "%s holds a DIFFERENT problem.\n  stored:  %s\n  wanted:  %s\n"
            "Ground-truth values are not comparable across these; use the "
            "folder that %s names, or move this one aside."
            % (os.path.join(path, META), doc.get("canonical"),
               canonical(identity), folder(identity)))
    return doc

def draw(nb_samples, d, seed):
    """Z exactly as the drivers draw it: the Chebyshev preconditioning measure
    cos(pi U). Defined here, next to the seed it depends on, so the generating
    scripts and the graphing scripts cannot draw different points. Prefix-stable
    in nb_samples, which is what lets nb_samples stay out of the identity --
    tests/test_gt_meta.py checks that rather than assuming it.
    """
    import numpy as np
    rng = np.random.default_rng(seed)
    return np.cos(np.pi * rng.uniform(-1.0, 1.0, (nb_samples, d)))

def identity_from_args(args, mesh):
    """The identity of what an argparse namespace describes.

    The generating scripts come in through here and the graphing scripts
    through identity_from_cfg; tests/test_gt_meta.py asserts the two agree on
    the same .ini, which is the only reason a figure can trust a folder.
    """
    idy = dict((k, getattr(args, k)) for k in IDENTITY if k not in ("d", "mesh"))
    idy.update(d=2 * args.nb_cosines, mesh=mesh)
    canonical(idy)
    return idy

def identity_from_cfg(cfg_path, mesh, seed, rhs=10.0):
    """The identity of the ground truth that belongs to an experiment .ini.

    The graphing scripts call this instead of spelling a folder name out, so
    "which reference goes with which experiment" is answered in one place. The
    swept parameter still has to be set by the caller: a d-sweep needs one
    reference per d, so the call site does identity["d"] = 2 * d afterwards.
    """
    from utilities.MLCSPG_param_parser import (define_MLCSPG_cli_parser,
                                               parse_cfg_file_and_args)
    main_cfg, pde_cfg, _ = parse_cfg_file_and_args(
        define_MLCSPG_cli_parser().parse_args(["--cfg", cfg_path]))
    section = {"main": main_cfg, "pde": pde_cfg}
    vals = {}
    for k, sec in CFG_FIELDS.items():
        if k not in section[sec]:
            raise KeyError("%s does not define %s, which the ground truth depends on"
                           % (cfg_path, k))
        vals[k] = CAST[k](section[sec][k])
    idy = dict((k, v) for k, v in vals.items() if k != "nb_cosines")
    idy.update(d=2 * vals["nb_cosines"], mesh=mesh, seed=int(seed), rhs=float(rhs))
    canonical(idy)                 # raises now rather than at hashing time
    return idy
