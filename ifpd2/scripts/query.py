#!/usr/bin/python3
# -*- coding: utf-8 -*-

# ------------------------------------------------------------------------------
#
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Date: 2019-10-03
#
# ------------------------------------------------------------------------------

import argparse
import copy
import ifpd2 as fp
import logging
import os
import re
import sys

parser = argparse.ArgumentParser(
    description="""
Lorem ipsum dolor sit amet, consectetur adipisicing elit. Soluta, aspernatur,
natus. Possimus recusandae distinctio, voluptatem fuga delectus laudantium ut,
inventore culpa sit amet ullam officiis, tenetur nobis eius vitae dolore.
""",
    formatter_class=argparse.RawDescriptionHelpFormatter,
)

parser.add_argument(
    "region",
    type=str,
    help="""Region in "chrom:start-end" format. Any string is accepted for the
    chrom field. In case of a full chromosome search, omit start and end:
    "chrom:-".""",
)

parser.add_argument("dbpath", type=str, help="""Path to database file.""")

parser.add_argument(
    "-X",
    metavar="probes",
    type=int,
    default=None,
    help="""Number of probes to query for.""",
)
parser.add_argument(
    "-N",
    metavar="oligos",
    type=int,
    default=48,
    help="""Number of oligos per probe. Default: 48.""",
)
parser.add_argument(
    "-W", metavar="size", type=int, default=None, help="""Window size."""
)
parser.add_argument(
    "-w",
    metavar="shift",
    type=float,
    default=0.1,
    help="""Window shift as a percentage of window size. Default: 0.1""",
)
parser.add_argument(
    "-R",
    metavar="size",
    type=float,
    default=8000,
    help="""Central focus region size, either in nt or as fraction of window
    size. Default: 8000""",
)
parser.add_argument(
    "-r",
    metavar="step",
    type=float,
    default=1000,
    help="""Central focus region step, either in nt or as fraction of central
    focus region size. Default: 1000""",
)

parser.add_argument(
    "-F",
    metavar="nOT",
    type=int,
    nargs=2,
    default=(0, 100),
    help="""Acceptable range of off-targets. Default: (0,100)""",
)
parser.add_argument(
    "-G",
    metavar="dG",
    type=float,
    nargs=2,
    default=(0.0, 0.5),
    help="""Acceptable range of secondary structure delta free energy. Either
    as absolute kcal/mol os as a fraction of delta free energy of hybridization.
    Default: (.0,.5)""",
)
parser.add_argument(
    "-o",
    metavar="step",
    type=float,
    default=0.1,
    help="""Step for oligo score relaxation. Default: .1""",
)

parser.add_argument(
    "-D",
    metavar="dist",
    type=float,
    default=2,
    help="""Minimum distance between consecutive oligos. Default: 2""",
)
parser.add_argument(
    "-t",
    metavar="temp",
    type=float,
    default=10.0,
    help="""Melting temperature range half-width. Default: 10.""",
)
parser.add_argument(
    "-P",
    metavar="size",
    type=float,
    default=10000,
    help="""Maximum probe size in nt. Default: 10000""",
)
parser.add_argument(
    "-H",
    metavar="size",
    type=float,
    default=0.1,
    help="""Maximum hole size in a probe as fraction of probe size.
    Default: .1""",
)
parser.add_argument(
    "-I",
    metavar="oligos",
    type=float,
    default=0.5,
    help="""Probe oligo intersection threshold as shared oligo fraction.
    Default: .5""",
)

parser.add_argument(
    "-k",
    metavar="length",
    type=int,
    help="""Oligo length in nt. This should be used when all oligos in the
    database have the same length.""",
)

parser.add_argument(
    "-O",
    metavar="outpath",
    type=str,
    default=".",
    help="Path to output directory. Default to current folder.",
)

parser.add_argument(
    "--threads",
    metavar="t",
    type=int,
    default=1,
    help="Number of threads for parallelization. Default: 1",
)
parser.add_argument(
    "--region-regexp",
    metavar="re",
    type=str,
    default="^(?P<chrom>[a-zA-Z0-9]+):(?P<start>[0-9]*)-(?P<end>[0-9]*)$",
    help='''Regular expression for region matching. Only for advanced users.
    Remember to use quotes. Default:
    "^(?P<chrom>[a-zA-Z0-9]+):(?P<start>[0-9]*)-(?P<end>[0-9]*)$"''',
)

parser.add_argument(
    "--non-binary",
    action="store_const",
    dest="binDB",
    const=False,
    default=True,
    help="When the database was not binarized.",
)
parser.add_argument(
    "--reuse",
    action="store_const",
    dest="reuse",
    const=True,
    default=False,
    help="Avoid overwriting previous database walk, load instead.",
)
parser.add_argument(
    "--single",
    action="store_const",
    dest="single",
    const=True,
    default=False,
    help="""Useful flag for
    designing single probes. Same as "-X 1 -w 1". It overrides -X and -w, and
    deactivate -W, if specified.""",
)

version = "0.0.1"
parser.add_argument("--version", action="version", version=f"{sys.argv[0]} v{version}")

args = parser.parse_args()

assert not all([isinstance(args.X, type(None)), isinstance(args.W, type(None))])

reg = re.compile(args.region_regexp)
assert_msg = "provided region '%s' does not match regular expression '%s'." % (
    args.region,
    args.region_regexp,
)
assert re.match(reg, args.region), assert_msg

regdict = re.search(reg, args.region).groupdict()
required_fields = ("chrom", "start", "end")
assert all([k in regdict.keys() for k in required_fields])
args.chrom, args.start, args.end = [regdict[k] for k in required_fields]
if 0 == len(args.start) and 0 == len(args.end):
    args.start = 0
    args.end = 0
else:
    args.start = int(args.start)
    args.end = int(args.end)

if args.R > 1:
    args.R = int(args.R)
if args.r > 1:
    args.r = int(args.r)

if args.single:
    args.X = 1
    args.w = 1.0
    args.W = None

# FUNCTIONS ====================================================================


def fparse(oligo, opb=None, *args, **kwargs):
    oligo.add_score(opb.F, opb.Gs)


def fprocess(oGroup, window, *args, **kwargs):
    opb = copy.copy(kwargs["opb"])
    assert isinstance(opb, fp.oligo.OligoProbeBuilder)
    logger = logging.getLogger(kwargs["loggerName"])
    probe_list = opb.start(oGroup, window, kwargs["cfr_step"], logger)
    reduced_probe_list = opb.reduce_probe_list(probe_list)
    logger.info(
        f"Reduced from {len(probe_list)} to " + f"{len(reduced_probe_list)} probes."
    )
    return reduced_probe_list


def fimport(path, *args, **kwargs):
    return fp.oligo.OligoProbeBuilder.import_probes(path)


def fpost(results, opath, *args, **kwargs):
    assert isinstance(kwargs["opb"], fp.oligo.OligoProbeBuilder)
    logger = logging.getLogger(kwargs["loggerName"])
    if 0 == len(results):
        logger.critical(f"Built {len(results)} oligo probe candidates")
        logger.handlers = logger.handlers[:-1]
    else:
        logger.info(f"Built {len(results)} oligo probe candidates")
        logger.handlers = logger.handlers[:-1]
        fp.oligo.OligoProbeBuilder.export_probes(results, opath)
    with open(os.path.join(opath, "builder.config"), "w+") as CPH:
        kwargs["opb"].config.write(CPH)
    return (True, results)


# RUN ==========================================================================

assert_msg = "output path should NOT direct towards an existing directory"
assert_msg += " or file.\nUse '--reuse' to load previous results."
assert_msg += f" Provided path '{args.O}'"
assert not os.path.isfile(args.O), assert_msg + " leads to a file"
if not args.reuse:
    assert not os.path.isdir(args.O), assert_msg + " leads to a directory."
if not os.path.isdir(args.O):
    os.mkdir(args.O)

logFormatter = logging.Formatter(
    fp.logging.Loggable.defaultfmt, datefmt=fp.logging.Loggable.datefmt
)
logger = logging.getLogger("ifpd2-main")
logger.setLevel(logging.DEBUG)

consoleHandler = logging.StreamHandler()
consoleHandler.setFormatter(logFormatter)
consoleHandler.setLevel(logging.DEBUG)
logger.addHandler(consoleHandler)
logger = fp.logging.Loggable(logger)

logPath = "{0}/{1}.log".format(args.O, "ifpd2-main")
logger.addFileHandler(logPath)

logger.log.info(f"This log is saved at '{logPath}'.")

opb = fp.oligo.OligoProbeBuilder()
opb.N = args.N
opb.D = args.D
opb.Tr = args.t
opb.Ps = args.P
opb.Ph = args.H
opb.Po = args.I
opb.k = args.k
opb.F = args.F
opb.Gs = tuple(args.G)
opb.Ot = args.o

opsb = fp.oligo.OligoProbeSetBuilder(os.path.join(args.O, "probe_sets"), logger.log)
logger.log.info(opb.get_prologue())

if args.binDB:
    ow = fp.database.WalkerBinary(args.dbpath, logger.log)
else:
    ow = fp.database.Walker(args.dbpath, logger.log)
ow.C = args.chrom
ow.S = args.start
ow.E = args.end
ow.X = args.X
ow.Ws = args.W
ow.Wh = args.w
ow.Rs = args.R
ow.Rt = args.r
ow.out_path = args.O
ow.reuse = args.reuse
ow.threads = args.threads

ow.start(fparse, fimport, fprocess, fpost, N=args.N, opb=opb, cfr_step=ow.Rt)

opsb.build(ow.walk_results)
opsb.export()

logging.shutdown()

# END ==========================================================================

################################################################################
