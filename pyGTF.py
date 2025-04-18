#!/usr/bin/env python
# -*- coding:utf-8 -*-

from __future__ import print_function, division
import sys
import bz2
import logging
import textwrap
import argparse
from io import open as iopen
from gzip import open as gopen
from collections import OrderedDict
from os.path import exists, basename

try:  # python2
    from string import maketrans
except ImportError:
    maketrans = str.maketrans

__package__ = "pyGTF"
__version__ = "1.0.0"
__licence__ = "MIT"
__author__ = "cheng cz."
__url__ = "https://github.com/chengcz/pyGTF"
__description__ = "python parser of Fastx, GTF, NCBI GFF files"

__all__ = [
    "Files",
    "Sequence",
    "FastaReader",
    "FastqReader",
    "Interval",
    "intervals",
    "Transcript",
    "GTFReader",
    "RefSeqGFFReader",
    "BedReader",
    "genePredReader",
]


class initLogging(object):
    """
    Init logging
    """

    def __init__(self):
        # self._logger_.setLevel(logging.DEBUG)
        self._logger_ = logging.getLogger("pyGTF")
        self._logger_.setLevel(logging.DEBUG)
        console = logging.StreamHandler(sys.stderr)
        console.setLevel(logging.INFO)
        formatter = logging.Formatter("%(levelname)-7s @ %(asctime)s | %(message)s", datefmt="%m/%d/%Y %H:%M:%S")
        console.setFormatter(formatter)
        self._logger_.addHandler(console)

    def addFileHandler(self, logFile):
        handler = logging.FileHandler(logFile, mode="w", encoding="utf-8")
        handler.setLevel(logging.DEBUG)
        formatter = logging.Formatter("%(levelname)-7s @ %(asctime)s | %(message)s")
        handler.setFormatter(formatter)
        self._logger_.addHandler(handler)

    def debug(self, msg, *args, **kwargs):
        self._logger_.debug(msg, *args, **kwargs)

    def info(self, msg, *args, **kwargs):
        self._logger_.info(msg, *args, **kwargs)

    def warning(self, msg, *args, **kwargs):
        self._logger_.warning(msg, *args, **kwargs)

    def error(self, msg, *args, **kwargs):
        self._logger_.error(msg, *args, **kwargs)

    def critical(self, msg, *args, **kwargs):
        self._logger_.critical(msg, *args, **kwargs)

    def log(self, level, msg, *args, **kwargs):
        self._logger_.log(level, msg, *args, **kwargs)

    def exception(self, msg, *args, **kwargs):
        self._logger_.exception(msg, *args, **kwargs)


logger = initLogging()


class Files(object):
    """
    open file as file handle for iterator
    """

    __slots__ = ("_filepath_", "_filehandle_")

    def __init__(self, File):
        self._filepath_ = File

    def __enter__(self):
        return self

    def __exit__(self, rtype, value, trace):
        self._filehandle_.close()
        logger.debug(rtype)
        logger.debug(value)
        logger.debug(trace)

    def __iter__(self):
        if self._filepath_.lower().endswith((".gz", ".gzip")):
            self._filehandle_ = gopen(self._filepath_, "rt")
        elif self._filepath_.lower().endswith(".bz2"):
            self._filehandle_ = bz2.BZ2File(self._filepath_, "rt")
        else:
            self._filehandle_ = iopen(self._filepath_, "r")

        logger.info("reading file: {}".format(basename(self._filepath_)))
        for index, line in enumerate(self._filehandle_):
            if index % 50000 == 0 and index != 0:
                logger.debug("  working on line NO. : {:,}".format(index))
            yield line
        # fp.close()


class SequenceError(Exception):
    def __init__(self, info):
        self.message = info


class Sequence(object):
    """
    parameter
    ---------
    name:    string,
    seq:     string, ATGCNatgcn for nucl
    descri:  string,
    qualstr: string
    """

    __slots__ = ("_name", "_seq", "_descr", "_qual")
    __transtab = maketrans("ATGCNatgcnRYMKSWBVDHrymkswbvdh", "TACGNtacgnYRKMSWVBHDyrkmswvbhd")
    # base     normal     RcBase
    # R        A/G        Y
    # Y        C/T        R
    # M        A/C        K
    # K        G/T        M
    # S        G/C        S
    # W        A/T        W
    # B        G/T/C      V
    # V        G/A/C      B
    # D        G/A/T      H
    # H        A/C/T      D

    def __init__(self, name, seq, descr=None, qual=None):
        self._name = name
        self._seq = seq.strip()
        self._descr = descr if descr else ""
        if qual and len(qual) != len(seq):
            raise SequenceError("string length of sequence and qualstr is inconsistent")
        self._qual = qual

    def __str__(self):
        return "<Sequence: {}, {} ... >".format(self._name, self._seq[:50])

    def __repr__(self):
        return self.__str__()

    def __len__(self):
        return len(self._seq)

    @property
    def name(self):
        return self._name

    @property
    def seq(self):
        return self._seq

    @property
    def descr(self):
        return self._descr

    @property
    def qual(self):
        return self._qual

    @property
    def sizes(self):
        return sum(map(len, [self._name, self._seq, self._descr, self._qual]))

    def is_nucl(self):
        if set(self._seq) - set("ATGCNatgcn"):
            return False
        return True

    def reverse_complement(self):
        if not self.is_nucl():
            return None
        rcseq = self._seq.translate(self.__transtab)[::-1]
        rcqual = self._qual[::-1] if self._qual else None
        return Sequence(self._name, rcseq, self._descr, rcqual)

    def __getitem__(self, item):    # Obj[idx], slice
        assert isinstance(item, (slice, int))
        return Sequence(
            self._name, self._seq[item], self._descr,
            self._qual[item] if self._qual else None
        )

    def write_to_fastq_file(self, fp):
        if self._descr:
            fp.write("@{} {}\n".format(self._name, self._descr))
        else:
            fp.write("@{}\n".format(self._name))
        fp.write("{}\n".format(self._seq))
        fp.write("+\n")
        qual = self._qual if self._qual else "I" * len(self._seq)
        fp.write("{}\n".format(qual))

    def write_to_fasta_file(self, fp):
        if self._descr:
            fp.write(">{} {}\n".format(self._name, self._descr))
        else:
            fp.write(">{}\n".format(self._name))
        fp.write("{}\n".format(self.__formatter__(self._seq)))

    @staticmethod
    def __formatter__(seq, length=80):
        fseq = ""
        for i in range(0, len(seq), length):
            fseq += seq[i : (i + length)] + "\n"
        return fseq.strip()

    def write_to_tab_file(self, fp):
        fp.write("{}\t{}\n".format(self._id, self._seq))


class FastaReader(Files):
    """
    File parser for fasta file of nucleic acid
    """

    __slots__ = ()

    def __init__(self, fasta):
        Files.__init__(self, fasta)

    def __str__(self):
        return '<FastaReader object: "{}">'.format(self._filepath_)

    def __repr__(self):
        return self.__str__()

    def __enter__(self):
        return self

    def __exit__(self, rtype, value, trace):
        logger.debug(rtype)
        logger.debug(value)
        logger.debug(trace)

    def __iter__(self):
        seq = None
        for line in Files.__iter__(self):
            if line.startswith(">"):
                if seq:
                    yield Sequence(seqid, seq, descr)
                seqid, _, descr = line.strip(">\t \n\r").partition(" ")
                seq = ""
            else:
                if seq is None:
                    raise FileFormatError("FASTA file does not start with '>'.")
                seq += line.strip()
        yield Sequence(seqid, seq, descr)


class FileFormatError(Exception):
    def __init__(self, info):
        self.message = info


class FastqReader(object):
    """
    File parser for fastq file of nucleic acid

    parameter
    ---------
    fq1, fq2:    string, file full path
    """

    __slots__ = ("_fastq", "_is_paired")

    def __init__(self, fq1, fq2=None):
        self._fastq = (fq1, fq2) if fq2 else (fq1,)
        self._is_paired = fq2 is not None

    def __str__(self):
        return "<FastqReader object, ... >"

    def __repr__(self):
        return self.__str__()

    def __enter__(self):
        return self

    def __exit__(self, rtype, value, trace):
        logger.debug(rtype)
        logger.debug(value)
        logger.debug(trace)

    def __iter__(self):
        return self.__pe_iter__() if self._is_paired else self.__se_iter__()

    def __se_iter__(self):
        iter1 = Files(self._fastq[0]).__iter__()
        while True:
            try:
                name1, _, descr1 = next(iter1)[1:-1].partition(" ")
            except StopIteration:
                break  # end of file
            try:
                seq1, _, qual1 = (next(iter1)[:-1], next(iter1), next(iter1)[:-1])
            except StopIteration:
                raise FileFormatError("FASTQ file is incompelete.")
            else:
                yield Sequence(name1, seq1, descr1, qual1)

    def __pe_iter__(self):
        iter1 = Files(self._fastq[0]).__iter__()
        iter2 = Files(self._fastq[1]).__iter__()
        while True:
            try:
                name1, _, descr1 = next(iter1)[1:-1].partition(" ")
                name2, _, descr2 = next(iter2)[1:-1].partition(" ")
            except StopIteration:
                break  # end of file
            try:
                seq1, _, qual1 = (next(iter1)[:-1], next(iter1), next(iter1)[:-1])
                seq2, _, qual2 = (next(iter2)[:-1], next(iter2), next(iter2)[:-1])
            except StopIteration:
                raise FileFormatError("FASTQ file is incompelete.")
            if name1 != name2:
                raise FileFormatError("Sequence order of paired fastq file is not match.")
            yield (Sequence(name1, seq1, descr1, qual1),
                   Sequence(name2, seq2, descr2, qual2))


class AtomicInterval(object):
    """ """

    __slots__ = ("_lower", "_upper")

    def __init__(self, start, end):
        self._lower = start
        self._upper = end
        if self.is_empty():
            self._lower = float("inf")
            self._upper = float("-inf")

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return f'AtomicInterval: [{self._lower}, {self._upper})'

    @property
    def lower(self):
        return self._lower

    @property
    def upper(self):
        return self._upper

    def is_empty(self):
        return self._lower >= self._upper

    def overlaps(self, other, adjacent=False):
        if not isinstance(other, AtomicInterval):
            raise TypeError("Only AtomicInterval instances are supported.")
        if other.is_empty() or self.is_empty():
            return False
        first, second = (self, other) if (self.lower <= other.lower) else (other, self)
        if adjacent:
            return first.upper >= second.lower
        return first.upper > second.lower

    def __len__(self):
        return self.upper - self.lower

    def __eq__(self, other):
        if isinstance(other, AtomicInterval):
            return (self.lower == other.lower) and (self.upper == other.upper)
        else:
            return NotImplemented

    def __ne__(self, other):
        return not self.__eq__(other)

    def __lt__(self, other):
        if isinstance(other, (AtomicInterval, Interval)):
            return self.upper <= other.lower
        elif isinstance(other, int):
            return self.upper <= other
        else:
            return NotImplemented

    def __le__(self, other):
        if isinstance(other, (AtomicInterval, Interval)):
            return self.upper <= other.upper
        elif isinstance(other, int):
            return self.upper <= other
        else:
            return NotImplemented

    def __gt__(self, other):
        if isinstance(other, (AtomicInterval, Interval)):
            return self.lower >= other.upper
        elif isinstance(other, int):
            return self.lower >= other
        else:
            return NotImplemented

    def __ge__(self, other):
        if isinstance(other, (AtomicInterval, Interval)):
            return self.lower >= other.lower
        elif isinstance(other, int):
            return self.lower >= other
        else:
            return NotImplemented

    def __and__(self, other):
        if isinstance(other, AtomicInterval):
            if self.overlaps(other):
                lower = max(self.lower, other.lower)
                upper = min(self.upper, other.upper)
                return AtomicInterval(lower, upper)
            else:
                return AtomicInterval(float("inf"), float("-inf"))
        else:
            return NotImplemented

    def __or__(self, other):
        if isinstance(other, AtomicInterval):
            if self.overlaps(other, adjacent=True):
                lower = min(self.lower, other.lower)
                upper = max(self.upper, other.upper)
                return AtomicInterval(lower, upper)
            else:
                return Interval(self, other)
        else:
            return NotImplemented

    def __contains__(self, item):
        if isinstance(item, AtomicInterval):
            left = self.lower <= item.lower
            right = item.upper <= self.upper
            return left and right
        elif isinstance(item, Interval):
            for interval in item:
                if interval not in self:
                    return False
            return True
        elif isinstance(item, int):
            return self.lower <= item < self.upper
        else:
            return NotImplemented

    def __invert__(self):
        if self.is_empty():
            return AtomicInterval(float("-inf"), float("inf"))
        return Interval(
            AtomicInterval(float("-inf"), self.lower),
            AtomicInterval(self.upper, float("inf")),
        )

    def __sub__(self, other):
        if isinstance(other, AtomicInterval):
            return self & ~other
        else:
            return NotImplemented


class Interval(object):
    """ """

    __slots__ = "_intervals"

    def __init__(self, *intervals):
        self._intervals = list()

        for x in intervals:
            if isinstance(x, Interval):
                self._intervals.extend(x)
            elif isinstance(x, AtomicInterval):
                # if not x.is_empty():
                self._intervals.append(x)
            else:
                raise TypeError("Parameters must be Interval or AtomicInterval instances")
        self.__deldup__()

    def __deldup__(self):
        self._intervals = [x for x in self._intervals if not x.is_empty()]
        if len(self._intervals) == 0:
            self._intervals.append(AtomicInterval(float("inf"), float("-inf")))
        else:
            self._intervals.sort(key=lambda x: (x.lower))

            x = 0
            while x < len(self._intervals) - 1:
                first = self._intervals[x]
                second = self._intervals[x + 1]
                if first.overlaps(second, adjacent=True):
                    interval = first | second
                    self._intervals.pop(x)  # pop first
                    self._intervals.pop(x)  # pop second
                    self._intervals.insert(x, interval)
                else:
                    x = x + 1

    def __getitem__(self, item):
        return self._intervals[item]

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        regions = ["[{}, {})".format(x.lower, x.upper) for x in self]
        return "Interval: {}".format(" / ".join(regions))

    @property
    def lower(self):
        return self._intervals[0].lower

    @property
    def upper(self):
        return self._intervals[-1].upper

    def overlaps(self, other):
        if isinstance(other, Interval):
            for x in other._intervals:
                if self.overlaps(x):
                    return True
            return False
        elif isinstance(other, AtomicInterval):
            for x in self._intervals:
                if x.overlaps(other):
                    return True
            return False
        elif isinstance(other, int):
            return other in self
        else:
            return NotImplemented

    def to_atomic(self):
        return Interval(AtomicInterval(self.lower, self.upper))

    def is_atomic(self):
        return len(self._intervals) == 1

    def is_empty(self):
        return self.is_atomic() and self._intervals[0].is_empty()

    def append(self, other):
        if isinstance(other, AtomicInterval):
            self._intervals.append(other)
        elif isinstance(other, Interval):
            self._intervals.extend(other._intervals)
            # self = self | other
        else:
            raise TypeError("Parameters must be Interval/AtomicInterval instances")
        self.__deldup__()

    @property
    def length(self):
        return sum([len(x) for x in self])

    def __len__(self):
        return len(self._intervals)

    def __iter__(self):
        return iter(self._intervals)

    def __eq__(self, other):
        if isinstance(other, Interval):
            return self._intervals == other._intervals
        elif isinstance(other, AtomicInterval):
            return Interval(other) == self
        else:
            return NotImplemented

    def __ne__(self, other):
        return not self.__eq__(other)

    def __lt__(self, other):
        if isinstance(other, (AtomicInterval, Interval)):
            return self.upper <= other.lower
        elif isinstance(other, int):
            return self.upper <= other
        else:
            return NotImplemented

    def __le__(self, other):
        if isinstance(other, (AtomicInterval, Interval)):
            return self.upper <= other.upper
        elif isinstance(other, int):
            return self.upper <= other
        else:
            return NotImplemented

    def __gt__(self, other):
        if isinstance(other, (AtomicInterval, Interval)):
            return self.lower >= other.upper
        elif isinstance(other, int):
            return self.lower >= other
        else:
            return NotImplemented

    def __ge__(self, other):
        if isinstance(other, (AtomicInterval, Interval)):
            return self.lower >= other.lower
        elif isinstance(other, int):
            return self.lower >= other
        else:
            return NotImplemented

    def __and__(self, other):
        if isinstance(other, (AtomicInterval, Interval)):
            if isinstance(other, AtomicInterval):
                lst = [other]
            else:
                lst = list(other._intervals)
            new_intervals = []
            for x in self._intervals:
                for y in lst:
                    new_intervals.append(x & y)
            return Interval(*new_intervals)
        else:
            return NotImplemented

    def __rand__(self, other):
        return self & other

    def __or__(self, other):
        if isinstance(other, AtomicInterval):
            return self | Interval(other)
        elif isinstance(other, Interval):
            return Interval(*(self._intervals + other._intervals))
        else:
            return NotImplemented

    def __ror__(self, other):
        return self | other

    def __sub__(self, other):
        if isinstance(other, (AtomicInterval, Interval)):
            return self & ~other
        else:
            return NotImplemented

    def __rsub__(self, other):
        return other & ~self

    def __invert__(self):
        complements = [~x for x in self._intervals]
        intersection = complements[0]
        for interval in complements:
            intersection = intersection & interval
        return intersection

    def __contains__(self, item):
        if isinstance(item, Interval):
            for x in item._intervals:
                if x not in self:
                    return False
            return True
        elif isinstance(item, AtomicInterval):
            for x in self._intervals:
                if item in x:
                    return True
            return False
        elif isinstance(item, int):
            for x in self._intervals:
                if item in x:
                    return True
            return False
        else:
            return NotImplemented


def intervals(start, end):
    """
    parameter
    ---------
    start:    int/list of int
    end:      int/list of int
    """
    if isinstance(start, list) and isinstance(end, list):
        return Interval(*[AtomicInterval(x, y) for x, y in zip(start, end)])
    elif isinstance(start, int) and isinstance(end, int):
        return Interval(AtomicInterval(start, end))
    else:
        raise Exception("wrong input")


class GenomicInterval(object):
    '''
    parameter
    ---------
    chrom:     string, [required]
    start:     int/list, [required]
    end:       int/list, [required]
    name:      string
    strand:    string, (-, +, None)
    '''
    __slots__ = ('_chrom', '_intervals', '_attr')
    def __init__(self, chrom, start, end, name=None, strand=None, **attr):
        self._chrom = chrom
        self._intervals = self._input_to_interval(start, end)
        self._name = name
        self._strand = strand
        self._attr = attr

    @staticmethod
    def _input_to_interval(start, end):
        region = []
        if isinstance(start, int) and isinstance(end, int):
            region.append(AtomicInterval(start, end))
        elif isinstance(start, list) and isinstance(end, list):
            if len(start) != len(end):
                raise Exception('length of start and end is unequal in Input interval')
            for x, y in zip(start, end):
                if not (isinstance(x, int) and isinstance(y, int)):
                    raise TypeError('Input interval data type is not int')
                region.append(AtomicInterval(x, y))
        else:
            raise TypeError('Unsupported input data type!')
        return Interval(region)

    def __str__(self):
        if self.strand:
            return f'GenomicInterval: {self.chrom}:{self.start}-{self.end}:{self.strand}'
        else:
            return f'GenomicInterval: {self.chrom}:{self.start}-{self.end}'

    def __repr__(self):
        return self.__str__()

    @property
    def chrom(self):
        return self._chrom

    @property
    def start(self):
        return self._intervals.lower

    @property
    def end(self):
        return self._intervals.upper

    @property
    def strand(self):
        return self._strand

    @property
    def name(self):
        return self._name

    @property
    def interval(self):
        return self._intervals

    def append_interval(self, start, end):
        region = AtomicInterval(start, end)
        self._intervals = self._intervals | region

    def mask_interval(self, start, end):
        region = AtomicInterval(start, end)
        self._intervals = self._intervals - region

    def overlaps(self, other):
        if isinstance(other, GenomicInterval):
            if self._chrom == other.chrom:
                if self._intervals.overlaps(self.region):
                    return True
            return False
        else:
            return NotImplemented

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        regions = [f'{x.lower} - {x.upper}' for x in self._intervals]
        return 'GenomicInterval: {}: {}'.format(self._chrom, " / ".join(regions))

    def __len__(self):
        return sum([len(x) for x in self._intervals])

    def __eq__(self, other):
        if isinstance(other, GenomicInterval):
            if (self._chrom == other.chrom) and (self._intervals == other.region):
                return True
            return False
        else:
            return NotImplemented

    def __ne__(self, other):
        return not self.__eq__(other)

    def __lt__(self, other):
        if isinstance(other, GenomicInterval):
            if (self._chrom == other.chrom):
                return self._intervals < other.region
            return False
        elif isinstance(other, int):
            return self._intervals.upper < other
        else:
            return NotImplemented

    def __le__(self, other):
        if isinstance(other, GenomicInterval):
            if (self._chrom == other.chrom):
                return self._intervals <= other.region
            return False
        elif isinstance(other, int):
            return self._intervals.upper <= other
        else:
            return NotImplemented

    def __gt__(self, other):
        if isinstance(other, GenomicInterval):
            if (self._chrom == other.chrom):
                return self._intervals > other.region
            return False
        elif isinstance(other, int):
            return self._intervals.lower > other
        else:
            return NotImplemented

    def __ge__(self, other):
        if isinstance(other, GenomicInterval):
            if (self._chrom == other.chrom):
                return self._intervals >= other.region
            return False
        elif isinstance(other, int):
            return self._intervals.lower >= other
        else:
            return NotImplemented

    def __and__(self, other):
        if isinstance(other, GenomicInterval):
            if self._chrom != other.chrom:
                return NotImplemented
            new_intervals = []
            for interval in self.region:
                for o_interval in other.region:
                    new_intervals.append(interval & o_interval)
            return GenomicInterval(self.chrom, Interval(*new_intervals))
        else:
            return NotImplemented

    def __rand__(self, other):
        return self & other

    def __or__(self, other):
        if isinstance(other, GenomicInterval):
            if self._chrom != other.chrom:
                return NotImplemented
            region = Interval(*(self.region._intervals + other.region._intervals))
            return GenomicInterval(self.chrom, region)
        else:
            return NotImplemented

    def __ror__(self, other):
        return self | other

    def __sub__(self, other):
        if isinstance(other, GenomicInterval):
            if self._chrom != other.chrom:
                return NotImplemented
            return self & ~other
        else:
            return NotImplemented

    def __rsub__(self, other):
        return other & ~self

    def __invert__(self):
        intervals = [~x for x in self.region._intervals]
        intersection = intervals[0]
        for interval in intervals:
            intersection = intersection & interval
        return GenomicInterval(self.chrom, Interval(*intersection))

    def __contains__(self, other):
        if isinstance(other, GenomicInterval):
            if self._chrom != other.chrom:
                return False
            for o_interval in other.region._intervals:
                if o_interval not in self.region:
                    return False
            return True
        elif isinstance(other, int):
            for interval in self.region._intervals:
                if other in interval:
                    return True
            return False
        else:
            return NotImplemented


class StructureInforError(Exception):
    def __init__(self, info):
        self.message = info


class Transcript(object):
    """
    Transcript structure annotation information

    parameter
    ---------
    Tid:        sring, Transcript IDs
    chro:       sting, chromesome ID
    interval:   Interval, 0-based location
    strand:     string, choices from {-, +}
    exon:       Interval,
    cds:        Interval,
    utr:        Interval,
    infor:      dict, attr of gtf file
    """

    __slots__ = (
        "_id",
        "_chrom",
        "_interval",
        "_start",
        "_end",
        "_strand",
        "_exon",
        "_cds",
        "_utr",
        "_attri",
        "_utr5",
        "_utr3",
        "_msg",
    )

    def __init__(
        self,
        Tid,
        chro,
        interval,
        strand,
        exon=None,
        cds=None,
        utr=None,
        infor=None,
        strict=True,
    ):
        self._id = Tid
        self._chrom = chro
        self._interval = interval
        self._start = interval.lower
        self._end = interval.upper
        self._strand = strand
        self._exon = exon if (exon is not None) else Interval()
        self._cds = cds if (cds is not None) else Interval()
        self._utr = utr if (utr is not None) else Interval()
        self._attri = infor if infor else {}
        ### adjustment and validation
        self.__Interval_self_adjustment()
        self._utr5 = Interval(*[x for x in self._utr if x < self.cds])
        self._utr3 = Interval(*[x for x in self._utr if x > self.cds])
        self._msg = "\n".join(
            [
                "-" * 50,
                "    Transcript ID: {}".format(Tid),
                "    Location:      {}:{}-{}:{}".format(chro, interval.lower, interval.upper, strand),
                "    Exon:          {}".format(self._exon),
                "    Cds:           {}".format(self._cds),
                "    Utr:           {}".format(self._utr),
                "    Attri:         {}".format(self._attri.items()),
                "-" * 50,
            ]
        )
        try:
            self.__data_validation(strict)
        except StructureInforError as e:
            if strict:
                logger.error("\n{}".format(self._msg))
                raise StructureInforError(e)
            logger.warning("\n{}".format(self._msg))
        if self._strand == "-":
            self._utr5, self._utr3 = self._utr3, self._utr5

    def __str__(self):
        return "<Transcript object>\n{}".format(self._msg)

    def __repr__(self):
        return self.__str__()

    def __Interval_self_adjustment(self):
        flag_exon = not self._exon.is_empty()
        flag_cds = not self._cds.is_empty()
        flag_utr = not self._utr.is_empty()
        if flag_exon and flag_cds:
            if not flag_utr:
                self._utr = self._exon - self._cds
        elif flag_exon and not flag_cds:
            if flag_utr:
                self._cds = self._exon - self._utr
            else:
                self._cds = Interval()
                self._utr = Interval()
        elif not flag_exon and flag_cds:
            if flag_utr:
                self._exon = self._cds | self._utr
            else:
                self._exon = self._cds
                self._utr = Interval()
        else:
            self._exon = intervals(self._start, self._end)
            logger.warning("transcript interval as Exon interval, {}".format(self._id))
            # raise StructureInforError("Exon & Cds interval both is missing.")

    def __data_validation(self, strict=True):
        if self._strand not in ("-", "+"):
            raise StructureInforError("unsupported strand Symbol, select from (-, +)")
        if self._interval.is_empty():
            raise StructureInforError("Transcript interval is wrong.")
        if self._exon.is_empty():
            raise StructureInforError("Exon interval is wrong.")

        # if self._exon.to_atomic() != self._interval:
        #     raise StructureInforError("Unmatched interval between Transcript with Exon.")
        if (not self._cds.is_empty()) or (not self._utr.is_empty()):
            if self._exon != (self._cds | self._utr):
                raise StructureInforError("Unmatched interval between Exon with CDS & UTR.")
        if (not self._cds.is_empty()) and (not self._utr.is_empty()):
            if not (self._cds & self._utr).is_empty():
                raise StructureInforError("CDS interval Overlap with UTR interval.")

    def __len__(self):
        return self.length

    @property
    def id(self):
        return self._id

    @property
    def chrom(self):
        return self._chrom

    @property
    def start(self):
        """
        1-based
        """
        return self._start + 1

    @property
    def end(self):
        return self._end

    @property
    def strand(self):
        return self._strand

    @property
    def exon(self):
        return self._exon

    @property
    def exon_count(self):
        return len(self._exon)

    @property
    def length(self):
        return self._exon.length

    @property
    def intron(self):
        TransRegion = AtomicInterval(self.exon.lower, self.exon.upper)
        # intron = self._interval - self._exon
        # intron = Interval(
        #     *[x for x in intron if x.overlaps(TransRegion, False)]
        # )
        intron = TransRegion - self.exon
        return intron

    @property
    def cds(self):
        return self._cds

    @property
    def cds_count(self):
        return len(self._cds)

    @property
    def length_cds(self):
        return self._cds.length

    @property
    def utr(self):
        return self._utr

    @property
    def utr5(self):
        return self._utr5

    @property
    def utr3(self):
        return self._utr3

    @property
    def name(self):
        return self._attri.get("transcript_name", self._id)

    @property
    def biotype(self):
        tmp = self._attri
        return tmp.get(
            "transcript_type", tmp.get("transcript_biotype", tmp.get("biotype", ""))
        )

    @property
    def gene_id(self):
        return self._attri.get("gene_id", self._id)

    @property
    def gene_name(self):
        return self._attri.get("gene_name", "")

    @property
    def gene_biotype(self):
        tmp = self._attri
        return tmp.get("gene_type", tmp.get("gene_biotype", tmp.get("biotype", "")))

    @property
    def bed(self):
        """
        basic bed:
            chro, start, end, name, score, strand
        """
        return (self._chrom, self._start, self._end, self._id, 0, self._strand)

    def del_version(self, sep="."):
        """
        remove version infor of transcript id, defaule version infor at '.' after
        """
        tid = self._id
        if sep in tid:
            self._id = tid[: tid.rindex(sep)]
        gid = self.gene_id
        if sep in gid:
            self._attri["gene_id"] = gid[: gid.rindex(sep)]

    def id_modifer(self, func):
        """
        modify transcript_id, gene_id

        parameter
        ---------
        func:     function, recommend use lambda
        """
        self._id = func(self._id)
        self._attri["gene_id"] = func(self.gene_id)

    def to_gtf(self, fp):
        """
        parameter
        ---------
        fp:    file handle for output standrand GTF file
        """
        attr = self.__attri_of_gtfline()

        transcript = (
            self._chrom,
            ".",
            "transcript",
            self._start + 1,
            self._end,
            ".",
            self._strand,
            ".",
            attr,
        )
        Record = []
        tmp = zip(
            ("exon", "CDS", "5UTR", "3UTR"),
            (self._exon, self._cds, self._utr5, self._utr3),
        )
        for label, regions in tmp:
            for region in regions:
                if region.is_empty():
                    continue
                start, end = region.lower, region.upper
                lines = (
                    self._chrom,
                    ".",
                    label,
                    start + 1,
                    end,
                    ".",
                    self._strand,
                    ".",
                    attr,
                )
                Record.append(lines)
        order = False if self._strand == "+" else True
        Record = sorted(Record, key=lambda x: x[3], reverse=order)

        Record.insert(0, transcript)
        for feature in Record:
            fp.write(self.__list2str(feature))

    def __attri_of_gtfline(self):
        skipkeys = {
            "gene_id",
            "gene_name",
            "gene_type",
            "gene_biotype",
            "transcript_id",
            "transcript_name",
            "transcript_type",
            "transcript_biotype",
            "Parent",
            "ID",
            "gene",
            "gbkey",
            "Name",
        }
        keepkeys = set(self._attri.keys()) - skipkeys
        attr = [
            'transcript_id "{}"; gene_id "{}"; '.format(self.id, self.gene_id),
        ]
        keys = ["transcript_name", "gene_name", "transcript_type", "gene_type"]
        for key in keys:
            value = self._attri.get(key, None)
            if value:
                attr.append('{} "{}"; '.format(key, value))
        attr += ['{} "{}"; '.format(i, self._attri[i]) for i in keepkeys]
        return "".join(attr)

    @staticmethod
    def __list2str(lst):
        return "\t".join([str(i) for i in lst]) + "\n"

    def to_bed(self, fp):
        """
        parameter
        ---------
        fp:    file handle for output 12 columns bed file
        """
        if not self._cds.is_empty():
            cstart, cend = self._cds.lower, self._cds.upper
        else:
            cstart, cend = (self.end, self.end)
        exon_num = self.exon_count
        exon_len = "".join(["{},".format(len(x)) for x in self._exon])
        exon_start = "".join(["{},".format(x.lower - self._start) for x in self._exon])

        Record = (
            self._chrom,
            self._start,
            self._end,
            self._id,
            0,
            self._strand,
            cstart,
            cend,
            0,
            exon_num,
            exon_len,
            exon_start,
        )
        fp.write(self.__list2str(Record))

    def to_genePred(self, fp, refFlat=True):
        """
        parameter
        ---------
        fp:     file handle for output GenePred file,
        refFlat:   bool
            whether create refFlat style genePred, {True, False}

        reference
        ---------
        url: https://genome.ucsc.edu/FAQ/FAQformat#format9
        """
        if not self._cds.is_empty():
            cstart, cend = self._cds.lower, self._cds.upper
        else:
            cstart, cend = (self.end, self.end)
        exon_num = self.exon_count
        estart = "".join(["{},".format(i.lower) for i in self._exon])
        eend = "".join(["{},".format(i.upper) for i in self._exon])

        Record = [
            self._id,
            self._chrom,
            self._strand,
            self._start,
            self._end,
            cstart,
            cend,
            exon_num,
            estart,
            eend,
        ]
        if refFlat:
            Record.insert(0, self.gene_id)
        fp.write(self.__list2str(Record))

    @property
    def summary(self):
        """
        return
        ------
        id, chro, start, end, strand, name, type, length, gene_id, gene_name, gene_type
        """
        summary = (
            self.id,
            self.chrom,
            self.start,
            self.end,
            self.strand,
            self.name,
            self.biotype,
            self.length,
            self.gene_id,
            self.gene_name,
            self.gene_biotype,
        )
        return summary

    def overlap_with(self, other, flag_strand=True):
        """ """
        flag = False
        if isinstance(other, Transcript):
            if self.chrom != other.chro:
                return False
            if flag_strand and self.strand != other.strand:
                return False
            if (self._exon & other.exon).is_empty():
                return False
            return True
        else:
            return NotImplemented

    @staticmethod
    def __extract_seq(faidx, chro, start, end):
        try:
            seq = faidx[chro][start:end]
        except KeyError:
            seq = ""
            logger.error(
                'Chromesome("{}") does not exist in fasta file, Skip...'.format(chro)
            )
        return seq

    def extract_genomic_seq(self, seqfp):
        """
        parameter
        ---------
        seqfp: dict
            {chro: seq, }

        return Sequence object
        """
        seq = self.__extract_seq(seqfp, self._chrom, self._start, self._end)
        obj = Sequence(
            "{}".format(self._id),
            seq,
            "genomic sequence, gene_id:{}".format(self.gene_id),
        )
        if self._strand == "-":
            obj = obj.reverse_complement()
        return obj

    def extract_transcript_seq(self, seqfp):
        """
        return Sequence object
        """
        seq = "".join(
            [
                self.__extract_seq(seqfp, self._chrom, x.lower, x.upper)
                for x in self._exon
            ]
        )
        obj = Sequence(
            self._id, seq, "transcript sequence, gene_id:{}".format(self.gene_id)
        )
        if self._strand == "-":
            obj = obj.reverse_complement()
        return obj

    def extract_cds_seq(self, seqfp):
        """
        return Sequence object
        """
        if self._cds.is_empty():
            return None
        seq = "".join(
            [
                self.__extract_seq(seqfp, self._chrom, x.lower, x.upper)
                for x in self._cds
            ]
        )
        obj = Sequence(
            self._id, seq, "coding sequence, gene_id:{}".format(self.gene_id)
        )
        if self._strand == "-":
            obj = obj.reverse_complement()
        return obj

    def extract_exon_seq(self, seqfp):
        """
        return tuple of Sequence
        """
        objlst = []
        exon_num = self.exon_count
        for index, x in enumerate(self._exon):
            order = index + 1 if self._strand == "+" else exon_num - index
            obj = Sequence(
                "{}_exon{}".format(self._id, order),
                self.__extract_seq(seqfp, self._chrom, x.lower, x.upper),
                "exonic sequence, gene_id:{}".format(self.gene_id),
            )
            if self._strand == "-":
                obj = obj.reverse_complement()
            objlst.append(obj)
        return tuple(objlst)

    def extract_intron_seq(self, seqfp):
        """
        return tuple of Sequence
        """
        if self.intron.is_empty():
            return tuple()
        intron_num = self.exon_count - 1  # len(self._exon)
        # assert len(self.intron) == intron_num, (self.intron, self._exon)
        objlst = []
        for index, x in enumerate(self.intron):
            # x, y = self._exon[index], self._exon[index + 1]
            order = index + 1 if self._strand == "+" else intron_num - index
            obj = Sequence(
                "{}_intron{}".format(self._id, order),
                self.__extract_seq(seqfp, self._chrom, x.lower, x.upper),
                "intronic sequence, gene_id:{}".format(self.gene_id),
            )
            if self._strand == "-":
                obj = obj.reverse_complement()
            objlst.append(obj)
        return tuple(objlst)

    def extract_utr5_seq(self, seqfp):
        """
        return Sequence object
        """
        if self._utr5.is_empty():
            return None
        seq = "".join(
            [
                self.__extract_seq(seqfp, self._chrom, x.lower, x.upper)
                for x in self._utr5
            ]
        )
        obj = Sequence(
            "{}_utr5".format(self._id),
            seq,
            "5' untranslated region sequence, gene_id:{}".format(self.gene_id),
        )
        if self._strand == "-":
            obj = obj.reverse_complement()
        return obj

    def extract_utr3_seq(self, seqfp):
        """
        return Sequence object
        """
        if self._utr3.is_empty():
            return None
        seq = "".join(
            [
                self.__extract_seq(seqfp, self._chrom, x.lower, x.upper)
                for x in self._utr3
            ]
        )
        obj = Sequence(
            "{}_utr3".format(self._id),
            seq,
            "3' untranslated region sequence, gene_id:{}".format(self.gene_id),
        )
        if self._strand == "-":
            obj = obj.reverse_complement()
        return obj

    def extract_utr_seq(self, seqfp):
        """
        return tuple of Sequence
        """
        utr5 = self.extract_utr5_seq(seqfp)
        utr3 = self.extract_utr3_seq(seqfp)
        return tuple([x for x in (utr5, utr3) if x is not None])

    def extract_promoter(self, seqfp, flank=2000):
        """
        """
        if self._strand == "+":
            seq = self.__extract_seq(seqfp, self._chrom, self._start - flank, self._start)
            seq = Sequence("{} {}:{}-{} promoter".format(self._id, self._chrom, self._start - flank, self._start), seq)
        elif self._strand == "-":
            seq = self.__extract_seq(seqfp, self._chrom, self._end, self._end + flank)
            seq = Sequence("{} {}:{}-{} promoter".format(self._id, self._chrom, self._end, self._end + flank), seq)
            seq = seq.reverse_complement()
        return seq

class GTFReader(Files):
    """
    File parser for GTF/GFF file of gene annotation

    parameter
    ---------
    gtf:            string, file of gtf format
    flag_stream:    bool, parse file style, stream style use less memery
    """

    __slots__ = (
        "_gene_feature",
        "_transcript_feature",
        "_utr_feature",
        "_keep_feature",
        "_skip_feature",
        "_novel_feature",
        "_flag_stream",
        "_drop_attr",
        "_strict",
    )

    def __init__(self, gtf, flag_stream=True, strict=True):
        self._gene_feature = {
            "gene",
            "ncRNA_gene",
            "pseudogene",
        }
        self._transcript_feature = {
            "mRNA",
            "transcript",
            "tRNA",
            "rRNA",
            "snRNA",
            "snoRNA",
            "pre_miRNA",
            "pseudogenic_transcript",
            "lnc_RNA",
            "SRP_RNA",
            "RNase_MRP_RNA",
        }
        self._utr_feature = {
            "five_prime_UTR",
            "five_prime_utr",
            "5UTR",
            "three_prime_UTR",
            "three_prime_utr",
            "3UTR",
            "UTR",
            "UTR_5",
            "UTR_3",
        }
        self._keep_feature = {
            "mRNA",
            "transcript",
            "exon",
            "CDS",
            "tRNA",
            "rRNA",
            "snRNA",
            "snoRNA",
            "pre_miRNA",
            "pseudogenic_transcript",
            "lnc_RNA",
            "SRP_RNA",
            "RNase_MRP_RNA",
            "five_prime_UTR",
            "five_prime_utr",
            "5UTR",
            "three_prime_UTR",
            "three_prime_utr",
            "3UTR",
            "UTR",
            "UTR_5",
            "UTR_3",
        }
        self._skip_feature = {
            "chromosome",
            "biological_region",
            "Selenocysteine",
            "start_codon",
            "stop_codon",
            "gene",
            "ncRNA_gene",
            "pseudogene",  # skip gene feature
        }
        self._novel_feature = []
        self._flag_stream = flag_stream
        self._strict = strict
        self._drop_attr = ("ID", "Parent", "Name")
        Files.__init__(self, gtf)

    def __str__(self):
        return '<GTFReader object: "{}">'.format(self._filepath_)

    def __repr__(self):
        return self.__str__()

    def __enter__(self):
        return self

    def __exit__(self, rtype, value, trace):
        logger.debug(rtype)
        logger.debug(value)
        logger.debug(trace)

    def __iter__(self):
        Iterator = self.__stream_mode() if self._flag_stream else self.__store_mode()
        for trans in Iterator:
            _transid, _chro, _interval, _strand, _exon, _cds, _utr, _attr = trans
            _attr = {x: y for x, y in _attr.items() if (x not in self._drop_attr) and y}
            yield Transcript(
                _transid,
                _chro,
                _interval,
                _strand,
                exon=_exon,
                cds=_cds,
                utr=_utr,
                infor=_attr,
                strict=self._strict,
            )

    def __stream_mode(self):
        """ """
        t_id, t_exon, t_cds, t_utr, t_info = (
            None,
            Interval(),
            Interval(),
            Interval(),
            {},
        )
        for line in self.__read_gtf():
            chro, _, feature, start, end, _, strand, _, attr, line_id = line

            if t_id and ((t_id != line_id) or ((t_id == line_id) and (t_chro != chro))):
                yield (t_id, t_chro, t_interval, t_strand, t_exon, t_cds, t_utr, t_info)
                t_exon, t_cds, t_utr, t_info = Interval(), Interval(), Interval(), {}

            if feature in self._transcript_feature:
                t_id, t_chro, t_interval, t_strand = (
                    line_id,
                    chro,
                    intervals(start, end),
                    strand,
                )
                try:  # gtf format
                    t_info["gene_id"] = attr.pop("gene_id")
                    t_info["gene_name"] = attr.get("gene_name", "")
                except KeyError:
                    t_info["gene_id"] = (
                        attr.pop("Parent") if ("Parent" in attr) else attr["ID"]
                    )
                    t_info["transcript_name"] = attr.get("Name", "")
                t_info.update(attr)
            elif feature == "exon":
                t_exon = t_exon | intervals(start, end)
            elif feature == "CDS":
                t_cds = t_cds | intervals(start, end)
            elif feature in self._utr_feature:
                t_utr = t_utr | intervals(start, end)
        yield (t_id, t_chro, t_interval, t_strand, t_exon, t_cds, t_utr, t_info)

    def __store_mode(self):
        """
        parse structure data of GTF/GFF file as python object to memory,
        store mode methods is Applicable to unsorted files,
        the same transcript annotation line is no longer in the same block
        """
        logger.info("Start Read gff file to python object...")

        annot = OrderedDict()
        for line in self.__read_gtf():
            chro, _, feature, start, end, _, strand, _, attr, t_id = line

            if feature in self._transcript_feature:
                try:
                    geneid = attr.pop("gene_id")
                except KeyError:
                    geneid = attr.pop("Parent")
                t_info = {"gene_id": geneid}
                t_info.update(attr)
                annot.setdefault((chro, strand, t_id), {})["summary"] = [
                    intervals(start, end),
                    t_info,
                ]
            elif feature == "exon":
                annot.setdefault((chro, strand, t_id), {}).setdefault(
                    "exon", Interval()
                ).append(intervals(start, end))
            elif feature == "CDS":
                annot.setdefault((chro, strand, t_id), {}).setdefault(
                    "cds", Interval()
                ).append(intervals(start, end))
            elif feature in self._utr_feature:
                annot.setdefault((chro, strand, t_id), {}).setdefault(
                    "utr", Interval()
                ).append(intervals(start, end))

        logger.info("Done of parse gff, sort transcript list...")
        for uniqx in annot:
            chro, strand, transid = uniqx
            t_interval, t_info = annot[uniqx]["summary"]
            t_exon = annot[uniqx].get("exon", Interval())
            t_cds = annot[uniqx].get("cds", Interval())
            t_utr = annot[uniqx].get("utr", Interval())
            yield (transid, chro, t_interval, strand, t_exon, t_cds, t_utr, t_info)

    def __read_gtf(self):
        logger.info(
            "Skip Known Annotation Feature: \n    ({})".format(", ".join(self._skip_feature))
        )
        for line in Files.__iter__(self):
            if line.startswith("#") or not line.strip():
                continue

            (
                chro,
                source,
                feature,
                start,
                end,
                score,
                strand,
                frame,
                attr,
            ) = line[:-1].split("\t")
            if feature not in self._keep_feature:
                if (feature not in self._skip_feature) and (
                    feature not in self._novel_feature
                ):
                    logger.warning("skip novel annotation feature: {}".format(feature))
                    self._novel_feature.append(feature)
                continue
            start, end = int(start) - 1, int(end)
            if "=" in attr:
                attr = [x.strip().partition("=") for x in attr.split(";") if x.strip()]
                attr = {x[0]: x[2].strip("'\t\n\r\"") for x in attr}
                line_id = (
                    attr["ID"]
                    if feature in self._transcript_feature
                    else attr["Parent"]
                )
            else:
                attr = [x.strip().partition(" ") for x in attr.split(";") if x.strip()]
                attr = {x[0]: x[2].strip("'\t\n\r\"") for x in attr}
                line_id = attr["transcript_id"]
            yield chro, source, feature, start, end, score, strand, frame, attr, line_id


class RefSeqGFFReader(Files):
    """
    File parser for GFF file of gene annotation from NCBI RefSeq or Genome database

    parameter
    ---------
    gtf:     string
    chrom:   string, the table of chrom symbol convert
    strict:  bool, whether to strictly verify the gene/exon interval
    """

    __slots__ = ("_strict", "_chrom")

    def __init__(self, gtf, chrom=None, strict=False):
        Files.__init__(self, gtf)
        self._strict = strict
        self._chrom = self.__convert_chrom_id(chrom) if chrom else {}

    def __str__(self):
        return '<RefSeqGFFReader object: "{}">'.format(self._filepath_)

    def __repr__(self):
        return self.__str__()

    def __enter__(self):
        return self

    def __exit__(self, rtype, value, trace):
        logger.debug(rtype)
        logger.debug(value)
        logger.debug(trace)

    def __convert_chrom_id(self, tab):
        with iopen(tab) as f:
            tab = [i.strip().split()[:2] for i in f if not i.startswith("#")]
        return {i[0]: i[1] for i in tab}

    def __iter__(self):
        for uniqx in self.__parse_refseq_gff():
            transid_, chro, interval_, strand, exon, cds, attr = uniqx
            yield Transcript(
                transid_,
                self._chrom.get(chro, chro),
                interval_,
                strand,
                exon=exon,
                cds=cds,
                infor=attr,
                strict=self._strict,
            )

    def __parse_refseq_gff(self):
        RegulateRegion = {
            "DNAseI_hypersensitive_site",
            "enhancer",
            "enhancer_blocking_element",
            "insulator",
            "promoter",
            "protein_binding_site",
            "replication_regulatory_region",
            "transcriptional_cis_regulatory_region",
        }
        logger.info(
            "skip annotation feature as Regulating region: \n    ({})".format(", ".join(RegulateRegion))
        )

        MotifRegion = {
            "centromere",
            "direct_repeat",
            "microsatellite",
            "tandem_repeat",
            "mobile_genetic_element",
            "nucleotide_motif",
            "repeat_instability_region",
        }
        logger.info(
            "skip annotation feature as Motif sequence region: \n    ({})".format(", ".join(MotifRegion))
        )

        UnknownRegion = {
            "cDNA_match",
            "repeat_region",
            "D_loop",
            "match",
            "region",
            "origin_of_replication",
            "sequence_feature",
            "biological_region",  # Igl
            "sequence_alteration",
            "CAGE_cluster",
            "meiotic_recombination_region",
            "mitotic_recombination_region",
        }
        logger.info(
            "skip annotation feature of unknown: \n    ({})".format(", ".join(UnknownRegion))
        )
        SkipRegion = UnknownRegion | MotifRegion | RegulateRegion

        GeneFeature = {"gene", "pseudogene"}
        TransFeature = {
            "mRNA",
            "lnc_RNA",
            "ncRNA",
            "antisense_RNA",
            "transcript",  # non coding
            "RNase_MRP_RNA",
            "RNase_P_RNA",
            "Y_RNA",
            "tRNA",
            "rRNA",
            "snoRNA",
            "snRNA",
            "miRNA",
            "primary_transcript",  # miRNA precursor seq
            "C_gene_segment",
            "D_gene_segment",
            "V_gene_segment",
            "J_gene_segment",
            "SRP_RNA",
            "telomerase_RNA",
            "vault_RNA",
            "guide_RNA",
            "scRNA",
        }
        logger.info(
            (
                "the following annotation feature as transcript: \n"
                "    ({})".format(", ".join(TransFeature))
            )
        )
        KeepFeature = GeneFeature | TransFeature | {"CDS", "exon"}
        NovelFeature = []

        logger.info("Start Read gff file...")

        GeneAnnot, TransAnnot = OrderedDict(), dict()
        # Gidlst, GAnnot, TAnnot, Tstructure = [], {}, {}, {}
        for line in Files.__iter__(self):
            if line.startswith("#") or not line.strip():
                continue
            chro, _, feature, start, end, _, strand, _, attr = line[:-1].split("\t")
            if feature not in KeepFeature:
                if (feature not in SkipRegion) and (feature not in NovelFeature):
                    logger.warning("skip novel annotation feature: {}".format(feature))
                    NovelFeature.append(feature)
                continue
            start, end = int(start) - 1, int(end)
            attr = [x.strip().partition("=") for x in attr.split(";") if x.strip()]
            attr = {x[0]: x[2].strip("'\t\n\r\"") for x in attr}
            # Dbxref = [i.partition(':') for i in attr.pop('Dbxref', '').split(',')]
            # Dbxref = {i[0]: i[2].strip('\'\t\n\r"') for i in Dbxref}
            # attr.update(Dbxref)

            biotype = attr.pop("gbkey", None)
            if feature in GeneFeature:
                gline_id = attr.pop("ID")
                attr["gene_name"] = attr.pop("gene")
                attr["gene_type"] = attr.pop("gene_biotype")
                GeneAnnot.setdefault((chro, strand, gline_id), {})["summary"] = [
                    intervals(start, end),
                    attr,
                ]
            elif feature in TransFeature:
                g_id = attr.pop("Parent")
                t_id = attr.pop("ID")
                if feature in ("miRNA", "tRNA"):
                    attr["transcript_id"] = attr.pop("product")
                    attr["transcript_name"] = attr.pop("gene")
                    attr["transcript_type"] = feature
                else:
                    attr["transcript_id"] = attr.pop("transcript_id", t_id)
                    attr["transcript_name"] = attr.pop("Name", t_id)
                    attr["transcript_type"] = (
                        "protein_coding" if biotype == "mRNA" else biotype
                    )
                if not GeneAnnot.get((chro, strand, g_id)):
                    GeneAnnot.setdefault((chro, strand, g_id), {})["summary"] = [
                        intervals(start, end),
                        attr,
                    ]
                GeneAnnot[(chro, strand, g_id)].setdefault("translst", []).append(t_id)
                TransAnnot.setdefault(t_id, {})["summary"] = [
                    intervals(start, end),
                    attr,
                ]
            elif feature == "exon":
                t_id = attr.pop("Parent")
                TransAnnot.setdefault(t_id, {}).setdefault("exon", Interval()).append(
                    intervals(start, end)
                )
            elif feature == "CDS":
                t_id = attr.pop("Parent")
                TransAnnot.setdefault(t_id, {}).setdefault("cds", Interval()).append(
                    intervals(start, end)
                )
        logger.info("Done of Read gff, parse transcript structure ...")

        drop_attr = (
            "ID",
            "Parent",
            "Name",
            "gene",
            "gbkey",
            "start_range",
            "pseudo",
            "Note",
            "description",
            "model_evidence",
            "standard_name",
        )
        index = 0
        for uniqx in GeneAnnot:
            chro, strand, geneid = uniqx
            g_interval, g_attr = GeneAnnot[uniqx]["summary"]

            translst = GeneAnnot[uniqx].get("translst", [])
            flag_pseudo = False
            if not translst:
                translst.append(geneid)
                TransAnnot.setdefault(geneid, {})["summary"] = [g_interval, g_attr]
                TransAnnot[geneid]["exon"] = g_interval
                flag_pseudo = True

            for transid in translst:
                if index % 5000 == 0 and index != 0:
                    logger.info("Already processed transcript NO. : {}".format(index))
                index += 1
                try:
                    t_interval, t_attr = TransAnnot[transid]["summary"]
                except:
                    msg = "missing transcript/gene feature in gff, tracking id: {}".format(transid)
                    raise FileFormatError(msg)
                t_exon = TransAnnot[transid].get("exon", Interval())
                t_cds = TransAnnot[transid].get("cds", Interval())
                assert not t_interval.is_empty()
                if not flag_pseudo:
                    t_attr["gene_id"] = geneid
                    t_attr["gene_name"] = g_attr.get("gene_name", None)
                    t_attr["gene_type"] = g_attr.get("gene_type", None)

                transid_ = t_attr.get("transcript_id", transid)
                t_attr = {x: y for x, y in t_attr.items() if (x not in drop_attr) and y}
                yield (transid_, chro, t_interval, strand, t_exon, t_cds, t_attr)


class BedReader(Files):
    """ """

    __slots__ = "_strict"

    def __init__(self, bed, strict=True):
        Files.__init__(self, bed)
        self._strict = strict

    def __str__(self):
        return '<BedReader object: "{}">'.format(self._filepath_)

    def __repr__(self):
        return self.__str__()

    def __enter__(self):
        return self

    def __exit__(self, rtype, value, trace):
        logger.debug(rtype)
        logger.debug(value)
        logger.debug(trace)

    def __iter__(self):
        for line in Files.__iter__(self):
            lines = line[:-1].split("\t")
            if len(lines) != 12:
                raise FileFormatError("uncomplete Bed formart, not is 12 cols")
            chro, start, end, name, _, strand, cstart, cend, _, _, elen, estart = lines
            start, end, cstart, cend = map(int, [start, end, cstart, cend])
            elen = [int(x) for x in elen.split(",") if x]
            estart = [int(x) + start for x in estart.split(",") if x]
            if not len(estart) == len(elen):
                raise StructureInforError("unmatched count of exon at line: {}".format(line))

            eend = [estart[x] + elen[x] for x, y in enumerate(elen)]
            _Exon = intervals(estart, eend)
            if cstart == cend:
                _Cds = Interval()
                attri = {"transcript_type": "noncoding"}
            else:
                _Cds = intervals(cstart, cend) & _Exon
                attri = {"transcript_type": "protein_coding"}
            yield Transcript(
                name,
                chro,
                intervals(start, end),
                strand,
                exon=_Exon,
                cds=_Cds,
                infor=attri,
                strict=self._strict,
            )


class genePredReader(Files):
    """ """

    __slots__ = "_strict"

    def __init__(self, refFlat, strict=True):
        Files.__init__(self, refFlat)
        self._strict = strict

    def __str__(self):
        return '<genePredReader object: "{}">'.format(self._filepath_)

    def __repr__(self):
        return self.__str__()

    def __enter__(self):
        return self

    def __exit__(self, rtype, value, trace):
        logger.debug(rtype)
        logger.debug(value)
        logger.debug(trace)

    def __iter__(self):
        for line in Files.__iter__(self):
            lines = line[:-1].split("\t")
            try:
                (
                    Tname,
                    chro,
                    strand,
                    start,
                    end,
                    cstart,
                    cend,
                    ecount,
                    estart,
                    eend,
                ) = lines
                Gname = None
            except:    # refFlat
                (
                    Gname,
                    Tname,
                    chro,
                    strand,
                    start,
                    end,
                    cstart,
                    cend,
                    ecount,
                    estart,
                    eend,
                ) = lines
            start, end, cstart, cend, ecount = map(int, [start, end, cstart, cend, ecount])
            estart = [int(x) for x in estart.split(",") if x]
            eend = [int(x) for x in eend.split(",") if x]
            if not len(estart) == len(eend) == ecount:
                raise StructureInforError("unmatched count of exon at line: {}".format(line))

            _Exon = intervals(estart, eend)
            if cstart == cend:
                _Cds = Interval()
                attri = {"transcript_type": "noncoding"}
            else:
                _Cds = intervals(cstart, cend) & _Exon
                attri = {"transcript_type": "protein_coding"}
            if Gname:
                attri["gene_id"] = Gname
                attri["gene_name"] = Gname
            yield Transcript(
                Tname,
                chro,
                intervals(start, end),
                strand,
                exon=_Exon,
                cds=_Cds,
                infor=attri,
                strict=self._strict,
            )


def _parse_mapping_file(mapping_file):
    """
    """
    chromName = {}
    if mapping_file and exists(mapping_file):
        with iopen(mapping_file) as fi:
            for line in fi:
                if line.startswith("#"):
                    continue
                lines = line[:-1].split()
                if len(lines) > 1:
                    oldname, newname = lines[:2]
                else:
                    oldname = newname = lines[0]
                chromName[oldname] = newname
    return chromName


def args_parser():
    parser = argparse.ArgumentParser(
        prog="pyGTF.py",
        usage="python %(prog)s [option]  ",
        description=textwrap.dedent(
        """
        Convert format of gene annotate file, Extract sequence from reference
        ----------------------------------------------------------------------
        """
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent(
        """
        ----------------------------------------------------------------------
        Example:
        ----------------------------------------------------------------------
          python %(prog)s -i example.gff -g example.gtf -b example.bed -f example.flat.txt
          python %(prog)s -i example.gff -r example.fa -cdna example.cdna.fa
        """
        )
    )
    parser.add_argument(
        "-i",
        "--input",
        dest="input_file",
        metavar="",
        required=True,
        help="[required], input file of gtf/gff/bed/genePred format",
    )
    parser.add_argument(
        "-ft",
        "--format",
        dest="input_format",
        metavar="",
        choices=("gtf", "gff", "refseqgff", "bed", "genePred"),
        help="the format of input file, default use input file suffix"
    )
    parser.add_argument(
        "-s",
        "--subset",
        dest="subset_mapping_file",
        metavar="",
        help="text file, subset of gene list",
    )
    parser.add_argument(
        "--strict",
        dest="strict_mode",
        action="store_true",
        help="strict validation of gene intervals, default: False",
    )

    aParams = parser.add_argument_group("[Output Annot Params]")
    aParams.add_argument(
        "-g", "--gtf",
        dest="gtf", metavar="", help="output gtf file"
    )
    aParams.add_argument(
        "-b", "--bed",
        dest="bed", metavar="", help="output bed file"
    )
    aParams.add_argument(
        "-f",
        "--flat",
        dest="genePred", metavar="",
        help="output genePred file, ref: https://genome.ucsc.edu/FAQ/FAQformat#format9",
    )

    sParams = parser.add_argument_group("[Sequence Params]")
    sParams.add_argument(
        "-r",
        "--reference",
        dest="reference_fasta_file",
        metavar="",
        help="input reference fasta file",
    )
    sParams.add_argument(
        "-m",
        "--mapping",
        dest="mapping_file",
        metavar="",
        help="chrom id mapping file, rename chrom id",
    )

    osParams = parser.add_argument_group("[Sequence Output Params]")
    osParams.add_argument(
        "-cdna",
        "--cdna",
        dest="cdna",
        metavar="",
        help="output cdna sequence",
    )
    osParams.add_argument(
        "-cds",
        "--cds",
        dest="cds",
        metavar="",
        help="output cds sequence",
    )
    osParams.add_argument(
        "-utr",
        "--utr",
        dest="utr",
        metavar="",
        help="output utr sequence",
    )
    osParams.add_argument(
        "-exon",
        "--exon",
        dest="exon",
        metavar="",
        help="output exon sequence",
    )
    osParams.add_argument(
        "-intron",
        "--intron",
        dest="intron",
        metavar="",
        help="output intron sequence",
    )
    osParams.add_argument(
        "-promoter",
        "--promoter",
        dest="promoter",
        metavar="",
        help="output promoter sequence",
    )
    args = parser.parse_args()

    output_fasta = [args.cdna, args.cds, args.utr, args.exon, args.intron, args.promoter]
    if any(map(lambda x: x is not None, output_fasta)) and (args.reference_fasta_file is None):
        parser.print_help()
        logger.error("required --reference parameter")
        exit(1)

    output_annot = [args.gtf, args.bed, args.genePred]
    if all(map(lambda x: x is None, output_annot + output_fasta)):
        parser.print_help()
        logger.error("Not provide output file")
        exit(1)

    return args


def _create_op_handle(fop):
    if not fop:
        return None
    if exists(fop):
        logger.warning('file "{}" is exists, overwrited!'.format(fop))
    return iopen(fop, "w")


def _guess_file_format(fop):
    """

    return:
        gtf: GTFReader
        gff: GTFReader
        bed: BedReader
        genePred: genePredReader
        refseqgff: RefSeqGFFReader
    """
    # file extend name
    filename, _, suffix = fop.lower().rpartition(".")
    if suffix in ("gz", "bz", "bz2"):
        filename, _, suffix = filename.rpartition(".")
    # file content
    for line in Files(fop):
        if line.startswith("#") or not line.strip():
            continue
    lines = line[:-1].split("\t")

    if (len(lines) == 9) or (suffix in ("gtf", "gff", "gff3")):
        if ("refseq" in filename):
            return "refseqgff"
        return "gtf"
    elif len(lines) == 12 or (suffix == "bed"):
        return "bed"
    elif (len(lines) in (10, 11)) or (suffix in ("flat", "genepred")):
        return "genePred"
    else:
        raise Exception("unknown file format, please input format")


def util():
    args = args_parser()
    inputfile, inputfmt = args.input_file, args.input_format

    subset_mapping_file, strict_mode = args.subset_mapping_file, args.strict_mode
    # output file
    opgtf, opbed, opFlat = args.gtf, args.bed, args.genePred
    # input sequence file
    genome_fasta_file, mapping_file = args.reference_fasta_file, args.mapping_file
    # output sequence file
    cdna, cds = args.cdna, args.cds
    utr, exon, intron, promoter = args.utr, args.exon, args.intron, args.promoter

    ##########
    ### starting
    subGeneList = _parse_mapping_file(subset_mapping_file)

    if genome_fasta_file:
        logger.info("parse genome fasta file: {} ... ".format(genome_fasta_file))
        chromName = _parse_mapping_file(mapping_file)
        genomeSeqDict = {
            chromName.get(i.name, i.name): i.seq
            for i in FastaReader(genome_fasta_file)
        }
        logger.info("done, parse fasta.")

    opgtf, opbed, opFlat, cdna, cds, utr, exon, intron, promoter = map(
        _create_op_handle, [opgtf, opbed, opFlat, cdna, cds, utr, exon, intron, promoter]
    )
    if not inputfmt:
        inputfmt = _guess_file_format(inputfile)
        logger.info("guess input file format: {}".format(inputfmt))

    formatparser = {
        "gtf": GTFReader,
        "gff": GTFReader,
        "bed": BedReader,
        "genePred": genePredReader,
        "refseqgff": RefSeqGFFReader,
    }
    fmtparser = formatparser[inputfmt]
    logger.info("parse gene structure file: {} ... ".format(inputfile))
    for trans in fmtparser(inputfile, strict=strict_mode):
        if subGeneList and (trans.id not in subGeneList):
            continue

        rid = trans.id
        trans._id = subGeneList.get(rid, rid)
        trans._attri["gene_id"] = subGeneList.get(rid, rid)

        if opgtf is not None:
            trans.to_gtf(opgtf)
        if opbed is not None:
            trans.to_bed(opbed)
        if opFlat is not None:
            trans.to_genePred(opFlat, refFlat=True)

        if cdna is not None:
            seq = trans.extract_transcript_seq(genomeSeqDict)
            seq.write_to_fasta_file(cdna)
        if cds is not None:
            seq = trans.extract_cds_seq(genomeSeqDict)
            seq.write_to_fasta_file(cds)
        if utr is not None:
            for seq in trans.extract_utr_seq(genomeSeqDict):
                seq.write_to_fasta_file(utr)
        if exon is not None:
            for seq in trans.extract_exon_seq(genomeSeqDict):
                seq.write_to_fasta_file(exon)
        if intron is not None:
            for seq in trans.extract_intron_seq(genomeSeqDict):
                seq.write_to_fasta_file(intron)
        if promoter is not None:
            seq = trans.extract_promoter(genomeSeqDict)
            seq.write_to_fasta_file(promoter)

    map(
        lambda x: x.close() if (x is not None) else None,
        [opgtf, opbed, opFlat, cdna, cds, utr, exon, intron, promoter],
    )


if __name__ == "__main__":
    util()
