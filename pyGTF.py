#!/usr/bin/env python
# -*- coding:utf-8 -*-

import os
import sys
import bz2
import gzip
import logging

__package__ = 'pyGTF'
__version__ = '0.1'
__licence__ = 'MIT'
__author__  = 'cheng cz.'
__url__     = 'https://github.com/chengcz/pyGTF'
__description__ = 'pure python parser of Fastx, GTF, NCBI GFF files'

__all__ = [
    'Files', 
    'Sequence', 'SequenceWithQual', 'FastaReader', 'FastqReader',
    'Transcript', 'GTFReader', 'RefSeqGFFReader', 'BedReader', 'genePredReader'
]


logging.basicConfig(
    level=logging.INFO,
    format='%(levelname)-5s @ %(filename)s:%(lineno)s, %(asctime)s, %(message)s',
    stream=sys.stderr)
logger = logging.getLogger(__name__)


class Files(object):
    '''
        open file as file handle for iterator
    '''

    def __init__(self, File):
        self._fos = File

    def __iter__(self):
        if self._fos.lower().endswith(('.gz', '.gzip')):
            fgz = True
            fp = gzip.open(self._fos)
        elif self._fos.lower().endswith('.bz2'):
            fgz = True
            fp = bz2.BZ2File(self._fos)
        else:
            fgz = False
            fp = open(self._fos)

        for index, line in enumerate(fp):
            if index % 50000 == 0 and index != 0:
                logger.info('working on line NO. of {}: {:,}'.format(
                    os.path.basename(self._fos), index
                ))
            try:
                if fgz:
                    line = line.decode()
            except:
                pass
            yield line
        fp.close()


class Sequence(object):
    '''
        nucleic acid sequence object, contain seq, name and seq description
    '''

    def __init__(self, name, seq, descr=''):
        self._id = name
        self._seq = seq.strip()
        self._descr = descr

    def __str__(self):
        return "<Sequence object: \"{}, {} ...\" >".format(self._id, self._seq[:30])

    # def __repr__(self):
    #     return "<Sequence: {}, {} ... >".format(self._id, self._seq[:50])

    def __len__(self):
        return len(self._seq)

    @property
    def name(self):
        return self._id

    @property
    def seq(self):
        return self._seq

    @property
    def descr(self):
        return self._descr

    @property
    def length(self):
        return len(self._seq)

    @property
    def is_empty(self):
        return False if self._seq else True

    def reverse_complement(self):
        paired_rc = {
            'A': 'T',
            'T': 'A',
            'C': 'G',
            'G': 'C',
            'N': 'N',
            'a': 't',
            't': 'a',
            'c': 'g',
            'g': 'c',
            'n': 'n',
        }
        try:
            seq = ''.join([paired_rc[i] for i in self._seq[::-1]])
        except KeyError:
            raise Exception('Illegal nucleotide character, out of "ATGCNatgcn".')
        return Sequence(self._id, seq, self._descr)

    def write_to_fasta_file(self, fp, chara=80):
        if self._descr:
            fp.write('>{} {}\n'.format(self._id, self._descr))
        else:
            fp.write('>{}\n'.format(self._id))
        fp.write('{}'.format(self._seq_formater(self._seq, chara)))

    def write_to_tab_file(self, fp):
        fp.write('{}\t{}\n'.format(self._id, self._seq))

    @staticmethod
    def _seq_formater(seq, length):
        tmp = ''
        for i in range(0, len(seq), length):
            tmp += seq[i : (i + length)] + '\n'
        return tmp


class SequenceWithQual(Sequence):
    '''
        nucleic acid sequence object, contain seq, name and qual
    '''

    def __init__(self, name, seq, qual):
        if len(seq) != len(qual):
            raise Exception('length is Inconsistent seq and qual string')
        self._id = name
        self._seq = seq
        self._qual = qual

    def __str__(self):
        return "<SequenceWithQual object: \"{}, {} ...\" >".format(self._id, self._seq[:30])

    def __len__(self):
        return len(self._seq)

    @property
    def name(self):
        return self._id

    @property
    def seq(self):
        return self._seq

    @property
    def qual(self):
        return [ord(i) for i in self._qual]

    @property
    def qualstr(self):
        return self._qual

    @property
    def length(self):
        return len(self._seq)

    def write_to_fastq_file(self, fp):
        fp.write('@{}\n'.format(self._id))
        fp.write('{}\n'.format(self._seq))
        fp.write('+\n{}\n'.format(self._qual))

    def write_to_fasta_file(self, fp):
        fp.write('>{}\n'.format(self._id))
        fp.write('{}\n'.format(self._seq))


class FastaReader(Files):
    '''
        File parser for fasta file of nucleic acid
    '''

    def __init__(self, fasta):
        Files.__init__(self, fasta)

    def __str__(self):
        return "<FastaReader object: \"{}\">".format(self._fos)

    def __iter__(self):
        seq = None
        for line in Files.__iter__(self):
            if line.startswith('>'):
                if seq:
                    yield Sequence(seqid, seq, descr)
                seqid, _, descr = line.strip('>\t \n').partition(' ')
                seq = ''
            else:
                assert seq is not None, 'FASTA file does not start with \'>\'.'
                seq += line.strip()
        yield Sequence(seqid, seq, descr)


class FastqFileError(Exception):
    pass


class FastqReader(Files):
    '''
        File parser for fastq file of nucleic acid
    '''

    def __init__(self, fastq):
        self._fastq = fastq
        Files.__init__(self, fastq)

    def __str__(self):
        return "<FastqReader object: \"{}\">".format(self._fos)

    def __iter__(self):
        iter = Files.__iter__(self)
        while True:
            try:
                seqid = next(iter)
            except StopIteration:
                break    # end of file
            try:
                seq, _, qual = next(iter), next(iter), next(iter)
            except StopIteration:
                raise FastqFileError('Incompelete Fastq file: {}\n'.format(self._fastq))
            seqid, seq, qual = [i.strip('\t \n') for i in [seqid[1:], seq, qual]]
            yield SequenceWithQual(seqid, seq, qual)

    def phred_judge(self):
        iter = Files.__iter__(self)
        _, _, _, qual = next(iter), next(iter), next(iter), next(iter)
        phred = min([ord(i) for i in qual.strip()])
        return 'phred+64' if phred > 64 else 'phred+33'


class Transcript(object):
    '''
        Transcript structure annotation information

        Parameter
        ---------
        Tid:   sring
            Transcript IDs
        chro:   sting
            chromesome ID
        start, end:   int
            0-based location
        strand:   string
            choices from {-, +}
        exon:  list
            [(exonstart1, exonend1), (exonstart2, exonend2), ]

        cds:   tuple
            (cds start site, cds end site)
        infor:    dict
            transcript_name, transcript_type:   string
            gene_id, gene_name, gene_type:   string
    '''

    def __init__(self, Tid, chro, start, end, strand, exon, cds=None, infor=None):
        logger.debug(('input info: ', Tid, (start, end), exon, cds))

        self._id = Tid
        self._chro = chro
        self._start = start
        self._end = end
        self._strand = strand
        self._exon = tuple([
            (x, y) for x, y in sorted(set(exon), key=lambda x: x[0], reverse=False)
        ])
        self._codingregion = cds
        self._attributes = infor if infor else {}
        self.__structure_data_validation()

        self._utr5, self._cds, self._utr3 = self.__transcriptional_region_group()
        logger.debug(
            ('output info:', Tid, (start, end), self._exon, self._utr5, self._cds, self._utr3)
        )

    def __str__(self):
        return "<Transcript object: \"{}, {}:{}-{} ...\" >".format(
            self._id, self._chro, self._start, self._end
        )

    def __len__(self):
        return self.length

    def __structure_data_validation(self):
        assert self._strand in ['-', '+'], 'unsupported strand Symbol, {-, +}'
        for each in self._exon:
            msg = '{}: Exon region Over range of the transcriptional region.'.format(self._id)
            assert self._start <= each[0] < each[1] <= self._end, msg
        if self._codingregion:
            msg = 'Error CDS postion, support format: (start, end)'
            assert len(self._codingregion) == 2, msg
            start, end = self._codingregion
            msg = '{}: CDS region Over range of the transcriptional region.'.format(self._id)
            assert self._start <= start < end <= self._end, msg

    def __transcriptional_region_group(self):
        if self._codingregion:
            cstart, cend = self._codingregion
            if cstart == cend:
                return (), (), ()
            utr5 = [(x, y) for x, y in self._exon if x < cstart]
            if utr5:
                start, end = utr5.pop(); utr5.append((start, cstart))
            cds = [
                (x, y) for x, y in self._exon
                if (x <= cstart < y) or (x < cend <= y) or \
                (cstart <= x < cend) or (cstart <= y < cend)
            ]
            if cds:
                start, end = cds.pop(); cds.append((start, cend))
                start, end = cds.pop(0); cds.insert(0, (cstart, end))
            utr3 = [(x, y) for x, y in self._exon if y > cend]
            if utr3:
                start, end = utr3.pop(0); utr3.insert(0, (cend, end))
            if self._strand == '-':
                utr5, utr3 = utr3, utr5
            return map(tuple, [utr5, cds, utr3])
        else:
            return (), (), ()

    @property
    def id(self):
        return self._id

    @property
    def chro(self):
        return self._chro

    @property
    def start(self):
        '''
            1-based
        '''
        return self._start + 1

    @property
    def end(self):
        return self._end

    @property
    def strand(self):
        return self._strand

    @property
    def bed(self):
        '''
        basic bed:
            chro, start, end, name, score, strand
        '''
        return (self._chro, self._start, self._end, self._id, 0, self._strand)

    @property
    def exon(self):
        return self._exon

    @property
    def exon_count(self):
        return len(self._exon)

    @property
    def length(self):
        return sum([(y - x) for x, y in self._exon])

    @property
    def cds(self):
        return self._cds

    @property
    def length_cds(self):
        return sum([(y - x) for x, y in self._cds])

    @property
    def name(self):
        return self._attributes.get('transcript_name', self._id)

    @property
    def type(self):
        tmp = self._attributes
        return tmp.get('transcript_type', tmp.get('transcript_biotype', tmp.get('biotype', "")))

    @property
    def gene_id(self):
        return self._attributes.get('gene_id', self._id)

    @property
    def gene_name(self):
        return self._attributes.get('gene_name', "")

    @property
    def gene_type(self):
        tmp = self._attributes
        return tmp.get('gene_type', tmp.get('gene_biotype', tmp.get('biotype', "")))

    def del_version(self, sep='.'):
        '''
        remove version infor of transcript id, defaule version infor at '.' after
        '''
        tid = self._id
        if sep in tid:
            self._id = tid[:tid.rindex(sep)]
        gid = self.gene_id
        if sep in gid:
            self._attributes['gene_id'] = gid[:gid.rindex(sep)]

    def id_modifer(self, func):
        '''
        modify transcript_id, gene_id

        parameter
        ---------
        func:     function, recommend use lambda
        '''
        self._id = func(self._id)
        self._attributes['gene_id'] = func(self.gene_id)

    def to_gtf(self, fp):
        '''
        parameter
        ---------
        fp:    file handle for output standrand GTF file
        '''
        attr = self.__attri_of_gtfline()

        transcript = (
            self._chro,
             '.',
             'transcript',
             self._start + 1,
             self._end,
             '.',
             self._strand,
             '.',
             attr,
        )
        Record = []
        tmp = zip(
            ('exon', 'CDS', '5UTR', '3UTR'),
            (self._exon, self._cds, self._utr5, self._utr3)
        )
        for label, region in tmp:
            for start, end in region:
                Record.append(
                    (self._chro, '.', label, start + 1, end, '.', self._strand, '.', attr)
                )
        order = False if self._strand == '+' else True
        Record = sorted(Record, key=lambda x: x[3], reverse=order)
        Record.insert(0, transcript)
        for feature in Record:
            fp.write(self.__list2str(feature))

    def __attri_of_gtfline(self):
        skipkeys = {
            'gene_id', 'gene_name', 'gene_type', 'gene_biotype',
            'transcript_id', 'transcript_name', 'transcript_type',
            'transcript_biotype', 'Parent', 'ID', 'gene', 'gbkey', 'Name',
        }
        keepkeys = set(self._attributes.keys()) - skipkeys
        attr = ['transcript_id "{}"; gene_id "{}"; '.format(self.id, self.gene_id), ]
        keys = ['transcript_name', 'gene_name', 'transcript_type', 'gene_type']
        for key in keys:
            value = self._attributes.get(key, None)
            if value:
                attr.append('{} "{}"; '.format(key, value))
        attr += ['{} "{}"; '.format(i, self._attributes[i]) for i in keepkeys]
        return ''.join(attr)

    @staticmethod
    def __list2str(lst):
        return '\t'.join([str(i) for i in lst]) + '\n'

    def to_bed(self, fp):
        '''
        parameter
        ---------
        fp:    file handle for output 12 columns bed file
        '''
        cstart, cend = self._codingregion if self._cds else (self._end, self._end)
        exon_num = self.exon_count
        exon_len = ''.join(['{},'.format(x[1] - x[0]) for x in self._exon])
        exon_start = ''.join(['{},'.format(x[0] - self._start) for x in self._exon])

        Record = (
            self._chro,
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

    def to_genePred(self, fp, refFlat=False):
        '''
        parameter
        ---------
        fp:     file handle for output GenePred file,
        refFlat:   bool
            whether create refFlat style genePred, {True, False}

        reference
        ---------
        url: https://genome.ucsc.edu/FAQ/FAQformat#format9
        '''
        if self._codingregion:
            cstart, cend = self._codingregion
        else:
            cstart, cend = self._end, self._end
        exon_num = self.exon_count
        estart = ''.join(['{},'.format(i[0]) for i in self._exon])
        eend = ''.join(['{},'.format(i[1]) for i in self._exon])

        Record = [
            self.name,
            self._chro,
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
        '''
        return
        ------
        id, chro, start, end, strand, name, type, length, gene_id, gene_name, gene_type
        '''
        summary = (
            self.id,
            self.chro,
            self.start,
            self.end,
            self.strand,
            self.name,
            self.type,
            self.length,
            self.gene_id,
            self.gene_name,
            self.gene_type
        )
        return summary

    def overlap_with(self, transcript):
        '''
        '''
        flag = False
        if self.chro == transcript.chro:
            if (self.start <=  transcript.start < self.end) or \
            (self.start < transcript.end <= self.end):
                flag = True
        return flag

    def __extract_seq(self, faidx, chro, start, end):
        try:
            seq = faidx[chro][start:end]
        except KeyError:
            seq = ''
            logger.error('Chromesome("{}") does not exist in fasta file, Skip...'.format(chro))
        return seq

    def extract_genomic_seq(self, seqfp):
        '''
        parameter
        ---------
        seqfp: dict
            {chro: seq, }
        '''
        seq = self.__extract_seq(seqfp, self._chro, self._start, self._end)
        obj = Sequence(
            '{}'.format(self._id),
            seq,
            'genomic sequence, gene_id:{}'.format(self.gene_id)
        )
        if self._strand == '-':
            obj = obj.reverse_complement()
        return (obj, )

    def extract_isoform_seq(self, seqfp):
        '''
        '''
        seq = ''.join([self.__extract_seq(seqfp, self._chro, x, y) for x, y in self._exon])
        obj = Sequence(
            self._id, seq, 'transcript sequence, gene_id:{}'.format(self.gene_id)
        )
        if self._strand == '-':
            obj = obj.reverse_complement()
        return (obj, )

    def extract_cds_seq(self, seqfp):
        '''
        '''
        seq = ''.join([self.__extract_seq(seqfp, self._chro, x, y) for x, y in self._cds])
        obj = Sequence(
            self._id, seq, 'coding sequence, gene_id:{}'.format(self.gene_id)
        )
        if self._strand == '-':
            obj = obj.reverse_complement()
        return (obj, )

    def extract_exon_seq(self, seqfp):
        '''
        '''
        objlst = []
        exon_num = self.exon_count
        for index, (x, y) in enumerate(self._exon):
            order = index + 1 if self._strand == '+' else exon_num - index
            obj = Sequence(
                '{}_exon{}'.format(self._id, order),
                self.__extract_seq(seqfp, self._chro, x, y),
                'exonic sequence, gene_id:{}'.format(self.gene_id),
            )
            if self._strand == '-':
                obj = obj.reverse_complement()
            objlst.append(obj)
        return tuple(objlst)

    def extract_intron_seq(self, seqfp):
        '''
        '''
        objlst = []
        intron_num = self.exon_count - 1
        if intron_num >= 1:
            for index in range(intron_num):
                _, x = self._exon[index]
                y, _ = self._exon[index + 1]
                order = index + 1 if self._strand == '+' else intron_num - index
                obj = Sequence(
                    '{}_intron{}'.format(self._id, order),
                    self.__extract_seq(seqfp, self._chro, x, y),
                    'intronic sequence, gene_id:{}'.format(self.gene_id),
                )
                if self._strand == '-':
                    obj = obj.reverse_complement()
                objlst.append(obj)
        return tuple(objlst)

    def extract_utr5_seq(self, seqfp):
        '''
        '''
        seq = ''.join([self.__extract_seq(seqfp, self._chro, x, y) for x, y in self._utr5])
        obj = Sequence(
            '{}_utr5'.format(self._id),
            seq,
            '5\' untranslated region sequence, gene_id:{}'.format(self.gene_id)
        )
        if self._strand == '-':
            obj = obj.reverse_complement()
        return (obj, )

    def extract_utr3_seq(self, seqfp):
        '''
        '''
        seq = ''.join([self.__extract_seq(seqfp, self._chro, x, y) for x, y in self._utr3])
        obj = Sequence(
            '{}_utr3'.format(self._id),
            seq,
            '3\' untranslated region sequence, gene_id:{}'.format(self.gene_id)
        )
        if self._strand == '-':
            obj = obj.reverse_complement()
        return (obj, )

    def extract_utr_seq(self, seqfp):
        '''
        '''
        utr5, = self.extract_utr5_seq(seqfp)
        utr3, = self.extract_utr3_seq(seqfp)
        return (utr5, utr3)


class GTFReader(Files):
    '''
        File parser for GTF/GFF file of gene annotation
    '''

    def __init__(self, gtf):
        self._transcript_feature = {
            'mRNA', 'transcript',
            'tRNA', 'rRNA', 'snRNA', 'snoRNA', 'pre_miRNA',
            'pseudogenic_transcript', 'lnc_RNA', 'SRP_RNA', 'RNase_MRP_RNA'
        }
        self._keep_feature = {
            'mRNA', 'transcript',
            'exon',
            'CDS',
            'tRNA', 'rRNA', 'snRNA', 'snoRNA', 'pre_miRNA',
            'pseudogenic_transcript', 'lnc_RNA', 'SRP_RNA', 'RNase_MRP_RNA'
        }
        self._skip_feature = {
            'chromosome',
            'biological_region',
            'Selenocysteine',
            'gene', 
            'start_codon', 'stop_codon',
            'five_prime_UTR', 'five_prime_utr', '5UTR',
            'three_prime_UTR', 'three_prime_utr', '3UTR',
            'UTR',
            'ncRNA_gene', 'pseudogene',
        }
        Files.__init__(self, gtf)

    def __str__(self):
        return "<GTFReader object: \"{}\">".format(self._fos)

    def __iter__(self):
        t_id, t_chro = None, None
        t_exon, t_cds, t_info = [], [], {}
        for line in self.__read_gtf():
            chro, _, feature, start, end, _, strand, _, attr, line_id = line

            if t_id and ((t_id != line_id) or ((t_id == line_id) and (t_chro != chro))):
                t_cds = sorted(t_cds, key=lambda x: x[0], reverse=False)
                t_cds = (t_cds[0][0], t_cds[-1][-1]) if t_cds else None
                yield Transcript(t_id, t_chro, t_start, t_end, t_strand, t_exon, t_cds, t_info)
                t_exon, t_cds, t_info = [], [], {}

            if feature in self._transcript_feature:
                t_id = line_id
                t_chro, t_start, t_end, t_strand = chro, start, end, strand
                try:
                    t_info['gene_id'] = attr.pop('gene_id')
                    t_info['gene_name'] = attr.get('gene_name', '')
                except KeyError:
                    t_info['gene_id'] = attr.pop('Parent')
                    t_info['transcript_name'] = attr.get('Name', '')
                t_info.update(attr)
            elif feature == 'exon':
                t_exon.append((start, end))
            elif feature == 'CDS':
                t_cds.append((start, end))
        t_cds = sorted(t_cds, key=lambda x: x[0], reverse=False)
        t_cds = (t_cds[0][0], t_cds[-1][-1]) if t_cds else None
        yield Transcript(t_id, t_chro, t_start, t_end, t_strand, t_exon, t_cds, t_info)

    def safe_mode(self):
        '''
            parse structure data of GTF/GFF file as python object to memory,
            safe_mode methods is Applicable to unsorted files,
            the same transcript annotation line is no longer in the same block
        '''
        logger.info('Start Read gff file to python object...')

        transcriptlst, annot = [], {}
        for line in self.__read_gtf():
            chro, _, feature, start, end, _, strand, _, attr, t_id = line
            transcriptlst.append((chro, t_id))

            if feature in self._transcript_feature:
                annot.setdefault((chro, t_id), {})['pos'] = [chro, start, end, strand]
                try:
                    geneid = attr.pop('gene_id')
                except KeyError:
                    geneid = attr.pop('Parent')
                t_info = {'gene_id': geneid}
                t_info.update(attr)
                annot[(chro, t_id)]['attri'] = t_info
            elif feature == 'exon':
                annot.setdefault((chro, t_id), {}).setdefault('exon', []).append((start, end))
            elif feature == 'CDS':
                annot.setdefault((chro, t_id), {}).setdefault('cds', []).append((start, end))

        logger.info('Done of parse gff, sort transcript list...')
        transcriptlst = sorted(
            set(transcriptlst),
            key=lambda x: transcriptlst.index(x),
            reverse=False
        )
        logger.info('Done of sort transcript list.')
        for CHR, iso in transcriptlst:
            chro, start, end, strand = annot[(CHR, iso)]['pos']
            assert CHR == chro, ''
            exon = annot[(CHR, iso)].get('exon', None)
            cds = annot[(CHR, iso)].get('cds', None)
            t_info = annot[(CHR, iso)].get('attri', None)
            if not any([exon, cds]):
                exon = [(start, end)]
            elif not exon:
                exon = cds.copy()
            if cds:
                cds = sorted(cds, key=lambda x: x[0], reverse=False)
                cds = (cds[0][0], cds[-1][-1])
            yield Transcript(iso, chro, start, end, strand, exon, cds, t_info)

    def __read_gtf(self):
        logger.info((
            'Skip known annotation feature: \n    ({})'.format(', '.join(self._skip_feature))
        ))
        for line in Files.__iter__(self):
            if line.startswith('#') or not line.strip():
                continue

            chro, source, feature, start, end, score, strand, frame, attr = line.strip().split('\t')
            if feature not in self._keep_feature:
                if feature not in self._skip_feature:
                    logger.warning('skip novel annotation feature: {}'.format(feature))
                continue
            start, end = int(start) - 1, int(end)
            if '=' in attr:
                attr = [i.strip().partition('=') for i in attr.split(';') if i.strip()]
                attr = {i[0]: i[2].strip('\'\t\n\r"') for i in attr}
                line_id = attr['ID'] if feature in self._transcript_feature else attr['Parent']
            else:
                attr = [i.strip().partition(' ') for i in attr.split(';') if i.strip()]
                attr = {i[0]: i[2].strip('\'\t\n\r"') for i in attr}
                line_id = attr['transcript_id']
            yield chro, source, feature, start, end, score, strand, frame, attr, line_id


class RefSeqGFFReader(Files):
    '''
        File parser for GFF file of gene annotation from NCBI RefSeq or Genome database
    '''

    def __init__(self, gtf, chrom=None):
        Files.__init__(self, gtf)
        self._chrom = self.__convert_chrom_id(chrom) if chrom else {}

    def __str__(self):
        return "<RefSeqGFFReader object: \"{}\">".format(self._fos)

    def __convert_chrom_id(self, tab):
        with open(tab) as f:
            tab = [i.strip().split()[:2] for i in f if not i.startswith('#')]
        return {i[0]: i[1] for i in tab}

    def __iter__(self):
        genelst, annot = self.__parse_gff()
        for CHR, gene in genelst:
            for iso in annot[(CHR, gene)]:
                chro, start, end, strand = annot[(CHR, gene)][iso]['infor']
                assert CHR == chro, ''
                chro = self._chrom.get(chro, chro)
                exon =  annot[(CHR, gene)][iso].get('exon', None)
                cds  =  annot[(CHR, gene)][iso].get('cds', None)
                attri = annot[(CHR, gene)][iso].get('attri', None)
                if not any([exon, cds]):
                    exon = [(start, end)]
                elif not exon:
                    exon = cds.copy()
                if cds:
                    cds = sorted(cds, key=lambda x: x[0], reverse=False)
                    cds = (cds[0][0], cds[-1][-1])
                yield Transcript(iso, chro, start, end, strand, exon, cds, attri)

    def __parse_gff(self):
        RegulateRegion = {
            'DNAseI_hypersensitive_site',
            'enhancer',
            'enhancer_blocking_element',
            'insulator',
            'promoter',
            'protein_binding_site',
            'replication_regulatory_region',
            'transcriptional_cis_regulatory_region',
        }
        logger.info((
            'skip annotation feature as Regulating region: \n'
            '    ({})'.format(', '.join(RegulateRegion))
        ))

        MotifRegion = {
            'centromere',
            'direct_repeat',
            'microsatellite',
            'tandem_repeat',
            'mobile_genetic_element',
            'nucleotide_motif',
            'repeat_instability_region',
        }
        logger.info((
            'skip annotation feature as Motif sequence region: \n'
            '    ({})'.format(', '.join(MotifRegion))
        ))

        UnknownRegion = {
            'cDNA_match',
            'repeat_region',
            'D_loop',
            'match',
            'region',
            'origin_of_replication',
            'sequence_feature',
            'biological_region',    # Igl
            'sequence_alteration',
            'CAGE_cluster',
            'meiotic_recombination_region',
            'mitotic_recombination_region',
        }
        logger.info((
            'skip annotation feature of unknown: \n'
            '    ({})'.format(', '.join(UnknownRegion))
        ))
        SkipRegion = UnknownRegion | MotifRegion | RegulateRegion

        TRANSCRIPT = {
            'mRNA', 'lnc_RNA', 'ncRNA',
            'antisense_RNA',
            'transcript',  # non coding
            'RNase_MRP_RNA',
            'RNase_P_RNA',
            'Y_RNA', 'tRNA', 'rRNA', 'snoRNA', 'snRNA',
            'primary_transcript',  # miRNA precursor seq
            'miRNA',
            'C_gene_segment',
            'D_gene_segment',
            'V_gene_segment',
            'J_gene_segment',
            'SRP_RNA',
            'telomerase_RNA',
            'vault_RNA',
            'guide_RNA',
            'scRNA',
        }
        logger.info((
            'the following annotation feature as transcript: \n'
            '    ({})'.format(', '.join(TRANSCRIPT))
        ))

        logger.info('Start Read gff file to python object...')
        Gidlst, GAnnot, TAnnot, Tstructure = [], {}, {}, {}
        for line in Files.__iter__(self):
            if line.startswith('#') or not line.strip():
                continue
            chro, _, feature, start, end, _, strand, _, attr = line.strip().split('\t')
            if feature in SkipRegion:
                continue

            start, end = int(start) - 1, int(end)
            attr = [i.strip().partition('=') for i in attr.split(';') if i.strip()]
            attr = {i[0]: i[2].strip('\'\t\n\r"') for i in attr}
            Dbxref = [i.partition(':') for i in attr.pop('Dbxref', '').split(',')]
            Dbxref = {i[0]: i[2].strip('\'\t\n\r"') for i in Dbxref}
            attr.update(Dbxref)

            if feature == 'gene' or feature == 'pseudogene':
                line_id = attr.pop('ID')
                Gidlst.append((chro, line_id))
                GAnnot.setdefault((chro, line_id), {})['infor'] = [chro, start, end, strand]
                attr['gene_name'] = attr.pop('Name')
                attr['gene_type'] = attr.pop('gene_biotype')
                GAnnot[(chro, line_id)]['attri'] = attr
            elif feature in TRANSCRIPT:
                g_id = attr.pop('Parent')
                t_id = attr.pop('ID')
                Gidlst.append((chro, g_id))
                TAnnot.setdefault((chro, g_id), {}).setdefault(t_id, {})['infor'] = [chro, start, end, strand]
                if feature == 'miRNA':
                    attr['transcript_id'] = attr.pop('product')
                    attr['transcript_name'] = attr.pop('gene')
                    attr['transcript_type'] = 'miRNA'
                elif feature == 'tRNA':
                    attr['transcript_id'] = t_id
                    attr['transcript_name'] = t_id
                    attr['transcript_type'] = 'tRNA'
                else:
                    attr['transcript_id'] = attr.get('transcript_id', t_id)
                    attr['transcript_name'] = attr.get('Name', t_id)
                    tmp = attr['gbkey']
                    attr['transcript_type'] = 'protein_coding' if tmp == 'mRNA' else tmp
                TAnnot[(chro, g_id)][t_id]['attri'] = attr
            elif feature == 'exon':
                t_id = attr.pop('Parent')
                Tstructure.setdefault((chro, t_id), {}).setdefault('exon', []).append((start, end))
            elif feature == 'CDS':
                t_id = attr.pop('Parent')
                Tstructure.setdefault((chro, t_id), {}).setdefault('cds', []).append((start, end))
            else:
                logger.warning('skip novel annotation feature: {}'.format(feature))

        logger.info('Done of Read gff, sort transcript list...')
        Gidlst = sorted(
            set(Gidlst), key=lambda x: Gidlst.index(x), reverse=False
        )
        logger.info((
            'Done of sort transcript list, convert raw annot infor to '
            'appropriate transcript structure annotation...'
        ))
        annot = {}
        for index, (CHR, gene) in enumerate(Gidlst):
            if index % 5000 == 0 and index != 0:
                logger.info('Already processed gene NO. : {}'.format(index))

            infor_gene = GAnnot.get((CHR, gene), None)
            infor_transcript = TAnnot.get((CHR, gene), None)
            if infor_transcript is None:
                logger.debug('Missing transcript feature line of gene: {}'.format(gene))
                annot.setdefault((CHR, gene), {}).setdefault(gene, {})['infor'] = infor_gene['infor']
                annot[(CHR, gene)][gene]['attri'] = infor_gene['attri']
                annot[(CHR, gene)][gene]['exon'] = Tstructure.get((CHR, gene), {}).get('exon', None)
                annot[(CHR, gene)][gene]['cds'] = Tstructure.get((CHR, gene), {}).get('cds', None)
            else:
                for iso in infor_transcript:
                    chro, start, end, strand = infor_transcript[iso]['infor']
                    attri = infor_transcript[iso]['attri']
                    if infor_gene:
                        G_chro, G_start, G_end, G_strand = infor_gene['infor']
                        msg = 'transcript({}) region over range of gene({}) region'.format(iso, gene)
                        assert G_start <= start < end <= G_end, msg
                        msg = 'chro is inconsistent of transcript({}) and gene({})'.format(iso, gene)
                        assert chro == G_chro, msg
                        msg = ('Transcription direction(strand) is inconsistent of '
                               'transcript({}) and gene({})').format(iso, gene)
                        assert strand == G_strand, msg
                        Gattr = infor_gene['attri']
                        attri.update(Gattr)
                    annot.setdefault((chro, gene), {}).setdefault(iso, {})['infor'] = infor_transcript[iso]['infor']
                    if attri.get('gene_type', '') == 'other':
                        attri['gene_type'] = attri['transcript_type']
                    annot[(chro, gene)][iso]['attri'] = attri
                    structure = Tstructure.get((chro, iso), Tstructure.get((chro, gene), {}))
                    annot[(chro, gene)][iso]['exon'] = structure.get('exon', None)
                    annot[(chro, gene)][iso]['cds'] = structure.get('cds', None)
        logger.info('Done of parse gff.')
        return Gidlst, annot


class BedReader(Files):
    '''
    '''
    def __init__(self, bed):
        Files.__init__(self, bed)

    def __str__(self):
        return "<BedReader object: \"{}\">".format(self._fos)

    def __iter__(self):
        for line in Files.__iter__(self):
            chro, start, end, name, _, strand, cstart, cend, _, _, elen, estart = line.strip().split('\t')
            start, end, cstart, cend = map(int, [start, end, cstart, cend])
            elen = [int(i) for i in elen.split(',') if i]
            estart = [int(i) + start for i in estart.split(',') if i]
            assert len(estart)==len(elen), 'Error line: {}'.format(line)

            eend = [estart[i] + elen[i] for i, j in enumerate(elen)]
            EXON = [(estart[i], eend[i]) for i, j in enumerate(estart)]
            if cstart == cend:
                CDS = []
                attri = {'transcript_type': 'noncoding'}
            else:
                CDS = (cstart, cend)
                attri = {'transcript_type': 'protein_coding'}
            yield Transcript(name, chro, start, end, strand, EXON, CDS, attri)


class genePredReader(Files):
    '''
    '''
    def __init__(self, refFlat):
        Files.__init__(self, refFlat)

    def __str__(self):
        return "<genePredReader object: \"{}\">".format(self._fos)

    def __iter__(self):
        for line in Files.__iter__(self):
            try:
                Tname, chro, strand, start, end, cstart, cend, ecount, estart, eend = line.strip().split('\t')
                Gname = None
            except:
                Gname, Tname, chro, strand, start, end, cstart, cend, ecount, estart, eend = line.strip().split('\t')
            start, end, cstart, cend, ecount = map(int, [start, end, cstart, cend, ecount])
            estart = [int(i) for i in estart.split(',') if i]
            eend = [int(i) for i in eend.split(',') if i]
            assert len(estart)==len(eend)==ecount, 'Error line: {}'.format(line)

            EXON = [(estart[i], eend[i]) for i, j in enumerate(estart)]
            if cstart == cend:
                CDS = None
                attri = {'transcript_type': 'noncoding'}
            else:
                CDS = (cstart, cend)
                attri = {'transcript_type': 'protein_coding'}
            if Gname:
                attri['gene_id'] = Gname
                attri['gene_name'] = Gname
            yield Transcript(Tname, chro, start, end, strand, EXON, CDS, attri)