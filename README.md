## Pure python parser of Fastx, GTF, NCBI GFF files

parse universal GTF/GFF file, return Transcript object, convert annotation infor as GTF, BED, GenePred format, and extract genome, transcript, CDS and UTR sequence with reference genome file

### download

```bash
git clone git@github.com:chengcz/pyGTF.git
cd pyGTF
```



### usage(script)

```bash
python pyGTF.py -h

# bed2gtf
python pyGTF.py --input input.bed --input-format bed --gtf output.gtf

# xxx
```



### usage(module)

##### install

```bash
python setup.py install
# done, enjoy
```



1. GTFReader

   ```python
   import sys
   from pyGTF import Transcript, GTFReader
   
   with GTFReader('data/ath.chr1.gff.gz', flag_stream=True) as fi:
     # flag_stream:    bool, parse stype of gtf, use less menory
       for i in fi:
     	# for i in GTFReader('data/ath.chr1.gff.gz', flag_stream=True):
           print(i.summary)
           print('format1: bed')
           i.to_bed(sys.stdout)
           print('format2: gtf')
           i.to_gtf(sys.stdout)
           break
   dir(t)    # check Transcript object attribution and methods
   ```

   

2. RefSeqGFFReader

   specific parse GFF file from NCBI, ref above

   ```python
   from pyGTF import Transcript, RefSeqGFFReader
   
   for i in RefSeqGFFReader('data/hg38.refseq.chr22.gff.gz', strict=False):
     	# strict： bool, whether to strictly verify the gene/exon interval
       print(i.summary)
      	break
   ```

   **important and Ignored feature is written inside the module script, adaptable changes that may be involved, please edit `RefSeqGFFReader.__parse_gff` function**

   

3. BedReader / genePredReader

   ```python
   from pyGTF import Transcript, BedReader, genePredReader
   
   for i in genePredReader('data/hg38.refFlat.chr22.txt.gz'):
       print(i.summary)
      	break
   ```

   

4. FastaReader / FastqReader

   ```python
   from pyGTF import FastqReader, FastqReader
   
   for i in FastaReader('data/ath.chr1.fa.gz'):
       print(i)
       break
   
   # for r1, r2 in FastqReader('fq1', 'fq2'):
   #   	pass
   ```
