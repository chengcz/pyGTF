## Pure python parser of Fastx, GTF, NCBI GFF files

parse universal GTF/GFF file, return Transcript object, convert annotation infor as GTF, BED, GenePred format, and extract genome, transcript, CDS and UTR sequence with reference genome file

### install

```bash
git clone git@github.com:chengcz/pyGTF.git
cd pyGTF

python setup.py install
# done, enjoy
```



### basic usage

1. GTFReader

   ```python
   import sys
   from pyGTF import Transcript, GTFReader
   
   for i in GTFReader('data/ath.chr1.gff.gz'):
       print(i.summary)
       print('format1: bed')
       i.to_bed(sys.stdout)
       print('format2: gtf')
       i.to_gtf(sys.stdout)
       break
   dir(t)    # check Transcript object attribution and methods
   
   ##########
   ### GTFReader.safe_mode accept unsort GTF/GFF file
   ##########
   fp = GTFReader('data/ath.chr1.gff.gz')
   for i in fp.safe_mode():
     print(i)
     break
   ```

   

2. RefSeqGFFReader

   specific parse GFF file from NCBI, ref above

   ```python
   from pyGTF import Transcript, RefSeqGFFReader
   
   for i in RefSeqGFFReader('data/hg38.refseq.chr22.gff.gz'):
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
   
   dir(i)    # check Sequence object attribution and methods
   ```
