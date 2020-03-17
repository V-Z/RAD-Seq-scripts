# References FASTA should be placed here

Processing RAD data require indexed FASTA reference files. See steps below.

## Required software

* BWA <https://github.com/lh3/bwa>
* Picard <https://broadinstitute.github.io/picard/>
* Samtools <http://www.htslib.org/>

## Indexing the reference

```sh
REF='file' # E.g. 'file.fasta' will be used as reference
# Create *.dict index
java -jar picard.jar CreateSequenceDictionary R="$REF".fasta O="$REF".dict
# Create *. fasta.amb, *. fasta.ann, *. fasta.bwt, *.fasta.pac, *.fasta.sa indices
bwa index "$REF".fasta
# Create *.fasta.fai index
samtools faidx "$REF".fasta
```

