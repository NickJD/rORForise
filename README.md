# rORForise
Read-based gene coverage evaluation

## Test data is available in the Genome_Processing directory.

### Menu - (rORForise -h):  
```bash
Running rORForise
usage: rORForise.py [-h] -d DIRECTORY -g GENOMES [GENOMES ...] --gc_prob GC_PROB -f FRAGMENTATION_TYPES [FRAGMENTATION_TYPES ...] -s SUBGROUPS [SUBGROUPS ...] -m METHODS
                    [METHODS ...]

Process genome evaluation parameters.

options:
  -h, --help            show this help message and exit
  -d DIRECTORY, --directory DIRECTORY
                        Path to the main directory containing genome processing subdirectories.
  -g GENOMES [GENOMES ...], --genomes GENOMES [GENOMES ...]
                        List of genome names to process. Example: 'Mycoplasma_genitalium_G37 Staphylococcus_aureus_502A'
  --gc_prob GC_PROB     GC probability of the genome being processed (e.g., 0.3169 for Mycoplasma genitalium).
  -f FRAGMENTATION_TYPES [FRAGMENTATION_TYPES ...], --fragmentation_types FRAGMENTATION_TYPES [FRAGMENTATION_TYPES ...]
                        List of fragmentation types to use (e.g., 'ART_errFree').
  -s SUBGROUPS [SUBGROUPS ...], --subgroups SUBGROUPS [SUBGROUPS ...]
                        List of subgroups for evaluation (e.g., 'Combined').
  -m METHODS [METHODS ...], --methods METHODS [METHODS ...]
                        List of methods for evaluation (e.g., 'FragGeneScan').
```
## Example Run:
```bash
rORForise -d "../../Genome_Processing" -g "Mycoplasma_genitalium_G37" -f "ART_errFree" -s "Combined" -m "FragGeneScan" --gc_prob 0.3169
```
```bash
Running rORForise
Mycoplasma_genitalium_G37
FragGeneScan
ART_errFree
../../Genome_Processing/Mycoplasma_genitalium_G37/FragGeneScan *_ART_errFree_Combined.gff.gz ['../../Genome_Processing/Mycoplasma_genitalium_G37/FragGeneScan/FragGeneScan_ART_errFree_Combined.gff.gz']
../../Genome_Processing/Mycoplasma_genitalium_G37/FragGeneScan/FragGeneScan_ART_errFree_Combined.gff.gz
Number of predictions made for reads    69709
Number of CDS-aligned reads     70073
Number of CDS-aligned reads without predictions 5929
Number of CDS-aligned reads with at least one pred      64440
Number of predictions   69430
Number of on-target preds       64695
correct start   1211
correct stop    2030
correct frame   57249
correct direction       60168
incorrect stop  13624
incorrect start 8015
incorrect frame 2919
incorrect direction     4527
prediction ends before cds starts       660
prediction starts after cds ends        713
alternative start       1837
middle or alternative start     1913
alternative stop        5201
middle or alternative stop      1166
middle  49700
... Below is the output of the Codon predictions ...
```
