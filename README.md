# sv-merge #
A python tool to perform sample-specific merging of SV vcf files by combining events with matching breakends given a slack parameter (e.g. 200bp) and reciprical overlap value (e.g. 0.8).  Creates combined .txt and .vcf files.
## Requirements
* python3
* python packages: pandas, pyranges, PyVCF

## Usage
* download this github repo
* to run: sv-merge/sv-merge.py --vcfs sample.manta.vcf,sample.svaba.vcf,sample.gridds.vcf --sample-name sample
* creates 3 types of files:
  1. caller.txt files: tab-delimited text file with all standard chromosome calls 
  2. match.txt files: file containing matched calls between different callers
  3. sv-merge.txt files: file containing the merged calls where calls by two or more callers have been combined

![image](https://github.com/bankhead3/sv-merge/assets/31142967/55ce5a82-5684-4890-86de-3d4f4c06cd81)

## Assumptions and caveats to be aware of...
* Currently: manta,svaba,gridss calls supported 
* Non-standard chromosome (e.g. alt, random, Unknown) are excluded
* Merged events require that both breakends match within 200bp and reciprical overlap of 0.8 if intra-chromosomal
* svs are matched using a simple greedy graph-based strategy fully connected calls or cliques are merged together
* manta/gridss del,ins,tandem dup events not merged due to missing alt orientation information.  All events are parsed and written to caller.txt files.
* If something doesn't work the way it should - let's discuss and fix

