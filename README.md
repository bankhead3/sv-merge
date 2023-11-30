# sv-merge #
A tool to perform sample-specific merging of SV vcf files by taking the intersection of variant events using a slack allowance for matching.  Generates combined table with all calls where matching calls have been combined.
## Requirements
* python3
* python packages: pandas, pyranges, PyVCF

## Usage
* download this github repo
* to run: sv-merge/sv-merge.py --vcfs sample.manta.vcf,sample.svaba.vcf --sample-name sample
* creates 3 types of files:
  1. caller.txt files: tab-delimited text file with all standard chromosome calls 
  2. comparison.txt files: file containing matched calls between different callers
  3. sv-merge.txt files: file containing the merged calls where calls by two or more callers have been combined

![image](https://github.com/bankhead3/sv-merge/assets/31142967/55ce5a82-5684-4890-86de-3d4f4c06cd81)

## Assumptions and caveats to be aware of...
* Currently: only manta and svaba calls are supported 
* Non-standard chromosome (e.g. alt, random, Unknown) are excluded
* Merged events require that both breakends match within 200bp
* If something doesn't work the way it should - let's discuss and fix

