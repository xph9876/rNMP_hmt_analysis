# rNMP_hmt_analysis
## Description
Analyze the rNMP embedment characteristics in human mtDNA.
## Citation
Please use the following citation if you use the code:

Enriched zones of embedded ribonucleotides are associated with DNA replication and coding sequences in the human mitochondrial genome,
Penghao Xu, Taehwan Yang, Deepali L. Kundnani, Mo Sun, Stefania Marsili, Alli L. Gombolay, Youngkyu Jeon, Gary Newnam, Sathya Balachander, Veronica Bazzani, Umberto Baccarani, Vivian S. Park, Sijia Tao, Adriana Lori, Raymond F. Schinazi, Baek Kim, Zachary F. Pursell, Gianluca Tell, Carlo Vascotto, Francesca Storici
_bioRxiv_ 2023.04.05.535745; doi: https://doi.org/10.1101/2023.04.05.535745

## Dependency
A conda environments containing all dependencies are available in __wrapper__ folder, to install:
```bash
conda create -f wrapper/env.yaml
conda activate rNMP_mtDNA
```
[RibosePreferenceAnalysis](https://github.com/xph9876/RibosePreferenceAnalysis) is required to run the scripts
```bash
git clone https://github.com/xph9876/RibosePreferenceAnalysis.git
cd RibosePreferenceAnalysis
pwd
```
copy the output path to the RibosePreferenceAnalysis variable in __wrapper/analyze.sh__

The curlyBrace function from Dr. Siyu Gao is used in control region analysis.

## Usage
Open the __wrapper/analyze.sh__, update the variables at the beginning, then run it
```bash
./wrapper/analyze.sh
```

## Contact
pxu64@gatech.edu
