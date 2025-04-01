# PupiPlot 
## Pipeline for Uncovering PileUp Information and plotting
PupiPlot is an improved version of genomealigner, designed to simplify the process of aligning sequencing reads to a reference genome and visualizing coverage depth. 

## Features
- Aligns FASTQ reads to a reference genome using `minimap2`
- Converts and sorts alignment files with `samtools`
- Generates a coverage depth plot from the BAM file
- Saves all output files (BAM, BAI, and PNG) in a user-defined directory
- Displays an adorable ASCII bunny upon successful completion

## Installation
Make sure you have the required dependencies installed:

```bash
pip install pysam matplotlib numpy
```

You'll also need `minimap2` and `samtools`, which can be installed via:

```bash
sudo apt install minimap2 samtools
```

## Usage
Run the script with:

```bash
python PupiPlot.py
```

Follow the prompts to provide the necessary input files and output folder. Pupi will handle the alignment for you! 🐰

## Output Files
All output files will be stored in the specified output folder:
- `alignment.sorted.bam` - The sorted BAM file
- `alignment.sorted.bam.bai` - The BAM index file
- `coverage_plot.png` - The coverage depth visualization

## Meet Pupi, the hardworking bunny!
```
  /\ /|
  \ V/
  | "")
  /  |
 /  \\
*(__\_\)
```
Pupi is the hardworking bunny who assembles your genome with precision and care. 
## License
This project is licensed under the MIT License.

## Acknowledgments
Special thanks to Pupi, the alignment bunny, for tirelessly processing genome data! 🐇✨
