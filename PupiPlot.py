import os
import subprocess
import pysam
import matplotlib.pyplot as plt
import numpy as np
import random

def run_command(command):
    """Function to execute a shell command"""
    print(f"Running command: {command}")
    subprocess.run(command, shell=True, check=True)

def get_coverage(bam_file):
    """Calculate coverage depth for the entire reference genome, including zero coverage."""
    samfile = pysam.AlignmentFile(bam_file, "rb")
    ref_length = samfile.lengths[0]  # length of the first (and assumed only) reference sequence
    coverage = np.zeros(ref_length, dtype=int)
    for pileupcolumn in samfile.pileup():
        pos = pileupcolumn.reference_pos
        if pos < ref_length:
            coverage[pos] = pileupcolumn.n
    samfile.close()
    return coverage

def smooth_data(coverage, window_size=50):
    """Function to smooth the data using a moving average"""
    if len(coverage) < window_size:
        return coverage  # no smoothing if data shorter than window
    return np.convolve(coverage, np.ones(window_size)/window_size, mode='valid')

def darken_color(color, factor=0.8):
    """Function to darken a hex color by a given factor"""
    r, g, b = [int(color[i:i+2], 16) for i in (1, 3, 5)]
    r, g, b = int(r * factor), int(g * factor), int(b * factor)
    return f'#{r:02x}{g:02x}{b:02x}'

def plot_coverage(coverage, smoothed_coverage, output_folder):
    """Function to generate and save a coverage plot"""
    colors = ['#264653', '#2A9D8F', '#F4A261', '#E76F51', '#0081A7', '#F77F00', '#FAA307', '#FF595E', '#9A031E', '#0F4C5C']
    selected_color = random.choice(colors)
    line_color = darken_color(selected_color, factor=0.8)
    
    plt.figure(figsize=(12, 6))
    plt.plot(range(len(smoothed_coverage)), smoothed_coverage, color=line_color, label='Smoothed Coverage')
    plt.fill_between(range(len(smoothed_coverage)), smoothed_coverage, color=selected_color, alpha=0.6)
    plt.title("Genome Coverage Depth", color='#333333')
    plt.xlabel("Position on Genome", color='#333333')
    plt.ylabel("Coverage Depth", color='#333333')
    plt.grid(True, color='#333333', alpha=0.3)
    plt.legend()
    
    output_path = os.path.join(output_folder, "coverage_plot.png")
    plt.savefig(output_path)
    plt.show()

def print_ascii_bunny():
    """Prints an ASCII art bunny to the terminal"""
    bunny = """
      /\\ /|
      \\ V/
      | "")
      /  |
     /  \\
    *(__\\_\\)
    """
    print("Pupi has assembled your genome for you!")
    print(bunny)

def main():
    # Ask the user for the output folder name
    output_folder = input("Enter the name of the output folder: ").strip()
    os.makedirs(output_folder, exist_ok=True)
    
    # Get user inputs
    reference_genome = input("Enter the path to the reference genome (reference.fasta): ").strip()
    input_fastqs = input("Enter the path(s) to the input FASTQ file(s) (comma separated if more than one): ").strip().split(',')
    combine_fastqs = input("Do you want to combine the FASTQ files? (yes/no): ").strip().lower()
    
    # Step 1: Index the reference genome
    reference_mmi = os.path.join(output_folder, "reference.mmi")
    run_command(f"minimap2 -d {reference_mmi} {reference_genome}")
    
    # Step 2: Combine FASTQ files if required
    combined_fastq = os.path.join(output_folder, "combined.fastq")
    if combine_fastqs == "yes":
        run_command(f"cat {' '.join([fastq.strip() for fastq in input_fastqs])} > {combined_fastq}")
    
    # Step 3: Align reads with loosened minimap2 parameters for better mapping
    if combine_fastqs == "yes":
        alignment_file = os.path.join(output_folder, "alignment.sam")
        run_command(
            f"minimap2 -x map-ont -a --secondary=yes -N 50 -p 0.5 -A1 -B2 -O2,32 -E1,0 -k11 -w5 "
            f"{reference_mmi} {combined_fastq} > {alignment_file}"
        )
    else:
        alignment_file = os.path.join(output_folder, "alignment.sam")
        run_command(
            f"minimap2 -x map-ont -a --secondary=yes -N 50 -p 0.5 -A1 -B2 -O2,32 -E1,0 -k11 -w5 "
            f"{reference_mmi} {input_fastqs[0].strip()} > {alignment_file}"
        )
    
    # Step 4: Convert SAM to sorted BAM
    sorted_bam_file = os.path.join(output_folder, "alignment.sorted.bam")
    run_command(f"samtools view -Sb {alignment_file} | samtools sort -o {sorted_bam_file}")
    
    # Step 5: Index the BAM file
    run_command(f"samtools index {sorted_bam_file}")
    
    # Step 6: Generate coverage plot for the whole genome
    coverage = get_coverage(sorted_bam_file)
    smoothed_coverage = smooth_data(coverage)
    plot_coverage(coverage, smoothed_coverage, output_folder)
    
    print(f"Pipeline completed successfully! All output files are stored in: {output_folder}")
    print_ascii_bunny()

if __name__ == "__main__":
    main()
