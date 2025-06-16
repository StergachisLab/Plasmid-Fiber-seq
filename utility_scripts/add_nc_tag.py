import pysam
import argparse

def add_nc_tag(input_bam_path, output_bam_path):
    # Open the input BAM file
    input_bam = pysam.AlignmentFile(input_bam_path, "rb")
    
    # Open the output BAM file
    output_bam = pysam.AlignmentFile(output_bam_path, "wb", template=input_bam)
    
    # Iterate through each alignment in the input BAM file
    for read in input_bam:
        if read.has_tag('ns'):
            nuc_count = len(read.get_tag('ns'))
        else:
            nuc_count = 0
        
        # Add a new tag to the alignment
        read.set_tag("nc", value=nuc_count, value_type="i")
        
        # Write the modified alignment to the output BAM file
        output_bam.write(read)
    
    # Close the BAM files
    input_bam.close()
    output_bam.close()

def main():
    parser = argparse.ArgumentParser(description="Add a nucleosome count tag ('nc') to alignments in a BAM file.")
    parser.add_argument("input_bam", help="Path to the input BAM file")
    parser.add_argument("output_bam", help="Path to the output BAM file")

    args = parser.parse_args()
    
    add_nc_tag(args.input_bam, args.output_bam)

if __name__ == "__main__":
    main()
