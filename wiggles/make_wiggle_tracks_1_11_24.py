import sys
import math

# Import chrM warning: If you want to study chromosome M, note that if there are reads in the first 0-50 nucleotides those are not displayed,
# and the last 50 nucleotides of reads are also not displayed. This is due to problems with visual cutoffs and chromosome lengths.
# If you need coverage in these regions, please adjust the script accordingly.

# Parameters:
#1. Bed12 file
#2. chrNameLength.txt file generated in the genome build step of STAR
#3. Header for output files to be displayed in genome browser
#4. Color of the wiggle track
#5. Whether to log10 normalize the data: y for log10 n for no normalization
#6. Bin size for the wiggle track
#7. Number of reads in the dataset (per wiggle track) for RPM standardization. If you don't wish to standardize by RPM, put the value 1 here.

# ex: python3 make_wiggle_tracks.py <bed12> <chrNameLength.txt> test-output blue n 50 1

# import variables
bedfile = sys.argv[1]
chrNameLength = sys.argv[2]
header = sys.argv[3]
color = sys.argv[4]
log10 = sys.argv[5]
bin_size = int(sys.argv[6])
reads_in_dataset = int(sys.argv[7])


if reads_in_dataset > 1:
    rpm = reads_in_dataset/1000000
else:
    rpm = 1


wig_dict = {}

with open(chrNameLength, "r") as infile:
    for line in infile:
        cols = line.split("\t")
        chrom = cols[0]
        if "GL" not in chrom and "JH" not in chrom and "KI" not in chrom:
            length = cols[1]
            wig_dict[chrom] = {}

# make the wiggle tracks
with open(bedfile, "r") as infile:
    for line in infile:
        cols = line.split("\t")
        chrom = cols[0]
        if chrom in wig_dict.keys():
            start = cols[1]
            strand = cols[5]
            block_sizes = cols[10].split(",")
            block_starts = cols[11].split(",")

            # find bin for each block of read
            my_covered_bins = []

            for i in range(len(block_sizes)):
                block_size = int(block_sizes[i])
                block_start = int(block_starts[i])
                block_end = block_start + block_size
                current_blockstart_bin = int(((int(start)+block_start)/bin_size))
                current_blockend_bin = int(((int(start)+block_end)/bin_size))

                for j in range(current_blockstart_bin, current_blockend_bin+1):
                    my_covered_bins.append(j)
                
            fractional_coverage = 1/len(my_covered_bins)
            for current_bin in my_covered_bins:
                if current_bin not in wig_dict[chrom]:
                    wig_dict[chrom][current_bin] = fractional_coverage
                else:
                    wig_dict[chrom][current_bin] += fractional_coverage

# make three number code for color
color_dict = {
    "blue":"0,0,255",
    "red":"255,0,0",
    "green":"0,255,0",
    "yellow":"255,255,0",
    "orange":"255,165,0",
    "purple":"128,0,128",
    "pink":"255,192,203",
    "black":"0,0,0",
}

color_code = color_dict[color]

# Need to sort the order of the bins for output
for chrom in wig_dict:
    wig_dict[chrom] = dict(sorted(wig_dict[chrom].items()))


# delete the first and last bin in the chrM wig_dict dictionary
if 0 in wig_dict["chrM"].keys():
    del wig_dict["chrM"][0]

if wig_dict["chrM"]:
    last_bin = max(wig_dict["chrM"])
    del wig_dict["chrM"][last_bin]

# write the wiggle track file
with open(header + ".wig", "w") as outfile:
    outfile.write("track type=wiggle_0 visibility=full name=" + header + " description=" + header + " color=" + color_code +  " maxHeightPixels=128:40:11  group=\"user\" graphType=bar priority=100 viewLimits=0:2.8 autoscale=on\n")
    for chrom in wig_dict:
        outfile.write("variableStep chrom=" + chrom + " span=" + str(bin_size) + "\n")

        for bin in wig_dict[chrom]:
            bin_name = bin*50
            if log10 == "y":
                outfile.write(str(bin_name) + "\t" + str(math.log10(wig_dict[chrom][bin])/rpm) + "\n")
            else:
                outfile.write(str(bin_name) + "\t" + str(wig_dict[chrom][bin]/rpm) + "\n")      