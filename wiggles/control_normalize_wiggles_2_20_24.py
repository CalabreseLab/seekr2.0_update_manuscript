import os

# Input received as serial cluster submission for each exp:control pairing
# python3 control_normalize_wiggles_2_20_24.py </path/to/exp.wig> </path/to/control.wig>
# Input for eCLIPs automatically developed with setup_normalize_wigs_2_20_24.sh


def read_in_wiggle_info(wiggle_file, initial_header):
    wiggle_dict = {}
    with open(wiggle_file, "r") as f:
        for i,line in enumerate(f):
            if i == 0:
                initial_header = line
            elif "variableStep" in line:
                chrom = line.split("chrom=")[1].split(" ")[0]
                span = line.split("span=")[1].strip()
                wiggle_dict[chrom] = {}
                continue
            else:
                position, signal = line.strip().split("\t")
                wiggle_dict[chrom][position] = signal
    return(wiggle_dict, initial_header, span)

def normalize_wiggles(exp_dict, control_dict):
    adj_dict = {}
    for chrom in exp_dict:
        adj_dict[chrom] = {}
        for position in exp_dict[chrom]:
            if position not in control_dict[chrom]:
                # this means that the control signal was zero
                adj_dict[chrom][position] = float(exp_dict[chrom][position])
            elif control_dict[chrom][position] != 0:
                adj_dict[chrom][position] = float(exp_dict[chrom][position]) - float(control_dict[chrom][position])
                if adj_dict[chrom][position] <= float(0):
                    del adj_dict[chrom][position]
                else:
                    adj_dict[chrom][position] = str(adj_dict[chrom][position])
    return(adj_dict)

def write_wiggle_file(outfile_name, adj_dict, initial_header, wig_span):
    with open(outfile_name, "w") as f:
        f.write(initial_header)
        for chrom in adj_dict:
            f.write("variableStep chrom=" + chrom + " span=" + wig_span + "\n")
            for position in adj_dict[chrom]:
                f.write(str(position) + "\t" + str(adj_dict[chrom][position]) + "\n")


if __name__ == "__main__":
    # check input validity
    if len(os.sys.argv) != 3:
        print("Usage: python control_normalize_wiggles_2_20_24.py <exp_file> <control_file>")
        print(os.sys.argv)
        os.sys.exit(1)

    # define inputs
    exp_file = os.sys.argv[1]
    control_file = os.sys.argv[2]
    control_map = {}

    # obtain paths for both file types
    exp_path = exp_file.split("/")[0:-1]
    exp_path = "/".join(exp_path)
    exp_path = exp_path + "/"

    # read in the wiggle files
    exp_dict, initial_header, wig_span = read_in_wiggle_info(exp_file, "")
    control_dict, null_var, wig_span_cont = read_in_wiggle_info(control_file, "")
    if wig_span != wig_span_cont:
        print("ERROR: The spans of the two wiggle files are not the same")
        os.sys.exit(1)

    # adjust wiggles by controls
    adj_dict = normalize_wiggles(exp_dict, control_dict)

    # check for output directory
    if not os.path.exists("control_normalized_wiggles"):
        os.makedirs("control_normalized_wiggles")

    # write out the adjusted wiggle file
    outfile_name = exp_file.strip(".wig")
    outfile_name = "control_normalized_wiggles/" + outfile_name.replace(exp_path, "")
    outfile_name = outfile_name + "_normalized.wig"
    write_wiggle_file(outfile_name, adj_dict, initial_header, wig_span)
