import argparse

parser = argparse.ArgumentParser(description = "Replace variant id in gen file with CHR:POS_REF_ALT.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--gen", help = "Path to the gen file.")
parser.add_argument("--chr", help = "Chromosome name.")
args = parser.parse_args()

for line in open(args.gen):
    line = line.rstrip()
    fields = line.split(" ")
    new_id = args.chr + ":" + fields[2] + "_" + fields[3] + "_" + fields[4]
    fields[1] = new_id
    new_line = " ".join(fields)
    print(new_line)
