import sys

sequence_file = sys.argv[1]

with open(sequence_file, 'r') as f:
    # print("tell(): {}".format(f.tell()))
    start = f.tell()
    line = f.readline()
    taxid = 0
    geneid = 0
    while line:
        # print("tell(): {}".format(f.tell()))
        line = line.rstrip()

        # print("Start {} end {}".format(start, f.tell()))

        if line.startswith('>'):
            tokens = line[1:].split('_')

            taxid = int(tokens[0])
            geneid = int(tokens[1])
        else:
            print("{}\t{}\t{}\t{}".format(taxid, geneid, start, f.tell()-1))

        start = f.tell()
        line = f.readline()

