# intersects two PAF files and determines the concordance of results, similarly to paftools mapeval

import sys
if len(sys.argv) < 3: 
    sys.stderr.write("arguments: [paf1] [paf2]\n")
    exit(1)
paf1_filename = sys.argv[1]
paf2_filename = sys.argv[2]


def parse_paf(filename):
    p = dict()
    for line in open(filename):
        #S1_1!chr1!224752794!224777027!+ 24299   0       24298   +       chr1    248387328       224752793       224777027       132     248387328       60
        ls = line.split()
        read = ls[0]
        chrm = ls[5]
        start = int(ls[6])
        end = int(ls[7])
        #print(read,chrm,start,end)
        p[read]=(chrm,start,end)
    return p


sys.stderr.write("Loading PAFs..\n")
paf1 = parse_paf(paf1_filename)
sys.stderr.write(f"Loaded {paf1_filename}\n")
paf2 = parse_paf(paf2_filename)
sys.stderr.write(f"Loaded {paf2_filename}\n")

#print(paf1,paf2)

nb_concordant = 0
nb_diff_chr   = 0
nb_discordant = 0
def intersect(read,coords1,coords2):
    global nb_concordant, nb_discordant, nb_diff_chr
    chr1, start1, end1 = coords1
    chr2, start2, end2 = coords2
    if chr1 != chr2:
        nb_diff_chr += 1
        nb_discordant += 1
    else:
        lowest  = min(start1,start2,end1,end2)
        highest = max(start1,start2,end1,end2)
        # same as mapeval as per https://github.com/lh3/minimap2/blob/master/misc/README.md#evaluating-mapping-accuracy-with-simulated-reads
        max1 = max(start1,end1)
        min1 = min(start1,end1)
        max2 = max(start2,end2)
        min2 = min(start2,end2)
        if max1 < max2:
            # min1          max1
            # ---------------
            #      ---------------
            #      min2          max2
            if max1 >= min2:
                o = max1-min2
            else:
                o = 0
        else:
            # min2          max2
            # ---------------
            #      ---------------
            #      min1          max1
            if max2 >= min1:
                o = max2-min1
            else:
                o = 0
        r = o/(highest-lowest)
        if r > 0.1:
            nb_concordant += 1
        else:
            nb_discordant += 1


paf1_missing_in_paf2 = set()
for read1 in paf1:
    coords1 = paf1[read1]
    if read1 in paf2:
        coords2 = paf2[read1]
        intersect(read1,coords1,coords2)
    else:
        paf1_missing_in_paf2.add(read1)


paf2_missing_in_paf1 = set()
for read2 in paf2:
    coords2 = paf2[read2]
    if read2 in paf1:
        # already processed
        #coords1 = paf1[read2]
        #intersect(read2,coords1,coords2)
        pass
    else:
        paf2_missing_in_paf1.add(read2)

print(f"Total number of mapped reads in {paf1_filename}: {len(paf1)}")
print(f"Total number of mapped reads in {paf2_filename}: {len(paf2)}")
print(f"Number of concordant mappings: {nb_concordant} ({nb_concordant/len(paf1)*100}% of {paf1_filename}, {nb_concordant/len(paf2)*100}% of {paf2_filename})")
print(f"Number of discordant mappings on same      chromosome: {nb_discordant} ({nb_discordant/len(paf1)*100}% of {paf1_filename}, {nb_discordant/len(paf2)*100}% of {paf2_filename})")
print(f"Number of discordant mappings on different chromosome: {nb_diff_chr}")
