# requires to know the strand of each read. for that, since wgsim doesn't output it, I just use minimap:
# minimap2 ecoli.genome.100k.fa perfect-ecoli.10X.100k.fa > minimap2.paf

import sys

if len(sys.argv) < 3:
    print("arguments: wgsim.fastq minimap2.paf")
    exit(1)

read_strand = dict()
for line in open(sys.argv[2]):
    #S1_1!NZ_CP027599.1!10587!12039!+        1456    3       1450    +       NZ_CP027599.1   99925   10590   12033   1408    1447    60      tp:A:P  cm:i:249        s1:i:1408       s2:i:0  dv:f:0.0036     rl:i:0
    strand=line.split()[4]
    read=line.split()[0]
    read_strand[read]=strand

counter=0
for line in open(sys.argv[1]):
    if line.startswith("@"):
       #@NZ_CP027599.1_88140_93139_0:0:0_0:0:0_2/1
       ref1, ref2,start,end = line.split('_')[:4]
       ref=ref1+"_"+ref2 # assumes the ref has a "_" (critical assumption)
       strand = read_strand[line.strip()[1:]]
       # >S1_5!NZ_CP027599.1!45234!51766!+
       print(f"@S1_{counter}!{ref}!{start}!{end}!{strand}")
       counter+=1
    else:
        print(line.strip())
