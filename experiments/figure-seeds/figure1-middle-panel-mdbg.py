import statistics, sys
for line in open(sys.argv[1]):
    abundances = list(map(int,line.split()[1:]))
    if len(abundances) == 0:
        print(0,0,0)
    else:
        print(statistics.median(abundances),statistics.mean(abundances),statistics.stdev(abundances))
