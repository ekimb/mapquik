import statistics, sys
for line in open(sys.argv[1]):
    if line.startswith('>'): continue
    abundances = list(map(int,line.split()))
    print(statistics.median(abundances),statistics.mean(abundances),statistics.stdev(abundances))
