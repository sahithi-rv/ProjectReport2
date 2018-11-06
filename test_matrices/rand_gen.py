import random
import sys

n = int(sys.argv[1])
filea = str(sys.argv[2])

f = open(filea, "w")
for i in range(n):
	a = random.random();
	f.write(str(a))
	f.write("\n")
