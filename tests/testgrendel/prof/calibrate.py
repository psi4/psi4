
import profile

def meanstdv(x):
    from math import sqrt
    n, mean, std = len(x), 0, 0
    for a in x:
        mean = mean + a
        mean = mean / float(n)
    for a in x:
        std = std + (a - mean)**2
        std = sqrt(std / float(n-1))
    return mean, std


results = []
pr = profile.Profile()
for i in range(10):
    results.append(pr.calibrate(50000))
mean, stdev = meanstdv(results)
print "Average: {}\nMax: {}\nMin: {}\nStdDev: {}".format(mean, max(*results), min(*results), stdev)

