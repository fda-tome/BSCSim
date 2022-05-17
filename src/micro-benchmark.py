import numpy as np
import time

n = 100
start = time.time()
while n <= 100000000:
    acc = 0
    A = np.random.rand(n)
    for k in A:
        acc += k
    n *= 10
    print(time.time() - start)
    start = time.time()
