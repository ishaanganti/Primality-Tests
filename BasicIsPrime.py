import numpy as np
import time, math, random

#Elementary isprime method
def BasicIsPrime(x):
    if x == 2: #Checking if the input is 2
        return True
    if x % 2 == 0: #Checking if the input is even
        return False
    for i in range(3, int(math.sqrt(x)), 2): #Checking divisibility of input by every odd number from 3 to sqrt(input) int(numpy.sqrt(x))
        if x % i == 0:
            return False
    else:                        #If prime
        return True

#Method to generate a random odd number    
def rng(bits):
    num = random.randrange(1 << bits, 1 << bits+1)
    while(num % 2 == 0):
        num = random.randrange(1 << bits, 1 << bits+1)
    #print(str(num) + "yo")
    return num

#Method to generate a prime using the basic primality test
def GenPrime(bits):
    num = 0
    state = False
    start = time.perf_counter()
    while(state == False):
        num = rng(bits)
        state = BasicIsPrime(num)
    end = time.perf_counter()
    return (end-start)

#Measures primality generation speed for different specified bit sizes, stores in array
def MethodTimer(bitRange):
    arr = [0]
    total = 0
    for i in range(3, bitRange):
        for k in range(20):
            total += GenPrime(i)
        arr.append(total/20)
        total = 0
    return arr

print(MethodTimer(48))
