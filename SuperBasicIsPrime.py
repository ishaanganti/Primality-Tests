import numpy as np
import time, math, random

#Elementary isprime method
def SuperBasicIsPrime(x):
    for i in range(2, x-1): #Checking divisibility of input by every odd number from 3 to sqrt(input) int(numpy.sqrt(x))
        if x == 2: return True
        if x % i == 0:
            return False                       
    return True       #If prime

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
        state = SuperBasicIsPrime(num)
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

print(MethodTimer(23))
