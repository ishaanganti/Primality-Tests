import time, math, random

def MillerRabin(n):
    return MR(n, 6)

def MR(n, k):
    if n < 2: return False
    if n < 4: return True
    if n % 2 == 0: return False
    #simple cases: n is less than 2, n is 3, or n is even.

    #n must be odd and at least 5
    s = 0 #will be number of 2s in prime factorization of n-1
    d = n-1 #product of remaining prime factors of n-1
    while d % 2 == 0:
        s += 1
        d //= 2
    # n-1 = 2^s * d, where d is odd

    for i in range(k): #run the test k times
        a = random.randrange(2, n-1)    # 2 <= a <= n-2 
        x = pow(a, d, n) #sets x to a^d mod n 
        if x == 1: continue #condition 1
        for j in range(s):
            if x == n-1: break #condition 2
            x = (x * x) % n
        else:
            return False
    return True
#Method to generate a random odd number    
def rng(bits):
    num = random.getrandbits(bits)
    while(num % 2 == 0):
        num = random.getrandbits(bits)
    return num

#Method to generate a prime using the basic primality test
def GenPrime(bits):
    num = 0
    state = False
    start = time.perf_counter()
    while(state == False):
        num = rng(bits)
        state = MillerRabin(num)
    end = time.perf_counter()
    print(num)
    return (end-start)

#Measures primality generation speed for different specified bit sizes, stores in array
def MethodTimer(bitRange):
    arr = [0]
    total = 0
    for i in range(2, bitRange):
        for k in range(25):
            total += GenPrime(i)
        arr.append(total/25)
        total = 0
    return arr

print(MethodTimer(20))
#Generally, 6 rounds of the Miller Rabin Test is good for numbers less than or equal to 1024 bits.
