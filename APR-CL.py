import copy
import numpy as np
import time, math, random
from math import gcd

#Implementation inspired by jacobi sum wikipedia page,
#“Implementation of a New Primality Test” by H. Cohen and A.K. Lenstra, and
#FLINT (Fast Library for Number Theory) source code
#(https://github.com/wbhart/flint2/tree/trunk/aprcl).

#standard primality checking by trial div, used to speed up processes within algorithm
def trialdiv(n):
    if n<2: return False
    elif n==2 or n==3: return True
    elif n%2==0: return False
    else:
        i = 3
        while i*i <= n:
            if n%i == 0: return False
            i += 2
    return True

# Jacobi sum class
class JacobiSum(object):
    def __init__(me, p, k, q):
        me.p = p
        me.k = k
        me.q = q
        me.pk = p**k
        me.m = (p-1)*p**(k-1)
        me.coefficient = [0]*me.m

    def one(me):
        me.coefficient[0] = 1
        for i in range(1,me.m):
            me.coefficient[i] = 0
        return me


    # product of JacobiSum, where "jac" refers to the jacobi sum. 
    def multiply(me, jac):
        jacobi_obj=JacobiSum(me.p, me.k, me.q)
        for i in range(me.m):
            for j in range(me.m):
                if (i+j)% me.pk < me.m:
                    jacobi_obj.coefficient[(i+j)% me.pk] += me.coefficient[i] * jac.coefficient[j]
                else:
                    subtract = (i+j) % me.pk - me.p ** (me.k-1)                    
                    while subtract>=0:
                        jacobi_obj.coefficient[subtract] -= me.coefficient[i] * jac.coefficient[j]
                        subtract-= me.p ** (me.k-1)

        return jacobi_obj


    def __mul__(me, right):
        if type(right) is int:
            jacobi_obj=JacobiSum(me.p, me.k, me.q)
            for i in range(me.m):
                jacobi_obj.coefficient[i] = me.coefficient[i]*right
            return jacobi_obj
        else: #JS case
            return me.multiply(right) 
    

    # power of JacobiSum
    def jacobiPowerFunc(me, power, mod):
        jacobi_obj=JacobiSum(me.p, me.k, me.q)
        jacobi_obj.coefficient[0]=1
        j_a = copy.deepcopy(me)
        while power>0:
            if power%2==1:
                jacobi_obj = (jacobi_obj * j_a).mod(mod)
            power //= 2
            j_a = j_a*j_a
            j_a.mod(mod)
        return jacobi_obj
    

    def mod(me, n):
        for i in range(me.m):
            me.coefficient[i] %= n
        return me
    

    def Sigma(me, x):
        jacobi_obj=JacobiSum(me.p, me.k, me.q)
        for i in range(me.m):
            if (i*x) % me.pk < me.m:
                jacobi_obj.coefficient[(i*x) % me.pk] += me.coefficient[i]
            else:
                subtract = (i*x) % me.pk - me.p ** (me.k-1)                    
                while subtract>=0:
                    jacobi_obj.coefficient[r] -= me.coefficient[i]
                    subtract=subtract- me.p ** (me.k-1)
        return jacobi_obj
    

    def OneOverSigma(me, x):
        jacobi_obj=JacobiSum(me.p, me.k, me.q)
        for i in range(me.pk):
            if i<me.m:
                if (i*x)%me.pk < me.m:
                    jacobi_obj.coefficient[i] += me.coefficient[(i*x)%me.pk]
            else:
                subtract = i - me.p ** (me.k-1)
                while subtract>=0:
                    if (i*x)%me.pk < me.m:
                        jacobi_obj.coefficient[subtract] -= me.coefficient[(i*x)%me.pk]
                    subtract=subtract- me.p ** (me.k-1)

        return jacobi_obj
    

    # If self/me is p^kth root of unity (mod N), return power where self/me is zeta^power
    def isRoU(me, N):
        m = me.m
        p = me.p
        k = me.k

        # power < m
        one = 0
        for i in range(m):
            if me.coefficient[i]==1:
                one += 1
                power = i
            elif me.coefficient[i] == 0:
                continue
            elif (me.coefficient[i] - (-1)) %N != 0:
                return False, None
        if one == 1:
            return True, power

        # power >= m
        for i in range(m):
            if me.coefficient[i]!=0:
                break
        r = i % (p**(k-1))
        for i in range(m):
            if i % (p**(k-1)) == r:
                if (me.coefficient[i] - (-1))%N != 0:
                    return False, None
            else:
                if me.coefficient[i] !=0:
                    return False, None

        return True, (p-1)*p**(k-1)+ r
            


#Timer Methods

    
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
        state = APRtest(num)
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

#Remaining Algo methods


# v_q(t)
def v(q, t):
    answer = 0
    while(t % q == 0):
        answer +=1
        t//= q #floor div
    return answer


def prime_factorize(n):
    listfac = []
    p=2
    while p*p <= n: 
        if n%p==0:
            num = 0
            while n%p==0:
                num+=1
                n//= p
            listfac.append((p,num))
        p+= 1

    if n!=1:
        listfac.append((n,1))

    return listfac


# calc e(t)
def e(t):
    s = 1
    q_list = []
    for q in range(2, t+2):
        if t % (q-1) == 0 and trialdiv(q):
            s *= q ** (1+v(q,t))
            q_list.append(q)
    return 2*s, q_list



#Primitive root finding
def spr(q):
    for r in range(2, q):
        lisst = set({})
        counter = 1
        for i in range(1, q):
            counter = (counter*r) % q
            lisst.add(counter)
        if len(lisst) == q-1:
            return r
    return None


#calc f_q(x)
def calc_f(q):
    g = spr(q)
    m = {}
    for x in range(1,q-1):
        m[pow(g,x,q)] = x
    f = {}
    for x in range(1,q-1):
        f[x] = m[ (1-pow(g,x,q))%q ]

    return f


#sum zeta^(ax+b[f(x)])
def calcJacobi_ab(p, k, q, a, b):
    jacobi_obj = JacobiSum(p,k,q)
    f = calc_f(q)
    for x in range(1,q-1):
        pk = p**k
        if (a*x+b*f[x]) % pk < jacobi_obj.m:
            jacobi_obj.coefficient[(a*x+b*f[x]) % pk] += 1
        else:
            r = (a*x+b*f[x]) % pk - p**(k-1)
            while r>=0:
                jacobi_obj.coefficient[r] -= 1
                r-= p**(k-1)
    return jacobi_obj


def calcJacobi(p, k, q):
    return calcJacobi_ab(p, k, q, 1, 1)

def calcJacobi3(p, k, q):
    j2q = calcJacobi(p, k, q)
    j21 = calcJacobi_ab(p, k, q, 2, 1)
    jacobi_obj = j2q * j21
    return jacobi_obj

def calcJacobi2(p, k, q):
    j31 = calcJacobi_ab(2, 3, q, 3, 1)
    j_conv = JacobiSum(p, k, q)
    for i in range(j31.m):
        j_conv.coefficient[i*(p**k)//8] = j31.coefficient[i]
    jacobi_obj = j_conv * j_conv
    return jacobi_obj

def APR_4_1(p, k, q, N):
    
    J = calcJacobi(p, k, q)
    s1 = JacobiSum(p,k,q).one()
    
    for x in range(p**k):
        if x % p == 0:
            continue
        t = J.OneOverSigma(x)
        t = t.jacobiPowerFunc(x, N)
        s1 = s1 * t
        s1.mod(N)

   
    r = N % (p**k)
    s2 = s1.jacobiPowerFunc(N//(p**k), N)
    J_alpha = JacobiSum(p,k,q).one()
    
    for x in range(p**k):
        if x % p == 0:
            continue
        t = J.OneOverSigma(x)
        t = t.jacobiPowerFunc((r*x)//(p**k), N)
        J_alpha = J_alpha * t
        J_alpha.mod(N)

    S = (s2 * J_alpha).mod(N)
    exist, h = S.isRoU(N)

    if not exist:
        # composite
        return False, None
    else:
        # possibly prime
        if h%p!=0:
            hold = 1
        else:
            hold = 0
        return True, hold



def APR_4_2(p, k, q, N):

    J = calcJacobi3(p, k, q)
    s1 = JacobiSum(p,k,q).one()
    for x in range(p**k):
        if x % 8 not in [1,3]:
            continue
        t = J.OneOverSigma(x)
        t = t.jacobiPowerFunc(x, N)
        s1 = s1 * t
        s1.mod(N)

    r = N % (p**k)
    s2 = s1.jacobiPowerFunc(N//(p**k), N)
    J_alpha = JacobiSum(p,k,q).one()
    
    for x in range(p**k):
        if x % 8 not in [1,3]:
            continue
        t = J.OneOverSigma(x)
        t = t.jacobiPowerFunc((r*x)//(p**k), N)
        J_alpha = J_alpha * t
        J_alpha.mod(N)

    if N%8 in [1,3]:
        S = (s2 * J_alpha ).mod(N)
    else:
        J2_delta = calcJacobi2(p,k,q)
        S = (s2 * J_alpha * J2_delta).mod(N)

    exist, h = S.isRoU(N)

    if not exist:
        # composite 
        return False, None
    else:
        # possible prime
        if h%p!=0 and (pow(q,(N-1)//2,N) + 1)%N==0:
            hold = 1
        else:
            hold = 0
        return True, hold

def APR_4_3(p, k, q, N):
    J2q = calcJacobi(p, k, q)
    s1 = (J2q * J2q * q).mod(N)
    s2 = s1.jacobiPowerFunc(N//4, N)

    if N%4 == 1:
        S = s2
    elif N%4 == 3:
        S = (s2 * J2q * J2q).mod(N)

    exist, h = S.isRoU(N)

    if not exist:
        # composite
        return False, None
    else:
        # possibly prime
        if h%p!=0 and (pow(q,(N-1)//2,N) + 1)%N==0:
            hold = 1
        else:
            hold = 0
        return True, hold

def APR_4_4(p, k, q, N):

    S2q = pow(-q, (N-1)//2, N)
    if (S2q-1)%N != 0 and (S2q+1)%N != 0:
        # composite
        return False, None
    else:
        # possible prime
        if (S2q + 1)%N == 0 and (N-1)%4==0:
            hold=1
        else:
            hold=0
        return True, hold


# Step 4
def APRtest_step4(p, k, q, N):

    if p>=3:
        result, hold = APR_4_1(p, k, q, N)
    elif p==2 and k>=3:
        result, hold = APR_4_2(p, k, q, N)
    elif p==2 and k==2:
        result, hold = APR_4_3(p, k, q, N)
    elif p==2 and k==1:
        result, hold = APR_4_4(p, k, q, N)

    return result, hold


def APRtest(N):
    t_list = [2,12,60,180,840,1260,1680,2520,5040,15120,55440,110880,720720,1441440,4324320,24504480,73513440]

    if N<=3:
        return False

    for t in t_list:
        et, q_list = e(t)
        if N < et*et:
            break
    else:
        return False

    g = gcd(t*et, N)
    if g > 1:
        return False

    # Step 2
    
    l = {}
    fac_t = prime_factorize(t)
    for p, k in fac_t:
        if p>=3 and pow(N, p-1, p*p)!=1:
            l[p] = 1
        else:
            l[p] = 0

    # Step 3, Step 4
    
    for q in q_list:
        if q == 2:
            continue
        fac = prime_factorize(q-1)
        for p,k in fac:
            # Step 4
            result, hold = APRtest_step4(p, k, q, N)
            if not result:
                # composite
                return False
            elif hold==1:
                l[p] = 1

    # Step 5          
    for p, value in l.items():
        if value==0:
            counter = 0
            i = 1
            found = False
            while counter < 30:
#found in online https://github.com/wacchoz's implementation, not exactly sure why 30
#is the number of times for the loop to be run still. Runtime works fine, though.
                q = p*i+1
                
                if N%q != 0 and trialdiv(q) and (q not in q_list):
                    counter += 1
                    k = v(p, q-1)
                    # Step 4
                    
                    result, hold = APRtest_step4(p, k, q, N)

                    if not result:
                        return False
                    elif hold == 1:
                        found = True
                        break
                i += 1
            if not found:
                return False

    # Step 6
    r = 1
    for t in range(t-1):
        r = (r*N) % et
        if r!=1 and r!= N and N % r == 0:
            return False
    return True


if __name__ == '__main__':

    print(APRtest(13))
    print(APRtest(876665514007))

