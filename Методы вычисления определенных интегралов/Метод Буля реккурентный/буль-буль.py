from numpy import sqrt,sum,abs
def f(x):
    return x

def T(N,f,a,b):
    if N==0:
        return (b-a)/2*(f(a)+f(b))
    else:
        g=2**N
        M=int(g/2)
        h=(b-a)/g
        def xk(k):
            return a+k*h
        return T(N-1,f,a,b)/2+h*sum([f(xk(2*k+1)) for k in range(1,M+1)])

def S(N,f,a,b):
    return (4*T(N,f,a,b)-T(N-1,f,a,b))/3

def B(N,f,a,b):
    return (16*S(N,f,a,b)-T(N-1,f,a,b))/15

def getIntegral(f,a,b,eps):
    N=2
    int1=B(1,f,a,b)
    int2=B(2,f,a,b)
    while abs(int2-int1)>eps:
        N+=1
        int1=int2
        int2=B(N,f,a,b)
    return [int2,N]

a=0
b=1
eps=0.01
print(getIntegral(f,a,b,eps))