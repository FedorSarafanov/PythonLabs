def f(x):
    return x

def boole(a,b):
    boole=(
        (b-a)/90*
            (    
             7*f(a)+32*f((3*a+b)/4)
             +12*f((a+b)/2)+
             32*f((a+3*b)/4)+
             7*f(b)    
             )
          )
    return boole

def integral(a,b,epsilon):
    n=1
    integral=0
    tmp=2
    while abs(integral-tmp)>epsilon:
        h=1/n
        print(h)
        tmp=integral
        integral=0
        for j in range(n):
            tmp1=boole(a+j*h,a+(j+1)*h)
            integral+=tmp1
        n*=2
    return integral

epsilon=0.01
# print(integral(2000,-0.5,epsilon))
print(boole(0,100))