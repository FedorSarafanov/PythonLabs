import scipy.integrate as spint
from numpy import cos,sin,exp, log as ln
import numpy as np
from numpy.linalg import inv,norm,cond,solve

# подинтегральная функция
def f(i,t):
    return exp(-i**2*t**2)/(1+i*t*cos(t)**2+t**2)

def Int(m):
	a=0
	b=1

	def F(t):
		return f(m,t)

	def f38(a,b):
		return (b-a)/8*(
				F(a)+
				3*F((2*a+b)/3)+
				3*F((a+2*b)/3)+
				F(b)
			)

	eps=0.01
	Iold=np.inf
	I=f38(a,b)
	n=2
	k=0

	while abs(I-Iold)/abs(I)>eps:
		k+=1
		n=n*2
		S=np.linspace(a,b,n) # эквидистантное разбиение
		h=S[1]-S[0] # шаг разбиения

		Iold=I
		I=np.sum([f38(ai,ai+h) for ai in S[:-1]])
		
	return [I,k]

n=5

# Заполнение матрицы A
A=np.zeros([n,n], dtype=np.float64)
for j in range(1,n-1):
	A[j][j]=-(3+sin(j+1)**2*cos(j+1)**5/(j+1+1))
for j in range(1,n-1):
	A[j][j-1]=1
	A[j][j+1]=(1+cos(j+1)**2)
A[0][0]=1
A[n-1][n-1]=1
A[n-1][n-2]=-0.9

# Заполнение матрицы b
b=np.zeros(n, dtype=np.float64)
for j in range(1,n-1):
	def F(t): 
		return f(j+1,t) 
	# print(j)
	if n<=1000:
		b[j]=-Int(j+1)[0]
	else:
		b[j]=-spint.quad(F,0,1)[0] 
b[n-1]=1
print('FiL matrix finished')


x9=solve(A,b)

d,c,a=np.diag(A,k=0),np.diag(A,k=1),np.diag(A,k=-1)
M,L,x=np.zeros(n),np.zeros(n),np.zeros(n)

# Прямая прогонка 
for i in range(0,n):
	if i==0:
		M[i]=-c[0]/d[0]
		L[i]=b[0]/d[0]
	else:
		M[i]=(b[i-1]-M[i-1]*a[i-2])/(a[i-2]*L[i-1]+d[i-1])
		L[i]=-c[i-1]/(a[i-2]*L[i-1]+d[i-1])

# Обратная прогонка 
x[n-1]=(b[n-1]-M[n-1]*a[n-2])/(a[n-2]*L[n-1]+d[n-1])

for i in reversed(range(0,n-1)):
	x[i]=x[i+1]*L[i+1]+M[i+1]

print(x)
print(x9)