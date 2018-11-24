import numpy as np
from numpy import cos,sin,pi,log
from matplotlib import pyplot as plt
from numpy.linalg import inv,norm,det,solve
from scipy.optimize import minimize

n=3
A=np.array([[1,2,2],[4,5,6],[7,8,9]])
b=np.array([2,3,4])
E=np.identity(n)

C=det(A)

B=A.T@A
v=A.T@b

# (B.T@B+alpha*E)z=B.T@v

# Псевдорешение
def z(alpha):
	return inv(B.T@B+alpha*E)@(B.T@v)

# Невязка псевдорешения
def r(alpha):
	return norm(A@z(alpha)-b)

Alpha=minimize(r,0.5,method='Nelder-Mead').x

# A=np.linspace()

# ответ 
x=z(Alpha)
print(x)

print(A@x-b)
print(r(Alpha))
if det(A)!=0:
	print(solve(A,b))