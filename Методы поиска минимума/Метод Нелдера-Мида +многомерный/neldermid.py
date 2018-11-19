# from matplotlib import pyplot as plt
import numpy as np
from numpy import sqrt,sum,hstack,ones,zeros,array,fill_diagonal,argsort,sin,cos
import matplotlib.animation as animation

# Вычислим минимум функции Розенброка. Известно, что он в точке [1,1]
# Вычислим минимум функции sin(x)*cos(y). известно, что один из минимумов в точке [-pi/2,0]
# Хитрое объявление, чтобы работать и с одним вектором [x1, y1], и набором векторов [[x1,y1],[x2,y2]]
def f(x):
	x=array(x).T
	return (sin(x[0])*cos(x[1])).T
	# return (100*(x[1]-x[0]**2)**2+(1-x[0])**2).T

# Размерность пространства, расстояние между вершинами начального симплекса, условие сходимости
n,t,eps=2,1,0.000001
# Параметры отражения, сжатия, растяжения
alpha,beta,gamma=1,0.5,2

# Построение регулярного симплекса
d1=(sqrt(n+1)+n-1)*t/sqrt(2)/n
d2=(d1-1/sqrt(2))*t
nn=ones([n,n])*d2
fill_diagonal(nn, d1)
simplex=hstack([zeros([n,1]),nn]).T

# Сортируем вершины симплекса, первая - наибольшая f
simplex = simplex[argsort(-f(simplex))]
# По всем точкам, кроме xh ищем центр тяжести
xc=1/n*sum(array([s for s in simplex[1:]]),axis=0)

simplex1=simplex2=simplex3=[]

N=0
while 1/(n+1)*sum(array([abs(f(s)-f(xc)) for s in simplex]))>eps: # Проверка сходимости
	N+=1
	if N>1e3:
		break
	# Выбираем три опорные точки
	xh,xg,xl=simplex[0],simplex[1],simplex[-1]

	# По всем точкам, кроме xh ищем центр тяжести
	xc=1/n*sum(array([s for s in simplex[1:]]),axis=0)

	# Генерируем отраженную точку
	xr=(1+alpha)*xc-alpha*xh

	Shrink=True
	if f(xr)<f(xl):
		xe=(1-gamma)*xc+gamma*xr
		if f(xe)<f(xl):
			xh=xe
			Shrink=False
		elif f(xe)>f(xl):
			xh=xr
			Shrink=False
	elif (f(xl)<f(xr))&(f(xr)<f(xg)):
		xh=xr
		Shrink=False
	elif (f(xh)>f(xr))&(f(xr)>f(xg)):
		xr,xh=xh,xr
		Shrink=True
	elif (f(xr)>f(xh)):
		Shrink=True

	if Shrink:
		xs=beta*xh+(1-beta)*xc
		if f(xs)<f(xh):
			xh=xs
		elif f(xs)>f(xh):
			# Обновляем симплекс и сжимаем его
			simplex[0],simplex[1],simplex[-1]=xh,xg,xl
			# Сжимаем симплекс, кроме xl
			simplex[:-1]=xl+(simplex[:-1]-xl)/2
			xh,xg,xl=simplex[0],simplex[1],simplex[-1]

	# Обновляем симплекс
	simplex[0],simplex[1],simplex[-1]=xh,xg,xl
	# Сортируем симплекс
	simplex = simplex[argsort(-f(simplex))]
	simplex1.append(simplex)

print('Сошлось за ',N,' шагов, конечная точка ',1/(n+1)*sum(array([s for s in simplex]),axis=0))
print(array(simplex1))
# print(arr)