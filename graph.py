#!/usr/bin/python

import random
import math
import copy
from OpenGL.GLUT import *
from OpenGL.GLU  import *
from OpenGL.GL   import *
import sys
from numpy import *


class Elem:
	def __init__(self,x,y,label):
		self.x = x
		self.y = y
		self.neighbours = set()
		self.label = label
	def dist(self, e):
		return math.sqrt((e.x-self.x)**2 + (e.y-self.y)**2)
	def cmp(self, a, b):
		if(self.dist(a) > self.dist(b)):
			return 1
		else:
			return -1


class Graph:
	def __init__(self):
		pass
	def __init__(self,elems, method="gaussian", sigma=1.0, k=5):
		self.elems = elems
		self.N = len(elems)
		self.W = zeros((self.N,self.N))
		self.subset = set()
		if(method == "KNN"):
			print "Generating KNN ..."
			for elem in elems:
				dist_list = sorted(elems,elem.cmp)
				for i in range(k):
					elem.neighbours.add(dist_list[i])

			print "Generating neighbour weight matrix ..."
			for i in range(self.N):
				for j in range(self.N):
					if elems[i] in elems[j].neighbours:
						self.W[i,j] = 1.0;
		elif(method == "gaussian"):
			print "Generating gaussian weight matrix ..."
			for i in range(self.N):
				for j in range(self.N):
					self.W[i,j] = exp( -(elems[i].dist(elems[j])**2) / (2*sigma**2) )
			print self.W
			print self.W.max()
		self.Y = zeros(self.N)
		for i in range(self.N):
			self.Y[i] = elems[i].label
		self.Yt = self.Y.view()
	def subgraph(self,subset):
		sg = Graph()
		elems = []
		for i in range(self.N):
			if i in subset:
				elems.append(self.elems[i])

	def algorithm_11_1(self,steps):
		D = zeros(self.W.shape)
		for i in arange(self.N):
			D[i,i] = self.W[i,:].sum()
		Di = linalg.pinv(D)
		for i in range(steps):
			self.Yt = dot(dot(Di,self.W),self.Yt)
			for i in range(self.N):
				if(self.Y[i] != 0.0):
					self.Yt[i] = self.Y[i]
	
	def select_subset(self,m):
		subset = set()
		#1 select max m elements from labelled data
		for i in range(self.N):
			if self.Yt[i] != 0.0 and len(subset) < m:
				subset.add(i)
		#2 select elements that are the furthest form the subset
		r = m - len(subset)
		for n in range(r):
			def subset_dist(i):
				mem = {}
				if i in mem:
					return mem[i]
				dist = 0.0
				for j in range(self.N):
					if j in subset:
						dist += self.W[i,j]
				mem[i] = dist
				return dist
			def cmp(a,b):
				if subset_dist(a) > subset_dist(b):
					return 1
				else:
					return -1
			candidates = []
			for i in range(self.N):
				if i not in subset:
					candidates.append(i)
			subset.add(sorted(candidates,cmp)[0])
			#todo outliers
		self.subset = subset

	def draw_graph(self):
		glLineWidth(0.75)
		for i in range(self.N):	
			for j in range(self.N):
				wij = self.W[i,j]
				if(wij > 0.0):
					wij *= 0.5
					if(self.Yt[i] > 0.0 and self.Yt[j] > 0.0):
						glColor3f(1.0,1.0-wij,1.0-wij)	
					elif(self.Yt[i] < 0.0 and self.Yt[j] < 0.0):
						glColor3f(1.0-wij,1.0-wij,1.0)
					else:
						glColor3f(1.0-wij,1.0-wij,1.0-wij)
					glBegin(GL_LINE)
					glVertex3f(self.elems[i].x,self.elems[i].y,wij)
					glVertex3f(self.elems[j].x,self.elems[j].y,wij)
					glEnd()
		glPointSize(5.0)
		for i in range(self.N):
			if i in self.subset:
				if(self.Yt[i] > 0.0):
					glColor3f(1.0,0.0,0.0)
				elif(self.Yt[i] < 0.0):
					glColor3f(0.0,0.0,1.0)
				else:
					glColor3f(0.0,0.0,0.0)
			else:
				if(self.Yt[i] > 0.0):
					glColor3f(1.0,0.5,0.5)
				elif(self.Yt[i] < 0.0):
					glColor3f(0.5,0.5,1.0)
				else:
					glColor3f(0.5,0.5,0.5)
			
			
			glBegin(GL_POINT)
			glVertex3f(self.elems[i].x,self.elems[i].y,1.1)
			glEnd()

	def draw_fun(self):
		glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT)
		glPushMatrix()
		self.draw_graph()
		glPopMatrix()
		glutSwapBuffers()
		

def draw_init(disp_fun):
	glutInit(sys.argv)
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH)
	glutInitWindowSize(400,400)
	glutCreateWindow("graph")
	glClearColor(1.0, 1.0, 1.0, 1.0)

	glutDisplayFunc(disp_fun)
	glMatrixMode(GL_PROJECTION)
	glOrtho(-10.0,110.0,-10.0,110.0,-10,10)
	glMatrixMode(GL_MODELVIEW)
	glLoadIdentity()
	glEnable(GL_DEPTH_TEST)
	glEnable(GL_POINT_SMOOTH)
	glEnable(GL_LINE_SMOOTH)
	glShadeModel(GL_SMOOTH)
	glEnable(GL_BLEND)
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
	#gluLookAt(0,0,10,
	#		  0,0,0,
	#		  0,1,0)
	#glPushMatrix()

def gen_uniform_circle_elems(px,py,r,N,Pl,klass):
	elems = []
	for i in range(N):
		x = px + (random.random()-0.5)*2.0*r
		y = py + (random.random()-0.5)*2.0*r
		while( (px-x)**2 + (py-y)**2 > r**2 ):
			x = px + (random.random()-0.5)*2.0*r
			y = py + (random.random()-0.5)*2.0*r
		if(random.random() < Pl):
			e = Elem(x,y,float(klass))
		else:
			e = Elem(x,y,0.0)
		elems.append(e)
	return elems

def main():
	S = 100.0
	N = 200
	Pl = 0.1
	elems = []
	print "Generating elems ..."
	#elems.extend(gen_uniform_circle_elems(25,50,25,100,0.5,1.0))
	#elems.extend(gen_uniform_circle_elems(75,50,25,100,0.5,-1.0))
	
	
	for i in range(N):
		x = random.random()*S
		y = random.random()*S
		if(random.random() < Pl):
			if(x < S/2):
				e = Elem(x,y,1.0)
			else:
				e = Elem(x,y,-1.0)
		else:
			e = Elem(x,y,0.0)
		elems.append(e)
	print "Generating graph ..."
	graph = Graph(elems,"gaussian",sigma=10)	
	#graph = Graph(elems,"KNN",k=20)	
	print "Subset Selection"
	graph.select_subset(60)
	#print "Class propagation ..."
	#graph.algorithm_11_1(100)
	print "Draw graph ... "
	draw_init(graph.draw_fun)
	glutMainLoop()


if __name__ == "__main__":
	main()
