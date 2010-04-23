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
	def __init__(self,elems = None, method="none", sigma=1.0, k=5, W = None):
		self.N  = 0
		self.W  = None
		self.Y  = None
		self.Yt = None
		self.elems = None

		self.subset       = set()
		self.graph_to_sub = {}
		self.sub_to_graph = {}
		self.sub_graph    = None

		if elems != None:
			self.N     = len(elems)
			self.elems = elems
			self.Y     = zeros(self.N)

			if method == "KNN" :
				self.W = self.KNN_W(elems,k)
			elif method == "gaussian" :
				self.W = self.gaussian_W(elems,sigma)
			elif W != None:
				self.W = W

			for i in range(self.N):
				self.Y[i] = elems[i].label
			self.Yt = self.Y.view()
	def KNN_W(self,elems,k):
		N = len(elems)
		W = zeros((N,N))
		print "Generating KNN ..."
		for elem in elems:
			dist_list = sorted(elems,elem.cmp)
			for i in range(k):
				elem.neighbours.add(dist_list[i])

		print "Generating neighbour weight matrix ..."
		for i in range(N):
			for j in range(N):
				if elems[i] in elems[j].neighbours:
					W[i,j] = 1.0
		return W
	def gaussian_W(self,elems,sigma):
		N = len(elems)
		W = zeros((N,N))
		print "Generating gaussian weight matrix ..."
		for i in range(N):
			for j in range(N):
				W[i,j] = exp( -(elems[i].dist(elems[j])**2) / (2*sigma**2) )
		return W
	def reorder(self,subset):
		S = len(subset)
		N = self.N
		order = []
		elems2 = range(N)
		W2 = zeros((N,N))
		subset2 = set(range(S))
		for i in range(N):
			if i in subset:
				order.append(i)
		for i in range(N):
			if i not in subset:
				order.append(i)
		for i in range(N):
			elems2[i] = self.elems[order[i]]
		for i in range(N):
			for j in range(N):
				W2[i,j] = self.W[order[i],order[j]]

		#test
		#for i in range(N):
		#	if i in subset:
		#		print self.elems[i].x
		#print "------"
		#for i in range(N):
		#	if i in subset2:
		#		print elems2[i].x

		G2 = Graph(elems=elems2,W=W2);
		G2.subset = subset2
		G2.order = order
		self.reordered  = G2
		return G2
	def algorithm_18_22(self, mu = 0.5, epsilon = 0.001):
		S = len(self.subset)
		N = self.N
		Sss = zeros((S,S))
		Wss = zeros((S,S))
		Wrs = zeros((N-S,S))
		Wsr = zeros((S,N-S))
		Dss = zeros((S,S))
		Id  = zeros((S,S))
		for i in range(S):
			if(self.Yt[i] < -0.99 or self.Yt[i] > 0.99):
				Sss[i,i] = 1
		for i in range(S):
			for j in range(S):
				Wss[i,j] = self.W[i,j]
		for i in range(S):
			for j in range(N-S):
				Wsr[i,j] = self.W[i,S+j]
		
		for I in range(N-S):
			i = I+S
			for j in range(S):
				sum = epsilon
				for k in range(S):
					sum += self.W[i,k]
				Wrs[I,j] = self.W[i,j] / sum
				
		for i in range(S):
			sum = 0.0
			for j in range(S):
				sum += self.W[i,j]
			D[i,i] = sum
		for i in range(S):
			Id[i,i] = 1
		
		A = (Sss + mu*(Dss - Wss - dot(Wsr,Wrs) + epsilon*Id))
		B = dot(Sss,self.Yt)
		Yt = linalg.solve(A,B)

	def subgraph(self,subset):
		M = len(subset)
		sub_elems = []
		sub_W = zeros((M,M))
		graph_to_sub = {}
		sub_to_graph = {}

		for i in range(self.N):
			if i in self.subset:
				sub_elems.append(self.elems[i])
				graph_to_sub[i] = len(sub_elems)-1
				sub_to_graph[len(sub_elems)-1] = i
		
		for i in range(M):
			for j in range(M):
				k = sub_to_graph[i]
				l = sub_to_graph[j]
				sub_W[i,j] = self.W[k,l]

		sg = Graph(elems = sub_elems, W = sub_W)
		sg.graph_to_sub  = graph_to_sub
		sg.sub_to_graph  = sub_to_graph
		sg.subset = set(range(M))
		self.sub_graph   = sg

		return sg
	def from_subgraph(self,sub_graph):
		for i in range(sub_graph.N):
			self.Yt[sub_graph.sub_to_graph[i]] = sub_graph.Yt[i]
	def estimate(self,i,subset):
		N = 0.0
		D = 0.00001
		for j in range(self.N):
			if j in subset:
				N += self.W[i,j]*self.Yt[j]
				D += self.W[i,j]
		return N/D
	def algorithm_11_1(self,steps = -1):
		if steps < 0:
			steps = self.N

		D = zeros(self.W.shape)
		for i in arange(self.N):
			D[i,i] = self.W[i,:].sum()
		Di = linalg.pinv(D)
		for i in range(steps):
			self.Yt = dot(dot(Di,self.W),self.Yt)
			for i in range(self.N):
				if(self.Y[i] != 0.0):
					self.Yt[i] = self.Y[i]
	def algorithm_18_2(self, mu = 0.5, epsilon = 0.001):
		N = self.N
		S = zeros((N,N))
		D = zeros((N,N))
		Id = zeros((N,N))
		for i in range(N):
			if(self.Yt[i] < -0.99 or self.Yt[i] > 0.99):
				S[i,i] = 1
		for i in range(N):
			sum = 0.0
			for j in range(N):
				sum += self.W[i,j]
			D[i,i] = sum
		for i in range(N):
			Id[i,i] = 1

		L = D - self.W
		
		A = S + mu*L + mu*epsilon*Id
		B = dot(S,self.Yt)
		self.Yt = linalg.solve(A,B)
	def select_subset(self,m,surf=0.5,sigma = 0.001):
		def furthest_from(subset, candidates):
			def subset_dist(subset,i):
				dist = 0.0
				for j in range(self.N):
					if j in subset:
						dist += self.W[i,j]
				return dist
			mindist = subset_dist(subset,list(candidates)[0])
			minc = list(candidates)[0]
			for c in candidates:
				newdist = subset_dist(subset,c)
				if newdist < mindist:
					mindist = newdist
					minc = c
			return minc
			
		subset = set()
		self.subset = subset

		print "step 1 select max m(",m,") elements from labelled data ... "
		for i in range(self.N):
			if self.Yt[i] != 0.0 and len(subset) < m:
				subset.add(i)
		print "... picked:",len(subset), " elements"

		print "step 2 select elements that are the furthest from the subset ... "
		r = m - len(subset)
		rest = set(range(self.N)) - subset
		for n in range(r):
			c = furthest_from(subset,rest)
			rest.remove(c)
			subset.add(c)
		print "...subset length :",len(subset)

		print "step 3 generate subgraph from subset ..."
		subgraph = self.subgraph(subset)

		print "step 4 estimate classes on subset ..."
		subgraph.algorithm_11_1()
		self.from_subgraph(subgraph)

		print "step 5 propagate subset estimation to R"
		for i in range(self.N):
			if i not in subset:
				self.Yt[i] = self.estimate(i,subset)

		print "step 6 pick least confident..."
		X = int(surf*len(subset))
		Sh = list(subset)
		Rl = []
		for i in range(self.N):
			if i not in subset:
				Rl.append(i)
		def confidence_cmp(i,j):
			if(abs(self.Yt[i]) > abs(self.Yt[j])):
				return 1
			else:
				return -1
		Sh.sort(confidence_cmp)
		Rl.sort(confidence_cmp)
		Sh = Sh[-X:]
		Rl = set(Rl[:X])
		for i in Sh:
			min = 10.0
			for j in range(self.N):
				s = 0.0
				for k in range(self.N):
					if k in subset:
						s += self.W[j,k]
				if s < min:
					min = s
			if min > sigma:
				subset.remove(i)
				c = furthest_from(subset,Rl)
				subset.add(c)
				Rl.remove(c)

		print "Generating new subgraph..."
		subgraph = self.subgraph(subset)
		
		return subset
	def draw_graph_confidence(self):
		glLineWidth(0.75)
		for i in range(self.N):	
			for j in range(self.N):
				wij = self.W[i,j]
				if(wij > 0.2):
					alpha = 1.0-min(abs(self.Yt[i]),abs(self.Yt[j]))
					glColor3f(1.0-(wij*alpha),1.0-wij,1.0-(wij*(1.0-alpha)))
					#wij *= 0.5
					#if(self.Yt[i] > 0.0 and self.Yt[j] > 0.0):
					#	glColor3f(1.0,1.0-wij,1.0-wij)	
					#elif(self.Yt[i] < 0.0 and self.Yt[j] < 0.0):
					#	glColor3f(1.0-wij,1.0-wij,1.0)
					#else:
					#	glColor3f(1.0-wij,1.0-wij,1.0-wij)
					glBegin(GL_LINE)
					glVertex3f(self.elems[i].x,self.elems[i].y,wij)
					glVertex3f(self.elems[j].x,self.elems[j].y,wij)
					glEnd()
		for i in range(self.N):
			glPointSize(5.0)
			alpha = float(abs(self.Yt[i]))
			glColor3f(alpha,0.0,1.0-alpha)
			if i in self.subset:
				glColor3f(0,0,0)
			else:
				continue
			#if i in self.Sh:
			#	glColor3f(0,0,0)
			#if i in self.Rl:
			#	glColor3f(1,0,1)
		#	if i in self.subset:
		#		glPointSize(5.0)
		#		if(self.Yt[i] > 0.0):
		#			glColor3f(1.0,0.2,0.0)
		#		elif(self.Yt[i] < 0.0):
		#			glColor3f(0.0,0.2,1.0)
		#		else:
		#			glColor3f(0.0,1.0,0.0)
		#	else:
		#		glPointSize(0.5)
		#		if(self.Yt[i] > 0.0):
		#			glColor3f(1.0,0.5,0.5)
		#		elif(self.Yt[i] < 0.0):
		#			glColor3f(0.5,0.5,1.0)
		#		else:
		#			glColor3f(0.5,0.5,0.5)
			
			
			glBegin(GL_POINT)
			glVertex3f(self.elems[i].x,self.elems[i].y,1.1)
			glEnd()
	def draw_graph_color(self):
		glLineWidth(0.75)
		for i in range(self.N):	
			for j in range(self.N):
				wij = self.W[i,j]
				if(wij < 0.2):
					continue
				#wij *= 0.5
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
		for i in range(self.N):
			glPointSize(5.0)
			#if i in self.subset:
			glPointSize(5.0)
			if(self.Yt[i] > 0.0):
				glColor3f(1.0,0.2,0.0)
			elif(self.Yt[i] < 0.0):
				glColor3f(0.0,0.2,1.0)
			else:
				glColor3f(0.0,1.0,0.0)
			#else:
			#	glPointSize(0.5)
			#	if(self.Yt[i] > 0.0):
			#		glColor3f(1.0,0.5,0.5)
			#	elif(self.Yt[i] < 0.0):
			#		glColor3f(0.5,0.5,1.0)
			#	else:
			#		glColor3f(0.5,0.5,0.5)
			
			
			glBegin(GL_POINT)
			glVertex3f(self.elems[i].x,self.elems[i].y,1.1)
			glEnd()
	def draw_fun(self):
		glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT)
		glPushMatrix()
		self.draw_graph_color()
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
def gen_uniform_rect(px,py,sx,sy,N,Pl,klass):
	elems = []
	labelled = 0
	for i in range(int(N)):
		x = px + random.random() * sx
		y = py + random.random() * sy
		if(random.random() < Pl):
			e = Elem(x,y,float(klass))
			labelled += 1
		else:
			e = Elem(x,y,0.0)
		elems.append(e)
	print "Labelled :",labelled
	return elems

def main():
	random.seed(0)
	N = 400
	Pl = 0.3
	elems = []
	print "Generating elems ..."
	#elems.extend(gen_uniform_circle_elems(25,50,25,N/2,Pl,1.0))
	#elems.extend(gen_uniform_circle_elems(75,50,25,N/2,Pl,-1.0))
	elems.extend(gen_uniform_rect(0,0,50,100,N/2,Pl,1.0))
	elems.extend(gen_uniform_rect(50,0,50,100,N/2,Pl,-1.0))
	
	print "Generating graph ..."
	graph = Graph(elems,"gaussian",sigma=10)	
	#graph = Graph(elems,"KNN",k=20)	
	print "Subset Selection ..."
	graph.select_subset(20)
	#graph.subgraph(graph.subset)
	#print "Class propagation ..."
	#graph.algorithm_11_1(100)
	#graph.algorithm_18_2()
	print "Reordering ..."
	graph.reorder(graph.subset)
	print "Draw graph ... "
	draw_init(graph.draw_fun)
	glutMainLoop()


if __name__ == "__main__":
	main()
