"""
:warning: deprecated
"""

from direct.showbase.ShowBase import ShowBase
import random
from math import pi, sin, cos, sqrt
import math
from direct.task import Task

import os
import re
import sys

from geometry import * #voro.geometry
import assets #voro.graphics
 
class TrajectoryGenerator(ShowBase):
 
	def __init__(self, repository, count, box_rad = 5.0):
		ShowBase.__init__(self)
		
		self.repository = repository
		self.count = count

		# Load the environment model.
		#self.environ = self.loader.loadModel("models/environment")
		# Reparent the model to render.
		#self.environ.reparentTo(self.render)
		# Apply scale and position transforms on the model.
		#self.environ.setScale(0.25, 0.25, 0.25)
		#self.environ.setPos(-8, 42, 0)
		
		self.rad = box_rad
		self.count = count	

		self.balls = []
		random.seed()
		for i in range(count):			
			x = random.random() * 2 * self.rad - self.rad
			y = random.random() * 2 * self.rad - self.rad
			z = random.random() * 2 * self.rad - self.rad
			ball  = assets.sphere(self.loader, x, y, z, 1.0)
			ball.reparentTo(self.render) 			
			self.balls.append(ball)
			
		self.oobe()
		self.camera.setPos(0, -10, 5)
		self.camera.setHpr(0, 0, 0)
			
		self.taskMgr.add(self.simulation, "simulation")		
		
		self.accept('a', self.startOutput)
		self.accept('s', self.stopOutput)
		self.isOutputStarted = False
		self.isOutputStopped = False
		
	def startOutput(self):
		if (not self.isOutputStarted):
			self.isOutputStarted = True
			self.current_frame = 0			
			if not os.path.isdir(self.repository):
				os.mkdir(self.repository)
			fileObj = open(os.path.join(self.repository, self.repository + ".desc"), "w")
			fileObj.write(" ".join(map(str, [2*self.rad, 2*self.rad, 2*self.rad, 0, 0, 0, "\n"])))
			i = 0
			for ball in self.balls:
				i += 1
				x = ball.getX() + self.rad
				y = ball.getY() + self.rad
				z = ball.getZ() + self.rad
				fileObj.write(" ".join(map(str, [i, "A", "1", "LIP", x, y, z, 1.0])))
				fileObj.write("\n")
			fileObj.close()
			print ("Output started")
			self.doOutput(); #doOutput to ensure that 0-th frame has the same positions as .desc file
	
	def stopOutput(self):
		if (self.isOutputStarted and not self.isOutputStopped):
			self.isOutputStopped = True
			print ("Output stopped")
		
	def dist(self, s, o):
		xx = s.getX() - o.getX()
		yy = s.getY() - o.getY()
		zz = s.getZ() - o.getZ()
		
		return sqrt(xx*xx + yy*yy + zz*zz)
		
			
	# Define a procedure to move the camera.
	def simulation(self, task):
		angleDegrees = task.time * 6.0
		angleRadians = angleDegrees * (pi / 180.0)
		#self.camera.setPos(20 * sin(angleRadians), -20.0 * cos(angleRadians), 3)
		#self.camera.setHpr(angleDegrees, 0, 0)
		
		N = len(self.balls) 
		v_1_x = {}       
		v_1_y = {}
		v_1_z = {}
		
		v_2_x = {}       
		v_2_y = {}
		v_2_z = {}
		for ball in self.balls:
			v_1_x[ball] = 0
			v_1_y[ball] = 0
			v_1_z[ball] = 0
			v_2_x[ball] = 0
			v_2_y[ball] = 0
			v_2_z[ball] = 0
			for other in self.balls:				
				if (ball != other):	
					v_1_x[ball] += -ball.getX() + other.getX()
					v_1_y[ball] += -ball.getY() + other.getY()
					v_1_z[ball] += -ball.getZ() + other.getZ()
				
					d = self.dist(ball, other)
					if (d < self.rad*0.8):
						if d < 0.001:
							d = 0.001
						v_2_x[ball] += (ball.getX() - other.getX()) / d
						v_2_y[ball] += (ball.getY() - other.getY()) / d
						v_2_z[ball] += (ball.getZ() - other.getZ()) / d
						
			v_1_x[ball] /= (N-1)
			v_1_y[ball] /= (N-1)
			v_1_z[ball] /= (N-1)
			
			v_1_x[ball] = (v_1_x[ball] - ball.getX())
			v_1_y[ball] = (v_1_y[ball] - ball.getY())
			v_1_z[ball] = (v_1_z[ball] - ball.getZ())									
		
		for ball in self.balls:	
			st1 = 0.001
			st2 = 0.005
			stz = random.random() * 15.0
			
			angleRadians = angleDegrees * (pi / 180.0)
			v_zewn_x = -0.00001 * ball.getY() + 0.01 * (random.random() - 0.5)
			v_zewn_y =  0.00001 * ball.getX() + 0.01 * (random.random() - 0.5)
			v_zewn_z =  0.01 * (random.random() - 0.5)
			
			v_x = st1*v_1_x[ball] + st2*v_2_x[ball] + stz*v_zewn_x
			v_y = st1*v_1_y[ball] + st2*v_2_y[ball] + stz*v_zewn_y
			v_z = st1*v_1_z[ball] + st2*v_2_z[ball] + stz*v_zewn_z
			ball.setPos(ball, min(0.5, v_x), min(0.5, v_y), min(0.5, v_z))	
			if (ball.getX() <= -self.rad):
				ball.setX(-self.rad + 0.0001);
			if (ball.getY() <= -self.rad):
				ball.setY(-self.rad + 0.0001);
			if (ball.getZ() <= -self.rad):
				ball.setZ(-self.rad + 0.0001);
			if (ball.getX() >= self.rad):
				ball.setX(self.rad - 0.0001);
			if (ball.getY() >= self.rad):
				ball.setY(self.rad - 0.0001);
			if (ball.getZ() >= self.rad):
				ball.setZ(self.rad - 0.0001);
						
		if (self.isOutputStarted and not self.isOutputStopped):					
			self.doOutput();	
		return Task.cont	
		
	def doOutput(self):		
		f_coord = open(os.path.join(self.repository, self.repository + ('%08d' % self.current_frame) + '.coord'), "w")
		for ball in self.balls:
			x = ball.getX() + self.rad
			y = ball.getY() + self.rad
			z = ball.getZ() + self.rad
			r = 1.0
			f_coord.write(" ".join(map(str, [x, y, z, r])))
			f_coord.write("\n")		
		f_coord.close()		
		f_voro = open(os.path.join(self.repository, self.repository + ('%08d' % self.current_frame) + '.voro'), "w")
		f_voro.close()		
		self.current_frame += 1	
 
# Below is the code of 'main' function
if __name__ == "__main__":

	if len(sys.argv) < 4:
		print "Trajectory generator for DMGa project\nUsage: voro_gen set_name atom_count box_radius"
		sys.exit(0)
		
	app = TrajectoryGenerator(sys.argv[1], int(sys.argv[2]), float(sys.argv[3]))
	app.run()