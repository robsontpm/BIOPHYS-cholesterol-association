"""
This module grants access to simple models located in models folder. Currently only sphere and cube are supported. 
"""

import sys,os
from panda3d.core import Filename
 
import os
mydir = os.path.dirname(__file__)
mydir = Filename.fromOsSpecific(mydir).getFullpath()

def sphere(loader, x, y, z, r):
	"""
	This loads a sphere object and return it as an object ready to be attached to renderer by add() function.
	samle code::

		sphere = assets.sphere(renderer, 0, 0, 0, 10.0)
		renderer.add(render.primitive(sphere))
	"""
	model = loader.loadModel(mydir + "/models/sphere")
	model.setScale(2 * r)
	model.setPos(x, y, z)
	model.setColor(1.0, 1.0, 1.0, 1.0)
	return model
	
def cube(loader, x, y, z, r):
	"""
	This loads a cube object and return it as an object ready to be attached to renderer by add() function.
	samle code::
	
		cube = assets.cube(renderer, 0, 0, 0, 10.0)
		renderer.add(render.primitive(cube))
	"""
	model = loader.loadModel(mydir + "/models/cube")
	model.setScale(r)
	model.setPos(x, y, z)
	model.setColor(1.0, 1.0, 1.0, 1.0)
	return model
