"""
This module uses Panda3D to allow for simple rendering in 3D. Panda3D may be downloaded and installed for free from http://www.panda3d.org/.

Basic usage consist of creating Render object and then add some objects to it. 
Finally, one executes run() function to display graphics. Output window
can be controlled with mouse in a standard way: scroll or move mouse with both buttons pressed to zoom, hold left button 
and move mouse to move up and down, hold right button and move to rotate the view around the origin. 

:warning: For now, only one invocation of run() per application is possible. We are working to allow for multiple windows. 

:warning: Interface of this submodule may slightly evolve during development. 

The simplest code to render a half-transparent cube with side = 10au at the origin looks like this::
	
	from pydmga import draw
	from pydmga.draw import render
	from pydmga.draw import assets

	r = draw.Render(render.Color.WHITE_OPAQUE)
	r.add(draw.Asset(assets.cube, (0.0, 0.0, 0.0), 5.0, render.Color.BLUE_QUARTER))
	r.run()
"""

import render
from render import *
import assets
import gen
import geometry

__all__ = ["geometry", "gen", "assets", "render"]

