"""
This is a Python interface to DMGAlpha library. It has more functions and it
is much simpler than the bare C++ counterpart. It uses dmga2py module to communicate
with C++ library dmga. 
"""

import container
import diagram
import geometry
import model
import shape

__all__ = ["container", "diagram", "geometry", "model", "shape"]