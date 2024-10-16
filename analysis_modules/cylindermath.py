"""
This module provides functions and a class for performing various mathematical and geometric operations.

The module contains the following functions:

- `moving_average`: Calculates the moving average of a numpy array with a specified window size and padding.
- `norm`: Calculates the Euclidean norm of a given vector.
- `distance_from_axis`: Calculates the distance from the axis of a cylinder to a given point.
- `distance_from_surface`: Calculates the distance from the surface of a cylinder to a given point.

The module also contains the `cylinder` class, which represents a cylinder in 3D space. The cylinder is defined by its radius and the position vectors of two end-points on its long axis.

This module requires the numpy library.

Example:
    To calculate the moving average of an array with a window size of 3 and constant padding:
    
    ```python
    import numpy as np
    import cylindermath

    array = np.array([1, 2, 3, 4, 5])
    window_size = 3
    moving_avg = cylindermath.moving_average(array, window_size, padding='constant')
    print(moving_avg)
    ```

    To calculate the distance from the surface of a cylinder to a point:
    
    ```python
    import numpy as np
    import cylindermath

    radius = 1.0
    rA = np.array([0.0, 0.0, 0.0])
    rB = np.array([0.0, 0.0, 1.0])
    cyl = cylindermath.cylinder(radius, rA, rB)

    rP = np.array([1.0, 1.0, 1.0])
    distance = cylindermath.distance_from_surface(cyl, rP)
    print(distance)
    ```
"""

import numpy as np
from numpy.linalg import norm

def moving_average(array, window_size, padding='constant'):
    """
    Calculates the moving average of an array with padding.

    Args:
      array (numpy.ndarray): The array to calculate the moving average of.
      window_size (int): The size of the moving window.
      padding (str, optional): The type of padding to use. Defaults to 'constant'.

    Returns:
      numpy.ndarray: The moving average of the array.
    """

    if padding not in ['constant', 'reflect', 'edge']:
        raise ValueError(
            'Padding must be one of "constant", "reflect", or "edge".')

    if window_size > len(array):
        raise ValueError(
            'Window size must be less than or equal to the length of the array.')

    padding_width = window_size // 2
    if padding == 'constant':
        padded_array = np.pad(array, padding_width, 'constant')
    elif padding == 'reflect':
        padded_array = np.pad(array, padding_width, 'reflect')
    elif padding == 'edge':
        padded_array = np.pad(array, padding_width, 'edge')

    moving_average = np.convolve(padded_array, np.ones(
        window_size), 'valid') / window_size

    return moving_average


class cylinder():
    """
    Represents a cylinder in 3D space.

    Attributes:
        radius (float): The radius of the cylinder.
        rA (numpy.ndarray): The position vector of end-point A on the long axis if the cylinder.
        rB (numpy.ndarray): The position vector of end-point B on the long axis if the cylinder.
    """

    def __init__(self, radius:float, rA:list, rB:list):
        """
        Initialize the CylinderMath class.

        Args:
            radius (float): The radius of the cylinder.
            rA (list): Position of end-point A.
            rB (list): Position of end-point B.
        """
        rA = np.array(rA)
        rB = np.array(rB)
        
        self.__radius = radius
        self.__rA = rA
        self.__rB = rB
    
    @property
    def radius(self):
        return self.__radius
    
    @property
    def rA(self):
        return self.__rA
    
    @property
    def rB(self):
        return self.__rB


# def norm(A):
#     """
#     Calculate the Euclidean norm of a given vector.

#     Parameters:
#     A (numpy.ndarray): The input vector.

#     Returns:
#     float: The Euclidean norm of the input vector.
#     """
#     sum = 0
#     for i in range(len(A)):
#         sum += A[i] * A[i]
#     return np.sqrt(sum)


def distance_from_axis(cyl:cylinder, rP):
    """
    Calculate the distance from the axis of a cylinder to a given point.

    Args:
        cyl (cylinder): The cylinder object with attributes rA and rB representing two points
             defining the axis of the cylinder.
        rP (numpy.ndarray): The point for which the distance from the cylinder axis is calculated.

    Returns:
        float: The distance from the cylinder axis to the given point.
    """
    
    if len(rP) != 3:
        raise ValueError("The point must be a 3D vector.")
    
    e = cyl.rA - cyl.rB
    d = norm((np.cross(e, rP - cyl.rA)))/norm(e)
    return d


def distance_from_surface(cyl:cylinder, rP):
    """
    Calculate the distance from the surface of a cylinder to a given point.

    Args:
        cyl (cylinder): The cylinder object.
        rP (numpy.ndarray): The point from which to calculate the distance.

    Returns:
        float: The distance from the surface of the cylinder to the given point.
    """
    
    if len(rP) != 3:
        raise ValueError("The point must be a 3D vector.")
    
    return cyl.radius - distance_from_axis(cyl, rP)
