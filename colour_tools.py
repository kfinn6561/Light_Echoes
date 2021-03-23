'''
Created on May 27, 2015

@author: kieran
'''
from numpy.random import random as rand


def hex_to_RGB(hx):
#    print hex, hex[0:2],hex[2:4],hex[4:]
    ''' "#FFFFFF" -> [255,255,255] '''
    # Pass 16 to the integer function for change of base
    return [int(hx[i:i+2], 16) for i in range(1,6,2)] 


def RGB_to_hex(RGB):
    ''' [255,255,255] -> "#FFFFFF" '''
    # Components need to be integers for hex to make sense
    RGB = [int(x) for x in RGB]
    return "#"+"".join(["0{0:x}".format(v) if v < 16 else 
                        "{0:x}".format(v) for v in RGB])

def color_dict(gradient):
    ''' Takes in a list of RGB sub-lists and returns dictionary of 
        colors in RGB and hex form for use in a graphing function 
        defined later on '''
    
    return {"hex":[RGB_to_hex(RGB) for RGB in gradient],
            "r":[RGB[0] for RGB in gradient],
            "g":[88 for RGB in gradient],
            "b":[RGB[2] for RGB in gradient]}


def linear_gradient(start_hex="#000000", finish_hex="#FFFFFF", n=10):
    ''' returns a gradient list of (n) colors between two hex colors.
        start_hex and finish_hex should be the full six-digit color string, 
        inlcuding the number sign ("#FFFFFF") 
    '''
    # Starting and ending colors in RGB form
    s = hex_to_RGB( start_hex)
    f = hex_to_RGB(finish_hex)
    # Initilize a list of the output colors with the starting color
    RGB_list = [s]
    # Calcuate a color at each evenly spaced value of t from 1 to n
    for t in range(1, n):
        # Interpolate RGB vector for color at the current value of t
        curr_vector = [ int(s[j] + (float(t)/(n-1))*(f[j]-s[j])) 
                        for j in range(3)]
        # Add it to our list of output colors
        RGB_list.append(curr_vector)

    return color_dict(RGB_list)

def rand_hex_color(num=1):
    ''' Generate random hex colors, default is one, returning a string.
        If num is greater than 1, an array of strings is returned. '''
    colors = [RGB_to_hex([x*255 for x in rand(3)]) for i in range(num)]
    if num == 1:
        return colors[0]
    else:
        return colors


def polylinear_gradient(colors, n):
    ''' returns a list of colors forming linear gradients between all 
        sequential pairs of colors. "n" specifies the total number of 
        desired output colors '''
    # The number of colors per individual linear gradient
    num_out = int(float(n) / (len(colors) - 1))
    gradient_dict = linear_gradient(colors[0], colors[1], num_out) 
    if len(colors) < 1:
        return gradient_dict
    
    for col in range(1, len(colors) - 1):
        nxt = linear_gradient(colors[col], colors[col+1], num_out)
        for k in ("hex", "r", "g", "b"): gradient_dict[k] += nxt[k][1:]
            # Exclude first point to avoid duplicates

    return gradient_dict

kelly_colors = [
    '#FFB300', # Vivid Yellow
    '#803E75', # Strong Purple
    '#007D34', # Vivid Green
    
    '#5D85B3',
    '#C10020', # Vivid Red
    '#CEA262', # Grayish Yellow

    # The following don't work well for people with defective color vision
    
    '#F6768E', # Strong Purplish Pink
    '#00538A', # Strong Blue
    '#FF7A5C', # Strong Yellowish Pink
    '#53377A', # Strong Violet
    '#FF8E00', # Vivid Orange Yellow
    '#B32851', # Strong Purplish Red
    '#F4C800', # Vivid Greenish Yellow
    '#7F180D', # Strong Reddish Brown
    '#93AA00', # Vivid Yellowish Green
    '#593315', # Deep Yellowish Brown
    '#F13A13', # Vivid Reddish Orange
    '#232C16', # Dark Olive Green
    '#817066', # Medium Gray
    '#A6BDD7', # Very Light Blue
    '#FF6800', # Vivid Orange
    ]