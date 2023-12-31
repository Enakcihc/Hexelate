import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path
from PIL import Image
from skimage.transform import downscale_local_mean as downscale

def squircle(a=1, b=1, p = 4, n=100):
    t = np.linspace(0, 0.5 * np.pi, n)
    x = a * np.cos(t) ** (2/p)
    y = b * np.sin(t) ** (2/p)

    xx = np.concatenate((x,np.flip(-x),-x,np.flip(x)))
    yy = np.concatenate((y,np.flip(y),-y,np.flip(-y)))
    
    verts = list(zip(xx,yy))
    codes = [Path.MOVETO] + [Path.LINETO]*(len(verts)-2)+[Path.CLOSEPOLY]
    
    sqcl = Path(verts,codes=codes)

    return sqcl

# size for fully filled 
canvas_width = 6 # inch
canvas_height = 6 # inch
monolisa_raw = np.array(Image.open('sample\Leonardo-Mona-Lisa.jpg'))
dim = monolisa_raw.shape
# print('variable: `monolisa_raw`')
# print(type(monolisa_raw))
# print(dim)
# print(monolisa_raw.dtype)
numX = 60 
factor = round(dim[1] / numX)

numY = int(np.ceil(dim[0] / factor))
numX = int(np.ceil(dim[1] / factor))
monolisa = downscale(monolisa_raw,(factor, factor, 1))

# print('variable: `monolisa`')
# print(type(monolisa))
# print(monolisa.shape)
# print(monolisa.dtype)

ch = 0;
im = np.uint8(monolisa[:,:,ch])

# print('variable: `im`')
# print(type(im))
# print(im.shape)
# print(im.dtype)

print('numX:'+str(numX)+' numY:'+str(numY))

# plt.imshow(im)
# plt.show()

width = canvas_height/numY #0.5 # inch; width of "pixel" or "grid"
dpi = 72
fullSize = (width*dpi) ** 2 

x = np.linspace(0,numX-1,numX)
y = np.linspace(0,numY-1,numY)

xv, yv = np.meshgrid(x, y)
# print(xv[0])
pred = ((xv + yv) % 2) == 1
xv[pred] = np.nan
# print(xv[0])
# print(xv.shape)
# print(im[0])

# -- Calculate the size of markers --
# s_area = (xv/numX)*fullSize
# s_lin = (xv/numX) ** 2 *fullSize
# s_area2 = 4*(xv/numX)*fullSize
s_monolisa = 4*((255-im)/255)*fullSize
#print(s_monolisa.shape)


minY = -0.5 
maxY = (numY - 0.5)
minX = -0.5 
maxX = (numX - 0.5) 

if (numX > numY) :
    margin = 0.5 * (numX - numY)
    minY = minY - margin
    maxY = maxY + margin
else :
    margin = 0.5 * (numY - numX)
    minX = minX - margin
    maxX = maxX + margin

plt.figure(figsize=(canvas_width,canvas_height))
# plt.xlim(0,numX-1)
# plt.ylim(0,numY-1)
plt.scatter(xv,yv,s=s_monolisa,c='k',marker=squircle(),linewidth=0)
plt.subplots_adjust(left=0, right=1, bottom=0, top=1)
plt.xlim(minX,maxX)
plt.ylim(minY,maxY)
plt.gca().invert_yaxis()
plt.savefig('monolisa_squircled.png')
plt.show()
