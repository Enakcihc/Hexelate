import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib.cm as cm


df = pd.read_csv('Hexalated.csv')
style = 1 # 1 = dot scatter; 

#print(df['px'][0:10])
canvas_width = 8 # inch
canvas_height = canvas_width # inch

# print("Min col: ",min(df["col"].values))
# print("Max col: ",max(df["col"].values))
# print("Min row ",min(df["row"].values))
# print("Max row: ",max(df["row"].values))

numX = (max(df["col"].values) - min(df["col"].values) )/2 + 1
numY = (max(df["row"].values) - min(df["row"].values) )*np.sqrt(3)/2 + 1
dP = (max(df["px"].values) - min(df["px"].values)) / (numX - 1)
minY = min(df["py"].values) - 0.5 * dP
maxY = max(df["py"].values) + 0.5 * dP
minX = min(df["px"].values) - 0.5 * dP
maxX = max(df["px"].values) + 0.5 * dP

print("dP :", dP)
print("numX :", numX)
print("deltaX / numX:", (maxX-minX) / numX)
print("numY :", numY)
print("deltaY / numY :", (maxY-minY)/numY)
print("num max :", max([numX,numY]))

if (numX > numY) :
    margin = 0.5 * dP * (numX - numY)
    minY = minY - margin
    maxY = maxY + margin
else :
    margin = 0.5 * dP * (numY - numX)
    minX = minX - margin
    maxX = maxX + margin

width = min([canvas_width/numX,canvas_height/numY]) #0.5 # inch; width of "pixel" or "grid"
dpi = 72
fullSize = (width*dpi) ** 2 
fullSize_hex = fullSize / (np.cos(np.pi/6)) ** 2
print(np.cos(np.pi/6))

plt.figure(figsize=(canvas_width,canvas_height))


if style == 1 :
    s_area = ((255-df["c_mean"])/255)*fullSize
    plt.scatter(df["px"].values,df["py"].values,s=s_area,c='k',marker='o',linewidth=0)
elif style == 2 :
    
    c_mean = 255-df["c_mean"].values
    plt.scatter(df["px"].values,df["py"].values,s=fullSize_hex,c=c_mean,cmap="Greys",marker=(6, 0, 0),linewidth=0)
elif style == 3 :
    c_mean = df["c_mean"].values
    plt.scatter(df["px"].values,df["py"].values,s=fullSize_hex,c=c_mean,cmap="jet",marker=(6, 0, 0),linewidth=0)




plt.subplots_adjust(left=0, right=1, bottom=0, top=1)
#plt.gca().invert_yaxis()
plt.xlim(minX,maxX)
plt.ylim(maxY,minY)
plt.show()