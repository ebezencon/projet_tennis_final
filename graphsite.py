#! /usr/bin/python3
import csv
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D as ax
import numpy as np

#####

"Entrez le path de votre fichier:"

file_path = "trajectoire.csv"
file_path_alea="impacts.csv"                #mettre le même path que celui de la simulation c. (même si la faisabilité n'a pas été calculée)
num=str(np.random.randint(1000,10000))      #l'indice de la simulation pour refresh les images sur site après simulation
print("Image index:"+num)                
####

with open(file_path, newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    x=[]
    y=[]
    z=[]
    for row in reader:
        x.append(float(row['X']))
        y.append(float(row['Y']))
        z.append(float(row['Z']))


def plottenniscourt2D():
    #plot tennis court en 2D
    x=[-5.49, -5.49,  5.49, 5.49,-5.49,-5.49, 5.49,  5.49,  4.12, 4.12,-4.12, -4.12,-4.12,4.12,   0,   0,4.12,-4.12]
    y=[11.89,-11.89,-11.89,11.89,11.89,    0,    0,-11.89,-11.89,11.89,11.89,-11.89,-6.40,-6.4,-6.4, 6.4, 6.4, 6.40]
    plt.plot(x,y,'tab:green')
    plt.plot([-5.49,5.49],[0,0],'tab:gray')
    plt.plot(12,12)
    plt.plot(-12,-12)
    
def plothalftenniscourt2D():
    #plot tennis court en 2D mais juste la partie du haut
        x=[-5.46,5.46,5.46,-5.46,-5.46,-4.12,-4.12,4.12 ,4.12,4.12,0,   0,  0 ,-4.12]
        y=[0    ,0   ,11.89,11.89,0   ,    0,11.89,11.89,0   ,6.40,6.40,0,6.40, 6.40]
        plt.plot(x,y,'tab:green')
        plt.plot([-5.49,5.49],[0,0],'tab:gray')

def plottenniscourt3D():
    x=[-5.49, -5.49,  5.49, 5.49,-5.49,-5.49, 5.49,  5.49,  4.12, 4.12,-4.12, -4.12,-4.12,4.12,   0,   0,4.12,-4.12]
    y=[11.89,-11.89,-11.89,11.89,11.89,    0,    0,-11.89,-11.89,11.89,11.89,-11.89,-6.40,-6.4,-6.4, 6.4, 6.4, 6.40]
    z=[0]*len(x)
    ax.plot3D(x, y, z, 'green')
    #On fait ça pour avoir un graph à l'echelle 1:1 car set_equal marche pas en 3d
    x=[12]
    y=[12]
    z=[12]
    ax.plot3D(x, y, z)  
    x=[-12]
    y=[-12]
    z=[-12]
    ax.plot3D(x, y, z)   

def plotnet3D():
    z=[0,0.91,0.91,0]
    y=[0,0,0,0]
    x=[-5.49,-5.49,5.49,5.49]
    ax.plot3D(x,y,z,'grey')
    

#Ligne qui évite problème de plot si la balle n'a jamais rebondi
if min(z)==0.00:
    stop=z.index(0.00)
else :
    stop = len(z)-1
    

#Vue du dessus 2D
plt.figure(figsize=(8,8))
plottenniscourt2D()
plt.xticks([-5.49,-4.12,0,4.12,5.49],['-5.49    ','    -4.12','0','4.12    ','    5.49'])
plt.yticks([-11.89,-6.40,0,6.40,11.89],['-11.89','-6.40','0','6.40','11.89'])
plt.plot(x[:stop],y[:stop],'#000000')
plt.plot(x[stop-1:],y[stop-1:],'-.',color='#000000')
plt.scatter(x[stop],y[stop],color='yellow',alpha=1,edgecolors='black',zorder=10)
plt.title("Vue du dessus du terrain de la trajectoire du coup sans aléas")
plt.savefig("image/vue_dessus"+num+".SVG")

#Plot 3D
plt.figure(figsize=(8,8))
ax = plt.axes(projection="3d")
ax.plot3D(x[:stop], y[:stop], z[:stop], '#000000',zorder=10)
ax.plot3D(x[stop-1:], y[stop-1:], z[stop-1:],'--',color='#000000',zorder=10)
ax.axis('off')
ax.margins(-0.49,-0.49,0)
ax.view_init(elev=20,azim=-65)
ax.dist=6
plottenniscourt3D()
plt.title("Plot 3D de la trajectoire du coup sans aléas")
plotnet3D()
plt.savefig("image/vue_plot_3D"+num+".SVG")


#Plot aléas
with open(file_path_alea, newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    x=[]
    y=[]
    z=[]
    i=1
    for row in reader:
        x.append(float(row['X']))
        y.append(float(row['Y']))
        z.append(float(row['Z']))
        i+=1
        if i >1250:                              #On récolte pas plus de points car on en plot max 1250 sinon c'est juste trop de points
            break


if x[0] != 0.0 or y[0] != 0.0 or z[0] != 0.0:    #On plot pas si on a les probas off
    plt.figure(figsize=(8,8))
    xpoints=[]
    ypoints=[]
    for i in range(1,min(len(x),1250)):          #On limite le plot de point a 1250 pour que ça reste visible
        if x[i] != -1:                           #on plote pas les points qui sont sortis des limites du terrains (caractérisés par x=-1)
            xpoints.append(x[i])
            ypoints.append(y[i])
    #plot tous les impacts maintenant
    plt.scatter(xpoints,ypoints,alpha=0.5,color='yellow',edgecolors='black')
    #Ploter le cercle d'incetitudes
    theta = np.linspace(0, 2*np.pi, 100)
    x1 = z[0]*np.cos(theta)+x[0]                 #On avait mis la taille du rayon du cercle d'incertitudes dans z[0] tandis que x[0] et y[0] sont les coords pt impact sans aléas.
    x2 = z[0]*np.sin(theta)+y[0]
    plt.plot(x1,x2,'--',color='red',label="Cercle de similitude (rayon=%.2fm)"%(z[0]))
    plt.scatter(x[0],y[0],color='red')
    plothalftenniscourt2D()
    ax=plt.gca()
    ax.set_aspect(1)            
    plt.legend(loc='upper right')
    plt.title("Plot des points d'impacts obtenus lors de la simulation avec aléas")
    plt.savefig("image/vue_impacts_imprécisions"+num+".SVG")
else:
    #Si on a pas de probas, on creer une figure extrèmement petite et invisible qui sera sur site mais invisible comme ça pas besoin de gerer le cas plus difficilement.
    plt.figure(figsize=(0.1,0.1))
    plt.savefig("image/vue_impacts_imprécisions"+num+".SVG")