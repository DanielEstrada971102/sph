from numpy import loadtxt, arange, pi, exp, sin, cos, array, sqrt
import matplotlib.pyplot as plt
from glob import glob
from statistics import mean

partID = 234
X = 1
Y = 2
VX = 3
VY = 4

Lx = 1
Ly = 1
dvc_eta = 8.9e-4;

Lref_x = Lx
Lref_y = Ly


def u(x, y, t, Re):
    b = -8 * pi * pi/ Re
    Vo = 5e-1
    return Vo * exp(b * t) * cos(2*pi*x / Lref_x) * sin(2*pi*y / Lref_y)

def v(t, x, y, Re):
    b = -8 * pi * pi/ Re
    Vo = 5e-1
    return -1 * Vo * exp(b * t) * cos(2*pi*y / Lref_y) * sin(2*pi*x / Lref_x)


def main():

    dt = .5e-5
    tTotal = 400*dt
    V_max = []
    t = arange(0, tTotal + dt, dt)

    x = []
    y = []
    velx = []
    vely = []


    for file_name in glob("output/*"):
        data = loadtxt(file_name, usecols=(0, X, Y, VX, VY))

        x.append(data[partID, X])
        y.append(data[partID, Y])
        velx.append(data[partID, VX])
        vely.append(data[partID, VY])

        V = sqrt( data[:,VX] * data[:, VX] +  data[:,VY] * data[:, VY])
        V_max.append(max(V))

    
    x , y = array(x), array(y)

    Re = mean(V_max) * Lx / dvc_eta
    print(Re)
    plt.plot(t, velx, 'ok')
    plt.plot(t, u(x,y,t, Re), '--k')
    #plt.plot(t, vely, 'or')
    #plt.plot(t, -1*u(t), '--r') 
    plt.title("%.2f"%Re)
    plt.show()


if __name__ == '__main__':
    main()