import matplotlib.animation as animation
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np

XN = 200; YN = 200; TN = 100; VN = 9
xl = 2.0*np.pi; yl = 2.0*np.pi; tl = np.pi
dx = xl/XN; dy = yl/YN; dt = tl/TN
x = np.arange(0.0, xl, dx); y = np.arange(0.0, yl, dy); t = np.arange(0.0, tl, dt)
fig = plt.figure()
var = ["r", "p", "u", "v", "w", "bx", "by", "bz", "ps"]
title = ["ρ", "p", "u", "v", "w", "Bx", "By", "Bz", "ψ"]
maxs = [8.0, 10.0,  1.5,  1.5,  1.0,  3.0,  2.5,  1.0,  0.2]
mins = [0.0,  0.0, -1.5, -1.5, -1.0, -3.0, -2.5, -1.0, -0.2]
cols = [cm.jet, cm.jet, cm.BrBG, cm.BrBG, cm.BrBG,
        cm.BrBG, cm.BrBG, cm.BrBG, cm.BrBG]
for i in range(VN):
    def animate(n):
        plt.gcf().clear()
        plt.axes().set_aspect("equal")
        data = np.loadtxt("../data/"+var[i]+"_{0:0>3}.csv".format(n), delimiter = ",")
        im = plt.pcolor(x, y, data, cmap = cols[i])
        im.set_clim(mins[i], maxs[i])
        fig.colorbar(im)
        plt.title(title[i]+ " (t = {0:.2f})".format(t[n]))
        plt.xlabel("X")
        plt.ylabel("Y")
        print(var[i]+": {0:>3}".format(n), end = "\b\b\b\b\b\b\b", flush = True)
    ani = animation.FuncAnimation(fig, animate, interval = 50, frames = TN)
    ani.save("../mov/"+var[i]+".gif", writer = "imagemagick")
