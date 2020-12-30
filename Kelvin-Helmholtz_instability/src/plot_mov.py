import matplotlib.animation as animation
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np

XN = 480; YN = 601; TN = 200; VN = 9
xl = 16.0; yl = 20.0; tl = 100.0
dx = xl/XN; dy = yl/(YN-1); dt = tl/TN
x = np.arange(0.0, xl, dx); y = np.arange(0.0, yl+dy, dy); t = np.arange(0.0, tl, dt)
fig = plt.figure(figsize = (8, 8))
var = ["r", "p", "u", "v", "w", "bx", "by", "bz", "ps"]
title = ["ρ", "p", "u", "v", "w", "Bx", "By", "Bz", "ψ"]
cols = [cm.jet, cm.jet, cm.BrBG, cm.BrBG, cm.BrBG,
        cm.BrBG, cm.BrBG, cm.BrBG, cm.BrBG]
data = np.zeros((VN, TN, YN, XN))
for i in range(VN):
    for n in range(TN):
        data[i, n, :, :] = np.loadtxt("../data/"+var[i]+"_{0:0>3}.csv".format(n), delimiter = ",")
    dmin = np.min(data[i, :, :, :])
    dmax = np.max(data[i, :, :, :])
    def animate(n):
        plt.gcf().clear()
        plt.axes().set_aspect("equal")
        im = plt.pcolor(x, y, data[i, n, :, :], cmap = cols[i])
        im.set_clim(dmin, dmax)
        fig.colorbar(im)
        plt.title(title[i]+ " (t = {0:.2f})".format(t[n]))
        plt.xlabel("X")
        plt.ylabel("Y")
        print(var[i]+": {0:>3}".format(n), end = "\b\b\b\b\b\b\b", flush = True)
    ani = animation.FuncAnimation(fig, animate, interval = 50, frames = TN)
    ani.save("../mov/"+var[i]+".gif", writer = "imagemagick")
