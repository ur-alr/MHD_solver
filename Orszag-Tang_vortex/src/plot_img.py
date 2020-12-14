import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np

XN = 200; YN = 200; TN = 100; VN = 9
xl = np.pi; yl = np.pi; tl = np.pi
dx = xl/XN; dy = yl/YN; dt = tl/TN
x = np.arange(0.0, xl, dx); y = np.arange(0.0, yl, dy); t = np.arange(0.0, tl, dt)
fig = plt.figure()
var = ["r", "p", "u", "v", "w", "bx", "by", "bz", "ps"]
title = ["ρ", "p", "u", "v", "w", "Bx", "By", "Bz", "ψ"]
maxs = [8.0, 10.0,  1.5,  1.5,  1.0,  3.0,  2.5,  1.0,  0.2]
mins = [0.0,  0.0, -1.5, -1.5, -1.0, -3.0, -2.5, -1.0, -0.2]
cols = [cm.hot, cm.hot, cm.BrBG, cm.BrBG, cm.BrBG,
        cm.BrBG, cm.BrBG, cm.BrBG, cm.BrBG]
for i in range(VN):
    plt.gcf().clear()
    plt.axes().set_aspect("equal")
    data = np.loadtxt("../data/"+var[i]+"_{0:0>3}.csv".format(TN-1), delimiter = ",")
    im = plt.pcolor(x, y, data, cmap = cols[i])
    im.set_clim(mins[i], maxs[i])
    fig.colorbar(im)
    plt.title(title[i]+ " (t = {0:.2f})".format(t[TN-1]))
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.savefig("../img/"+var[i]+".png")
