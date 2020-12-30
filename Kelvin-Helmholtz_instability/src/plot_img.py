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
for i in range(VN):
    plt.gcf().clear()
    plt.axes().set_aspect("equal")
    data = np.loadtxt("../data/"+var[i]+"_{0:0>3}.csv".format(TN-1), delimiter = ",")
    im = plt.pcolor(x, y, data, cmap = cols[i])
    fig.colorbar(im)
    plt.title(title[i]+ " (t = {0:.2f})".format(t[TN-1]))
    plt.xlabel("X")
    plt.ylabel("Y")
    #plt.show()
    plt.savefig("../img/"+var[i]+".png")
