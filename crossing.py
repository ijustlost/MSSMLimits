def crossing(x1, x2, y1, y2, y1p, y2p):
    dx = x2-x1
    dy = y2-y1
    dyp = y2p - y1p
    xc = x1 + float(y1-y1p)*dx/(dyp-dy)
    yc = y1 + float(y1-y1p)/(dyp-dy) * dy / dx
    print x1, x2, y1, y2, y1p, y2p, "xing:", xc, yc


crossing(1, 2, 1, 2, 1.5, 1.5)
crossing(1, 2, 1, 2, 2, 1)
crossing(1, 2, 1, 2, 3, 3)
