# pylint: disable=I0011,C0103,C0111,W0401,W0614,W0622
from __future__ import division, print_function
from visual import *
from visual.graph import *
scene.width = 800

X_momentum = gdisplay(x=0, y=400, height=200)
px_Alpha = gcurve(color=color.red)
px_Au = gcurve(color=color.blue)
px_total = gcurve(color=color.cyan)

Y_momentum = gdisplay(x=0, y=600, height=200)
py_Alpha = gcurve(color=color.red)
py_Au = gcurve(color=color.blue)
py_total = gcurve(color=color.cyan)


# Constants
massAu = (79 + 118) * 1.7e-27
massAlpha = (2 + 2) * 1.7e-27
qAu = 2 * 1.6e-19
qAlpha = 79 * 1.6e-19
oofpez = 9e9  # one-over-four-pi-epsilon-zero
deltat = 1e-23
initx = -1e-13

Au = sphere(pos=vector(0, 0, 0), radius=8e-15,
            color=color.yellow, make_trail=True)
b = 5e-15
Alpha = sphere(pos=vector(initx, b, 0), radius=2e-15,
               color=color.red, make_trail=True)
Alpha.visible = True
pAu = massAu * vector(0, 0, 0)
pAlpha = vector(1.043e-19, 0, 0)
t = 0
while t < 2e-20:
    rate(10000)
    Alpha.pos = Alpha.pos + (pAlpha / massAlpha) * deltat
    Au.pos = Au.pos + (pAu / massAu) * deltat
    f = (oofpez * qAu * qAlpha) / \
        ((Alpha.pos - pAu).mag2) * (Alpha.pos - pAu).norm()
    pAlpha += f * deltat
    pAu -= f * deltat
    px_Alpha.plot(pos=(t, pAlpha.x))
    px_Au.plot(pos=(t, pAu.x))
    px_total.plot(pos=(t, (pAu + pAlpha).x))
    py_Alpha.plot(pos=(t, pAlpha.y))
    py_Au.plot(pos=(t, pAu.y))
    py_total.plot(pos=(t, (pAu + pAlpha).y))
    t = t + deltat

qry = [5e-15, 10e-15, 100e-15]
for b in qry:
    Au = sphere(pos=vector(0, 0, 0), radius=8e-15,
            color=color.yellow)
    Au.visible=False
    Alpha = sphere(pos=vector(initx, b, 0), radius=2e-15)
    Alpha.visible = False
    pAu = massAu * vector(0, 0, 0)
    pAlpha = vector(1.043e-19, 0, 0)
    t = 0
    while t < 2e-20:
        Alpha.pos = Alpha.pos + (pAlpha / massAlpha) * deltat
        Au.pos = Au.pos + (pAu / massAu) * deltat
        f = (oofpez * qAu * qAlpha) / \
            ((Alpha.pos - pAu).mag2) * (Alpha.pos - pAu).norm()
        pAlpha += f * deltat
        pAu -= f * deltat
        t = t + deltat
    theta = acos(pAlpha.norm().x)
    print(degrees(theta))
    del Alpha

qry = [45, 90, 135]
for q in qry:
    lb = 1e-13
    rb = 0
    for loop in range(20):
        b = (lb + rb) / 2
        Au = sphere(pos=vector(0, 0, 0), radius=8e-15,
            color=color.yellow)
        Au.visible=False
        Alpha = sphere(pos=vector(initx, b, 0), radius=2e-15)
        Alpha.visible = False
        pAu = massAu * vector(0, 0, 0)
        pAlpha = vector(1.043e-19, 0, 0)
        t = 0
        while t < 2e-20:
            Alpha.pos = Alpha.pos + (pAlpha / massAlpha) * deltat
            Au.pos = Au.pos + (pAu / massAu) * deltat
            f = (oofpez * qAu * qAlpha) / \
                ((Alpha.pos - pAu).mag2) * (Alpha.pos - pAu).norm()
            pAlpha += f * deltat
            pAu -= f * deltat
            t = t + deltat
        theta = acos(pAlpha.norm().x)
        if theta > radians(q):
            rb = b
        else:
            lb = b
        del Alpha
    print((lb + rb) / 2)
