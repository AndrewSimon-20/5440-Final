import sympy as sym
import matplotlib.pyplot as plt
import math

kf = sym.Symbol('kf')
rs = sym.Symbol('Rs')
kr = sym.Symbol('kr')
rsStar = sym.Symbol('RsStar')
ke = sym.Symbol('ke')
vs = sym.Symbol('vs')
keStar = sym.Symbol('keStar')
q = sym.Symbol('q')
nc = sym.Symbol('nc')
km = sym.Symbol('km')
riStar = sym.Symbol('RiStar')
ri = sym.Symbol('Ri')
kdeg = sym.Symbol('kdeg')
# Lc = sym.Symbol('Lc')
Lc = (nc*(q+kr*rsStar)/(km+kf*rs*nc))

# Equations 1 -4 from the notes with Lc subbed in.
solution = sym.solve((-kf*Lc*rs+kr*rsStar-ke*rs+vs, kf*Lc*rs-kr*rsStar-keStar*rsStar, ke*rs+keStar*rsStar-kdeg*(ri+riStar), keStar*rsStar-kdeg*riStar), (rsStar, rs, riStar, ri))


# To be physically possible we must take the positive root.
print('\n'+'The Solution for RsStar is:')
print(solution[0][0])
print('\n'+'The Solution for RsStar is:')
print(solution[0][1])
print('\n'+'The Solution for RiStar is:')
print(solution[0][2])
print('\n'+'The Solution for Ri is:')
print(solution[0][3])
print('\n'+'Thus Lc is:')
print(sym.simplify((nc*(q+kr*solution[0][0])/(km+kf*solution[0][1]*nc))))







xAxis = []
yAxis = []
for i in range(1,1000,1):
    j = i/1000
    xAxis.append(j)
    z = j

    # Plugging in variables

    ke = 10 ** -4
    keStar = 5 * 10 ** -3
    kf = 5.14 * 10 ** -21
    kr = 2.5 * 10 ** -2
    kdeg = 8 * 10 ** -4
    vs = 18
    q = 1000
    nc = 3 * 10 ** 8
    gammadot = 100
    dl = 10 ** -10
    km = ((gammadot * dl ** 2) / z) ** (1 / 3)

# I just copied and pasted the positive roots in here.
    rsStar = (ke*keStar*km + ke*km*kr + keStar*kf*nc*q + keStar*kf*nc*vs - math.sqrt(ke**2*keStar**2*km**2 + 2*ke**2*keStar*km**2*kr + ke**2*km**2*kr**2 + 2*ke*keStar**2*kf*km*nc*q + 2*ke*keStar**2*kf*km*nc*vs + 2*ke*keStar*kf*km*kr*nc*q + 2*ke*keStar*kf*km*kr*nc*vs + keStar**2*kf**2*nc**2*q**2 - 2*keStar**2*kf**2*nc**2*q*vs + keStar**2*kf**2*nc**2*vs**2))/(2*keStar**2*kf*nc)
    rs = (-ke*keStar*km - ke*km*kr - keStar*kf*nc*q + keStar*kf*nc*vs + math.sqrt(ke**2*keStar**2*km**2 + 2*ke**2*keStar*km**2*kr + ke**2*km**2*kr**2 + 2*ke*keStar**2*kf*km*nc*q + 2*ke*keStar**2*kf*km*nc*vs + 2*ke*keStar*kf*km*kr*nc*q + 2*ke*keStar*kf*km*kr*nc*vs + keStar**2*kf**2*nc**2*q**2 - 2*keStar**2*kf**2*nc**2*q*vs + keStar**2*kf**2*nc**2*vs**2))/(2*ke*keStar*kf*nc)
    riStar = (ke*keStar*km + ke*km*kr + keStar*kf*nc*q + keStar*kf*nc*vs - math.sqrt(ke**2*keStar**2*km**2 + 2*ke**2*keStar*km**2*kr + ke**2*km**2*kr**2 + 2*ke*keStar**2*kf*km*nc*q + 2*ke*keStar**2*kf*km*nc*vs + 2*ke*keStar*kf*km*kr*nc*q + 2*ke*keStar*kf*km*kr*nc*vs + keStar**2*kf**2*nc**2*q**2 - 2*keStar**2*kf**2*nc**2*q*vs + keStar**2*kf**2*nc**2*vs**2))/(2*kdeg*keStar*kf*nc)
    ri = (-ke*keStar*km - ke*km*kr - keStar*kf*nc*q + keStar*kf*nc*vs + math.sqrt(ke**2*keStar**2*km**2 + 2*ke**2*keStar*km**2*kr + ke**2*km**2*kr**2 + 2*ke*keStar**2*kf*km*nc*q + 2*ke*keStar**2*kf*km*nc*vs + 2*ke*keStar*kf*km*kr*nc*q + 2*ke*keStar*kf*km*kr*nc*vs + keStar**2*kf**2*nc**2*q**2 - 2*keStar**2*kf**2*nc**2*q*vs + keStar**2*kf**2*nc**2*vs**2))/(2*kdeg*keStar*kf*nc)

    kss = keStar * kf / (ke * (kr + keStar))

# Also note I didn't know how to incorporate the approximation that LcKss << 1 for these long expressions so I just plotted it as is.
# NOTE BECAUSE I DON'T KNOW THE MITOGENIC SIGNAL I AM GOING TO JUST ADD A CONSTANT WHICH NORMALIZES MY FUNCTION TO PEAK AT 100% AS I THINK IT WOULD LOGICALLY.
# I'M STILL GOING TO CALL THIS CONSTANT GAMMA.
    gamma = 1/(1525*12.3)
    yAxis.append(gamma*(rsStar+riStar))

plt.scatter(xAxis, yAxis)
plt.xlabel('Z (meters)')
plt.ylabel('mitotic Response')
plt.show()
