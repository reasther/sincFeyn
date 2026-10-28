# g2_02, translated from the generated C by c2mma.py
def pt(x):
    x0 = x[0]
    x1 = x[1]
    x2 = x[2]
    x3 = x[3]
    x4 = x[4]
    x5 = x[5]
    # node F2: seg (0, 4) kids ()  L=2
    n_F2_0_4_x_W = [[[0.0 for _ in range(4)] for _ in range(1)] for _ in range(1)]
    U = x0*x1 + x0*x2 + x0*x3 + x0*x5 + x1*x3 + x1*x4 + x1*x5 + x2*x3 + x2*x4 + x2*x5 + x3*x4 + x4*x5
    iU = 1/U
    V = ((x0*x0)*x1 + (x0*x0)*x2 + (x0*x0)*x3 + (x0*x0)*x5 + x0*(x1*x1) + 2*x0*x1*x2 + 2*x0*x1*x3 + 2*x0*x1*x5 + x0*(x2*x2) + 2*x0*x2*x3 + 2*x0*x2*x5 + x0*(x3*x3) + (x1*x1)*x3 + (x1*x1)*x4 + (x1*x1)*x5 + 2*x1*x2*x3 + 2*x1*x2*x4 + 2*x1*x2*x5 + x1*(x3*x3) + 2*x1*x3*x4 + (x2*x2)*x3 + (x2*x2)*x4 + (x2*x2)*x5 + x2*(x3*x3) + 2*x2*x3*x4 + (x3*x3)*x4)*iU
    cp0 = (x1*x4 + x2*x3 + x2*x4 + x2*x5 + x3*x4 + x4*x5)*iU
    cpp0 = (x1*x3 - x2*x5)*iU
    cp1 = (x0*x2 + x2*x3 + x2*x4 + x2*x5 + x3*x4 + x4*x5)*iU
    cpp1 = (-x0*x2 - x0*x3 - x2*x3 - x2*x4 - x2*x5 - x3*x4)*iU
    cp2 = (-x0*x1 - x0*x3 - x0*x5 - x1*x3 - x1*x4 - x1*x5)*iU
    cpp2 = (x0*x1 + x0*x5 + x1*x3 + x1*x4 + x1*x5 + x4*x5)*iU
    cp3 = (x0*x2 - x1*x4)*iU
    cpp3 = (x0*x1 + x0*x5 + x1*x4 + x1*x5 + x2*x5 + x4*x5)*iU
    w0_1 = (x3 + x5)*iU
    w0_2 = (x3 + x5)*iU
    w0_3 = (-x1 - x2)*iU
    w1_2 = (x0 + x3 + x4 + x5)*iU
    w1_3 = (x0 + x4)*iU
    w2_3 = (x0 + x4)*iU
    n_F2_0_4_x_W[0][0][0] += cpp2*(256+cp3*(-256+cpp1*(-128+cpp0*(256)+cp0*(256))+cp1*(-256+cpp0*(512)+cp0*(512))+cp0*(128))+cpp3*(-128+cpp0*(-128+cpp1*(256)+cp1*(512))+cp1*(-128+cp0*(512))+cpp1*(cp0*(256)))+cpp1*(-128+cpp0*(-128))+cp1*(cpp0*(-256)+cp0*(-128))+cp0*(-128))+cp2*(cp1*(-128+cpp0*(-128+cpp3*(256)+cp3*(256))+cp3*(-128+cp0*(256))+cpp3*(cp0*(256)))+cpp3*(cpp1*(128)+cpp0*(-256)+cp0*(-128))+cpp1*(-256+cp0*(128))+cpp0*(128+cp3*(-128))+cp3*(-128))+cp1*(256+cpp3*(-128+cpp0*(128))+cp0*(-128+cp3*(-128))+cpp0*(-256))+cpp1*(cp3*(128+cpp0*(-128)+cp0*(-256))+cpp3*(cp0*(-128))+cpp0*(-128))+cpp3*(-128)+cp0*(-128)
    n_F2_0_4_x_W[0][0][1] += cpp2*(w1_3*(128+cpp0*(-128)+cp0*(-256))+w0_3*(128+cpp1*(-512)+cp1*(-1024))+w0_1*(256+cpp3*(-128)+cp3*(-256)))+cp1*(w2_3*(256+cpp0*(-256)+cp0*(-128))+w0_2*(128+cpp3*(-256)+cp3*(-128))+w0_3*(128+cp2*(-512)))+cpp3*(w0_2*(256+cpp1*(-128))+w1_2*(-128)+w0_1*(-128))+cp2*(w1_3*(128+cp0*(-128))+w0_3*(128)+w0_1*(cp3*(-128)))+cpp1*(w2_3*(cpp0*(-128))+w0_3*(128)+w0_2*(128))+cp0*(w2_3*(-128)+w1_3*(256)+w1_2*(-128))
    iU2 = iU*iU
    n_F2_0_4_x_W = [[[v*iU2 for v in r] for r in p] for p in n_F2_0_4_x_W]
    n_F2_0_4_x_V = V
    tot = 0
    t1 = [0.0 for _ in range(4)]
    t1 = [a+b for a, b in zip(t1, n_F2_0_4_x_W[0][0])]
    iV = 1/(n_F2_0_4_x_V)
    r = 0
    ivp = 1
    ivp = ivp*iV
    r = r+1*t1[1]*ivp
    ivp = ivp*iV
    r = r+1*t1[0]*ivp
    tot = tot+r
    return tot
GLX = [0.1127016653792583,0.5,0.8872983346207417]
GLW = [0.27777777777777785,0.4444444444444444,0.27777777777777785]
def f(sf):
    acc = 0.0
    for u, wg in zip(GLX, GLW):
        s = [1.0]*6
        s[0] = sf[1]/1
        s[3] = sf[2]/1
        s[5] = sf[3]/1
        s[1] = sf[0]*u/1
        s[2] = sf[0]*(1-u)/1
        sa, sb = sf[0]*(1-u), sf[0]*u
        acc += wg*1*pt(s)
    return 1.0/256*sf[0]*acc
