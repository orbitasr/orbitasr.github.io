import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from scipy.integrate import quad
import math
from matplotlib.patches import Circle

def v(u, l):
    v = -u + (l ** 2) * (u ** 2) / 2 - (l ** 2) * (u ** 3)
    return v

def gerarPotencial1():
    
    from js import m, mode

    momento = float(m)

    if momento == 0:
        momento = 0.000001
    else:
        momento = momento

    
    global vmax, vmin, l

    l = momento

    if mode == 'light':
        eixos = '#5e6469'
        borda = '#f7f7f7'
    else:
        eixos = 'white'
        borda = '#0E1117'

    if l > np.sqrt(12):
        coef = [0, - 3 * (l ** 2), l ** 2, -1]
        a = np.roots(coef)
        umax = sorted(a)[1].real
        umin = sorted(a)[0].real
        vmin = v(umin, l)
        vmax = v(umax, l)
        vlim = vmax
        rmax = 2 / umin
        thereIsValue = True
        # t = 'O mínimo e o máximo da energia potencial efetiva são ' + str(vmin) + ' e ' + str(vmax) + ' (ver pontos no gráfico)'
        # display(t, target="graph", append=True)
        r = np.arange(2, rmax, rmax / 30000)
        u = 1 / r

        #PLOT DO POTENCIAL

        fig1 = plt.figure()
        plt.subplot(1, 1, 1)
        plt.axhline(0, linewidth=0.3, color='white')
        plt.plot(1.477 * r, v(u, l), color='white')
        plt.plot([1.477 / umin, 1.477 / umax], [vmin, vmax], 'bo', color="gold")
        plt.xlabel("r [km]")
        plt.axis([-1, 1.2 * rmax, -0.5, vlim + 0.1])
        # plt.axis([1.477/umax*0.05, 1.477/umin*1.2, -0.5, vlim*1.2])
        ax = plt.gca()
        ax.spines['bottom'].set_color(eixos)
        ax.tick_params(axis='x', colors=eixos)
        ax.tick_params(axis='y', colors=eixos)
        ax.spines['top'].set_color(eixos)
        ax.spines['right'].set_color(eixos)
        ax.spines['left'].set_color(eixos)
        ax.xaxis.label.set_color(eixos)
        ax.yaxis.label.set_color(eixos)
        fig1.patch.set_facecolor(borda)
        ax.set_facecolor('black')
        display(plt, target="graph", append=True)

    else:
        vlim = 0
        rmax = 80
        r = np.arange(2, rmax, rmax / 300)
        u = 1 / r

        #PLOT POTENCIAL

        fig1 = plt.figure()
        plt.subplot(1, 1, 1)
        plt.plot(1.477 * r, v(u, l), color="white")
        plt.axhline(0, linewidth=0.3, color='white')
        plt.xlabel("r [km]")
        plt.axis([0, 1.477 * rmax, -0.5, vlim + 0.1])
        ax = plt.gca()
        ax.spines['bottom'].set_color(eixos)
        ax.tick_params(axis='x', colors=eixos)
        ax.tick_params(axis='y', colors=eixos)
        ax.spines['top'].set_color(eixos)
        ax.spines['right'].set_color(eixos)
        ax.spines['left'].set_color(eixos)
        ax.xaxis.label.set_color(eixos)
        ax.yaxis.label.set_color(eixos)
        fig1.patch.set_facecolor(borda)
        ax.set_facecolor("black")
        display(plt, target="graph", append=True)

def gerarOrbita1():
     
    from js import e, o, mode

    if mode == 'light':
        eixos = '#5e6469'
        borda = '#f7f7f7'
    else:
        eixos = 'white'
        borda = '#0E1117'

    E = float(e)

    if E < vmin:
        display("ATENÇÃO: Você escolheu um valor de parâmetro de energia menor que o valor mínimo permitido. Escolha novamente um valor dentro do limite especificado antes de continuar.", target="graph2", append=True)
    else:
        norbit = int(o)
        if E == 0:
            E = E + 1e-10
        coef = [- (l ** 2), (l ** 2) / 2, -1, -E]
        roots = np.roots(coef)
        tp1 = roots[2]
        tp2 = roots[1]
        tp3 = roots[0]

        eps = 0.00000001
        rst = 10 * l / (E + 0.5)
        ust = 1 / rst
        correction = 1.477

        if l > math.sqrt(12):
            if E < 0 and ust < tp2.real:
                u1 = tp1.real * (1 + eps)
                u2 = tp2.real * (1 - eps)
            elif 0 < E < vmax and ust < tp2.real:
                u1 = ust
                u2 = tp2.real * (1 - eps)
                norbit = 1
            elif E < vmax and ust > tp3.real:
                u1 = 0.5
                u2 = tp3.real * (1 + eps)
                norbit = 0.5
            elif E > vmax:
                u1 = ust
                u2 = 0.5
                norbit = 0.5
        else:
            if E >= 0:
                u1 = ust
                u2 = 0.5
                norbit = 0.5
            else:
                u1 = tp1.real * (1 + eps)
                u2 = 0.5
                norbit = 0.5


        def theta(w):
            theta = (l / (2 ** (1 / 2))) * ((E - v(w, l)) ** (-1 / 2))
            return theta


        delphi, erro = quad(theta, u1, u2)

        n = 1000
        uc = np.arange(u1, u2, (u2 - u1) / n)
        ud = np.arange(u2, u1, (u1 - u2) / n)

        phi1 = []
        for i in range(len(uc)):
            a = quad(theta, u1, uc[i])
            phi1.append(abs(a[0]))

        phi2 = []
        for j in range(len(ud)):
            b = quad(theta, u2, ud[j])
            phi2.append(abs(b[0]))

        if norbit == 0.5:
            utotal = uc
        else:
            utotal = np.concatenate([uc, ud] * (norbit))

        accphi = [0] * (len(utotal))

        if norbit == 0.5:
            accphi = phi1
            x = [0] * (len(uc))
            y = [0] * (len(uc))
            for i in range(len(uc)):
                x[i] = (math.cos(accphi[i])) / utotal[i] * correction
                y[i] = (math.sin(accphi[i])) / utotal[i] * correction
        else:
            for i in range(norbit):
                for j in range(n):
                    accphi[j + (2 * i * n)] = 2 * i * delphi + phi1[j]
                    accphi[j + ((2 * i + 1) * n)] = ((2 * i) + 1) * delphi + phi2[j]
            x = [0] * (2 * norbit * n)
            y = [0] * (2 * norbit * n)
            for i in range(2 * norbit * n):
                x[i] = ((math.cos(accphi[i])) / utotal[i]) * correction
                y[i] = ((math.sin(accphi[i])) / utotal[i]) * correction

        #PLOT DA ÓRBITA

        fig2 = plt.figure()
        plt.plot(x, y, color="gold")
        plt.xlabel("x [km]")
        plt.ylabel("y [km]")
        circle = Circle((0, 0), 2 * 1.477, color='dimgrey')
        plt.gca().add_patch(circle)
        plt.gca().set_aspect('equal')
        plt.axis(
            [(-1 / u1 - 1) * 1.477, (1 / u1 + 1) * 1.477, (-1 / u1 - 1) * 1.477, (1 / u1 + 1) * 1.477])
        ax = plt.gca()
        ax.spines['bottom'].set_color(eixos)
        ax.tick_params(axis='x', colors=eixos)
        ax.tick_params(axis='y', colors=eixos)
        ax.spines['top'].set_color(eixos)
        ax.spines['right'].set_color(eixos)
        ax.spines['left'].set_color(eixos)
        ax.xaxis.label.set_color(eixos)
        ax.yaxis.label.set_color(eixos)
        fig2.patch.set_facecolor(borda)
        ax.set_facecolor("black")
        display(plt, target="graph2", append=True)

def gerarOrbita2(choice):
    
    from js import D1, D2, mode

    if choice == 1: d = float(D1)
    if choice == 2: d = float(D2)

    if mode == 'light':
        eixos = '#5e6469'
        borda = '#f7f7f7'
    else:
        eixos = 'white'
        borda = '#0E1117'
    
    if d == 0:
        d = 0.00001
    else:
        d = d

    b = 2 * d

    rst = 50
    norbit = 10

    r = np.arange(0.1, 80, 0.1)
    u = 1 / r
    k = 1 / (b ** 2)


    def w(u):
        w = u ** 2 - 2 * (u ** 3)
        return w


    umax = 1 / 3
    wmax = 1 / 27
    ust = 1 / rst

    coef = [-2, 1, 0, -k]
    roots = np.roots(coef)
    tp2 = roots[1]
    tp3 = roots[0]

    eps = 0.000000001

    if k == wmax:
        phi = np.arange(0, 2 * math.pi, math.pi / 500)
        rmax = [1 / umax] * len(phi)
        xymax = [[0] * (len(phi)), [0] * (len(phi))]
        xymin = [[0] * (len(phi)), [0] * (len(phi))]
        for i in range(len(phi)):
            xymax[0][i] = rmax[i] * math.cos(phi[i]) / 2  # Divisão por 2 para colocar em unidades de rg, não de M.
            xymax[1][i] = rmax[i] * math.sin(phi[i]) / 2

        #PLOT ORBITA 1

        plt.figure()
        plt.subplot(1, 1, 1)
        plt.xlabel("x / rg")
        plt.ylabel("y / rg")
        circle = Circle((0, 0), 2 * 1.477, color='black')
        plt.gca().add_patch(circle)
        plt.gca().set_aspect('equal')
        plt.plot(xymax[0], xymax[1], 'k--', color='black')
        plt.axis([1 / umax + 1, -1 / umax - 1, 1 / umax + 1, -1 / umax - 1])

    else:
        if k < wmax and ust < umax:
            uint = ust
            uext = tp2 * (1 - eps)
            norbit = 1
        elif k > wmax:
            uint = ust
            uext = 0.5 * (1 - eps)
            norbit = 0.5
        elif k < wmax and ust > umax:
            uint = 0.5
            uext = tp3 * (1 + eps)
            norbit = 0.5
        else:
            display("Há uma incoerência entre os parametros fornecidos", target="graph3", append=True)

        v = sp.Symbol('v')


        def theta(v):
            theta = (k - w(v)) ** (-1 / 2)
            return theta


        delphi, erro = quad(theta, uint, uext)

        n = 500
        uc = np.arange(uint, uext, (uext - uint) / n)
        ud = np.arange(uext, uint, (uint - uext) / n)

        phi1 = []
        for i in range(len(uc)):
            a = quad(theta, uint, uc[i])
            phi1.append(abs(a[0]))

        phi2 = []
        for j in range(len(ud)):
            b = quad(theta, uext, ud[j])
            phi2.append(abs(b[0]))

        if norbit == 0.5:
            utotal = uc
        else:
            utotal = utotal = np.concatenate([uc, ud] * (norbit))

        accphi = [0] * (len(utotal))

        if norbit == 0.5:
            accphi = phi1
            x = [0] * (len(uc))
            y = [0] * (len(uc))
            for i in range(len(uc)):
                x[i] = (math.cos(accphi[i])) / utotal[
                    i] / 2  # Divisão por 2 para colocar em unidades de rg, não de M.
                y[i] = (math.sin(accphi[i])) / utotal[i] / 2
        else:
            for i in range(norbit):
                for j in range(n):
                    accphi[j + (2 * i * n)] = 2 * i * delphi + phi1[j]
                    accphi[j + ((2 * i + 1) * n)] = ((2 * i) + 1) * delphi + phi2[j]
            x = [0] * (2 * norbit * n)
            y = [0] * (2 * norbit * n)
            for i in range(2 * norbit * n):
                x[i] = (math.cos(accphi[i])) / utotal[i] / 2
                y[i] = (math.sin(accphi[i])) / utotal[
                    i] / 2  # Divisão por 2 para colocar em unidades de rg, não de M.

        #PLOT ORBITA 2

        fig2 = plt.figure()
        plt.subplot(1, 1, 1)
        plt.plot(x, y, 'k--', color='gold')
        plt.xlabel("x [rg]")
        plt.ylabel("y [rg]")
        circle = Circle((0, 0), 1, color='dimgrey')
        plt.gca().add_patch(circle)
        plt.gca().set_aspect('equal')
        plt.axis([(-1 / uint + 10) / 2, (1 / uint - 10) / 2, (-1 / uint + 10) / 2, (1 / uint - 10) / 2])
        ax = plt.gca()
        ax.spines['bottom'].set_color(eixos)
        ax.tick_params(axis='x', colors=eixos)
        ax.tick_params(axis='y', colors=eixos)
        ax.spines['top'].set_color(eixos)
        ax.spines['right'].set_color(eixos)
        ax.spines['left'].set_color(eixos)
        ax.xaxis.label.set_color(eixos)
        ax.yaxis.label.set_color(eixos)
        fig2.patch.set_facecolor(borda)
        ax.set_facecolor("black")
        display(plt, target="graph3", append=True)
