import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
import math
from matplotlib.patches import Circle

#Muda cor do gráfico de acordo com tema do site ~ A
def mudarCor():
    from js import mode
    if mode == 'light':
        eixos = '#5e6469'
        borda = '#f7f7f7'
    else:
        eixos = 'white'
        borda = '#181F2F'
    return eixos, borda

def isNumeric(input): # Função genérica para checar se a natureza do input é numérica
    try:
        float(input)
        return True
    except ValueError:
        return False

#elif tipo_orbita == "L":
def gerarOrbitaL():
    global msg
    msg = False
    # Definições básicas
    r_escala = 1/2 # Todas as expressões estão em unidades de M, mas o input é em unidades do raio gravitacional (2*M).
    eps = 1e-8
    
    def w(u):
        w = u**2 - 2 * (u**3)
        return w
    
    def theta(v, k):
        theta = (k - w(v)) ** (-1 / 2)
        return theta
    
    eixos, borda = mudarCor()
    
    umax = 1 / 3
    wmax = w(1/3) #=1 / 27

    """
    # GRÁFICO DA ENERGIA POTENCIAL EFETIVA
    #print("O potencial efetivo para raios de luz tem um máximo local de 0.037037037... em r = 1.5rg.")
    display("O potencial efetivo para raios de luz tem um máximo local de 0.037037037... em r = 1.5rg.", target="info2", append=True)
    r = np.linspace(0.1, 10, 100)
    u = 1/r
    fig3 = plt.figure()
    plt.subplot(1, 1, 1)
    plt.plot(r_escala * r, w(u), color="white")
    plt.plot([r_escala * 3], [w(1/3)], 'bo', color="gold")
    plt.axhline(0, linewidth=0.3, color='white')
    plt.xlabel("r [rg]")
    plt.axis([0, r_escala * 10, -0.05, 0.05])
    ax = plt.gca()
    ax.spines['bottom'].set_color(eixos)
    ax.tick_params(axis='x', colors=eixos)
    ax.tick_params(axis='y', colors=eixos)
    ax.spines['top'].set_color(eixos)
    ax.spines['right'].set_color(eixos)
    ax.spines['left'].set_color(eixos)
    ax.xaxis.label.set_color(eixos)
    ax.yaxis.label.set_color(eixos)
    fig3.patch.set_facecolor(borda)
    ax.set_facecolor("black")
    #plt.show()
    display(plt, target="graph3", append=True)
    """
    
    #INPUT (opção 1): Parâmetro de impacto
   # dinput = input("Escolha o valor do parâmetro de impacto (d), em unidades de $r_g$: ")
   # if isNumeric(dinput) and eval(dinput)>=0:
   #     d = eval(dinput)
   #     b = d / r_escala
   #     k = 1 / (b ** 2)

    #INPUT (opção 2): Parâmetro de energia
    from js import d
    kinput = d
    # kinput = input("Escolha o valor do parâmetro de energia (k): ")
    if isNumeric(kinput) and eval(kinput) > 0:
        k = eval(kinput) / 4
        
        rst = 50
        ust = 1 / rst
        r = np.linspace(0.1, rst, 100)
        u = 1 /r 
        
        coef = [-2, 1, 0, -k]
        roots = np.roots(coef)
        positiverealroots = [x.real for x in roots if (x.imag == 0 and x.real > 0)]
        tps = sorted(positiverealroots)
        
        if k == wmax:
            k = k + eps
        if k < wmax:
            if ust < tps[0] / 2:
                u1 = ust
            else: 
                u1 = tps[0] / 2
            u2 = tps[0] * (1 - eps)
            rlim = 1 / u1 - 20
            norbit = 1
            display("Para essa escolha de parâmetros, dois tipos de órbitas são possíveis. A órbita de espalhamento será mostrada.", target="info3", append=True)
            msg = True
        else:
            u1 = ust
            u2 = 0.5 # * (1 - eps)
            rlim = 1 / u1 - 20
            norbit = 0.5
        
        delphi, erro = quad(theta, u1, u2, args=(k))
    
        uc1 = np.linspace(u1, 0.96*u2, 1000)
        uc2 = np.linspace(0.961*u2, u2, 2000)
        uc = np.concatenate((uc1,uc2))
        ud = uc[::-1]
        n = len(uc)
    
        phi1 = []
        for i in range(n):
            a = quad(theta, u1, uc[i], args=(k))
            phi1.append(abs(a[0]))

        phi2 = []
        for j in range(n):
            b = quad(theta, u2, ud[j], args=(k))
            phi2.append(abs(b[0]))

        if norbit == 0.5:
            utotal = uc
        else:
            utotal = utotal = np.concatenate([uc, ud] * (norbit))

        accphi = [0] * (len(utotal))

        if norbit == 0.5:
            accphi = phi1
            x = [0] * (n)
            y = [0] * (n)
            for i in range(n):
                x[i] = (math.cos(accphi[i])) / utotal[i] * r_escala
                y[i] = (math.sin(accphi[i])) / utotal[i] * r_escala
        elif norbit == 1:
            for j in range(n):
                accphi[j] = phi1[j]
                accphi[j + n] = delphi + phi2[j]
            x = [0] * (2 * n)
            y = [0] * (2 * n)
            for i in range(2 * n):
                x[i] = (math.cos(accphi[i])) / utotal[i] * r_escala
                y[i] = (math.sin(accphi[i])) / utotal[i] * r_escala
        
        #Output: GRÁFICO DA ÓRBITA

        fig4 = plt.figure()
        plt.subplot(1, 1, 1)
        plt.plot(x, y, 'k--', color='gold')
        plt.xlabel("x [rs]")
        plt.ylabel("y [rs]")
        circle = Circle((0, 0), 1, color='dimgrey')
        plt.gca().add_patch(circle)
        plt.gca().set_aspect('equal')
        plt.axis([- rlim * r_escala, rlim * r_escala, -rlim * r_escala, rlim * r_escala])
        ax = plt.gca()
        ax.spines['bottom'].set_color(eixos)
        ax.tick_params(axis='x', colors=eixos)
        ax.tick_params(axis='y', colors=eixos)
        ax.spines['top'].set_color(eixos)
        ax.spines['right'].set_color(eixos)
        ax.spines['left'].set_color(eixos)
        ax.xaxis.label.set_color(eixos)
        ax.yaxis.label.set_color(eixos)
        fig4.patch.set_facecolor(borda)
        ax.set_facecolor("black")
        #plt.show()
        display(plt, target="graph4", append=True)
        
    else:
        #print("ATENÇÃO: O parâmetro de impacto deve ser um número positivo.")
        #print("ATENÇÃO: O parâmetro de energia deve ser um número positivo (>1e-8) para garantir que a órbita tenha início fora do horizonte de eventos.")
        display("ATENÇÃO: O parâmetro de energia deve ser um número positivo (>1e-8) para garantir que a órbita tenha início fora do horizonte de eventos.", target="info3", append=True)
        msg = True

#else:
#    print("ATENÇÃO: Apenas 'M' e 'L' podem ser usados para identificar os tipos de órbitas.")
