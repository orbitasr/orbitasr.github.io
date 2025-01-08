import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
import math
from matplotlib.patches import Circle

global firstRun
firstRun = True

def reiniciar():
    global orbitas
    orbitas = False
    js.orbitas = orbitas
    global firstRun
    firstRun = True

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

def gerarPotencial():
    def v(u, l):
        v = -u + (l ** 2) * (u ** 2) / 2 - (l ** 2) * (u ** 3)
        return v

    # Muda o escopo das variáveis abaixo para global:
    global l, vmin, vmax, eps, r_escala, rg
    # Definições básicas
    rg = 2.953 # Raio de Schwarzachild para o Sol (2*G*M_sun/c^2), em km. Pode ser alterado caso se queira mudar a massa do buraco negro.
    r_escala = rg/2 # Fator de recuperação de unidades.
    eps = 1e-8
    
    eixos, borda = mudarCor()
    
    #Input: valor do momento angular adimensional (l>0)
    #linput = input("Insira o valor do momento angular adimensional: ")
    from js import m
    linput = m
    if isNumeric(linput):
        l = float(linput)
    else:
        l = -1 # Se o input não é numérico, atribui um valor negativo a l para gerar mensagem de erro abaixo.

    #Output: GRÁFICO DO POTENCIAL. Distingue os casos em que o potencial possui máximo/mínimo local ou não.
    if l > np.sqrt(12):
        coef = [- 3 * (l ** 2), l ** 2, -1]
        a = np.roots(coef)
        umax = sorted(a)[1].real
        umin = sorted(a)[0].real
        vmin = v(umin, l)
        vmax = v(umax, l)
        #print("O mínimo e o máximo da energia potencial efetiva são", vmin, "e", vmax,
        #         " (ver pontos no gráfico).")
        t = 'O mínimo e o máximo da energia potencial efetiva são ' + str(round(vmin,8)+eps) + ' e ' + str(round(vmax,8)) + ' (ver pontos no gráfico).'
        display(t, target="valores", append=True)
        
        vlim = vmax
        rmax = 2 / umin
        r = np.arange(2, rmax, rmax / 30000)
        u = 1 / r

        fig1 = plt.figure()
        plt.subplot(1, 1, 1)
        plt.axhline(0, linewidth=0.3, color='white')
        plt.plot(r_escala * r, v(u, l), color="white")
        plt.plot([r_escala / umin, r_escala / umax], [vmin, vmax], 'bo', color="gold")
        plt.xlabel("r [km]")
        plt.ylabel(r"$\bar{U}_{efetiva}^{(RG)}$")
        plt.axis([0, r_escala * rmax, -0.5, vlim + 0.1])
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
        #plt.show()
        display(plt, target="graph", append=True)
        
    elif l>=0 and l <= np.sqrt(12):
        vlim = 0
        rmax = 30
        #print("Para esse valor de momento angular, a energia potencial efetiva não possui um mínimo ou máximo local.")
        t = "Para esse valor de momento angular, a energia potencial efetiva não possui um mínimo ou máximo local."
        display(t, target="valores", append=True)
        r = np.arange(2, rmax, rmax / 1000)
        u = 1 / r

        fig1 = plt.figure()
        plt.subplot(1, 1, 1)
        plt.plot(r_escala * r, v(u, l), color="white")
        plt.axhline(0, linewidth=0.3, color='white')
        plt.xlabel("r [km]")
        plt.axis([0, r_escala * rmax, -0.5, vlim + 0.1])
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
        #plt.show()
        display(plt, target="graph", append=True)
        
       
def gerarOrbitaM(): 

        global circular
        circular = False

        def theta(w, l, E):
            theta = (l / (2 ** (1 / 2))) / ((E - v(w, l)) ** (1 / 2))
            return theta

        def v(u, l):
            v = -u + (l ** 2) * (u ** 2) / 2 - (l ** 2) * (u ** 3)
            return v

        global mensagem
        mensagem = False
        js.mensagem = mensagem
        
        eixos, borda = mudarCor()

        #Input: valor do parâmetro de energia (E>-0.5)
        #Einput = input("Insira o valor do parâmetro de energia: ")
        from js import e
        Einput = e
        if isNumeric(Einput):
            E = float(Einput)
        else:
            E = -1 # Se o input não é numérico, atribui um valor menor que -1/2 a E para gerar mensagem de erro abaixo.
        
        rst = 60 # Raio de início da órbita no caso de órbita de espalhamento ou captura.
        ust = 1/rst
        Emin = v(1/2,l) # Valor mínimo permitido para o parâmetro de energia, igual a -0.5.
        
        # Define os pontos de retorno (tps)
        if l == 0:
            tps = [-E]
        else:
            coef = [- (l ** 2), (l ** 2) / 2, -1, -E]
            roots = np.roots(coef)           
            positiverealroots = [x.real for x in roots if (x.imag == 0 and x.real > 0)]
            tps = sorted(positiverealroots)
        
        # Define os pontos u1 e u2 (os pontos de u = 1/r máximo e mínimo da órbita).
        if l > math.sqrt(12):
            if E == vmin or (0 < E - vmin < eps):
                E = vmin
                rcirc = l**2 / 2 * (1 - math.sqrt(1 - 12/l**2))
                circular = True
            if E == vmax:
                E = E + 10*eps #Força uma órbita de captura.
            if E < 0 and E > vmin and E < vmax:
                u1 = tps[0] * (1 + eps)
                u2 = tps[1] * (1 - eps)
                rlim = 1 / u1 + 2 # Limite para os gráficos
                #INPUT: Número de voltas para órbita ligada
                #Faz aparecer input no site:
                global orbitas
                orbitas = True
                js.orbitas = orbitas
                
                global ocultar
                ocultar = True
                js.ocultar = ocultar

                global firstRun
                if firstRun:
                    firstRun = False
                    return

                ocultar = False
                js.ocultar = ocultar

                #norbitinput = input("Para essa escolha de parâmetros, dois tipos de órbitas são possíveis. A órbita ligada será mostrada. Escolha o número de voltas que deseja traçar: ")
                from js import o
                norbitinput = o
                if isNumeric(norbitinput) and eval(norbitinput)>0:
                    norbit = int(eval(norbitinput))
                else:
                    #print("ATENÇÃO: Valor inválido para o número de voltas. Mostrando 5 voltas.")
                    display("ATENÇÃO: Valor inválido para o número de voltas. Mostrando 5 voltas.", target="infos", append=True)
                    mensagem = True
                    js.mensagem = mensagem
                    norbit = 5
            elif (E < 0 and E > vmax) or E < vmin: #Órbita de captura começando do ponto de retorno.
                u1 = tps[0] * (1 + eps)
                u2 = 0.5
                rlim = 1 / u1 + 1
                norbit = 0.5
            elif 0 <= E < vmax:
                u1 = ust
                u2 = tps[0] * (1 - eps)
                rlim = 1 / u1 - 10
                norbit = 1
                display("Para essa escolha de parâmetros, dois tipos de órbitas são possíveis. A órbita de espalhamento está sendo mostrada.", target="infos", append=True)
                mensagem = True
                js.mensagem = mensagem
            elif E >=0 and E > vmax:
                u1 = ust
                u2 = 0.5
                rlim = 1 / u1 - 5
                norbit = 0.5
                
        else: # l < sqrt(12)
            if E >= 0:
                u1 = ust
                u2 = 0.5
                rlim = 1 / u1 - 5
                norbit = 0.5
            elif Emin < E < 0:
                u1 = tps[0] * (1 + eps)
                u2 = 0.5
                rlim = 1 / u1 + 1
                norbit = 0.5

        # Calcula a órbita se E>Emin.
        if E <= Emin:
            display("ATENÇÃO: O parâmetro de energia deve ser maior que Emin = -0.5 para garantir que a órbita tenha início fora do horizonte de eventos.", target="infos", append=True)
            mensagem = True
            js.mensagem = mensagem
            
        else:          
            if not circular:
                delphi, erro = quad(theta, u1, u2, args=(l, E))
                
                if 1/u1 > 1000:
                    #print("ATENÇÃO: Órbitas muito excêntricas podem não ser representadas adequadamente.")
                    display("ATENÇÃO: Órbitas muito excêntricas podem não ser representadas adequadamente.", target="infos", append=True)
                    mensagem = True
                    js.mensagem = mensagem
                if 1.1*u1 < 0.9*u2:
                    uc1 = np.linspace(u1, 1.1*u1, 300)
                    uc2 = np.linspace(1.11*u1, 0.9*u2, 1000)
                    uc3 = np.linspace(0.91*u2, u2, 300)
                    uc = np.concatenate((uc1,uc2,uc3))
                else:
                    uc = np.linspace(u1,u2,500)
                ud = uc[::-1]
                n = len(uc)
    
                phi1 = []
                for i in range(n):
                    a = quad(theta, u1, uc[i], args=(l, E))
                    phi1.append(abs(a[0]))
    
                phi2 = []
                for j in range(n):
                    b = quad(theta, u2, ud[j], args=(l, E))
                    phi2.append(abs(b[0]))
    
                if norbit == 0.5:
                    utotal = uc
                else:
                    utotal = np.concatenate([uc, ud] * (norbit))
    
                accphi = [0] * (len(utotal))
    
                if norbit == 0.5:
                    accphi = phi1
                    x = [0] * (n)
                    y = [0] * (n)
                    for i in range(n):
                        x[i] = (math.cos(accphi[i])) / utotal[i] * r_escala
                        y[i] = (math.sin(accphi[i])) / utotal[i] * r_escala
                else:
                    for i in range(norbit):
                        for j in range(n):
                            accphi[j + (2 * i * n)] = 2 * i * delphi + phi1[j]
                            accphi[j + ((2 * i + 1) * n)] = ((2 * i) + 1) * delphi + phi2[j]
                    x = [0] * (2 * norbit * n)
                    y = [0] * (2 * norbit * n)
                    for i in range(2 * norbit * n):
                        x[i] = ((math.cos(accphi[i])) / utotal[i]) * r_escala
                        y[i] = ((math.sin(accphi[i])) / utotal[i]) * r_escala
            
            else: #Caso circular
                n = 500
                x = [0] * (n)
                y = [0] * (n)
                for i in range(n):
                    x[i] = rcirc * r_escala * math.cos(2 * math.pi * i / n)
                    y[i] = rcirc * r_escala * math.sin(2 * math.pi * i / n)
                rlim = rcirc + 1
                    
            #GRÁFICO DA ÓRBITA
            fig2 = plt.figure()
            plt.plot(x, y, color="gold")
            plt.xlabel("x [km]")
            plt.ylabel("y [km]")
            circle = Circle((0, 0), 2 * r_escala, color='dimgrey')
            plt.gca().add_patch(circle)
            plt.gca().set_aspect('equal')
            plt.axis(
                [- rlim * r_escala, rlim * r_escala, - rlim * r_escala, rlim * r_escala])
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
            #plt.show()
            display(plt, target="graph2", append=True)
