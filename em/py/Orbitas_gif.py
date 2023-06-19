import numpy as np
import matplotlib.pyplot as plt
#import sympy as sp
from scipy.integrate import quad, solve_ivp
import math
#from PIL import Image
import matplotlib.animation as animation
from matplotlib.patches import Circle

#import warnings
#warnings.filterwarnings('ignore')

#from matplotlib.animation import FuncAnimation, writers
# matplotlib.rc('animation', html='html5')
# from IPython.display import HTML

def isNumeric(input): # Checa se a natureza do input é numérica.
    try:
        float(input)
        return True
    except ValueError:
        return False
    
def closest(lst, value): # Encontra o valor mais próximo a "value" na lista "lst".
    lst = np.asarray(lst)
    idx = (np.abs(lst - value)).argmin()
    return idx, lst[idx]

# Definições básicas
rg = 2.953 # Raio de Schwarzachild para o Sol (2*G*M_sun/c^2), em km. Pode ser alterado caso se queira mudar a massa do buraco negro.
r_escala = rg/2 # Fator de recuperação de unidades.
eps = 1e-8

test_mode = False # Se True, gráfico do potencial é mostrado, a fim de testar consistência dos resultados. Se False, só o gif é exibido.

# Input: trajetória tipo-tempo ('M') ou tipo-luz ('L')
tipo_orbita = input("Escolha 'M' para órbita de corpos celestes e 'L' órbitas de raios de luz: ")

if tipo_orbita == "M":
    
    # Definindo funções importantes

    def v(u, l):
        v = -u + (l ** 2) * (u ** 2) / 2 - (l ** 2) * (u ** 3)
        return v    
    
    def theta(w, l, E):
        theta = l * (2.0 * (E - v(w, l))) ** (-1 / 2)
        return theta
    
    def tau_integrand(w, l, E):
        tau_integrand = w ** (-2) * (2.0 * (E - v(w, l))) ** (-1 / 2)
        return tau_integrand
    
    # Input: posição inicial
    x0input = input("Escolha o valor da posição inicial (em km): ")
    if isNumeric(x0input):
        x0 = float(x0input)
    else:
        x0 = 0 # Se o input não é numérico, atribui um valor nulo para gerar mensagem de erro abaixo.

    if x0 <= rg:
        print("ATENÇÃO: Escolha um valor da posição inicial maior que 3 km, caso contrário, a trajetória já começa dentro do buraco negro.")    
    else:
        # Input: módulo da velocidade inicial
        v0input = input("Escolha o valor da velocidade inicial (em unidades da velocidade da luz): ")
 
        if isNumeric(v0input) and eval(v0input) < 1:
            v0 = abs(float(v0input)) #Força um valor positivo.
            
            rst = x0 / r_escala
            ust = 1 / rst
        
            l = rst * v0 # Momento angular por unidade de massa, adimensional.
            E = v(ust, l) + eps  # Hipótese: dr/dt = 0 inicialmente.
            
            if test_mode == True:      
                print("Energia: ",E,". Momento angular: ", l)
                r = np.arange(2, 10, 1 / 10000)
                u = 1 / r
                fig1 = plt.figure()
                plt.subplot(1, 1, 1)
                plt.axhline(0, linewidth=0.3, color='white')
                plt.plot(r_escala * r, v(u, l), color="white")
                plt.plot(np.linspace(0, max(r) * r_escala, 50), np.ones(50) * E, color="green")
                plt.xlabel("r [km]")
                plt.axis([0, r_escala * 10, -0.5,  1.2*max(v(u, l))])
                ax = plt.gca()
                ax.spines['bottom'].set_color('white')
                ax.tick_params(axis='x', colors='white')
                ax.tick_params(axis='y', colors='white')
                ax.spines['top'].set_color('white')
                ax.spines['right'].set_color('white')
                ax.spines['left'].set_color('white')
                ax.xaxis.label.set_color('white')
                ax.yaxis.label.set_color('white')
                fig1.patch.set_facecolor('#0E1117')
                ax.set_facecolor("black")
                plt.show()
                
            # Calcula os pontos de retorno.
            if l == 0:
                tps = [-E]
            else:
                coef = [- (l ** 2), (l ** 2) / 2, -1, -E]
                roots = np.roots(coef)           
                positiverealroots = [x.real for x in roots if (x.imag == 0 and x.real > 0)]
                tps = sorted(positiverealroots)
        
            # Define os pontos u1 e u2 (os pontos de u = 1/r máximo e mínimo da órbita).
            if l > math.sqrt(12):                
                
                ist = closest(tps, ust)[0]
                if E >= 0:
                    if ist == 0:
                        print("Para essa escolha de parâmetros, a partícula escapa do buraco negro.") # Órbita de espalhamento
                        u1 = ust * (1 - eps)
                        u2 = ust / 10
                        norbit = 0.5
                    elif ist == 1:
                        print("Para essa escolha de parâmetros, a partícula cai no buraco negro.") # Órbita de captura
                        u1 = ust * (1 + eps)
                        u2 = 0.55
                        norbit = 0.5
                    else:
                        print("ERRO 1")
                else: # E < 0
                    if len(tps) == 3:
                        if ist == 0:
                            u1 = ust * (1 + eps)
                            u2 = tps[1] * (1 - eps)
                        elif ist == 1:
                            u1 = ust * (1 - eps)
                            u2 = tps[0] * (1 + eps)
                        elif ist == 2:
                            print("Para essa escolha de parâmetros, a partícula cai no buraco negro.") # Órbita de captura
                            u1 = ust * (1 + eps)
                            u2 = 0.55
                            norbit = 0.5
                        
                        #INPUT: Número de voltas para órbita ligada
                        if ist == 0 or ist == 1:
                            norbitinput = input("Para essa escolha de parâmetros, a partícula orbita o buraco negro. Escolha o número de voltas que deseja traçar: ") # Órbita ligada
                            if isNumeric(norbitinput) and eval(norbitinput)>0:
                                norbit = int(eval(norbitinput))
                            else:
                                print("ATENÇÃO: Valor inválido para o número de voltas. Mostrando 5 voltas.")
                                norbit = 5                            
                        
                    elif len(tps) == 1:
                        print("Para essa escolha de parâmetros, a partícula cai no buraco negro.") # Órbita de captura
                        u1 = ust * (1 + eps)
                        u2 = 0.55
                        norbit = 0.5
                    
                    else:
                        print("ERRO 2")
            else: # l < sqrt(12)
                print("Para essa escolha de parâmetros, a partícula cai no buraco negro.") # Órbita de captura
                u1 = ust * (1 + eps)
                u2 = 0.55
                norbit = 0.5
            
            umin = min(u1, u2)
            umax = max(u1, u2)
            if 1/umin > 1000:
                print("ATENÇÃO: Órbitas muito excêntricas podem não ser representadas adequadamente.")

            n = 500
            Ttotal, erroT = quad(tau_integrand, umin, umax, args=(l, E))  # Tempo total para ir de u1 para u2.
            teval = np.linspace(0, Ttotal, n)
            dt = teval[1] - teval[0]
            
            f = lambda t, u: u ** 2 * (2.0 * abs(E - v(u, l))) ** (1 / 2.0)
            sol = solve_ivp(f, [0, Ttotal], [umin], t_eval = teval)
            uc = np.append(sol.y[0][:-1] , [umax])
            ud = uc[::-1]
            
            delphi, erro = quad(theta, umin, umax, args=(l, E))
            
            if umin == u1:
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

            else: # umin == u2, 
                phi1 = []
                for j in range(n):
                    a = quad(theta, u1, ud[j], args=(l, E))
                    phi1.append(abs(a[0]))
                    
                phi2 = []
                for i in range(n):
                    b = quad(theta, u2, uc[i], args=(l, E))
                    phi2.append(abs(b[0]))

                if norbit == 0.5:
                    utotal = ud
                else:
                    utotal = np.concatenate([ud, uc] * (norbit))                

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
                    
  
            fig = plt.figure()
            
            plt.xlabel("x (km)")
            plt.ylabel("y (km)")
            plt.gca().set_aspect('equal')
            ax = plt.gca()
            ax.spines['bottom'].set_color('white')
            ax.tick_params(axis='x', colors='white')
            ax.tick_params(axis='y', colors='white')
            ax.spines['top'].set_color('white')
            ax.spines['right'].set_color('white')
            ax.spines['left'].set_color('white')
            ax.xaxis.label.set_color('white')
            ax.yaxis.label.set_color('white')
            fig.patch.set_facecolor('#0E1117')
            ax.set_facecolor("black")
            circle = Circle((0, 0), rg, color='dimgrey', linewidth=0)
            plt.gca().add_patch(circle)
            if u1 < u2 or norbit > 0.5:
                plt.axis([- r_escala * 1.1 / umin, r_escala * 1.1/ umin, - r_escala * 1.1/ umin, r_escala * 1.1/ umin])
            else:
                plt.axis([- r_escala / umin, r_escala / umin, - r_escala / umin, r_escala / umin])
            
            # Montagem do gif
            
            graph, = plt.plot([], [], color="gold", markersize=3, label='Tempo: 0 s')
            #L = plt.legend(loc=1) #Legenda
            
    #        plt.close()  # Não mostra a imagem de fundo [Raissa: Comentei porque com isso o gif não aparecia.]
                        
            def animate(i):
                 #lab = 'Tempo: ' + str(round(teval[i] * r_escala / (3e2), -int(math.floor(math.log10(dt * r_escala /(3e2)))))) + ' ms'
                 graph.set_data(x[:i], y[:i])
                 #L.get_texts()[0].set_text(lab)  # Atualiza a legenda a cada frame
                 
                 return graph,
            
            skipframes = int(len(uc) / 150)
            if skipframes == 0:
                skipframes = 1
            
     #       for i in range(len[x]):
      #          animate(i)
         #       time.sleep(0.01)
            #
            
            ani1 = animation.FuncAnimation(fig, animate, frames=range(0, len(x), skipframes), interval=30, blit=True, repeat=False)
            
            plt.show()
        
        else:
            print("ATENÇÃO: Valor inválido para a velocidade inicial, que deve ser um número menor que 1.") 

    # HTML(ani1.to_jshtml())
    # components.html(ani1.to_jshtml(),height=800)


if tipo_orbita == "L":
        # Definindo funções importantes

    def w(u):
        w = u ** 2 - 2 * (u ** 3)
        return w
    
    def lambda_integrand(v, k):
        lambda_integrand = 1 / (v ** 2) * (k - w(v)) ** (-1 / 2.0)
        return lambda_integrand
    
    def theta(v, k):
        theta = (k - w(v)) ** (-1 / 2.0)
        return theta
    
    umax = 1 / 3
    wmax = w(umax) #1 / 27

    # Input: parâmetro de impacto
    dinput = input("Escolha o valor do parâmetro de impacto (em km): ")
    if isNumeric(dinput):
        d = float(dinput)
    else:
        d = -1 # Se o input não é numérico, atribui um valor nulo para gerar mensagem de erro abaixo.

    if d >= 0:
        param_imp = d / r_escala
        k = 1 / param_imp**2
        
        rst = 50 # Esse valor pode ser alterado, mas deve ser maior que 3.
        ust = 1 / rst

        coef = [-2, 1, 0, -k]
        roots = np.roots(coef)
        positiverealroots = [x.real for x in roots if (x.imag == 0 and x.real > 0)]
        tps = sorted(positiverealroots)
        
        if k < wmax:
            print("Para essa escolha de parâmetros, o raio de luz escapa do buraco negro.")
            if ust < tps[0] / 2:
                u1 = ust
            else: 
                u1 = tps[0] / 2
            u2 = tps[0] * (1 - eps)
            norbit = 1
        else:
            print("Para essa escolha de parâmetros, o raio de luz é capturado pelo buraco negro.")
            u1 = ust
            u2 = 0.5 
            norbit = 0.5    
    
        npoints = 400
        lambdatotal, errolambda = quad(lambda_integrand, u1, u2, args=(k))  # Calcula parâmetro afim total para ir de u1 a u2     
        dlambda = lambdatotal/npoints

#        lambda_eval = np.linspace(0, lambdatotal, npoints)        
#        f = lambda t, v: (v ** 2) * (k - w(v)) ** (1 / 2.0)
 #       sol = solve_ivp(f, [0, lambdatotal], [u1], t_eval = lambda_eval)
  #      uc = np.append(sol.y[0][:-1] , [u2])
   #     ud = uc[::-1]
        
        ud = [u2]
        for i in range(npoints + 20):
            ud.append(ud[i] - dlambda * ud[i] ** 2 * (k - w(ud[i])) ** (1 / 2.0))
            if ud[-1].imag != 0 or math.isnan(ud[-1]):
                ud = ud[:-2]
                break
        uc = ud[::-1]
        n = len(uc)

        delphi, erro = quad(theta, u1, u2, args=(k))
    
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
        phi0 = np.arcsin(param_imp / rst)
    
        if norbit == 0.5:
            accphi = phi1
            x = [0] * (len(uc))
            y = [0] * (len(uc))
            for i in range(len(uc)):
                x[i] = (math.cos(phi0 + accphi[i])) / utotal[i] * r_escala
                y[i] = (math.sin(phi0 + accphi[i])) / utotal[i] * r_escala
        else:
            for i in range(norbit):
                for j in range(n):
                    accphi[j + (2 * i * n)] = 2 * i * delphi + phi1[j]
                    accphi[j + ((2 * i + 1) * n)] = ((2 * i) + 1) * delphi + phi2[j]
            x = [0] * (2 * norbit * n)
            y = [0] * (2 * norbit * n)
            for i in range(2 * norbit * n):
                x[i] = (math.cos(phi0 + accphi[i])) / utotal[i] * r_escala
                y[i] = (math.sin(phi0 + accphi[i])) / utotal[i] * r_escala
    
        fig = plt.figure()
    
        plt.xlabel("x (km)")
        plt.ylabel("y (km)")
        plt.gca().set_aspect('equal')
        ax = plt.gca()
        ax.spines['bottom'].set_color('white')
        ax.tick_params(axis='x', colors='white')
        ax.tick_params(axis='y', colors='white')
        ax.spines['top'].set_color('white')
        ax.spines['right'].set_color('white')
        ax.spines['left'].set_color('white')
        ax.xaxis.label.set_color('white')
        ax.yaxis.label.set_color('white')
        fig.patch.set_facecolor('#0E1117')
        ax.set_facecolor("black")
        circle = Circle((0, 0), rg, color='dimgrey')
        plt.gca().add_patch(circle)
        plt.axis([- r_escala * 0.85 / u1, r_escala * 0.85 / u1, - r_escala * 0.85 / u1, r_escala * 0.85 / u1])
    
        # Montagem do gif
        graph, = plt.plot([], [], 'k--', color="gold", markersize=3)
     #   plt.close()  # Não mostra a imagem de fundo
    
        def animate(i):
            graph.set_data(x[:i], y[:i])
            return graph,
    
        skipframes = int(len(uc) / 200)
        if skipframes == 0:
            skipframes = 1
    
        ani2 = animation.FuncAnimation(fig, animate, frames=range(0, len(x), skipframes), interval=10, blit=True, repeat=False)
        plt.show()
        
    else:
        print("ATENÇÃO: Valor inválido para o parâmetro de impacto.")    

