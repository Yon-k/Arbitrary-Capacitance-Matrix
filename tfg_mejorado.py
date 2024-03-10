import numpy as np
import random
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import copy
from mpl_toolkits.mplot3d import Axes3D
import warnings
import matplotlib.cbook
from primeras_vecinas import *
warnings.filterwarnings("ignore",category=matplotlib.cbook.mplDeprecation)


#-------------------------------META-FUNCTIONS----------------------------------

def escoge_radio(ratio, mu1, sigma1, mu2, sigma2):
    cutoff = random.uniform(0,1)
    if cutoff > ratio:
        R = np.random.lognormal(mu1,sigma1)
    else:
        R = np.random.lognormal(mu2,sigma2)
    return R

#-----------------

def calcula_Vfake(Vcaja, coord):
    column = 2
    sumaR = sum(row[column]**2 for row in coord)                        #Lo guardo y recalculo la 3a columna de nuevo
    Vfake = (np.pi*sumaR)/Vcaja
    return Vfake

#-----------------

def dist_centros(vector1, vector2):
    return np.sqrt((vector1[0]-vector2[0])**2 + (vector1[1]-vector2[1])**2)

#-----------------

def dist_surface(vector1, vector2):
    d = np.sqrt((vector1[0]-vector2[0])**2 + (vector1[1]-vector2[1])**2)
    return d - vector1[2] - vector2[2]

#----------------------------------PLOTS----------------------------------------

def plot(coord):
    """Me gustaria pasarle un particula.index y que me recorra el for de abajo a
    partir de ese indice.
    """
    xs = []
    ys = []
    rs = []
    area = []
    for i in range(len(coord)):
        xs.append(coord[i][0])
        ys.append(coord[i][1])
        rs.append(coord[i][2])
        area.append(np.pi*(coord[i][2]**2))

    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i,particula_i in enumerate(coord):
        ax.text(particula_i[0], particula_i[1], i+1)
    #plt.axis([-7, 7, -10, 45])

    plt.grid(linestyle='--')
    plt.scatter(xs, ys, s=area, marker='o',c='blue', edgecolors='r')#facecolors='none')
    plt.xlabel('x (nm)')
    plt.ylabel('y (nm)')
    plt.savefig('Plot1')
    plt.show()

#-----------------

def plot_2(coord):
    distS = []
    distS_2 = []
    for i in range(len(coord)):
        for j in range(len(coord)):
            if i!=j:
                ds = dist_surface(coord[i], coord[j])
                if ds < 5:
                    distS.append(ds)
                if ds < 1:
                    distS_2.append(ds)


    print('Numero de particulas con distancia menor a 1 nm=', len(distS_2))
    print('Numero de particulas con distancia menor a 1 nm=', len(distS))
    n, bins, patches = plt.hist(distS, 32, density=True,color = "skyblue", ec='black')
    sigma = np.std(distS)
    mu = np.mean(distS)
    y = mlab.normpdf( bins, mu, sigma)
    l = plt.plot(bins, y, 'r--', linewidth=2)

    plt.xlabel(r'$d \ (nm)$')
    plt.ylabel(r'$Probability$')
    plt.title(r'Histogram of surface distances')
    plt.text(301, 0.003, r'$\mu=247.7,\ \sigma=124.4$')
    #plt.axis([0, 700, 0, 0.0035])
    #ax.yaxis.set_major_formatter(PercentFormatter(xmax=1))
    plt.grid(linestyle='--')
    #sns.distplot(distS, color="dodgerblue", label="Compact")
    plt.show()
    print('mu_2=',mu,'sigma_2=',sigma)

#-----------------

def plot_3(coord):
    rs = []
    for i in range(len(coord)):
        rs.append(coord[i][2])

    n, bins, patches = plt.hist(rs, 32, density=True,color='green',ec='black')

    plt.xlabel('$R \ (nm)$')
    plt.ylabel(r'$Probability$')
    plt.title('Radial distribution of particles')
    #plt.axis([0, 10, 0, 0.5])
    plt.grid(linestyle='--')
    #ax.yaxis.set_major_formatter(PercentFormatter(xmax=1))
    plt.show()

#-----------------

def plot_4(coord,Lx,Ly):
    xs = []
    ys = []
    rs = []
    area = []
    for i in range(len(coord)):
        xs.append(coord[i][0])
        ys.append(coord[i][1])
        rs.append(coord[i][2])
        area.append(np.pi*(coord[i][2]**2))

    fig= plt.figure()
    ax = fig.add_axes([0.1,0.1,0.8, 0.8])
    ax.grid(True,linestyle='--')
    ax.scatter(xs, ys, s=area, c='blue', marker='o')
    ax.set_xlabel('x (nm)')
    ax.set_ylabel('y (nm)')
    ax.axis([0, Lx, 0, Ly])
    ax.vlines(20,0,Ly,linestyles='dashed',colors='red',label='Left electrode')
    ax.vlines(Lx-20,0,Ly,linestyles='dashed',colors='red',label='Right electrode')
    fig.savefig('Plot1')

    ax2 = fig.add_axes([0.6,0.6,0.25, 0.25])
    n, bins, patches = ax2.hist(rs, 32, density=True,color='green',ec='black')

    ax2.set_xlabel('$R \ (nm)$')
    ax2.set_ylabel(r'$Probability$')
    ax2.set_facecolor('#eafff5')
    ax2.grid(True, linestyle='--')
    vals = ax2.get_yticks()
    ax2.set_yticklabels(['%1.0f%%' %(i*100) for i in vals])

    plt.show()


def plot_5(coord):

    d_tol = 5
    nn = [0]*len(coord)
    for i,particula_i in enumerate(coord):
        for j,particula_j in enumerate(coord):
            d = dist_surface(particula_i, particula_j)
            if d <= d_tol:
                nn[i] += 1
    print(nn)

    n, bins, patches = plt.hist(nn, 32, density=True,color = "skyblue", ec='black')
    plt.xlabel('Nº de primeros vecinos')
    plt.ylabel('Probability')
    #plt.axis([0, 700, 0, 0.0035])
    #ax.yaxis.set_major_formatter(PercentFormatter(xmax=1))
    plt.grid(linestyle='--')
    #sns.distplot(distS, color="dodgerblue", label="Compact")
    plt.show()


    mat = np.matrix(nn)
    with open('primeros_vecinos.txt','w+') as datafile_id:
        for line in mat:
            np.savetxt(datafile_id, line, delimiter = '\n' ,fmt='%.10f')

#-------------------------------------------------------------------------------

def solapa(coord, particula):
    for particula_j in coord[:-1]:
        if dist_surface(particula_j, particula) <= 0:
            return True
    return False

#-----------------

def es_salvable(coord, particula):
    """d tiene que ser mayor que los radios de las dos particulas.
    """
    for particula_j in coord[:-1]:
        d = dist_centros(particula, particula_j)
        if d<particula_j[2] or d<particula[2]:
            return False
    return True

#-----------------

def solapa_0(coord, particula):
    for i,particula_i in enumerate(coord):
        if i==coord.index(particula):
            pass
        elif dist_surface(particula_i, particula) <= 0:
            return True
    return False

#-----------------

def es_salvable_0(coord, particula):
    """d tiene que ser mayor que los radios de las dos particulas.
    """
    for i,particula_i in enumerate(coord):
        d = dist_centros(particula, particula_i)
        if i==coord.index(particula):
            pass
        elif d<particula_i[2] or d<particula[2]:
            return False
    return True

#-----------------

def toca_borde(particula,Lx,Ly):
    """Si toca el borde == True
    """
    condicion_1 = particula[0] > particula[2] and particula[0] < Lx-particula[2]
    condicion_2 = particula[1] > particula[2] and particula[1] < Ly-particula[2]

    if condicion_1 and condicion_2:
        return False
    return True

#-------------------------------------------------------------------------------
#---------------------------------FUNCTIONS-------------------------------------
#-------------------------------------------------------------------------------

def crea_particula(ratio, Lx, Ly, mu1, sigma1, mu2, sigma2, coord):
    """ Me crea una x,y,R siguiendo la distribucion que toca y me lo apendea en coord.
    """
    cutoff = random.uniform(0,1)
    if cutoff > ratio:
        x = random.uniform(0,Lx)
        y = random.uniform(0,Ly)
        R = np.random.lognormal(mu1,sigma1)                                     #Siempre trabajo en nanometros!!!
        coord.append([x,y,R])
    else:
        x = random.uniform(0,Lx)
        y = random.uniform(0,Ly)
        R = np.random.lognormal(mu2,sigma2)
        coord.append([x,y,R])
    return coord

#-----------------

def pre_desolapamiento(Lx, Ly, coord):
    max_iter = 10**5
    max_iter_2 = 10**5
    for i,particula_i in enumerate(coord):
        possibles = []
        for contador in range(max_iter):
            if toca_borde(coord[i],Lx,Ly)==False:
                if solapa_0(coord, coord[i]):
                    if es_salvable_0(coord, coord[i]):
                        possibles.append(copy.deepcopy(coord[i]))
                    coord[i][0] = random.uniform(0,Lx)
                    coord[i][1] = random.uniform(0,Ly)
                else:
                    break
            else:
                coord[i][0] = random.uniform(0,Lx)
                coord[i][1] = random.uniform(0,Ly)
                coord = pre_desolapamiento(Lx, Ly, coord)                       #RECURSIVIDAD POWERRR


        else:                                                                   #Si salgo por el break no se ejecuta el else
            if len(possibles) == 0:
                coord[i][2] =  escoge_radio(ratio, mu1, sigma1, mu2, sigma2)
                return pre_desolapamiento(Lx, Ly, coord)
            else:
                coord[i] = possibles[0]
                for j,particula_j in enumerate(coord):
                    ds = dist_surface(coord[i], coord[j])
                    if ds <= 0:
                        for contador in range(max_iter_2):
                            coord[j][2] = escoge_radio(ratio, mu1, sigma1, mu2, sigma2)
                            ds = dist_surface(coord[i], coord[j])
                            print(ds,'R_i=',coord[i][2],'R_{j}=',coord[j][2])
                            if ds>0:
                                break
                        else:                                                   #Si salgo por el break no se ejecuta el else
                            R = coord[j][2]-abs(ds)
                            R -= random.uniform(0,R)
                            coord[j][2] = R
    return coord

#-----------------

def desolapamiento(Lx, Ly, coord):
    print('reubicando')
    max_iter = 10**5
    max_iter_2 = 10**5
    possibles = []
    for contador in range(max_iter):
        possibles = []
        if solapa(coord, coord[-1]) or toca_borde(coord[-1],Lx,Ly):
            if len(possibles) == 0:
                if es_salvable(coord, coord[-1]):
                    possibles.append(copy.deepcopy(coord[-1]))
            coord[-1][0] = random.uniform(0,Lx)
            coord[-1][1] = random.uniform(0,Ly)
        else:
            break

    else:
        if len(possibles) == 0:
            coord = crea_particula(ratio, Lx, Ly, mu1, sigma1, mu2, sigma2, coord[:-1])
            print("Reintentando reubicar")
            return desolapamiento(Lx, Ly, coord)
        else:
            coord[-1] = possibles[0]
            for i,particula_i in enumerate(coord[:-1]):
                ds = dist_surface(coord[i], coord[-1])
                if ds <= 0:
                    for contador in range(max_iter_2):
                        coord[-1][2] = escoge_radio(ratio, mu1, sigma1, mu2, sigma2)
                        ds = dist_surface(coord[i], coord[-1])
                        print(ds,'R_i=',coord[i][2],'R_{-1}=',coord[-1][2])
                        if ds>0:
                            break
                    else:                                                       #Si salgo por el break no se ejecuta el else
                        R = coord[-1][2]-abs(ds)
                        R -= random.uniform(0,R)
                        coord[-1][2] = R
                        ds = dist_surface(coord[i], coord[-1])
                        print(contador, ds,'R_i=',coord[i][2],'R_{-1}=',coord[-1][2])
    return coord

#-----------------

def corrige_n(calcula_Vfake, crea_particula, desolapamiento, Lx, Ly, n_min, n_max, coord):
    """Hace que el nº de particulas sea el deseado. Esto se consigue haciendo la caja
    mas grande o pequeña segun convenga.
    """
    while len(coord) < n_min or len(coord) > n_max:
        if n_max < len(coord):
            coord.pop(random.randint(0,len(coord)-1))

        elif n_min > len(coord):
            crea_particula(ratio, Lx, Ly, mu1, sigma1, mu2, sigma2, coord)
            coord = desolapamiento(Lx, Ly, coord)

    Vfake = calcula_Vfake(Lx*Ly, coord)
    return coord

#-----------------

def electrodos(coord, Lx):
    """Ojo que empiezo a contar desde el 1 y no el 0.
    """

    lim_izq = 20
    lim_der = Lx-20

    electro_izquierdo = []
    electro_derecho = []
    for i,particula_i in enumerate(coord):
        if particula_i[0] < lim_izq:
            electro_izquierdo.append(i+1)
        elif particula_i[0] > lim_der:
            electro_derecho.append(i+1)

    return electro_izquierdo,electro_derecho

#----------------------------------Main-----------------------------------------
if __name__== "__main__":

    #Parametros del sistema
    L = 350
    Lx = L
    Ly = L

    n_ini = 50
    n_min = 300
    n_max = 1000                                                                #Nº de partículas
    Vreal = 0.20                                                                # % en volumen ocupado que quiero

    Smallpart = 83
    Largepart = 17
    ratio = Largepart/Smallpart

    mu1 = np.log(5.7/2)                                                         #Small part.
    sigma1 = 0.24
    mu2 = np.log(17.0/2)                                                        #Big part.
    sigma2 = 0.20

    tol = 10**(-4)
    lim_inf = Vreal-tol
    lim_sup = Vreal+tol

    seed = 46482619+1                                                           #CAMBIA LA SEED CUANDO CAMBIES DE CONCENTRACION
    np.random.seed(seed)
    random.seed(seed)


    #---------------------------------------------------------------------------
    #---------------------------------------------------------------------------

    coord = []
    for i in range(n_ini):
        crea_particula(ratio, Lx, Ly, mu1, sigma1, mu2, sigma2, coord)
    coord = pre_desolapamiento(Lx, Ly, coord)

    Vfake = calcula_Vfake(Lx*Ly, coord)
    print('hanga=', Vfake)

    while Vfake < lim_inf or Vfake > lim_sup:
        if Vfake < lim_inf:
            crea_particula(ratio, Lx, Ly, mu1, sigma1, mu2, sigma2, coord)
            coord = desolapamiento(Lx, Ly, coord)
        else:
            coord.pop(random.randint(0,len(coord)-1))
        coord = corrige_n(calcula_Vfake, crea_particula, desolapamiento, Lx, Ly, n_min, n_max, coord)
        Vfake = calcula_Vfake(Lx*Ly, coord)
        print('Vfake=',Vfake,'Nº de part=',len(coord))


    print('Nº de elementos de la matriz coord; Nº de partículas =',len(coord))
    print('% en volumnen ocupado por las partículas; Vfake =',calcula_Vfake(Lx*Ly,coord))
    print('Lx=',Lx)
    print('Ly=',Ly)
    print('Area de la caja = ',Lx*Ly)


    #Plot1
    plot(coord)

    #plot_2
    plot_2(coord)

    #Plot3
    plot_3(coord)

    #Plot4
    plot_4(coord,Lx,Ly)

    #Plot5
    plot_5(coord)


    mat = np.matrix(coord)
    with open('coord.txt','w+') as datafile_id:
        for line in mat:
            np.savetxt(datafile_id, line, delimiter="              ", fmt='%.10f')



#-------------------------------------------------------------------------------
#-------------------------------RESISTENCIAS------------------------------------
#-------------------------------------------------------------------------------
    d_0 = 10
    ds_max = 5*d_0                                                              #cut-off 2 ordenes de magnitud
    resist_0 = (10**(16))
    """
    h_bar = 1.054571818*(10**(-34))
    me = 9.1093897*(10**(-31))
    gap = 1*1.60217733*(10**(-19))
    alpha = np.sqrt(2*me*gap)/h_bar
    """

    electro_izquierdo, electro_derecho = electrodos(coord, Lx)
    resistencias = []
    to_write = ""
    for i,particula_i in enumerate(coord):
        for j,particula_j in enumerate(coord):
            if i<j:
                ds = dist_surface(particula_i, particula_j)
                if ds<ds_max:
                    resist = resist_0*np.exp(ds/d_0)
                    resistencias.append(resist)
                    if (i+1 in electro_izquierdo) and (j+1 in electro_izquierdo):
                        to_write += "R_iz_"+str(i+1)+"_"+str(j+1) +" "+str(0)+" "+str(j+1)+" "+str(resist)+"\n"
                    elif (i+1 in electro_derecho) and (j+1 in electro_derecho):
                        if i+1 == electro_derecho[0]:
                            pass
                        else:
                            to_write += "R_der_"+str(i+1)+"_"+str(j+1) +" "+str(i+1)+" "+str(electro_derecho[0])+" "+str(resist)+"\n"
                    else:
                        to_write += "R_"+str(i+1)+"_"+str(j+1) +" "+str(i+1)+" "+str(j+1)+" "+str(resist)+"\n"

    with open('FicheroPalSpice_Resistencias.txt','w+') as file:
        file.write(to_write)
    print('Resistencia maxima=',max(resistencias))
    print('Resistencia minima=',min(resistencias))

#-------------------------------------------------------------------------------
#--------------------------------CAPACIDADES------------------------------------
#-------------------------------------------------------------------------------

    particulas = []
    for i in coord:
        particulas.append(Particula([i[0],i[1]], i[2]))

    FinalCapacitances = []
    for index,particula in enumerate(particulas):
        vecinas = vecinasf(particula, particulas)
        q0=Carga(valor=1, dist_radial=0, coordenadas=particula.coordenadas, index_particula=index)
        particula.asignar_carga(q0)
        q0s=[q0]
        #first iter
        qprimas = calcular_qprimas(q0s, particulas, vecinas)
        calcular_coeficientes(particulas)
        print(1)
        #second iter
        qprimas = calcular_qprimas(qprimas, particulas, vecinas)
        calcular_coeficientes(particulas)
        print(2)

        for contador in range(3,12):
            qprimas = calcular_qprimas(qprimas, particulas, vecinas)
            calcular_coeficientes(particulas)
            print(contador)

        FinalCapacitances.append(obtener_coeficientes_finales(particulas, particula.radio, particula.coordenadas, vecinas))
        for particula in particulas:
            particula.reset()


    for i,fila in enumerate(FinalCapacitances):
        for j,coeficiente_ij in enumerate(FinalCapacitances):
            if i<j:
                FinalCapacitances[i][j] = (abs(FinalCapacitances[j][i])+abs(FinalCapacitances[i][j]))/2


#-----------------.--Escritura en fichero de capacidades------------------------

    electro_izquierdo, electro_derecho = electrodos(coord, Lx)
    to_write = ""
    for index_i,row in zip(range(1,len(FinalCapacitances)+1),FinalCapacitances):
        for index_j,capacitancia in zip(range(1,len(row)+1),row):
            if index_i < index_j:
                index_izq = index_i
                if index_izq in electro_derecho:
                    index_izq = -1
                elif index_izq in electro_izquierdo:
                    index_izq  = 0

                index_der = index_j
                if index_der in electro_derecho:
                    index_der = -1
                elif index_der in electro_izquierdo:
                    index_der  = 0

                if (index_izq in (0,-1)) and (index_der in (0,-1)):
                    if index_izq != index_der:
                        pass
                    else:
                        if index_izq == 0:
                            to_write+="C_izq"+str(index_i)+"_"+str(index_j) +" "+str(0)+" "+str(index_j)+" "+str(abs(row[index_j-1]))+"\n"
                        else:
                            if index_i == electro_derecho[0]:                   #CHECK THIS
                                pass
                            else:
                                to_write+="C_der"+str(index_i)+"_"+str(index_j) +" "+str(index_i)+" "+str(electro_derecho[0])+" "+str(abs(row[index_j-1]))+"\n"
                else:
                    to_write+="C_"+str(index_i)+"_"+str(index_j) +" "+str(index_i)+" "+str(index_j)+" "+str(abs(row[index_j-1]))+"\n"


    with open('FicheroPalSpiceCapacidades.txt','w+') as file:
        file.write(to_write)
    print('Electrodo derecho=',electro_derecho)
