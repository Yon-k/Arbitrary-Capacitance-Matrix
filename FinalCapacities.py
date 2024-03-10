import numpy as np

class Carga(object):
    def __init__(self, valor, dist_radial, coordenadas, index_particula):
        self.valor = valor #valor de la carga
        self.dist_radial = dist_radial #distancia des de el centro de la particula
        self.coordenadas = coordenadas #x,y,z
        self.index_particula = index_particula #indice de la particula a la que pertenece

#----------------------------------------------------------------------------------

class Particula(object):
    def __init__(self, coordenadas, radio):
        self.coordenadas = np.array(coordenadas) #x,y,z
        self.radio = radio
        self.cargas = 0
        self.coeficientes_influencia = []


    def asignar_carga(self, carga):
        """Este metodo me guarda una carga en el valor de cargas de arriba.
        """
        self.cargas += carga.valor


    def coeficiente_influencia(self):
        return self.cargas


    def coeficiente_final(self):
        return self.coeficientes_influencia[-1]


    def Qprimas(self, carga0, particulas):
        """Calculo todas las Qprimas dada una carga0 concreta.
        """
        qprimas = []
        for j,particle in enumerate(particulas):
            if self==particle:
                continue
            else:
                #1a parte: parametros de la carga.
                d = np.linalg.norm(self.coordenadas - particle.coordenadas)
                Q_prima = -(carga0.valor*particle.radio)/(d-carga0.dist_radial)
                dist_radial_carga = (particle.radio**2)/(d-carga0.dist_radial)
                coordenadasQ_prima = particle.coordenadas-((particle.radio/d)**2)*(particle.coordenadas - self.coordenadas)
                #2a parte: creacion objeto carga.
                carga_prima = Carga(Q_prima, dist_radial_carga, coordenadasQ_prima,  j)
                particle.asignar_carga(carga_prima)
                qprimas.append(carga_prima)
        return qprimas


    def reset(self):
        self.cargas = 0
        self.coeficientes_influencia = []

#-------------------------------------------------------------------------------
#---------------------------------FUNCTIONS-------------------------------------
#-------------------------------------------------------------------------------

def calcular_qprimas(q0s, particulas):
    qprimas = []
    for carga in q0s:
        for index,particula in enumerate(particulas):
            if carga.index_particula==index:
                qprimas += particula.Qprimas(carga, particulas)
    return qprimas


def calcular_coeficientes(particulas):
    for particula in particulas:
        coef = particula.coeficiente_influencia()
        particula.coeficientes_influencia.append(coef)


def converge(particulas):
    tolerancia = 10**(-7)
    errores_relativos = []
    for particula in particulas:
        errores_relativos.append(abs(particula.coeficientes_influencia[-1]-particula.coeficientes_influencia[-2]))
    return max(errores_relativos) < tolerancia


def obtener_coeficientes_finales(particulas, radius):
    eps_r = 1
    eps= 8.8541878176*(10**(-12))*eps_r
    coeficientes = []
    for particula in particulas:
        coeficientes.append(4*np.pi*eps*radius*particula.coeficiente_final())
        print(coeficientes)
    return coeficientes

#----------------------------------MAIN---------------------------------------

if __name__== "__main__":
    #Inicializacion del sistema
    R1 = 10**(-9)
    R2 = 10**(-9)
    R3 = 10**(-9)
    #particulas = [Particula([0,0,0], 1), Particula([5,0,0],0.5)]
    #particulas = [Particula([-3*10**(-9),0,0], R1), Particula([0,0,0], R2), Particula([3*10**(-9),0,0], R3)]
    #particulas = [Particula([-6*10**(-9),0,0], R1), Particula([-3*10**(-9),0,0], R1), Particula([0,0,0], R2), Particula([3*10**(-9),0,0], R3),]
    #Particula([6*10**(-9),0,0], R3)]
    particulas = [Particula([0,0,0], R1), Particula([3*10**(-9),0,0], R2)]


    FinalCapacitances = []
    for index,particula in enumerate(particulas):
        q0=Carga(valor=1, dist_radial=0, coordenadas=particula.coordenadas, index_particula=index)
        particula.asignar_carga(q0)
        q0s=[q0]
        #first iter
        qprimas = calcular_qprimas(q0s, particulas)
        calcular_coeficientes(particulas)
        print(1)
        #second iter
        qprimas = calcular_qprimas(qprimas, particulas)
        calcular_coeficientes(particulas)
        print(2)

        contador = 2
        while(not converge(particulas)):
            qprimas = calcular_qprimas(qprimas, particulas)
            calcular_coeficientes(particulas)
            contador += 1
            print(contador)
            if contador==100:
                break

        FinalCapacitances.append(obtener_coeficientes_finales(particulas, particula.radio))
        for particula in particulas:
            particula.reset()


    for index_i,row in enumerate(FinalCapacitances):
        for index_j,elemento in enumerate(row):
            if index_i==index_j:
                new_coef = elemento-(sum(np.abs(row))-elemento)
                row.pop(index_j)
                row.insert(index_j, new_coef)

    capacitances = np.reshape(np.abs(FinalCapacitances),(len(particulas),len(particulas)))
    print(capacitances)
