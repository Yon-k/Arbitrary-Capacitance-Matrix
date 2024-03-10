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

def vecinasf(particula_i, particulas):
    distancias = np.array([np.linalg.norm(particula_i.coordenadas - particulas[i].coordenadas) for i in range(len(particulas))])
    primeros_indices = np.argsort(distancias)[:4] #Las N primeras , mas ella misma
    return np.sort(primeros_indices)


def calcular_qprimas(q0s, particulas,vecinasf):
    qprimas = []
    vecinas = [particulas[i] for i in vecinasf]
    for carga in q0s:
        for index,particula in enumerate(particulas):
            if carga.index_particula==index:
                qprimas += particula.Qprimas(carga, vecinas)
    return qprimas


def calcular_coeficientes(particulas):
    for i, particula in enumerate(particulas):
        coef = particula.coeficiente_influencia()
        particula.coeficientes_influencia.append(coef)


def obtener_coeficientes_finales(particulas, radius, coords, vecinas):
    eps_r = 23
    eps= 8.8541878176*(10**(-12))*eps_r
    coeficientes = []
    for i, particula in enumerate(particulas):
        if i in vecinas:
            coeficientes.append(4*np.pi*eps*radius*(10**(-9))*particula.coeficiente_final())
        else:
            d = np.linalg.norm(particula.coordenadas - coords)
            if d != 0:
                coeficientes.append(4*np.pi*eps*particula.radio*radius*(10**(-9))/d)
            else:
                coeficientes.append(0)

    return coeficientes

#----------------------------------MAIN---------------------------------------

if __name__== "__main__":

    coord = [[3,4,5],[6,38,4],[10,4,0.3],[20,4,1],[30,6,1]]
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

        for contador in range(3,10):
            qprimas = calcular_qprimas(qprimas, particulas, vecinas)
            calcular_coeficientes(particulas)
            print(contador)

        FinalCapacitances.append(obtener_coeficientes_finales(particulas, particula.radio, particula.coordenadas, vecinas))
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
