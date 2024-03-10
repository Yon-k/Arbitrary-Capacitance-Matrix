import numpy as np
from pony.orm import *
import os
db = Database()


class Carga(db.Entity):
    valor = Required(float) #valor de la carga
    dist_radial = Required(float) #distancia des de el centro de la particula
    coordenadas = Required(FloatArray) #x,y,z
    index_particula = Required(int) #indice de la particula a la que pertenece
    iteracion = Required(int) #Numero que marca la iteracion a la que corresponde

#----------------------------------------------------------------------------------
db.bind(provider='sqlite', filename='database.sqlite', create_db=True)
db.generate_mapping(create_tables=True)

class Particula(object):
    def __init__(self, coordenadas, radio):
        self.coordenadas = np.array(coordenadas) #x,y,z
        self.radio = radio
        self.cargas = 0
        self.coeficientes_influencia = []

    def asignar_carga(self, valor):
        """Este metodo me guarda una carga en el valor de cargas de arriba.
        """
        self.cargas += valor

    def coeficiente_influencia(self):
        return self.cargas

    def coeficiente_final(self):
        return self.coeficientes_influencia[-1]

    @db_session
    def Qprimas(self, carga0, particulas):
        """Calculo todas las Qprimas dada una carga0 concreta. Me lo guardo todo
        en la database.
        """
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
                carga_prima = Carga(valor=Q_prima, dist_radial=dist_radial_carga, coordenadas=coordenadasQ_prima, index_particula=j, iteracion=carga0.iteracion+1)
                particle.asignar_carga(carga_prima.valor)

    def reset(self):
        self.cargas = 0
        self.coeficientes_influencia = []

################################################################################
################################################################################

#Inicializacion del sistema
R1 = 10**(-9)
R2 = 10**(-9)
R3 = 10**(-9)
eps_r = 1   #sera 23 en un futuro
eps= 8.8541878176*(10**(-12))*eps_r
tolerancia = 10**(-7)
particulas = [Particula([-3*10**(-9),0,0], R1), Particula([0,0,0], R2), Particula([3*10**(-9),0,0], R3)]
#particulas = [Particula([-6*10**(-9),0,0], R1), Particula([-3*10**(-9),0,0], R1), Particula([0,0,0], R2), Particula([3*10**(-9),0,0], R3),
#Particula([6*10**(-9),0,0], R3)]
#particulas = [Particula([0,0,0], R1), Particula([3*10**(-9),0,0], R2)]

#-----------------------------------------------------------------------------
@db_session
def calcular_qprimas(particulas, iteracion):                                    #iteracion de la anterior
    medida_paquete = int(10e4)
    numQ = select(q for q in Carga if q.iteracion == iteracion).count()         #numero de cargas de la iteracion iteracion
    paquete = (numQ // medida_paquete) + 1
    for i in range(paquete):
        q0s = select(q for q in Carga if q.iteracion == iteracion).sort_by(Carga.valor).page(i+1, pagesize=medida_paquete)
        for carga in q0s:
            for index,particula in enumerate(particulas):
                if carga.index_particula==index:
                    particula.Qprimas(carga, particulas)


def calcular_coeficientes(particulas):
    for particula in particulas:
        coef = particula.coeficiente_influencia()
        particula.coeficientes_influencia.append(coef)


def converge(particulas):
    errores_relativos = []
    for particula in particulas:
        errores_relativos.append(abs(particula.coeficientes_influencia[-1]-particula.coeficientes_influencia[-2]))
    Convergencia_metodo.append((particulas[0].coeficientes_influencia[-1], particulas[1].coeficientes_influencia[-1]))
    return max(errores_relativos) < tolerancia


def obtener_coeficientes_finales(particulas, radius):
    coeficientes = []
    for particula in particulas:
        coeficientes.append(4*np.pi*eps*radius*particula.coeficiente_final())
    return coeficientes


@db_session
def carga_inicial():
    """La primera q0 me la guardo en la database.
    """
    q0=Carga(valor=1, dist_radial=0, coordenadas=particula.coordenadas, index_particula=index, iteracion=0)
    particula.asignar_carga(q0.valor)


@db_session
def clear_database():
    delete(q for q in Carga)


#----------------------------------MAIN---------------------------------------

FinalCapacitances = []
Convergencia_metodo = []

for index,particula in enumerate(particulas):
    carga_inicial()
    #first iter
    calcular_qprimas(particulas, 0)
    calcular_coeficientes(particulas)
    print(1)
    #second iter
    calcular_qprimas(particulas, 1)
    calcular_coeficientes(particulas)
    print(2)
    #iters
    contador = 2
    while(not converge(particulas)):
        calcular_qprimas(particulas, contador)
        calcular_coeficientes(particulas)
        contador += 1
        print(contador)
        if contador==10:
            break

    FinalCapacitances.append(obtener_coeficientes_finales(particulas, particula.radio))
    for particula in particulas:
        particula.reset()
    clear_database()

to_write = ""
for index_i,row in enumerate(FinalCapacitances):
    for index_j,elemento in enumerate(row):
        if index_i==index_j:
            new_coef = elemento-(sum(np.abs(row))-elemento)
            row.pop(index_j)
            row.insert(index_j, new_coef)
        to_write+="C_"+str(index_i)+"_"+str(index_j) +" "+ str(abs(row[index_j]))+"\n"

with open('FicheroPalSpice_Capacidades.txt','w+') as file:
    file.write(to_write)
with open('Coeficiente_convergencia.csv','w+') as file:
    for line in Convergencia_metodo:
        file.write(str(line)[1:-1]+"\n")                                        #Esto me quita los parentesis

capacitances = np.reshape(np.abs(FinalCapacitances),(len(particulas),len(particulas)))
print(capacitances)




db.disconnect()
os.remove("database.sqlite")
