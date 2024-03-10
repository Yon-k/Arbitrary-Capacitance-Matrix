
particles = []
for i in coord:
    particles.append(Particle([i[0],i[1]], i[2]))

final_capacitances = []
for index, particle in enumerate(particles):
    neighbors = find_neighbors(particle, particles)
    q0 = Charge(value=1, radial_dist=0, coordinates=particle.coordinates, particle_index=index)
    particle.assign_charge(q0)
    q0s = [q0]
    #first iter
    q_primes = calculate_q_primes(q0s, particles, neighbors)
    calculate_coefficients(particles)
    print(1)
    #second iter
    q_primes = calculate_q_primes(q_primes, particles, neighbors)
    calculate_coefficients(particles)
    print(2)

    for counter in range(3,12):
        q_primes = calculate_q_primes(q_primes, particles, neighbors)
        calculate_coefficients(particles)
        print(counter)

    final_capacitances.append(get_final_coefficients(particles, particle.radius, particle.coordinates, neighbors))
    for particle in particles:
        particle.reset()

for i, row in enumerate(final_capacitances):
    for j, coefficient_ij in enumerate(final_capacitances):
        if i<j:
            final_capacitances[i][j] = (abs(final_capacitances[j][i])+abs(final_capacitances[i][j]))/2

#-----------------.--Writing capacitances to file------------------------

left_electrode, right_electrode = electrodes(coord, Lx)
to_write = ""
for index_i, row in zip(range(1, len(final_capacitances)+1), final_capacitances):
    for index_j, capacitance in zip(range(1, len(row)+1), row):
        if index_i < index_j:
            left_index = index_i
            if left_index in right_electrode:
                left_index = -1
            elif left_index in left_electrode:
                left_index  = 0

            right_index = index_j
            if right_index in right_electrode:
                right_index = -1
            elif right_index in left_electrode:
                right_index  = 0

            if (left_index in (0,-1)) and (right_index in (0,-1)):
                if left_index != right_index:
                    pass
                else:
                    if left_index == 0:
                        to_write += "C_left"+str(index_i)+"_"+str(index_j) +" "+str(0)+" "+str(index_j)+" "+str(abs(row[index_j-1]))+"\n"
                    else:
                        if index_i == right_electrode[0]:                   #CHECK THIS
                            pass
                        else:
                            to_write += "C_right"+str(index_i)+"_"+str(index_j) +" "+str(index_i)+" "+str(right_electrode[0])+" "+str(abs(row[index_j-1]))+"\n"
            else:
                to_write += "C_"+str(index_i)+"_"+str(index_j) +" "+str(index_i)+" "+str(index_j)+" "+str(abs(row[index_j-1]))+"\n"

with open('SpiceCapacitancesFile.txt','w+') as file:
    file.write(to_write)
print('Right electrode=', right_electrode)