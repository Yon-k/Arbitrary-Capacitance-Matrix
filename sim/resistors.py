# Constants
d_0 = 10
ds_max = 5 * d_0
resistance_0 = 10 ** 16

def calculate_resistances(coordinates, Lx):
    left_electrode, right_electrode = electrodes(coordinates, Lx)
    resistances = []
    resistances_info = ""

    for i, particle_i in enumerate(coordinates):
        for j, particle_j in enumerate(coordinates):
            if i < j:
                ds = surface_distance(particle_i, particle_j)
                if ds < ds_max:
                    resistance = resistance_0 * np.exp(ds / d_0)
                    resistances.append(resistance)
                    resistances_info += format_resistance_info(i, j, resistance, left_electrode, right_electrode)

    write_to_file('SpiceResistancesFile.txt', resistances_info)
    print('Maximum resistance=', max(resistances))
    print('Minimum resistance=', min(resistances))

def format_resistance_info(i, j, resistance, left_electrode, right_electrode):
    if (i + 1 in left_electrode) and (j + 1 in left_electrode):
        return f"R_left_{i + 1}_{j + 1} 0 {j + 1} {resistance}\n"
    elif (i + 1 in right_electrode) and (j + 1 in right_electrode):
        if i + 1 != right_electrode[0]:
            return f"R_right_{i + 1}_{j + 1} {i + 1} {right_electrode[0]} {resistance}\n"
    else:
        return f"R_{i + 1}_{j + 1} {i + 1} {j + 1} {resistance}\n"

def write_to_file(filename, content):
    with open(filename, 'w+') as file:
        file.write(content)