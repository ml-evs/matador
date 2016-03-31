#!/usr/bin/python

def castep2py(seed):

    castep = dict()
    
    # read .castep file
    with open(seed+'.castep', 'r') as f:
        flines = f.readlines()

    if 'CollCode' in seed:
        castep['icsd_ref'] = seed.split('CollCode')[-1] 

    # wrangle castep file for basic parameters
    for line_no, line in enumerate(flines):
        if 'type of calculation' in line:
            castep['task'] = line.split(':')[-1].strip()
        elif 'functional' in line:
            castep['xc_functional'] = line.split(':')[-1].strip()
        elif 'plane wave basis set' in line:
            castep['cut_off_energy'] = float(line.split(':')[-1].split()[0])
        elif 'finite basis set correction  ' in line:
            castep['finite_basis_corr'] = line.split(':')[-1].strip()
        elif 'MP grid size for SCF' in line:
            castep['kpoints_mp_grid'] = map(int, list(line.split('is')[-1].split()))
        elif 'Number of kpoints used' in line:
            castep['kpoints_calculated'] = int(line.split('=')[-1])
        elif 'Space group of crystal' in line:
            castep['space_group'] = line.split(':')[-1].split(',')[0].strip()
        elif 'External pressure/stress' in line:
            castep['external_pressure'] = list()
            castep['external_pressure'].append(map(float, flines[line_no+1].split()))
            castep['external_pressure'].append(map(float, flines[line_no+2].split()))
            castep['external_pressure'].append(map(float, flines[line_no+3].split()))
        elif 'Cell Contents' in line:
            castep['atomic_types'] = list()
            i = 1
            atoms = False
            while True:
                if atoms:
                    if 'xxxxxxxxx' in flines[line_no+i]:
                        atoms = False
                        break
                    else:
                        castep['atomic_types'].append(flines[line_no+i].split()[1])
                if 'x------' in flines[line_no+i]:
                    atoms= True
                i += 1
        elif 'Mass of species in AMU' in line:
            i = 1
            atomic_masses = list()
            while True:
                if len(flines[line_no+i].strip())==0:
                    break
                else:
                    atomic_masses.append([float(flines[line_no+i].split()[1]), flines[line_no+i].split()[0]])
                i += 1
            atomic_masses.sort(reverse=True)
        elif 'Files used for pseudopotentials' in line:
            castep['species_pot'] = dict()
            i = 1
            while True:
                if len(flines[line_no+i].strip())==0:
                    break
                else:
                    castep['species_pot'][flines[line_no+i].split()[0].strip()] = flines[line_no+i].split()[1].strip()
                    i += 1

    # generate mass-ordered chemical composition string
    castep['composition'] = ''
    count = list()
    for ind, species in enumerate(atomic_masses):
        count_tmp = 0
        for atom in castep['atomic_types']:
            if species[1]==atom:
                count_tmp += 1
        count.append(count_tmp)
    reducible = True
    for num in count:
        if num % min(count) != 0:
            reducible = False
    min_count = min(count)
    if reducible:
        count = [num/min_count for num in count]
    for ind, species in enumerate(atomic_masses):
        castep['composition'] += species[1]
        if count[ind] != 1: 
            castep['composition'] += str(count[ind])
    # task specific options
    if castep['task'] == 'geometry optimization':
        final = False
        for line_no, line in enumerate(flines):
            if 'WARNING - Geometry optimization failed to converge' in line:
                castep['optimised'] = False
            if 'Final Configuration' in line:
                final = True
            if final:
                if 'Real Lattice' in line:
                    castep['lattice_cart'] = list()
                    i = 1
                    while True:
                        if len(flines[line_no+i].strip()) == 0:
                            break
                        else:
                            castep['lattice_cart'].append((map(float, (flines[line_no+i].split()[0:3]))))
                        i += 1
                elif 'Lattice parameters' in line :
                    castep['lattice_abc'] = list()
                    i = 1
                    if 'Current cell volume' in flines[line_no+i]:
                        break
                    else:
                        castep['lattice_abc'].append([flines[line_no+i].split('=')[1].strip().split(' ')[0], 
                                                     flines[line_no+i+1].split('=')[1].strip().split(' ')[0],
                                                     flines[line_no+i+2].split('=')[2].strip().split(' ')[0]])
                        castep['lattice_abc'].append([flines[line_no+i].split('=')[-1],
                                                     flines[line_no+i+1].split('=')[-1],
                                                     flines[line_no+i+2].split('=')[-1]])
                elif 'Cell Contents' in line :
                    castep['positions_frac'] = list()
                    i = 1
                    atoms = False
                    while True:
                        if atoms:
                            if 'xxxxxxxxx' in flines[line_no+i]:
                                atoms = False
                                break
                            else:
                                castep['positions_frac'].append((map(float, (flines[line_no+i].split()[3:6]))))
                        if 'x------' in flines[line_no+i]:
                            atoms = True
                        i += 1
                elif 'Forces' in line:
                    i = 1
                    max_force = 0
                    forces = False
                    while True:
                        if forces:
                            if '*' in flines[line_no+i].split()[1]:
                                forces = False
                                break
                            else:
                                force_on_atom = 0
                                for j in range(1):
                                    force_on_atom += float(flines[line_no+i].split()[3+j])**2
                                if force_on_atom > max_force:
                                    max_force = force_on_atom
                        elif 'x' in flines[line_no+i]:
                            i += 1                      # skip next blank line
                            forces = True
                        i += 1
                    castep['max_force_on_atom'] = pow(max_force, 0.5)
                elif 'Stress Tensor' in line:
                    i = 1
                    while True and i < 20:
                        if 'Pressure' in flines[line_no+i]:
                            castep['pressure_on_cell'] = float(flines[line_no+i].split()[-2])
                        i += 1
                elif 'Final Enthalpy' in line:
                    castep['final_enthalpy'] = float(line.split('=')[-1].split()[0])
                elif 'Final bulk modulus' in line:
                    castep['final_bulk_modulus'] = float(line.split('=')[-1].split()[0])

    # computing metadata, i.e. parallelism, time, memory, version
    for line in flines:
        if 'Release CASTEP version' in line:
            castep['castep_version'] = line.split()[-2]
        elif 'Total time' in line:
            castep['total_time_hrs'] = float(line.split()[-2])/3600
        elif 'Peak Memory Use' in line:
            castep['peak_mem_MB'] = int(float(line.split()[-2])/1000)
    print castep['castep_version']
    
    import pymongo as pm
    client = pm.MongoClient()
    db = client.crystals
    lcmo = db.lcmo
    struct_id = lcmo.insert_one(castep).inserted_id
    

if __name__ == '__main__':
    castep2py('La2CoMnO6CollCode195072')
