from polymer_maker import polymer_maker
import os, sys

def oss(fun):
    os.system(fun)

def gmx(*fun):
    gmxrun = '/usr/local/gromacs_gmx/bin/gmx'
    for arg in fun:
        gmxrun += ' '+arg
    oss(gmxrun)

if __name__ == "__main__":
    oss('export PATH=$PATH:/usr/local/openmpi/bin')
    oss('export OMP_NUM_THREADS=4')

    cation_type = sys.argv[1] # Cation type to be used
    anion_type = sys.argv[2] # Anion type to be used
    ion_concentration = sys.argv[3] # Concentration of ions in system
    monomer = sys.argv[4]
    n_polymers = sys.argv[5] # Number of polymers in system
    dop = sys.argv[6]
    trial_number = sys.argv[7]

    # Create polymer_formula file to describe polymer
    with open('polymer_formula','w') as polyf:
        polyf.write("%s\nA %s" % ('A'*dop, monomer))

    # Create linear polymers with 120 degree rotation angle
    polymer_maker(n_polymers, 120)

    # Generate multi-polymer simulation system
    gmx('insert-molecules','-ci newmol.gro','-nmol '+str(n_polymers),'-box 20 20 20','-o box.gro','-rot xyz','-try 200')

    # Minimize polymers in vacuum
    gmx('grompp','-f source/minim.mdp','-c box.gro','-p topol.top','-o em.tpr')
    gmx('mdrun','-deffnm em')

    # Solvate the system and minimize energy
    # SPC/E (spc216) water model chosen for its accurate representation of the dipole of water
    gmx('solvate','-cp em.gro','-cs spc216.gro','-p topol.top','-o wetbox.gro')
    gmx('grompp','-f source/minim.mdp','-c wetbox.gro','-p topol.top','-o em_wet.tpr')
    gmx('mdrun','-deffnm em_wet')

    # Run isothermal-isobaric (NPT) simulation to equilibrate the system
    gmx('grompp','-f source/npt.mdp','-c em_wet.gro','-p topol.top','-o npt.tpr')
    gmx('mdrun','-deffnm npt')

    # Check that the NPT simulation completed correctly
    nptflag = 0
    while not os.path.isfile('npt.gro'):
        gmx('grompp','-f source/npt.mdp','-c em_wet.gro','-p topol.top','-o npt.tpr')
        gmx('mdrun','-deffnm npt')
        nptflag += 1
        if nptflag == 10:
            print('ERROR: NPT simulation failed')
            break

    # Add ions to the system
    gmx('grompp','-f source/ions.mdp','-c npt.gro','-p topol.top','-o ions.tpr')
    gmx('echo "\SOL\" | /usr/local/gromacs_gmx/bin/gmx genion -s ions.tpr -pname '+cation_type+' -nname '+anion_type+' -conc '+conc+' -o ionwetbox.gro -p topol.top')

    # Run isothermal-isochoric (NVT) simulation for data production
    gmx('grompp','-f source/nvt.mdp','-c ionwetbox.gro','-p topol.top','-o nvt.tpr')
    gmx('mdrun','-deffnm nvt')

    # Make indices for ions
    index = 'printf \"r '+cation_type+' TO_USE \n r '+anion_type+' TO_USE \n q \n\"'
    oss(index+' | /usr/local/gromacs_gmx/bin/gmx make_ndx -f nvt.gro')

    # Calculate cation and anion mean-squared displacement
    oss('echo \"'+cation_type+'_TO_USE\" | gmx msd -f nvt.trr -s nvt.gro -n index.ndx -o msd_cation.xvg')
    oss('echo \"'+anion_type+'_TO_USE\" | gmx msd -f nvt.trr -s nvt.gro -n index.ndx -o msd_anion.xvg')

    ### Calculate ionic conductivity ###
    # Read simulation structure file to calculate box volume and number of ions
    cation_n = 0
    anion_n = 0
    with open('nvt.gro','r') as boxf:
        doc = boxf.readlines()
        for line in doc:
            if line.split()[1] == cation_type:
                cation_n += 1
            elif line.split()[1] == anion_type:
                anion_n += 1
        box_x, box_y, box_z = [float(i) for i in doc[-1].split()]

     # Grab cation and anion diffusion coefficients
     with open('msd_cation.xvg','r') as f:
         doc = f.readlines()
         for line in doc:
             if '# D[  '+cation_type+'_TO_USE]' in line:
                 cation_d = float(line.split()[4]+line.split()[-2][2:])
    with open('msd_anion.xvg','r') as f:
        doc = f.readlines()
        for line in doc:
            if '# D[  '+anion_type+'_TO_USE]' in line:
                anion_d = float(line.split()[4]+line.split()[-2][2:])

     volume = box_x*box_y*box_z # System volume
     k_b = 1.380649E-23 # Boltzmann constant, J/K [=] (kg*m^2)/(K*s^2)
     T = 298.0 # Temperature in K
     e = 1.602176634E-19 # Elementary Charge, C [=] A*s

     # Calculated from Nernst-Einstein equation
     conductivity = ()(e**2)/(volume*k_b*T))*(cation_n*cation_d+anion_n*anion_d)
     with open('conductivity','w') as condf:
         condf.write('Conductivity: %f' % (conductivity))
