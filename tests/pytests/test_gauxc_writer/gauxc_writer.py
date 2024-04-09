import h5py
import numpy as np
import os
import psi4

def write_results(run_name, wfn):
    try:
        #== create reference ==#
        with h5py.File(os.path.join(os.path.dirname(__file__), f'{run_name}.hdf5'), "w") as f:
            #== write molecule info ==#
            wfn_molecule = wfn.molecule()

            molecule_dtype = {
                'names': ['Atomic Number', 'X Coordinate', 'Y Coordinate', 'Z Coordinate'], 
                'formats': ['<i4', '<f8', '<f8', '<f8'], 
                'offsets': [0, 8, 16, 24], 
                'itemsize': 32
            }
           
            molecule_dset = f.create_dataset("/MOLECULE", shape=(wfn_molecule.natom(),), dtype=molecule_dtype)  
          
            for iatom in np.arange(0, wfn_molecule.natom()):
                atom_data = ( int(wfn_molecule.Z(iatom)), wfn_molecule.x(iatom), wfn_molecule.y(iatom), wfn_molecule.z(iatom) )
                molecule_dset[iatom] = atom_data 
          
            #== write basis info ==#
            basis = wfn.basisset()
            basis_dtype = {
                'names': ['NPRIM', 'L', 'PURE', 'ALPHA', 'COEFF', 'ORIGIN'], 
                'formats': ['<i4', '<i4', '<i4', ('<f8', (16,)), ('<f8', (16,)), ('<f8', (3,))], 
                'offsets': [0, 4, 8, 16, 272, 528], 
                'itemsize': 552
            }
            
            basis_dset = f.create_dataset("/BASIS", shape=(basis.nshell(),), dtype=basis_dtype)  
          
            for ishell in np.arange(0, basis.nshell()):
                basis_data = [0, 0, 0, 
                  np.zeros(16, dtype=np.float64), 
                  np.zeros(16, dtype=np.float64), 
                  np.zeros(3, dtype=np.float64)
                ]
                
                psi4_shell = basis.shell(ishell)
          
                basis_data[0] = psi4_shell.nprimitive
                basis_data[1] = psi4_shell.am
                basis_data[2] = 0 if any([ psi4_shell.is_cartesian() ]) else 1
           
                for iprim in np.arange(0, psi4_shell.nprimitive): 
                  basis_data[3][iprim] = psi4_shell.exp(iprim)
                  basis_data[4][iprim] = psi4_shell.coef(iprim)
          
                for icoord in np.arange(0,3):
                  basis_data[5][icoord] = psi4_shell.coord(icoord)
          
                basis_dset[ishell] = tuple(basis_data)
        
            initial_index = np.arange(basis.nbf())
            final_index = np.zeros(basis.nbf(), dtype=int)
            
            max_am = 7;
            cca_integral_order = [ [] for _ in np.arange(max_am) ] 
        
            # s shell, easy
            cca_integral_order[0] = [ 0 ]
        
            # p shell, just CCA ordering 
            cca_integral_order[1] = [ 0, 1, -1 ]
         
            # d shells or larger
            for l in np.arange(2, max_am):
                cca_integral_order[l] = [ 0 ]
                val = 1
                for idx in np.arange(1, 2*l + 1, 2): 
                    cca_integral_order[l].append(val);
                    cca_integral_order[l].append(-val);
                    val += 1
        
            ibf = 0 
            for ish in np.arange(0, basis.nshell()):
                sh = basis.shell(ish);
                am = sh.am;
        
                ibf_base = ibf;
                for ishbf in np.arange(0, 2*am + 1):
                   final_index[ibf] = ibf_base + cca_integral_order[am][ishbf] + am;
                   ibf += 1
        
            permutation_matrix = np.zeros((basis.nbf(), basis.nbf()), dtype=np.float64)
            permutation_matrix[final_index, initial_index] = 1.0
        
            #== write density info ==#
            D_ref = wfn.jk().D()[0].to_array()
            D_ref = np.matmul(permutation_matrix, D_ref) 
            D_ref = np.matmul(D_ref, np.transpose(permutation_matrix)) 
            D_dset = f.create_dataset("/DENSITY", data=D_ref)  
          
            #== write exchange info ==#
            K_ref = wfn.jk().K()[0].to_array()
            K_ref = np.matmul(permutation_matrix, K_ref) 
            K_ref = np.matmul(K_ref, np.transpose(permutation_matrix)) 
            K_dset = f.create_dataset("/K", data=K_ref)
        
        return True # write was successful
    except Exception as e:
        raise 

def validate_results(run_name, tol, *, write_output=False):
    #== define comparisons to make ==#
    mol_same = False # are the molecules the same?
    basis_same = False # are the basis sets the same?
    D_same = -1.0 # is the difference density RMS smaller than 1E(-tol)? 
    K_same = -1.0 # is the difference Exchange RMS smaller than 1E(-tol)?
   
    #== fill in above data ==# 
    try:
       #== open files ==#
       run = h5py.File(os.path.join(os.path.dirname(__file__), f'{run_name}.hdf5'), "r")
       ref = h5py.File(os.path.join(os.path.dirname(__file__), f'{run_name}-ref.hdf5'), "r")
  
       outfile = open(os.path.join(os.path.dirname(__file__), f'{run_name}.out'), "w") if write_output else None
 
       #== process molecules ==#
       molrun_dset = run["/MOLECULE"] 
       molrun = molrun_dset[()]
       
       molref_dset = ref["/MOLECULE"] 
       molref = molref_dset[()]
       
       dmol = molrun == molref
       mol_same = all(dmol.flatten())
       
       if write_output:
           outfile.write( "                     MOLECULE                     \n")
           outfile.write(str(mol_same) + "\n")
     
       #== process basis sets ==#
       basisrun_dset = run["/BASIS"] 
       basisrun = basisrun_dset[()]
     
       basisref_dset = ref["/BASIS"] 
       basisref = basisref_dset[()]
     
       dbasis = basisrun == basisref 
       basis_same = all(dbasis.flatten())

       if write_output:
           outfile.write("                      BASIS                       \n")
           outfile.write(str(basis_same) + "\n")
     
       #== process densities ==#
       Drun_dset = run["/DENSITY"] 
       Drun = Drun_dset[()]
     
       Dref_dset = ref["/DENSITY"] 
       Dref = Dref_dset[()]

       dD = Drun - Dref 
       dD_norm = np.linalg.norm(dD)
       dD_rms = np.sqrt(np.square(dD).mean())
       
       D_same = abs(dD_rms) <= pow(10, -tol)
 
       if write_output:
           outfile.write("                     DENSITY                      \n")

           outfile.write("Norm of dD: " + str(dD_norm) + "\n")
           dD_max_index = np.unravel_index(np.argmax(np.absolute(dD)), dD.shape)
           outfile.write("Max of dD: " + str(np.absolute(dD).max()) + " at (" + str(dD_max_index[0]) + ", " + str(dD_max_index[1]) + ")\n")
           outfile.write("RMSD of dD: " + str(dD_rms) + "\n\n")

           outfile.write("Norm of Drun - Drun^T: " + str(np.linalg.norm(Drun - np.transpose(Drun))) + "\n")
           outfile.write("Norm of Dref - Dref^T: " + str(np.linalg.norm(Dref - np.transpose(Dref))) + "\n")
           outfile.write("Norm of dD - dD^T: " + str(np.linalg.norm(dD - np.transpose(dD))) + "\n")
           outfile.write("\n\n")
 
       #== process Exchange matrices ==#
       Krun_dset = run["/K"] 
       Krun = Krun_dset[()]
       
       Kref_dset = ref["/K"] 
       Kref = Kref_dset[()]
      
       dK = Krun - Kref 
       dK_norm = np.linalg.norm(dK) 
       dK_rms = np.sqrt(np.square(dK).mean())
       
       K_same = abs(dK_rms) <= pow(10, -tol)
 
       if write_output:
           outfile.write("                    EXCHANGE                     \n")

           outfile.write("Norm of Krun: " + str(np.linalg.norm(Krun)) + "\n")
           outfile.write("Norm of Kref: " + str(np.linalg.norm(Kref)) + "\n")
           outfile.write("Norm of dK: " + str(dK_norm) + "\n\n")

           dK_max_index = np.unravel_index(np.argmax(np.absolute(dK)), dK.shape)
           outfile.write("Max of dK: " + str(np.absolute(dK).max()) + " at (" + str(dK_max_index[0]) + ", " + str(dK_max_index[1])  + ")\n")
           outfile.write("RMS of dK: " + str(dK_rms) + "\n\n")

           outfile.write("Norm of Krun - Krun^T: " + str(np.linalg.norm(Krun - np.transpose(Krun))) + "\n")
           outfile.write("Norm of Kref - Kref^T: " + str(np.linalg.norm(Kref - np.transpose(Kref))) + "\n")
           outfile.write("Norm of dK - dK^T: " + str(np.linalg.norm(dK - np.transpose(dK))) + "\n")
           outfile.write("\n")
   
           outfile.close()
 
    except:
        raise
        
    return mol_same, basis_same, D_same, K_same 
