#include <psi4-dec.h>
#include <libciomr/libciomr.h>
#include <libmints/molecule.h>
#include <libmints/pointgrp.h>
#include <libmints/petitelist.h>
#include <libmints/cdsalclist.h>

#include <algorithm>

using namespace boost;

namespace psi {

bool cdsalc_sort_predicate(const CdSalc& d1, const CdSalc& d2)
{
    return d1.irrep() < d2.irrep();
}

void CdSalc::print() const
{
    fprintf(outfile, "\tirrep = %d, ncomponent = %ld\n", irrep_, ncomponent());
    for (int i=0; i<ncomponent(); ++i) {
        fprintf(outfile, "\t\t%d: atom %d, direction %d, coef %lf\n",
                i,
                components_[i].atom,
                components_[i].xyz,
                components_[i].coef);
    }
}

CdSalcList::CdSalcList(const boost::shared_ptr<Molecule>& mol,
                       char needed_irreps,
                       bool project_out_translations,
                       bool project_out_rotations)
    : molecule_(mol), needed_irreps_(needed_irreps),
      project_out_translations_(project_out_translations),
      project_out_rotations_(project_out_rotations)
{
    // Ensure point group has been set.
    if (!molecule_->point_group()) {
        throw PSIEXCEPTION("CdSalcList::CdSalcList: Molecule point group has not been set.");
    }

    // Ideally this could be 3n-5 or 3n-6.
    ncd_ = 3 * molecule_->natom();
    double *salc = new double[ncd_];

    // Obtain handy reference to point group.
    PointGroup& pg = *molecule_->point_group().get();
    CharacterTable char_table = pg.char_table();
    nirrep_ = char_table.nirrep();

    salc_symblock_ = new double**[nirrep_];
    for (int i=0; i<nirrep_; ++i)
        salc_symblock_[i] = new double*[ncd_];

    cdsalc2cd_ = block_matrix(ncd_, ncd_);

    // Clear the counting array
    memset(cdsalcpi_, 0, sizeof(int) * 8);

    // Obtain atom mapping of atom * symm op to atom
    int **atom_map = compute_atom_map(molecule_);

    print_int_mat(atom_map, molecule_->natom(), nirrep_, outfile);

    for (int uatom=0; uatom < molecule_->nunique(); ++uatom) {
        int atom = molecule_->unique(uatom);

        // Project each displacement
        for (int xyz=0; xyz<3; ++xyz) {

            // on each irrep
            for (int irrep=0; irrep<nirrep_; ++irrep) {
                IrreducibleRepresentation gamma = char_table.gamma(irrep);
                memset(salc, 0, sizeof(double)*ncd_);

                // This is the order of the atom stabilizer
                int stab_order = 0;

                // Apply the projection
                for (int G=0; G<nirrep_; ++G) {
                    SymmetryOperation so = char_table.symm_operation(G);
                    int Gatom = atom_map[atom][G];
                    if (Gatom == atom)
                        ++stab_order;

                    int Gcd = 3*Gatom + xyz;
                    double coeff = so(xyz, xyz) * gamma.character(G);

                    salc[Gcd] += coeff;
                }

//                fprintf(outfile, "salc atom %d, xyz %d, h %d = [", atom, xyz, irrep);
                int nonzero=0;
                for (int cd=0; cd<ncd_; ++cd) {
                    // Normalize the salc
                    salc[cd] /= sqrt((double)nirrep_*stab_order);

                    // Count number of nonzeros
                    if (fabs(salc[cd]) > 1e-10)
                        ++nonzero;

//                    fprintf(outfile, " %lf,", salc[cd]);
                }
//                fprintf(outfile, "]\n");

                // We're only interested in doing the following if there are nonzeros
                // AND the irrep that we're on is what the user wants.
                if (nonzero && (1 << irrep) & needed_irreps) {
                    CdSalc new_salc(irrep);

                    for (int cd=0; cd<ncd_; ++cd) {
                        if (fabs(salc[cd]) > 1e-10)
                            new_salc.add(salc[cd], cd/3, xyz);
                    }

                    // Save this salc.
                    salcs_.push_back(new_salc);
                }

                for (int i=0; i<ncd_; ++i) {
                    if (fabs(salc[i]) > 1e-10) {
                        salc_symblock_[irrep][cdsalcpi_[irrep]] = new double[ncd_];
                        memcpy(salc_symblock_[irrep][cdsalcpi_[irrep]],salc,sizeof(double)*ncd_);
                        ++cdsalcpi_[irrep];
                        break;
                    }
                }
            }
        }
    }

    // Sort the salcs based on irrep
    std::sort(salcs_.begin(), salcs_.end(), cdsalc_sort_predicate);

    /* copy salc_symblk to cdsalc2cd */
    {
        int c = 0;
        for(int irrep=0; irrep<nirrep_; irrep++) {
            int num_per_irrep = cdsalcpi_[irrep];
            for(int i=0; i<num_per_irrep; i++,c++) {
                for(int j=0; j<ncd_; j++) {
                    cdsalc2cd_[j][c] = salc_symblock_[irrep][i][j];
                }
            }
        }
    }

    fprintf(outfile,"    -Cartesian displacement SALCs per irrep:\n");
    fprintf(outfile,"    Irrep  #SALCs\n");
    fprintf(outfile,"    -----  ------\n");
    for (int irrep=0; irrep<nirrep_; irrep++) {
        fprintf(outfile,"    %3d    %4d\n", irrep, cdsalcpi_[irrep]);
    }
    fprintf(outfile,"\n");

    fprintf(outfile,"    -Cartesian displacement SALCs:\n");
    print_mat(cdsalc2cd_, ncd_, ncd_, outfile);
    fprintf(outfile,"\n");

    // Free memory.
    delete[] salc;
    delete_atom_map(atom_map, molecule_);
}

CdSalcList::~CdSalcList()
{
    /* copy salc_symblk to cdsalc2cd */
    {
        int c = 0;
        for(int irrep=0; irrep<nirrep_; irrep++) {
            int num_per_irrep = cdsalcpi_[irrep];
            for(int i=0; i<num_per_irrep; i++,c++) {
                delete[] salc_symblock_[irrep][i];
            }
            delete[] salc_symblock_[irrep];
        }
        delete[] salc_symblock_;
    }
}

void CdSalcList::print() const
{
    fprintf(outfile, "Number of SALCs %ld, irreps = %d\n"
            "Project out translations %s (not programmed)\n"
            "Project out rotations %s (not programmed)\n",
            salcs_.size(), needed_irreps_,
            project_out_translations_ ? "true" : "false",
            project_out_rotations_ ? "true" : "false");

    for (int i=0; i<salcs_.size(); ++i)
        salcs_[i].print();
}

}
