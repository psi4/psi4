/*! \file
    \ingroup OEPROP
    \brief Enter brief description of file here 
*/
#define EXTERN
#include <libiwl/iwl.h>
#include <libqt/qt.h>
#include "includes.h"
#include "globals.h"
#include "prototypes.h"
#include <psifiles.h>
#include <ccfiles.h>
#define TOL 1E-14

namespace psi { namespace oeprop {

void read_density()
{ 
  int i,j,k,l,dim_i,count,ibf;
  int *locs;
  double **psq_so, *tmp_arr, **tmp_mat, **psq_ao;
  int irrep, mo_offset, so_offset, i_ci, j_ci, max_opi, errcod;
  int fzc, populated_orbs;
  int *docc, *socc, *frozen_docc, *frozen_uocc, *reorder, **ras_opi; 
  int *rstr_docc, *rstr_uocc;
  int *reorder_a, *reorder_b;
  double **scfvec, **opdm_blk, **onepdm;
  char opdm_key[80];
  int nso=0;
  int maxmopi=0;  /* Maximum number of MOs per irrep */
  int *mopi;      /* MOs per irrep */
  int *doccpi;    /* DOCC per irrep */
  int *soccpi;    /* SOCC per irrep */
  int *fzdoccpi;  /* Frozen DOCC per irrep */
  int *fzvirtpi;  /* Frozen VIRT per irrep */
  double **aopdm;  /* Alpha OPDM */
  double **bopdm;  /* Beta OPDM */
  double **aopdm_blk;  /* Symmetry Blocked Alpha OPDM */
  double **bopdm_blk;  /* Symmetry Blocked Beta OPDM */
  double **ascfvec; 
  double **bscfvec;
  double **aPsq_so;
  double **bPsq_so;
  double **P_so_tot;
  double **P_mo_tot;
  double **tmat;
  double *Stri;
  double **Smat;
  double *P_eigvals;
  double **P_eigvecs;
  double **P_mo_block;
  double **P_so_block;
  double **scfvec_irrep;
  double *eval_tmp;
  int ntri;
  int irrep_dim;
  std::string id;
  
  Ptot = init_array(natri);
  nso = nbfso;
  ntri = (nmo*(nmo+1))/2;

  if(opdm_basis == "MO" && ref == "RHF") {
    id == "RHF";
  }
  else if(opdm_basis == "MO" && ref == "ROHF") {
    /* Methods using semi-canonical orbitals */
    if(wfn == "MP2" || wfn == "CC2" ||wfn == "EOM_CC2" ||wfn == "CCSD(T)" ||wfn == "CC3" ||wfn == "EOM_CC3")
      id == "UHF";
    else id == "RHF";
  }
  else if(opdm_basis == "MO" && ref == "UHF") {
    id == "UHF";
  }

  /* Read OPDM in MO-basis */

  if(id == "RHF") {
    psq_so = block_matrix(nbfso,nbfso);

    max_opi = 0;
    for (irrep=0; irrep<nirreps; irrep++)
      if (sopi[irrep] > max_opi) max_opi = sopi[irrep];
    opdm_blk = block_matrix(max_opi, max_opi); 
    tmp_mat = block_matrix(max_opi, max_opi);

    docc = init_int_array(nirreps);
    socc = init_int_array(nirreps);
    ras_opi = init_int_matrix(4,nirreps);
    reorder = init_int_array(nmo);
    rstr_docc = init_int_array(nirreps);
    rstr_uocc = init_int_array(nirreps);

    fzc = 1;
    fzc = options.get_int("FREEZE_CORE");
 
    if (ci_wfn(wfn)) {
        frozen_docc = init_int_array(nirreps);
        frozen_uocc = init_int_array(nirreps);
      if (!ras_set2(nirreps, nmo, fzc, 1, orbspi, docc, socc,
		   frozen_docc, frozen_uocc, rstr_docc, rstr_uocc,
                   ras_opi, reorder, 1, 0) )
        throw PsiException("Error in ras_set()", __FILE__, __LINE__);
      /* treat restricted vir as frozen vir for now */
      for (i=0; i<nirreps; i++) {
        frozen_uocc[i] += rstr_uocc[i];
      }
    }
    else { /* CC densities */
      if (chkpt_rd_override_occ()) { /* ignore input occupations */
        docc = chkpt_rd_clsdpi();
        socc = chkpt_rd_openpi();
        frozen_docc = chkpt_rd_frzcpi();
        frozen_uocc = chkpt_rd_frzvpi();
      }
      else { /* try to read input occupations if you can */
        if(options["DOCC"].has_changed()) 
          docc = options.get_int_array("DOCC");
        else {
          free(docc);
          docc = chkpt_rd_clsdpi();
        }
                
        if(options["SOCC"].has_changed()) 
          socc = options.get_int_array("SOCC");
        else {
          free(socc);
          socc = chkpt_rd_clsdpi();
        }



        frozen_docc = get_frzcpi();
        frozen_uocc = get_frzvpi();
        }

      reorder_qt(docc, socc, frozen_docc, frozen_uocc,
               reorder, orbspi, nirreps);
    }

    if (print_lvl >= PRINTOPDMLEVEL) {
      fprintf(outfile, "Reorder array:\n");
      for (i=0; i<nmo; i++) fprintf(outfile, "%d ", reorder[i]);
      fprintf(outfile, "\n");
    }

    populated_orbs = nmo;
    for (irrep=0; irrep<nirreps; irrep++) {
      populated_orbs -= frozen_uocc[irrep];
    }
    onepdm = block_matrix(populated_orbs, populated_orbs);

    psio_open(opdm_file, PSIO_OPEN_OLD);
    psio_read_entry(opdm_file, opdm_lbl[irho], (char *) onepdm[0],
    populated_orbs * populated_orbs * sizeof(double));
    /* psio_read_entry(opdm_file, "MO-basis OPDM", (char *) onepdm[0],
    populated_orbs * populated_orbs * sizeof(double)); */
    psio_close(opdm_file, 1);

    if (print_lvl > 2) { 
      fprintf(outfile, "\n  Density matrix read");
      fprintf(outfile, " for label %s:\n", opdm_lbl[irho]);
      print_mat(onepdm,populated_orbs,populated_orbs,outfile);
      fprintf(outfile, "\n");
    }

    mo_offset = 0;
    so_offset = 0;
    for (irrep=0; irrep<nirreps; irrep++) {
      if (orbspi[irrep] == 0) continue;
      for (i=0; i<orbspi[irrep]-frozen_uocc[irrep]; i++) {
        for (j=0; j<orbspi[irrep]-frozen_uocc[irrep]; j++) {
          i_ci = reorder[i+mo_offset];
          j_ci = reorder[j+mo_offset];
          opdm_blk[i][j] = onepdm[i_ci][j_ci];
        }
      }

      /*if (print_lvl >= PRINTOPDMLEVEL) {
        fprintf(outfile, "Irrep %d (MO basis)\n", irrep);
        print_mat(opdm_blk,orbspi[irrep]-frozen_uocc[irrep], 
                  orbspi[irrep]-frozen_uocc[irrep],outfile);
      }*/

      scfvec = chkpt_rd_scf_irrep(irrep);

      mmult(opdm_blk,0,scfvec,1,tmp_mat,0,orbspi[irrep]-frozen_uocc[irrep],
            orbspi[irrep]-frozen_uocc[irrep],sopi[irrep],0);
      mmult(scfvec,0,tmp_mat,0,opdm_blk,0,sopi[irrep],
            orbspi[irrep]-frozen_uocc[irrep],sopi[irrep],0);
      /*if (print_lvl >= PRINTOPDMLEVEL) {
        fprintf(outfile,"Irrep %d (SO basis)\n", irrep);
        print_mat(opdm_blk,sopi[irrep], sopi[irrep],outfile);
      }*/
      for (i=0; i<sopi[irrep]; i++)
        for (j=0; j<sopi[irrep]; j++)
          psq_so[i+so_offset][j+so_offset] = opdm_blk[i][j];
      mo_offset += orbspi[irrep];
      so_offset += sopi[irrep];
    }

    /*if (print_lvl >= PRINTOPDMLEVEL) {
      fprintf(outfile,"  Total density matrix in SO basis :\n");
      print_mat(psq_so,nbfso,nbfso,outfile);
      fprintf(outfile,"\n");
    }*/
    
    free_block(onepdm);
    free_block(opdm_blk);
    free_block(tmp_mat);
    free_block(scfvec);
    free(docc); free(socc); free(frozen_docc); free(frozen_uocc);
    free(reorder);
    free_int_matrix(ras_opi);

    tmp_mat = init_matrix(nbfso,nbfao);
    psq_ao = init_matrix(nbfao,nbfao);
    mmult(psq_so,0,usotao,0,tmp_mat,0,nbfso,nbfso,nbfao,0);
    mmult(usotao,1,tmp_mat,0,psq_ao,0,nbfao,nbfso,nbfao,0);
    free_matrix(tmp_mat,nbfso);
    free_block(psq_so);

    /* Symmetrize an Asymmetric OPDM */
    if(asymm_opdm && opdm_format == "SQUARE") {
      for(i=0;i<nbfao;i++) {
        for(j=0;j<=i;j++) {
          Ptot[ioff[i]+j] = 0.5 * (psq_ao[i][j] + psq_ao[j][i]);
	}
      }
    }
    else {
      sq_to_tri(psq_ao,Ptot,nbfao);
    }

    free_matrix(psq_ao,nbfao);
  } /* end read RHF MO-basis case */
  else if(id == "UHF") {
    eval_tmp = chkpt_rd_alpha_evals();
    chkpt_wt_evals(eval_tmp);
    free(eval_tmp);

    /* One symmetry block at a time */
    maxmopi = 0;
    for (irrep=0; irrep<nirreps; irrep++) {
      if (sopi[irrep] > maxmopi) {
        maxmopi = sopi[irrep];
      }
    }
    
    reorder_a = init_int_array(nmo);
    reorder_b = init_int_array(nmo);

    mopi = chkpt_rd_orbspi();
    doccpi = chkpt_rd_clsdpi();
    soccpi = chkpt_rd_openpi();
    fzdoccpi = get_frzcpi();
    fzvirtpi = get_frzvpi();
    
    /* Pitzer to QTS ordering */
    reorder_qt_uhf(doccpi,soccpi,fzdoccpi,fzvirtpi,reorder_a,reorder_b,
      mopi,nirreps);
    
    aopdm = block_matrix(nmo,nmo);
    bopdm = block_matrix(nmo,nmo);
    
    psio_open(PSIF_MO_OPDM, PSIO_OPEN_OLD);
		psio_read_entry(PSIF_MO_OPDM, opdm_a_lbl[irho], (char *)aopdm[0],
                    sizeof(double)*nmo*nmo);
    psio_read_entry(PSIF_MO_OPDM, opdm_b_lbl[irho], (char *)bopdm[0],
                    sizeof(double)*nmo*nmo);
		/* psio_read_entry(PSIF_MO_OPDM, "MO-basis Alpha OPDM", (char *)aopdm[0],
                    sizeof(double)*nmo*nmo);
    psio_read_entry(PSIF_MO_OPDM, "MO-basis Beta OPDM", (char *)bopdm[0],
                    sizeof(double)*nmo*nmo); */
    psio_close(PSIF_MO_OPDM, 1);
    
    if (print_lvl >= PRINTOPDMLEVEL) {
      fprintf(outfile, "Alpha OPDM (MO)\n");
      print_mat(aopdm,nmo,nmo,outfile);
      fprintf(outfile, "Beta OPDM (MO)\n");
      print_mat(bopdm,nmo,nmo,outfile);
    }
 
    aPsq_so = block_matrix(nso,nso);
    bPsq_so = block_matrix(nso,nso);
    
    mo_offset = 0;
    so_offset = 0;
    
    aopdm_blk = block_matrix(maxmopi,maxmopi);
    bopdm_blk = block_matrix(maxmopi,maxmopi);
    tmp_mat = block_matrix(maxmopi,maxmopi);

    for (irrep=0; irrep<nirreps; irrep++) {
      zero_mat(aopdm_blk,maxmopi,maxmopi);
      zero_mat(bopdm_blk,maxmopi,maxmopi);
      zero_mat(tmp_mat,maxmopi,maxmopi);

      if (mopi[irrep] == 0) {
        continue;
      }
      
      for (i=0; i<mopi[irrep]; i++) {
        for (j=0; j<mopi[irrep]; j++) {
	  aopdm_blk[i][j] = 
            aopdm[reorder_a[i+mo_offset]][reorder_a[j+mo_offset]];
	  bopdm_blk[i][j] = 
            bopdm[reorder_b[i+mo_offset]][reorder_b[j+mo_offset]];
	}
      }
      
      /*fprintf(outfile,"Alpha OPDM Irrep %d (MO)\n",irrep);
      print_mat(aopdm_blk,mopi[irrep],mopi[irrep],outfile);

      fprintf(outfile,"Beta ODPM Irrep %d (MO)\n",irrep);
      print_mat(bopdm_blk,mopi[irrep],mopi[irrep],outfile);*/
      
      ascfvec = chkpt_rd_alpha_scf_irrep(irrep);
      bscfvec = chkpt_rd_beta_scf_irrep(irrep);
      
      if (print_lvl >= PRINTOPDMLEVEL) {
        fprintf(outfile, "Alpha SCF Vector\n");
        print_mat(ascfvec, sopi[irrep], mopi[irrep], outfile);
      
        fprintf(outfile, "Beta SCF Vector\n");
        print_mat(bscfvec, sopi[irrep], mopi[irrep], outfile);
      }
 
      mmult(aopdm_blk,0,ascfvec,1,tmp_mat,0,
            mopi[irrep],mopi[irrep],sopi[irrep],0);
      mmult(ascfvec,0,tmp_mat,0,aopdm_blk,0,
            sopi[irrep],mopi[irrep],sopi[irrep],0);
      
      zero_mat(tmp_mat,maxmopi,maxmopi);

      mmult(bopdm_blk,0,bscfvec,1,tmp_mat,0,
            mopi[irrep],mopi[irrep],sopi[irrep],0);
      mmult(bscfvec,0,tmp_mat,0,bopdm_blk,0,
            sopi[irrep],mopi[irrep],sopi[irrep],0);
     
      /*fprintf(outfile,"Alpha OPDM Irrep %d (SO)\n",irrep);
      print_mat(aopdm_blk,sopi[irrep],sopi[irrep],outfile);
      
      fprintf(outfile,"Beta ODPM Irrep %d (SO)\n",irrep);
      print_mat(bopdm_blk,sopi[irrep],sopi[irrep],outfile);*/
      
      for (i=0; i<sopi[irrep]; i++) {
        for (j=0; j<sopi[irrep]; j++) {
	  aPsq_so[i+so_offset][j+so_offset] = aopdm_blk[i][j];
	  bPsq_so[i+so_offset][j+so_offset] = bopdm_blk[i][j];
	}
      }
      
      mo_offset += mopi[irrep];
      so_offset += sopi[irrep];
      
      free_block(ascfvec);
      free_block(bscfvec);
    }

    free_block(tmp_mat);
    
    if (print_lvl >= PRINTOPDMLEVEL) {
      fprintf(outfile,"Alpha OPDM (SO)\n");
      print_mat(aPsq_so,nso,nso,outfile);

      fprintf(outfile,"Beta ODPM (SO)\n");
      print_mat(bPsq_so,nso,nso,outfile);
    }
   
    free(reorder_a);
    free(reorder_b);
    free_block(aopdm);
    free_block(aopdm_blk);
    free_block(bopdm);
    free_block(bopdm_blk);
    
    P_so_tot = block_matrix(nso,nso);
    
    for (i=0; i<nso; i++) {
      for (j=0; j<nso; j++) {
        P_so_tot[i][j] = aPsq_so[i][j] + bPsq_so[i][j];
      }
    }
    
    /* 
       Two Cases:
       1.  Properties (AO)
       2.  Natural Orbitals (MO)
    */
    
    if(wrtnos) {
      fprintf(outfile,"OPDM (SO)\n");
      print_mat(P_so_tot,nso,nso,outfile);

      P_mo_tot = block_matrix(nmo, nmo);
      tmat = block_matrix(nmo, nmo);
      scfvec = chkpt_rd_alpha_scf();
  
      Stri = init_array(ntri);
      Smat = block_matrix(nmo, nmo);
      iwl_rdone(PSIF_OEI, PSIF_SO_S, Stri, ntri, 0, 0, outfile);
      tri_to_sq(Stri, Smat, nmo);
      free(Stri);
 
      /*  What is this? */

      mmult(P_so_tot, 0, Smat, 0, tmat, 0, nmo, nmo, nmo, 0);
      mmult(tmat, 0, scfvec, 0, P_so_tot, 0, nmo, nmo, nmo, 0);
      mmult(Smat, 0, P_so_tot, 0, tmat, 0, nmo, nmo, nmo, 0);
      mmult(scfvec, 1, tmat, 0, P_mo_tot, 0, nmo, nmo, nmo, 0);
   
      fprintf(outfile, "OPM (alpha MO) \n");
      print_mat(P_mo_tot, nmo, nmo, outfile); 

      P_eigvals = init_array(nmo);
      P_eigvecs = block_matrix(nmo,nmo);
      P_mo_block = block_matrix(nmo,nmo);
      P_so_block = block_matrix(nmo,nmo);

      chkpt_wt_scf(P_mo_block);

      for (irrep=0,mo_offset=0; irrep<nirreps; irrep++) {
        irrep_dim = mopi[irrep];
        if (irrep_dim == 0) {
          continue;
        }
        scfvec_irrep = chkpt_rd_alpha_scf_irrep(irrep);
        for (i=0; i<irrep_dim; i++) {
          for (j=0; j<irrep_dim; j++) {
            i_ci = i+mo_offset;
	    j_ci = j+mo_offset;
	    P_mo_block[i][j] = P_mo_tot[i_ci][j_ci];
          }
        }
    
        fprintf(outfile, "%d Block of OPDM\n",irrep);
        print_mat(P_mo_block, irrep_dim, irrep_dim, outfile);
        fprintf(outfile, "\n");

        sq_rsp(irrep_dim, irrep_dim, P_mo_block, P_eigvals, 3, P_eigvecs, TOL);

        fprintf(outfile, "%d NOs in terms of molecular orbitals\n",irrep);
        eivout(P_eigvecs, P_eigvals, irrep_dim, irrep_dim, outfile);

        mmult(scfvec_irrep, 0, P_eigvecs, 0, P_so_block, 0, irrep_dim, 
              irrep_dim, irrep_dim, 0);

        fprintf(outfile, "%d NOs in terms of symmetry orbitals\n",irrep);
        eivout(P_so_block, P_eigvals, irrep_dim, irrep_dim, outfile);

        chkpt_wt_scf_irrep(P_so_block,irrep);
        free_block(scfvec_irrep);
        mo_offset += mopi[irrep];
      }
      
      fprintf(outfile, "NOs written to checkpoint file\n\n");
      chkpt_close();

      free_block(tmat);
      free_block(Smat);
      free_block(scfvec);
      free_block(P_mo_tot);
      free_block(P_mo_block);
      free_block(P_so_tot);
      free_block(P_so_block);
      free_block(P_eigvecs);
      free(P_eigvals);
    }
    else if(wrtnos==0) {
      tmp_mat = init_matrix(nbfso,nbfao);
      psq_ao = init_matrix(nbfao,nbfao);
      mmult(P_so_tot,0,usotao,0,tmp_mat,0,nbfso,nbfso,nbfao,0);
      mmult(usotao,1,tmp_mat,0,psq_ao,0,nbfao,nbfso,nbfao,0);
      free_matrix(tmp_mat,nbfso);
      free_block(P_so_tot);

      if(asymm_opdm && opdm_format == "SQUARE") {
        for(i=0;i<nbfao;i++) {
          for(j=0;j<=i;j++) {
            Ptot[ioff[i]+j] = 0.5 * (psq_ao[i][j] + psq_ao[j][i]);
	  }
        }
      }
      else {
        sq_to_tri(psq_ao,Ptot,nbfao);
      }

      free_matrix(psq_ao,nbfao);
    }
  }

  /* Read OPDM in SO-basis */
  if(opdm_basis == "SO") {
    if(opdm_format == "SQUARE") {
      psq_so = block_matrix(nbfso,nbfso);
      psio_open(opdm_file, PSIO_OPEN_OLD);
      psio_read_entry(opdm_file, "SO-basis OPDM", (char *) psq_so[0],
	  nbfso*nbfso*sizeof(double));
      psio_close(opdm_file, 1);
    }
    else {
      tmp_arr = init_array(nstri);
      psq_so = init_matrix(nbfso,nbfso);
      psio_open(opdm_file, PSIO_OPEN_OLD);
      psio_read_entry(opdm_file, "SO-basis OPDM", (char *) tmp_arr,
	  nstri*sizeof(double));
      psio_close(opdm_file, 1);
      tri_to_sq(tmp_arr,psq_so,nbfso);
      free(tmp_arr);
    }
          
    tmp_mat = init_matrix(nbfso,nbfao);
    psq_ao = init_matrix(nbfao,nbfao);
    mmult(psq_so,0,usotao,0,tmp_mat,0,nbfso,nbfso,nbfao,0);
    mmult(usotao,1,tmp_mat,0,psq_ao,0,nbfao,nbfso,nbfao,0);
    free_matrix(tmp_mat,nbfso);

    if(asymm_opdm && opdm_format == "SQUARE") {
      for(i=0; i<nbfao; i++) {
        for(j=0; j<=i; j++) {
          Ptot[ioff[i]+j] = 0.5 * (psq_ao[i][j] + psq_ao[j][i]);
	}
      }
    }
    else {
      sq_to_tri(psq_ao,Ptot,nbfao);
    }

    free_matrix(psq_ao,nbfao);
    if (print_lvl >= PRINTOPDMLEVEL) {
      fprintf(outfile,"  Total density matrix in SO basis :\n");
      print_mat(psq_so,nbfso,nbfso,outfile);
      fprintf(outfile,"\n");
    }
    free_block(psq_so);
  }
  /* Read in OPDM in AO-basis */
  else if(opdm_basis == "AO") {
    if(opdm_format == "SQUARE") {
      psq_ao = block_matrix(nbfao,nbfao);
      psio_open(opdm_file, PSIO_OPEN_OLD);
      psio_read_entry(opdm_file, "AO-basis OPDM", (char *) psq_ao[0],
	  nbfao*nbfao*sizeof(double));
      psio_close(opdm_file, 1);
      /* Antisymmetric OPDM */
      if (asymm_opdm) {
        for(i=0; i<nbfao; i++) {
          for(j=0; j<=i; j++) {
            Ptot[ioff[i]+j] = (psq_ao[i][j] + psq_ao[j][i])/2;
	  }
	}
      }
      /* Symmetric OPDM */
      else {
        sq_to_tri(psq_ao,Ptot,nbfao);
      }
      
      free_matrix(psq_ao,nbfao);
    }
    /* There is nothing for AO-basis in Diagonal form */
    else {
      fprintf(outfile, "Unable to read the OPDM\n");
      abort();
    }
  }

  if (print_lvl >= PRINTOPDMLEVEL) {
    fprintf(outfile,"  Total density matrix in AO basis :\n");
    print_array(Ptot,nbfao,outfile);
    fprintf(outfile,"\n\n");
  }

  return;
}

}} // namespace psi::oeprop
