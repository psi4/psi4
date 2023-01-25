/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2023 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include "ccsd.h"
#include "blas.h"

#include "psi4/psi4-dec.h"
#include "psi4/psifiles.h"
#include "psi4/libiwl/iwl.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libqt/qt.h"
#include "psi4/libpsi4util/process.h"

namespace psi {
namespace fnocc {
typedef long int size_t;
struct integral {
    size_t ind;
    double val;
};
void SortAllIntegrals(iwlbuf *Buf, int nfzc, int nfzv, int norbs, int ndoccact, int nvirt, Options &options);
void klcd_terms_incore(double val, size_t pq, size_t rs, size_t p, size_t q, size_t r, size_t s, size_t o, size_t v,
                       double *klcd);
void ijkl_terms(double val, size_t pq, size_t rs, size_t p, size_t q, size_t r, size_t s, size_t o, size_t &nijkl,
                struct integral *ijkl);
void ijak_terms(double val, size_t p, size_t q, size_t r, size_t s, size_t o, size_t v, size_t &nijak,
                struct integral *ijak);
void ijak2_terms(double val, size_t p, size_t q, size_t r, size_t s, size_t o, size_t v, size_t &nijak2,
                 struct integral *ijak2);
void klcd_terms(double val, size_t pq, size_t rs, size_t p, size_t q, size_t r, size_t s, size_t o, size_t v,
                size_t &nklcd, struct integral *klcd);
void akjc_terms(double val, size_t p, size_t q, size_t r, size_t s, size_t o, size_t v, size_t &nklcd,
                struct integral *klcd);
void abci1_terms(double val, size_t p, size_t q, size_t r, size_t s, size_t o, size_t v, size_t &nabci1,
                 struct integral *abci1);
void abci3_terms(double val, size_t p, size_t q, size_t r, size_t s, size_t o, size_t v, size_t &nabci3,
                 struct integral *abci3);
void abci4_terms(double val, size_t p, size_t q, size_t r, size_t s, size_t o, size_t v, size_t &nabci4,
                 struct integral *abci4);
void abci5_terms(double val, size_t p, size_t q, size_t r, size_t s, size_t o, size_t v, size_t &nabci5,
                 struct integral *abci5);
void abcd1_terms(double val, size_t pq, size_t rs, size_t p, size_t q, size_t r, size_t s, size_t o, size_t v,
                 size_t &nabcd1, struct integral *abcd1);
void abcd2_terms(double val, size_t pq, size_t rs, size_t p, size_t q, size_t r, size_t s, size_t o, size_t v,
                 size_t &nabcd2, struct integral *abcd2);

void abcd1_terms_new(double val, size_t pq, size_t rs, size_t p, size_t q, size_t r, size_t s, size_t o, size_t v,
                     size_t *nabcd1, size_t *totalnabcd1, struct integral **abcd1, size_t binsize, size_t bucketsize,
                     psio_address *addr, size_t nfiles);
void abcd2_terms_new(double val, size_t pq, size_t rs, size_t p, size_t q, size_t r, size_t s, size_t o, size_t v,
                     size_t *nabcd2, size_t *totalnabcd2, struct integral **abcd2, size_t binsize, size_t bucketsize,
                     psio_address *addr, size_t nfiles);

void abci1_terms_new(double val, size_t p, size_t q, size_t r, size_t s, size_t o, size_t v, size_t *nabci1,
                     size_t *totalnabci1, struct integral **abci1, size_t binsize, size_t bucketsize,
                     psio_address *addr, size_t filestart, size_t nfiles);
void abci3_terms_new(double val, size_t p, size_t q, size_t r, size_t s, size_t o, size_t v, size_t *nabci3,
                     size_t *totalnabci3, struct integral **abci3, size_t binsize, size_t bucketsize,
                     psio_address *addr, size_t filestart, size_t nfiles);
void abci5_terms_new(double val, size_t p, size_t q, size_t r, size_t s, size_t o, size_t v, size_t *nabci5,
                     size_t *totalnabci5, struct integral **abci5, size_t binsize, size_t bucketsize,
                     psio_address *addr, size_t filestart, size_t nfiles);

void abcd3_terms(double val, size_t pq, size_t rs, size_t p, size_t q, size_t r, size_t s, size_t o, size_t v,
                 size_t &nabcd1, struct integral *abcd1);
void SortBlock(size_t nelem, size_t blockdim, struct integral *buffer, double *tmp, size_t PSIFILE, const char *string,
               size_t maxdim);
void SortBlockNew(size_t nelem, size_t blockdim, struct integral *buffer, double *tmp, size_t PSIFILE,
                  const char *string, size_t maxdim);
void SortBlockNewNew(size_t *nelem, size_t blockdim, struct integral *buffer, double *tmp, size_t PSIFILE,
                     const char *string, size_t maxdim, size_t filestart, size_t nfiles);
}
}

namespace psi {
namespace fnocc {
void SortIntegrals(int nfzc, int nfzv, int norbs, int ndoccact, int nvirt, Options &options) {
    struct iwlbuf Buf;
    iwl_buf_init(&Buf, PSIF_MO_TEI, 0.0, 1, 1);

    outfile->Printf("\n");
    outfile->Printf("        **********************************************************\n");
    outfile->Printf("        *                                                        *\n");
    outfile->Printf("        *                   CCSD Integral Sort                   *\n");
    outfile->Printf("        *                                                        *\n");
    outfile->Printf("        **********************************************************\n");
    outfile->Printf("\n");
    outfile->Printf("\n");

    SortAllIntegrals(&Buf, nfzc, nfzv, norbs, ndoccact, nvirt, options);

    iwl_buf_close(&Buf, 1);
}
void SortAllIntegrals(iwlbuf *Buf, int nfzc, int nfzv, int norbs, int ndoccact, int nvirt, Options &options) {
    double val;
    size_t o = ndoccact;
    size_t v = nvirt;
    size_t fstact = nfzc;
    size_t lstact = norbs - nfzv;

    size_t lastbuf;
    Label *lblptr;
    Value *valptr;
    size_t nocc, idx, p, q, r, s, pq, rs, pqrs;

    lblptr = Buf->labels;
    valptr = Buf->values;

    lastbuf = Buf->lastbuf;

    // buckets for integrals:
    struct integral *ijkl, *klcd, *akjc, **abci1, **abci3;
    struct integral **abci5, **abcd1, **abcd2, *ijak, *ijak2;

    // available memory:
    size_t memory = Process::environment.get_memory();

    // 8 bytes for tmp, 16 for integral struct
    size_t maxelem = memory / (sizeof(double) + sizeof(struct integral));

    // what is the biggest block?
    size_t maxblock = o * o * o * o;
    if (maxblock < o * o * o * v) maxblock = o * o * o * v;
    if (maxblock < o * o * v * v) maxblock = o * o * v * v;
    if (maxblock < o * v * v * v) maxblock = o * v * v * v;
    if (maxblock < 2 * v * (v + 1) / 2 * v * (v + 1) / 2) maxblock = 2 * v * (v + 1) / 2 * v * (v + 1) / 2;

    // maxelem should be, at max, the size of the biggest block
    if (maxelem > maxblock) maxelem = maxblock;

    outfile->Printf("        CC integral sort will use                   %7.2lf mb\n",
                    maxelem * (sizeof(double) + sizeof(struct integral)) / 1024. / 1024.);
    if (maxelem < 2 * v * (v + 1) / 2 * v * (v + 1) / 2) {
        outfile->Printf("       (for most efficient sort, increase memory by %7.2lf mb)\n",
                        (2 * v * (v + 1) / 2 * v * (v + 1) / 2 - maxelem) * (sizeof(double) + sizeof(struct integral)) /
                            1024. / 1024.);
    }
    outfile->Printf("\n");

    // how many files does (ac|bd) need to be?
    size_t filesize;
    size_t vtri = v * (v + 1L) / 2L;
    size_t nfiles = 0;
    for (size_t i = 1; i <= vtri * vtri; i++) {
        if (maxelem >= vtri * vtri / i) {
            filesize = vtri * vtri / i;
            if (i * filesize < vtri * vtri) filesize++;
            nfiles = i;
            break;
        }
    }
    if (nfiles == 0) throw PsiException("how is this possible? (ab|cd)", __FILE__, __LINE__);

    // how many (ab|ci) files?
    size_t ov3filesize;
    size_t ov3nfiles = 0;
    size_t ov3 = o * v * v * v;
    for (size_t i = 1; i <= ov3; i++) {
        if (maxelem >= ov3 / i) {
            ov3filesize = ov3 / i;
            if (i * ov3filesize < ov3) ov3filesize++;
            ov3nfiles = i;
            break;
        }
    }
    if (ov3nfiles == 0) throw PsiException("how is this possible (ab|ci)?", __FILE__, __LINE__);

    struct integral *integralbuffer;
    size_t nelem = maxelem / (5 + 3 * ov3nfiles + 2 * nfiles) - 20;
    if ((nelem + 20) * (5 + 3 * ov3nfiles + 2 * nfiles) > maxelem)
        integralbuffer = new integral[(nelem + 20) * (5 + 3 * ov3nfiles + 2 * nfiles)];
    else
        integralbuffer = new integral[maxelem];

    outfile->Printf("        Number of (ab|cd) temporary files:            %5li\n", 2 * nfiles);
    outfile->Printf("        Number of (ab|ci) temporary files:            %5li\n", 3 * ov3nfiles);
    outfile->Printf("        Starting temporary file number:               %5i\n", PSIF_DCC_SORT_START);
    outfile->Printf("\n");

    // buckets:
    size_t bucketsize = nelem;

    ijkl = integralbuffer;
    ijak = integralbuffer + (nelem + 20);
    klcd = integralbuffer + (nelem + 20) * 2;
    akjc = integralbuffer + (nelem + 20) * 3;
    ijak2 = integralbuffer + (nelem + 20) * 4;

    // abci1 = integralbuffer+(nelem+20)*5;
    // abci3 = integralbuffer+(nelem+20)*6;
    // abci5 = integralbuffer+(nelem+20)*7;

    abci1 = (struct integral **)malloc(ov3nfiles * sizeof(struct integral *));
    abci3 = (struct integral **)malloc(ov3nfiles * sizeof(struct integral *));
    abci5 = (struct integral **)malloc(ov3nfiles * sizeof(struct integral *));
    for (size_t k = 0; k < ov3nfiles; k++) {
        abci1[k] = integralbuffer + (nelem + 20L) * (5L + k);
        abci3[k] = integralbuffer + (nelem + 20L) * (5L + 1L * ov3nfiles + k);
        abci5[k] = integralbuffer + (nelem + 20L) * (5L + 2L * ov3nfiles + k);
    }

    abcd1 = (struct integral **)malloc(nfiles * sizeof(struct integral *));
    abcd2 = (struct integral **)malloc(nfiles * sizeof(struct integral *));
    for (size_t k = 0; k < nfiles; k++) {
        abcd1[k] = integralbuffer + (nelem + 20L) * (5L + 3L * ov3nfiles + k);
        abcd2[k] = integralbuffer + (nelem + 20L) * (nfiles + 5L + 3L * ov3nfiles + k);
    }

    auto psio = std::make_shared<PSIO>();

    psio_address ijkl_addr = PSIO_ZERO;
    psio_address klcd_addr = PSIO_ZERO;
    psio_address akjc_addr = PSIO_ZERO;
    psio_address ijak_addr = PSIO_ZERO;
    psio_address ijak2_addr = PSIO_ZERO;
    psio_address abci2_addr = PSIO_ZERO;

    auto *abci1_addr = new psio_address[ov3nfiles];
    auto *abci3_addr = new psio_address[ov3nfiles];
    auto *abci5_addr = new psio_address[ov3nfiles];
    auto *abcd1_addr = new psio_address[nfiles];
    auto *abcd2_addr = new psio_address[nfiles];
    for (size_t k = 0; k < ov3nfiles; k++) {
        abci1_addr[k] = PSIO_ZERO;
        abci3_addr[k] = PSIO_ZERO;
        abci5_addr[k] = PSIO_ZERO;
    }
    for (size_t k = 0; k < nfiles; k++) {
        abcd1_addr[k] = PSIO_ZERO;
        abcd2_addr[k] = PSIO_ZERO;
    }

    psio->open(PSIF_DCC_IJKL, PSIO_OPEN_NEW);
    psio->close(PSIF_DCC_IJKL, 1);
    psio->open(PSIF_DCC_IAJB, PSIO_OPEN_NEW);
    psio->close(PSIF_DCC_IAJB, 1);
    psio->open(PSIF_DCC_IJAK, PSIO_OPEN_NEW);
    psio->close(PSIF_DCC_IJAK, 1);
    psio->open(PSIF_DCC_IJAK2, PSIO_OPEN_NEW);
    psio->close(PSIF_DCC_IJAK2, 1);
    psio->open(PSIF_DCC_ABCI, PSIO_OPEN_NEW);
    psio->close(PSIF_DCC_ABCI, 1);
    psio->open(PSIF_DCC_ABCI2, PSIO_OPEN_NEW);
    psio->close(PSIF_DCC_ABCI2, 1);
    psio->open(PSIF_DCC_ABCI3, PSIO_OPEN_NEW);
    psio->close(PSIF_DCC_ABCI3, 1);
    psio->open(PSIF_DCC_ABCD1, PSIO_OPEN_NEW);
    psio->close(PSIF_DCC_ABCD1, 1);
    psio->open(PSIF_DCC_ABCD2, PSIO_OPEN_NEW);
    psio->close(PSIF_DCC_ABCD2, 1);
    psio->open(PSIF_DCC_IJAB, PSIO_OPEN_NEW);
    psio->close(PSIF_DCC_IJAB, 1);
    size_t nijkl = 0;
    size_t totalnijkl = 0;
    size_t nijak2 = 0;
    size_t totalnijak2 = 0;
    size_t nijak = 0;
    size_t totalnijak = 0;
    size_t nklcd = 0;
    size_t totalnklcd = 0;
    size_t nakjc = 0;
    size_t totalnakjc = 0;
    auto *nabci1 = new size_t[ov3nfiles];
    auto *totalnabci1 = new size_t[ov3nfiles];
    auto *nabci3 = new size_t[ov3nfiles];
    auto *totalnabci3 = new size_t[ov3nfiles];
    auto *nabci5 = new size_t[ov3nfiles];
    auto *totalnabci5 = new size_t[ov3nfiles];
    auto *nabcd1 = new size_t[nfiles];
    auto *totalnabcd1 = new size_t[nfiles];
    auto *nabcd2 = new size_t[nfiles];
    auto *totalnabcd2 = new size_t[nfiles];
    for (size_t k = 0; k < ov3nfiles; k++) {
        nabci1[k] = 0;
        nabci3[k] = 0;
        nabci5[k] = 0;
        totalnabci1[k] = 0;
        totalnabci3[k] = 0;
        totalnabci5[k] = 0;
    }
    for (size_t k = 0; k < nfiles; k++) {
        nabcd1[k] = 0;
        nabcd2[k] = 0;
        totalnabcd1[k] = 0;
        totalnabcd2[k] = 0;
    }

    // (ab|cd) files:
    for (size_t k = 0; k < nfiles; k++) {
        psio->open(PSIF_DCC_SORT_START + k, PSIO_OPEN_NEW);
        psio->open(PSIF_DCC_SORT_START + k + nfiles, PSIO_OPEN_NEW);
        psio->close(PSIF_DCC_SORT_START + k, 1);
        psio->close(PSIF_DCC_SORT_START + k + nfiles, 1);
    }
    // (ab|ci) files come after (ab|cd) files:
    for (size_t k = 0; k < nfiles; k++) {
        psio->open(PSIF_DCC_SORT_START + k + 2 * nfiles, PSIO_OPEN_NEW);
        psio->open(PSIF_DCC_SORT_START + k + 2 * nfiles + ov3nfiles, PSIO_OPEN_NEW);
        psio->open(PSIF_DCC_SORT_START + k + 2 * nfiles + 2 * ov3nfiles, PSIO_OPEN_NEW);
        psio->close(PSIF_DCC_SORT_START + k + 2 * nfiles, 1);
        psio->close(PSIF_DCC_SORT_START + k + 2 * nfiles + ov3nfiles, 1);
        psio->close(PSIF_DCC_SORT_START + k + 2 * nfiles + 2 * ov3nfiles, 1);
    }

    outfile->Printf("        Initial sort........");
    /**
      * first buffer (read in when Buf was initialized)
      */
    for (idx = 4 * Buf->idx; Buf->idx < Buf->inbuf; Buf->idx++) {
        p = (size_t)lblptr[idx++];
        q = (size_t)lblptr[idx++];
        r = (size_t)lblptr[idx++];
        s = (size_t)lblptr[idx++];

        if (p < fstact || q < fstact || r < fstact || s < fstact) continue;
        if (p > lstact || q > lstact || r > lstact || s > lstact) continue;
        p -= fstact;
        q -= fstact;
        r -= fstact;
        s -= fstact;

        pq = Position(p, q);
        rs = Position(r, s);
        pqrs = Position(pq, rs);

        nocc = 0;
        if (p < o) nocc++;
        if (q < o) nocc++;
        if (r < o) nocc++;
        if (s < o) nocc++;

        // which type of integral?

        if (nocc == 4) {
            val = (double)valptr[Buf->idx];
            ijkl_terms(val, pq, rs, p, q, r, s, o, nijkl, ijkl);

            if (nijkl >= nelem) {
                psio->open(PSIF_DCC_IJKL, PSIO_OPEN_OLD);
                psio->write(PSIF_DCC_IJKL, "E2ijkl", (char *)&ijkl[0], nijkl * sizeof(struct integral), ijkl_addr,
                            &ijkl_addr);
                psio->close(PSIF_DCC_IJKL, 1);
                totalnijkl += nijkl;
                nijkl = 0;
            }
        } else if (nocc == 3) {
            val = (double)valptr[Buf->idx];
            ijak_terms(val, p, q, r, s, o, v, nijak, ijak);
            if (nijak >= nelem) {
                psio->open(PSIF_DCC_IJAK, PSIO_OPEN_OLD);
                psio->write(PSIF_DCC_IJAK, "E2ijak", (char *)&ijak[0], nijak * sizeof(struct integral), ijak_addr,
                            &ijak_addr);
                psio->close(PSIF_DCC_IJAK, 1);
                totalnijak += nijak;
                nijak = 0;
            }
            ijak2_terms(val, p, q, r, s, o, v, nijak2, ijak2);
            if (nijak2 >= nelem) {
                psio->open(PSIF_DCC_IJAK2, PSIO_OPEN_OLD);
                psio->write(PSIF_DCC_IJAK2, "E2ijak2", (char *)&ijak2[0], nijak2 * sizeof(struct integral), ijak2_addr,
                            &ijak2_addr);
                psio->close(PSIF_DCC_IJAK2, 1);
                totalnijak2 += nijak2;
                nijak2 = 0;
            }
        } else if (nocc == 2) {
            val = (double)valptr[Buf->idx];

            if ((p < o && q >= o) || (p >= o && q < o)) {
                klcd_terms(val, pq, rs, p, q, r, s, o, v, nklcd, klcd);

                if (nklcd >= nelem) {
                    psio->open(PSIF_DCC_IAJB, PSIO_OPEN_OLD);
                    psio->write(PSIF_DCC_IAJB, "E2iajb", (char *)&klcd[0], nklcd * sizeof(struct integral), klcd_addr,
                                &klcd_addr);
                    psio->close(PSIF_DCC_IAJB, 1);
                    totalnklcd += nklcd;
                    nklcd = 0;
                }
            } else {
                akjc_terms(val, p, q, r, s, o, v, nakjc, akjc);
                if (nakjc >= nelem) {
                    psio->open(PSIF_DCC_IJAB, PSIO_OPEN_OLD);
                    psio->write(PSIF_DCC_IJAB, "E2ijab", (char *)&akjc[0], nakjc * sizeof(struct integral), akjc_addr,
                                &akjc_addr);
                    psio->close(PSIF_DCC_IJAB, 1);
                    totalnakjc += nakjc;
                    nakjc = 0;
                }
            }
        } else if (nocc == 1) {
            val = (double)valptr[Buf->idx];
            abci1_terms_new(val, p, q, r, s, o, v, nabci1, totalnabci1, abci1, ov3filesize, bucketsize, abci1_addr,
                            PSIF_DCC_SORT_START + 2 * nfiles, ov3nfiles);
            abci3_terms_new(val, p, q, r, s, o, v, nabci3, totalnabci3, abci3, ov3filesize, bucketsize, abci3_addr,
                            PSIF_DCC_SORT_START + 2 * nfiles + ov3nfiles, ov3nfiles);
            abci5_terms_new(val, p, q, r, s, o, v, nabci5, totalnabci5, abci5, ov3filesize, bucketsize, abci5_addr,
                            PSIF_DCC_SORT_START + 2 * nfiles + 2 * ov3nfiles, ov3nfiles);
        } else if (nocc == 0) {
            val = (double)valptr[Buf->idx];
            abcd1_terms_new(val, pq, rs, p, q, r, s, o, v, nabcd1, totalnabcd1, abcd1, filesize, bucketsize, abcd1_addr,
                            nfiles);
            abcd2_terms_new(val, pq, rs, p, q, r, s, o, v, nabcd2, totalnabcd2, abcd2, filesize, bucketsize, abcd2_addr,
                            nfiles);
        }
    }

    /**
      * now do the same for the rest of the buffers
      */
    while (!lastbuf) {
        iwl_buf_fetch(Buf);
        lastbuf = Buf->lastbuf;
        for (idx = 4 * Buf->idx; Buf->idx < Buf->inbuf; Buf->idx++) {
            p = (size_t)lblptr[idx++];
            q = (size_t)lblptr[idx++];
            r = (size_t)lblptr[idx++];
            s = (size_t)lblptr[idx++];

            if (p < fstact || q < fstact || r < fstact || s < fstact) continue;
            if (p > lstact || q > lstact || r > lstact || s > lstact) continue;
            p -= fstact;
            q -= fstact;
            r -= fstact;
            s -= fstact;

            pq = Position(p, q);
            rs = Position(r, s);
            pqrs = Position(pq, rs);

            nocc = 0;
            if (p < o) nocc++;
            if (q < o) nocc++;
            if (r < o) nocc++;
            if (s < o) nocc++;

            // which type of integral?

            if (nocc == 4) {
                val = (double)valptr[Buf->idx];
                ijkl_terms(val, pq, rs, p, q, r, s, o, nijkl, ijkl);

                if (nijkl >= nelem) {
                    psio->open(PSIF_DCC_IJKL, PSIO_OPEN_OLD);
                    psio->write(PSIF_DCC_IJKL, "E2ijkl", (char *)&ijkl[0], nijkl * sizeof(struct integral), ijkl_addr,
                                &ijkl_addr);
                    psio->close(PSIF_DCC_IJKL, 1);
                    totalnijkl += nijkl;
                    nijkl = 0;
                }
            } else if (nocc == 3) {
                val = (double)valptr[Buf->idx];
                ijak_terms(val, p, q, r, s, o, v, nijak, ijak);
                if (nijak >= nelem) {
                    psio->open(PSIF_DCC_IJAK, PSIO_OPEN_OLD);
                    psio->write(PSIF_DCC_IJAK, "E2ijak", (char *)&ijak[0], nijak * sizeof(struct integral), ijak_addr,
                                &ijak_addr);
                    psio->close(PSIF_DCC_IJAK, 1);
                    totalnijak += nijak;
                    nijak = 0;
                }
                ijak2_terms(val, p, q, r, s, o, v, nijak2, ijak2);
                if (nijak2 >= nelem) {
                    psio->open(PSIF_DCC_IJAK2, PSIO_OPEN_OLD);
                    psio->write(PSIF_DCC_IJAK2, "E2ijak2", (char *)&ijak2[0], nijak2 * sizeof(struct integral),
                                ijak2_addr, &ijak2_addr);
                    psio->close(PSIF_DCC_IJAK2, 1);
                    totalnijak2 += nijak2;
                    nijak2 = 0;
                }
            } else if (nocc == 2) {
                val = (double)valptr[Buf->idx];

                if ((p < o && q >= o) || (p >= o && q < o)) {
                    klcd_terms(val, pq, rs, p, q, r, s, o, v, nklcd, klcd);

                    if (nklcd >= nelem) {
                        psio->open(PSIF_DCC_IAJB, PSIO_OPEN_OLD);
                        psio->write(PSIF_DCC_IAJB, "E2iajb", (char *)&klcd[0], nklcd * sizeof(struct integral),
                                    klcd_addr, &klcd_addr);
                        psio->close(PSIF_DCC_IAJB, 1);
                        totalnklcd += nklcd;
                        nklcd = 0;
                    }
                } else {
                    akjc_terms(val, p, q, r, s, o, v, nakjc, akjc);

                    if (nakjc >= nelem) {
                        psio->open(PSIF_DCC_IJAB, PSIO_OPEN_OLD);
                        psio->write(PSIF_DCC_IJAB, "E2ijab", (char *)&akjc[0], nakjc * sizeof(struct integral),
                                    akjc_addr, &akjc_addr);
                        psio->close(PSIF_DCC_IJAB, 1);
                        totalnakjc += nakjc;
                        nakjc = 0;
                    }
                }
            } else if (nocc == 1) {
                val = (double)valptr[Buf->idx];
                abci1_terms_new(val, p, q, r, s, o, v, nabci1, totalnabci1, abci1, ov3filesize, bucketsize, abci1_addr,
                                PSIF_DCC_SORT_START + 2 * nfiles, ov3nfiles);
                abci3_terms_new(val, p, q, r, s, o, v, nabci3, totalnabci3, abci3, ov3filesize, bucketsize, abci3_addr,
                                PSIF_DCC_SORT_START + 2 * nfiles + ov3nfiles, ov3nfiles);
                abci5_terms_new(val, p, q, r, s, o, v, nabci5, totalnabci5, abci5, ov3filesize, bucketsize, abci5_addr,
                                PSIF_DCC_SORT_START + 2 * nfiles + 2 * ov3nfiles, ov3nfiles);
            } else if (nocc == 0) {
                val = (double)valptr[Buf->idx];
                abcd1_terms_new(val, pq, rs, p, q, r, s, o, v, nabcd1, totalnabcd1, abcd1, filesize, bucketsize,
                                abcd1_addr, nfiles);
                abcd2_terms_new(val, pq, rs, p, q, r, s, o, v, nabcd2, totalnabcd2, abcd2, filesize, bucketsize,
                                abcd2_addr, nfiles);
            }
        }
    }
    outfile->Printf("done.\n\n");
    /**
      * write any leftover bits that might not have been dumped to disk
      */
    for (size_t k = 0; k < nfiles; k++) {
        if (nabcd2[k] > 0) {
            psio->open(PSIF_DCC_SORT_START + k + nfiles, PSIO_OPEN_OLD);
            psio->write(PSIF_DCC_SORT_START + k + nfiles, "E2abcd2", (char *)&abcd2[k][0],
                        nabcd2[k] * sizeof(struct integral), abcd2_addr[k], &abcd2_addr[k]);
            psio->close(PSIF_DCC_SORT_START + k + nfiles, 1);
            totalnabcd2[k] += nabcd2[k];
            nabcd2[k] = 0;
        }
        if (nabcd1[k] > 0) {
            psio->open(PSIF_DCC_SORT_START + k, PSIO_OPEN_OLD);
            psio->write(PSIF_DCC_SORT_START + k, "E2abcd1", (char *)&abcd1[k][0], nabcd1[k] * sizeof(struct integral),
                        abcd1_addr[k], &abcd1_addr[k]);
            psio->close(PSIF_DCC_SORT_START + k, 1);
            totalnabcd1[k] += nabcd1[k];
            nabcd1[k] = 0;
        }
    }
    for (size_t k = 0; k < ov3nfiles; k++) {
        if (nabci5[k] > 0) {
            psio->open(PSIF_DCC_SORT_START + k + 2 * nfiles + 2 * ov3nfiles, PSIO_OPEN_OLD);
            psio->write(PSIF_DCC_SORT_START + k + 2 * nfiles + 2 * ov3nfiles, "E2abci2", (char *)&abci5[k][0],
                        nabci5[k] * sizeof(struct integral), abci5_addr[k], &abci5_addr[k]);
            psio->close(PSIF_DCC_SORT_START + k + 2 * nfiles + 2 * ov3nfiles, 1);
            totalnabci5[k] += nabci5[k];
            nabci5[k] = 0;
        }
        if (nabci3[k] > 0) {
            psio->open(PSIF_DCC_SORT_START + k + 2 * nfiles + ov3nfiles, PSIO_OPEN_OLD);
            psio->write(PSIF_DCC_SORT_START + k + 2 * nfiles + ov3nfiles, "E2abci3", (char *)&abci3[k][0],
                        nabci3[k] * sizeof(struct integral), abci3_addr[k], &abci3_addr[k]);
            psio->close(PSIF_DCC_SORT_START + k + 2 * nfiles + ov3nfiles, 1);
            totalnabci3[k] += nabci3[k];
            nabci3[k] = 0;
        }
        if (nabci1[k] > 0) {
            psio->open(PSIF_DCC_SORT_START + k + 2 * nfiles, PSIO_OPEN_OLD);
            psio->write(PSIF_DCC_SORT_START + k + 2 * nfiles, "E2abci", (char *)&abci1[k][0],
                        nabci1[k] * sizeof(struct integral), abci1_addr[k], &abci1_addr[k]);
            psio->close(PSIF_DCC_SORT_START + k + 2 * nfiles, 1);
            totalnabci1[k] += nabci1[k];
            nabci1[k] = 0;
        }
    }
    if (nakjc != 0) {
        psio->open(PSIF_DCC_IJAB, PSIO_OPEN_OLD);
        psio->write(PSIF_DCC_IJAB, "E2ijab", (char *)&akjc[0], nakjc * sizeof(struct integral), akjc_addr, &akjc_addr);
        psio->close(PSIF_DCC_IJAB, 1);
        totalnakjc += nakjc;
        nakjc = 0;
    }
    if (nklcd != 0) {
        psio->open(PSIF_DCC_IAJB, PSIO_OPEN_OLD);
        psio->write(PSIF_DCC_IAJB, "E2iajb", (char *)&klcd[0], nklcd * sizeof(struct integral), klcd_addr, &klcd_addr);
        psio->close(PSIF_DCC_IAJB, 1);
        totalnklcd += nklcd;
        nklcd = 0;
    }
    if (nijkl != 0) {
        psio->open(PSIF_DCC_IJKL, PSIO_OPEN_OLD);
        psio->write(PSIF_DCC_IJKL, "E2ijkl", (char *)&ijkl[0], nijkl * sizeof(struct integral), ijkl_addr, &ijkl_addr);
        psio->close(PSIF_DCC_IJKL, 1);
        totalnijkl += nijkl;
        nijkl = 0;
    }
    if (nijak != 0) {
        psio->open(PSIF_DCC_IJAK, PSIO_OPEN_OLD);
        psio->write(PSIF_DCC_IJAK, "E2ijak", (char *)&ijak[0], nijak * sizeof(struct integral), ijak_addr, &ijak_addr);
        psio->close(PSIF_DCC_IJAK, 1);
        totalnijak += nijak;
        nijak = 0;
    }
    if (nijak2 != 0) {
        psio->open(PSIF_DCC_IJAK2, PSIO_OPEN_OLD);
        psio->write(PSIF_DCC_IJAK2, "E2ijak2", (char *)&ijak2[0], nijak2 * sizeof(struct integral), ijak2_addr,
                    &ijak2_addr);
        psio->close(PSIF_DCC_IJAK2, 1);
        totalnijak2 += nijak2;
        nijak2 = 0;
    }

    /**
      * sort values in each of the files
      */
    double *tmp;
    tmp = new double[maxelem];

    outfile->Printf("        Sort (IJ|KL)........");
    SortBlock(totalnijkl, o * o * o * o, integralbuffer, tmp, PSIF_DCC_IJKL, "E2ijkl", maxelem);
    outfile->Printf("done.\n");
    outfile->Printf("        Sort (IJ|KA) 1/2....");
    SortBlock(totalnijak, o * o * o * v, integralbuffer, tmp, PSIF_DCC_IJAK, "E2ijak", maxelem);
    outfile->Printf("done.\n");
    outfile->Printf("        Sort (IJ|KA) 2/2....");
    SortBlock(totalnijak2, o * o * o * v, integralbuffer, tmp, PSIF_DCC_IJAK2, "E2ijak2", maxelem);
    outfile->Printf("done.\n");
    outfile->Printf("        Sort (IA|JB)........");
    SortBlock(totalnklcd, o * o * v * v, integralbuffer, tmp, PSIF_DCC_IAJB, "E2iajb", maxelem);
    outfile->Printf("done.\n");
    outfile->Printf("        Sort (IJ|AB)........");
    SortBlock(totalnakjc, o * o * v * v, integralbuffer, tmp, PSIF_DCC_IJAB, "E2ijab", maxelem);
    outfile->Printf("done.\n");

    delete[] integralbuffer;

    auto *integralbuffer2 = new integral[maxelem];

    outfile->Printf("        Sort (IA|BC) 1/3....");
    // SortBlock(totalnabci1,o*v*v*v,integralbuffer,tmp,PSIF_DCC_ABCI,"E2abci",maxelem);
    SortBlockNewNew(totalnabci1, ov3, integralbuffer2, tmp, PSIF_DCC_ABCI, "E2abci", maxelem,
                    PSIF_DCC_SORT_START + 2 * nfiles, ov3nfiles);
    outfile->Printf("done.\n");
    outfile->Printf("        Sort (IA|BC) 2/3....");
    SortBlockNewNew(totalnabci3, ov3, integralbuffer2, tmp, PSIF_DCC_ABCI3, "E2abci3", maxelem,
                    PSIF_DCC_SORT_START + 2 * nfiles + ov3nfiles, ov3nfiles);
    // SortBlock(totalnabci3,o*v*v*v,integralbuffer,tmp,PSIF_DCC_ABCI3,"E2abci3",maxelem);
    outfile->Printf("done.\n");
    outfile->Printf("        Sort (IA|BC) 3/3....");
    SortBlockNewNew(totalnabci5, ov3, integralbuffer2, tmp, PSIF_DCC_ABCI2, "E2abci2", maxelem,
                    PSIF_DCC_SORT_START + 2 * nfiles + 2 * ov3nfiles, ov3nfiles);
    // SortBlock(totalnabci5,o*v*v*v,integralbuffer,tmp,PSIF_DCC_ABCI2,"E2abci2",maxelem);
    outfile->Printf("done.\n");

    outfile->Printf("        Sort (AB|CD) 1/2....");
    SortBlockNewNew(totalnabcd1, v * (v + 1) / 2 * v * (v + 1) / 2, integralbuffer2, tmp, PSIF_DCC_ABCD1, "E2abcd1",
                    maxelem, PSIF_DCC_SORT_START, nfiles);
    outfile->Printf("done.\n");
    outfile->Printf("        Sort (AB|CD) 2/2....");
    SortBlockNewNew(totalnabcd2, v * (v + 1) / 2 * v * (v + 1) / 2, integralbuffer2, tmp, PSIF_DCC_ABCD2, "E2abcd2",
                    maxelem, PSIF_DCC_SORT_START + nfiles, nfiles);
    outfile->Printf("done.\n");
    outfile->Printf("\n");

    delete[] integralbuffer2;

    double *tmp2;
    tmp2 = new double[maxelem];
    /**
      *  Sort ABCI2 integrals (actually, just ABCI2-2*ABCI3)
      */
    size_t nbins = 0;
    size_t binsize, lastbin;
    for (size_t i = 1; i <= o * v * v * v; i++) {
        if (maxelem >= (double)o * v * v * v / i) {
            binsize = o * v * v * v / i;
            if (i * binsize < o * v * v * v) binsize++;
            nbins = i;
            break;
        }
    }
    lastbin = o * v * v * v - (nbins - 1) * binsize;
    psio->open(PSIF_DCC_ABCI3, PSIO_OPEN_OLD);
    psio->open(PSIF_DCC_ABCI2, PSIO_OPEN_OLD);
    psio->open(PSIF_DCC_ABCI5, PSIO_OPEN_NEW);

    abci2_addr = PSIO_ZERO;
    abci3_addr[0] = PSIO_ZERO;
    abci5_addr[0] = PSIO_ZERO;
    psio_address abci4_addr = PSIO_ZERO;
    psio_address abci5a_addr = PSIO_ZERO;

    for (size_t i = 0; i < nbins - 1; i++) {
        psio->read(PSIF_DCC_ABCI3, "E2abci3", (char *)&tmp[0], binsize * sizeof(double), abci3_addr[0], &abci3_addr[0]);
        psio->read(PSIF_DCC_ABCI2, "E2abci2", (char *)&tmp2[0], binsize * sizeof(double), abci5_addr[0],
                   &abci5_addr[0]);
        psio->write(PSIF_DCC_ABCI5, "E2abci5", (char *)&tmp2[0], binsize * sizeof(double), abci5a_addr, &abci5a_addr);
        C_DAXPY(binsize, -2.0, tmp, 1, tmp2, 1);
        psio->write(PSIF_DCC_ABCI2, "E2abci2", (char *)&tmp2[0], binsize * sizeof(double), abci2_addr, &abci2_addr);
    }
    psio->read(PSIF_DCC_ABCI3, "E2abci3", (char *)&tmp[0], lastbin * sizeof(double), abci3_addr[0], &abci3_addr[0]);
    psio->read(PSIF_DCC_ABCI2, "E2abci2", (char *)&tmp2[0], lastbin * sizeof(double), abci5_addr[0], &abci5_addr[0]);
    psio->write(PSIF_DCC_ABCI5, "E2abci5", (char *)&tmp2[0], lastbin * sizeof(double), abci5a_addr, &abci5a_addr);
    C_DAXPY(lastbin, -2.0, tmp, 1, tmp2, 1);
    psio->write(PSIF_DCC_ABCI2, "E2abci2", (char *)&tmp2[0], lastbin * sizeof(double), abci2_addr, &abci2_addr);
    psio->close(PSIF_DCC_ABCI2, 1);
    psio->close(PSIF_DCC_ABCI3, 1);
    psio->close(PSIF_DCC_ABCI5, 1);

    /**
      *  Combine ABCD1 and ABCD2 integrals if SJS packing
      */
    for (size_t i = 1; i <= v * (v + 1) / 2 * v * (v + 1) / 2; i++) {
        if (maxelem >= (double)v * (v + 1) / 2 * v * (v + 1) / 2 / i) {
            binsize = v * (v + 1) / 2 * v * (v + 1) / 2 / i;
            if (i * binsize < v * (v + 1) / 2 * v * (v + 1) / 2) binsize++;
            nbins = i;
            break;
        }
    }
    lastbin = v * (v + 1) / 2 * v * (v + 1) / 2 - (nbins - 1) * binsize;
    psio->open(PSIF_DCC_ABCD1, PSIO_OPEN_OLD);
    psio->open(PSIF_DCC_ABCD2, PSIO_OPEN_OLD);
    psio_address abcd1_again = PSIO_ZERO;
    psio_address abcd1_new = PSIO_ZERO;
    psio_address abcd2_new = PSIO_ZERO;
    abcd1_addr[0] = abcd2_addr[0] = PSIO_ZERO;
    for (size_t i = 0; i < nbins - 1; i++) {
        psio->read(PSIF_DCC_ABCD1, "E2abcd1", (char *)&tmp[0], binsize * sizeof(double), abcd1_addr[0], &abcd1_addr[0]);
        psio->read(PSIF_DCC_ABCD2, "E2abcd2", (char *)&tmp2[0], binsize * sizeof(double), abcd2_addr[0],
                   &abcd2_addr[0]);
        C_DAXPY(binsize, -1.0, tmp2, 1, tmp, 1);
        psio->write(PSIF_DCC_ABCD2, "E2abcd2", (char *)&tmp[0], binsize * sizeof(double), abcd2_new, &abcd2_new);
        psio->read(PSIF_DCC_ABCD1, "E2abcd1", (char *)&tmp[0], binsize * sizeof(double), abcd1_again, &abcd1_again);
        C_DAXPY(binsize, 1.0, tmp2, 1, tmp, 1);
        psio->write(PSIF_DCC_ABCD1, "E2abcd1", (char *)&tmp[0], binsize * sizeof(double), abcd1_new, &abcd1_new);
    }
    psio->read(PSIF_DCC_ABCD1, "E2abcd1", (char *)&tmp[0], lastbin * sizeof(double), abcd1_addr[0], &abcd1_addr[0]);
    psio->read(PSIF_DCC_ABCD2, "E2abcd2", (char *)&tmp2[0], lastbin * sizeof(double), abcd2_addr[0], &abcd2_addr[0]);
    C_DAXPY(lastbin, -1.0, tmp2, 1, tmp, 1);
    psio->write(PSIF_DCC_ABCD2, "E2abcd2", (char *)&tmp[0], lastbin * sizeof(double), abcd2_new, &abcd2_new);
    psio->read(PSIF_DCC_ABCD1, "E2abcd1", (char *)&tmp[0], lastbin * sizeof(double), abcd1_again, &abcd1_again);
    C_DAXPY(lastbin, 1.0, tmp2, 1, tmp, 1);
    psio->write(PSIF_DCC_ABCD1, "E2abcd1", (char *)&tmp[0], lastbin * sizeof(double), abcd1_new, &abcd1_new);
    psio->close(PSIF_DCC_ABCD1, 1);
    psio->close(PSIF_DCC_ABCD2, 1);

    delete[] tmp;
    delete[] tmp2;
}
void klcd_terms_incore(double val, size_t pq, size_t rs, size_t p, size_t q, size_t r, size_t s, size_t o, size_t v,
                       double *klcd) {
    size_t k, l, c, d;
    long int ind;
    if (p < o) {
        k = p;
        c = q - o;
        if (r < o) {
            l = r;
            d = s - o;
        } else {
            d = r - o;
            l = s;
        }
    } else {
        c = p - o;
        k = q;
        if (r < o) {
            l = r;
            d = s - o;
        } else {
            d = r - o;
            l = s;
        }
    }
    ind = k * o * v * v + c * o * v + l * v + d;
    klcd[ind] = val;
    ind = l * o * v * v + d * o * v + k * v + c;
    klcd[ind] = val;
}
void klcd_terms(double val, size_t pq, size_t rs, size_t p, size_t q, size_t r, size_t s, size_t o, size_t v,
                size_t &nklcd, struct integral *klcd) {
    size_t k, l, c, d;

    if (p < o) {
        k = p;
        c = q - o;
        if (r < o) {
            l = r;
            d = s - o;
        } else {
            d = r - o;
            l = s;
        }
    } else {
        c = p - o;
        k = q;
        if (r < o) {
            l = r;
            d = s - o;
        } else {
            d = r - o;
            l = s;
        }
    }
    klcd[nklcd].ind = k * o * v * v + c * o * v + l * v + d;
    klcd[nklcd++].val = val;
    if (pq != rs) {
        klcd[nklcd].ind = l * o * v * v + d * o * v + k * v + c;
        klcd[nklcd++].val = val;
    }
}
void ijkl_terms(double val, size_t pq, size_t rs, size_t p, size_t q, size_t r, size_t s, size_t o, size_t &nijkl,
                struct integral *ijkl) {
    if (p == q) {
        if (r == s) {
            ijkl[nijkl].ind = p * o * o * o + r * o * o + q * o + s;
            ijkl[nijkl++].val = val;
            if (pq != rs) {
                ijkl[nijkl].ind = r * o * o * o + p * o * o + s * o + q;
                ijkl[nijkl++].val = val;
            }
        } else {
            ijkl[nijkl].ind = p * o * o * o + r * o * o + q * o + s;
            ijkl[nijkl++].val = val;
            ijkl[nijkl].ind = p * o * o * o + s * o * o + q * o + r;
            ijkl[nijkl++].val = val;
            if (pq != rs) {
                ijkl[nijkl].ind = r * o * o * o + p * o * o + s * o + q;
                ijkl[nijkl++].val = val;
                ijkl[nijkl].ind = s * o * o * o + p * o * o + r * o + q;
                ijkl[nijkl++].val = val;
            }
        }
    } else {
        if (r == s) {
            ijkl[nijkl].ind = p * o * o * o + r * o * o + q * o + s;
            ijkl[nijkl++].val = val;
            ijkl[nijkl].ind = q * o * o * o + r * o * o + p * o + s;
            ijkl[nijkl++].val = val;
            if (pq != rs) {
                ijkl[nijkl].ind = r * o * o * o + p * o * o + s * o + q;
                ijkl[nijkl++].val = val;
                ijkl[nijkl].ind = r * o * o * o + q * o * o + s * o + p;
                ijkl[nijkl++].val = val;
            }
        } else {
            ijkl[nijkl].ind = p * o * o * o + r * o * o + q * o + s;
            ijkl[nijkl++].val = val;
            ijkl[nijkl].ind = q * o * o * o + r * o * o + p * o + s;
            ijkl[nijkl++].val = val;
            ijkl[nijkl].ind = p * o * o * o + s * o * o + q * o + r;
            ijkl[nijkl++].val = val;
            ijkl[nijkl].ind = q * o * o * o + s * o * o + p * o + r;
            ijkl[nijkl++].val = val;
            if (pq != rs) {
                ijkl[nijkl].ind = r * o * o * o + p * o * o + s * o + q;
                ijkl[nijkl++].val = val;
                ijkl[nijkl].ind = r * o * o * o + q * o * o + s * o + p;
                ijkl[nijkl++].val = val;
                ijkl[nijkl].ind = s * o * o * o + p * o * o + r * o + q;
                ijkl[nijkl++].val = val;
                ijkl[nijkl].ind = s * o * o * o + q * o * o + r * o + p;
                ijkl[nijkl++].val = val;
            }
        }
    }
}

void akjc_terms(double val, size_t p, size_t q, size_t r, size_t s, size_t o, size_t v, size_t &nakjc,
                struct integral *akjc) {
    size_t a, k, j, c;
    if (p >= o) {
        a = p - o;
        c = q - o;
        j = r;
        k = s;
    } else {
        j = p;
        k = q;
        a = r - o;
        c = s - o;
    }
    if (j == k) {
        if (a == c) {
            akjc[nakjc].ind = k * o * v * v + c * o * v + j * v + a;
            akjc[nakjc++].val = val;
        } else {
            akjc[nakjc].ind = k * o * v * v + c * o * v + j * v + a;
            akjc[nakjc++].val = val;
            akjc[nakjc].ind = k * o * v * v + a * o * v + j * v + c;
            akjc[nakjc++].val = val;
        }
    } else {
        if (a == c) {
            akjc[nakjc].ind = k * o * v * v + c * o * v + j * v + a;
            akjc[nakjc++].val = val;
            akjc[nakjc].ind = j * o * v * v + c * o * v + k * v + a;
            akjc[nakjc++].val = val;
        } else {
            akjc[nakjc].ind = k * o * v * v + c * o * v + j * v + a;
            akjc[nakjc++].val = val;
            akjc[nakjc].ind = j * o * v * v + c * o * v + k * v + a;
            akjc[nakjc++].val = val;
            akjc[nakjc].ind = k * o * v * v + a * o * v + j * v + c;
            akjc[nakjc++].val = val;
            akjc[nakjc].ind = j * o * v * v + a * o * v + k * v + c;
            akjc[nakjc++].val = val;
        }
    }
}
void ijak2_terms(double val, size_t p, size_t q, size_t r, size_t s, size_t o, size_t v, size_t &nijak2,
                 struct integral *ijak2) {
    size_t i = 0, j = 0, a = 0, k = 0;
    if (p >= o) {
        a = p - o;
        i = q;
        j = r;
        k = s;
    } else if (q >= o) {
        i = p;
        a = q - o;
        j = r;
        k = s;
    } else if (r >= o) {
        a = r - o;
        i = s;
        j = p;
        k = q;
    } else if (s >= o) {
        i = r;
        a = s - o;
        j = p;
        k = q;
    }
    ijak2[nijak2].ind = j * o * o * v + a * o * o + k * o + i;
    ijak2[nijak2++].val = val;
    if (k != j) {
        ijak2[nijak2].ind = k * o * o * v + a * o * o + j * o + i;
        ijak2[nijak2++].val = val;
    }
}
void ijak_terms(double val, size_t p, size_t q, size_t r, size_t s, size_t o, size_t v, size_t &nijak,
                struct integral *ijak) {
    size_t i = 0, j = 0, a = 0, k = 0;
    if (p >= o) {
        a = p - o;
        i = q;
        j = r;
        k = s;
    } else if (q >= o) {
        i = p;
        a = q - o;
        j = r;
        k = s;
    } else if (r >= o) {
        a = r - o;
        i = s;
        j = p;
        k = q;
    } else if (s >= o) {
        i = r;
        a = s - o;
        j = p;
        k = q;
    }
    ijak[nijak].ind = j * o * o * v + i * o * v + k * v + a;
    ijak[nijak++].val = val;
    if (k != j) {
        ijak[nijak].ind = k * o * o * v + i * o * v + j * v + a;
        ijak[nijak++].val = val;
    }
}
void abci5_terms_new(double val, size_t p, size_t q, size_t r, size_t s, size_t o, size_t v, size_t *nabci5,
                     size_t *totalnabci5, struct integral **abci5, size_t binsize, size_t bucketsize,
                     psio_address *addr, size_t filestart, size_t nfiles) {
    size_t k, rem;
    size_t i = 0, a = 0, b = 0, c = 0;
    if (p < o) {
        i = p;
        b = q - o;
        a = r - o;
        c = s - o;
    } else if (q < o) {
        i = q;
        b = p - o;
        a = r - o;
        c = s - o;
    } else if (r < o) {
        i = r;
        b = s - o;
        a = p - o;
        c = q - o;
    } else if (s < o) {
        i = s;
        b = r - o;
        a = p - o;
        c = q - o;
    }
    size_t ind = a * v * v * o + b * v * o + i * v + c;
    rem = ind % binsize;
    k = (ind - rem) / binsize;
    abci5[k][nabci5[k]].ind = ind;
    abci5[k][nabci5[k]++].val = val;
    if (a != c) {
        ind = c * v * v * o + b * v * o + i * v + a;
        rem = ind % binsize;
        k = (ind - rem) / binsize;
        abci5[k][nabci5[k]].ind = ind;
        abci5[k][nabci5[k]++].val = val;
    }
    for (k = 0; k < nfiles; k++) {
        if (nabci5[k] >= bucketsize) {
            auto psio = std::make_shared<PSIO>();
            // write
            psio->open(filestart + k, PSIO_OPEN_OLD);
            psio->write(filestart + k, "E2abci2", (char *)&abci5[k][0], nabci5[k] * sizeof(struct integral), addr[k],
                        &addr[k]);
            psio->close(filestart + k, 1);
            totalnabci5[k] += nabci5[k];
            nabci5[k] = 0;
        }
    }
}
void abci5_terms(double val, size_t p, size_t q, size_t r, size_t s, size_t o, size_t v, size_t &nabci5,
                 struct integral *abci5) {
    size_t i = 0, a = 0, b = 0, c = 0;
    if (p < o) {
        i = p;
        b = q - o;
        a = r - o;
        c = s - o;
    } else if (q < o) {
        i = q;
        b = p - o;
        a = r - o;
        c = s - o;
    } else if (r < o) {
        i = r;
        b = s - o;
        a = p - o;
        c = q - o;
    } else if (s < o) {
        i = s;
        b = r - o;
        a = p - o;
        c = q - o;
    }
    abci5[nabci5].ind = a * v * v * o + b * v * o + i * v + c;
    abci5[nabci5++].val = val;
    if (a != c) {
        abci5[nabci5].ind = c * v * v * o + b * v * o + i * v + a;
        abci5[nabci5++].val = val;
    }
}
void abci3_terms_new(double val, size_t p, size_t q, size_t r, size_t s, size_t o, size_t v, size_t *nabci3,
                     size_t *totalnabci3, struct integral **abci3, size_t binsize, size_t bucketsize,
                     psio_address *addr, size_t filestart, size_t nfiles) {
    size_t k, rem;
    size_t a = 0, f = 0, m = 0, e = 0;
    if (p < o) {
        m = p;
        e = q - o;
        a = r - o;
        f = s - o;
    } else if (q < o) {
        m = q;
        e = p - o;
        a = r - o;
        f = s - o;
    } else if (r < o) {
        m = r;
        e = s - o;
        a = p - o;
        f = q - o;
    } else if (s < o) {
        m = s;
        e = r - o;
        a = p - o;
        f = q - o;
    }
    size_t ind = a * v * v * o + f * v * o + m * v + e;
    rem = ind % binsize;
    k = (ind - rem) / binsize;
    abci3[k][nabci3[k]].ind = ind;
    abci3[k][nabci3[k]++].val = val;
    if (a != f) {
        ind = f * v * v * o + a * v * o + m * v + e;
        rem = ind % binsize;
        k = (ind - rem) / binsize;
        abci3[k][nabci3[k]].ind = ind;
        abci3[k][nabci3[k]++].val = val;
    }
    for (k = 0; k < nfiles; k++) {
        if (nabci3[k] >= bucketsize) {
            auto psio = std::make_shared<PSIO>();
            // write
            psio->open(filestart + k, PSIO_OPEN_OLD);
            psio->write(filestart + k, "E2abci3", (char *)&abci3[k][0], nabci3[k] * sizeof(struct integral), addr[k],
                        &addr[k]);
            psio->close(filestart + k, 1);
            totalnabci3[k] += nabci3[k];
            nabci3[k] = 0;
        }
    }
}
void abci3_terms(double val, size_t p, size_t q, size_t r, size_t s, size_t o, size_t v, size_t &nabci3,
                 struct integral *abci3) {
    size_t a, f, m, e;
    if (p < o) {
        m = p;
        e = q - o;
        a = r - o;
        f = s - o;
    } else if (q < o) {
        m = q;
        e = p - o;
        a = r - o;
        f = s - o;
    } else if (r < o) {
        m = r;
        e = s - o;
        a = p - o;
        f = q - o;
    } else if (s < o) {
        m = s;
        e = r - o;
        a = p - o;
        f = q - o;
    }
    abci3[nabci3].ind = a * v * v * o + f * v * o + m * v + e;
    abci3[nabci3++].val = val;
    if (a != f) {
        abci3[nabci3].ind = f * v * v * o + a * v * o + m * v + e;
        abci3[nabci3++].val = val;
    }
}
void abci1_terms_new(double val, size_t p, size_t q, size_t r, size_t s, size_t o, size_t v, size_t *nabci1,
                     size_t *totalnabci1, struct integral **abci1, size_t binsize, size_t bucketsize,
                     psio_address *addr, size_t filestart, size_t nfiles) {
    size_t k, rem;
    size_t i, a, b, c;
    if (p < o) {
        i = p;
        b = q - o;
        a = r - o;
        c = s - o;
    } else if (q < o) {
        i = q;
        b = p - o;
        a = r - o;
        c = s - o;
    } else if (r < o) {
        i = r;
        b = s - o;
        a = p - o;
        c = q - o;
    } else if (s < o) {
        i = s;
        b = r - o;
        a = p - o;
        c = q - o;
    }
    size_t ind = i * v * v * v + a * v * v + b * v + c;
    rem = ind % binsize;
    k = (ind - rem) / binsize;
    abci1[k][nabci1[k]].ind = ind;
    abci1[k][nabci1[k]++].val = val;
    if (a != c) {
        size_t ind = i * v * v * v + c * v * v + b * v + a;
        rem = ind % binsize;
        k = (ind - rem) / binsize;
        abci1[k][nabci1[k]].ind = ind;
        abci1[k][nabci1[k]++].val = val;
    }
    for (k = 0; k < nfiles; k++) {
        if (nabci1[k] >= bucketsize) {
            auto psio = std::make_shared<PSIO>();
            // write
            psio->open(filestart + k, PSIO_OPEN_OLD);
            psio->write(filestart + k, "E2abci", (char *)&abci1[k][0], nabci1[k] * sizeof(struct integral), addr[k],
                        &addr[k]);
            psio->close(filestart + k, 1);
            totalnabci1[k] += nabci1[k];
            nabci1[k] = 0;
        }
    }
}
void abci1_terms(double val, size_t p, size_t q, size_t r, size_t s, size_t o, size_t v, size_t &nabci1,
                 struct integral *abci1) {
    size_t i, a, b, c;
    if (p < o) {
        i = p;
        b = q - o;
        a = r - o;
        c = s - o;
    } else if (q < o) {
        i = q;
        b = p - o;
        a = r - o;
        c = s - o;
    } else if (r < o) {
        i = r;
        b = s - o;
        a = p - o;
        c = q - o;
    } else if (s < o) {
        i = s;
        b = r - o;
        a = p - o;
        c = q - o;
    }
    abci1[nabci1].ind = i * v * v * v + a * v * v + b * v + c;
    abci1[nabci1++].val = val;
    if (a != c) {
        abci1[nabci1].ind = i * v * v * v + c * v * v + b * v + a;
        abci1[nabci1++].val = val;
    }
}
/**
  * ABCD-type integrals, because of weird SJS packing, are really
  * confusing to sort.  I couldn't think of an analytic way to do
  * this, so I resorted to brute force.
  */
void abcd2_terms(double val, size_t pq, size_t rs, size_t p, size_t q, size_t r, size_t s, size_t o, size_t v,
                 size_t &nabcd2, struct integral *abcd2) {
    size_t ind3, a, b, c, d, ind1, ind2, index, flag;
    size_t nvals, vals[16];
    nvals = 0;

    a = p - o;
    d = q - o;
    b = r - o;
    c = s - o;
    ind1 = Position(a, b);
    ind2 = Position(c, d);
    if ((a <= b && c <= d) || (b <= a && d <= c)) {
        ind3 = ind1 * v * (v + 1) / 2 + ind2;
        flag = 1;
        for (index = 0; index < nvals; index++) {
            if (vals[index] == ind3) {
                flag = 0;
                break;
            }
        }
        if (flag) {
            abcd2[nabcd2].ind = ind3;
            abcd2[nabcd2++].val = val;
            vals[nvals++] = ind3;
        }
        if (ind1 != ind2) {
            ind3 = ind2 * v * (v + 1) / 2 + ind1;
            flag = 1;
            for (index = 0; index < nvals; index++) {
                if (vals[index] == ind3) {
                    flag = 0;
                    break;
                }
            }
            if (flag) {
                abcd2[nabcd2].ind = ind3;
                abcd2[nabcd2++].val = val;
                vals[nvals++] = ind3;
            }
        }
    }
    d = p - o;
    a = q - o;
    b = r - o;
    c = s - o;
    ind1 = Position(a, b);
    ind2 = Position(c, d);
    if ((a <= b && c <= d) || (b <= a && d <= c)) {
        ind3 = ind1 * v * (v + 1) / 2 + ind2;
        flag = 1;
        for (index = 0; index < nvals; index++) {
            if (vals[index] == ind3) {
                flag = 0;
                break;
            }
        }
        if (flag) {
            abcd2[nabcd2].ind = ind3;
            abcd2[nabcd2++].val = val;
            vals[nvals++] = ind3;
        }
        if (ind1 != ind2) {
            ind3 = ind2 * v * (v + 1) / 2 + ind1;
            flag = 1;
            for (index = 0; index < nvals; index++) {
                if (vals[index] == ind3) {
                    flag = 0;
                    break;
                }
            }
            if (flag) {
                abcd2[nabcd2].ind = ind3;
                abcd2[nabcd2++].val = val;
                vals[nvals++] = ind3;
            }
        }
    }
    a = p - o;
    d = q - o;
    c = r - o;
    b = s - o;
    ind1 = Position(a, b);
    ind2 = Position(c, d);
    if ((a <= b && c <= d) || (b <= a && d <= c)) {
        ind3 = ind1 * v * (v + 1) / 2 + ind2;
        flag = 1;
        for (index = 0; index < nvals; index++) {
            if (vals[index] == ind3) {
                flag = 0;
                break;
            }
        }
        if (flag) {
            abcd2[nabcd2].ind = ind3;
            abcd2[nabcd2++].val = val;
            vals[nvals++] = ind3;
        }
        if (ind1 != ind2) {
            ind3 = ind2 * v * (v + 1) / 2 + ind1;
            flag = 1;
            for (index = 0; index < nvals; index++) {
                if (vals[index] == ind3) {
                    flag = 0;
                    break;
                }
            }
            if (flag) {
                abcd2[nabcd2].ind = ind3;
                abcd2[nabcd2++].val = val;
                vals[nvals++] = ind3;
            }
        }
    }
    d = p - o;
    a = q - o;
    c = r - o;
    b = s - o;
    ind1 = Position(a, b);
    ind2 = Position(c, d);
    if ((a <= b && c <= d) || (b <= a && d <= c)) {
        ind3 = ind1 * v * (v + 1) / 2 + ind2;
        flag = 1;
        for (index = 0; index < nvals; index++) {
            if (vals[index] == ind3) {
                flag = 0;
                break;
            }
        }
        if (flag) {
            abcd2[nabcd2].ind = ind3;
            abcd2[nabcd2++].val = val;
            vals[nvals++] = ind3;
        }
        if (ind1 != ind2) {
            ind3 = ind2 * v * (v + 1) / 2 + ind1;
            flag = 1;
            for (index = 0; index < nvals; index++) {
                if (vals[index] == ind3) {
                    flag = 0;
                    break;
                }
            }
            if (flag) {
                abcd2[nabcd2].ind = ind3;
                abcd2[nabcd2++].val = val;
                vals[nvals++] = ind3;
            }
        }
    }
}
void abcd2_terms_new(double val, size_t pq, size_t rs, size_t p, size_t q, size_t r, size_t s, size_t o, size_t v,
                     size_t *nabcd2, size_t *totalnabcd2, struct integral **abcd2, size_t binsize, size_t bucketsize,
                     psio_address *addr, size_t nfiles) {
    size_t ind3, a, b, c, d, ind1, ind2, index, flag;
    size_t nvals, vals[16];
    size_t k;
    nvals = 0;

    a = p - o;
    d = q - o;
    b = r - o;
    c = s - o;
    ind1 = Position(a, b);
    ind2 = Position(c, d);
    if ((a <= b && c <= d) || (b <= a && d <= c)) {
        ind3 = ind1 * v * (v + 1) / 2 + ind2;
        flag = 1;
        for (index = 0; index < nvals; index++) {
            if (vals[index] == ind3) {
                flag = 0;
                break;
            }
        }
        if (flag) {
            size_t rem = ind3 % binsize;
            k = (ind3 - rem) / binsize;
            abcd2[k][nabcd2[k]].ind = ind3;
            abcd2[k][nabcd2[k]++].val = val;
            vals[nvals++] = ind3;
        }
        if (ind1 != ind2) {
            ind3 = ind2 * v * (v + 1) / 2 + ind1;
            flag = 1;
            for (index = 0; index < nvals; index++) {
                if (vals[index] == ind3) {
                    flag = 0;
                    break;
                }
            }
            if (flag) {
                size_t rem = ind3 % binsize;
                k = (ind3 - rem) / binsize;
                abcd2[k][nabcd2[k]].ind = ind3;
                abcd2[k][nabcd2[k]++].val = val;
                vals[nvals++] = ind3;
            }
        }
    }
    d = p - o;
    a = q - o;
    b = r - o;
    c = s - o;
    ind1 = Position(a, b);
    ind2 = Position(c, d);
    if ((a <= b && c <= d) || (b <= a && d <= c)) {
        ind3 = ind1 * v * (v + 1) / 2 + ind2;
        flag = 1;
        for (index = 0; index < nvals; index++) {
            if (vals[index] == ind3) {
                flag = 0;
                break;
            }
        }
        if (flag) {
            size_t rem = ind3 % binsize;
            k = (ind3 - rem) / binsize;
            abcd2[k][nabcd2[k]].ind = ind3;
            abcd2[k][nabcd2[k]++].val = val;
            vals[nvals++] = ind3;
        }
        if (ind1 != ind2) {
            ind3 = ind2 * v * (v + 1) / 2 + ind1;
            flag = 1;
            for (index = 0; index < nvals; index++) {
                if (vals[index] == ind3) {
                    flag = 0;
                    break;
                }
            }
            if (flag) {
                size_t rem = ind3 % binsize;
                k = (ind3 - rem) / binsize;
                abcd2[k][nabcd2[k]].ind = ind3;
                abcd2[k][nabcd2[k]++].val = val;
                vals[nvals++] = ind3;
            }
        }
    }
    a = p - o;
    d = q - o;
    c = r - o;
    b = s - o;
    ind1 = Position(a, b);
    ind2 = Position(c, d);
    if ((a <= b && c <= d) || (b <= a && d <= c)) {
        ind3 = ind1 * v * (v + 1) / 2 + ind2;
        flag = 1;
        for (index = 0; index < nvals; index++) {
            if (vals[index] == ind3) {
                flag = 0;
                break;
            }
        }
        if (flag) {
            size_t rem = ind3 % binsize;
            k = (ind3 - rem) / binsize;
            abcd2[k][nabcd2[k]].ind = ind3;
            abcd2[k][nabcd2[k]++].val = val;
            vals[nvals++] = ind3;
        }
        if (ind1 != ind2) {
            ind3 = ind2 * v * (v + 1) / 2 + ind1;
            flag = 1;
            for (index = 0; index < nvals; index++) {
                if (vals[index] == ind3) {
                    flag = 0;
                    break;
                }
            }
            if (flag) {
                size_t rem = ind3 % binsize;
                k = (ind3 - rem) / binsize;
                abcd2[k][nabcd2[k]].ind = ind3;
                abcd2[k][nabcd2[k]++].val = val;
                vals[nvals++] = ind3;
            }
        }
    }
    d = p - o;
    a = q - o;
    c = r - o;
    b = s - o;
    ind1 = Position(a, b);
    ind2 = Position(c, d);
    if ((a <= b && c <= d) || (b <= a && d <= c)) {
        ind3 = ind1 * v * (v + 1) / 2 + ind2;
        flag = 1;
        for (index = 0; index < nvals; index++) {
            if (vals[index] == ind3) {
                flag = 0;
                break;
            }
        }
        if (flag) {
            size_t rem = ind3 % binsize;
            k = (ind3 - rem) / binsize;
            abcd2[k][nabcd2[k]].ind = ind3;
            abcd2[k][nabcd2[k]++].val = val;
            vals[nvals++] = ind3;
        }
        if (ind1 != ind2) {
            ind3 = ind2 * v * (v + 1) / 2 + ind1;
            flag = 1;
            for (index = 0; index < nvals; index++) {
                if (vals[index] == ind3) {
                    flag = 0;
                    break;
                }
            }
            if (flag) {
                size_t rem = ind3 % binsize;
                k = (ind3 - rem) / binsize;
                abcd2[k][nabcd2[k]].ind = ind3;
                abcd2[k][nabcd2[k]++].val = val;
                vals[nvals++] = ind3;
            }
        }
    }
    for (k = 0; k < nfiles; k++) {
        if (nabcd2[k] >= bucketsize) {
            auto psio = std::make_shared<PSIO>();
            // write
            psio->open(PSIF_DCC_SORT_START + k + nfiles, PSIO_OPEN_OLD);
            psio->write(PSIF_DCC_SORT_START + k + nfiles, "E2abcd2", (char *)&abcd2[k][0],
                        nabcd2[k] * sizeof(struct integral), addr[k], &addr[k]);
            psio->close(PSIF_DCC_SORT_START + k + nfiles, 1);
            totalnabcd2[k] += nabcd2[k];
            nabcd2[k] = 0;
            // memset((void*)abcd2[k],'\0',bucketsize*sizeof(struct integral));
        }
    }
}
void abcd3_terms(double val, size_t pq, size_t rs, size_t p, size_t q, size_t r, size_t s, size_t o, size_t v,
                 size_t &nabcd3, struct integral *abcd3) {
    p -= o;
    q -= o;
    r -= o;
    s -= o;
    if (p == q) {
        if (r == s) {
            abcd3[nabcd3].ind = p * v * v * v + r * v * v + q * v + s;
            abcd3[nabcd3++].val = val;
            if (pq != rs) {
                abcd3[nabcd3].ind = r * v * v * v + p * v * v + s * v + q;
                abcd3[nabcd3++].val = val;
            }
        } else {
            abcd3[nabcd3].ind = p * v * v * v + r * v * v + q * v + s;
            abcd3[nabcd3++].val = val;
            abcd3[nabcd3].ind = p * v * v * v + s * v * v + q * v + r;
            abcd3[nabcd3++].val = val;
            if (pq != rs) {
                abcd3[nabcd3].ind = r * v * v * v + p * v * v + s * v + q;
                abcd3[nabcd3++].val = val;
                abcd3[nabcd3].ind = s * v * v * v + p * v * v + r * v + q;
                abcd3[nabcd3++].val = val;
            }
        }
    } else {
        if (r == s) {
            abcd3[nabcd3].ind = p * v * v * v + r * v * v + q * v + s;
            abcd3[nabcd3++].val = val;
            abcd3[nabcd3].ind = q * v * v * v + r * v * v + p * v + s;
            abcd3[nabcd3++].val = val;
            if (pq != rs) {
                abcd3[nabcd3].ind = r * v * v * v + p * v * v + s * v + q;
                abcd3[nabcd3++].val = val;
                abcd3[nabcd3].ind = r * v * v * v + q * v * v + s * v + p;
                abcd3[nabcd3++].val = val;
            }
        } else {
            abcd3[nabcd3].ind = p * v * v * v + r * v * v + q * v + s;
            abcd3[nabcd3++].val = val;
            abcd3[nabcd3].ind = q * v * v * v + r * v * v + p * v + s;
            abcd3[nabcd3++].val = val;
            abcd3[nabcd3].ind = p * v * v * v + s * v * v + q * v + r;
            abcd3[nabcd3++].val = val;
            abcd3[nabcd3].ind = q * v * v * v + s * v * v + p * v + r;
            abcd3[nabcd3++].val = val;
            if (pq != rs) {
                abcd3[nabcd3].ind = r * v * v * v + p * v * v + s * v + q;
                abcd3[nabcd3++].val = val;
                abcd3[nabcd3].ind = r * v * v * v + q * v * v + s * v + p;
                abcd3[nabcd3++].val = val;
                abcd3[nabcd3].ind = s * v * v * v + p * v * v + r * v + q;
                abcd3[nabcd3++].val = val;
                abcd3[nabcd3].ind = s * v * v * v + q * v * v + r * v + p;
                abcd3[nabcd3++].val = val;
            }
        }
    }
}
void abcd1_terms_new(double val, size_t pq, size_t rs, size_t p, size_t q, size_t r, size_t s, size_t o, size_t v,
                     size_t *nabcd1, size_t *totalnabcd1, struct integral **abcd1, size_t binsize, size_t bucketsize,
                     psio_address *addr, size_t nfiles) {
    size_t ind3, a, b, c, d, ind1, ind2, index, flag;
    size_t nvals, vals[16];
    size_t k;
    nvals = 0;

    a = p - o;
    c = q - o;
    b = r - o;
    d = s - o;
    ind1 = Position(a, b);
    ind2 = Position(c, d);
    if ((a <= b && c <= d) || (b <= a && d <= c)) {
        ind3 = ind1 * v * (v + 1) / 2 + ind2;
        flag = 1;
        for (index = 0; index < nvals; index++) {
            if (vals[index] == ind3) {
                flag = 0;
                break;
            }
        }
        if (flag) {
            size_t rem = ind3 % binsize;
            k = (ind3 - rem) / binsize;
            abcd1[k][nabcd1[k]].ind = ind3;
            abcd1[k][nabcd1[k]++].val = val;
            vals[nvals++] = ind3;
        }
        if (ind1 != ind2) {
            ind3 = ind2 * v * (v + 1) / 2 + ind1;
            flag = 1;
            for (index = 0; index < nvals; index++) {
                if (vals[index] == ind3) {
                    flag = 0;
                    break;
                }
            }
            if (flag) {
                size_t rem = ind3 % binsize;
                k = (ind3 - rem) / binsize;
                abcd1[k][nabcd1[k]].ind = ind3;
                abcd1[k][nabcd1[k]++].val = val;
                vals[nvals++] = ind3;
            }
        }
    }
    c = p - o;
    a = q - o;
    b = r - o;
    d = s - o;
    ind1 = Position(a, b);
    ind2 = Position(c, d);
    if ((a <= b && c <= d) || (b <= a && d <= c)) {
        ind3 = ind1 * v * (v + 1) / 2 + ind2;
        flag = 1;
        for (index = 0; index < nvals; index++) {
            if (vals[index] == ind3) {
                flag = 0;
                break;
            }
        }
        if (flag) {
            size_t rem = ind3 % binsize;
            k = (ind3 - rem) / binsize;
            abcd1[k][nabcd1[k]].ind = ind3;
            abcd1[k][nabcd1[k]++].val = val;
            vals[nvals++] = ind3;
        }
        if (ind1 != ind2) {
            ind3 = ind2 * v * (v + 1) / 2 + ind1;
            flag = 1;
            for (index = 0; index < nvals; index++) {
                if (vals[index] == ind3) {
                    flag = 0;
                    break;
                }
            }
            if (flag) {
                size_t rem = ind3 % binsize;
                k = (ind3 - rem) / binsize;
                abcd1[k][nabcd1[k]].ind = ind3;
                abcd1[k][nabcd1[k]++].val = val;
                vals[nvals++] = ind3;
            }
        }
    }
    a = p - o;
    c = q - o;
    d = r - o;
    b = s - o;
    ind1 = Position(a, b);
    ind2 = Position(c, d);
    if ((a <= b && c <= d) || (b <= a && d <= c)) {
        ind3 = ind1 * v * (v + 1) / 2 + ind2;
        flag = 1;
        for (index = 0; index < nvals; index++) {
            if (vals[index] == ind3) {
                flag = 0;
                break;
            }
        }
        if (flag) {
            size_t rem = ind3 % binsize;
            k = (ind3 - rem) / binsize;
            abcd1[k][nabcd1[k]].ind = ind3;
            abcd1[k][nabcd1[k]++].val = val;
            vals[nvals++] = ind3;
        }
        if (ind1 != ind2) {
            ind3 = ind2 * v * (v + 1) / 2 + ind1;
            flag = 1;
            for (index = 0; index < nvals; index++) {
                if (vals[index] == ind3) {
                    flag = 0;
                    break;
                }
            }
            if (flag) {
                size_t rem = ind3 % binsize;
                k = (ind3 - rem) / binsize;
                abcd1[k][nabcd1[k]].ind = ind3;
                abcd1[k][nabcd1[k]++].val = val;
                vals[nvals++] = ind3;
            }
        }
    }
    c = p - o;
    a = q - o;
    d = r - o;
    b = s - o;
    ind1 = Position(a, b);
    ind2 = Position(c, d);
    if ((a <= b && c <= d) || (b <= a && d <= c)) {
        ind3 = ind1 * v * (v + 1) / 2 + ind2;
        flag = 1;
        for (index = 0; index < nvals; index++) {
            if (vals[index] == ind3) {
                flag = 0;
                break;
            }
        }
        if (flag) {
            size_t rem = ind3 % binsize;
            k = (ind3 - rem) / binsize;
            abcd1[k][nabcd1[k]].ind = ind3;
            abcd1[k][nabcd1[k]++].val = val;
            vals[nvals++] = ind3;
        }
        if (ind1 != ind2) {
            ind3 = ind2 * v * (v + 1) / 2 + ind1;
            flag = 1;
            for (index = 0; index < nvals; index++) {
                if (vals[index] == ind3) {
                    flag = 0;
                    break;
                }
            }
            if (flag) {
                size_t rem = ind3 % binsize;
                k = (ind3 - rem) / binsize;
                abcd1[k][nabcd1[k]].ind = ind3;
                abcd1[k][nabcd1[k]++].val = val;
                vals[nvals++] = ind3;
            }
        }
    }

    for (k = 0; k < nfiles; k++) {
        if (nabcd1[k] >= bucketsize) {
            // write
            auto psio = std::make_shared<PSIO>();
            psio->open(PSIF_DCC_SORT_START + k, PSIO_OPEN_OLD);
            psio->write(PSIF_DCC_SORT_START + k, "E2abcd1", (char *)&abcd1[k][0], nabcd1[k] * sizeof(struct integral),
                        addr[k], &addr[k]);
            psio->close(PSIF_DCC_SORT_START + k, 1);
            totalnabcd1[k] += nabcd1[k];
            nabcd1[k] = 0;
            // memset((void*)abcd1[k],'\0',bucketsize*sizeof(struct integral));
        }
    }
}

void abcd1_terms(double val, size_t pq, size_t rs, size_t p, size_t q, size_t r, size_t s, size_t o, size_t v,
                 size_t &nabcd1, struct integral *abcd1) {
    size_t ind3, a, b, c, d, ind1, ind2, index, flag;
    size_t nvals, vals[16];
    nvals = 0;

    a = p - o;
    c = q - o;
    b = r - o;
    d = s - o;
    ind1 = Position(a, b);
    ind2 = Position(c, d);
    if ((a <= b && c <= d) || (b <= a && d <= c)) {
        ind3 = ind1 * v * (v + 1) / 2 + ind2;
        flag = 1;
        for (index = 0; index < nvals; index++) {
            if (vals[index] == ind3) {
                flag = 0;
                break;
            }
        }
        if (flag) {
            abcd1[nabcd1].ind = ind3;
            abcd1[nabcd1++].val = val;
            vals[nvals++] = ind3;
        }
        if (ind1 != ind2) {
            ind3 = ind2 * v * (v + 1) / 2 + ind1;
            flag = 1;
            for (index = 0; index < nvals; index++) {
                if (vals[index] == ind3) {
                    flag = 0;
                    break;
                }
            }
            if (flag) {
                abcd1[nabcd1].ind = ind3;
                abcd1[nabcd1++].val = val;
                vals[nvals++] = ind3;
            }
        }
    }
    c = p - o;
    a = q - o;
    b = r - o;
    d = s - o;
    ind1 = Position(a, b);
    ind2 = Position(c, d);
    if ((a <= b && c <= d) || (b <= a && d <= c)) {
        ind3 = ind1 * v * (v + 1) / 2 + ind2;
        flag = 1;
        for (index = 0; index < nvals; index++) {
            if (vals[index] == ind3) {
                flag = 0;
                break;
            }
        }
        if (flag) {
            abcd1[nabcd1].ind = ind3;
            abcd1[nabcd1++].val = val;
            vals[nvals++] = ind3;
        }
        if (ind1 != ind2) {
            ind3 = ind2 * v * (v + 1) / 2 + ind1;
            flag = 1;
            for (index = 0; index < nvals; index++) {
                if (vals[index] == ind3) {
                    flag = 0;
                    break;
                }
            }
            if (flag) {
                abcd1[nabcd1].ind = ind3;
                abcd1[nabcd1++].val = val;
                vals[nvals++] = ind3;
            }
        }
    }
    a = p - o;
    c = q - o;
    d = r - o;
    b = s - o;
    ind1 = Position(a, b);
    ind2 = Position(c, d);
    if ((a <= b && c <= d) || (b <= a && d <= c)) {
        ind3 = ind1 * v * (v + 1) / 2 + ind2;
        flag = 1;
        for (index = 0; index < nvals; index++) {
            if (vals[index] == ind3) {
                flag = 0;
                break;
            }
        }
        if (flag) {
            abcd1[nabcd1].ind = ind3;
            abcd1[nabcd1++].val = val;
            vals[nvals++] = ind3;
        }
        if (ind1 != ind2) {
            ind3 = ind2 * v * (v + 1) / 2 + ind1;
            flag = 1;
            for (index = 0; index < nvals; index++) {
                if (vals[index] == ind3) {
                    flag = 0;
                    break;
                }
            }
            if (flag) {
                abcd1[nabcd1].ind = ind3;
                abcd1[nabcd1++].val = val;
                vals[nvals++] = ind3;
            }
        }
    }
    c = p - o;
    a = q - o;
    d = r - o;
    b = s - o;
    ind1 = Position(a, b);
    ind2 = Position(c, d);
    if ((a <= b && c <= d) || (b <= a && d <= c)) {
        ind3 = ind1 * v * (v + 1) / 2 + ind2;
        flag = 1;
        for (index = 0; index < nvals; index++) {
            if (vals[index] == ind3) {
                flag = 0;
                break;
            }
        }
        if (flag) {
            abcd1[nabcd1].ind = ind3;
            abcd1[nabcd1++].val = val;
            vals[nvals++] = ind3;
        }
        if (ind1 != ind2) {
            ind3 = ind2 * v * (v + 1) / 2 + ind1;
            flag = 1;
            for (index = 0; index < nvals; index++) {
                if (vals[index] == ind3) {
                    flag = 0;
                    break;
                }
            }
            if (flag) {
                abcd1[nabcd1].ind = ind3;
                abcd1[nabcd1++].val = val;
                vals[nvals++] = ind3;
            }
        }
    }
}

void SortBlockNewNew(size_t *nelem, size_t blockdim, struct integral *buffer, double *tmp, size_t PSIFILE,
                     const char *string, size_t maxdim, size_t filestart, size_t nfiles) {
    auto psio = std::make_shared<PSIO>();
    // bins are coreloads
    size_t nbins = 0;
    size_t binsize, lastbin;
    for (size_t i = 1; i <= blockdim; i++) {
        if (maxdim >= (double)blockdim / i) {
            binsize = blockdim / i;
            if (i * binsize < blockdim) binsize++;
            nbins = i;
            break;
        }
    }
    lastbin = blockdim - (nbins - 1) * binsize;
    // open temporary files:
    for (size_t i = 0; i < nbins; i++) {
        psio->open(filestart + i, PSIO_OPEN_OLD);
    }

    // loop over temporary files:
    psio_address addrwrite = PSIO_ZERO;
    psio->open(PSIF_DCC_TEMP, PSIO_OPEN_NEW);
    for (size_t k = 0; k < nbins; k++) {
        memset((void *)tmp, '\0', binsize * sizeof(double));
        // read coreload from this file
        psio->read_entry(filestart + k, string, (char *)&buffer[0], nelem[k] * sizeof(struct integral));
        // loop over elements and sort into tmp
        for (size_t j = 0; j < nelem[k]; j++) {
            tmp[buffer[j].ind - k * binsize] = buffer[j].val;
        }
        // write this guy to PSIF_DCC_TEMP
        psio->write(PSIF_DCC_TEMP, string, (char *)&tmp[0], binsize * sizeof(double), addrwrite, &addrwrite);
    }
    psio->close(PSIF_DCC_TEMP, 1);
    psio->rename_file(PSIF_DCC_TEMP, PSIFILE);

    // close temporary files:
    for (size_t i = 0; i < nbins; i++) {
        psio->close(filestart + i, 0);
    }
}

void SortBlockNew(size_t nelem, size_t blockdim, struct integral *buffer, double *tmp, size_t PSIFILE,
                  const char *string, size_t maxdim) {
    auto psio = std::make_shared<PSIO>();
    // does the block fit in core?
    if (nelem <= maxdim && blockdim <= maxdim) {
        psio->open(PSIFILE, PSIO_OPEN_OLD);
        psio->read_entry(PSIFILE, string, (char *)&buffer[0], nelem * sizeof(struct integral));
        psio->close(PSIFILE, 0);

        memset((void *)tmp, '\0', blockdim * sizeof(double));
        for (size_t j = 0; j < nelem; j++) {
            tmp[buffer[j].ind] = buffer[j].val;
        }

        psio->open(PSIFILE, PSIO_OPEN_NEW);
        psio->write_entry(PSIFILE, string, (char *)&tmp[0], blockdim * sizeof(double));
        psio->close(PSIFILE, 1);
    } else {
        // bins are coreloads
        size_t nbins, binsize, lastbin;
        for (size_t i = 1; i <= blockdim; i++) {
            if (maxdim >= (double)blockdim / i) {
                binsize = blockdim / i;
                if (i * binsize < blockdim) binsize++;
                nbins = i;
                break;
            }
        }
        lastbin = blockdim - (nbins - 1) * binsize;
        outfile->Printf("  %5li files, %5li %5li\n", nbins, blockdim, binsize);

        // new integral buffers to hold rough sorted ints
        struct integral **buffer2;
        buffer2 = (struct integral **)malloc(nbins * sizeof(struct integral *));
        // open temporary files:
        for (size_t i = 0; i < nbins; i++) {
            psio->open(PSIF_DCC_SORT_START + i, PSIO_OPEN_NEW);
        }
        // count elements in buffer2
        size_t *buffer2_count = (size_t *)malloc(nbins * sizeof(long int));
        // total elements from buffer2 written to disk
        size_t *buffer2_total = (size_t *)malloc(nbins * sizeof(long int));
        // array of addresses
        auto *addr = new psio_address[nbins];
        // bucketsize is binsize / nbins
        size_t bucketsize = binsize / nbins;
        if (bucketsize * nbins < binsize) bucketsize++;
        printf("%5li %5li\n", binsize, bucketsize);
        for (size_t i = 0; i < nbins; i++) {
            addr[i] = PSIO_ZERO;
            buffer2_count[i] = 0;
            buffer2_total[i] = 0;
            buffer2[i] = (struct integral *)malloc(bucketsize * sizeof(struct integral));
            memset((void *)buffer2[i], '\0', bucketsize * sizeof(struct integral));
        }

        size_t initialnbins = 0;
        size_t initialbinsize, initiallastbin;
        for (size_t i = 1; i <= nelem; i++) {
            if (maxdim >= (double)nelem / i) {
                initialbinsize = nelem / i;
                if (i * initialbinsize < nelem) initialbinsize++;
                initialnbins = i;
                break;
            }
        }
        initiallastbin = nelem - (initialnbins - 1) * initialbinsize;

        psio_address addr1, addr2;
        addr1 = PSIO_ZERO;

        // loop over data on disk:
        psio->open(PSIFILE, PSIO_OPEN_OLD);
        for (size_t i = 0; i < initialnbins - 1; i++) {
            // read buffer
            psio->read(PSIFILE, string, (char *)&buffer[0], initialbinsize * sizeof(struct integral), addr1, &addr1);
            // loop over elements:
            for (size_t j = 0; j < initialbinsize; j++) {
                // which bin to i belong to?
                size_t myk;
                for (size_t k = 0; k < nbins; k++) {
                    if (buffer[j].ind < (k + 1) * binsize && buffer[j].ind >= k * binsize) {
                        myk = k;
                        break;
                    }
                }
                buffer2[myk][buffer2_count[myk]].ind = buffer[j].ind;
                buffer2[myk][buffer2_count[myk]++].val = buffer[j].val;
                if (buffer2_count[myk] == bucketsize) {
                    psio->write(PSIF_DCC_SORT_START + myk, string, (char *)&buffer2[myk][0],
                                bucketsize * sizeof(struct integral), addr[myk], &addr[myk]);
                    buffer2_total[myk] += bucketsize;
                    buffer2_count[myk] = 0;
                    memset((void *)buffer2[myk], '\0', bucketsize * sizeof(struct integral));
                }
            }
        }
        // lastbin
        // read buffer
        psio->read(PSIFILE, string, (char *)&buffer[0], initiallastbin * sizeof(struct integral), addr1, &addr1);
        // loop over elements:
        for (size_t j = 0; j < initiallastbin; j++) {
            // which bin to i belong to?
            size_t myk;
            for (size_t k = 0; k < nbins; k++) {
                if (buffer[j].ind < (k + 1) * binsize && buffer[j].ind >= k * binsize) {
                    myk = k;
                    break;
                }
            }
            buffer2[myk][buffer2_count[myk]].ind = buffer[j].ind;
            buffer2[myk][buffer2_count[myk]++].val = buffer[j].val;
            if (buffer2_count[myk] == bucketsize) {
                psio->write(PSIF_DCC_SORT_START + myk, string, (char *)&buffer2[myk][0],
                            bucketsize * sizeof(struct integral), addr[myk], &addr[myk]);
                buffer2_total[myk] += bucketsize;
                buffer2_count[myk] = 0;
                memset((void *)buffer2[myk], '\0', bucketsize * sizeof(struct integral));
            }
        }
        psio->close(PSIFILE, 1);

        // write leftovers
        for (size_t k = 0; k < nbins; k++) {
            if (buffer2_count[k] > 0) {
                psio->write(PSIF_DCC_SORT_START + k, string, (char *)&buffer2[k][0],
                            buffer2_count[k] * sizeof(struct integral), addr[k], &addr[k]);
                buffer2_total[k] += buffer2_count[k];
                buffer2_count[k] = 0;
                memset((void *)buffer2[k], '\0', bucketsize * sizeof(struct integral));
            }
        }

        // loop over temporary files:
        psio_address addrwrite = PSIO_ZERO;
        psio->open(PSIF_DCC_TEMP, PSIO_OPEN_NEW);
        for (size_t k = 0; k < nbins; k++) {
            memset((void *)tmp, '\0', binsize * sizeof(double));
            // read coreload from this file
            psio->read_entry(PSIF_DCC_SORT_START + k, string, (char *)&buffer[0],
                             buffer2_total[k] * sizeof(struct integral));
            // loop over elements and sort into tmp
            for (size_t j = 0; j < buffer2_total[k]; j++) {
                tmp[buffer[j].ind - k * binsize] = buffer[j].val;
            }
            // write this guy to PSIF_DCC_TEMP
            psio->write(PSIF_DCC_TEMP, string, (char *)&tmp[0], binsize * sizeof(double), addrwrite, &addrwrite);
        }
        psio->close(PSIF_DCC_TEMP, 1);
        psio->rename_file(PSIF_DCC_TEMP, PSIFILE);

        /*psio->open(PSIF_DCC_TEMP,PSIO_OPEN_NEW);
        addr2=PSIO_ZERO;
        for (size_t k=0; k<nbins; k++){
            addr1=PSIO_ZERO;
            memset((void*)tmp,'\0',binsize*sizeof(double));
            for (size_t i=0; i<initialnbins-1; i++){
                psio->read(PSIFILE,string,(char*)&buffer[0],initialbinsize*sizeof(struct integral),addr1,&addr1);
                for (size_t j=0; j<initialbinsize; j++){
                    if (buffer[j].ind < (k+1)*binsize && buffer[j].ind>=k*binsize){
                        tmp[buffer[j].ind-k*binsize] = buffer[j].val;
                    }
                }
            }
            psio->read(PSIFILE,string,(char*)&buffer[0],initiallastbin*sizeof(struct integral),addr1,&addr1);
            for (size_t j=0; j<initiallastbin; j++){
                if (buffer[j].ind < (k+1)*binsize && buffer[j].ind>=k*binsize){
                    tmp[buffer[j].ind-k*binsize] = buffer[j].val;
                }
            }
            psio->write(PSIF_DCC_TEMP,string,(char*)&tmp[0],binsize*sizeof(double),addr2,&addr2);
        }
        psio->close(PSIFILE,1);
        psio->close(PSIF_DCC_TEMP,1);

        psio->rename_file(PSIF_DCC_TEMP,PSIFILE);*/

        delete[] addr;

        // close temporary files:
        for (size_t i = 0; i < nbins; i++) {
            psio->close(PSIF_DCC_SORT_START + i, 0);
            free(buffer2[i]);
        }
        // free extra memory
        free(buffer2_count);
        free(buffer2_total);
        free(buffer2);
    }
}
void SortBlock(size_t nelem, size_t blockdim, struct integral *buffer, double *tmp, size_t PSIFILE, const char *string,
               size_t maxdim) {
    auto psio = std::make_shared<PSIO>();
    // does the block fit in core?
    if (nelem <= maxdim && blockdim <= maxdim) {
        psio->open(PSIFILE, PSIO_OPEN_OLD);
        psio->read_entry(PSIFILE, string, (char *)&buffer[0], nelem * sizeof(struct integral));
        psio->close(PSIFILE, 0);

        memset((void *)tmp, '\0', blockdim * sizeof(double));
        for (size_t j = 0; j < nelem; j++) {
            tmp[buffer[j].ind] = buffer[j].val;
        }

        psio->open(PSIFILE, PSIO_OPEN_NEW);
        psio->write_entry(PSIFILE, string, (char *)&tmp[0], blockdim * sizeof(double));
        psio->close(PSIFILE, 1);
    } else {
        size_t nbins, binsize, lastbin;
        for (size_t i = 1; i <= blockdim; i++) {
            if (maxdim >= (double)blockdim / i) {
                binsize = blockdim / i;
                if (i * binsize < blockdim) binsize++;
                nbins = i;
                break;
            }
        }
        lastbin = blockdim - (nbins - 1) * binsize;

        size_t initialnbins = 0;
        size_t initialbinsize, initiallastbin;
        for (size_t i = 1; i <= nelem; i++) {
            if (maxdim >= (double)nelem / i) {
                initialbinsize = nelem / i;
                if (i * initialbinsize < nelem) initialbinsize++;
                initialnbins = i;
                break;
            }
        }
        initiallastbin = nelem - (initialnbins - 1) * initialbinsize;

        psio_address *addr, addr1, addr2;
        addr = new psio_address[nbins];
        addr1 = PSIO_ZERO;

        psio->open(PSIFILE, PSIO_OPEN_OLD);

        psio->open(PSIF_DCC_TEMP, PSIO_OPEN_NEW);
        addr2 = PSIO_ZERO;
        for (size_t k = 0; k < nbins; k++) {
            addr1 = PSIO_ZERO;
            memset((void *)tmp, '\0', binsize * sizeof(double));
            for (size_t i = 0; i < initialnbins - 1; i++) {
                psio->read(PSIFILE, string, (char *)&buffer[0], initialbinsize * sizeof(struct integral), addr1,
                           &addr1);
                for (size_t j = 0; j < initialbinsize; j++) {
                    if (buffer[j].ind < (k + 1) * binsize && buffer[j].ind >= k * binsize) {
                        tmp[buffer[j].ind - k * binsize] = buffer[j].val;
                    }
                }
            }
            psio->read(PSIFILE, string, (char *)&buffer[0], initiallastbin * sizeof(struct integral), addr1, &addr1);
            for (size_t j = 0; j < initiallastbin; j++) {
                if (buffer[j].ind < (k + 1) * binsize && buffer[j].ind >= k * binsize) {
                    tmp[buffer[j].ind - k * binsize] = buffer[j].val;
                }
            }
            psio->write(PSIF_DCC_TEMP, string, (char *)&tmp[0], binsize * sizeof(double), addr2, &addr2);
        }
        psio->close(PSIFILE, 1);
        psio->close(PSIF_DCC_TEMP, 1);

        psio->rename_file(PSIF_DCC_TEMP, PSIFILE);

        delete[] addr;
    }
}

/*
 * low memory and cim (t) routines need the ov^3 integrals in a different
 * ordering than the efficient default routine.  this function writes the
 * integrals to disk in the proper ordering.
 */
void Sort_OV3_LowMemory(long int memory, long int o, long int v) {
    outfile->Printf("\n");
    outfile->Printf("\n");
    outfile->Printf("        ==> Resort (ov|vv) integrals for low-memory (T) computation <==\n");
    outfile->Printf("\n");

    long int maxelem = memory / 8 / 2;
    auto *tmp = new double[maxelem];
    auto *tmp2 = new double[maxelem];

    // in the interest of diskspace, get rid of PSIF_DCC_ABCI
    auto psio = std::make_shared<PSIO>();
    psio->open(PSIF_DCC_ABCI, PSIO_OPEN_NEW);
    psio->close(PSIF_DCC_ABCI, 0);

    size_t nbins, binsize, lastbin;
    for (size_t i = 1; i <= o * v * v * v; i++) {
        if (maxelem >= (double)o * v * v * v / i) {
            binsize = o * v * v * v / i;
            if (i * binsize < o * v * v * v) binsize++;
            nbins = i;
            break;
        }
    }
    lastbin = o * v * v * v - (nbins - 1) * binsize;

    psio->open(PSIF_DCC_ABCI3, PSIO_OPEN_OLD);
    psio->open(PSIF_DCC_ABCI2, PSIO_OPEN_OLD);
    psio->open(PSIF_DCC_ABCI4, PSIO_OPEN_NEW);

    psio_address abci2_addr = PSIO_ZERO;
    psio_address abci3_addr = PSIO_ZERO;
    psio_address abci5_addr = PSIO_ZERO;
    psio_address abci4_addr = PSIO_ZERO;
    for (size_t i = 0; i < nbins - 1; i++) {
        psio->read(PSIF_DCC_ABCI3, "E2abci3", (char *)&tmp[0], binsize * sizeof(double), abci3_addr, &abci3_addr);
        psio->read(PSIF_DCC_ABCI2, "E2abci2", (char *)&tmp2[0], binsize * sizeof(double), abci5_addr, &abci5_addr);
        C_DAXPY(binsize, 2.0, tmp, 1, tmp2, 1);
        // this is for the local triples
        psio->write(PSIF_DCC_ABCI4, "E2abci4", (char *)&tmp2[0], binsize * sizeof(double), abci4_addr, &abci4_addr);
    }
    psio->read(PSIF_DCC_ABCI3, "E2abci3", (char *)&tmp[0], lastbin * sizeof(double), abci3_addr, &abci3_addr);
    psio->read(PSIF_DCC_ABCI2, "E2abci2", (char *)&tmp2[0], lastbin * sizeof(double), abci5_addr, &abci5_addr);
    C_DAXPY(lastbin, 2.0, tmp, 1, tmp2, 1);
    psio->write(PSIF_DCC_ABCI4, "E2abci4", (char *)&tmp2[0], lastbin * sizeof(double), abci4_addr, &abci4_addr);
    psio->close(PSIF_DCC_ABCI2, 0);
    psio->close(PSIF_DCC_ABCI3, 1);
    psio->close(PSIF_DCC_ABCI4, 1);

    delete[] tmp;
    delete[] tmp2;
}
/**
  * OVOV in-core integral sort.  requires o^2v^2 doubles
  */
void SortOVOV(struct iwlbuf *Buf, int nfzc, int nfzv, int norbs, int ndoccact, int nvirt) {
    double val;
    size_t o = ndoccact;
    size_t v = nvirt;
    size_t fstact = nfzc;
    size_t lstact = norbs - nfzv;

    size_t lastbuf;
    Label *lblptr;
    Value *valptr;
    size_t nocc, idx, p, q, r, s, pq, rs;

    lblptr = Buf->labels;
    valptr = Buf->values;

    lastbuf = Buf->lastbuf;

    // available memory:
    size_t memory = Process::environment.get_memory();

    // 8 bytes for integrals
    size_t maxelem = memory / (sizeof(double));
    if (maxelem > o * o * v * v) maxelem = o * o * v * v;

    outfile->Printf("        CC integral sort will use %7.2lf mb\n", maxelem * (sizeof(double)) / 1024. / 1024.);
    if (maxelem < o * o * v * v) {
        throw PsiException("out of memory: o^2v^2 won't fit in core.", __FILE__, __LINE__);
    }

    auto *klcd = new double[o * o * v * v];
    memset((void *)klcd, '\0', o * o * v * v * sizeof(double));

    outfile->Printf("        Sort (IA|JB)........");
    /**
      * first buffer (read in when Buf was initialized)
      */
    for (idx = 4 * Buf->idx; Buf->idx < Buf->inbuf; Buf->idx++) {
        p = (size_t)lblptr[idx++];
        q = (size_t)lblptr[idx++];
        r = (size_t)lblptr[idx++];
        s = (size_t)lblptr[idx++];

        if (p < fstact || q < fstact || r < fstact || s < fstact) continue;
        if (p > lstact || q > lstact || r > lstact || s > lstact) continue;
        p -= fstact;
        q -= fstact;
        r -= fstact;
        s -= fstact;

        pq = Position(p, q);
        rs = Position(r, s);

        if (pq > rs) continue;

        val = (double)valptr[Buf->idx];
        klcd_terms_incore(val, pq, rs, p, q, r, s, o, v, klcd);
    }

    /**
      * now do the same for the rest of the buffers
      */
    while (!lastbuf) {
        iwl_buf_fetch(Buf);
        lastbuf = Buf->lastbuf;
        for (idx = 4 * Buf->idx; Buf->idx < Buf->inbuf; Buf->idx++) {
            p = (size_t)lblptr[idx++];
            q = (size_t)lblptr[idx++];
            r = (size_t)lblptr[idx++];
            s = (size_t)lblptr[idx++];

            if (p < fstact || q < fstact || r < fstact || s < fstact) continue;
            if (p > lstact || q > lstact || r > lstact || s > lstact) continue;
            p -= fstact;
            q -= fstact;
            r -= fstact;
            s -= fstact;

            pq = Position(p, q);
            rs = Position(r, s);

            if (pq > rs) continue;

            val = (double)valptr[Buf->idx];
            klcd_terms_incore(val, pq, rs, p, q, r, s, o, v, klcd);
        }
    }

    /**
      * write to disk
      */
    auto psio = std::make_shared<PSIO>();
    psio->open(PSIF_DCC_IAJB, PSIO_OPEN_NEW);
    psio->write_entry(PSIF_DCC_IAJB, "E2iajb", (char *)&klcd[0], o * o * v * v * sizeof(double));
    psio->close(PSIF_DCC_IAJB, 1);

    delete[] klcd;

    outfile->Printf("done.\n");
    outfile->Printf("\n");
}
}
}  // end of namespaces
