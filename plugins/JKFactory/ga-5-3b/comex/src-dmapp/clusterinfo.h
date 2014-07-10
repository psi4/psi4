#ifndef CLUSTERINFO_H_
#define CLUSTERINFO_H_

/* consider up to HOSTNAME_LEN characters in host name */
#define HOSTNAME_LEN 64

typedef struct {
    int master;
    int nslave;
    char hostname[HOSTNAME_LEN];
} armci_clus_t;

extern armci_clus_t *armci_clus_info;
extern int armci_me;
extern int armci_nproc;
extern int armci_nclus;
extern int armci_clus_me;
extern int armci_master;
extern int armci_clus_first;
extern int armci_clus_last;

extern int armci_clus_id(int p);
extern void armci_init_clusinfo();

#endif /* CLUSTERINFO_H_ */
