int sim_sp(struct efp *,
	   const struct config *,
	   struct sys *);

int sim_grad(struct efp *,
	     const struct config *,
	     struct sys *);

int sim_cg(struct efp *,
	   const struct config *,
	   struct sys *);

int sim_nve(struct efp *,
	    const struct config *,
	    struct sys *);

int sim_nvt(struct efp *,
	    const struct config *,
	    struct sys *);
