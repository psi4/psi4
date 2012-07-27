char * efp_strndup(const char *s, size_t n);

int parse_config(const char *,
		 struct config *,
		 struct sys *);

void free_config(struct config *,
		 struct sys *);
