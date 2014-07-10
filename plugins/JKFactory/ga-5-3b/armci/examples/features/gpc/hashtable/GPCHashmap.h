/* $Id:  */
#ifndef GPCHASHMAP_H_
#define GPCHASHMAP_H_

#include <string>
using std::string;

#include "Hash_common.h"

class GPCHashmap {

    public:
        /**
         * Constructor
         */
        GPCHashmap();

        /**
         * Default Destructor
         */
        ~GPCHashmap();

        /**
         * creates a new hashmap (local)
         */
        void create();

        /**
         * destroys this hashmap
         */
        void destroy();

        /**
         * inserts elements into hashmap
         */
        void insert(const char *buf, size_t size);

        /**
         * get the global term IDs
         */
        void getGlobalIds(const char *buf, size_t bufsize,
                int *globalTermIds);

        /**
         * prints the local hashmap
         */
        void print();

        void rehash(int *size);

        /**
         * returns true if a hashmap already exists
         */
        static bool isCreated();

    private:
        VocabIntMap *mVocabMap;
        static short int sm_initialized;
};

#endif /* GPCHASHMAP_H_ */
