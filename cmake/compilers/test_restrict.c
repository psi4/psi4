int foo(int * TESTRESTRICTDEF i, int * TESTRESTRICTDEF j){return *i+*j;}
int main(int argc, char *argv[]){ int i=0; int j=0; return foo(&i,&j); }
