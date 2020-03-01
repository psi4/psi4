* There are currently three very similar families of dpdfiles that look like Fock matrices:
  Fock <X|X> - An intermediate used during construction of the Fock matrix. Just the J/K part.
  F <X|X>    - The Fock matrix without diagonal elements
  FD <X|X>   - The Fock matrix
* Ideally, those can be condensed. The F uses in the amplitude conditions will change to FD upon changing the conditions from = to += form
  As for Fock uses, reconsider if that's the best way to construct the Fock matrix. If it is, then rename it.
