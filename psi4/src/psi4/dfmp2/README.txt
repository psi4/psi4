General Notes About the Code

1. There are several references to previous gradient papers in the comments.
  DiStasio : doi:10.1002/jcc.20604, J. Comput. Chem. 28, 839–856, (2007)
    "An Improved Algorithm for Analytical Gradient
    Evaluation in Resolution-of-the-Identity Second-Order
    Møller-Plesset Perturbation Theory: Application to
    Alanine Tetrapeptide Conformational Analysis"
      Most of our equations come from here.
  Wang : doi:10.1063/1.5100175; J. Chem. Phys. 151, 044118 (2019)
    "Analytic gradients for the single-reference
    driven similarity renormalization group
    second-order perturbation theory"
      This theory gives gradients for a family of methods that add
      another parameter to MP2, so these formulas apply to MP2 as well.
      Table II of this paper will be crucial in explaining a trick that
      does not appear in DiStasio. (See Sec. 3 here.)
  Aikens : doi:10.1007/s00214-003-0453-3; Theor. Chem. Acc. 110, 233-253 (2003)
    "A derivation of the frozen-orbital unrestricted open-shell
    and restricted closed-shell second-order perturbation theory
    analytic gradient expressions"
      This is the conventional-integral paper that DiStasio worked from to
      get their DF gradients. This sometimes clarifies obscure points in
      DiStasio, but they are not clear about one crucial point:
      If one block of a quantity is not given, but its adjoint is, is the first
      block zero, or the adjoint of the first block?
2. Conventions
    OPDM = One-particle density matrix
    EWDM = Energy-weighted density matrix
    All intermediates used in other functions are defined by comments starting DEFINITION.
    When a previous intermediate is updated, it is defined by comments starting UPDATE.
3. The Trick
    The RHF gradient code does NOT use the formulas of the DiStasio or Aikens papers.
    I can't find where the formulas come from, but they're very close to the DiStasio and
    Aikens formulas, but there are consistent deviations from the formulas that are
    accounted for by the addition of a "magic" term later on. In the comments, I refer to
    this as the "subtle trick." Grep for that if you want to see where it appears.
4. Cholesky Decomposition
    Several times, the code needs to contract a two-index MO quantity by the four-index electron
    repulsion integrals and end with a two-index AO quantity. If the two-index quantity were
    separable, this would be a perfect opportunity to use our efficient libfock/jk.cc technology.
    This two-index quantity is the one-particle density matrix, which is NOT separable.
    To remedy this, we Cholesky decompose the OPDM, which gets it in the desired form.
