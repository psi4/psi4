When I load Q back, I only get one nonzero element.
This disagrees with the behavior of my own code.
I suspect that there is a bug in saving or loading the quantity.
I've traced this back to the set function. The current loop structure assumes rowspi = colspi, for each irrep.
This is likely a bad assumption, but fixing that may degrade performance when you can assume square.
Further, I don't understand why this set step is necessary. Why not use the block_matrix you already have?

CONVENTIONS:
Much of the DCT code in Psi was written in ways that made perfect sense at the same time, but no longer
do now. Here, we document some more esoteric conventions:
* The amplitudes are often written as Lambda (the 2RDM cumulant). This is perfectly reasonable for most
  implemented theories, but the amplitudes no longer match Lambda for ODC-13. Whether the amplitudes should
  match lambda depends on your ansatz. Choosing the ansatz is a matter of ongoing research.
* In the -06 theories, the intermediates d and tau are identical. Some variable names say one when they
  (properly speaking) mean the other.
