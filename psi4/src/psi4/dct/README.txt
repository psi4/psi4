When I load Q back, I only get one nonzero element.
This disagrees with the behavior of my own code.
I suspect that there is a bug in saving or loading the quantity.
I've traced this back to the set function. The current loop structure assumes rowspi = colspi, for each irrep.
This is likely a bad assumption, but fixing that may degrade performance when you can assume square.
Further, I don't understand why this set step is necessary. Why not use the block_matrix you already have?
