# FuzzyC -- A simple implementation of Fuzzy logic in C

Originally just the product of a friendly competition with a friend to develop the fastest fuzzy inference system.
It later developed into a useful tool for control optimization for simple dynamic systems in my course work (combined
with a genetic fuzzy addition). It is now mostly a plaything, but still altogether useful. The code is mostly stable,
but changes may still occur. New features may be added, but they will be slow to appear if at all.

I should note that this project also represents my attempts to learn C, so there are many sins here. That being said,
I did take care to try to be as correct as I knew how. I may get the motivation in future months to go back and fix
many of my missteps, but I likely won't. A poor design choice made early on was to rely heavily on VLAs, which
is a GCC compiler extension. There is also a quadruple pointer in at least one place... that can't lead to
anything good. That being said, no memory leaks have been detected and except for the known bugs (see below),
the code does what is expected.

Known bugs:

- If the alpha exploration/mutation factor is anything but 0 in the blended crossover function, than crossover
  will eventually produce invalid individuals (individuals that exceed the assumed bounds). To fix this,
  bounds need to be checked during crossover, either as (lo, hi) params during crossover or as clamps directly
  after

