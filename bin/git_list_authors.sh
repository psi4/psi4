#!/usr/bin/env tcsh

# I saw something analogous to this
#   in a presentation.  It doesn't 
#   tell you how much of the work
#   people did, but it gives you
#   a general idea.  Feel free to 
#   make changes or add username 
#   conversions.

git log --format=format:%an | sed \
    -e "s/jturney/Justin Turney/" \
    -e "s/andysim/Andy Simmonett/" \
    -e "s/rparrish/Rob Parrish/" \
    -e "s/loriab/Lori Burns/" \
    -e "s/rking/Rollin King/" \
    -e "s/Dummy (RAK)/Rollin King/" \
    -e "s/dhollman/David Hollman/" \
    -e "s/sherrill/David Sherrill/" \
    -e "s/crawdad/Daniel Crawford/" \
    -e "s/ugur/Ugur Bozkaya/" \
    -e "s/deprince/Eugene DePrince/" \
    -e "s/bmintz/Benjamin Mintz/" \
    -e "s/hohenstein/Ed Hohenstein/" \
    -e "s/Jet/Justin Turney/" \
    -e "s/marshall/Michael Marshall/" \
    -e "s/aevaughn/Alexander Vaughn/" \
    -e "s/Robert Parrish/Rob Parrish/" \
    -e "s/evaleev/Ed Valeev/" \
    -e "s/ssh2/Alexander Sokolov/" \
    -e "s/katierose/Katie Compaan/" \
    -e "s/magers/Brandon Magers/" \
    -e "s/avcopan/Andreas Copan/" \
    -e "s/kennedy/Matthew Kennedy/" \
    -e "s/massimo.malagoli/Massimo Malagoli/" \
    -e "s/jjwilke/Jeremy Wilke/" \
    -e "s/shrymc/Shane McNew/" \
    -e 's/Justin$/Justin Turney/' \
| sort | uniq -c | sort -r
