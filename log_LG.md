# Zfel development

Lipi Gupta started working on this on Sept. 16, 2019:

Initial findings (Sept. 17th, 2019): 

**Bunch loading is done incorrectly; should be 'quiet start' or Hammersley
- beamlets are wrong -- result in initial bunching of ~ 10E-4 (should be ~ 1/(number_of_electrons)

**Perhaps some issue with RF bucket vs. pondermotive bucket?
- particles should be able to travel from one p. bucket to the other (check periodic BC's)

**Need to output final theta and dgamma

**Tapering can be done in one of many ways; people use different methods:
- und_K as an array of values at different z steps (z is the coord along undulators)
- und_K as a continuous function of z
- und_K is recalculated to maintain a resonant phase (enter a phase, program calcs und_K for you)

**What was zfel checked against?

**Order of integration is-- First along z, then along s:
- Is the same as Genesis1.3 but the slippage is calculated very carefully in Genesis. 
Is this the case for zfel? Not sure..
- For a 1-D code, because it is easy to store a single E(z) per z step, could do it the other way to make accounting for slippage easier


******** Adding tapering is easy. Fixing the rest is harder/needs to be done first*********** 
