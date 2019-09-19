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


Sept 18th, 2019: 

Conclusions from yesterday's meeting with Chris is that the loading seems somewhat wrong but not entirely wrong. 

Need to determine the difference between bucket and slice!!

Genesis1.3 Manual, "Time Dependent Effects" pg 9: 
Bucket: electrons in 1 radiation wavelength within the electron bunch. USED FOR STEADY STATE
Slice: discritiztion of electron bunch in time. USED FOR TIME DEPENDENT

Loading should be changed to not be beamlets but Hammersley.
^^ Started working on this
 
Talked to Claudio - Definitions are correct, a slice is 1 rad length worth of electrons, electrons only rotate within
a slice, do not move slices.

To do:
1) Change loading to quick/random because beamlets is wrong. Update to Hammersley in future
2) Tapering?
3) Plot separatrix option 


Sept. 19, 2019: 

Making progress on adding tapering. 

Changing all "kai" to chi_1

Note: need to calculate gain length not just set it (can be calculated ~~~~ \lambda_u/4\pi\rho) need to find full formula

For benchmarking, need to get K's from Claudio's code to run both and check if results are similar? 
His code creates and saves outputs so that's good. Need to go back through those and figure out what is what. 

Have not yet run benchmark program in zfel - do that!

Update: Tapering is now included (basic). Need to: 
- test it
- include error messages for wrong tapering arrays??
- include other ways of inputing tapering (function??)

