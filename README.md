# decouplingArtifact

These are my GNU Octave scripts used for to simulate the 1H decoupling artifact 
in 15N CPMG experiments as described in the paper:  

*Removal of slow-pulsing artifacts in in-phase 15N relaxation dispersion experiments using broadband 1H decoupling*  

Soumya Deep Chatterjee, Marcellus Ubbink, Hugo van Ingen  
Journal of Biomolecular NMR, 2018, [paper](https://link.springer.com/article/10.1007/s10858-018-0193-2)

I put them now on Github so that anyone can make us of it.  

# How to use

1. Install GNU Octave from  
	a. *Windows users* [octave.org](https://octave.org)  
	b. *Mac users* [octave-app.org](https://octave-app.org)  
	c. *Linux users* using your package manager  

2. Download the simulations scripts 

3. Start octave  
	a. from the command line in the directory where you saved the simulation scripts  
	b. as gui, then navigate within the app to the directory where you saved the simulation scripts  

4. start the simulation by typing the name of the particular main simulation script "sim*xxxx*"

I verified that the scripts work using Mac OS 14 and Octave-app version 6.2.0.
It should also work under Windows with Octave version 7.2.0.s

# Short description of scripts

The simulations are based on numerical evalation of the density operator for a 2-spin NH system in Liouville space,  
as described in the Allard et al [1998 JMR paper](https://doi.org/10.1006/jmre.1998.1509) with calculation of relaxation rates as described in the Helgstrand et al [2000 JBNMR paper](https://doi.org/10.1023/A:1008309220156).
The simulations do not include the effect of exchange.  

+ simFigure1A | *main script to simulate the results shown in Figure 1A*
  | *comparison of matched 1H CW decoupling and single 1H power value, single train decoupling*.
+ simFigure1B | *main script to simulate the results shown in Figure 1B*
  | *predicted maximum artifact as function of protein size*.
+ simFigure1B | *main script to simulate the results shown in Figure 1C*
  | *predicted maximum artifact as function of magnetic field*.
+ simCPD | *main script to compare different composite pulse decoupling schemes*
  | *not used in the paper*.

+ ini | *initialization script w/ plotting control options*.
+ definePars* | *definition spin system and CPMG paramaters*.
+ makePlot* |*plotting file*.
+ buildRelaxationMatrix | *calculation of relaxation rates based on spin system, based on Helgstrand paper*.
+ LVM | *Liouvillian propagator matrix, based on Allard paper*.
+ LVM_CPD | *Liouvillian propagator during 1H CPD decoupling*
  | *w/ extensive checks on proper execution of CPD blocks*.
+ CPMG_CW | *in-phase 15N CPMG dispersion experiment w/ matched 1H decoupling*
  | *as decribed in Flemming et al 2008 paper*.
+ CPMG_ST_CW | *in-phase single-train 15N CPMG dispersion experiment w/ constant 1H CW decoupling*
  | *as decribed in Jiang et al 2008 paper*
+ CPMG_ST_CPD | *in-phase single-train 15N CPMG dispersion experiment w/ constant 1H CDP decoupling*
  | *as decribed in Chatterjee et al 2018 paper*
  | *CPD can be set to CW, 90-240-90, WALTZ, MLEV, SPA*
+ mymod | *custom modulo function used for Zuiderweg phase cycle CPMG*
+ Zuidersel | *custom modulo function used for Zuiderweg phase cycle CPMG*

