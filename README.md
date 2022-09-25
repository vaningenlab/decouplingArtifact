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

# Short description of scripts

The simulations are based on numerical evalation of the density operator for a 2-spin NH system in Liouville space,  
as described in the Allard et al [1998 JMR paper](https://doi.org/10.1006/jmre.1998.1509) with calculation of relaxation rates as described in the Helgstrand et al [2000 JBNMR paper](https://doi.org/10.1023/A:1008309220156)
The simulations do not include the effect of exchange.  



