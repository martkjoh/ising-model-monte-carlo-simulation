# ising-model-monte-carlo-simulation
Simulating the Ising model with Monte Carlo simulation and Metropolis sampling

Equilibration of 200x200 matrices at different temperatures

![pic](/figs/equilibration.png)


The python version runs the Metropolis arlogrithm on large chuncks of the spin-matrix, about half, simultaniously. This seemed to work, as it produces good graphs for the heat capacity and magnetic susceptibility. It causes a large speed-up, as it allows for the use of function on numpy-arrays, which ar *much* faster than for loops. The julia-algorithm only works on one spin at a time, which is tolerable as it is a compiled language. 

<p float="left">
  <img src="/figs/quantities/energy_density.png" width="400" />
  <img src="/ising-mod-data/energy.png" width="400" /> 
</p>

<p float="left">
  <img src="/figs/quantities/absolute_magnetization.png" width="400" />
  <img src="/ising-mod-data/abs_mag.png" width="400" /> 
</p>



<p float="left">
  <img src="/figs/quantities/heat_capacity.png" width="400" />
  <img src="/ising-mod-data/heat_cap.png" width="400" /> 
</p>


<p float="left">
  <img src="/figs/quantities/susceptibility.png" width="400" />
  <img src="/ising-mod-data/susc.png" width="400" /> 
</p>

These are pretty similar, but when comparing results for the interfacial tension, the differences are much larger

<p float="left">
  <img src="/figs/quantities/tension.png" width="400" />
  <img src="/ising-mod-data/tension.png" width="400" /> 
</p>
