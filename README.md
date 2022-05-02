# TimingModels

Code to simulate models of how input timing affects velocity tuning of motion detectors.



## Circuit model notes

The circuit model code in this repository is based on code from [our previous work (Zavatone-Veth, Badwan, and Clark, *Journal of Vision*, 2020)](https://doi.org/10.1167/jov.20.2.2), which is available [here](https://github.com/ClarkLabCode/SynapticModel). However, it is not backward-compatible with that previously-released code. 

# BLM

Code to simulate responses of a pure delay model of motion detection: one that considers only relative delays between the two inputs, with no low pass filtering. Code also includes a version of the Barlow-Levick model, which considers different decay dynamics for the excitatory and inhibitory inputs. Simulations are for periodic and single points in space. 



