# IDC+S - Interacting Discrete Continua + Stokes
## A physics-based computational method for the solid and fluid mechanics modelling in fractures. 

This is a clean example of the solid and fluid mechanics simulator to model fracture compaction deformation with permeability calculation at each simulation step. The simulator is called Interacting Discrete Contunua and Stokes (IDC+S) and uses Virtual Element Method (VEM) and Finite Volume Method (FVM). 

Numerical grid can be selected from a lot of pre-defined grids using the natural Camel mudrock. Usually hard-coded grids are used in the corresponding publications.  

# Requirements

- MATLAB 2017 or above
- MRST-2017a or above (MATLAB Reservoir Simulation Toolbox https://www.sintef.no/projectweb/mrst/)

# Installation 

Once installed mrst, please place this github code to the "vem" module inside mrst. It needs some of the functions of it. 

# Examples

The main example file to run is compaction.m using MATLAB programming language. This file shows the simulation process in a (hopefully) comprehensible way. 

# References

Please cite the following papers if you are using IDC+S: 

1) 
``` 
  @article{kubeyev2022digital,
    title={Digital Image-Based Stress--Permeability Relationships of Rough Fractures Using Numerical Contact Mechanics and Stokes Equation},
    author={Kubeyev, Amanzhol and Forbes Inskip, Nathaniel and Phillips, Tomos and Zhang, Yihuai and Maier, Christine and Bisdom, Kevin and Busch, Andreas and Doster, Florian},
    journal={Transport in Porous Media},
    volume={141},
    number={2},
    pages={295--330},
    year={2022},
    publisher={Springer}
  }
``` 
2) 

``` 
@inproceedings{kubeyev2019geomechanics,
  title={Geomechanics numerical code for modelling contact in fractures and other discontinuities using virtual element method},
  author={Kubeyev, A and Maier, C and Doster, F},
  booktitle={53rd US Rock Mechanics/Geomechanics Symposium},
  year={2019},
  organization={OnePetro}
}

``` 
