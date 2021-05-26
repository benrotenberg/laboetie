# laboetie

---

laboetie is a fluid dynamics code for chemical applications.  
It is based on the Lattice-Boltzmann algorithm to solve fluid dynamics equations.  
[more to come]

It is written by Maximilien Levesque¹ and Benjamin Rotenberg²  
¹ École Normale Supérieure and CNRS, Paris, France  
² CNRS, Sorbonne Université, UMR 8234, PHENIX, F-75005 Paris, France

---

## Please cite us!

We are researchers: our work is evaluated on the basis of citations to our publications. If you use Laboetie, please cite us (see list below)

------

## Science done with laboetie

  1. **Accounting for adsorption and desorption in lattice Boltzmann simulations**  
      Maximilien Levesque, Magali Duvail, Ignacio Pagonabarraga, Daan Frenkel and Benjamin Rotenberg  
      Phys. Rev. B 88, 013308 (2013)  
      http://dx.doi.org/10.1103/PhysRevE.88.013308    
  2. **Unexpected coupling between flow and adsorption in porous media**  
     Jean-Mathieu Vanson, François-Xavier Coudert, Benjamin Rotenberg, Maximilien Levesque, Caroline Tardivat, Michaela Klotz and Anne Boutin  
     Soft Matter 11, 6125-6133 (2015)  
     http://dx.doi.org/10.1039/C5SM01348H  
  3. **Transport and adsorption under liquid flow: the role of pore geometry**  
      Jean-Mathieu Vanson, Anne Boutin, Michaela Klotz and François-Xavier Coudert  
      Soft Matter 13, 875-885 (2017)  
     http://dx.doi.org/10.1039/C6SM02414A
  4. **Kinetic Accessibility of Porous Material Adsorption Sites Studied through the Lattice Boltzmann Method**  
      Jean-Mathieu Vanson, François-Xavier Coudert, Michaela Klotz and Anne Boutin  
      Langmuir 33, 1405-1411 (2017)  
      http://dx.doi.org/10.1021/acs.langmuir.6b04472  
  5. **Transient hydrodynamic finite size effects in simulations under periodic boundary conditions**  
      Adelchi J. Asta, Maximilien Levesque, Rodolphe Vuilleumier and Benjamin Rotenberg  
      Phys. Rev. E 95, 061301 (2017)  
      https://dx.doi.org/10.1103/PhysRevE.95.061301
  6. **Lattice Boltzmann electrokinetics simulation of nanocapacitors**  
      Adelchi J. Asta, Ivan Palaia, Emmanuel Trizac, Maximilien Levesque and Benjamin Rotenberg  
      J. Chem. Phys., 151, 114104 (2019)  
      https://dx.doi.org/10.1063/1.5119341  

---

## Installation instructions

The prefered way is to use CMake.  If you have CMake on your computer, then: 
```sh
git clone https://github.com/benrotenberg/laboetie
cd laboetie
mkdir build
cd build
cmake ..
make -j
```

If you don't have CMake:
```sh
git clone https://github.com/benrotenberg/laboetie
cd laboetie
make -j
```

## Git, github, issues etc.

[A successful git branching model](http://nvie.com/posts/a-successful-git-branching-model/)

[How to reports bugs effectively](http://www.chiark.greenend.org.uk/~sgtatham/bugs.html)

[How to write a git commit message](http://chris.beams.io/posts/git-commit/)


### Parallelism

The moment propagation is parallelized. It uses the OPENMP API. It is enabled by default (see `-fopenmp` in the Makefile).  
To disable openmp parallelism, remove `-fopenmp` from line 13 of Makefile.

By default, laboetie will use all the threads of your computer.

How do you control the number of threads used by OPENMP?  
You should export OMP_NUM_THREADS=3 in your terminal before executing laboetie if you want laboetie to use 3 threads:  
Thus, to compile and execute laboetie limited to 8 threads on my computer, I type in my terminal:
```bash
$ export OMP_NUM_THREADS=8
$ ./laboetie
```

Please note that the run time, user time, system time and cpu time printed to you by laboetie are not usable anymore when OMP_NUM_THREADS > 1. This is an issue that needs to be corrected.

Parallelism is implemented by dividing the system along almost independent slices along the z direction: if your system has {100,100,1} nodes along the x, y and z directions, respectively, then it is absolutely useless to multithread your job: use `export OMP_NUM_THREADS=1`.

For now, OPENMP in laboetie is memory bound. Preliminary speed-ups for systems of few hundred of nodes are as follow:
```
# number of threads, speed-up
1, 1
2, 1.9
3, 2.6
4, 3.1
6, 2.7
8, 2.7
```

I would recommand to use 2 or 4 threads only.

## Inputs (lb.in)

* `lx` Number of nodes in x direction
* `ly` Number of nodes in y direction
* `lz` Number of nodes in z direction
* `lbmodel` Lattice Boltzmann geometric model for velocities, e.g., D3Q19
* `timestepmax_for_PoissonNernstPlanck` Maximum number of timesteps in trying to find the equilibrium distribution of charged solutes. Should be 1 salt-free systems.
* `D_equil` number of steps for equilibrating charges (finding PB solution). 1 if charge (salt) free fluid.
* `t_equil` number of steps for equilibrating the flux without constraints
* `tmom` number of steps for equilibrating the flux without constraints
* `tmax` between tmom and tmax, moment propagation is done
* `geometryLabel`
                  #0 for a custom cell (written in geom.in)
                  #1 for a slit, i.e. two walls at z=zmin and z=zmax
                  #2 for a cylinder along Z. lx have to be equal to ly.
                  #3 for a body centered cubic cell with solid spheres in contact
                  # Other options are possible, check source code (supercell_definition.f90) 
* `initialSolventDensity` fluid density in LB units
* `f_ext` external force
* `charge_distrib` sol) charge distributed in the whole solid. int) charge distributed on interfacial nodes only
* `sigma` = 0.0 # charge distributed in solid
* `bjl` = 0.4 # bjerum length
* `lambda_D` = -2.0 # old debye_l
* `D_plus` = 0.05
* `D_minus` = 0.05
* `D_iter` = 1
* `tracer_Db` = 0.01   # bulk diffusion coefficient of the tracer
* `tracer_Ds` = 0.0   # surface diffusion coefficient of the tracer
* `tracer_z` = 0.0     # charge of the tracer
* `tracer_ka` = 0.1    # adsorption coefficient of the tracer
* `tracer_kd` = 0.01    # desorption coefficient of the tracer
* `c_tracer_0` = 0.1
* `left wall vel` = 0.0
* `right wall vel` = 0.0
* `Tprint_eq` = 1000
* `Tprint_run` = 10000
* `elec_slope` = 0.0 0.0 0.0 # external electric field
* `lncb_slope` = 0.0 0.0 0.0 # external gradient of salt concentration
* `f_gen` = 0.0 0.0 0.0



## Outputs

All outputs files are found in `output/`.

### supercell.xsf

A 3-dimensional representation of the supercell in [xsf format](http://www.xcrysden.org/doc/XSF.html).
*supercelf.xsf* can be opened with [VMD](http://www.ks.uiuc.edu/Research/vmd/): ```vmd -xsf output/supercell.xsf```.
The color code is:
* `pink` Solid nodes
* `green` Interfacial fluid nodes, i.e., fluid nodes close to a solid node
* `white` Non-interfacial fluid nodes
