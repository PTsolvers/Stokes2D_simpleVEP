# Stokes2D visco-elasto-plastic (simple)
Visco-elastic rheology for 2D Stokes solvers and continuum mechanics.

![](docs/output_ve_vm_phi.png)

This repository contains 2D iterative visco-elastic incompressible and single phase Stokes solvers to (1) resolve pressure, velocity and visco-elastic stress distribution around a buoyant ductile spherical inclusion and (2) capture the visco-elastic stress build-up in a homogeneous sample.

## Content
* [Julia Codes](#julia-codes)
* [Experiment results](#experiment-results)
* [Extra material](#extra-material)

## Julia codes
The Julia codes implementing 2D Stokes equations and visco-elastic shear rheology:
- [`Stokes2D_ve_grav.jl`](Stokes2D_ve_grav.jl) resolves a visco-elastic buoyant inclusion setup;
- [`Stokes2D_ve_bench.jl`](Stokes2D_ve_bench.jl) captures the visco-elastic stress build-up shearing a homogenous bloc.
- [`Stokes2D_ve_pureshear.jl`](Stokes2D_ve_pureshear.jl) captures the visco-elastic stress build-up shearing a visco-elastic inclusion.
- [`Stokes2D_ve_vm.jl`](Stokes2D_ve_vm.jl) resolve brittle failure of a bloc containing a visco-elastic inclusion.

## Experiment results
The rise of a buoyant and ductile inclusion generates, among others, pressure deviation from the hydrostatic gradient, vertical (y) velocity field and vertical normal stress as depicted in the following figure:

![](docs/output_ve_grav.png)

The visco-elastic stress build-up benchmark captures stress build up while applying pure shear on a homogeneous visco-elastic body. The current non-dimensional configuration is expected to reach a maximal stress level of 2.0 once the elastic build-up is completed, recovering the viscous solution. The figure depicts the horizontal and vertical velocity fields, and the stress build-up curve as function of time, matching the analytical solution (red line):

![](docs/output_ve_bench.png)

---

🚧 **On-work section**

Shearing a visco-elastic inclusion results:

![](docs/output_ve_pureshear.png)

Shearing a visco-elasto-plastic inclusion results.

- Von Mises plasticity:

![](docs/output_ve_vm.png)

- Adding friction angel 👼:

![](docs/output_ve_vm_phi.png)


## Extra material
- A succinct [intro to continuum mechanics](docs/intro_continuum_mechanics.pdf) as written up by a former colleague from the Uni Lausanne.
- Some [slides](docs/visco-elast_schmalholz_unil.pdf) from a mechanics course given in Earth sciences at Uni Lausanne by Prof. S. Schmalholz.
