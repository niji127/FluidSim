# FluidSim
SPH solver based on CUDA
### Features:   
- CUDA 12.1 based
- Fluid-structure interaction
- Weakly compressible SPH for fluid
- Total lagrangian SPH for solid
- Coupling force based on the Riemann-SPH
- Spatial partitioning for neighborhood search
- Cell lists data structure
- Multi-resolution time stepping
- .vtk output
### Configuration
- physics: gravity
- fluid: density, sound_speed, gamma, viscosity
- solid: density, youngs_modulus, poisson_ratio
- domain: domain_min, domain_max
- time: dt, start_step, end_step, solid_sub_step
- kernel: smoothing length, particle diameter, impact_length
### Example
#### Fluid-Structure Interaction (fluid: 4M particles, solid: 27k particles) 
surface reconstruction using splashsurf (https://github.com/InteractiveComputerGraphics/splashsurf)
![output_video](https://github.com/niji127/FluidSim/assets/152270816/5fe22c65-3e71-4f5f-ab66-5a7a3029a2fa)
