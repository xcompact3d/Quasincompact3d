***********************************************************************
*                   incompact3d-lmn Test case                         *
*                                                                     *
*         TEST: T2a                                                   *
* WORK PACKAGE: WP1-TG2                                               *
*      SUMMARY: Test the conservation of density in a (somewhat)      *
*               complex flow with only advection.                     *
***********************************************************************

Description of test case:

If density is not included in the momentum/pressure equations, it should have no effect on the
result. The flow is the 2D mixing layer of, e.g. Golanski2005, for an initially constant density
field, the density should not increase over time.

Again a passive scalar is used for comparison purposes.

Setup:

- Disable diffusion in density equation (advection only)
- Domain: 30x60x1
- Mesh: 128x257x8
- BC: 0-1-0
- U1: 0.5
- U2: -0.5
- rho: 1.0
- Re: 400.0
- Pr: 1.0
- noise: 0.025
- dt: 0.05
- nsteps: 12,000
- iscalar: 1

Passing criterion:

Status: 
