***********************************************************************
*                   incompact3d-lmn Test case                         *
*                                                                     *
*         TEST: T1                                                    *
* WORK PACKAGE: WP1-TG2                                               *
*      SUMMARY: Test the advection of density by comparing with the   *
*               already implemented passive scalars.                  *
***********************************************************************

Description of test case:

If density is not included in the momentum/pressure equations, then it should behave as a passive
scalar. To test this a stepped profile of rho/phi is advected (without diffusion) by a uniform
velocity field. A potential concern is that the convection term of the continuity equation in LMN is
split using the chain rule to obtain an advective component and a diffusion-like component. This
therefore acts as a basic test that conservation isn't affected (for uniform velocity it shouldn't
be).

Setup:

- Disable diffusion in density and scalar subroutines
- Disable all subroutines updating velocity and pressure (ensures velocity field remains constant
  and speeds up test)
- Domain: 1x1x1
- Mesh: 128x9x8
- BC: 0-1-0
- U: 1.0
- rho/phi (t=0): 2.0, 0.25 <= x <= 0.75
                 1.0, otherwise
- dt: 0.00078125
- nsteps: 128,000

Passing criterion:

- Visual comparison of density/scalar fields: PASS
- Integral error: PASS
	int rho dV = 1.31253
	int phi dV = 1.31253

Status: PASS (13-OCT-2017)
