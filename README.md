# NS-solver
Laminar flow Navier-Stokes (NS) solver using the NS equations in vorticity form and finite differences.
This implementation is also in cylindrical coordinates (2D) to simulate fluid flow on a pipe 
(which was checked against a simple analytic solution).

The code could also be adapted to other cases with simple geometries cases such as the [lid driven cavity]
(http://www.cfd-online.com/Wiki/Lid-driven_cavity_problem).

The vorticity formulation of the NS equations was used to prevent the pressure checkerboarding issue that happens
when trying to apply finite differences to the NS equations.

The solution for each iteration until convergence is carried out as follows:
* The vorticity is diffused
* The radial velocity is implicitly calculated from the vorticity using the [Bi-CGSTAB algorithm]
(https://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method)
* Finaly the axial velocity and new vorticity is calculated from continuity

The quivers.m code used to graph the velocity plot is a slightly modified version of Bertran Dano's [code] (http://www.mathworks.com/matlabcentral/fileexchange/24482-quivers/content/quivers.m).
