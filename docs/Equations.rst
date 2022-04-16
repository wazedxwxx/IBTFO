
 .. role:: cpp(code)
    :language: c++

 .. role:: f(code)
    :language: fortran

 
.. _Equations:



Equations
=========

Euler Equations
-------------------

IBTFO advances the following set of 2D Euler equations for the conserved state vector: :math:`\mathbf{U} = (\rho, \rho u, \rho v,\rho E):`

.. math::

`\mathbf{U}_t+\mathbf{F}_x+\mathbf{G}_y = 0
\mathbf{F} = (\rho u, \rho u^2 + p, \rho u v,u ( E + p))
 mathbf{G} = (\rho v, \rho u v, \rho v^2 + p ,v ( E + p)`
  


Here :math:`\rho, u, v`, and :math:`p` are the density, x-velocity, y-velocity and pressure, respectively. 
:math:`E = e + (u^2+v^2) / 2` is the total energy with :math:`e` representing the internal energy.

In the code, We define the number of equations in 'EQdefine.H'. We set it in the variable :math:`U(i,j,k)` , where k represents different equation.


Some notes:

* Regardless of the dimensionality of the problem, we always carry
  all 3 components of the velocity. You should always initialize all velocity components to zero, and
  always construct the kinetic energy with all three velocity components.

* There are ``NADV`` advected quantities, which range from :math:`{\tt
  UFA: UFA+nadv-1}`.  The advected quantities have no effect at all on
  the rest of the solution but can be useful as tracer quantities.

* There are ``NSPECIES`` species defined in the chemistry model, which range from :math:`{\tt UFS: UFS+nspecies-1}`.

* There are ``NAUX`` auxiliary variables, from :math:`{\tt UFX:UFX+naux-1}`. The auxiliary variables are passed into the equation
  of state routines along with the species.




