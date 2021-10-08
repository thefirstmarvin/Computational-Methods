COMPUTATIONAL METHODS

A C++ program which solves the the first order wave equation below on a uniform grid with the prescribed initial and boundary conditions 

Wave Equation:
∂f/∂t + u ∂f/∂x = 0

Assume that a disturbance is introduced at t = 0 in a one-dimensional long tube of length L = 100m with x ∈ [−50,50] with both ends closed. The sets of imposed initial / boundary conditions are : 

SET 1 
t = 0 f(x,0) = 1/2(sign(x) + 1)
x = 0 f(−L/2, t) = 0
x = L f(L/2, t) = 1

SET 2
t = 0 f(x,0) = 1/2 exp(−x2)
x = 0 f(−L/2, t) = 0
x = L f(L/2, t) = 0

Using the following methods:
•Explicit Upwind FTBS (Forward time, Backward space)
•Implicit Upwind FTBS (Forward time, Backward space)
•Lax-Wendroff
•Richtmyer multi-step 

The analytical solution of this problem, subject to the imposed initial and boundary conditions, is respectively:

SET 1: f(x, t) =1/2(sign(x−1.75t) + 1)
SET 2: f(x, t) =1/2exp(−(x−1.75t)2)

where u, the speed of sound, is 1.75m/s.   
