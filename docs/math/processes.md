## Linear Elasticity

**Strong Form.**
Conservation of linear and angular momentum must be satisfied.
For smooth fields with applied body force field $f$, the equations are:

\begin{align}
  \newcommand{\Div}{\mathrm{div}\ }
  \label{eq:momentum}
  \Div \sigma + f &= 0 \\
  \sigma &= \sigma^T
\end{align}

For linear elasticity, the second order stress and strain tensors,
$ \sigma$ and $\epsilon$
are related by the fourth order stiffness tensor $\mathcal{C}$.
\begin{equation}
  \label{eq:stiffness}
  \sigma = \mathcal{C}:\epsilon
\end{equation}

**Boundary Conditions.**
In practice, boundary conditions are applied on sections of the boundary
either as displacements or tractions, or a combination of both. For example,
a box mesh has six natural boundary sections, corresponding to the min and max value
in each component direction.  You may want to apply zero displacements on the
$z=0$ section and apply a displacement only in the $z$-direction on the top
section, while the $x$ and $y$ components of applied traction are zero on the
same surface.

To make this more precise, the boundary is divided into a number ($n$) of
sections, $\Gamma_i$. Over each section, a displacement field or a traction
field can be applied for each component. On $\Gamma_i$
\begin{align}
  \label{eq:bcs}
	u \cdot e_j &= d_{ij}, \text{ or} \\
        t \cdot e_j &= g_{ij}, \text{ where $t=\sigma\ n$ is the traction} \\
\end{align}

**Null Space and Consistency Conditions**

Note that in the case of pure traction boundary conditions, there are consistency conditions.
\begin{align}
   \int f dx = \int_\Gamma t ds
\end{align}
Also, in that case, any rigid motion (6 degrees of freedom) is a solution when
$t=0$.

In fact, if any displacement component is unspecified over the whole boundary,
there is a consistency condition, as above, but for that component.

If your displacement boundary conditions eliminate rigid motions, you should
be good.  See the Implementation section for more details.

**Weak Form.**
Notationally, ${\cal H}$ is the function space of 3D distributions with
square integrable weak derivatives with zero displacements as implied by the
boundary conditions described above.
The surface is denoted by $\Gamma$ and $t$  represents the traction field
applied there.
For simplicity, $t$ is defined to have value 0 where displacement boundary condtions are in effect.

For the homogeneous problem, We seek the material displacement field
$u \in {\cal H}$ that satisfies the weak form of linear elasticity
(see equations above):
\begin{equation}
  \newcommand{\symm}{\mathrm{symm}\ }

  \int  \mathcal{C}:\symm\nabla u \cdot  \symm\nabla v\ dx = \int_\Gamma t \cdot v\ ds,
  \text{ for all test functions $v \in {\cal H}$}
\end{equation}
Here $dx$ is the volume integration element, and $ds$ is the surface integration element.
The test function $v$ is an arbitrary member of ${\cal H}$.
Note that due to the discontinuities in the stiffness tensor,
$\sigma$ will have jump discontinuities at grain boundaries, and
so the smooth form of the equilibrium equations does not apply,
but the weak form does.

**Discrete Form.**
We obtain the finite element formulation when we restrict the solution
space ${\cal H}$ to
a finite dimensional subspace ${\cal H}_h$ associated with the mesh.
The discrete equations are:
\begin{equation}
\newcommand{\symm}{\mathrm{symm}\ }
  \int  \mathcal{C}:\symm\nabla u_h \cdot  \symm\nabla v_h\ dx = \int_\Gamma t \cdot v_h\ ds,
  \text{ for all test functions $v_h \in {\cal H}_h$}
\end{equation}
where $u_h \in {\cal H}_h$ is the discrete displacement field.
The discrete test function $v_h$ is an arbitrary member of ${\cal H}_h$.


**Implementation.**
Our implementation calls for these coefficients:
* orientation field
* stiffness field
* traction boundary conditions
* displacement boundary conditions

We do not check for the consistency condition (at this point).  As for the
null space (rigid motions), the Krylov iterations give a solution with zero
projection on the null space. See [Bochev et al., 2005](https://doi.org/10.1137/S0036144503426074). That should probably be OK for problems with zero body
forces and zero tractions on the same components.  Point boundary conditions
have not been implemented.


### References

Bochev, Pavel, and R. B. Lehoucq. “On the Finite Element Solution of the Pure Neumann Problem.” SIAM Review 47, no. 1 (January 2005): 50–66. https://doi.org/10.1137/S0036144503426074.



## Thermal Conductivity

**Strong Form.**
The equations for steady heat flow are:
\begin{align}
  \label{eq:thermal}
  -\Div q + h &= 0 \\
\end{align}
where $q = -K \nabla \theta $ is the flux vector field,
$K$ is the (symmetric) thermal conductivity tensor
and $\theta$ is the temperature field, $h$ and is a volumetric heat source.

**Boundary Conditions.**
Boundary conditions will be some combination of fixed temperature ($\theta$) or
normal flux ($q_n = q \cdot n$).

**Weak Form.**
\begin{equation}
  \int  K \nabla \theta \cdot  \nabla v\ dx =
   \int h v\ dx + \int_\Gamma q_n\  v\ ds
\end{equation}
where $q_b$ is flux on the boundary, and $v$ is a test function.

### References

Lienhard, J. H., V and Lienhard, J. H., IV. A Heat Transfer Textbook, 6th ed. Cambridge MA: Phlogiston Press, 2024.


Ohta, Kenji, Yu Nishihara, Yuki Sato, Kei Hirose, Takashi Yagi, Saori I. Kawaguchi, Naohisa Hirao, and Yasuo Ohishi. “An Experimental Examination of Thermal Conductivity Anisotropy in Hcp Iron.” Frontiers in Earth Science 6 (November 6, 2018). https://doi.org/10.3389/feart.2018.00176.
