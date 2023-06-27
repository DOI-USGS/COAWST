# **Regional Ocean Modeling System (ROMS)**

![roms-jedi logo](https://www.myroms.org/trac/roms_src_600px.png)

The dynamical kernel of **ROMS** is comprised of four separate models, including
the nonlinear (**NLM**), perturbation tangent linear (**TLM**), finite amplitude
tangent linear (**RPM**), and adjoint (**ADM**). They are located in the
**Nonlinear**, **Tangent**, **Representer**, and **Adjoint** sub-directories,
respectively. The **TLM** and **ADM** were hand-coded from the discrete **NLM**
code using the recipes of Giering and Kaminski (1998). Therefore, any change to
its dynamical and numerical kernels will affect the symmetry of the **TLM** and
**ADM** operators. The discrete adjoint is exact and is defined relative to
the inner product that prescribes the L2-norm.

**ROMS** solves the free-surface, hydrostatic, flux form of the primitive
equations over variable bathymetry using stretched terrain-following in the
vertical and orthogonal curvilinear coordinates in the horizontal. The finite
volume grid is discretized on a staggered Arakawa C-grid. Detailed information
about its governing equation, numerical discretization, algorithms, usage, and
tutorials is available in the **WikiROMS** documentation portal at
**`www.myroms.org/wiki`**.

The official community version of **ROMS** is developed and maintained at Rutgers,
The State University of New Jersey, New Brunswick, New Jersey, USA. The scarlet
flag in the above logo indicates the location of our institution.

The **doxygen** version of **ROMS** is available at:
```
https://www.myroms.org/doxygen
```

Registered users of **ROMS** have access to:

- Official **svn** and **git** repositories
  ```
  svn checkout --username joe_roms https://www.myroms.org/svn/src/trunk <source_dir>

  git clone https://joe_roms@www.myroms.org/git/src <source_dir>
  ```

- Idealized and Realistic test cases repository
  ```
  svn checkout https://www.myroms.org/svn/src/test <test_dir>
  ```

- **ROMS** User's Forum:
  ```
  https://www.myroms.org/forum
  ```
