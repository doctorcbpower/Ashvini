"""
Physics components of the general chain (milestone M2 of docs/SUPERSET_DESIGN.md). Each module holds one
piece of the model, with the equation it implements in its docstring, so it can be read, replaced or
extended on its own. `integrators.grid.GridIntegrator` composes them into the step-by-step integration of
`main.run_forest`.

    ode            closed-form update of dy/dt = forcing - decay * y
    star_formation free-fall star formation rate coefficient
    inflow         gas accretion rate on to the galaxy (includes the UV suppression)
    delay          which star formation drives the supernova winds: none, present or delayed
    wind model        supernova mass loading (how much gas is blown out per unit star formation)
    chemistry      metals and dust
    blackholes     seeding, growth and AGN wind of the central black hole
"""
