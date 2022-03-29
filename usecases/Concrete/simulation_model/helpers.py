import dolfin as df


class Parameters(dict):
    """
    Dict that also allows to access the parameter
        p["parameter"]
    via the matching attribute
        p.parameter
    to make access shorter
    """
    def __getattr__(self, key):
        return self[key]

    def __setattr__(self, key, value):
        assert key in self
        self[key] = value

    def __add__(self, other):
        if other == None:
            dic = self
        else:
            dic = Parameters({**self, **other})
        return dic


class Sensors(dict):
    """
    Dict that also allows to access the parameter
        p["parameter"]
    via the matching attribute
        p.parameter
    to make access shorter
    """
    # TESTING a sensor dictionary
    def __getattr__(self, key):
        return self[key]

    def __setattr__(self, key, value):
        assert key in self
        self[key] = value



# helper functions for quadrature spaces
def set_q(q, values):
    """
    q:
        quadrature function space
    values:
        entries for `q`
    """
    v = q.vector()
    v.zero()
    v.add_local(values.flat)
    v.apply("insert")


class LocalProjector:
    def __init__(self, expr, V, dxm):
        """
        expr:
            expression to project
        V:
            quadrature function space
        dxm:
            dolfin.Measure("dx") that matches V
        """
        dv = df.TrialFunction(V)
        v_ = df.TestFunction(V)
        a_proj = df.inner(dv, v_) * dxm
        b_proj = df.inner(expr, v_) * dxm
        self.solver = df.LocalSolver(a_proj, b_proj)
        self.solver.factorize()

    def __call__(self, u):
        """
        u:
            function that is filled with the solution of the projection
        """
        self.solver.solve_local_rhs(u)
