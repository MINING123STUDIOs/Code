def newton_solve(F, x0, tol=1e-9, max_iter=20):
    """
    Solve F(x) = 0 for vector x using Newton's method
    with numerical Jacobian approximation.

    Parameters
    ----------
    F : function(x) -> ndarray
        Function returning residual vector.
    x0 : ndarray
        Initial guess.
    tol : float
        Convergence tolerance on ||F(x)||.
    max_iter : int
        Maximum Newton iterations.

    Returns
    -------
    x : ndarray
        Approximate root of F(x).
    """

    x = x0.astype(float).copy()
    n = len(x)

    for _ in range(max_iter):
        Fx = F(x)
        if np.linalg.norm(Fx, ord=2) < tol:
            return x

        # Numerical Jacobian
        J = np.zeros((n, n), dtype=float)
        eps = 1e-8

        for i in range(n):
            x_pert = x.copy()
            x_pert[i] += eps
            J[:, i] = (F(x_pert) - Fx) / eps

        # Solve J Δ = -F
        delta = np.linalg.solve(J, -Fx)
        x += delta

        if np.linalg.norm(delta, ord=2) < tol:
            return x

    return x   # last iterate (same behavior as many solvers)
