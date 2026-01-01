For a quick example, see below. For detailed derivations, numerical method descriptions, and full API reference, see docs.pdf (coming soon).

# Non-Linear ODE Solver API

Solve non-linear first-order vector ODEs of the form:

dy/dx = f(x, y) 

and higher-order ODEs via reduction of order. Supports multiple solvers including **explicit/implicit Euler**, **Runge–Kutta**, and **Gauss–Legendre Runge–Kutta** (orders 2, 4, and 6).

**Status:** Work in progress — functions and API may change.

---

## 1. Why use this API?

- Solve systems of ODEs easily without writing solver code from scratch  
- Supports multiple solver types with a consistent interface  
- Designed for clarity, simplicity, and reproducibility

---

