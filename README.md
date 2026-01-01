For high level documentation, see below. For detailed derivations, numerical method descriptions, and full API reference, see docs.pdf (coming soon).

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

## 2. Quick Installation

1. Copy the following files into your project folder:

```

sim_API.py
math_API.py
special_API.py
config.ini

````

2. Import the functions you need. Example:

```
from math_API import integrate, deriv
from sim_API import run_sim
````

That’s it — you’re ready to use the API.

---


## 3. Working Example
A working example is provided in `sim.py`.

---

## 4. Running Tests

A partial test suite is included:

```
tests.py
```

This verifies that core functions behave correctly.

---

## 5. Further Reading

For detailed function documentation, equations, examples, and derivations, see `docs.pdf`. (coming soon.) 

---

**Enjoy solving ODEs without headaches!**
