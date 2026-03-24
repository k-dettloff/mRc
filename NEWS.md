# mRc 0.1.2 (2026-03-24)

* Made further speed improvements using matrix multiplication and optimized memory allocation.
* Added validation checks to ensure integer-valued inputs.

# mRc 0.1.1 (2026-03-23)

* Optimized `closedCI()` using vectorized matrix operations for significantly faster performance.
* Added correction to remove excess decimal value from adjusted Schumacher-Eschmeyer estimator.
* Applied clamping logic to estimates and confidence bounds to prevent impossible values in low-sample scenarios.
* Refactored code to follow R best practices (spacing, assignment operators, and `stats` namespace imports).
