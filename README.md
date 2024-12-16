# Fortran Programming Lab II (MTHL-3206)

## Introduction

This repository contains a collection of numerical methods implemented in modern Fortran, which were learned during the course Fortran Programming Lab II (MTHL-3206). These methods are widely used in scientific computing, engineering simulations, and mathematical problem-solving. Each method is implemented in Fortran 90/95/2003/2008 style (modern Fortran).

## Index

### 1. Solution of Algebraic Equations in Single Variables

- [Bisection Method](./01-bisection-method/main.f90)
- [Method of False Position](./02-false-position/main.f90)
- [Fixed-Point Iteration Method](./03-fixed-point-iteration/main.f90)
- [Newton-Raphson Method](./04-newton-raphson/main.f90)

### 2. Interpolation and Polynomial Approximation

- [Taylor's Polynomial Approximation](./05-taylor-polynomial/main.f90)
- [Newton's Forward Difference Interpolation](./06-newton-forward-difference/main.f90)
- [Newton's Backward Difference Interpolation](./07-newton-backward-difference/main.f90)
- [Newton's Divided Difference Interpolation](./08-newton-divided-difference/main.f90)
- [Lagrange's Interpolation](./09-lagrange-interpolation/main.f90)
- [Stirling's Central Difference Interpolation](./10-stirling-central-difference/main.f90)
- [Bessel's Central Difference Interpolation](./11-bessel-central-difference/main.f90)

### 3. Differentiation and Integration

- [Forward Difference](./12-forward-difference/main.f90)
- [Backward Difference](./13-backward-difference/main.f90)
- [Central Divided Difference](./14-central-divided-difference/main.f90)
- [Richardson Extrapolation](./15-richardson-extrapolation/main.f90)
- [Trapezoidal Rule](./16-trapezoidal-rule/main.f90)
- [Simpson's 1/3 Rule](./17-simpson-1-3-rule/main.f90)
- [Simpson's 3/8 Rule](./18-simpson-3-8-rule/main.f90)

## Reference used for algorithm

- Numerical Analysis (9th edition) by Richard L. Burden and J. Douglas Faires
- [Engineering Mathematics (Stirling's Central Difference Interpolation)](https://theengineeringmaths.com/wp-content/uploads/2017/11/interpolation-web.pdf)
- [LibreTexts Physics (Bessel's Central Difference Interpolation)](<https://phys.libretexts.org/Bookshelves/Astronomy__Cosmology/Celestial_Mechanics_(Tatum)/01%3A_Numerical_Methods/1.10%3A_1.10-Besselian_Interpolation>)
- [LibreTexts Mathematics (Richardson Extrapolation)](<https://math.libretexts.org/Bookshelves/Calculus/CLP-2_Integral_Calculus_(Feldman_Rechnitzer_and_Yeager)/04%3A_Appendices/4.03%3A_C%3A_More_About_Numerical_Integration/4.3.01%3A_C.1%3A_Richardson_Extrapolation>)
