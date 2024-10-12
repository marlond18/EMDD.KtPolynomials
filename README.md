&nbsp; [![Nuget](https://img.shields.io/nuget/v/EMDD.KtPolynomials)](https://www.nuget.org/packages/EMDD.KtPolynomials/)[![Nuget](https://img.shields.io/nuget/dt/EMDD.KtPolynomials)](https://www.nuget.org/stats/packages/EMDD.KtPolynomials?groupby=Version&groupby=ClientName&groupby=ClientVersion)[![Deploy to Nuget](https://github.com/marlond18/EMDD.KtPolynomials/actions/workflows/nuget.yml/badge.svg)](https://github.com/marlond18/EMDD.KtPolynomials/actions/workflows/nuget.yml)[![Run Tests](https://github.com/marlond18/EMDD.KtPolynomials/actions/workflows/runTest.yml/badge.svg)](https://github.com/marlond18/EMDD.KtPolynomials/actions/workflows/runTest.yml)
&nbsp; 

----------------
# EMDD.KtPolynomials
a library for Polynomial equation manipulations

## Requirements

[.Net 5.0.102 sdk](https://dotnet.microsoft.com/download/dotnet/5.0) or greater

## Nuget Package Usage

https://www.nuget.org/packages/EMDD.KtPolynomials/

`<PackageReference Include="EMDD.KtPolynomials" Version="*.*.*" />`

## How it works
- this library was created in order to make representation of the polynomial equations. ```KtPolynomial``` is a class that represents a polynomial of any degree.


## Usage
A polynomial can be represented by ```KtPolynomial``` through the ```KtPolynomial.Create``` method, with an array of [```EMDD.KtNumeric.Number```](https://github.com/marlond18/EMDD.KtNumerics) as a parameter, representing the coefficient of the polynomial.

So, `x^3+2x^2-1=0` is instantiated as

```c#
KtPolynomial.Create(3,2,-1) or KtPolynomial.Create(new double[] { 3,2,-1 })
```

### ```KtPolynomials``` Methods
#### Basic Math Ops
- addition, subtraction, multiplications with other polynomials, multiplications with constants, division by a constant
#### Differentiation and Integration
- ```KtPolynomials``` can be differentiated and integrated
#### Root Extraction/ Factorization
- ```KtPolynomials``` has methods for root extraction and factorization using:
  * JenkinsTraub
  * Laguerre
  * Weirstrass
