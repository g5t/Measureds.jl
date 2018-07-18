# Measureds.jl

A package for Julia to handle measured values, which invariably have some uncertainty.

In general any measurement has some uncertainty due to statistics and/or instrumental precision.
This package defines a type which carries its uncertainty (as its variance) so that calculations 
based on the measured values can deal with the uncertainty as well.
In the case of simple calculations, this package has overloaded most simple mathematical functions to
correctly handle the uncertainty. For more complicated functions it is likely necessary to define a function
which calculates the uncertainty in the function based on the uncertainty in the input.

Measured values can be entered through a simple constructor as their value and variance, e.g., `Measured(100,10^2)`
for the value and uncertainty of counting 100 of something (uncertainty in counts is the square root of the number of counts, 10, and the variance is always the uncertainty squared).

Alternatively, the Julia non-standard string litteral functionality can be used, e.g., `m"100(10)"` where the
value and its uncertainty are given. This precise syntax is somewhat limiting in the case of non-integer values and uncertainties
but ten distinct forms of syntax are supported for string litteral input.

## String litteral input
In such a string, a relative uncertainty is indicated by one or two numbers in
parentheses following a value, e.g., `1.234(5)` indicates an uncertainty of plus or minus `5` in the
last digit of `1.234`; alternatively the same uncertainty can be specified as an absolute value as
`1.234±0.005`.

The ten valid syntax available for use with `@m_str` are:

| syntax | sym | rel |
|:---:|:---:|:---:|
| `val`(`err`)                                   | symmetric  | relative |
| `val`(`err`)×10`pow`                           | symmetric  | relative |
| `val`±`err`                                    | symmetric  | absolute |
| (`val`±`err`)×10`pow`                          | symmetric  | absolute |
| `val`×10`vpow`±`err`×10`epow`                  | symmetric  | absolute |
| `val`(`pos`,`neg`)                             | asymmetric | relative |
| `val`(`pos`,`neg`)×10`pow`                     | asymmetric | relative |
| `val`±(`pos`,`neg`)                            | asymmetric | absolute |
| (`val`±(`pos`,`neg`))×10`pow`                  | asymmetric | absolute |
| `val`×10`vpow`±(`pos`×10`ppow`,`neg`×10`npow`) | asymmetric | absolute |

where `val` represents its value, `err` symmetric uncertainty, `pow` overall exponent,
`vpow` value-specific exponent, `epow` uncertainty-exponent, `pos` positive uncertainty,
`neg` negative uncertainty, `ppow` postive-uncertainty-exponent, and `npow` negative-uncertainty
exponent. In the above table, each `±` can alternatively be `+/-` or `+-`, and each `×10` can
alternatively be `x10`, `e`, `E`, or `ᴇ`. The exponents `pow`, `vpow`, `epow`, `ppow`, and `npow`
must be integers (unsigned integers are taken to represent positive exponents) and can optionally
be raised characters, e.g., `×10⁻³`.


