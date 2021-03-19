Augmenting basis sets by extra functions
========================================

In many cases, one may want to add more basis function to a given
basis set in order to make it more suitable for properties it was not
designed to reproduce. Common examples are the improved description of
excited states or of very weakly bound electrons by adding more
diffuse functions into augmented basis sets, or the improved
description of the wave function close to the nucleus by adding more
steep functions into a basis set.

Manual extrapolation
--------------------

The traditional, well-known approach to augmenting basis sets with
extra functions is to rely on a rule of thumb to obtain the extra
functions. For instance, one can generate extra steep functions by
taking the largest exponent in the basis and multiplying it by a
factor of e.g. 3. If this function results in noticeable changes to
computed properties, the next step is to add another function by
multiplying the largest exponent in the original basis by 3*3 = 9, and
so on. The same procedure can also be used for diffuse exponents;
dividing the smallest exponent by e.g. 3 yields a new trial exponent.

This procedure has also been used by Zheng and coworkers to define
"minimally augmented" Karlsruhe basis sets
(http://doi.org/10.1007/s00214-010-0846-z), where one diffuse s and p
exponent are added to an original def2 basis set.

The manual extrapolation functionality is implemented in the
:func:`basis_set_exchange.manip.manual_augmentation` function.


Dunning-style extrapolation
---------------------------

The extrapolation works by taking the two outermost primitives in the
basis, :math:`X` and :math:`Y`. We assume :math:`X` to be the steepest
or most diffuse, and :math:`Y` is the second-steepest or second-most
diffuse. (If you want to do steep augmentation, you should also fully
decontract your basis set.)  Now, we can just add more primitives as
:math:`X(X/Y)`, :math:`X(X/Y)^2`, ... The procedure is the best
justified for even-tempered basis sets, where all exponents follow the
same pattern, but the procedure is used often also for fully optimized
basis sets.

Adding more diffuse functions to an already augmented basis set is
generally known as multiple augmentation, as described by Woon &
Dunning (https://doi.org/10.1063/1.466439). Adding a single
function per angular momentum shell to an aug set leads to a doubly
augmented (d-aug) set; adding two leads to a triply augmented (t-aug)
set, and so on.

This functionality is implemented in the :func:`basis_set_exchange.manip.geometric_augmentation` function.

Take, for example, He/aug-cc-pVTZ::

    BASIS "ao basis" SPHERICAL PRINT
    #BASIS SET: (7s,3p,2d) -> [4s,3p,2d]
    He    S
          2.340000E+02           2.587000E-03           0.000000E+00           0.000000E+00           0.00000000
          3.516000E+01           1.953300E-02           0.000000E+00           0.000000E+00           0.00000000
          7.989000E+00           9.099800E-02           0.000000E+00           0.000000E+00           0.00000000
          2.212000E+00           2.720500E-01           0.000000E+00           0.000000E+00           0.00000000
          6.669000E-01           4.780650E-01           1.000000E+00           0.000000E+00           0.00000000
          2.089000E-01           3.077370E-01           0.000000E+00           1.000000E+00           0.00000000
          0.0513800              0.00000000             0.00000000             0.00000000             1.0000000
    He    P
          3.044000E+00           1.000000E+00           0.000000E+00           0.00000000
          7.580000E-01           0.000000E+00           1.000000E+00           0.00000000
          0.1993000              0.00000000             0.00000000             1.0000000
    He    D
          1.965000E+00           1.0000000              0.00000000
          0.4592000              0.00000000             1.0000000
    END


The parameters are generated as below:

+----------------+--------------------+--------------------+--------------------+
|                | **S**              | **P**              | **D**              |
+----------------+--------------------+--------------------+--------------------+
| Smallest       | 2.089000E-01       | 7.580000E-01       | 1.965000E+00       |
| Exponents      +--------------------+--------------------+--------------------+
|                | 0.0513800          | 0.1993000          | 0.4592000          |
+----------------+--------------------+--------------------+--------------------+
| X              | 5.1380E-02         | 1.9930E-01         | 4.5920E-01         |
+----------------+--------------------+--------------------+--------------------+
| Y              | 2.089E-01          | 7.58E-01           | 1.965E+00          |
+----------------+--------------------+--------------------+--------------------+
| d-aug exponent | 1.2637E-02         | 5.2402E-02         | 1.0731E-01         |
+----------------+--------------------+--------------------+--------------------+
| t-aug exponent | 3.1082E-03         | 1.3778E-02         | 2.5077E-02         |
+----------------+--------------------+--------------------+--------------------+
| q-aug exponent | 7.6447E-04         | 3.6226E-03         | 5.8603E-03         |
+----------------+--------------------+--------------------+--------------------+

Therefore, He/q-aug-cc-pVTZ will be::

    BASIS "ao basis" SPHERICAL PRINT
    #BASIS SET: (7s,3p,2d) -> [4s,3p,2d]
    He    S
          2.340000E+02           2.587000E-03           0.000000E+00           0.000000E+00           0.00000000
          3.516000E+01           1.953300E-02           0.000000E+00           0.000000E+00           0.00000000
          7.989000E+00           9.099800E-02           0.000000E+00           0.000000E+00           0.00000000
          2.212000E+00           2.720500E-01           0.000000E+00           0.000000E+00           0.00000000
          6.669000E-01           4.780650E-01           1.000000E+00           0.000000E+00           0.00000000
          2.089000E-01           3.077370E-01           0.000000E+00           1.000000E+00           0.00000000
          0.0513800              0.00000000             0.00000000             0.00000000             1.0000000
    He    S
          1.264E-02              1.0
    He    S
          3.108E-03              1.0
    He    S
          7.645E-04              1.0
    He    P
          3.044000E+00           1.000000E+00           0.000000E+00           0.00000000
          7.580000E-01           0.000000E+00           1.000000E+00           0.00000000
          0.1993000              0.00000000             0.00000000             1.0000000
    He    P
          5.240E-02              1.0
    He    P
          1.378E-02              1.0
    He    P
          3.623E-03              1.0
    He    D
          1.965000E+00           1.0000000              0.00000000
          0.4592000              0.00000000             1.0000000
    He    D
          1.073E-01              1.0
    He    D
          2.508E-02              1.0
    He    D
          5.860E-03              1.

