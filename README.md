##### Description

Fits a non-linear dose-response curve using the 4PL or 5Pl models as implemented by either `nplr` or `drda` packages in R.

##### Usage

Input projection|.
---|---
`x-axis`        | numeric (log10), Concentration
`y-axis`        | numeric, Dose response
`row`           | 'Assay ID' (mandatory)
`colors`        | char, 'Sample type', Binary, QC or Standard


Input parameters|.
---|---
`Library`        | Either nplr or drda: library used to perform curve fitting
`Model`          | Either 4PL or 5PL: Number of parameters in the model used to fit

Output relations|.
---|---
`output_var`        | output relation
`Operator view`        | view of the Shiny application

##### Details

Details on the computation.

##### See Also

[nonlinear_regression_operator](https://github.com/tercen/nonlinear_regression_operator)


