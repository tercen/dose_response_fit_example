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

Models used to fit the data|.
---|---
`nplr 4PL`      | yfit = B + (T-B) / (1 + 10^(b*(xmid-x)))
`nplr 5PL`      | yfit = B + (T-B) / (1 + 10^(b*(xmid-x)))^s
`drda 4PL`      | yfit = &alpha; + (&beta; - &alpha;) / ( 1 * exp(-&eta; * (x - &phi;)))
`drda 5PL`      | yfit = &alpha; + (&beta; - &alpha;) / ( 1 + &nu; * exp(-&eta; * (x - &phi;)))^(1/&nu;)

Where:
 `nplr`  
 B: Bottom asymptote  
 T: Top asymptote  
 b: Hill slope  
 xmid: x-coordinate at inflection point  
 s: asymmetric coefficient  
 
 `drda`  
  &alpha; : lower asymptote  
  &beta; : upper asymptote  
  &eta; : Steepness of the curve (growth rate)  
  &phi; : Related to the value of the function at x = 0  
  &nu; : Affects asymptote near which maximum growth occurs  

Further reading:

[Introduction to non-linear regression](https://www.statforbiology.com/nonlinearregression/usefulequations).  
[Nonlinear regression on Wikipedia](https://en.wikipedia.org/wiki/Nonlinear_regression).  
[nplr vignette](https://cran.r-project.org/web/packages/nplr/vignettes/nplr.pdf).  
[drda vignette](https://cran.r-project.org/web/packages/drda/vignettes/drda.pdf).  


##### See Also

[nonlinear_regression_operator](https://github.com/tercen/nonlinear_regression_operator)


