{smcl}
{* *! version 1.0.0  14January2022}{...}
{cmd:help wentropy}
{hline}

{title:Title}

{p2colset 5 24 26 2}{...}
{p2col :{cmd:wentropy} {hline 1} Implements cross-entropy weight calibration}{p_end}
{p2colreset}{...}

{title:Syntax}

{p 8 23 2}
{opt wentropy} {varlist} {ifin} {cmd:,} 
{opt NEWweight(varlist max=1 numeric)}
{opt constraints(matname)}
[{opt OLDweight(varlist max=1 numeric)}
{opt POPtotal(real 1)}
]


{title:Description}

{pstd}
{cmd:wentropy} Creates a vector of calibrated weights using cross-entropy. The weights are calibrated to ensure that the weighted means of variables match the constraints provided by the user. Under cross-entropy, the weights are
calibrated by the smallest amount possible given the constraints.


{title:Options}

{phang}
{opt NEWweight} New variable name that will be contain the cross-entropy calibrated weights.

{phang}
{opt constraints} Matrix K x 1 matrix containing the target means for the re-weighting. Matrix must have the same number of rows as the number of variables in varlist.

{phang}
{opt OLDweight} Variable containing the current set of weights which will be used as priors in the algorithm.

{phang}
{opt POPtotal} Population total to rescale weights.



{title:Example}

sysuse auto, clear
//constraints for foreign and price
mat A = (0.235\6500)
wentropy foreign price, old(weight) newweight(nw) constraints(A)
//Old means
sum foreign price [aw=weight]
//New means
sum foreign price [aw=nw]



{title:Authors:}

{pstd}
Paul Corral{break}
The World Bank - Poverty and Equity Global Practice {break}
Washington, DC{break}
pcorralrodas@worldbank.org{break}

{pstd}
Rodrigo Salcedo Du Bois{break}
rsalcedo@gmail.com{p_end}

{title:References:}

Golan, A., Judge, G. G., & Miller, D. (1996). Maximum entropy econometrics: robust estimation with limited data. Chichester [England], Wiley.




{pstd}
Any error or omission is the author's responsibility alone.


