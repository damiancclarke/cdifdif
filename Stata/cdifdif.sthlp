{smcl}
{* 14 September 2017}{...}
{  hline}
help for {hi:cdifdif}
{hline}

{title:Title}

{p 8 20 2}
    {hi:cdifdif} {hline 2} Routine to estimate difference-in-differences models in the presence of treatment spillovers

{title:Syntax}

{p 8 20 2}
{cmdab:cdifdif} yvar xvars [if] [in] [weight]{cmd:,} {cmd:distance(}{it:varname}{cmd:)} {cmd:maxdist(}{it:real}{cmd:)} [{it:options}]

{synoptset 20 tabbed}{...}
{synopthdr}
{synoptline}
{syntab :Options}
{synopt :{cmd:delta(#)}}Defines fineness of grid to be searched when testing for optimal distance bandwidth; should be set based on units of measurement in the variable specified in {cmd:distance(}{it:varname}{cmd:)}{p_end}
{synopt :{cmd:tlimit(#)}}Sets level of significance required for marginal distance bin to be considered statistically significant; default is {cmd:tlimit(1.96)}{p_end}
{synopt :{cmd:kfold(#)}}Defines number of folds in k-fold cross-validation routine when determining RMSE associated with each distance bandwidth; default is {cmd:kfold(10)}{p_end}
{synopt :{opth stub(string)}}Specifies names for returned close to treatment variables.  By default, a series of variables are returned beginning with _close{p_end}
{synopt :{cmd:nogenerate}}Requests for no close to treatment variables to be returned, simply printing out the optimal regression model{p_end}
{synopt :{cmd:plotrmse}}Request that a graph is produced showing RMSE values for each spillover bandwidth tested{p_end} 
{synopt :{opth regtype(string)}}Specifies the regression model that is estimated with yvar and xvars.  This allows any type of model which is used in the baseline diff-in-diff model, including {help areg}, {help regress}, and so forth.  When not specified, regress is assumed {p_end}
{synopt :{opt *}}Any other options to be passed to the {help regress:regression}.  This allows for clustering standard errors, absorbing fixed effects (when areg is specified), and any other options permitted in {cmd:regtype(string)}  {p_end}
{synopt :{cmd:loocv}} Requests for leave-one-out cross-validation to be used rather than k-fold cross validation.  This is potentially slow when there are many observations in the model.{p_end}
{synopt :{cmd:verbose}}Requests additional output.  This is suggested when specifying {cmd:loocv} as it provides more details related to the state of advance of the program.{p_end}
{synopt :{cmd:nonoptimal}}Requests that the optimal bandwidth search process should not be used.  This is not a suggested option.  When specified, the spillover distances must be set by the user.{p_end}
{synopt :{cmd:h(#)}}Specifies the bandwidth for use in regressions.  This is only relevant when {cmd:nonoptimal} is specified, as otherwise {cmd:cdifdif} determines bandwidth optimally.{p_end}
{...}
{synoptline}
{p2colreset}


{title:Description}

	{p 6 6 2}
{hi:cdifdif} Implements the spillover-robust difference-in-differences (DD) estimator
described in Clarke (2017).  This routine is based on a baseline DD model of the
following form:

{p 8 12 2}y(i,t) = alpha + beta*Treat(i,t) + mu_i + phi_t + u_{i,t}      (1)

{p 6 6 2}
where Treat(i,t) is a variable of interest varying by time (t) and
units (i) or some other double difference setting, and mu and phi are
time and unit fixed effects.

{p 6 6 2}
The {hi:cdifdif} routine examines whether units which are untreated, but
spatially proximate to treated units were impacted by treatment spillovers.
Proximity must be measured by some distance to treatment variable.  As
described in Clarke (2017), {hi:cdifdif} tests for the presence of
"close to treatment" effects, augmenting equation (1) with a series of
Close variables:

{p 8 12 2}y(i,t) = alpha + beta*Treat(i,t) + gamma*Close(i,t) + mu_i + phi_t + u_{i,t}      (2)

{p 6 6 2}
This search is based on an optimal bandwidth procedure, in which the
distance variable is partiationed into mutually exclusive Close(i,t)
variables, which are then included in (2) up until treatment spillovers
are no longer observed.  The size of Close dummies are determined
using Leave-One-Out of k-fold Cross-Validation.  For further details regarding the functionality of {cmd:cdifdif} or estimating a
difference-in-differences model in the presence of spillovers more generally,
refer to {it: Difference-in-Differences in the Presence of Spillovers} available
at: {browse "https://sites.google.com/site/damiancclarke/research"}.

{p 6 6 2}
Estimating using {hi:cdifdif} simply requires that the user pass a
baseline difference-in-difference model as yvar and xvars, where
this model can be estimated using any of Stata's {help regress:regression}
commands, including {help regress}, {help areg}, {help xtreg}, and
so forth.  It is assumed that a full set of geographic baseline fixed
effects are included in the original model (ie area dummy variables, perhaps
absorbed using {help areg}) so
that when adding additional Close(i,t) variables, their initial fixed-effect
has already been included. Along with the baseline DD model, the name of the distance
variable should be passed, along with a maximal distance to consider
in tests of optimal spillover bandwidths when creating the matrix
of Close(i,t) variables.  The {cmd:distance(}{it:varname}{cmd:)} variable
should be equal to 0 if an area is itself treated, and the distance to the
nearest treated unit when an area is not itself treated.  In time periods
in which no treatment is available anywhere, distance can be missing, or
zero.

{p 6 6 2}
After running, {hi:cdifdif} displays the optimal model including any
spillovers (if present), and returns indicators for close to treatment
dummies (unless these are disabled using the {cmd:nogenerate} option).
{hi:cdifdif} is an e-class command, and also returns in the e(list)
the maximum spillover distance, optimal bandwidth size, Root Mean
Squared Error associated with optimal model and the name of all spillover
variables if present.

{title:Options}

{phang}
{cmd:distance(}{it:varname}{cmd:)} Specifies the variable which measures
distance to the nearest treatment unit for a given observation in a given
time-period.  This variable should be equal to zero when a unit
is treated.  In pre-treatment periods the variable can either be equal to
zero, or set as missing.

{phang}
{cmd:maxdist(}{it:real}{cmd:)} Scalar value which describes the maximum
possible spillover distance to test for a particular spillover bin.  
{hi:cdifdif} tests spillover bandwidths ranging from {cmd:delta(#)},
increasing in units of {cmd:delta(#)}, until reaching {cmd:maxdist}.
Note that {cmd:maxdist} does not imply that spillovers cannot be estimated
beyond the value set in {cmd:maxdist} via the iterative estimation procedure,
but rather that optimal bandwidths will not be considered beyond {cmd:maxdist}.

{phang}
{cmd:delta(#)} Changing {cmd:delta(#)} allows for the fineness of the grid to be
searched when testing for optimal distance bandwidth to be varied.  Indicating a
small value of {cmd:delta(#)} suggests that subsequent bandwidths tested should be
closely spaced, resulting in more bandwidth options considered when determining the
RMSE optimal bandwidth.  The units for {cmd:delta(#)} should be set based on units
of measurement in the variable specified in {cmd:distance(}{it:varname}{cmd:)},
as {cmd:delta(#)} partitions the distance variable into blocks of this size.

{phang}
{cmd:tlimit(#)} The spillover robust DD methodology considers models up until the
marginal distance bin is considered statistically significant. {cmd:tlimit(#)} sets
the level of significance required for marginal distance bin to be considered
statistically significant.  By default it is set at {cmd:tlimit(1.96)}.	

{phang}
{cmd:kfold(#)} When determining RMSE optimal models, k-fold cross-validation is
used by default (unless {cmd:loocv} is specified). Specifying {cmd:kfold(#)} allows
for the number of folds used when predicting y-hat to be varied.  By default this
is set at {cmd:kfold(10)}.

{phang}
{opth stub(string)} By default {hi:cdifcdif} returns the "close to treatment"
variables used in the optimal model starting as _close_dist1_dist2 (where dist1 and
dist2 refer to scalar values of distance cut-offs).  Setting stub replaces _close
with an alternative name for returned variables.

{phang}
{cmd:nogenerate} Requests for no close to treatment variables to be returned, simply
printing out the optimal regression model

{phang}
{cmd:plotrmse} Specifying {cmd:plotrmse} results in a graph showing RMSE values for
each spillover bandwidth tested between {cmd:delta(#)} and {cmd:maxdist(}{it:real}{cmd:)}.
This shows the RMSE for each model, with the lowest RMSE model providing the
otpimal bandwidth.

{phang}
{opth regtype(string)} Specifies the regression model that is estimated with yvar and
xvars.  This should provide the basic DD model used and how it is estimated in Stata.
This allows any type of model which is accepted as a Stata {help regress:regression},
including {help areg}, {help regress}, and so forth.  When not specified, regress is assumed.

{phang}
{opt *} Any other options passed to the {hi:cdifdif} command are assumed to be from
the original {help regress:regression}.  This allows for the clustering of standard
errors, absorbing fixed effects (when areg is specified), and the use of any other
desired options permitted in {cmd:regtype(string)}.

{phang}
{cmd:loocv} Requests for leave-one-out cross-validation to be used rather than k-fold
cross validation in determining the RMSE optimal spillover bandwidth.  This is
potentially slow when there are many observations in the model.

{phang}
cmd:verbose} Requests for additional output to be printed while the program progresses.
This is suggested when specifying {cmd:loocv} as it provides more details related to
the state of advance of the program, which may be slow when many observations are used.

{phang}
{cmd:nonoptimal} This option is generally not suggested.  It requests that the optimal
bandwidth search process should not be used.  When specified, the spillover distances
must be set by the user, and {hi:cdifdif} simply partitions the distance variable
into spillover bins, and tests for local spillovers..

{phang}
{cmd:h(#)} This Specifies the bandwidth for use in regressions.  This is only relevant
when {cmd:nonoptimal} is specified, as otherwise {cmd:cdifdif} determines bandwidth
optimally.



{marker examples}{...}
{title: Examples}

{pstd}Setup{p_end}
{phang2}{cmd:. webuse set http://www.damianclarke.net/data/}{p_end}
{phang2}{cmd:. webuse spilloverDGPs}{p_end}

{pstd} Run basic model for spillovers by geographic distance with k-fold cross-validation{p_end}
{phang2}{cmd:. cdifdif y1 treat i.time, distance(distance) maxdist(25) regtype(areg) abs(id) cluster(id)}{p_end}

{pstd} Replicate estimation using leave one out cross-validation (likely to be slow){p_end}
{phang2}{cmd:. cdifdif y1 treat i.time, distance(distance) maxdist(25) regtype(areg) abs(id) cluster(id) loocv verbose}{p_end}

{pstd} Estimate with k-fold cross-validation and produce graph of RMSE at search points {p_end}
{phang2}{cmd:. cdifdif y1 treat i.time, distance(distance) maxdist(25) regtype(areg) abs(id) cluster(id) plotrmse}{p_end}



{marker results}{...}
{title:Saved results}

{pstd}
{cmd:cdifdif} saves the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(max)}}Maximum distance up to which spillovers are identified{p_end}
{synopt:{cmd:e(h)}}Optimal search bandwidth{p_end}
{synopt:{cmd:e(rmse)}}Root Mean Squared Error associated with optimal bandwidth model{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(closevars)}}List of close to treatment spillover variables in optimal specification{p_end}
{p2colreset}{...}


{marker references}{...}
{title:References}

{marker Clarke2015}{...}
{phang}
Clarke, Damian (2017).
{browse "https://sites.google.com/site/damiancclarke/research":{it:Estimating Difference-in-Differences in the Presence of Spillovers}.}
Manuscript.


{title:Acknowledgements}
{pstd}
This project was supported by FONDECYT grant 11160200 from the Government of Chile.
{p_end}


{title:Author}

{pstd}
Damian Clarke, Department of Economics, Universidad de Santiago de Chile. {browse "mailto:damian.clarke@usach.cl":damian.clarke@usach.cl}
{p_end}
