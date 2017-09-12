{smcl}
{* 16 January 2015}{...}
{hline}
help for {hi:cdifdif}
{hline}

{title:Title}

{p 8 20 2}
    {hi:cdifdif} {hline 2} Difference-in-difference estimation with local spillovers

{title:Syntax}

{p 8 20 2}
{cmdab:cdifdif} yvar xvars [if] [in] [weight] {cmd:,} {cmd:close(}{it:varname}{cmd:)} {cmd:bandwidth(}{it:real}{cmd:)} [{it:options}]

{synoptset 20 tabbed}{...}
{synopthdr}
{synoptline}
{syntab :Options}
{synopt :{cmd:close(}{it:varname}{cmd:)}}
{p_end}
{...}
{synoptline}
{p2colreset}


{title:Description}

{p 6 6 2}
{hi:cdifdif} tests difference-in-differences regressions for local spillovers.  

{p 6 6 2}
For further details regarding the functionality of {cmd:cdifdif} or estimating a difference-in-differences model in the presence of spillovers more generally, refer to {it: Difference-in-Differences in the Presence of Spillovers} available at: {browse "https://sites.google.com/site/damiancclarke/research"}.


{marker examples}{...}
{title:Examples}

    {hline}
{pstd}Basic model for spillovers by geographic distance {break}

{phang2}{cmd:. cdifdif y i.year i.state postTreat, close(distance) bandwidth(2.5) vce(cluster state)}{p_end}

    {hline}



{marker results}{...}
{title:Saved results}

{pstd}
{cmd:cdifdif} saves the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(fit)}}Bayesian Information Criterion of final specification {p_end}
	

{marker references}{...}
{title:References}

{marker Clarke2015}{...}
{phang}
Clarke D.C., 2015.
{browse "https://sites.google.com/site/damiancclarke/research":{it:Difference-in-Differences in the Presence of Spillovers}.}
Manuscript.


{title:Also see}

{psee}
Online:  {manhelp regress_postestimation R} {manhelp regress_postestimation_time_series R}, {manhelp xtreg_postestimation XT}



{title:Author}

{pstd}
Damian C. Clarke, Department of Economics, University of Oxford. {browse "mailto:damian.clarke@economics.ox.ac.uk"}
{p_end}
