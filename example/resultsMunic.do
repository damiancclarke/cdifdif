/* resultsMunic.do v1.00         damiancclarke             yyyy-mm-dd:2017-07-23
----|----1----|----2----|----3----|----4----|----5----|----6----|----7----|----8

   This file replicates the results from Abouk and Adams (2013)' analysis of the  
impact of text messaging bans on fatal single-vehicle-single occupant accidents.
First the file replicates the original state-level results using their data and
specifications, and then extends to estimate municipal-level results with ident-
ical specifications, but municipal level measures.  Finally, it estimates using 
the spillover-robust DD methods described in Clarke (2017). In order to replica-
te the results shown in Clarke (2017), tables 2, 3, A2, A3 and A4, as well as f-
igure 2, this file should be run. The remaining figures and tables are from sim-
ulation results whose replication code is available in the "Monte Carlo" folder.

This file should be run in Stata, and the only lines that need to be changed are
the globals defined on lines 37 and 38 of this file, which point to the location
of data, as well as desired result output on the users' computer.

This file requires use of the estout and spmap ados made available on the SSC f-
from Boston College.  If these are not installed on the users' machine, they can
be installed by typing:
   ssc install estout
   ssc install spmap
in Stata prior to running this do file.

Questions or comments should be directed to damian.clarke@usach.cl.
*/

vers 11
clear all
set more off
cap log close
set matsize 5000

*-------------------------------------------------------------------------------
*--- (1) globals
*-------------------------------------------------------------------------------    
global DAT "~/investigacion/2014/Spillovers/data/AEJ/2013/AboukAdams/replicate"
global OUT "~/investigacion/2014/Spillovers/data/AEJ/2013/AboukAdams/replicate/results"

do cdifdif.ado
cap mkdir "$OUT"
log using "$OUT/resultsMunic.txt", replace text

*-------------------------------------------------------------------------------
*--- (2) Original Abouk and Adams analysis for comparison (state-level)
*-------------------------------------------------------------------------------    
use "$DAT/estdata"
#delimit ;
local conds treated==0&state!=2&txmsban==0
            treated==1&state!=2
            treated==1&state!=2&txmsban==0
            treated==1&state!=2&txmsban==1;
local names `" "Number of single-vehicle single-occupant accidents"
               "Population (annual)"
               "Unemployment rate (monthly)"
               "Proportion male (monthly)"
               "Real gas tax in 1983 cents (monthly)"
               "Sample size""';
#delimit cr

file open sums using "$OUT/StateSums.tex", write replace
replace rgastax=rgastax*100*100
tokenize accidentsvso pop unemp permale2 rgastax rgastax
local j = 1

foreach var of local names {
    dis "`var'"
    local fmt "%5.2f"
    if `j'==2  local fmt "%12.0gc"
    file write sums "`var'"
    foreach cond of local conds {
        sum ``j'' if `cond'
        local vo = string(r(mean),"`fmt'")
        if `j'==6 local vo = string(r(N), "%12.0gc")
        file write sums "& `vo'"
    }

    local ++j
    file write sums "\\" _n
}
file close sums

gen strongban   = second!=1 & agelimit!=1 & txmsban==1
gen weakban     = strongban==0 & txmsban==1
gen handheldban = HHBAN==1 & txmsban==1

local bans  strongban weakban handheldban
local conts lpop lunemp permale2 lrgastax
local wt    [aweight=pop]
local se    cluster(state) abs(state)
local yvar  laccidentsvso2


lab var strongban "Strong Ban"
lab var weakban   "Weak Ban"
lab var handheld  "Handheld Ban"
eststo: areg `yvar' strongban   `conts' laccidentmv2 i.time  `wt', `se'
eststo: areg `yvar' weakban     `conts' laccidentmv2 i.time  `wt', `se'
eststo: areg `yvar' handheldban `conts' laccidentmv2 i.time  `wt', `se'

#delimit ;
esttab est1 est2 est3 using "$OUT/StateDifDif.tex",
b(%-9.3f) se(%-9.3f) brackets nonumbers keep(strongban weakban handheldban)
nonotes mlabels(, none) style(tex) fragment replace noline label stats
(N r2, fmt(%9.0g %5.3f) label(Observations R-Squared))
starlevel ("*" 0.10 "**" 0.05 "***" 0.01) ;
#delimit cr
estimates clear


*-------------------------------------------------------------------------------
*--- (3) open Municipal data, summary statistics
*-------------------------------------------------------------------------------    
use "$DAT/estdata_munic", clear

file open sums using "$OUT/CountySums.tex", write replace
replace rgastax=rgastax*100*100
tokenize accidentsvso countypop unemp permale2 rgastax rgastax
local j = 1

foreach var of local names {
    dis "`var'"
    local fmt "%5.2f"
    if `j'==1  local fmt "%5.4f"
    if `j'==2  local fmt "%9.0gc"
    file write sums "`var'"
    foreach cond of local conds {
        sum ``j'' if `cond'
        local vo = string(r(mean),"`fmt'")
        if `j'==6 local vo = string(r(N), "%12.0gc")
        file write sums "& `vo'"
    }
    local ++j
    file write sums "\\" _n
}
file close sums


*-------------------------------------------------------------------------------
*--- (4a) Original table 3 at municipal level
*-------------------------------------------------------------------------------    
egen statecounty=concat(stateFIPS countyFIPS)
gen lcpop=log(countypop)
local conts lpop lunemp permale2 lrgastax
local wt    [aweight=countypop]
local se    cluster(state) abs(statecounty)
local yvar  laccidentsvso2
local bans  strongban weakban

eststo: areg `yvar' txmsban             `conts' i.time                `wt', `se'
eststo: areg `yvar' `bans'              `conts' i.time                `wt', `se'
eststo: areg `yvar' `bans' laccidentmv2 `conts' i.time                `wt', `se'
eststo: areg `yvar' `bans' laccidentmv2 `conts'        c.time#i.state `wt', `se'
eststo: areg `yvar' `bans' laccidentmv2 `conts' i.time c.time#i.state `wt', `se' 

lab var txmsban      "Texting ban in place"
lab var strongban    "$\times$ universally applied, primarily enforced"
lab var weakban      "$\times$ limited coverage/enforcement"
lab var lpop         "Log of population"
lab var lunemp       "Log of unemployment"
lab var permale2     "Percent male"
lab var lrgastax     "Log of gas tax"
lab var laccidentmv2 "Other accidents"

#delimit ;
esttab est1 est2 est3 est4 est5 using "$OUT/MunicipalReplicate.tex", replace
keep(txmsban `bans' `conts' laccidentmv2) booktabs b(%-9.4f) se(%-9.4f)
order(txmsban `bans' `conts' laccidentmv2) noobs label
title("Determinants of Fatal, Single-Vehicle, Single-Occupant Crashes
       (Municipal data)" \label{baselineDD}) nomtitles nodepvars
postfoot("Including 48 month fixed effects & Yes& Yes& Yes& No & Yes\\"
         "Including differential monthly trend & No& No& No&Yes& Yes\\"
         "\ \ \ for all states&&&&&\\"
         "\bottomrule\multicolumn{6}{p{16.4cm}}{\begin{footnotesize}    "
         "\textsc{Notes to Table \ref{baselineDD}}: Coefficients are    "
         "reported from weighted least squares regressions, where       "
         "weights are based on county population for 3,111 counties     "
         "over 48 months (following original specifications, Alaska     "
         "and Hawaii are not included). The dependent variable is the   "
         "natural logarithm of fatal accidents + 1. All specifications  "
         "include county by state fixed effects, and standard errors    "
         "are clustered by state. Specifications follow precisely those "
         "of \citet{AboukAdams2013}, however with county-level data in  "
         "place of state-level data."
         "*** Significant at the 1 percent level."
         "** Significant at the 5 percent level."
         "* Significant at the 10 percent level."
         "\end{footnotesize}}\end{tabular}\end{table}") style(tex);

#delimit cr
estimates clear

*-------------------------------------------------------------------------------
*--- (4b) Simple diff-in-diff models at county level
*-------------------------------------------------------------------------------
lab var strongban "Strong Ban"
lab var weakban   "Weak Ban"
lab var handheld  "Handheld Ban"
foreach ban in strongban weakban handheldban {
    areg `yvar' `ban'              `conts' i.time                `wt', `se'
    eststo: areg `yvar' `ban' laccidentmv2 `conts' i.time        `wt', `se'
    areg `yvar' `ban' laccidentmv2 `conts'        c.time#i.state `wt', `se'
    areg `yvar' `ban' laccidentmv2 `conts' i.time c.time#i.state `wt', `se' 
}

#delimit ;
esttab est1 est2 est3 using "$OUT/CountyDifDif.tex",
b(%-9.3f) se(%-9.3f) brackets nonumbers keep(strongban weakban handheldban)
nonotes mlabels(, none) style(tex) fragment replace noline label stats
(N r2, fmt(%9.0g %5.3f) label(Observations R-Squared))
starlevel ("*" 0.10 "**" 0.05 "***" 0.01) ;
#delimit cr
estimates clear

*-------------------------------------------------------------------------------
*--- (5a) Spillovers 'by hand'
*-------------------------------------------------------------------------------
local cont2 lpop lunemp permale2 lrgastax
local cdname Strong Weak HH
tokenize `cdname'

foreach ban in strong weak handheld {
    local noobs
    if `"`ban'"'!="handheld" local noobs noobs
    local spillovers
    eststo: areg `yvar' `ban'ban `spillovers' `cont2' i.time `wt', `se'
    local j=1
    foreach sp of numlis 10 20 30 40 {
        local spm5 = `sp'-10
        gen c_`ban'_`j'=closeDistance`1'>`spm5'&closeDistance`1'<=`sp'
        lab var c_`ban'_`j' "Close to `ban' ban [`spm5'-`sp')km"
        local spillovers `spillovers' c_`ban'_`j'
        eststo: areg `yvar' `ban'ban `spillovers' `conts' i.time `wt', `se'
        local ++j
    }
    #delimit ;
    esttab est1 est2 est3 est4 est5 using "$OUT/Spillover_manual_`ban'.tex",
    b(%-9.3f) se(%-9.3f) brackets nonumbers keep(`ban'ban `spillovers')
    nonotes mlabels(, none) style(tex) fragment replace noline label
    starlevel ("*" 0.10 "**" 0.05 "***" 0.01) `noobs';
    #delimit cr
    estimates clear
    drop `spillovers'
    macro shift
}


*-------------------------------------------------------------------------------
*--- (5b) Spillovers: cdifcdif min RMSE
*-------------------------------------------------------------------------------    
local cont2 lpop lunemp permale2 lrgastax i.time
local mdist 40
local delta 2
set seed 2711

#delimit ;
cdifdif `yvar' `cont2' handheldban `wt', distance(closeDistanceHH) `se'
h(11) regtype(areg) tlimit(1.64) maxdist(100) ;

eststo: cdifdif `yvar' handheldban `cont2' `wt', distance(closeDistanceHH) `se'
maxdist(`mdist') delta(`delta') regtype(areg) optimal crossfold kfold(10)
tlimit(1.64) verbose stub(hhclose);
graph export "$OUT/RMSE_handheld.eps", replace;
local handhelddist = e(distmax);
local kvars1=e(closevars);
dis length(`"`kvars1'"');

eststo: cdifdif `yvar' strongban `cont2' `wt', distance(closeDistanceSt) `se'
maxdist(`mdist') delta(`delta') regtype(areg) optimal crossfold kfold(10)
tlimit(1.64) verbose stub(strongclose);
graph export "$OUT/RMSE_strong.eps", replace;
local strongdist = e(distmax);

eststo: cdifdif `yvar' weakban `cont2' `wt', distance(closeDistanceWeak) `se'
maxdist(`mdist') delta(`delta') regtype(areg) optimal crossfold kfold(10)
tlimit(1.64) verbose stub(weakclose);
graph export "$OUT/RMSE_weak.eps", replace;
local weakdist = e(distmax);
local kvars3=e(closevars);
sum `kvars3';

lab var handheldban    "Treated (Handheld ban)";
lab var strongban      "Treated (Strongly enforced ban)";
lab var weakban        "Treated (Weakly enforced ban)";

esttab est2 est3 est1 using "$OUT/SpilloverDifDif.tex",
b(%-9.3f) se(%-9.3f) brackets nonumbers nolines
keep(strongban weakban handheldban `kvars1' `kvars2' `kvars3')
nonotes mlabels(, none) style(tex) fragment replace label stats
(N r2 h max rmseOpt, fmt(%9.0gc %5.3f %5.2f %5.2f %5.3f)
    label("\midrule Observations" "R-Squared" "Optimal Bandwidth (h)"
          "Maximum Spillover" "RMSE CV(h)" ))
starlevel ("*" 0.10 "**" 0.05 "***" 0.01) ;
#delimit cr
estimates clear


*-------------------------------------------------------------------------------
*--- (6) Basic Visualisation
*-------------------------------------------------------------------------------
preserve
replace closeDistanceHH     = 0 if handheldban==1
replace closeDistanceStrong = 0 if strongban==1
replace closeDistanceWeak   = 0 if weakban==1
replace closeDistanceAll    = 0 if txmsban==1


collapse (min) closeDistance*, by(countyFIPS stateFIPS)
rename countyFIPS COUNTYFP
rename stateFIPS STATEFP

merge 1:1 COUNTYFP STATEFP using "$DAT/c_data"

drop if STATEFP=="02"|STATEFP=="15"|STATEFP=="60"|STATEFP=="72"|STATEFP=="78"
drop if STATEFP=="66"|STATEFP=="69"
foreach vtype in All Strong Weak HH {
    sum closeDistance`vtype' if GEOID=="46047"
    replace closeDistance`vtype'=r(mean) if GEOID=="46102"
    replace closeDistance`vtype'=501 if closeDistance`vtype'==.
}


foreach vtype in All Strong Weak HH {
    format closeDistance`vtype' %5.2f

    #delimit ;
    spmap closeDistance`vtype' using "$DAT/c_coords_mercator", id(_ID)
    legstyle(2) osize(*0.01 ..) ocolor(none ..)
    legend(symy(*0.8) symx(*0.8) size(*1) rowgap(1) title("Distance"))
    clmethod(custom) clbreaks(0 0.1 100 200 300 400 500 5000) fcolor(Spectral)
    legend(label(2 "0.00 - 0.00") label(3 "0.00 - 100.00") label(8 ">500.00"));
    #delimit cr
    graph export "$OUT/minDistance`vtype'.eps", replace
}
restore


*-------------------------------------------------------------------------------
*--- (X) End
*-------------------------------------------------------------------------------
log close
