*! cdifdif: Estimating difference-in-differences in the presence of local spillovers
*! Version 1.0.0 august 14, 2017 @ 16:25:00
*! Author: Damian Clarke
*! Department of Economics
*! Universidad de Santiago de Chile
*! damian.clarke@usach.cl

cap program drop cdifdif
program cdifdif, eclass


vers 11.0
#delimit ;
syntax anything(id="Basic Model") [if] [in] [pweight fweight aweight iweight],
maxdist(real) distance(varlist numeric min=1 max=1)
[
 delta(real 1)
 tlimit(real 1.96)
 h(numlist >0)
 optimal
 crossfold
 kfold(integer 10)
 regtype(string)
 verbose
 stub(string)
 nogenerate
 *
]
;
#delimit cr

if length(`"`regtype'"')==0  local regtype regress
if length(`"`stub'"')==0 local stub _close


qui `regtype' `anything' [`weight' `exp'], `options'
tokenize `anything'
local y `1'

*-------------------------------------------------------------------------------
*--- (1) Using user-defined bandwidth
*-------------------------------------------------------------------------------
if length(`"`optimal'"')==0 {
    local l=0
    local u=`h'
    local vars
    local t=100
    qui {
        local nn=1
        while `t'>`tlimit' {
            local vars `vars' `proposal'
            if length(`"`nogenerate'"')!=0 {
                tempvar pd`nn'
                gen `pd`nn''=`distance'>(`l')&`distance'<=(`u')
                local proposal `pd`nn''
            }
            else {
                gen `stub'_`l'_`u'=`distance'>(`l')&`distance'<=(`u')
                local proposal `stub'_`l'_`u'
            }
            `regtype' `anything' `vars' `proposal' [`weight' `exp'], `options'
            local t = abs(_b[`proposal']/_se[`proposal'])
            local l=`l'+`h'
            local u=`u'+`h'
            local ++nn
        }
    }
    `regtype' `anything' `vars' [`weight' `exp'], `options'
    drop `proposal' `vars'
    ereturn local cvars `vars'
}

*-------------------------------------------------------------------------------
*--- (2) Using optimal bandwidth based on LOOCV or k-fold CV
*-------------------------------------------------------------------------------
if length(`"`optimal'"')!=0 {
    local RMSEmin = .
    local distopt = 0
    tempvar rmseVal
    qui gen `rmseVal' = .
    tempvar rmseNum
    qui gen `rmseNum' = .
    
    local iter=1
    local NnoS=0
    foreach k of numlist `delta'(`delta')`maxdist' {
        dis "Iteration `k'..."
        local l = 0
        local u = `k'
        local vars
        local proposal
    
        local t=100
        local nn=1
        while `t'>`tlimit' {
            local vars `vars' `proposal'
            if length(`"`nogenerate'"')!=0 {
                tempvar pd`nn'
                qui gen `pd`nn''=`distance'>(`l')&`distance'<=(`u')
                qui sum `pd`nn''
                local vcover = r(mean)
                while `vcover'==0 {
                    dis "No observations found between `l' and `u' from treatment"
                    local l=`l'+`k'
                    local u=`u'+`k'
                    local ++nn
                    qui gen `pd`nn''=`distance'>(`l')&`distance'<=(`u')
                    qui sum `pd`nn''
                    local vcover = r(mean)
                }
                local proposal `pd`nn''
            }
            else {
                qui gen `stub'_`l'_`u'=`distance'>(`l')&`distance'<=(`u')
                qui sum `stub'_`l'_`u'
                local vcover = r(mean)
                while `vcover'==0 {
                    dis "No observations found between `l' and `u' from treatment"
                    drop `stub'_`l'_`u'
                    local l=`l'+`k'
                    local u=`u'+`k'
                    qui gen `stub'_`l'_`u'=`distance'>(`l')&`distance'<=(`u')
                    qui sum `stub'_`l'_`u'
                    local vcover = r(mean)
                }
                local proposal `stub'_`l'_`u'
            }
            qui `regtype' `anything' `vars' `proposal' [`weight' `exp'], `options'
            local t = abs(_b[`proposal']/_se[`proposal'])
            local l=`l'+`k'
            local u=`u'+`k'
            local ++nn
        }
        if length(`"`vars'"')==0&`NnoS'>0 {
            qui replace `rmseVal' = `rmseZ' in `iter'
            qui replace `rmseNum' = `k'     in `iter'
            local ++iter
            local ++NnoS
            dis "RMSE for `k' is `rmseZ'" 
            drop `vars' `proposal' 
            continue
        }
        
        ****Leave-one-out Cross-Validation
        if length(`"`crossfold'"')==0 {
            qui count
            local max = r(N)
            tempvar diff
            qui gen `diff' = .
            tempvar yhat
            if length(`"`verbose'"')!= 0 dis "Calculating LOO CV"
            forval i = 1/`max' {
                if length(`"`verbose'"')!= 0 dis `i'
                qui `regtype' `anything' `vars' if _n!=`i' [`weight' `exp'], `options'
                qui predict `yhat' if _n==`i'
                qui replace `diff' = `y'-`yhat' if _n==`i'
                drop `yhat'
            }
        }
        ****k-fold Cross-Validation
        if length(`"`crossfold'"')!=0 {
            tempvar ksplit
            qui gen `ksplit' = ceil(runiform()*`kfold')        
            tempvar diff
            qui gen `diff' = .
            tempvar yhat
            if length(`"`verbose'"')!= 0  dis "Calculating k-fold CV"
            forval i = 1/`kfold' {
                if length(`"`verbose'"')!= 0 dis `i'
                qui `regtype' `anything' `vars' if `ksplit'!=`i' [`weight' `exp'], `options'
                qui predict `yhat' if `ksplit'==`i'
                qui replace `diff' = `y'-`yhat' if `ksplit'==`i'
                drop `yhat'
            }
        }
       
        tempvar rmse
        qui gen `rmse' = `diff'^2
        qui sum `rmse'
        dis "RMSE for `k' is " sqrt(r(mean))
        local rmseV = sqrt(r(mean))
        if sqrt(r(mean))<`RMSEmin' {
            local RMSEmin = sqrt(r(mean))
            local distopt = `k'
            local optcon treat `vars'
            local distmax = `l'-`k'
            local numband : word count `vars'
            local cvars `"`vars'"'
        }
        drop `vars' `proposal' `diff' `rmse'
        qui replace `rmseVal' = `rmseV' in `iter'
        qui replace `rmseNum' = `k'     in `iter'
        if length(`"`vars'"')==0&`NnoS'==0 {
            local rmseZ=`rmseV'
            local ++NnoS
        }
        local ++iter
    }

    local spill
    if length(`"`nogenerate'"')==0&`distmax'>0 {
        foreach u of numlist `distopt'(`distopt')`distmax' {
            local l = `u'-`distopt'
            qui gen `stub'_`l'_`u'=`distance'>(`l')&`distance'<=(`u')
            qui sum `stub'_`l'_`u'
            local vcover = r(mean)
            if `vcover'==0  drop `stub'_`l'_`u'
            if `vcover'!=0  local spill `spill' `stub'_`l'_`u'
        }
    }
    local nn=1
    if length(`"`nogenerate'"')!=0&`distmax'>0 {
        foreach u of numlist `distopt'(`distopt')`distmax' {
            local l = `u'-`distopt'
            tempvar pd`nn'
            qui gen `pd`nn''=`distance'>(`l')&`distance'<=(`u')
            qui sum `pd`nn''
            local vcover = r(mean)
            if `vcover'==0  drop `pd`nn''
            if `vcover'!=0  local spill `spill' `pd`nn''
            local ++nn
        }
    }
    `regtype' `anything' `spill' [`weight' `exp'], `options'
    dis "Optimal Distance calculated is `distopt'"
    dis "Maximum spillover distance is `distmax'"
    
    ereturn scalar max = `distmax'
    ereturn scalar h   = `distopt'
    ereturn local  closevars `"`cvars'"'
    ereturn scalar rmseOpt = `RMSEmin'

    
    #delimit ;
    twoway line `rmseVal' `rmseNum', scheme(s1mono) ytitle("RMSE")
    xtitle("Bandwidth (h)") lwidth(medthick) lcolor(cranberry);
    #delimit cr
}


end
