*! blp 2.0 04april2015
*! David Vincent, email: davivincent@deloitte.co.uk
program blp, eclass
version 10.1
set more off

if replay() {
	if "`e(cmd)'"!="blp" {
		error 301 
	}
	ereturn display
}
else {

syntax varlist [if] [in], 	///
	stochastic(string) 	///
	endog(string) 		///
	markets(name) [ 		///
	draws(integer 200)	///
	Optinsta			///
	Optinstb(string)		///
	Itera				///
	Iterb(string)		///	
	demofile(string)		///
      tolin(real 10E-15) 	///
	tolout(real 10E-12)	///
	initsd(string)		///
	initdemo(string)		///
	elast(string)		///
	NOIsily			///
	NOCons 			///
	ROBUSTWeight 		///
	ROBUST			///
	RANdom			///
	Burn(string) ]		///
			
	
marksample touse								
gettoken y_var exog_vars: varlist					
local exog_vars2 `exog_vars'						

//indicators					
local use_optinst = "`optinsta'"!="" |"`optinstb'"!=""	//use optimal instruments 
local optinst = cond(`use_optinst', "`optinstb'", "")		//optimal instruments (varlist) for endogenous vars
local robustweight = "`robustweight'"!=""				//use robustweighting matrix (i.e. permits E[u_jt,u_kt]!=sig_jk,t)
local robust = "`robust'"!=""						//use robust standard errors 
local noisily = "`noisily'"!=""					//display contraction mapping by market 
local random = "`random'"!=""						//random instead of halton draws								


//exit if share variable outside interval (0,1)
qui sum `y_var'
if r(min)<=0 | r(max)>=1{
	di as error "share variable must be within the interval 0,1"
	exit 498
} 

//exit if iterative estimator specified when no iteration required
if ("`itera'"!="" | "`iterb'"!="") & !`robustweight' & !`use_optinst' {
	di as err "may only specify {opt iterate()} for options {opt robustweight} or {opt optinst()}" 
	exit 198
}

//iterative estimator options
if "`itera'"!=""{
	local numiter=1000	//option iter specified - set to 1000 iterations
}
else if "`iterb'"==""{
	local numiter=1		//iter & iter(#) not specified - set to 1 iteration
}
else{
	//exit if not real
	if missing(real("`iterb'")){
		di as err "{opt iter()} must be a positive integer"	 
		exit 198
	}
	//exit if not an integer or negative
	if `iterb'!=int(`iterb') | `iterb'<1 {
		di as err "{opt iter()} must be a positive integer"	
		exit 198
	}
	local numiter=`iterb'	//set to number specified
}

//exit if draws in monte-carlo integration not a positive integer
if `draws'<=0{
	di as err "{opt draws()} must be an positive integer"
	exit 198
}

//exit if random draws and burns (initial elements dropped in Halton sequences) specified
if `random' & "`burn'"!=""{
	di as err "cannot combine {opt burn()} with {opt random}"
	exit 198
}

//burns options
if "`burn'"==""{
	local burn=15		//defualt value for Halton draws 
}
else{
	//exit if not real
	if missing(real("`burn'")){
		di as err "{opt burn()} must be a non-negative integer"
		exit 198
	}
	//exit if not an integer or negative
	if `burn'!=int(`burn') | `burn'<0 {
		di as err "{opt burn()} must be a non-negative integer"
		exit 198
	}
}

//create constant & include in exogenous varlist
if "`nocons'"=="" {
	tempvar cons
	gen `cons'=1
	local exog_vars `cons' `exog_vars'
}

//parse stochastic input
gettoken eqn_1 stochastic: stochastic, parse(",") 
gettoken eqn_0 stochastic: stochastic, parse(",") 
gettoken lhs rhs_1: eqn_1, parse("=")		//lhs is the first variable with a stochastic parameter
gettoken rhs_0 rhs_1: rhs_1, parse("=")		//rhs_1 is the list of demographic variables in first stochastic parameter 
local stoch_vars `stoch_vars' `lhs'			//variables with stochastic parameters
local rhs `rhs' `rhs_1'					//all demographic variables specified

//continue parsing for remaining varnames in stochastic
local c=1
while "`eqn_0'"==","{
	local c `++c'
	gettoken eqn_`c' stochastic: stochastic, parse(",")
	gettoken eqn_0 stochastic: stochastic, parse(",")
	gettoken lhs rhs_`c': eqn_`c', parse("=")
	gettoken rhs_0 rhs_`c': rhs_`c', parse("=")
	local stoch_vars `stoch_vars' `lhs'
	local rhs `rhs' `rhs_`c''
}

local rows: word count `stoch_vars'			//number of variables with stochastic parameters
local uniq_rhs: list uniq rhs				//unique demographic variables 
local cols: word count `uniq_rhs'			//number of unique demographic variables
local num_rhs_vars: word count `rhs'		//total number of demographic variables in all equations


//initial values for random parameter standard deviations
if "`initsd'"==""{					
	matrix A_INIT_SD=J(1,`rows',0.5)		//default values 0.5
}
else{
	matrix A_INIT_SD=`initsd'			//user specified
	
	//exit if insufficient supplied
	if colsof(A_INIT_SD)<`rows'{
		di as err "insufficient initial values specified in {opt initsd()} for standard deviations of stochastic coefficients"
		exit 498	
	}
	//exit if too many supplied
	else if colsof(A_INIT_SD)>`rows'{	
		di as err "too many initial values specified in {opt initsd()} for standard deviations of stochastic coefficients"
		exit 498	
	}
}

//exit if negative initial values
mata: st_numscalar("neginitsd",sum(st_matrix("A_INIT_SD"):<=0))			
if neginitsd>0{
	di as err "initial values in {opt initsd()} must be positive"
	exit 498
}

//exit if variables in stochastic are not found
foreach var of local stoch_vars {
	capture confirm variable `var'
	if _rc {
		di as err "stochastic coefficient variable " "`var'"  " not found"
		exit 111
	}
}

//exit if market varname not in dataset
capture confirm variable `markets'
if _rc {
	di as err "market identifier " "`markets'"  " not found in share data"
	exit 111
}


//parse demand elasticity input
gettoken elast_x elast_t: elast, parse(",")	
gettoken elast_0 elast_t: elast_t, parse(",")	
gettoken elast_t elast_0: elast_t, parse(",")	
gettoken elast_0 elast_n: elast_0, parse(",")	
local elast_x_vars: word count `elast_x'		//elasticity variable
local elast_t_vars: word count `elast_t'		//markets(numeric) for elasticities to be computed
local elast_n_vars: word count `elast_n'		//product string labels variable (optional)


//option elast specified
if "`elast'"!=""{

	//exit if multiple elasticity variables specified
	if `elast_x_vars'>1{
		di as err "cannot include more than one variable in option {opt elast()}"
		exit 198
	}

	//exit if multiple markets specified 
	if `elast_t_vars'>1{
		di as err "cannot include more than one market in option {opt elast()}"
		exit 198
	}

	//exit if multiple product label variables specified
	if `elast_n_vars'>1{
		di as err "cannot include more than one product name variable in option {opt elast()}"
		exit 198
	}

	//exit if elasticity variable not found
	capture confirm variable `elast_x'
	if _rc{
		di as err "variable " "`elast_x'"  " in option {opt elast()} not found"
		exit 111
	}

	//exit if product label variable not found
	capture confirm variable `elast_n'
	if _rc & "`elast_n'"!=""{
		di as err "variable " "`elast_n'"  " in option {opt elast()} not found"
		exit 111
	}

	local c=0
	foreach var of local stoch_vars {
		local c=`c'+strmatch("`elast_x'","`var'")
	}

	//exit if elasticity variable does not appear in stochastic
	if !`c'{
		di as error "variable " "`elast_x'"  " in option {opt elast()} does not appear in {opt stochastic()} "
		exit 111
	}

	//exit if elasticity market not specified 
	if "`elast_t'"==""{
		di as error "must specify a market to compute elasticities"
		exit 198
	}
	else{
		//exit if elasticity market not observed in data
		qui sum `markets'
		if `elast_t'<r(min) | `elast_t'>r(max){
			di as error "selected market for computing elasticities not found"
			exit 111
		}	
	}
}

//exit if demographic variables not included when demographic variable specified
if "`demofile'"=="" & "`initdemo'"!=""{
		di as err "initial values specified in {opt initdemo()} but no demographic file provided in {opt demofile()}"
		exit 198
}


//demographic file specified

if "`demofile'"!=""{

	//exit when observations are excluded 
	qui count if `touse'
	if r(N)!=_N{
		di as err "cannot combine if or in when using demographic data"
		exit 198
	}

	//exit if model does not include any demographic variables
	if !`cols' {
		di as err "no demographic variables entered but demographic file specified in {opt demofile()}"
		exit 198
	}

	preserve
	use "`demofile'",clear

	//exit if market varname does not match varname in demofile
	capture confirm variable `markets'
	if _rc {
		di as err "market identifier " "`markets'"  " not found in demographic data"
		exit 111
	}

	//exit if demographic variables not found
	foreach var of local uniq_rhs {
		capture confirm variable `var'
			if _rc {
			di as err "demographic variable " "`var'"  " not found in demographic data"
			exit 111
		}
	}

	//create copy of demographic data
	mata: demo_data_in("`uniq_rhs'","`markets'")	

	//exit if the number of demographic draws differs between markets
	if rowsof(uniq_D)!=1{
		di as err "number of demographic draws must be equal between markets, check data"
		exit 498
	}
	else{	
		local draws=uniq_D[1,1]
		di as txt "number of draws set to number in demographic data file"
		di as txt "draws per market is:  " "`draws'"	
	}
	restore

	//initial values for demographic variable coefficients
	if "`initdemo'"==""{
		matrix A_INIT_DEMO=J(1,`num_rhs_vars',0.5)	//default values 0.5
	}
	else{
		matrix A_INIT_DEMO=`initdemo'				//user specified

		//exit if insufficient supplied
		if colsof(A_INIT_DEMO)<`num_rhs_vars'{
			di as err "insufficient initial values specified in {opt initdemo()} for impacts of demographic-variables in stochastic coefficient equations"
			exit 498	
		}

		//exit if too many supplied
		else if colsof(A_INIT_DEMO)>`num_rhs_vars'{
			di as err "too many initial values specified in {opt initdemo()} for impacts of demographic-variables in stochastic coefficient equations"
			exit 498	
		}

	}
	
	//populate the parameter matrix PI in: beta(i)~N(PI*D(i),Sigma)
	local index=0
	matrix PI=J(`rows',`cols',.)					
	forvalues i=1/`rows' {
		forvalues j=1/`cols' {
			local ij=0
			local m: word `j' of `uniq_rhs'	
			foreach v of local rhs_`i' {					
				local ij=`ij'+strmatch("`v'","`m'")			//indicator for inclusion of demographic variable m in ith-stochastic equation
			}
			local index=`index'+`ij'					//location of parameter in initial values vector
			if `ij'{
				matrix PI[`i',`j']=A_INIT_DEMO[1,`index']		//assign parameter value
			}
		}
	}

	//display included demographic variables and initial values
	di as txt "Initial values for included demographic-variables in stochastic coefficient equations"
	matrix rownames PI = `stoch_vars'
	matrix colnames PI = `uniq_rhs'
	matrix list PI, noheader
	
	local a=0
	while !`a' {
		display as txt "Do you wish to continue?: 1=yes, 0=no?" _request(ans)		//ask for user confirmation
		if "$ans"=="0"{
			di as txt "estimation aborted"
			exit
		}
		else if "$ans"=="1" {
			di as txt "estimation continuing"
			local a=1
		}
		else {
			di as err "incorrect selection"
		} 
	}
	
}

//no demographic file specified

else{
	mata: D=.	

	//exit if demographic variables specified
	if `cols' {
		di as err "cannot include demographic variables without specifying {opt demofile()}"
		exit 198
	}
}


//parse endog input
gettoken endog inst: endog, parse("=") 	//endog contains endogenous variables
gettoken inst_0 inst: inst, parse("=")	//inst contains instruments		

//model does not contain any endogenous variables

if "`inst_0'"=="" {				
	
	local inst `endog'			//endog are all instruments
	local endog					//set to empty
	
	foreach var_1 of local exog_vars {
		foreach var_2 of local inst {	
			capture confirm variable `var_2'

			//exit if instrument is not found 
			if _rc {
				di as err "instrument " "`var_2'"  " not found"
				exit 111
			}

			//exit if instrument appears in the exogenous varlist
			if strmatch("`var_1'","`var_2'") {
				di as err "`var_2'" " is exogenous and already used as an instrument"
				exit 498
			}				
		}
	}
} 

//model contains endogenous variables

else{
	
	foreach var_1 of local exog_vars {	
		foreach var_2 of local endog {
			capture confirm variable `var_2'

			//exit if endogenous variable not found
			if _rc {
				di as err "endogenous variable " "`var_2'"  " not found"
				exit 111
			}

			//exit if endogenous variable appears in the exogenous varlist
			if strmatch("`var_1'","`var_2'") {
				di as err "`var_2'" " is included in both the exogenous and endogenous variable list"
				exit 498
			}
		}
	}	

	foreach var_1 of local exog_vars  {	
		foreach var_2 of local inst {	
			capture confirm variable `var_2'

			//exit if instrument not found
			if _rc {
				di as err "instrument " "`var_2'"  " not found"
				exit 111
			}

			//exit if instrument appears in the exogenous varlist
			if strmatch("`var_1'","`var_2'") {
				di as err "`var_2'" " is exogenous and already used as an instrument"
				exit 498
			}
		}
	}
}


local inst_all `exog_vars' `inst'	//all instruments 
local vars `exog_vars' `endog'	//all variables

//exit if more parameters than instruments
if wordcount("`vars' `rhs' `stoch_vars'")>wordcount("`inst_all'"){
		di as err "insufficient instruments"
		exit 481
}

//set to one step estimator if instruments equal parameters
if wordcount("`vars' `rhs' `stoch_vars'")==wordcount("`inst_all'") & `robustweight'{
		di as err "parameters exactly identified, switching from {opt robustweight} to one step GMM-estimator"
		local robustweight=0
}


if "`nocons'"==""{
	local vars2 cons `exog_vars2' `endog'	//model includes a constant 
}
else{
	local vars2 `exog_vars2' `endog'		//model excludes a constant
}

foreach var_1 of local stoch_vars{
	local c=0
	foreach var_2 of local vars2{
		local c=`c'+strmatch("`var_1'","`var_2'")
	}

	//exit if variable stochastic does snot appear in exogenous or endogenous varlist
	if !`c'{
		di as err "stochastic coefficient variable " "`var_1'"  " does not appear as an exogenous or endogenous variable in the model."
		di as err "To include a constant in {opt stochastic()} generate a variable named cons and ensure that {opt nocons()} is not included"
		exit 498
	}
}


//use chamberlain optimal instruments

if `use_optinst' {

	//exit if robust weighting matrix specified
	if `robustweight'{
		di as err "cannot combine {opt robustweight} with { opt optinst()}"
		exit 198
	}


	//exit if no endogenous variable specified but instruments included in optinst()
	if "`endog'"=="" & "`optinst'"!="" {
		di as err "instruments specified in {opt optinst()} but no endogenous variables included in {opt endog()} . 
            di as err "Type {opt optinst}  to use optimal instruments for model without endogenous variables"
		exit 498
	}

	//exit if more endogenous variables specified than instruments included in optinst()
	if wordcount("`endog'")>wordcount("`optinst'"){
		di as err "insufficient instruments specified in {opt optinst()} for endogneous variables"
		exit 481
	}

	foreach var_1 of local exog_vars {
		foreach var_2 of local optinst {
			capture confirm variable `var_2'

			//exit if instrument in optinst() not found 
			if _rc {
				di as err "variable " "`var_2'"  " in {opt optinst()} not found"
				exit 111
			}

			//exit if instrument in optinst() included in the exogenous varlist (already included by default)
			if strmatch("`var_1'","`var_2'") {
				di as err "`var_2'" " in {opt optinst()} is exogenous and already included in the instrument set for " "`endog'""
				exit 498
			}
		}
	}

	foreach var_1 of local endog{
		foreach var_2 of local optinst{
			if strmatch("`var_1'","`var_2'"){

				//exit if instrument in optinst() included in endogenous varlist 
				di as err "`var_2'" " in {opt optinst()} is endogenous and cannot be included"
				exit 498
			}
		}
	}	
}


//call the mata subroutine that estimates the BLP-model
mata: blp("`y_var'","`exog_vars'","`endog'","`inst'",`use_optinst',"`optinst'","`stoch_vars'","`markets'",`robustweight',`robust',`noisily',"`touse'","`elast_x'","`elast_t'","`elast_n'","`nocons'",`draws',`tolin',`tolout',`numiter',`random',`burn')

//post the estimation results
ereturn post b VAR, esample(`touse')
ereturn scalar N=r(N)				//number of observations
ereturn scalar T=r(T)				//number of markets
ereturn scalar conv=conv			//indicator for convergence
ereturn scalar draws=`draws'			//number of random draws

if !`random'{
	ereturn scalar burn=`burn'		//number of burns if halton sequence used
}

ereturn scalar obj=r(obj)			//value of GMM-objective function
ereturn local cmd "blp"				
ereturn matrix initsd=A_INIT_SD		//initial values for random parameter standard deviations

if "`demofile'"!=""{
	ereturn matrix initd=A_INIT_DEMO	//initial values for demographic variable coefficients
}

if "`elast'"!=""{
	ereturn matrix elast=elast		//elasticity matrix
}
	

//display the parameter estimates and options specified
display _newline as txt "GMM estimator of BLP-model" 

if `robustweight'{
	display _newline as txt "GMM weight matrix: robust" _col(45)"Number of obs " _col(70)"=  " as result e(N)
}
else if `use_optinst'{
		display  _newline as txt "Instruments: Chamberlain optimal" _col(45)"Number of obs " _col(70)"=  " as result e(N)
}
else{
	display _newline as txt "GMM weight matrix: unadjusted" _col(45)"Number of obs " _col(70)"=  " as result e(N)
}

display as txt _col(45)"Number of markets" _col(70)"=  " as result e(T)

if `random'{
	display as txt _col(45)"Number of random draws " _col(70)"=  " as result `draws'
}
else{
	display as txt _col(45)"Number of Halton draws " _col(70)"=  " as result `draws'

}
	
if `robust'{
	ereturn local vcetype "Robust"
}
ereturn display
}
end