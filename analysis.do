log using analysis.log, text replace

********************************************************************************
*                                                                                     
*	Filename		:	analysis.do 
*
*   Project         :   Statistical analysis for the publication "Heat stress 
*                       reveals the existence of a specialized variant of the 
*                       pachytene checkpoint in meiosis of Arabidopsis 
*                       thaliana" by Joke De Jaeger-Braet, Linda Krause, 
*                       Anika Buchholz and Arp Schnittger
*   Program author  :   Anika Buchholz                                                                    
*   Purpose         :   Prediction of median time in state (and differences
*                       therein between treatments) by interval censored 
*                       parametric survival models                                                                    
*   Date            :   05.10.2021                                                                   
*   Stata Version   :   SE 16.1                                                                    
*                                                                                                                                                        
*
*   Input data files  :	   sum_values.dta,
*                          sum_values_mutations.dta                                                           
*   Output data files :    ---                                                           
*                                                                                     
*                                                                                     
*	Required packages (non-standard ado's):	
*      - 'labelsof' by Ben Jann (distribution date 2007-04-13)
*      - 'mplotoffset' by Nicholas Winter (distribution date 2017-05-14)
*      both packages are available on SSC and may be installed via the 
*      'ssc install' command
*
*
*   Associated Programs:                                                              
*    - Include      :                       ---                                                     
*    - Run the present program before  :	---                                                                  
*    - Run the present program after   :    data_preparation.R (by Linda Krause)                                                 
*                                                                    
********************************************************************************



* Set Stata version
version 16.1



********************************************************************************



*********************************************
*****                                   *****
*****     Generate header of report     *****
*****                                   *****
*********************************************

putdocx clear
putdocx begin

putdocx paragraph, style("Title") halign(center)
putdocx text ("Statistical analysis")
putdocx paragraph, style("Subtitle") halign(center)
putdocx text ("for"), linebreak
putdocx text ("Heat stress reveals the existence of a specialized variant of the pachytene checkpoint in meiosis of Arabidopsis thaliana"), font(italic) linebreak
putdocx text (" by Joke De Jaeger-Braet, Linda Krause, Anika Buchholz and Arp Schnittger")



**********************************************
***   Description of statistical methods   ***
**********************************************

putdocx paragraph, style("Heading1")
putdocx text ("Statistical methods")
putdocx paragraph
putdocx text ("Because modelling the issues of clustered data, left, right and/or interval censoring as well as multistate data in one combined model is not possible, we decided to reduce the complexity of the analysis. Since the lengths of the different states do not depend on each other, we will build one separate model for each state. This simultaneously allows us to simplify the mixture of left, right and/or interval censored data w.r.t. duration of each state to interval (and right) censoring. Thus, we apply parametric models for interval-censored survival time data with a clustered sandwich estimator of variance to address the clustering of meiocytes within anther-sacs, including effects of the heat treatment, genotype and their interactions. The underlying distribution of each parametric model will be chosen based on the Akaike Information Criterion (AIC) with exponential, Gompertz, log-logistic, Weibull and log-normal distribution as candidates."), linebreak
putdocx text ("The state-specific models will be based on information from all cells which were captured through the microscope for at least one time point in that specific state. The event of interest is the transition of a cell from one state to the next. Each cell for which the exact beginning and end of the state is known is modelled as having an event, with the event time calculated as the difference between the start of the next state and the end of the previous state. Cells where the exact time points of either the transition from the previous state to the state of interest or to the next state are not known are modelled as interval-censored data points with the lower limit of the interval being the time where the cell was observed in this specific state and the upper limit of the interval being one time unit after / before when the cell was observed in the previous / next state, respectively. If for a cell a certain state either at the beginning or the end was not observed, the cell was modelled as right-censored with the censoring time being the minimum observed time for this cell in the state of interest."), linebreak 
putdocx text ("In addition to the individual states we calculated a model for the duration from MT array state 2 to 13 in an analogous fashion."), linebreak 
putdocx text ("Since for wildtype cells more growth conditions (21°C, HS30°C, HS34°C, LT30°C ) have been considered than in mutated cells (21°C and HS34°C), the interaction is restricted to those combinations available. In addition, some combinations had to be excluded from specific models because no (or hardly any) events had been observed."), linebreak 
putdocx text ("Estimation results are presented as predicted marginal median times (or corresponding contrasts), together with 95% confidence intervals. Since the analysis is of exploratory nature, no adjustment for multiplicity was applied."), linebreak   





********************************************************************************
********************************************************************************
********************************************************************************



****************************************************************************
*****                                                                  *****
*****     Prepare the data and                                         *****
*****     combine wildtype and mutations into one data set             *****
*****                                                                  *****
****************************************************************************



*************************************
***   Prepare mutation data set   ***
*************************************


* Read data
use sum_values_mutations.dta, clear


* Extract genotype (mutation) and actual treatment, both stored in "treatment" 
* indicator string
gen tmp = subinstr(treatment, "_", " ", 1)

gen mutation = word(tmp, 1)
rename treatment treatment_original

gen treatment = word(tmp, 2)
	* CHECK
	tab2 treatment_original treatment mutation, first
codebook treatment
replace treatment = "HS34_A" if treatment=="HS34"

drop treatment_original
drop tmp
save tmp, replace



*************************************
***   Prepare wildtype data set   ***
*************************************

* Read data
use sum_values.dta, clear


* Anther and cell ID are not unique, but are repeated over the
* different therapies -> generate a unique cell identifier
gen mutation = "Wildtype"
codebook treatment



*********************************
***   Append both data sets   ***
*********************************

append using tmp
erase tmp.dta



**************************************
***   Recode and label variables   ***
**************************************

* Generate / encode ID, genotype and treatment variables
codebook treatment
egen tmp = concat(mutation treatment ID), punct("_")
encode tmp, gen(id)
drop tmp ID
egen tmp = concat(treatment anther), punct("_")
encode tmp, gen(anther_id)
drop tmp

encode treatment, gen(treat) label(l_treat)
label values treat l_treat
	* CHECK: 
	codebook treat
lab define lab_treat 1 "HS30°C" 2 "HS30°C late" 3 "HS34°C" 4 "HS34°C late" 5 "21°C" 6 "LT30°C"
lab val treat lab_treat
	* CHECK:
	codebook treat
drop treatment
	
encode mutation, gen(mut) label(l_mut)
label values mut l_mut
	* CHECK:
	codebook mut
label def lab_mut 1 "atm" 2 "dmc1" 3 "msh4" 4 "spo11" 5 "WT"
lab val mut lab_mut
drop anther mutation
	* CHECK:
	codebook mut
	

* Drop variables that are not required for the analyses
drop P_time* PT_time* HT_time* FT_time* MM2_time*

	
* Introduce appropriate variable labels
label var id "Cell identifier"
label var anther_id "Anther identifier"
label var treat "Treatment (heat condition)"
label var mut "Genotype"

/*
Remark:
	State names in the data set and their 'labels' used in the manuscript
	H:  MT array states 2-3-4 / late leptotene - early pachytene
	F:  MT array states 5-6   / pachytene - diakinesis
	M:  MT array states 7-8-9 / metaphase I - anaphase ID
	I:  MT array states 10-11 / telophase I - interkinesis
	M2: MT array states 12-13 / metaphase II - anaphase II
	T:  MT array state 14     / telophase II
*/	
label var H_time_l	"Time in MT array states 2-3-4 (lower limit)"
label var H_time_u	"Time in MT array states 2-3-4 (upper limit)"	
label var F_time_l	"Time in MT array states 5-6 (lower limit)"
label var F_time_u	"Time in MT array states 5-6 (upper limit)"	
label var M_time_l	"Time in MT array states 7-8-9 (lower limit)"
label var M_time_u	"Time in MT array states 7-8-9 (upper limit)"	
label var I_time_l	"Time in MT array states 10-11 (lower limit)"
label var I_time_u	"Time in MT array states 10-11 (upper limit)"	
label var M2_time_l	"Time in MT array states 12-13 (lower limit)"
label var M2_time_u	"Time in MT array states 12-13 (upper limit)"	
label var T_time_l	"Time in MT array state 14 (lower limit)"
label var T_time_u	"Time in MT array state 14 (upper limit)"	
label var HM2_time_l	"Time in MT array states 2-13 (lower limit)"
label var HM2_time_u	"Time in MT array states 2-13 (upper limit)"	
/* 
Remark:
   The time variables for interval-censored time-to-event analysis are coded as
   described in the Stata help for 'stintreg':
   The interval time variables time_l and time_u have the following form

             Type of data                        t_l    t_u
             --------------------------------------------------
             uncensored data         a = [a,a]    a      a 
             interval-censored data      (a,b]    a      b
             left-censored data          (0,b]    .      b
             left-censored data          (0,b]    0      b
             right-censored data      [a,+inf)    a      . 
             missing                              .      .
             missing                              0      .
             --------------------------------------------------
*/


* Define local macros with state labels for the report document
local H_name "MT array states 2-3-4"
local F_name "MT array states 5-6"
local M_name "MT array states 7-8-9"
local I_name "MT array states 10-11"
local M2_name "MT array states 12-13"
local T_name "MT array state 14"
local HM2_name "MT array states 2-13"
	
	
	
	
	
********************************************************************************
********************************************************************************
********************************************************************************	





*******************************************************
*****                                             *****
*****     Provide general information on data     *****
*****                                             *****
*******************************************************


******************************************************************
***   No. of cells that are under observation "in principle"   ***
******************************************************************

putdocx paragraph, style("Heading1")
putdocx text ("General information")
putdocx paragraph, style("Heading2")
putdocx text ("No. of cells under observation 'in principle'")

tab mut treat, matcell(tab1_cell)
putdocx table tbl = matrix(tab1_cell)

putdocx table tbl(1,.), addrows(2,before)
putdocx table tbl(1,1) = ("Heat condition"), colspan(6) halign(center)
qui levelsof treat
forvalues i = 1/6 {
	local lab: label (treat) `i'
	putdocx table tbl(2,`i') = ("`lab'")
}

putdocx table tbl(.,1), addcols(1,before)
qui levelsof mut
forvalues i = 1/5 {
	local row = `i'+2
	local lab: label (mut) `i'
	putdocx table tbl(`row',1) = ("`lab'")
}

putdocx table tbl(1,1) = (""), rowspan(2)



**************************************************************************
***   No. of cells and no. of events that are available for analysis   ***
***   in the different states                                          ***
**************************************************************************

putdocx paragraph, style("Heading2")
putdocx text ("No. of cells and no. of events that are available for analysis in the different states")
putdocx paragraph
putdocx text ("In survival analysis, the information actually originates from the events rather than the sample size. The former is also called 'effective sample size'. Here, an 'event' means that the time in state (duration) has been observed exactly (both transition into the current and next state have been observed). In case the transition into the current and/or next state has not been observed exactly, the time in state is 'interval-censored'. Yet, in this case we can incorporate the information that the 'time in state' lies within a specified time span if we have observations in the previous and next state. If we do not have information on both the previous and next state, we only know that the minimum time in the current state, i.e. this is defined as right censoring here.")


* Generate event indicator variable
lab def lab_ev_c 0 "only minimum known" 1 "event known exactly" 2 "time span known"
foreach state in "H" "F" "M" "I" "M2" "T" "HM2" {
	gen ev_type_`state' = cond(missing(`state'_time_l) & missing(`state'_time_l), ., cond(!missing(`state'_time_l) & !missing(`state'_time_u) & `state'_time_l==`state'_time_u, 1, cond(!missing(`state'_time_l) & !missing(`state'_time_u) & `state'_time_l!=`state'_time_u, 2, cond(!missing(`state'_time_l) & missing(`state'_time_u), 0, 99))))
	lab val ev_type_`state' lab_ev_c
	lab var ev_type_`state' "Event status for ``state'_name'"
	bysort mut: tab treat ev_type_`state' 
	* This version will be used throughout the code 
}


* Use a differently coded version for the report (censoring shall by default be 
* given in the third column, not in the first)
lab def lab_ev 1 "event known exactly" 2 "time span known" 3 "only minimum known"
foreach state in "H" "F" "M" "I" "M2" "T" "HM2" {
	
	gen ev_type_`state'_out = cond(ev_type_`state'==0, 3, ev_type_`state')
	lab val ev_type_`state'_out lab_ev
	lab var ev_type_`state'_out "Event status for ``state'_name' (recoded ev_type_`state' - for reporting only)"
	
	putdocx paragraph, style("Heading3")
	putdocx text ("``state'_name'")
	
	* Loop over all mutations
	forvalues i = 1/5 {
				
		local lab: label (mut) `i'
		putdocx paragraph
		putdocx text ("`lab'"), bold
		
		tab treat ev_type_`state'_out if mut==`i', matcell(tab_cell) matrow(tab_row) 
		putdocx table tbl = matrix(tab_cell)
		local nrow = rowsof(tab_cell)
		local ncol = colsof(tab_cell)
		* Extract available column categories:
		tab ev_type_`state'_out if mut==`i', matrow(tab_col)

		* Not all event types are present in all combinations 
		* -> use dynamic range and labelling
		putdocx table tbl(1,.), addrows(2,before)
		forvalues j = 1/`ncol' {
			local ev`j' = tab_col[`j',1]
		}
		forvalues j = 1/`ncol' {
			local lab: label (ev_type_`state'_out) `j'
			putdocx table tbl(2,`j') = ("`lab'")
		}

		putdocx table tbl(.,1), addcols(2,before)
		
		* Not all treatments are present for all mutations -> Identify those 
		* present and use lables accordingly
		forvalues j = 1/`nrow' {
			local val`j' = tab_row[`j',1]
		}
		forvalues j = 1/`nrow' {
			local row = `j'+2
			local lab: label (treat) `val`j''
			putdocx table tbl(`row',2) = ("`lab'")
		}

		
		* Add 'total' column
		local col = `ncol'+2
		putdocx table tbl(.,`col'), addcols(1,after) 
			* Here we have to use index column '3', since the merged first row 
			* counts as 1 column
		local col = `ncol'+3
		tab treat if mut==`i' & !missing(ev_type_`state'), matcell(tab_total) 
		forvalues j = 1/`nrow' {
			local row = `j'+2
			putdocx table tbl(`row',`col') = (tab_total[`j', 1])
		}
		
		* Merge header columns/rows
		putdocx table tbl(1,`col') = ("Total"), rowspan(2)
		putdocx table tbl(1,3) = ("Duration in ``state'_name'"), colspan(`ncol') halign(center)
		putdocx table tbl(3,1) = ("Treatment"), rowspan(`nrow') valign(center)
		putdocx table tbl(1,1) = (""), span(2,2)	
	}
}





********************************************************************************
********************************************************************************
********************************************************************************





***********************************************************************
*****                                                             *****
*****     Initial inspection of potential problematic aspects     *****
*****                                                             *****
***********************************************************************



/* 
Since in analyses of previous versions of the data, the models with generalized 
gamma distribution caused problems, it is now checked in advance whether it is 
sensible to include this distribution as a candidate for model fitting. 
*/



********************************************
* Checks on generalized gamma distribution *
********************************************

foreach state in "H" "F" "M" "I" "M2" "T" "HM2" {
	di ""
	di "**********     State `state'     **********"
	di ""
	if ("`state'"=="H") local ifexclude = "if !(inlist(treat, 2, 4) & mut==5)"
	else local ifexclude ""
	cap noi stintreg b5.treat##b5.mut `ifexclude', interval(`state'_time_l `state'_time_u) vce(cluster anther_id) dist(ggamma) base difficult 
		* Remark: 
		* Control treatment is coded as '5' -> use 'b5.treat'
		* Wildtype cells are coded as '5' too -> use 'b5.mut' 
	di ""
	di "**********     END (State `state')     **********"
	di ""
	pause
}

/*
There are convergence problems in many states which also cannot be resolved /
improved by using the 'difficult' option (using a different stepping 
algorithm in nonconcave regions). 

Consequence: Do not use ggamma as candidate distribution, since the estimates
             will most probably not be reliable.

Remaining candidate distributions for model building:
exponential, Gompertz, log-logistic, Weibull and lognormal
*/





********************************************************************************
********************************************************************************
********************************************************************************





******************************************************************************
*****                                                                    *****
*****     Analysis of MT array states with prediction of median time     *****   
*****                                                                    *****
******************************************************************************



**********************************************************************
***                                                                ***
***                      Investigate state H                       ***
***   (MT array states 2-3-4 / late leptotene - early pachytene)   ***
***                                                                ***
**********************************************************************

local state "H"

putdocx pagebreak
putdocx paragraph, style("Heading1")
putdocx text ("``state'_name'")



************************************************************************
* Check whether there is sufficient information on all combinations of *
* genotype and heat condition                                          *
************************************************************************

* No. of events per combination of genotype and heat condition
bysort mut: tab treat ev_type_`state'
* -> For wildtype cells, the treatments HS30_B and HS34_B have hardly any
*    events. Therefore they are excluded from the analysis. 

codebook treat mut
local ifexclude = "if !(inlist(treat, 2, 4) & mut==5)"
putdocx paragraph
putdocx text ("Remarks:"), bold linebreak
putdocx text ("- 'HS30°C late' and 'HS34°C late' in wildtype cells have to be excluded from the analysis of ``state'_name', because there are hardly any completely ('event known exactly') or closely ('time span known') observed durations for these two conditions."), linebreak




********************************************
* Select best distribution in terms of AIC *
********************************************
	
* Prepare matrix for collecting AIC values
mat aicmat_`state' = J(1,5,.)
mat colnames aicmat_`state' = exponential gompertz loglogistic weibull lognormal
	
* Identify distributions that cannot be estimated and write them to report
foreach dist in "exponential" "gompertz" "loglogistic" "weibull" "lognormal" {
	di ""
	di "**********     State `state', `dist' distribution     **********"
	di ""
	cap noi stintreg b5.treat##b5.mut `ifexclude', interval(`state'_time_l `state'_time_u) vce(cluster anther_id) dist(`dist') base difficult 
		* Remark: 
		* Control treatment is coded as '5' -> use 'b5.treat'
		* Wildtype cells are coded as '5' too -> use 'b5.mut' 
	cap noi estat ic
	mat ic = r(S)
	local aic = ic[1,colnumb(ic, "AIC")]
	mat aicmat_`state'[1,colnumb(aicmat_`state', "`dist'")] = `aic'
	local aic .
	di ""
	di "**********     END (State `state', `dist' distribution)     **********"
	di ""
	pause
}
mat list aicmat_`state'


* Identify distribution with best AIC and write it to report
qui mat list aicmat_`state'
local mcn : coln(aicmat_`state')

mata mata clear
mata x = st_matrix("aicmat_`state'")
mata xs = sort(x',1)
mata v = xs[(1), .]
mata min_aic_col = selectindex(x:==rowmin(v))
mata st_numscalar("min_aic", min_aic_col)
	
local mincol = scalar(min_aic) 
local mincolname: word `mincol' of `mcn'
local best_`state' = "`mincolname'"

local best_`state'_name = "`best_`state''"
if ("`best_`state'_name'"=="weibull") local best_`state'_name = "Weibull"

putdocx text ("- The distribution yielding the best AIC is `best_`state'_name'.")
		
		
* Clean up			
scalar drop _all   



	
***************
* Final model *
***************

* Specify best distribution chosen above
local dist = "`best_`state''"

* Mixed effects model
cap noi stintreg b5.treat##b5.mut `ifexclude', interval(`state'_time_l `state'_time_u) vce(cluster anther_id) dist(`dist') noempty difficult  
	* Remarks: 
	* - Control treatment is coded as '5' -> use 'b5.treat'
	* - Wildtype cells are coded as '5' too -> use 'b5.mut' 
est sto `state'_model	

putdocx paragraph, style("Heading2")
local d = "`dist'"
if ("`dist'"=="weibull") local d = "Weibull"
putdocx text ("Interval censored parametric survival model with `d' distribution")
putdocx table tbl = etable	
putdocx table tbl(.,2/4), nformat(%6.2f)
putdocx table tbl(.,5), nformat(%5.3f)
putdocx table tbl(.,6/7), nformat(%6.2f)


* Median survival time (i.e. time in state)
putdocx pagebreak
putdocx paragraph, style("Heading2")
putdocx text ("Predicted marginal median time in state `state' based on interval censored `d' survival model")

putdocx paragraph, style("Heading3")
putdocx text ("Predicted marginal median times by genotype and heat condition")
margins mut#treat, predict(median) noempty post coeflegend
est sto `state'_margins

est resto `state'_model	
margins mut#treat, predict(median) noempty
cap putdocx table tbl = etable

mplotoffset, recast(scatter) graphregion(color(white) margin(r+5)) plotregion(margin(l+5 r+5)) xtitle("Genotype") ytitle("Predicted median time") title("``state'_name'") legend(col(4)) offset(0.1)
graph export MarginalMedians_`state'.png, replace
putdocx paragraph
putdocx image MarginalMedians_`state'.png, width(15cm)
putdocx table tbl(.,2/4), nformat(%6.2f)
putdocx table tbl(.,5), nformat(%5.3f)
putdocx table tbl(.,6/7), nformat(%6.2f)
			
			
putdocx pagebreak
putdocx paragraph, style("Heading3")
putdocx text ("Pairwise differences in predicted marginal medians")
putdocx paragraph

est resto `state'_model
margins, predict(median) at(treat=(1 3 5 6) mut=(1/5)) pwcompare noempty
	* This works fine, but does not give self-explaining labels.
* -> Calculate margins manually via lincom to have nicer labels.
est resto `state'_margins


* Pairwise comparisons of all 6 treatments within wildtype cells
local tmp_mut = 5
local row_counter = 1
putdocx table tbl = (1, 4)
putdocx table tbl(1,1) = ("Contrast") 
putdocx table tbl(1,2) = ("Estimated difference in median survival times") 
putdocx table tbl(1,3) = ("95% Confidence Interval"), colspan(2)
forvalues tmp_treat1 = 1/6 {
	forvalues tmp_treat2 = `tmp_treat1'/6 {
		cap lincom (_b[`tmp_treat1'.treat#`tmp_mut'.mut]-_b[`tmp_treat2'.treat#`tmp_mut'.mut])
		if (_rc==0 & !(missing(r(estimate)) | r(estimate)==0)) {
			putdocx table tbl(`row_counter',.), addrows(1,after)
			local row_counter = `row_counter'+1
			local lab_m: label (mut) `tmp_mut'
			local lab_t1: label (treat) `tmp_treat1'
			local lab_t2: label (treat) `tmp_treat2'
			putdocx table tbl(`row_counter',1) = ("`lab_m'#`lab_t1' - `lab_m'#`lab_t2'")						
			putdocx table tbl(`row_counter',2) = (`r(estimate)')
			putdocx table tbl(`row_counter',3) = (`r(lb)')
			putdocx table tbl(`row_counter',4) = (`r(ub)')
		}
		putdocx table tbl(.,2/4), nformat(%6.2f)
	}
}


* Pairwise comparison of the 2 available treatments within mutated cells
forvalues tmp_mut = 1/4 {
	local tmp_treat1 3
	local tmp_treat2 5
	cap lincom (_b[`tmp_treat1'.treat#`tmp_mut'.mut]-_b[`tmp_treat2'.treat#`tmp_mut'.mut])
	if (_rc==0 & !(missing(r(estimate)) | r(estimate)==0)) {
		local row_counter = 1
		putdocx table tbl = (1, 4)
		putdocx table tbl(1,1) = ("Contrast") 
		putdocx table tbl(1,2) = ("Estimated difference in median survival times") 
		putdocx table tbl(1,3) = ("95% Confidence Interval"), colspan(2)
		putdocx table tbl(`row_counter',.), addrows(1,after)
		local row_counter = `row_counter'+1
		local lab_m: label (mut) `tmp_mut'
		local lab_t1: label (treat) `tmp_treat1'
		local lab_t2: label (treat) `tmp_treat2'
		putdocx table tbl(`row_counter',1) = ("`lab_m'#`lab_t1' - `lab_m'#`lab_t2'")						
		putdocx table tbl(`row_counter',2) = (`r(estimate)')
		putdocx table tbl(`row_counter',3) = (`r(lb)')
		putdocx table tbl(`row_counter',4) = (`r(ub)')
	}
	putdocx table tbl(.,2/4), nformat(%6.2f)
}

* Pairwise comparisons of genotypes for the 2 treatments available in all cell lines
foreach tmp_treat in 3 5 {
	local row_counter = 1
	putdocx table tbl = (1, 4)
	putdocx table tbl(1,1) = ("Contrast") 
	putdocx table tbl(1,2) = ("Estimated difference in median survival times") 
	putdocx table tbl(1,3) = ("95% Confidence Interval"), colspan(2)
	forvalues tmp_mut1 = 1/5 {
		forvalues tmp_mut2 = `tmp_mut1'/5 {
			cap lincom (_b[`tmp_treat'.treat#`tmp_mut1'.mut]-_b[`tmp_treat'.treat#`tmp_mut2'.mut])
			if (_rc==0 & !(missing(r(estimate)) | r(estimate)==0)) {
				putdocx table tbl(`row_counter',.), addrows(1,after)
				local row_counter = `row_counter'+1
				local lab_m1: label (mut) `tmp_mut1'
				local lab_m2: label (mut) `tmp_mut2'
				local lab_t: label (treat) `tmp_treat'
				putdocx table tbl(`row_counter',1) = ("`lab_m1'#`lab_t' - `lab_m2'#`lab_t'")						
				putdocx table tbl(`row_counter',2) = (`r(estimate)')
				putdocx table tbl(`row_counter',3) = (`r(lb)')
				putdocx table tbl(`row_counter',4) = (`r(ub)')
			}
		}
	}
	putdocx table tbl(.,2/4), nformat(%6.2f)
}


putdocx paragraph, style("Heading3")
putdocx text ("Comparison of HS34°C vs. 21°C in spo11 with HS34°C vs. 21°C in WT")
est resto `state'_margins
lincom (_b[3.treat#5.mut]-_b[5.treat#5.mut])-(_b[3.treat#4.mut]-_b[5.treat#4.mut])

putdocx table tbl = (2, 4)
putdocx table tbl(1,1) = ("Contrast") 
putdocx table tbl(1,2) = ("Estimated difference between differences in median surv. times") 
putdocx table tbl(1,3) = ("95% Confidence Interval"), colspan(2)
local lab_m5: label (mut) 5
local lab_m4: label (mut) 4
local lab_t3: label (treat) 3
local lab_t5: label (treat) 5
putdocx table tbl(2,1) = ("(`lab_m5'#`lab_t3' - `lab_m5'#`lab_t5') - (`lab_m4'#`lab_t3' - `lab_m4'#`lab_t5')")						
putdocx table tbl(2,2) = (`r(estimate)')
putdocx table tbl(2,3) = (`r(lb)')
putdocx table tbl(2,4) = (`r(ub)')
putdocx table tbl(.,2/4), nformat(%6.2f)





********************************************************************************
********************************************************************************
********************************************************************************





**********************************************************
***                                                    ***
***                Investigate state F                 ***
***   (MT array states 5-6 / pachytene - diakinesis)   ***
***                                                    ***
**********************************************************

local state "F"

putdocx pagebreak
putdocx paragraph, style("Heading1")
putdocx text ("``state'_name'")



************************************************************************
* Check whether there is sufficient information on all combinations of *
* genotype and heat condition                                          *
************************************************************************

* No. of events per combination of genotype and heat condition
bysort mut: tab treat ev_type_`state'
* -> Available information seems to be sufficient.

local ifexclude = ""
putdocx paragraph
putdocx text ("Remark:"), bold linebreak



********************************************
* Select best distribution in terms of AIC *
********************************************
	
* Prepare matrix for collecting AIC values
mat aicmat_`state' = J(1,5,.)
mat colnames aicmat_`state' = exponential gompertz loglogistic weibull lognormal
	
* Identify distributions that could not be estimated and write them to report
foreach dist in "exponential" "gompertz" "loglogistic" "weibull" "lognormal" {
	di ""
	di "**********     State `state', `dist' distribution     **********"
	di ""
	cap noi stintreg b5.treat##b5.mut `ifexclude', interval(`state'_time_l `state'_time_u) vce(cluster anther_id) dist(`dist') base difficult 
		* Remark: 
		* Control treatment is coded as '5' -> use 'b5.treat'
		* Wildtype cells are coded as '5' too -> use 'b5.mut' 
	cap noi estat ic
	mat ic = r(S)
	local aic = ic[1,colnumb(ic, "AIC")]
	mat aicmat_`state'[1,colnumb(aicmat_`state', "`dist'")] = `aic'
	local aic .
	di ""
	di "**********     END (State `state', `dist' distribution)     **********"
	di ""
	pause
}
mat list aicmat_`state'


* Identify distribution with best AIC and write it to report
qui mat list aicmat_`state'
local mcn : coln(aicmat_`state')

mata mata clear
mata x = st_matrix("aicmat_`state'")
mata xs = sort(x',1)
mata v = xs[(1), .]
mata min_aic_col = selectindex(x:==rowmin(v))
mata st_numscalar("min_aic", min_aic_col)
	
local mincol = scalar(min_aic) 
local mincolname: word `mincol' of `mcn'
local best_`state' = "`mincolname'"

local best_`state'_name = "`best_`state''"
if ("`best_`state'_name'"=="weibull") local best_`state'_name = "Weibull"

putdocx text ("- The distribution yielding the best AIC is `best_`state'_name'.")
			
* Clean up			
scalar drop _all   



	
***************
* Final model *
***************

* Specify best distribution chosen above
local dist = "`best_`state''"

* Mixed effects model
cap noi stintreg b5.treat##b5.mut `ifexclude', interval(`state'_time_l `state'_time_u) vce(cluster anther_id) dist(`dist') noempty difficult  
	* Remarks: 
	* - Control treatment is coded as '5' -> use 'b5.treat'
	*   Wildtype cells are coded as '5' too -> use 'b5.mut' 
est sto `state'_model	
putdocx paragraph, style("Heading2")
local d = "`dist'"
if ("`dist'"=="weibull") local d = "Weibull"
putdocx text ("Interval censored parametric survival model with `d' distribution")
putdocx table tbl = etable	
putdocx table tbl(.,2/4), nformat(%6.2f)
putdocx table tbl(.,5), nformat(%5.3f)
putdocx table tbl(.,6/7), nformat(%6.2f)


* Median survival time (i.e. time in state)
putdocx pagebreak
putdocx paragraph, style("Heading2")
putdocx text ("Predicted marginal median time in state `state' based on interval censored `d' survival model")

putdocx paragraph, style("Heading3")
putdocx text ("Predicted marginal median times by genotype and heat condition")
margins mut#treat, predict(median) noempty post coeflegend
est sto `state'_margins

est resto `state'_model	
margins mut#treat, predict(median) noempty
cap putdocx table tbl = etable

mplotoffset, recast(scatter) graphregion(color(white) margin(r+5)) plotregion(margin(l+5 r+5)) xtitle("Genotype") ytitle("Predicted median time") title("``state'_name'") legend(col(4)) offset(0.1)
graph export MarginalMedians_`state'.png, replace
putdocx paragraph
putdocx image MarginalMedians_`state'.png, width(15cm)
putdocx table tbl(.,2/4), nformat(%6.2f)
putdocx table tbl(.,5), nformat(%5.3f)
putdocx table tbl(.,6/7), nformat(%6.2f)
			

putdocx pagebreak
putdocx paragraph, style("Heading3")
putdocx text ("Pairwise differences in predicted marginal medians")
putdocx paragraph

est resto `state'_model
margins, predict(median) at(treat=(1 3 5 6) mut=(1/5)) pwcompare noempty
	* This works fine, but does not give self-explaining labels.
* -> Calculate margins manually via lincom to have nicer labels.
est resto `state'_margins

* Pairwise comparisons of all 6 treatments within wildtype cells
local tmp_mut = 5
local row_counter = 1
putdocx table tbl = (1, 4)
putdocx table tbl(1,1) = ("Contrast") 
putdocx table tbl(1,2) = ("Estimated difference in median survival times") 
putdocx table tbl(1,3) = ("95% Confidence Interval"), colspan(2)
forvalues tmp_treat1 = 1/6 {
	forvalues tmp_treat2 = `tmp_treat1'/6 {
		cap lincom (_b[`tmp_treat1'.treat#`tmp_mut'.mut]-_b[`tmp_treat2'.treat#`tmp_mut'.mut])
		if (_rc==0 & !(missing(r(estimate)) | r(estimate)==0)) {
			putdocx table tbl(`row_counter',.), addrows(1,after)
			local row_counter = `row_counter'+1
			local lab_m: label (mut) `tmp_mut'
			local lab_t1: label (treat) `tmp_treat1'
			local lab_t2: label (treat) `tmp_treat2'
			putdocx table tbl(`row_counter',1) = ("`lab_m'#`lab_t1' - `lab_m'#`lab_t2'")						
			putdocx table tbl(`row_counter',2) = (`r(estimate)')
			putdocx table tbl(`row_counter',3) = (`r(lb)')
			putdocx table tbl(`row_counter',4) = (`r(ub)')
		}
		putdocx table tbl(.,2/4), nformat(%6.2f)
	}
}


* Pairwise comparison of the 2 available treatments within mutated cells
forvalues tmp_mut = 1/4 {
	local tmp_treat1 3
	local tmp_treat2 5
	cap lincom (_b[`tmp_treat1'.treat#`tmp_mut'.mut]-_b[`tmp_treat2'.treat#`tmp_mut'.mut])
	if (_rc==0 & !(missing(r(estimate)) | r(estimate)==0)) {
		local row_counter = 1
		putdocx table tbl = (1, 4)
		putdocx table tbl(1,1) = ("Contrast") 
		putdocx table tbl(1,2) = ("Estimated difference in median survival times") 
		putdocx table tbl(1,3) = ("95% Confidence Interval"), colspan(2)
		putdocx table tbl(`row_counter',.), addrows(1,after)
		local row_counter = `row_counter'+1
		local lab_m: label (mut) `tmp_mut'
		local lab_t1: label (treat) `tmp_treat1'
		local lab_t2: label (treat) `tmp_treat2'
		putdocx table tbl(`row_counter',1) = ("`lab_m'#`lab_t1' - `lab_m'#`lab_t2'")						
		putdocx table tbl(`row_counter',2) = (`r(estimate)')
		putdocx table tbl(`row_counter',3) = (`r(lb)')
		putdocx table tbl(`row_counter',4) = (`r(ub)')
	}
	putdocx table tbl(.,2/4), nformat(%6.2f)
}


* Pairwise comparisons of genotypes for the 2 treatments available in all cell lines
foreach tmp_treat in 3 5 {
	local row_counter = 1
	putdocx table tbl = (1, 4)
	putdocx table tbl(1,1) = ("Contrast") 
	putdocx table tbl(1,2) = ("Estimated difference in median survival times") 
	putdocx table tbl(1,3) = ("95% Confidence Interval"), colspan(2)
	forvalues tmp_mut1 = 1/5 {
		forvalues tmp_mut2 = `tmp_mut1'/5 {
			cap lincom (_b[`tmp_treat'.treat#`tmp_mut1'.mut]-_b[`tmp_treat'.treat#`tmp_mut2'.mut])
			if (_rc==0 & !(missing(r(estimate)) | r(estimate)==0)) {
				putdocx table tbl(`row_counter',.), addrows(1,after)
				local row_counter = `row_counter'+1
				local lab_m1: label (mut) `tmp_mut1'
				local lab_m2: label (mut) `tmp_mut2'
				local lab_t: label (treat) `tmp_treat'
				putdocx table tbl(`row_counter',1) = ("`lab_m1'#`lab_t' - `lab_m2'#`lab_t'")						
				putdocx table tbl(`row_counter',2) = (`r(estimate)')
				putdocx table tbl(`row_counter',3) = (`r(lb)')
				putdocx table tbl(`row_counter',4) = (`r(ub)')
			}
		}
	}
	putdocx table tbl(.,2/4), nformat(%6.2f)
}







********************************************************************************
********************************************************************************
********************************************************************************





***************************************************************
***                                                         ***
***                     Investigate state M                 ***
***   (MT array states 7-8-9 / metaphase I - anaphase ID)   ***
***                                                         ***
***************************************************************

local state "M"

putdocx pagebreak
putdocx paragraph, style("Heading1")
putdocx text ("``state'_name'")



************************************************************************
* Check whether there is sufficient information on all combinations of *
* genotype and heat condition                                          *
************************************************************************

* No. of events per combination of genotype and heat condition
bysort mut: tab treat ev_type_`state'
* -> Available information seems to be sufficient.

local ifexclude = ""
putdocx paragraph
putdocx text ("Remark:"), bold linebreak



********************************************
* Select best distribution in terms of AIC *
********************************************

* Prepare matrix for collecting AIC values
mat aicmat_`state' = J(1,5,.)
mat colnames aicmat_`state' = exponential gompertz loglogistic weibull lognormal
	
* Identify distributions that could not be estimated and write them to report
foreach dist in "exponential" "gompertz" "loglogistic" "weibull" "lognormal" {
	di ""
	di "**********     State `state', `dist' distribution     **********"
	di ""
	cap noi stintreg b5.treat##b5.mut `ifexclude', interval(`state'_time_l `state'_time_u) vce(cluster anther_id) dist(`dist') base difficult 
		* Remark: 
		* Control treatment is coded as '5' -> use 'b5.treat'
		* Wildtype cells are coded as '5' too -> use 'b5.mut' 
	cap noi estat ic
	mat ic = r(S)
	local aic = ic[1,colnumb(ic, "AIC")]
	mat aicmat_`state'[1,colnumb(aicmat_`state', "`dist'")] = `aic'
	local aic .
	di ""
	di "**********     END (State `state', `dist' distribution)     **********"
	di ""
	pause
}
mat list aicmat_`state'


* Identify distribution with best AIC and write it to report
qui mat list aicmat_`state'
local mcn : coln(aicmat_`state')

mata mata clear
mata x = st_matrix("aicmat_`state'")
mata xs = sort(x',1)
mata v = xs[(1), .]
mata min_aic_col = selectindex(x:==rowmin(v))
mata st_numscalar("min_aic", min_aic_col)
	
local mincol = scalar(min_aic) 
local mincolname: word `mincol' of `mcn'
local best_`state' = "`mincolname'"

local best_`state'_name = "`best_`state''"
if ("`best_`state'_name'"=="weibull") local best_`state'_name = "Weibull"

putdocx text ("- The distribution yielding the best AIC is `best_`state'_name'.")
			
* Clean up			
scalar drop _all   



	
***************
* Final model *
***************

* Specify best distribution chosen above
local dist = "`best_`state''"

* Mixed effects model
cap noi stintreg b5.treat##b5.mut `ifexclude', interval(`state'_time_l `state'_time_u) vce(cluster anther_id) dist(`dist') noempty difficult  
	* Remarks: 
	* - Control treatment is coded as '5' -> use 'b5.treat'
	*   Wildtype cells are coded as '5' too -> use 'b5.mut' 
est sto `state'_model
putdocx paragraph, style("Heading2")
local d = "`dist'"
if ("`dist'"=="weibull") local d = "Weibull"
putdocx text ("Interval censored parametric survival model with `d' distribution")
putdocx table tbl = etable	
putdocx table tbl(.,2/4), nformat(%6.2f)
putdocx table tbl(.,5), nformat(%5.3f)
putdocx table tbl(.,6/7), nformat(%6.2f)


* Median survival time (i.e. time in state)
putdocx pagebreak
putdocx paragraph, style("Heading2")
putdocx text ("Predicted marginal median time in state `state' based on interval censored `d' survival model")

putdocx paragraph, style("Heading3")
putdocx text ("Predicted marginal median times by genotype and heat condition")
margins mut#treat, predict(median) noempty post coeflegend
est sto `state'_margins

est resto `state'_model	
margins mut#treat, predict(median) noempty
cap putdocx table tbl = etable

mplotoffset, recast(scatter) graphregion(color(white) margin(r+5)) plotregion(margin(l+5 r+5)) xtitle("Genotype") ytitle("Predicted median time") title("``state'_name'") legend(col(4)) offset(0.1)
graph export MarginalMedians_`state'.png, replace
putdocx paragraph
putdocx image MarginalMedians_`state'.png, width(15cm)
putdocx table tbl(.,2/4), nformat(%6.2f)
putdocx table tbl(.,5), nformat(%5.3f)
putdocx table tbl(.,6/7), nformat(%6.2f)
			

putdocx pagebreak
putdocx paragraph, style("Heading3")
putdocx text ("Pairwise differences in predicted marginal medians")
putdocx paragraph

est resto `state'_model
margins, predict(median) at(treat=(1 3 5 6) mut=(1/5)) pwcompare noempty
	* This works fine, but does not give self-explaining labels.
* -> Calculate margins manually via lincom to have nicer labels.
est resto `state'_margins

* Pairwise comparisons of all 6 treatments within wildtype cells
local tmp_mut = 5
local row_counter = 1
putdocx table tbl = (1, 4)
putdocx table tbl(1,1) = ("Contrast") 
putdocx table tbl(1,2) = ("Estimated difference in median survival times") 
putdocx table tbl(1,3) = ("95% Confidence Interval"), colspan(2)
forvalues tmp_treat1 = 1/6 {
	forvalues tmp_treat2 = `tmp_treat1'/6 {
		cap lincom (_b[`tmp_treat1'.treat#`tmp_mut'.mut]-_b[`tmp_treat2'.treat#`tmp_mut'.mut])
		if (_rc==0 & !(missing(r(estimate)) | r(estimate)==0)) {
			putdocx table tbl(`row_counter',.), addrows(1,after)
			local row_counter = `row_counter'+1
			local lab_m: label (mut) `tmp_mut'
			local lab_t1: label (treat) `tmp_treat1'
			local lab_t2: label (treat) `tmp_treat2'
			putdocx table tbl(`row_counter',1) = ("`lab_m'#`lab_t1' - `lab_m'#`lab_t2'")						
			putdocx table tbl(`row_counter',2) = (`r(estimate)')
			putdocx table tbl(`row_counter',3) = (`r(lb)')
			putdocx table tbl(`row_counter',4) = (`r(ub)')
		}
		putdocx table tbl(.,2/4), nformat(%6.2f)
	}
}


* Pairwise comparison of the 2 available treatments within mutated cells
forvalues tmp_mut = 1/4 {
	local tmp_treat1 3
	local tmp_treat2 5
	cap lincom (_b[`tmp_treat1'.treat#`tmp_mut'.mut]-_b[`tmp_treat2'.treat#`tmp_mut'.mut])
	if (_rc==0 & !(missing(r(estimate)) | r(estimate)==0)) {
		local row_counter = 1
		putdocx table tbl = (1, 4)
		putdocx table tbl(1,1) = ("Contrast") 
		putdocx table tbl(1,2) = ("Estimated difference in median survival times") 
		putdocx table tbl(1,3) = ("95% Confidence Interval"), colspan(2)
		putdocx table tbl(`row_counter',.), addrows(1,after)
		local row_counter = `row_counter'+1
		local lab_m: label (mut) `tmp_mut'
		local lab_t1: label (treat) `tmp_treat1'
		local lab_t2: label (treat) `tmp_treat2'
		putdocx table tbl(`row_counter',1) = ("`lab_m'#`lab_t1' - `lab_m'#`lab_t2'")						
		putdocx table tbl(`row_counter',2) = (`r(estimate)')
		putdocx table tbl(`row_counter',3) = (`r(lb)')
		putdocx table tbl(`row_counter',4) = (`r(ub)')
	}
	putdocx table tbl(.,2/4), nformat(%6.2f)
}


* Pairwise comparisons of genotypes for the 2 treatments available in all cell lines
foreach tmp_treat in 3 5 {
	local row_counter = 1
	putdocx table tbl = (1, 4)
	putdocx table tbl(1,1) = ("Contrast") 
	putdocx table tbl(1,2) = ("Estimated difference in median survival times") 
	putdocx table tbl(1,3) = ("95% Confidence Interval"), colspan(2)
	forvalues tmp_mut1 = 1/5 {
		forvalues tmp_mut2 = `tmp_mut1'/5 {
			cap lincom (_b[`tmp_treat'.treat#`tmp_mut1'.mut]-_b[`tmp_treat'.treat#`tmp_mut2'.mut])
			if (_rc==0 & !(missing(r(estimate)) | r(estimate)==0)) {
				putdocx table tbl(`row_counter',.), addrows(1,after)
				local row_counter = `row_counter'+1
				local lab_m1: label (mut) `tmp_mut1'
				local lab_m2: label (mut) `tmp_mut2'
				local lab_t: label (treat) `tmp_treat'
				putdocx table tbl(`row_counter',1) = ("`lab_m1'#`lab_t' - `lab_m2'#`lab_t'")						
				putdocx table tbl(`row_counter',2) = (`r(estimate)')
				putdocx table tbl(`row_counter',3) = (`r(lb)')
				putdocx table tbl(`row_counter',4) = (`r(ub)')
			}
		}
	}
	putdocx table tbl(.,2/4), nformat(%6.2f)
}






********************************************************************************
********************************************************************************
********************************************************************************





****************************************************************
***                                                          ***
***                   Investigate state I                    ***
***   (MT array states 10-11 / telophase I - interkinesis)   ***
***                                                          ***
****************************************************************

local state "I"

putdocx pagebreak
putdocx paragraph, style("Heading1")
putdocx text ("``state'_name'")



************************************************************************
* Check whether there is sufficient information on all combinations of *
* genotype and heat condition                                          *
************************************************************************

* No. of events per combination of genotype and heat condition
bysort mut: tab treat ev_type_`state'
* -> Available information seems to be sufficient.

local ifexclude = ""
putdocx paragraph
putdocx text ("Remark:"), bold linebreak



********************************************
* Select best distribution in terms of AIC *
********************************************

* Prepare matrix for collecting AIC values
mat aicmat_`state' = J(1,5,.)
mat colnames aicmat_`state' = exponential gompertz loglogistic weibull lognormal
	
* Identify distributions that could not be estimated and write them to report
foreach dist in "exponential" "gompertz" "loglogistic" "weibull" "lognormal" {
	di ""
	di "**********     State `state', `dist' distribution     **********"
	di ""
	cap noi stintreg b5.treat##b5.mut `ifexclude', interval(`state'_time_l `state'_time_u) vce(cluster anther_id) dist(`dist') base difficult 
		* Remark: 
		* Control treatment is coded as '5' -> use 'b5.treat'
		* Wildtype cells are coded as '5' too -> use 'b5.mut' 
	cap noi estat ic
	mat ic = r(S)
	local aic = ic[1,colnumb(ic, "AIC")]
	mat aicmat_`state'[1,colnumb(aicmat_`state', "`dist'")] = `aic'
	local aic .
	di ""
	di "**********     END (State `state', `dist' distribution)     **********"
	di ""
	pause
}
mat list aicmat_`state'


* Identify distribution with best AIC and write it to report
qui mat list aicmat_`state'
local mcn : coln(aicmat_`state')

mata mata clear
mata x = st_matrix("aicmat_`state'")
mata xs = sort(x',1)
mata v = xs[(1), .]
mata min_aic_col = selectindex(x:==rowmin(v))
mata st_numscalar("min_aic", min_aic_col)
	
local mincol = scalar(min_aic) 
local mincolname: word `mincol' of `mcn'
local best_`state' = "`mincolname'"

local best_`state'_name = "`best_`state''"
if ("`best_`state'_name'"=="weibull") local best_`state'_name = "Weibull"

putdocx text ("- The distribution yielding the best AIC is `best_`state'_name'.")
			
* Clean up			
scalar drop _all   



	
***************
* Final model *
***************

* Specify best distribution chosen above
local dist = "`best_`state''"

* Mixed effects model
cap noi stintreg b5.treat##b5.mut `ifexclude', interval(`state'_time_l `state'_time_u) vce(cluster anther_id) dist(`dist') noempty difficult  
	* Remarks: 
	* - Control treatment is coded as '5' -> use 'b5.treat'
	*   Wildtype cells are coded as '5' too -> use 'b5.mut' 
est sto `state'_model
putdocx paragraph, style("Heading2")
local d = "`dist'"
if ("`dist'"=="weibull") local d = "Weibull"
putdocx text ("Interval censored parametric survival model with `d' distribution")
putdocx table tbl = etable	
putdocx table tbl(.,2/4), nformat(%6.2f)
putdocx table tbl(.,5), nformat(%5.3f)
putdocx table tbl(.,6/7), nformat(%6.2f)


* Median survival time (i.e. time in state)
putdocx pagebreak
putdocx paragraph, style("Heading2")
putdocx text ("Predicted marginal median time in state `state' based on interval censored `d' survival model")

putdocx paragraph, style("Heading3")
putdocx text ("Predicted marginal median times by genotype and heat condition")
margins mut#treat, predict(median) noempty post coeflegend
est sto `state'_margins

est resto `state'_model	
margins mut#treat, predict(median) noempty
cap putdocx table tbl = etable

mplotoffset, recast(scatter) graphregion(color(white) margin(r+5)) plotregion(margin(l+5 r+5)) xtitle("Genotype") ytitle("Predicted median time") title("``state'_name'") legend(col(4)) offset(0.1)
graph export MarginalMedians_`state'.png, replace
putdocx paragraph
putdocx image MarginalMedians_`state'.png, width(15cm)
putdocx table tbl(.,2/4), nformat(%6.2f)
putdocx table tbl(.,5), nformat(%5.3f)
putdocx table tbl(.,6/7), nformat(%6.2f)
			

putdocx pagebreak
putdocx paragraph, style("Heading3")
putdocx text ("Pairwise differences in predicted marginal medians")
putdocx paragraph

est resto `state'_model
margins, predict(median) at(treat=(1 3 5 6) mut=(1/5)) pwcompare noempty
	* This works fine, but does not give self-explaining labels.
* -> Calculate margins manually via lincom to have nicer labels.
est resto `state'_margins

* Pairwise comparisons of all 6 treatments within wildtype cells
local tmp_mut = 5
local row_counter = 1
putdocx table tbl = (1, 4)
putdocx table tbl(1,1) = ("Contrast") 
putdocx table tbl(1,2) = ("Estimated difference in median survival times") 
putdocx table tbl(1,3) = ("95% Confidence Interval"), colspan(2)
forvalues tmp_treat1 = 1/6 {
	forvalues tmp_treat2 = `tmp_treat1'/6 {
		cap lincom (_b[`tmp_treat1'.treat#`tmp_mut'.mut]-_b[`tmp_treat2'.treat#`tmp_mut'.mut])
		if (_rc==0 & !(missing(r(estimate)) | r(estimate)==0)) {
			putdocx table tbl(`row_counter',.), addrows(1,after)
			local row_counter = `row_counter'+1
			local lab_m: label (mut) `tmp_mut'
			local lab_t1: label (treat) `tmp_treat1'
			local lab_t2: label (treat) `tmp_treat2'
			putdocx table tbl(`row_counter',1) = ("`lab_m'#`lab_t1' - `lab_m'#`lab_t2'")						
			putdocx table tbl(`row_counter',2) = (`r(estimate)')
			putdocx table tbl(`row_counter',3) = (`r(lb)')
			putdocx table tbl(`row_counter',4) = (`r(ub)')
		}
		putdocx table tbl(.,2/4), nformat(%6.2f)
	}
}


* Pairwise comparison of the 2 available treatments within mutated cells
forvalues tmp_mut = 1/4 {
	local tmp_treat1 3
	local tmp_treat2 5
	cap lincom (_b[`tmp_treat1'.treat#`tmp_mut'.mut]-_b[`tmp_treat2'.treat#`tmp_mut'.mut])
	if (_rc==0 & !(missing(r(estimate)) | r(estimate)==0)) {
		local row_counter = 1
		putdocx table tbl = (1, 4)
		putdocx table tbl(1,1) = ("Contrast") 
		putdocx table tbl(1,2) = ("Estimated difference in median survival times") 
		putdocx table tbl(1,3) = ("95% Confidence Interval"), colspan(2)
		putdocx table tbl(`row_counter',.), addrows(1,after)
		local row_counter = `row_counter'+1
		local lab_m: label (mut) `tmp_mut'
		local lab_t1: label (treat) `tmp_treat1'
		local lab_t2: label (treat) `tmp_treat2'
		putdocx table tbl(`row_counter',1) = ("`lab_m'#`lab_t1' - `lab_m'#`lab_t2'")						
		putdocx table tbl(`row_counter',2) = (`r(estimate)')
		putdocx table tbl(`row_counter',3) = (`r(lb)')
		putdocx table tbl(`row_counter',4) = (`r(ub)')
	}
	putdocx table tbl(.,2/4), nformat(%6.2f)
}


* Pairwise comparisons of genotypes for the 2 treatments available in all cell lines
foreach tmp_treat in 3 5 {
	local row_counter = 1
	putdocx table tbl = (1, 4)
	putdocx table tbl(1,1) = ("Contrast") 
	putdocx table tbl(1,2) = ("Estimated difference in median survival times") 
	putdocx table tbl(1,3) = ("95% Confidence Interval"), colspan(2)
	forvalues tmp_mut1 = 1/5 {
		forvalues tmp_mut2 = `tmp_mut1'/5 {
			cap lincom (_b[`tmp_treat'.treat#`tmp_mut1'.mut]-_b[`tmp_treat'.treat#`tmp_mut2'.mut])
			if (_rc==0 & !(missing(r(estimate)) | r(estimate)==0)) {
				putdocx table tbl(`row_counter',.), addrows(1,after)
				local row_counter = `row_counter'+1
				local lab_m1: label (mut) `tmp_mut1'
				local lab_m2: label (mut) `tmp_mut2'
				local lab_t: label (treat) `tmp_treat'
				putdocx table tbl(`row_counter',1) = ("`lab_m1'#`lab_t' - `lab_m2'#`lab_t'")						
				putdocx table tbl(`row_counter',2) = (`r(estimate)')
				putdocx table tbl(`row_counter',3) = (`r(lb)')
				putdocx table tbl(`row_counter',4) = (`r(ub)')
			}
		}
	}
	putdocx table tbl(.,2/4), nformat(%6.2f)
}






********************************************************************************
********************************************************************************
********************************************************************************





****************************************************************
***                                                          ***
***                   Investigate state M2                   ***
***   (MT array states 12-13 / metaphase II - anaphase II)   ***
***                                                          ***
****************************************************************

local state "M2"

putdocx pagebreak
putdocx paragraph, style("Heading1")
putdocx text ("``state'_name'")



************************************************************************
* Check whether there is sufficient information on all combinations of *
* genotype and heat condition                                          *
************************************************************************

* No. of events per combination of genotype and heat condition
bysort mut: tab treat ev_type_`state'
* -> Available information seems to be sufficient.

local ifexclude = ""
putdocx paragraph
putdocx text ("Remarks:"), bold linebreak



********************************************
* Select best distribution in terms of AIC *
********************************************

* Prepare matrix for collecting AIC values
mat aicmat_`state' = J(1,5,.)
mat colnames aicmat_`state' = exponential gompertz loglogistic weibull lognormal
	
* Identify distributions that could not be estimated and write them to report
foreach dist in "exponential" "gompertz" "loglogistic" "weibull" "lognormal" {
	di ""
	di "**********     State `state', `dist' distribution     **********"
	di ""
	cap noi stintreg b5.treat##b5.mut `ifexclude', interval(`state'_time_l `state'_time_u) vce(cluster anther_id) dist(`dist') base difficult
		* Remark: 
		* Control treatment is coded as '5' -> use 'b5.treat'
		* Wildtype cells are coded as '5' too -> use 'b5.mut' 
	cap noi estat ic
	mat ic = r(S)
	local aic = ic[1,colnumb(ic, "AIC")]
	mat aicmat_`state'[1,colnumb(aicmat_`state', "`dist'")] = `aic'
	local aic .
	di ""
	di "**********     END (State `state', `dist' distribution)     **********"
	di ""
	pause
}
mat list aicmat_`state'


* Identify distribution with best AIC and write it to report
qui mat list aicmat_`state'
local mcn : coln(aicmat_`state')

mata mata clear
mata x = st_matrix("aicmat_`state'")
mata xs = sort(x',1)
mata v = xs[(1), .]
mata min_aic_col = selectindex(x:==rowmin(v))
mata st_numscalar("min_aic", min_aic_col)
	
local mincol = scalar(min_aic) 
local mincolname: word `mincol' of `mcn'
local best_`state' = "`mincolname'"

local best_`state'_name = "`best_`state''"
if ("`best_`state'_name'"=="weibull") local best_`state'_name = "Weibull"

putdocx text ("- The distribution yielding the best AIC is `best_`state'_name'.")
			
* Clean up			
scalar drop _all   



	
***************
* Final model *
***************

* Specify best distribution chosen above
local dist = "`best_`state''"

* Mixed effects model
cap noi stintreg b5.treat##b5.mut `ifexclude', interval(`state'_time_l `state'_time_u) vce(cluster anther_id) dist(`dist') noempty difficult 
	* Remarks: 
	* - Control treatment is coded as '5' -> use 'b5.treat'
	*   Wildtype cells are coded as '5' too -> use 'b5.mut' 
est sto `state'_model
putdocx paragraph, style("Heading2")
local d = "`dist'"
if ("`dist'"=="weibull") local d = "Weibull"
putdocx text ("Interval censored parametric survival model with `d' distribution")
putdocx table tbl = etable	
putdocx table tbl(.,2/4), nformat(%6.2f)
putdocx table tbl(.,5), nformat(%5.3f)
putdocx table tbl(.,6/7), nformat(%6.2f)


* Median survival time (i.e. time in state)
putdocx pagebreak
putdocx paragraph, style("Heading2")
putdocx text ("Predicted marginal median time in state `state' based on interval censored `d' survival model")

putdocx paragraph, style("Heading3")
putdocx text ("Predicted marginal median times by genotype and heat condition")
margins mut#treat, predict(median) noempty post coeflegend
est sto `state'_margins

est resto `state'_model	
margins mut#treat, predict(median) noempty
cap putdocx table tbl = etable

mplotoffset, recast(scatter) graphregion(color(white) margin(r+5)) plotregion(margin(l+5 r+5)) xtitle("Genotype") ytitle("Predicted median time") title("``state'_name'") legend(col(4)) offset(0.1)
graph export MarginalMedians_`state'.png, replace
putdocx paragraph
putdocx image MarginalMedians_`state'.png, width(15cm)
putdocx table tbl(.,2/4), nformat(%6.2f)
putdocx table tbl(.,5), nformat(%5.3f)
putdocx table tbl(.,6/7), nformat(%6.2f)
			

putdocx pagebreak
putdocx paragraph, style("Heading3")
putdocx text ("Pairwise differences in predicted marginal medians")
putdocx paragraph

est resto `state'_model
margins, predict(median) at(treat=(1 3 5 6) mut=(1/5)) pwcompare noempty
	* This works fine, but does not give self-explaining labels.
* -> Calculate margins manually via lincom to have nicer labels.
est resto `state'_margins

* Pairwise comparisons of all 6 treatments within wildtype cells
local tmp_mut = 5
local row_counter = 1
putdocx table tbl = (1, 4)
putdocx table tbl(1,1) = ("Contrast") 
putdocx table tbl(1,2) = ("Estimated difference in median survival times") 
putdocx table tbl(1,3) = ("95% Confidence Interval"), colspan(2)
forvalues tmp_treat1 = 1/6 {
	forvalues tmp_treat2 = `tmp_treat1'/6 {
		cap lincom (_b[`tmp_treat1'.treat#`tmp_mut'.mut]-_b[`tmp_treat2'.treat#`tmp_mut'.mut])
		if (_rc==0 & !(missing(r(estimate)) | r(estimate)==0)) {
			putdocx table tbl(`row_counter',.), addrows(1,after)
			local row_counter = `row_counter'+1
			local lab_m: label (mut) `tmp_mut'
			local lab_t1: label (treat) `tmp_treat1'
			local lab_t2: label (treat) `tmp_treat2'
			putdocx table tbl(`row_counter',1) = ("`lab_m'#`lab_t1' - `lab_m'#`lab_t2'")						
			putdocx table tbl(`row_counter',2) = (`r(estimate)')
			putdocx table tbl(`row_counter',3) = (`r(lb)')
			putdocx table tbl(`row_counter',4) = (`r(ub)')
		}
		putdocx table tbl(.,2/4), nformat(%6.2f)
	}
}


* Pairwise comparison of the 2 available treatments within mutated cells
forvalues tmp_mut = 1/4 {
	local tmp_treat1 3
	local tmp_treat2 5
	cap lincom (_b[`tmp_treat1'.treat#`tmp_mut'.mut]-_b[`tmp_treat2'.treat#`tmp_mut'.mut])
	if (_rc==0 & !(missing(r(estimate)) | r(estimate)==0)) {
		local row_counter = 1
		putdocx table tbl = (1, 4)
		putdocx table tbl(1,1) = ("Contrast") 
		putdocx table tbl(1,2) = ("Estimated difference in median survival times") 
		putdocx table tbl(1,3) = ("95% Confidence Interval"), colspan(2)
		putdocx table tbl(`row_counter',.), addrows(1,after)
		local row_counter = `row_counter'+1
		local lab_m: label (mut) `tmp_mut'
		local lab_t1: label (treat) `tmp_treat1'
		local lab_t2: label (treat) `tmp_treat2'
		putdocx table tbl(`row_counter',1) = ("`lab_m'#`lab_t1' - `lab_m'#`lab_t2'")						
		putdocx table tbl(`row_counter',2) = (`r(estimate)')
		putdocx table tbl(`row_counter',3) = (`r(lb)')
		putdocx table tbl(`row_counter',4) = (`r(ub)')
	}
	putdocx table tbl(.,2/4), nformat(%6.2f)
}


* Pairwise comparisons of genotypes for the 2 treatments available in all cell lines
foreach tmp_treat in 3 5 {
	local row_counter = 1
	putdocx table tbl = (1, 4)
	putdocx table tbl(1,1) = ("Contrast") 
	putdocx table tbl(1,2) = ("Estimated difference in median survival times") 
	putdocx table tbl(1,3) = ("95% Confidence Interval"), colspan(2)
	forvalues tmp_mut1 = 1/5 {
		forvalues tmp_mut2 = `tmp_mut1'/5 {
			cap lincom (_b[`tmp_treat'.treat#`tmp_mut1'.mut]-_b[`tmp_treat'.treat#`tmp_mut2'.mut])
			if (_rc==0 & !(missing(r(estimate)) | r(estimate)==0)) {
				putdocx table tbl(`row_counter',.), addrows(1,after)
				local row_counter = `row_counter'+1
				local lab_m1: label (mut) `tmp_mut1'
				local lab_m2: label (mut) `tmp_mut2'
				local lab_t: label (treat) `tmp_treat'
				putdocx table tbl(`row_counter',1) = ("`lab_m1'#`lab_t' - `lab_m2'#`lab_t'")						
				putdocx table tbl(`row_counter',2) = (`r(estimate)')
				putdocx table tbl(`row_counter',3) = (`r(lb)')
				putdocx table tbl(`row_counter',4) = (`r(ub)')
			}
		}
	}
	putdocx table tbl(.,2/4), nformat(%6.2f)
}






********************************************************************************
********************************************************************************
********************************************************************************





**********************************************
***                                        ***
***           Investigate state T          ***
***   (MT array state 14 / telophase II)   ***
***                                        ***
**********************************************

local state "T"

putdocx pagebreak
putdocx paragraph, style("Heading1")
putdocx text ("``state'_name'")



************************************************************************
* Check whether there is sufficient information on all combinations of *
* genotype and heat condition                                          *
************************************************************************

* No. of events per genotype x heat condition combination
bysort mut: tab treat ev_type_`state'
* -> There is insufficient information (no events) on HS34_A in SPO11. 

local ifexclude = "if !(treat==3 & mut==4)"
putdocx paragraph
putdocx text ("Remarks:"), bold linebreak
putdocx text ("- HS34_A in SPO11 cells has to be excluded from the analysis of state `state', because there are no completely (or closely) observed durations for this condition."), linebreak



********************************************
* Select best distribution in terms of AIC *
********************************************

* Prepare matrix for collecting AIC values
mat aicmat_`state' = J(1,5,.)
mat colnames aicmat_`state' = exponential gompertz loglogistic weibull lognormal
	
* Identify distributions that could not be estimated and write them to report
foreach dist in "exponential" "gompertz" "loglogistic" "weibull" "lognormal" {
	di ""
	di "**********     State `state', `dist' distribution     **********"
	di ""
	cap noi stintreg b5.treat##b5.mut `ifexclude', interval(`state'_time_l `state'_time_u) vce(cluster anther_id) dist(`dist') base difficult
		* Remark: 
		* Control treatment is coded as '5' -> use 'b5.treat'
		* Wildtype cells are coded as '5' too -> use 'b5.mut' 
	cap noi estat ic
	mat ic = r(S)
	local aic = ic[1,colnumb(ic, "AIC")]
	mat aicmat_`state'[1,colnumb(aicmat_`state', "`dist'")] = `aic'
	local aic .
	di ""
	di "**********     END (State `state', `dist' distribution)     **********"
	di ""
	pause
}
mat list aicmat_`state'


* Identify distribution with best AIC and write it to report
qui mat list aicmat_`state'
local mcn : coln(aicmat_`state')

mata mata clear
mata x = st_matrix("aicmat_`state'")
mata xs = sort(x',1)
mata v = xs[(1), .]
mata min_aic_col = selectindex(x:==rowmin(v))
mata st_numscalar("min_aic", min_aic_col)
	
local mincol = scalar(min_aic) 
local mincolname: word `mincol' of `mcn'
local best_`state' = "`mincolname'"

local best_`state'_name = "`best_`state''"
if ("`best_`state'_name'"=="weibull") local best_`state'_name = "Weibull"

putdocx text ("- The distribution yielding the best AIC is `best_`state'_name'.")
			
* Clean up			
scalar drop _all   



	
***************
* Final model *
***************

* Specify best distribution chosen above
local dist = "`best_`state''"

* Mixed effects model
cap noi stintreg b5.treat##b5.mut `ifexclude', interval(`state'_time_l `state'_time_u) vce(cluster anther_id) dist(`dist') noempty difficult 
	* Remarks: 
	* - Control treatment is coded as '5' -> use 'b5.treat'
	*   Wildtype cells are coded as '5' too -> use 'b5.mut' 
est sto `state'_model
putdocx paragraph, style("Heading2")
local d = "`dist'"
if ("`dist'"=="weibull") local d = "Weibull"
putdocx text ("Interval censored parametric survival model with `d' distribution")
putdocx table tbl = etable	
putdocx table tbl(.,2/4), nformat(%6.2f)
putdocx table tbl(.,5), nformat(%5.3f)
putdocx table tbl(.,6/7), nformat(%6.2f)


* Median survival time (i.e. time in state)
putdocx pagebreak
putdocx paragraph, style("Heading2")
putdocx text ("Predicted marginal median time in state `state' based on interval censored `d' survival model")

putdocx paragraph, style("Heading3")
putdocx text ("Predicted marginal median times by genotype and heat condition")
margins mut#treat, predict(median) noempty post coeflegend
est sto `state'_margins

est resto `state'_model	
margins mut#treat, predict(median) noempty
cap putdocx table tbl = etable

mplotoffset, recast(scatter) graphregion(color(white) margin(r+5)) plotregion(margin(l+5 r+5)) xtitle("Genotype") ytitle("Predicted median time") title("``state'_name'") legend(col(4)) offset(0.1)
graph export MarginalMedians_`state'.png, replace
putdocx paragraph
putdocx image MarginalMedians_`state'.png, width(15cm)
putdocx table tbl(.,2/4), nformat(%6.2f)
putdocx table tbl(.,5), nformat(%5.3f)
putdocx table tbl(.,6/7), nformat(%6.2f)

			
putdocx pagebreak
putdocx paragraph, style("Heading3")
putdocx text ("Pairwise differences in predicted marginal medians")
putdocx paragraph

est resto `state'_model
margins, predict(median) at(treat=(1 3 5 6) mut=(1/5)) pwcompare noempty
	* This works fine, but does not give self-explaining labels.
* -> Calculate margins manually via lincom to have nicer labels.
est resto `state'_margins


* Pairwise comparisons of all treatments within wildtype cells
local tmp_mut = 5
local row_counter = 1
putdocx table tbl = (1, 4)
putdocx table tbl(1,1) = ("Contrast") 
putdocx table tbl(1,2) = ("Estimated difference in median survival times") 
putdocx table tbl(1,3) = ("95% Confidence Interval"), colspan(2)
forvalues tmp_treat1 = 1/6 {
	forvalues tmp_treat2 = `tmp_treat1'/6 {
		cap lincom (_b[`tmp_treat1'.treat#`tmp_mut'.mut]-_b[`tmp_treat2'.treat#`tmp_mut'.mut])
		if (_rc==0 & !(missing(r(estimate)) | r(estimate)==0)) {
			putdocx table tbl(`row_counter',.), addrows(1,after)
			local row_counter = `row_counter'+1
			local lab_m: label (mut) `tmp_mut'
			local lab_t1: label (treat) `tmp_treat1'
			local lab_t2: label (treat) `tmp_treat2'
			putdocx table tbl(`row_counter',1) = ("`lab_m'#`lab_t1' - `lab_m'#`lab_t2'")						
			putdocx table tbl(`row_counter',2) = (`r(estimate)')
			putdocx table tbl(`row_counter',3) = (`r(lb)')
			putdocx table tbl(`row_counter',4) = (`r(ub)')
		}
		putdocx table tbl(.,2/4), nformat(%6.2f)
	}
}


* Pairwise comparison of the 2 available treatments within mutated cells
forvalues tmp_mut = 1/4 {
	local tmp_treat1 3
	local tmp_treat2 5
	cap lincom (_b[`tmp_treat1'.treat#`tmp_mut'.mut]-_b[`tmp_treat2'.treat#`tmp_mut'.mut])
	if (_rc==0 & !(missing(r(estimate)) | r(estimate)==0)) {
		local row_counter = 1
		putdocx table tbl = (1, 4)
		putdocx table tbl(1,1) = ("Contrast") 
		putdocx table tbl(1,2) = ("Estimated difference in median survival times") 
		putdocx table tbl(1,3) = ("95% Confidence Interval"), colspan(2)
		putdocx table tbl(`row_counter',.), addrows(1,after)
		local row_counter = `row_counter'+1
		local lab_m: label (mut) `tmp_mut'
		local lab_t1: label (treat) `tmp_treat1'
		local lab_t2: label (treat) `tmp_treat2'
		putdocx table tbl(`row_counter',1) = ("`lab_m'#`lab_t1' - `lab_m'#`lab_t2'")						
		putdocx table tbl(`row_counter',2) = (`r(estimate)')
		putdocx table tbl(`row_counter',3) = (`r(lb)')
		putdocx table tbl(`row_counter',4) = (`r(ub)')
	}
	putdocx table tbl(.,2/4), nformat(%6.2f)
}


* Pairwise comparisons of genotypes for the 2 treatments available in all cell lines
foreach tmp_treat in 3 5 {
	local row_counter = 1
	putdocx table tbl = (1, 4)
	putdocx table tbl(1,1) = ("Contrast") 
	putdocx table tbl(1,2) = ("Estimated difference in median survival times") 
	putdocx table tbl(1,3) = ("95% Confidence Interval"), colspan(2)
	forvalues tmp_mut1 = 1/5 {
		forvalues tmp_mut2 = `tmp_mut1'/5 {
			cap lincom (_b[`tmp_treat'.treat#`tmp_mut1'.mut]-_b[`tmp_treat'.treat#`tmp_mut2'.mut])
			if (_rc==0 & !(missing(r(estimate)) | r(estimate)==0)) {
				putdocx table tbl(`row_counter',.), addrows(1,after)
				local row_counter = `row_counter'+1
				local lab_m1: label (mut) `tmp_mut1'
				local lab_m2: label (mut) `tmp_mut2'
				local lab_t: label (treat) `tmp_treat'
				putdocx table tbl(`row_counter',1) = ("`lab_m1'#`lab_t' - `lab_m2'#`lab_t'")						
				putdocx table tbl(`row_counter',2) = (`r(estimate)')
				putdocx table tbl(`row_counter',3) = (`r(lb)')
				putdocx table tbl(`row_counter',4) = (`r(ub)')
			}
		}
	}
	putdocx table tbl(.,2/4), nformat(%6.2f)
}






********************************************************************************
********************************************************************************
********************************************************************************





*****************************************************************
***                                                           ***
***            Investigate state range H to M2                ***
***   (MT array states 2-13 / late leptotene - anaphase II)   ***
***                                                           ***
*****************************************************************

local state "HM2"

putdocx pagebreak
putdocx paragraph, style("Heading1")
putdocx text ("``state'_name'")



************************************************************************
* Check whether there is sufficient information on all combinations of *
* genotype and heat condition                                          *
************************************************************************

* No. of events per combination of genotype and heat condition
bysort mut: tab treat ev_type_`state'
* -> There is insufficient information (no or hardly any events) on 
*    control in DMC as well as HS30_B and HS34_B in wildtype. 

local ifexclude = "if !((treat==5 & mut==2) | (inlist(treat, 2, 4) & mut==5))"
putdocx paragraph
putdocx text ("Remarks:"), bold linebreak
putdocx text ("- For the following conditions we have no or hardly any completely observed durations, so that we cannot use them for estimation:"), linebreak
putdocx text ("  * control in DMC cells"), linebreak
putdocx text ("  * HS30_B and HS34_B in wildtype cells"), linebreak




********************************************
* Select best distribution in terms of AIC *
********************************************

* Prepare matrix for collecting AIC values
mat aicmat_`state' = J(1,5,.)
mat colnames aicmat_`state' = exponential gompertz loglogistic weibull lognormal
	
* Identify distributions that could not be estimated and write them to report
foreach dist in "exponential" "gompertz" "loglogistic" "weibull" "lognormal" {
	di ""
	di "**********     State `state', `dist' distribution     **********"
	di ""
	cap noi stintreg b5.treat##b5.mut `ifexclude', interval(`state'_time_l `state'_time_u) vce(cluster anther_id) dist(`dist') base difficult
		* Remark: 
		* Control treatment is coded as '5' -> use 'b5.treat'
		* Wildtype cells are coded as '5' too -> use 'b5.mut' 
	cap noi estat ic
	mat ic = r(S)
	local aic = ic[1,colnumb(ic, "AIC")]
	mat aicmat_`state'[1,colnumb(aicmat_`state', "`dist'")] = `aic'
	local aic .
	di ""
	di "**********     END (State `state', `dist' distribution)     **********"
	di ""
	pause
}
mat list aicmat_`state'


* Identify distribution with best AIC and write it to report
qui mat list aicmat_`state'
local mcn : coln(aicmat_`state')

mata mata clear
mata x = st_matrix("aicmat_`state'")
mata xs = sort(x',1)
mata v = xs[(1), .]
mata min_aic_col = selectindex(x:==rowmin(v))
mata st_numscalar("min_aic", min_aic_col)
	
local mincol = scalar(min_aic) 
local mincolname: word `mincol' of `mcn'
local best_`state' = "`mincolname'"

local best_`state'_name = "`best_`state''"
if ("`best_`state'_name'"=="weibull") local best_`state'_name = "Weibull"

putdocx text ("- The distribution yielding the best AIC was `best_`state'_name'."), linebreak
			
* Clean up			
scalar drop _all   



	
***************
* Final model *
***************

* Specify best distribution chosen above
local dist = "`best_`state''"

* Mixed effects model
cap noi stintreg b5.treat##b5.mut `ifexclude', interval(`state'_time_l `state'_time_u) vce(cluster anther_id) dist(`dist') noempty difficult 
	* Remarks: 
	* - Control treatment is coded as '5' -> use 'b5.treat'
	*   Wildtype cells are coded as '5' too -> use 'b5.mut' 
est sto `state'_model
putdocx text ("- HS34_A in DMC was omitted from the model automatically due to collinearity.")

putdocx paragraph, style("Heading2")
local d = "`dist'"
if ("`dist'"=="weibull") local d = "Weibull"
putdocx text ("Interval censored parametric survival model with `d' distribution")
putdocx table tbl = etable	
putdocx table tbl(.,2/4), nformat(%6.2f)
putdocx table tbl(.,5), nformat(%5.3f)
putdocx table tbl(.,6/7), nformat(%6.2f)


* Median survival time (i.e. time in state)
putdocx pagebreak
putdocx paragraph, style("Heading2")
putdocx text ("Predicted marginal median time in state `state' based on interval censored `d' survival model")

putdocx paragraph, style("Heading3")
putdocx text ("Predicted marginal median times by genotype and heat condition")
margins mut#treat, predict(median) noempty post coeflegend
est sto `state'_margins

est resto `state'_model	
margins mut#treat, predict(median) noempty
cap putdocx table tbl = etable

mplotoffset, recast(scatter) graphregion(color(white) margin(r+5)) plotregion(margin(l+5 r+5)) xtitle("Genotype") ytitle("Predicted median time") title("``state'_name'") legend(col(4)) offset(0.1)
graph export MarginalMedians_`state'.png, replace
putdocx paragraph
putdocx image MarginalMedians_`state'.png, width(15cm)
putdocx table tbl(.,2/4), nformat(%6.2f)
putdocx table tbl(.,5), nformat(%5.3f)
putdocx table tbl(.,6/7), nformat(%6.2f)





********************************************************************************
********************************************************************************
********************************************************************************



putdocx save Report_analysis.docx, replace
log close








