*! wentropy v1 January 14, 2022
* Paul Corral - WBG - Poverty and Equity GP
* Rodirgo Salcedo Dubois

cap prog drop wentropy
cap set matastrict off
program define wentropy, eclass
	version 11.2
	#delimit ;
	syntax varlist (min=1 numeric) [if] [in],
	NEWweight(string)
	constraints(name)
	[OLDweight(varlist max=1 numeric) 
	POPtotal(real 1)
	iter(integer 10000)
	missok];
	
	#delimit cr		
marksample touse

local state = c(rngstate)

// Remove observations if prior is missing
if ("`oldweight'"!=""){
	replace `touse' = 0 if missing(`oldweight')
}
	
	//Send data to MATA	
	//Data
	mata: st_view(X=.,.,"`varlist'","`touse'")
	//Prior weights
	if ("`oldweight'"=="") mata: q = J(rows(X),1,1):/rows(X)
	else{
		mata: q = st_data(.,"`oldweight'","`touse'")
		mata: q = q:/quadsum(q)
	}
	//Constraint matrix
	mata: y = st_matrix("`constraints'")
	
	*===========================================================================
	// Data checks
	*===========================================================================
	//Ensure the constraints are the same number as the Xs and a column vector
		//First column vector
		if (colsof(`constraints')>1){
			dis as error "Constraints matrix must be a column vector"
			error 198
			exit
		}
		local numx: list sizeof varlist
		if (rowsof(`constraints')!=`numx'){
			dis as error "Constraints must match the number of variables"
			error 198
			exit			
		}
	// Ensure constraints are not out of bounds
	local a = 1
	local thereiserror = 0
	foreach x of local varlist{
		qui:sum `x'
		if (~inrange(`constraints'[`a',1],r(min),r(max))){
			local thereiserror = 1
			dis as error "`x''s constraint is out of bounds'"
		}
		local ++a
	}
	if (`thereiserror'==1){
		error 198 
		exit
	}
	if ("`missok'"!=""){
		local miss = 1
	}
	else{
		local miss = 0
	}
	
	qui:gen double `newweight' = .
		lab var `newweight' "wentropy calibrated weights"
	mata: st_store(., tokens("`newweight'"), "`touse'",_wentropy(y,q,X, `iter', `miss'))
	if (~missing("`poptotal'")) qui:replace `newweight' = `poptotal'*`newweight'

set rngstate `state'	
end

mata

function _wentropy(y,q,X, iter, miss){
	display("Newton-Rhapson - Own")
	p = _myNR(y,q,X,10)
	//ALternative solution if _myNR didn't work
	if (hasmissing(p) & miss==0){
		quadsum((p:==.))
		display("Broyden–Fletcher–Goldfarb–Shanno")
		s=optimize_init()
		optimize_init_evaluator(s,&_other_solver())
		optimize_init_evaluatortype(s,"d2")
		optimize_init_which(s,"min")
		optimize_init_technique(s, "bfgs")
		optimize_init_singularHmethod(s, "hybrid")
		optimize_init_argument(s,1,y)
		optimize_init_argument(s,2,X)
		optimize_init_argument(s,3,q)
		optimize_init_conv_maxiter(s,iter)
		optimize_init_params(s,J(1,rows(y),0))
		g = optimize_result_gradient(s)
		B=optimize(s)'
		eXl		= exp(-quadcross(X',B))
		Omega	= mean(eXl, q)
		p		= q :* eXl / Omega
		if (quadsum(abs(g))>0.05){
			p = J(rows(p),1,.)
			display("Sum of gradient over 0.5")
		}
	}
	if (hasmissing(p) & miss==0){
		display("Davidon–Fletcher–Powell")
		s=optimize_init()
		optimize_init_evaluator(s,&_other_solver())
		optimize_init_evaluatortype(s,"d2")
		optimize_init_which(s,"min")
		optimize_init_technique(s, "dfp")
		optimize_init_singularHmethod(s, "hybrid")
		optimize_init_argument(s,1,y)
		optimize_init_argument(s,2,X)
		optimize_init_argument(s,3,q)
		optimize_init_conv_maxiter(s,iter)
		optimize_init_params(s,J(1,rows(y),0))
		g = optimize_result_gradient(s)
		B=optimize(s)'
		eXl		= exp(-quadcross(X',B))
		Omega	= mean(eXl, q)
		p		= q :* eXl / Omega
		if (quadsum(abs(g))>0.05){
			p = J(rows(p),1,.)
			display("Sum of gradient over 0.5")
		}
	}
	if (hasmissing(p) & miss==0){
		display("modified Newton–Raphson")
		s=optimize_init()
		optimize_init_evaluator(s,&_other_solver())
		optimize_init_evaluatortype(s,"d2")
		optimize_init_which(s,"min")
		optimize_init_technique(s, "nr")
		optimize_init_singularHmethod(s, "hybrid")
		optimize_init_argument(s,1,y)
		optimize_init_argument(s,2,X)
		optimize_init_argument(s,3,q)
		optimize_init_conv_maxiter(s,iter)
		optimize_init_params(s,J(1,rows(y),0))
		g = optimize_result_gradient(s)
		B=optimize(s)'
		eXl		= exp(-quadcross(X',B))
		Omega	= mean(eXl, q)
		p		= q :* eXl / Omega
		if (quadsum(abs(g))>0.05){
			p = J(rows(p),1,.)
			display("Sum of gradient over 0.5")
		}
	}	
	if (hasmissing(p) & miss==0) _error(999, "Convergence was not achieved")
	return(p)
	
}

function _myNR(y,q,X, iter){
	change = 1
	B      = J(rows(y),1,-1)
	tries  = 0
	rseed(23736)
	//Newton Rhapson
	while (abs(change)>1e-10 & tries<iter){	
		eXl		= exp(-quadcross(X',B))
		Omega	= mean(eXl, q)
		p		= q :* eXl / Omega
		mX	    = mean(X, p)
		g	    = y - mX'	        // gradient		
		H	    = ((quadcross(p:* X, X) - quadcross(mX, mX)))	// Hessian	
		//iteration counter
		tries    = tries+1
		//Old lambdas (called here B)
		bold    = B
		//New lambdas (called here B)
		_invsym(H)
		B       = bold - (quadcross(H,g))
		//Change to determine convergence
		change=quadcross((bold-B),(bold-B))/(quadcross(bold,bold))
		//Change starting values for B in case of problems
		//if (hasmissing(H) | allof(H,0) | hasmissing(g) ){
		//	tries = tries+1
		//	change = 10
		//	if (tries==1)  B = J(rows(y),1,1)
		//	if (tries==2)  B = J(rows(y),1,0)
		//	if (tries>=3){
		//		if (mod(tries,2)==0) S = 1
		//		else S = -1
		//		B = S*runiform(rows(y),1):*((tries-2)*10)
		//	}
		//	display("Tries :")
		//	tries
		//}	
	}
	eXl		= exp(-quadcross(X',B))
	Omega	= mean(eXl, q)
	p		= q :* eXl / Omega
	if (quadsum(abs(g))>0.05){
		p = J(rows(p),1,.)
		display("Sum of gradient over 0.05")
	}
	return(p)
}

function _other_solver(todo, B,y,X,q,L,g,H){
	B       = B'
	eXl		= exp(-quadcross(X',B))
	Omega	= mean(eXl, q)
	p		= q :* eXl / Omega
	mX	    = mean(X, p)
	L = quadcross(y,B) + ln(Omega)
	if (todo>=1) g = (y - mX')'
	if (todo==2) H	    = ((quadcross(p:* X, X) - quadcross(mX, mX)))	
	B = B'
}



end




