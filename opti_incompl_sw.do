mata: mata set matastrict on

mata:


/* GREEDY ALGORITHM TO SEARCH FOR OPTIMAL DESIGN */

/* seq  is a matrix of the permissible sequences */
/* initweight  is a vector of initial weights */
/* r  is the correlation between cluster-period means */
/* fineness is the fineness of changes per iteration */
/* (large value for fineness = small changes) */
/* maxiter  is the maximum number of iterations */
/* (set maxiter to missing for no limit) */

real colvector opti_incompl_sw(real matrix seq, ///
		real colvector initweight, real scalar r, ///
		real scalar fineness, real scalar maxiter)
{
	real colvector ncompl, weight, altwt, nextwt
	real scalar nper, nseq, inc, inc1, i, j, iter, ///
			var, altvar, nextvar, reweight
	
	nper = cols(seq)
	nseq = rows(seq)
	ncompl = rownonmissing(seq)
	reweight = sum(initweight :* ncompl) / nper
	inc = 1 / fineness
	
	weight = initweight :/ reweight
	var = glsvar(seq, weight, r)
	nextvar = var
	
	for (iter=1; iter<maxiter; iter++) {
		altwt = weight
		for (i=1; i<=nseq; i++) {
			if (weight[i]>0) {
				if (weight[i]*ncompl[i] < inc) {
					altwt[i] = 0
					inc1 = weight[i]*ncompl[i]
				}
				else {
					altwt[i] = altwt[i] - inc/ncompl[i]
					inc1 = inc
				}
				for (j=1; j<=nseq; j++) {
					if (j!=i) {
						altwt[j] = altwt[j] + ///
								inc1/ncompl[j]
						altvar = glsvar(seq, ///
								altwt, r)
						if (altvar < nextvar) {
							nextvar = altvar
							nextwt = altwt
						}
						altwt[j] = weight[j]
					}
				}
				altwt[i] = weight[i]
			}
		}
		if (nextvar < var) {
			weight = nextwt
			var = nextvar
		}
		else {
			iter=maxiter
		}
	}

	return(weight)
}


/* VARIANCE OF GLS ESTIMATOR */

real scalar glsvar(real matrix seq, ///
		real colvector weight, real scalar r)
{
	real matrix x, y, x1, y1
	real colvector ncompl
	real scalar nper, nseq, n, i, k, var
	
	nper = cols(seq)
	nseq = rows(seq)
	ncompl = rownonmissing(seq)
	n = sum( weight :* ncompl )
	
	y = J(nper+1, 0, .)
	x = J(0, nper+1, .)
	
	for (i=1; i<=nseq; i++) {
		if (weight[i]>0) {
			k = ncompl[i]
			x1 = select((seq[i,.]' , I(nper)), ///
					seq[i,.]' :!= .)
			y1 = x1' * invsym((I(k) :* (1-r)) + ///
					J(k, k, r)) :* weight[i]
			x = x \ x1
			y = y , y1
		}
	}
	
	var = invsym(y*x)[1,1] * n
	return(var)
}


/* ACCESSORIES */

/* Matrix of all incomplete sequences */

real matrix incompleteseq(real scalar nper)
{
	real scalar i, i1, j, nj, k1, k2
	real rowvector jlist
	real matrix seq, newseq
	
	seq = J(0, nper, .)
	for (i=1; i<=2^nper-1; i++) {
		i1 = i
		jlist = J(1, 0, .)
		for (j=1; j<=nper; j++) {
			if ( mod(i1,2) == 1 ) {
				jlist = jlist, j
			}
			i1 = floor(i1/2)
		}
		nj = cols(jlist)
		newseq = J(nj+1, nper, .)
		for (k1=1; k1<=nj+1; k1++) {
			for (k2=1; k2<=nj; k2++) {
				newseq[k1, jlist[k2]] = (k2>=k1)
			}
		}
		seq = seq \ newseq
	}
	return(seq)
}

/* Vector of weights for staircase design */
/* Could be used as starting weights */

real colvector staircaseweight(real matrix seq)
{
	real matrix diffseq
	real colvector weight
	real scalar nper, nseq, i
	
	nper = cols(seq)
	nseq = rows(seq)
	diffseq = J(nseq, 0, .)
	
	for (i=2; i<=nper; i++) {
		diffseq = diffseq, ( seq[.,i] - seq[.,i-1] )
	}
	
	weight = (( rownonmissing(seq) :== 2 ) :& ///
			( rowsum(diffseq) :== 1 )) :| ///
			(( rownonmissing(seq) :== 1 ) :& ///
			(( seq[.,1] :== 1 ) :| ///
			( seq[.,nper] :== 0 )))
		
	return(weight)
}

/* Vector of weights for classic stepped wedge */
/* Could be used as starting weights */

real colvector classicswweight(real matrix seq)
{
	real scalar nper
	real colvector weight
	
	nper = cols(seq)
	weight = ( rownonmissing(seq) :== nper ) :& ///
			( rowmin(seq) :!= rowmax(seq) )
			
	return(weight)
}

end


/* EXAMPLE */

/* 5-period incomplete designs */
/* Starts with staircase design */
/* r=0.8 */

mata:

seq = incompleteseq(5)
initweight = staircaseweight(seq)
weight = opti_incompl_sw(seq, initweight, 0.8, 1000, .)

seq, weight


/* ADD OTHER CALLS TO opti_incompl_sw AS NEEDED */


end

