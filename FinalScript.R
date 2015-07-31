####################################################
#
# R script for final exam. 
# by Dan Jin. 
# last update: Dec. 12, 2012
# 
# To run this script, please do the following step in R.
#   1. source("FinalScript.r")
#   2. result.ls=list() # it's better to sign the result to a list, since it will return large data set.
#   3. result.ls=control.fun()
#
# Result of this script:
#   1. Return one list, result.ls, that includes all results required.
#   2. Create 7 plots in png format. (See below for more details)
#
#
# Summary of the content in result.ls:
#   1. $covar: beta.MLE, means of the absolute value of beta and p values calculated with covariance Xz 
#   2. $nocovar: beta.MLE, means of the absolute value of beta and p values calculated without covariance Xz
#   3. $hits.ls: all genotype marker (marker number and their p value) significant than Bonferroni correction;
#   4. $hit.set.ls: distinct sets;
#   5. $sig.hit.ls: the most significant marker in each set (marker number and the associated p value);
#   6. $epi.X: a list containing the nine dummy variables for the first five individuals for epistasis analysis.
#
#
#Structure of result.ls:
# result.ls
# |-$covar
# | |-beta.MLE
# | |-beta.a.mean
# | |-beta.d.mean
# | |-beta.z.mean
# | |-pval
# |
# |-$nocovar
# | |-beta.minus.MLE
# | |-beta.a.minus.mean
# | |-beta.d.minus.mean
# | |-pval
# |
# |-$hits.ls
# | |-$position
# | |-$pval
# |
# |-$hit.set.ls
# |
# |-$sig.hit.ls
# | |-$position
# | |-$pval
# |
# |-$epi.X
#
#
# This script will also create 7 plots in png format.
#   1. Q1_Phenotype distribution.png
#   2. Q3_PCA.png
#   3. Q6.1_Histogram of p value (covariance).png
#   3. Q6.2_Histogram of p value (no covariance).png
#   5. Q7.1_QQ Plot with covariate.png
#   6. Q7.2_QQ Plot without covariate.png
#	7. Q8_Manhattan Plot with Bonferroni correction (covariance).png
#
####################################################



#:::::::::::::::::Import & convert data:::::::::::::::::

import.sample.fun <- function() {
	cat("Importing data.\n")
	phenodata.df = read.table("final_phenotype_fall12.txt", header = FALSE, sep = "\n")
	genodata.df = read.table("final_genotype_fall12.txt", header = FALSE, sep = "")

	size = dim(genodata.df)
	sample.size = size[1]/2 # = 499
	num.sample = size[2] # = 2480
	cat("The sample size (n) =",sample.size,"\n")
	cat("The number of genotypes (N) =",num.sample,"\n\n")

	sample.ls = list()
	sample.ls$y = as.matrix(phenodata.df)
	cat("Q1: Check phenotype distribution.\n")
	png("Q1_Phenotype distribution.png", width = 500, height = 500)
	hist(sample.ls$y, main = "Phenotype distribution", ylab = "Number of individuals", xlab = "Phenotype (0=control;1=case)")
	dev.off()
	cat("\"Q1_Phenotype distribution.png\" created.\n\n")
	sample.ls$g = genodata.df
	cat("Q1: Check genotype data.\n")
	if ((sum(genodata.df != "A" & genodata.df != "T" & genodata.df != "G" & genodata.df != "C") == 0) & length(sample.ls$y) == 
		sample.size) {
		cat("No missing genotype data.\n\n")
	}
	return(sample.ls)
}


#:::::::::::::::::Filter data:::::::::::::::::

MAF.fun <- function(genodata.df) {
	threshold.MAF = 0.05
	cat("Q2: Remove genotype with MAF<", threshold.MAF * 100, "%.\n")
	MAF.vec = c()
	MAF.fun <- function(col) {
		unique.allele = unique(col) # find the unique allele
		unique.allele.non0 = unique.allele[which(unique(col) != 0)] # remove zero
		if (length(unique.allele.non0) == 1) {
			return(0)
		} else {
			allele.freq = c(sum(col == unique.allele.non0[1]), sum(col == unique.allele.non0[2])) # calculate allele frequence
			return(min(allele.freq)/length(col))
		}
	}
	MAF.vec = apply(genodata.df, 2, MAF.fun)
	MAF.bool = MAF.vec >= threshold.MAF
	cat(sum(MAF.bool == 1), "SNPs (", sum(MAF.bool == 1)/length(MAF.bool) * 100, "%) have MAF >=", threshold.MAF * 
		100, "%.\n\n")
	filtered.geno = genodata.df[, which(MAF.bool == 1, arr.ind = TRUE)]
	return(filtered.geno)
}


#:::::::::::::::::Convert genotype to Xa and Xd:::::::::::::::::

Xa.coding.fun <- function(filtered.geno) {
	cat("Coding Xa and Xd...\n\n")
	size = dim(filtered.geno)
	sample.size = size[1]/2
	num.sample = size[2]

	Xa = matrix(0, sample.size, num.sample)
	for (col in 1:num.sample) {
		unique.allele = unique(filtered.geno[, col]) # find the unique allele
		allele.freq = c(sum(filtered.geno[, col] == unique.allele[1]), sum(filtered.geno[, col] == unique.allele[2])) # calculate allele frequence
		unique.allele = unique.allele[order(allele.freq)]

		for (row in 1:sample.size) {
			if (filtered.geno[2 * row - 1, col] == (filtered.geno[2 * row, col])) {
				if (filtered.geno[2 * row - 1, col] == unique.allele[2]) {
					Xa[row, col] = 1
				} else if (filtered.geno[2 * row - 1, col] == unique.allele[1]) {
					Xa[row, col] = -1
				}
			}

		}
	}
	return(Xa)
}


#:::::::::::::::::PCA:::::::::::::::::

PCA.fun <- function(Xa) {
	#scale the genotypes by their mean and standard deviation
	X = t(Xa)
	W <- (X - rowMeans(X))/sqrt(diag((cov(X))))

	#perform a PCA
	geno.pc <- princomp(W)

	#plot the loadings of the first two PCs
	cat("Q3: PCA.\n")
	png("Q3_PCA.png", width = 1000, height = 500)
	par(mfrow = c(1, 2))
	plot(geno.pc$loadings[, c(1, 2)], main = "PCA")
	hist(geno.pc$loadings[, 1], breaks = dim(Xa)[1], main = "Distribution of Comp.1", xlab = "Comp.1", ylab = "Frequency")
	dev.off()
	cat("\"Q3_PCA.png\" created.\n")

	return(geno.pc)
}


################################## GWAS with Covariance ##################################
result.ls = list()

#::::::::::::::::::::IRLS function::::::::::::::::::::
logistic.IRLS <- function(Y, X, beta.0 = c(0, 0, 0, 0), D = 1e-06, it.max = 100) {
	sample.size = length(Y)

	beta.t1 = c() # beta[t]
	beta.t2 = c() # beta[t+1]
	gamma.t1 = c() #gamma.t1=gamma^(-1)(Xbeta.t1)
	gamma.t2 = c() #gamma.t2=gamma^(-1)(Xbeta.t2)

	# starting value
	beta.t1 = beta.0
	gamma.t1 = 1/(1 + exp(-X %*% beta.0))

	W = matrix(0, sample.size, sample.size)
	for (i in 1:sample.size) {
		W[i, i] = gamma.t1[i] %*% (1 - gamma.t1[i])
	}

	beta.t2 = beta.t1 + solve(t(X) %*% W %*% X) %*% t(X) %*% (Y - gamma.t1)
	gamma.t2 = exp(X %*% beta.t2)/(1 + exp(X %*% beta.t2))

	D1 = D2 = 0
	for (i in 1:sample.size) {
		if (Y[i]) {
			D1 = D1 + 2 * Y[i] * log(Y[i]/gamma.t1[i])
			D2 = D2 + 2 * Y[i] * log(Y[i]/gamma.t2[i])
		} else if (1 - Y[i]) {
			D1 = D1 + 2 * (1 - Y[i]) * log((1 - Y[i])/(1 - gamma.t1[i]))
			D2 = D2 + 2 * (1 - Y[i]) * log((1 - Y[i])/(1 - gamma.t2[i]))
		}

	}
	deltaD = abs(D2 - D1)

	# loop until deltaD is smaller than 1e-06 or the given D
	j = 1
	while ((deltaD >= D) & (j <= it.max)) {
		beta.t1 = beta.t2
		gamma.t1 = 1/(1 + exp(-X %*% beta.t1))
		W = matrix(0, sample.size, sample.size)
		for (i in 1:sample.size) {
			W[i, i] = gamma.t1[i] %*% (1 - gamma.t1[i])
		}

		beta.t2 = beta.t1 + solve(t(X) %*% W %*% X) %*% t(X) %*% (Y - gamma.t1)
		gamma.t2 = 1/(1 + exp(-X %*% beta.t2))

		D1 = D2 = 0
		for (i in 1:sample.size) {
			if (Y[i]) {
				D1 = D1 + 2 * Y[i] * log(Y[i]/gamma.t1[i])
				D2 = D2 + 2 * Y[i] * log(Y[i]/gamma.t2[i])
			} else if (1 - Y[i]) {
				D1 = D1 + 2 * (1 - Y[i]) * log((1 - Y[i])/(1 - gamma.t1[i]))
				D2 = D2 + 2 * (1 - Y[i]) * log((1 - Y[i])/(1 - gamma.t2[i]))
			}

		}
		deltaD = abs(D2 - D1)
		j = j + 1
	}
	beta.MLE = beta.t2

	return(beta.MLE)
}


#::::::::::::::::::::LRT function::::::::::::::::::::

logistic.LRT <- function(Y, Xa, Xd, Xz, beta.0 = c(0, 0, 0, 0), D = 1e-06, it.max = 100) {
	LRT.ls = list()
	beta1.hat.vec = c()

	sample.size = length(Y)
	ones = rep(1, sample.size)
	X0 = cbind(ones, Xz)
	XA = cbind(ones, Xa, Xd, Xz)

	#beta0.hat	
	beta.mu.z.hat = logistic.IRLS(Y, X0, c(0, 0), D, it.max)
	beta0.hat = c(beta.mu.z.hat[1], 0, 0, beta.mu.z.hat[2])
	gamma0 = 1/(1 + exp(-XA %*% beta0.hat))

	#beta1.hat
	beta1.hat = logistic.IRLS(Y, XA, beta.0, D, it.max)
	beta1.hat.vec = c(beta1.hat.vec, beta1.hat)
	gamma1 = 1/(1 + exp(-XA %*% beta1.hat))

	#LRT
	l.beta0.hat = 0
	l.beta1.hat = 0
	for (i in 1:length(Y)) {
		l.beta0.hat = l.beta0.hat + (Y[i] * log(gamma0[i]) + (1 - Y[i]) * log(1 - gamma0[i]))
		l.beta1.hat = l.beta1.hat + (Y[i] * log(gamma1[i]) + (1 - Y[i]) * log(1 - gamma1[i]))
	}
	LRT = 2 * (l.beta1.hat - l.beta0.hat)

	#p value
	LRT.ls$pval = pchisq(LRT, 2, lower.tail = FALSE, log.p = FALSE)
	LRT.ls$beta.MLE = beta1.hat.vec

	return(LRT.ls)
}



################################## GWAS without Covariance ##################################

#::::::::::::::::::::IRLS function::::::::::::::::::::
logistic.IRLS.minus <- function(Y, X, beta.0 = c(0, 0, 0), D = 1e-06, it.max = 100) {
	sample.size = length(Y)

	beta.t1 = c() # beta[t]
	beta.t2 = c() # beta[t+1]
	gamma.t1 = c() #gamma.t1=gamma^(-1)(Xbeta.t1)
	gamma.t2 = c() #gamma.t2=gamma^(-1)(Xbeta.t2)

	# starting value
	beta.t1 = beta.0
	gamma.t1 = 1/(1 + exp(-X %*% beta.0))

	W = matrix(0, sample.size, sample.size)
	for (i in 1:sample.size) {
		W[i, i] = gamma.t1[i] %*% (1 - gamma.t1[i])
	}

	beta.t2 = beta.t1 + solve(t(X) %*% W %*% X) %*% t(X) %*% (Y - gamma.t1)
	gamma.t2 = exp(X %*% beta.t2)/(1 + exp(X %*% beta.t2))

	D1 = D2 = 0
	for (i in 1:sample.size) {
		if (Y[i]) {
			D1 = D1 + 2 * Y[i] * log(Y[i]/gamma.t1[i])
			D2 = D2 + 2 * Y[i] * log(Y[i]/gamma.t2[i])
		} else if (1 - Y[i]) {
			D1 = D1 + 2 * (1 - Y[i]) * log((1 - Y[i])/(1 - gamma.t1[i]))
			D2 = D2 + 2 * (1 - Y[i]) * log((1 - Y[i])/(1 - gamma.t2[i]))
		}

	}
	deltaD = abs(D2 - D1)

	# loop until deltaD is smaller than 1e-06 or the given D
	j = 1
	while ((deltaD >= D) & (j <= it.max)) {
		beta.t1 = beta.t2
		gamma.t1 = 1/(1 + exp(-X %*% beta.t1))
		W = matrix(0, sample.size, sample.size)
		for (i in 1:sample.size) {
			W[i, i] = gamma.t1[i] %*% (1 - gamma.t1[i])
		}

		beta.t2 = beta.t1 + solve(t(X) %*% W %*% X) %*% t(X) %*% (Y - gamma.t1)
		gamma.t2 = 1/(1 + exp(-X %*% beta.t2))

		D1 = D2 = 0
		for (i in 1:sample.size) {
			if (Y[i]) {
				D1 = D1 + 2 * Y[i] * log(Y[i]/gamma.t1[i])
				D2 = D2 + 2 * Y[i] * log(Y[i]/gamma.t2[i])
			} else if (1 - Y[i]) {
				D1 = D1 + 2 * (1 - Y[i]) * log((1 - Y[i])/(1 - gamma.t1[i]))
				D2 = D2 + 2 * (1 - Y[i]) * log((1 - Y[i])/(1 - gamma.t2[i]))
			}

		}
		deltaD = abs(D2 - D1)
		j = j + 1
	}
	beta.MLE = beta.t2

	return(beta.MLE)
}


#::::::::::::::::::::LRT function::::::::::::::::::::

logistic.LRT.minus <- function(Y, Xa, Xd, beta.0 = c(0, 0, 0), D = 1e-06, it.max = 100) {
	LRT.ls = list()
	beta1.hat.vec = c()

	sample.size = length(Y)
	ones = rep(1, sample.size)
	X = cbind(ones, Xa, Xd)

	beta0.hat = c(mean(Y), 0, 0)
	gamma0 = 1/(1 + exp(-X %*% beta0.hat))

	beta1.hat = logistic.IRLS.minus(Y, X, beta.0, D, it.max)
	beta1.hat.vec = c(beta1.hat.vec, beta1.hat)
	gamma1 = 1/(1 + exp(-X %*% beta1.hat))

	#LRT
	l.beta0.hat = 0
	l.beta1.hat = 0
	for (i in 1:length(Y)) {
		l.beta0.hat = l.beta0.hat + (Y[i] * log(gamma0[i]) + (1 - Y[i]) * log(1 - gamma0[i]))
		l.beta1.hat = l.beta1.hat + (Y[i] * log(gamma1[i]) + (1 - Y[i]) * log(1 - gamma1[i]))
	}
	LRT = 2 * (l.beta1.hat - l.beta0.hat)

	#p value
	LRT.ls$pval = pchisq(LRT, 2, lower.tail = FALSE, log.p = FALSE)
	LRT.ls$beta.MLE = beta1.hat.vec

	return(LRT.ls)
}



################################## Control Function #################################

control.fun <- function() {

	result.ls = list() # store all the results together in one list! Please refer to the top of this file for a summary and the structure of result.ls.

	#:::::::::::::::::::: Import, filter and code genotype ::::::::::::::::::::
	
	#--------------Import sample--------------
	sample.ls = list()
	sample.ls = import.sample.fun()


	#--------------Remove genotype with MAF <0.05--------------
	sample.ls$filtered.geno = MAF.fun(sample.ls$g)


	#--------------Code Xa and Xd--------------
	sample.ls$Xa = Xa.coding.fun(sample.ls$filtered.geno)
	sample.ls$Xd = 1 - 2 * abs(sample.ls$Xa)
	size = dim(sample.ls$Xa)
	sample.size = size[1]
	num.sample = size[2]
	

	#--------------PCA and Xz--------------
	geno.pc = PCA.fun(sample.ls$Xa)

	# Create covariate matrix Xz
	Xz = NULL
	threshold.comp1 <- -0.03 #X-axis value which appears to split group 0 from group 1
	group0 <- which(geno.pc$loadings[, 1] < threshold.comp1)
	Xz <- array(1, sample.size)
	Xz[group0] <- 0
	num.group0 = sum(Xz == 0) #147
	num.group1 = sum(Xz == 1) # 352
	cat("There are", num.group0, "individuals in one group and", num.group1, "individuals in the other one.\n\n")

	sample.ls$Xz = Xz



	#:::::::::::::::::::: GWAS with Covariance ::::::::::::::::::::
	cat("::: GWAS with covariance :::\n")
	
	#--------------Calculate p-value--------------
	cat("Q5: Estimating beta parameters and calculating p values...\n")
	LRT.ls = list()
	beta.MLE.matrix = matrix(0, num.sample, 4)
	for (i in 1:num.sample) {
		logistic.LRT.ls <- logistic.LRT(sample.ls$y, sample.ls$Xa[, i], sample.ls$Xd[, i], sample.ls$Xz, beta.0 = c(0, 
			0, 0, 0), D = 1e-06, it.max = 100)
		LRT.ls$pval[i] = logistic.LRT.ls$pval
		LRT.ls$log10pval[i] = -log10(logistic.LRT.ls$pval)
		beta.MLE.matrix[i, ] = logistic.LRT.ls$beta.MLE
	}

	covar = list()
	covar$beta.MLE = beta.MLE.matrix
	covar$beta.a.mean = mean(abs(beta.MLE.matrix[, 2]))
	covar$beta.d.mean = mean(abs(beta.MLE.matrix[, 3]))
	covar$beta.z.mean = mean(abs(beta.MLE.matrix[, 4]))
	covar$pval = LRT.ls$pval
	result.ls$covar = covar
	result.ls$s0 = "############################################"
	cat("beta parameters and p values are stored in sample.ls$covar\n")
	cat("The mean of the absolute values of beta.a.hat is", result.ls$covar$beta.a.mean,"\n")
	cat("The mean of the absolute values of beta.d.hat is", result.ls$covar$beta.d.mean,"\n")
	cat("The mean of the absolute values of beta.z.hat is", result.ls$covar$beta.z.mean,"\n\n")
	
	#--------------Histogram of p value--------------
	
	cat("Q6: Histogram of p values.\n")
	png("Q6.1_Histogram of p value (covariance).png", width = 500, height = 500)
	hist(LRT.ls$pval, main = "Histogram of p value (covariance)", xlab = "p value")
	dev.off()
	cat("\"Q6.1_Histogram of p value (covariance).png\" created.\n\n")


	#--------------QQ plot--------------
	
	cat("Q7: QQ plot.\n")
	log10pval.expected = -log10((1:ncol(sample.ls$Xa))/ncol(sample.ls$Xa))
	ranklog10pval.expected = log10pval.expected[order(log10pval.expected)]

	log10pval.covar = -log10(covar$pval)
	ranklog10pval.covar = log10pval.covar[order(log10pval.covar)]

	png("Q7.1_QQ Plot with covariate.png", width = 500, height = 500)
	plot(ranklog10pval.expected, ranklog10pval.covar, main = "QQ Plot with covariate", xlab = "Expected", ylab = "pval without covariate")
	dev.off()
	cat("\"Q7.1_QQ Plot with covariate.png\" created.\n\n")


	#--------------Manhattan plot with Bonferroni correction--------------
	
	cat("Q8: Manhattan plot.\n")
	position = seq(1, num.sample)
	y.bonferroni.0.05 = rep(-log10(0.05/num.sample), num.sample)
	y.bonferroni.0.1 = rep(-log10(0.1/num.sample), num.sample)
	png("Q8_Manhattan Plot with Bonferroni correction (covariance).png", width = 500, height = 500)
	plot(position, LRT.ls$log10pval, main = "Manhattan Plot with Bonferroni Correction (covariance)", ylab = "-log(p)", 
		xlab = "position in genotype markers")
	lines(position, y.bonferroni.0.05, type = "l", col = "red")
	lines(position, y.bonferroni.0.1, type = "l", col = "blue")
	dev.off()
	cat("\"Q8_Manhattan Plot with Bonferroni correction (covariance).png\" created.\n")

	cat(sum(LRT.ls$log10pval > y.bonferroni.0.05), "hits are indicatd at a (single test) alpha=0.05.\n")
	cat(sum(LRT.ls$log10pval > y.bonferroni.0.1), "hits are indicatd at a (single test) alpha=0.1.\n\n")



	#:::::::::::::::::::: GWAS without Covariance ::::::::::::::::::::
	cat("::: GWAS without covariance :::\n")
	
	#--------------Calculate p-value--------------
	cat("Q5: Estimating beta parameters and calculating p values...\n")
	LRT.minus.ls = list()
	beta.minus.MLE.matrix = matrix(0, num.sample, 3)
	for (i in 1:num.sample) {
		logistic.LRT.minus.ls <- logistic.LRT.minus(sample.ls$y, sample.ls$Xa[, i], sample.ls$Xd[, i], beta.0 = c(0, 
			0, 0), D = 1e-06, it.max = 100)
		LRT.minus.ls$pval[i] = logistic.LRT.minus.ls$pval
		LRT.minus.ls$log10pval[i] = -log10(logistic.LRT.minus.ls$pval)
		beta.minus.MLE.matrix[i, ] = logistic.LRT.minus.ls$beta.MLE
	}
	nocovar = list()
	nocovar$beta.minus.MLE = beta.minus.MLE.matrix
	nocovar$beta.a.minus.mean = mean(abs(beta.minus.MLE.matrix[, 2]))
	nocovar$beta.d.minus.mean = mean(abs(beta.minus.MLE.matrix[, 3]))
	nocovar$pval = LRT.minus.ls$pval
	result.ls$nocovar = nocovar
	result.ls$s1 = "############################################"
	cat("beta parameters and p values are stored in sample.ls$covar\n")
	cat("The mean of the absolute values of beta.a.minus.hat is", result.ls$nocovar$beta.a.minus.mean,"\n")
	cat("The mean of the absolute values of beta.d.minus.hat is", result.ls$nocovar$beta.d.minus.mean,"\n\n")


	#--------------Histogram of p value--------------	
	
	cat("Q6: Histogram of p values.\n")
	png("Q6.2_Histogram of p value (no covariance).png", width = 500, height = 500)
	hist(LRT.minus.ls$pval, main = "Histogram of p value (no covariance)", xlab = "p value")
	dev.off()
	cat("\"Q6.2_Histogram of p value (no covariance).png\" created.\n\n")
	

	#--------------QQ plot--------------
	
	cat("Q7: QQ plot.\n")
	log10pval.expected = -log10((1:ncol(sample.ls$Xa))/ncol(sample.ls$Xa))
	ranklog10pval.expected = log10pval.expected[order(log10pval.expected)]

	log10pval.nocovar = -log10(nocovar$pval)
	ranklog10pval.nocovar = log10pval.nocovar[order(log10pval.nocovar)]


	png("Q7.2_QQ Plot without covariate.png", width = 500, height = 500)
	plot(ranklog10pval.expected, ranklog10pval.nocovar, main = "QQ Plot without covariate", xlab = "Expected", ylab = "pval without covariate")
	dev.off()
	cat("\"Q7.2_QQ Plot without covariate.png\" created.\n\n")


	#:::::::::::::::::::: finding hits ::::::::::::::::::::
	
	cat("Q9: Find the two most significant hits indicating different locations in the genome.\n")
		
	#--------------List of hits--------------
	cat("Finding potential hits...\n")
	alpha = 0.05
	bonferroni = alpha/num.sample
	bonferroni.log10 = -log10(bonferroni)

	hits.ls = list()
	h = 1
	for (p in 1:num.sample) {
		if (-log10(LRT.ls$pval[p]) > bonferroni.log10) {
			hits.ls$positon[h] = p
			hits.ls$pval[h] = -log(LRT.ls$pval[p])
			h = h + 1
		}
	}

	result.ls$hits.ls = hits.ls
	result.ls$s2 = "############################################"


	#--------------Cluster hits into sets--------------
	
	# measure distance between two adjacent hits.
	# set a threshold for clustering by calculating the mean distance. (this method depends on the following assumptions: 1.there should not be too many hits AND 2.the distance between two hits from two sets is much larger than the distance from two hits within a set)
	cat("Clustering hits into sets...\n")
	hit.position = c(hits.ls$positon)
	hit.dist = hit.position[-1] - hit.position[-length(hit.position)] # calculate the distance between two adjacent hits
	threshold = mean(hit.dist) # calculate the threshold of distance for clustering

	hit.set.ls = list() # list to store different sets of hits
	k = 1 # index for hit.set.ls
	hit.set = c() # vector to store the position of hits in the same set
	h = 1 # index for hit within a set
	hit.set[1] = hit.position[h]
	h = h + 1
	while (h <= length(hit.position)) {
		while (((hit.position[h] - hit.position[h - 1])) < threshold & (h <= length(hit.position))) {
			hit.set = c(hit.set, hit.position[h])
			h = h + 1
		}
		hit.set.ls[[k]] = hit.set
		k = k + 1
		hit.set = c() # empty his.set
		hit.set[1] = hit.position[h] # sign the current hit.position as the first element in the new hit.set
		h = h + 1
	}

	result.ls$hit.set.ls = hit.set.ls
	result.ls$s3 = "############################################"


	#--------------Find the most significant marker in each sets--------------
	cat("Finding the most significant hits in each sets...\n")
	num.set = length(hit.set.ls) # the number of distinct sets
	sig.hit.position = c()
	sig.hit.pval = c()
	sig.hit.MAF = c()
	for (i in 1:num.set) {
		# find the position with the most significant(smallest) p value in a set
		sig.pval.index = match(min(LRT.ls$pval[hit.set.ls[[i]]]), LRT.ls$pval[hit.set.ls[[i]]]) # index in hit.set.ls[[i]]
		sig.pval.position = hit.set.ls[[i]][sig.pval.index] # real marker position in genome

		sig.hit.position = c(sig.hit.position, sig.pval.position) # store the position of most sig p value from each set
		sig.hit.pval = c(sig.hit.pval, LRT.ls$pval[sig.pval.position]) # store the most sig p value from each set

		unique.allele = unique(sample.ls$filtere.geno[, sig.pval.position]) # find the unique allele
		allele.freq = c(sum(sample.ls$filtere.geno[, sig.pval.position] == unique.allele[1]), sum(sample.ls$filtere.geno[, sig.pval.position] == 
			unique.allele[2])) # calculate allele frequence
	}
	sig.hit.ls = list()
	sig.hit.ls$position = sig.hit.position
	sig.hit.ls$pval = sig.hit.pval

	result.ls$sig.hit.ls = sig.hit.ls
	result.ls$s4 = "############################################"
	cat("The positions and p values of the two most significant hits from different locations are stored in result.ls$sig.hit.ls\n")
	cat("Position:", result.ls$sig.hit.ls$position,"\n")
	cat("Corresponding p values:", result.ls$sig.hit.ls$pval,"\n\n")


	#:::::::::::::::::::: Epistasis ::::::::::::::::::::
	cat("Q9: The nine dummy variables for the first five individuals for epistatic analysis.\n")
	position.1 = sig.hit.ls$position[1]
	position.2 = sig.hit.ls$position[2]

	epi.X = list()
	for (i in 1:5) {
		Xa.1 = sample.ls$Xa[i, position.1]
		Xd.1 = sample.ls$Xd[i, position.1]
		Xa.2 = sample.ls$Xa[i, position.2]
		Xd.2 = sample.ls$Xd[i, position.2]
		Xa.1_Xa.2 = Xa.1 * Xa.2
		Xa.1_Xd.2 = Xa.1 * Xd.2
		Xd.1_Xa.2 = Xd.1 * Xa.2
		Xd.1_Xd.2 = Xd.1 * Xd.2
		epi.X[[i]] = c(1, Xa.1, Xd.1, Xa.2, Xd.2, Xa.1_Xa.2, Xa.1_Xd.2, Xd.1_Xa.2, Xd.1_Xd.2)
	}
	result.ls$epi.X = epi.X
	cat("The dummy variables are stored in result.ls$epi.X\n")
	for (i in 1:5) {
		cat("Individual",i,":", result.ls$epi.X[[i]],"\n")
	}
	cat("\nFinalScript.r is done.")
	return(result.ls)

}
