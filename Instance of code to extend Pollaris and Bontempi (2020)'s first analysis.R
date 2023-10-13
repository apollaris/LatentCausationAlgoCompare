#Main author : A. Pollaris
#Date : October 13, 2023
#Concern: implementations in the R language to generate data and perform causal inference. This code can help providing results to compare 3 causal inference algorithms:

#1. The factor scores computation followed by the Direction Dependence Analysis independence component (FSC+DDA)
#2. The Latent Causation Algorithm (LC)
#3. The Hypothesis Test Algorithm (HT)

#Welcome in this instance of code implementing causal inference algorithms for latent pairs of variables, tools for performing simulations and still other features. 
#For every statistical analysis, the author did his best to make sure this code is working fine. However, if some problems or inconsistencies are found using the code below,
#please inform the author, A. Pollaris (Arnaud.Pollaris@ulb.ac.be), to help him to improve this work.  

#By using this code, you accept to take all responsibility for any negative consequences that might result. This code is provided "as is", without warranty of any kind,
#express or implied. In no event shall the author or copyright holders be liable for any claim, damages or other liability.

#N.B.: 
#1. Some parts of the code below come from the DDA project (available online : https://www.ddaproject.com/ ; see also Wolfgang Wiedermann, Xintong Li and Alexander von Eye)
#   and were not written by A. Pollaris.
#2. Some functions below were already available in https://github.com/apollaris/LatentCausation
#3. Note that some parts of the current implementation are still in a beta version. Parts of the code (as for example the function "causal.interpret()"
#   used with only one specified indicator for the "assumed_cause_indicators" parameter, for the "assumed_consequence_indicators" parameter or even for both) are not validated yet.

###############################################################################################################################################################################
###############################################################################################################################################################################

#Imports:

library(FactoMineR)
library(energy)
library(dHSIC)

############################################################################################
############################################################################################

## ----------------------------------------------------------------------------------------------------------------------- ##
## File:    boot_hsic_test.R 
## Purpose: Computes Sen & Sen's (2014) HSIC test. This is essentially a wrapper for the source code provided 
##          by Sen & Sen (2014). Source code available from http://www.stat.columbia.edu/~bodhi/Bodhi/Publications.html
## Author:  WW / Bodhisattva Sen 
## Date:    Dec 7, 2018 
## ----------------------------------------------------------------------------------------------------------------------- ##


hsic.test <- function(model, x = NULL, hx = 1, hy = 1, B = 1000, parallelize = FALSE, cores = 2)
{

## --------------------------------------------------------------------------------------------------
## model:       An "lm" object representing the model to be tested. 
## x:           A vector specifying the target predictor that is used to test independence. 
## hx:          Numeric value specifying the bandwidth parameter for "x". 
## hy:          Numeric value specifying the bandwidth parameter for "y". 
## B:           Numeric value specifying the number of resamples used to compute the HSIC p-value. 
## parallelize: A logical value indicating whether boostrapping is performed on multiple cores. 
## cores:       A numeric value indicating the number of cores. Only used if parallelize = TRUE.
## Example:     m <- lm(y ~ x + z) 
##              hsic.test(m, x, B = 500, parallelize = TRUE, cores = 4)     
## -------------------------------------------------------------------------------------------------- 

    X <- model.matrix(model)
	b <- as.matrix(coef(model))
	n <- nobs(model) 
	e <- resid(model)
	
	if(is.null(x)) stop("Predictor 'x' is missing.")
	if(!is.matrix(x)) x <- as.matrix(x)

	T_hat <- HSIC(x, e, 2*hx*hx, 2*hy*hy)   ## Compute the test-statistic

    # Implementing the bootstrap procedure (parallelized if specified)
    
	e_0 <- e - mean(e)      ## centered residuals

    if(parallelize){

      require(doParallel)
	  cl <- makeCluster(cores)
      registerDoParallel(cl)

         T_hat_B <- foreach(icount(B), .combine = c, .export = "HSIC") %dopar% {

               idx <- sample(n, n, replace=TRUE)    ## with replacement samples from the errors
        	   e_B <- e_0[idx]
    
        	   idx2 <- sample(n, n, replace=TRUE)   ## with replacement samples from the predictors
        	   x_B <- x[idx2,]
        	   X_B <- X[idx2,]       

        	   yhat_B  <- X_B %*% b + e_B                              ## Create the bootstrap response values
        	   bhat_B  <- solve(t(X_B) %*% X_B) %*% t(X_B) %*% yhat_B  ## Get new slope estimates 
        	   e_hat_B <- yhat_B - X_B %*% bhat_B                      ## Bootstrap residuals

        	   hx_B <- hx; hy_B <- hy
        	   HSIC(x_B, e_hat_B, 2*hx_B*hx_B, 2*hy_B*hy_B)
               }
    
	   stopCluster(cl)  
    
	} else {

      T_hat_B <- matrix(0, B, 1)
      for(j in 1:B){
      	    idx <- sample(n, n, replace=TRUE)                ## with replacement samples from the errors
        	e_B <- e_0[idx]
    
        	idx2 <- sample(n, n, replace=TRUE)               ## with replacement samples from the predictors
        	x_B <- x[idx2,]
        	X_B <- X[idx2,]       

        	yhat_B  <- X_B %*% b + e_B                               ## Create the bootstrap response values
        	bhat_B  <- solve(t(X_B) %*% X_B) %*% t(X_B) %*% yhat_B   ## Get new slope estimates 
        	e_hat_B <- yhat_B - X_B %*% bhat_B                       ## Bootstrap residuals

        	hx_B <- hx; hy_B <- hy
        	T_hat_B[j] <- HSIC(x_B, e_hat_B, 2*hx_B*hx_B, 2*hy_B*hy_B)

 	    }
	}
	
	pval <- mean(T_hat_B >= T_hat)
	names(T_hat) <- "HSIC"
	#hist(T_hat_B, col = "grey"); abline(v = T_hat)
	output <- list(statistic = T_hat, p.value = pval, 
	               method = "Hilbert-Schmidt Independence Test", data.name = deparse(formula(model)))
	class(output) <- "htest"
	return(output)
}

## --- Helper function to calculate HSIC value 

HSIC <- function(x, y, hx, hy)
{
   n <- length(y)
   if(n != length(x)) stop("Variables must have same length.")   
   H <- diag(n) - matrix(1, n, n)/n    
   K <- exp(-as.matrix(dist(x, method = "euclidean", diag = TRUE, upper = TRUE))^2/hx)
   L <- exp(-as.matrix(dist(y, method = "euclidean", diag = TRUE, upper = TRUE))^2/hy)
   hsic <- sum(diag(K %*% H %*% L %*% H))/n
   return(hsic)
}



############################################################################################
############################################################################################

#Auxiliary functions to generate datasets:

sd_resid_std=function(beta1,beta2=0,corr=0){
	# in a regression model (Y~beta1*X1+beta2*X2+err), compute the theoretical standard deviation of the disturbance "err" (noise term) 
	# such that the total standard deviation of the target variable (called here Y) is equal to 1 (with a mean equal to 0).
	# beta1 is mandatory: it is the theoretical regression coefficient related to a given explanatory variable X1.
	# beta2 is optional :  it is the theoretical regression coefficient related to a other given explanatory variable X2.
	# corr is the theoretical Pearson's correlation between X1 and X2. By default, explanatory variables are uncorrelated. 
	# The numerical value of corr has no influence if beta2 is equal to 0. 
	# return sd_err : the theoretical standard deviation of the disturbance.

	sd_err=sqrt(1-beta1^2-beta2^2-2*beta1*beta2*corr)
	return(sd_err)
}

runif_sd=function(n, sd=1){
	#generate n samples from a uniform distribution with mean equal to 0 and standard deviation equal to sd.
	# by default, sd is equal to 1 (i.e., standardized uniform distribution)
	rand_nbrs=runif(n, min = -0.5, max = 0.5)/sqrt(1/12)*sd
	return(rand_nbrs)
}

rchisq_sd=function(n, sd=1, df=1){
	#generate n samples from a chi-square distribution with degree of freedom equal to df (by default df=1), mean equal to 0 and standard deviation equal to sd.
	# by default, sd is equal to 1 (i.e., standardized chi-square distribution)
	rand_nbrs=(rchisq(n, df=df)-df)/sqrt(2*df)*sd
	return(rand_nbrs)
}

causal_sim=function(model, n, n_samples, obs_var="all", seed=4321){
	# model is a string like:'
	#	#Exogenous variables (come first)
	#	U ~ unif
	#	#Latent Endogeneous variables
	#	X ~ 0.2*U + chisq
	#	Y ~ 0.3*U + 0.5*X + unif
	#
	#	#Observed indicators
	#	
	#	x1 ~ 0.7*X + norm
	#	x2 ~ 0.8*X + norm
	#	x3 ~ 0.9*X + norm
	#	
	#	y1 ~ 0.7*Y + norm
	#	y2 ~ 0.8*Y + norm
	#	y3 ~ 0.9*Y + norm	
	#'
	
	#n is the sample ssize for each sample to simulate
	#n_samples is the number of samples to generate
	#obs_var indicates for which variables to keep records in the generated dataframes
	#seed: to indicate a seed at the beginning of the simulations

	# Returns a list of datasets "datasets_list"

	if (is.numeric(seed))	set.seed(seed)

	#pretreatement:
	prep_model_lines=gsub(" ", "", model)
	prep_model_lines=gsub("\t", "", prep_model_lines)
	prep_model_lines=unlist(strsplit(prep_model_lines,"\n",fixed=TRUE))
	model_lines=c()
	for (i in 1:length(prep_model_lines)){
		if (prep_model_lines[i]!="" & substr(prep_model_lines[i],1,1)!="#") {
			if (length(grep("#", prep_model_lines[i]))){
				model_lines=c(model_lines, unlist(strsplit(prep_model_lines[i],"#", fixed=TRUE))[1])
			}
			else {
				model_lines=c(model_lines, prep_model_lines[i])
			}	
		}
	}

	coefficients=matrix(nrow=1,ncol=0)
	resid_shapes=matrix(nrow=1,ncol=0)
	predict_part=matrix(nrow=1,ncol=0)

	for (i in 1:length(model_lines)){
		if (0==length(grep("~", model_lines[i]))) {
			stop(paste("There is no '~' in line '", model_lines[i],"'"))
		}
		if (substr(model_lines[i],1,1)=="~" | substr(model_lines[i],nchar(model_lines[i]),nchar(model_lines[i]))=="~") {
			stop(paste("The character '~' in line '", model_lines[i],"' is in a bad position !"))
		}

		VD_VI=unlist(strsplit(model_lines[i],"~",fixed=TRUE))
		VD_name=VD_VI[1]
		VI=VD_VI[2]

		if (substr(VI,1,1)=="+" | substr(VI,nchar(model_lines[i]),nchar(model_lines[i]))=="+") {
			stop(paste("A character '+' in line '", model_lines[i],"' is in a bad position !"))
		}

		VI=unlist(strsplit(VI,"+",fixed=TRUE))
		resid_shape=VI[length(VI)]

		if (resid_shape!="norm" & resid_shape!="unif" & resid_shape!="chisq"){
			stop(paste("The distribution of the residuals in '", model_lines[i],"' was neither 'norm' nor 'unif' nor 'chisq' !"))
		}

		resid_shapes=cbind(resid_shapes, matrix(resid_shape, nrow=1))
		colnames(resid_shapes)[dim(resid_shapes)[2]]=VD_name

		VI=VI[-length(VI)]

		if (length(VI)>0){

			element_predict_part=""

			for (j in 1:length(VI)){
				coef_variable=unlist(strsplit(VI[j],"*",fixed=TRUE))				
				if (length(coef_variable)!=2){
					stop(paste("The term '", VI[j],"' does not contains a '*' or can not be splitted in exactly 2 elements (one coefficient and one explanatory variable)."))					
				}

				variable=coef_variable[2]
				if (is.na(variable) | !is.na(suppressWarnings(as.numeric(variable))) ){
					stop(paste("The second factor in '", coef_variable,"' is not an acceptable variable name !"))					
				}
				coef_name=paste("beta",VD_name,variable,sep="_")
				coef_value=as.numeric(coef_variable[1])
				if (is.na(coef_value)){
					stop(paste("The first factor in '", coef_variable,"' is not a numerical coefficient !"))
				}
				if (1<abs(coef_value)){
					stop(paste("The coefficient in '", coef_variable,"' is not between -1 and 1 (included) !"))
				}

				coefficients=cbind(coefficients, matrix(coef_value, nrow=1))
				colnames(coefficients)[dim(coefficients)[2]]=coef_name

				coef_predict_part_variable=paste(coef_value , "*data_mat[,\"" , variable,"\"]",sep="")
				element_predict_part=paste(element_predict_part,coef_predict_part_variable,sep="+")

			}

			predict_part=cbind(predict_part, matrix(element_predict_part, nrow=1))
			colnames(predict_part)[dim(predict_part)[2]]=VD_name

		}
		
	}

	#Data Generation:

	#Use:
	#coefficients
	#resid_shapes
	#predict_part

	datasets_list=list()

	for (k in 1:n_samples) {
	
		data_mat=matrix(nrow=n,ncol=dim(resid_shapes)[2])
		colnames(data_mat)=colnames(resid_shapes)
	
		if ("all"==obs_var[1]) {
			obs_var=colnames(data_mat)
		}
	
		var_index_pred=1
	
		for (var_index in 1:(dim(resid_shapes)[2])){
	
			current_VD_name=colnames(resid_shapes)[var_index]
			coeff_current_VD_name=paste("beta_",current_VD_name, sep="")
	
			#identification of the relevant beta values+ if needed the correlation between the 2 related explained variables
	
			coef_relevant_colnum=grep(coeff_current_VD_name,colnames(coefficients))
	
			if (length(coef_relevant_colnum)>2){
				stop(paste("The function causal_sim does not support an explained variable '", current_VD_name,"' with more than 2 explanatory variables."))
			}
	
			if (length(coef_relevant_colnum)==2){
	
				beta1=coefficients[,coef_relevant_colnum[1]]
				beta2=coefficients[,coef_relevant_colnum[2]]
	
				coef_relevant_colnames=colnames(coefficients)[coef_relevant_colnum]
				coef_relevant_colnames_splitted=strsplit(coef_relevant_colnames,"_", fixed=TRUE)
	
				v_expl1=coef_relevant_colnames_splitted[[1]][[3]]
				v_expl2=coef_relevant_colnames_splitted[[2]][[3]]
				test=c(paste("beta_", v_expl1,"_",v_expl2, sep=""),paste("beta_", v_expl2,"_",v_expl1, sep=""))
	
				corr_v_expl1_2=coefficients[, which( colnames(coefficients)==test[1] | colnames(coefficients)==test[2] ) ]
	
			}
			if (length(coef_relevant_colnum)==1){
				beta1=coefficients[,coef_relevant_colnum[1]]
				beta2=0
				corr_v_expl1_2=0
			}
	
	
			if (current_VD_name==colnames(predict_part)[var_index_pred]){#if the variable is not exogeneous
	
				if (resid_shapes[var_index]=="norm"){
					#eval(parse(text=predict_part[1]))
					data_mat[,current_VD_name]=eval(parse(text=predict_part[var_index_pred])) + rnorm(n, sd=sd_resid_std(beta1=beta1,beta2=beta2,corr=corr_v_expl1_2))
				}
				else {
					if (resid_shapes[var_index]=="unif"){
						data_mat[,current_VD_name]=eval(parse(text=predict_part[var_index_pred])) + runif_sd(n, sd=sd_resid_std(beta1=beta1,beta2=beta2,corr=corr_v_expl1_2))
					}
					else { #if "chisq"
						data_mat[,current_VD_name]=eval(parse(text=predict_part[var_index_pred])) + rchisq_sd(n, sd=sd_resid_std(beta1=beta1,beta2=beta2,corr=corr_v_expl1_2))
					}
				}
	
	
				var_index_pred=var_index_pred+1
			} else { #variable exogeneous
	
				if (resid_shapes[var_index]=="norm"){
					data_mat[,current_VD_name]=rnorm(n, sd=1)
				}
				else {
					if (resid_shapes[var_index]=="unif"){
						data_mat[,current_VD_name]=runif_sd(n, sd=1)
					}
					else { #if "chisq"
						data_mat[,current_VD_name]=rchisq_sd(n, sd=1)
					}
				}
				
			}
	
		}
		datasets_list[[k]]=as.data.frame(data_mat[,obs_var])
	}

	return(datasets_list)

}


############################################################################################
############################################################################################

causation_stat=function(data, assumed_cause_indicators, assumed_consequence_indicators, indep.measures="all", dcor.pvalues=TRUE, HSIC.pvalues=c("gamma","bootstrap"), dcor.R=1000, HSIC.B=500){

	#The goal of this function is providing some statistics helping for the causal inference. After a principal component analysis for each group of indicators (where only the first factor is kept for each PCA),
	# 2 resgressions are made with the 2 factors (one using the first factor as the predictor and the second one as the dependant variable and vice evrsa for the second regression model).
	# Then, independency statistics are computed for each regression model between the residuals and the predictor. Available methods are "spearman", "dcor", "hsic".
   
	#data is a data.frame
	#"assumed_cause_indicators" is a vector of characters giving the names of the indicators for the assumed (latent) cause.
	#"assumed_consequence_indicators" is a vector of characters giving the names of the indicators for the assumed (latent) consequence.
	#indep.measures: a vector containing either : "spearman", "dcor", "hsic" or "all" in order to specify for which method relevant causation statistics are required.
	#dcor.pvalues: logical, should p.values for dcor tests be computed? (it is a permutation bootstrap, with dcor.R replicates)
	#HSIC.pvalues: a vector indicating which method to use to compute p.values for hsic tests. Choices are : "gamma","bootstrap" (or both or "none"). Bootstrap calculus is made according Sen & Sen's (2014) HSIC test and can take much time to run.
	#dcor.R is the number of replicates and is only used when dcor.pvalues is TRUE
	#HSIC.B number of resamples used to compute the bootstrap's HSIC p-value. Only used when HSIC.pvalues contains "bootstrap".

	#returns "output.indep.measures", a list() referencing all the asked indep.measures and their statistics. 


	#Performs a PCA for each latent variable (assumed cause Xi and assumed consequence Eta):
	resXi=PCA(data[,assumed_cause_indicators], graph=FALSE)
	Fx=resXi$ind$coord[,1]
	resEta=PCA(data[,assumed_consequence_indicators], graph=FALSE)
	Fy=resEta$ind$coord[,1]

	#Performs a linear regression Fy~Fx and storage of the residuals:
	modelFyFx=lm(Fy~Fx)
	residFy=residuals(modelFyFx)
	#Performs a linear regression Fx~Fy and storage of the residuals:
	modelFxFy=lm(Fx~Fy)
	residFx=residuals(modelFxFy)

	output.indep.measures=list() #empty list referencing all the asked indep.measures 
	
	if ( (sum(indep.measures=="spearman")>0) | (sum(indep.measures=="all")>0) ) {

		#computing the Spearman correlation between the residuals and the factors scores of the "predictor"
		# for each possible causal direction:
		rs.y=cor(Fx, residFy, method="spearman")
		rs.x=cor(Fy, residFx, method="spearman")

		#computing a difference score between both statistics (a positive difference is assumed to mean Xi->Eta)
		diff.rs=abs(rs.x)-abs(rs.y)

		rs.mat=matrix(c(rs.x,rs.y,diff.rs),nrow=1)
		colnames(rs.mat)=c("rs.x","rs.y","diff.rs")
		output.indep.measures=c(output.indep.measures, list("spearman"=rs.mat))

	}

	if ( (sum(indep.measures=="dcor")>0) | (sum(indep.measures=="all")>0) ) {

		if (dcor.pvalues) {
	
			#computing the Distance correlation between the residuals and the factors scores of the "predictor"
			# for each possible causal direction:
	
			dcor.y=dcor.test(Fx, residFy, index = 1.0, R=dcor.R)
			dcor.x=dcor.test(Fy, residFx, index = 1.0, R=dcor.R)
			
	
			#computing a difference score between both statistics (a positive difference is assumed to mean Xi->Eta)
			diff.dcor=dcor.x$statistic-dcor.y$statistic
	
			#computing a difference score between both p-values (a negative difference is assumed to mean Xi->Eta)
			diff.dcor.pval=dcor.x$p.value-dcor.y$p.value
	
			dcor.mat=matrix(c(dcor.x$statistic,dcor.x$p.value,dcor.y$statistic,dcor.y$p.value,diff.dcor,diff.dcor.pval),nrow=2,ncol=3)
			colnames(dcor.mat)=c("dcor.x","dcor.y","diff.dcor")
			rownames(dcor.mat)=c("statistic","p.value")
			output.indep.measures=c(output.indep.measures, list("dcor"=dcor.mat))
		}
		else {#only descriptives statistics:
		
			#computing the Distance correlation between the residuals and the factors scores of the "predictor"
			# for each possible causal direction:
	
			dcor.y=dcor(Fx, residFy)
			dcor.x=dcor(Fy, residFx)
	
			#computing a difference score between both statistics (a positive difference is assumed to mean Xi->Eta)
			diff.dcor=dcor.x-dcor.y
	
			dcor.mat=matrix(c(dcor.x,dcor.y,diff.dcor),nrow=1,ncol=3)
			colnames(dcor.mat)=c("dcor.x","dcor.y","diff.dcor")
			rownames(dcor.mat)=c("statistic")
			output.indep.measures=c(output.indep.measures, list("dcor"=dcor.mat))			
		
		}

	}

	if ( (sum(indep.measures=="hsic")>0) | (sum(indep.measures=="all")>0) ) {
		#computing the Hilbert-Schmidt Independence Criterion between the residuals and the factors scores of the "predictor"
		# for each possible causal direction:

		### Using "boot_hsic_test.R" written by : WW / Bodhisattva Sen 

		save_seed=.Random.seed

		hsic.gamma.y = dhsic.test(Fx, residFy, method = "gamma", kernel = "gaussian")	 
		hsic.gamma.x = dhsic.test(Fy, residFx, method = "gamma", kernel = "gaussian")

		.Random.seed<<-save_seed #assignation to a global variable ".Random.seed"

		#computing a difference score between both statistics (a positive difference is assumed to mean Xi->Eta)
		diff.hsic=abs(hsic.gamma.x$statistic)-abs(hsic.gamma.y$statistic)

		hsic.mat=matrix(c(hsic.gamma.x$statistic,hsic.gamma.y$statistic,diff.hsic),nrow=1,ncol=3)    
		colnames(hsic.mat)=c("hsic.x","hsic.y","diff.hsic")
		rownames(hsic.mat)="statistic"

		if ( (sum(HSIC.pvalues=="gamma")>0) ) {

			#computing a difference score between both p-values (a negative difference is assumed to mean Xi->Eta)
			diff.hsic.pval.gamma=hsic.gamma.x$p.value-hsic.gamma.y$p.value

			hsic.mat=rbind(hsic.mat, matrix(c(hsic.gamma.x$p.value,hsic.gamma.y$p.value,diff.hsic.pval.gamma),nrow=1,ncol=3) )
			rownames(hsic.mat)[dim(hsic.mat)[1]]="p.value - gamma"
		}

		if ( (sum(HSIC.pvalues=="bootstrap")>0) ) {

		     	bw.y <- hsic.gamma.y$bandwidth
			#print( HSIC(Fx, residFy, hx=2*bw.y[1]*bw.y[1], hy=2*bw.y[2]*bw.y[2]) )
			 hsic.y <- hsic.test(modelFyFx, Fx, hx = bw.y[1], hy = bw.y[2], B = HSIC.B, parallelize = FALSE, cores = 2)
	
		     	bw.x <- hsic.gamma.x$bandwidth
			#print( HSIC(Fy, residFx, hx=2*bw.x[1]*bw.x[1], hy=2*bw.x[2]*bw.x[2]) )
			 hsic.x <- hsic.test(modelFxFy, Fy, hx = bw.x[1], hy = bw.x[2], B = HSIC.B, parallelize = FALSE, cores = 2)

			#computing a difference score between both p-values (a negative difference is assumed to mean Xi->Eta)
			diff.hsic.pval.boot=hsic.x$p.value-hsic.y$p.value

			hsic.mat=rbind(hsic.mat, matrix(c(hsic.x$p.value,hsic.y$p.value,diff.hsic.pval.boot),nrow=1,ncol=3) )
			rownames(hsic.mat)[dim(hsic.mat)[1]]="p.value - bootstrap"
		}

		output.indep.measures=c(output.indep.measures, list("hsic"=hsic.mat))


	}
	return(output.indep.measures)
}

LC=function(data, assumed_cause_indicators, assumed_consequence_indicators, indep.measures="all", CI=c(0.025,0.975), boot.pval=FALSE, dcor.R=1000, HSIC.B=500, B=1000, seed=1234)
{
	#This function helps to compute some percentiles bootstrap confidence intervals for statistics produced in "causation_stat()". It can include "rs.x.stat", "rs.y.stat", "diff.rs", "dcor.x.stat", "dcor.y.stat", "diff.dcor.stat", "diff.dcor.pvalue", "hsic.x.stat", "hsic.y.stat", "diff.hsic.stat", "diff.hsic.gamma_pval".
	#data is a data.frame
	#"assumed_cause_indicators" is a vector of characters giving the names of the indicators for the assumed (latent) cause.
	#"assumed_consequence_indicators" is a vector of characters giving the names of the indicators for the assumed (latent) consequence.
	#indep.measures: a vector containing either : "spearman", "dcor", "hsic" or "all" in order to specify for which method relevant causation statistics are required.
	#CI is a vector (length equal to 2) giving the bounds for the percentiles confidence intervals to compute by bootstrap. Default is c(0.025,0.975).
	#boot.pval: logical, do we need additional bootstrap CI for the difference of dcor's p.values and for the difference of hsic's gamma_p.value. Default is FALSE (quicker).
	#dcor.R is the number of replicates used for estimating each dcor.pvalue
	#HSIC.B is the number of resamples used to compute each bootstrap's HSIC p-value.
	#B is the number of bootstrap samples to draw before computing the additional percentiles bootstrap confidence intervals.
	#seed is a value used to initialize the randoms generators before computing first statistics on the original dataset and before generating the bootstrapped samples (and computing the related statistics).

	#The return value "output" is a list giving [1] statistics from the original sample (uses "causation_stat()" function); 
	#	[2] statistics computed with "causation_stat()" for each bootstrap sample; 
	#	[3] information about the number of bootstrapped sample where some statistics could not be computed. 
	#	[4] the asked additional bootstrapped CI.  

	NB_MAX_STATISTICS=11 # rs.x.stat, rs.y.stat, diff.rs, dcor.x.stat, dcor.y.stat, diff.dcor.stat, diff.dcor.pvalue, hsic.x.stat, hsic.y.stat, diff.hsic.stat, diff.hsic.gamma_pval

	if (is.numeric(seed))	set.seed(seed)
	
	n=dim(data)[1]

	output=list()

	#Analysis of the original sample:
	original_sample=causation_stat(data=data, assumed_cause_indicators=assumed_cause_indicators, assumed_consequence_indicators=assumed_consequence_indicators, indep.measures=indep.measures, dcor.pvalues=TRUE, HSIC.pvalues=c("gamma","bootstrap"), dcor.R=dcor.R, HSIC.B=HSIC.B)
	output=c(output, list("original_sample"=original_sample))

	if (is.numeric(seed))	set.seed(seed)

	#creation of a logical vector giving the asked statistics to bootstrap 
	selected_stat=rep(FALSE,NB_MAX_STATISTICS)
	if ( (sum(indep.measures=="spearman")>0) | (sum(indep.measures=="all")>0) ) {
		selected_stat[1]=TRUE
		selected_stat[2]=TRUE
		selected_stat[3]=TRUE
	}
	if ( (sum(indep.measures=="dcor")>0) | (sum(indep.measures=="all")>0) ) {
		selected_stat[4]=TRUE
		selected_stat[5]=TRUE
		selected_stat[6]=TRUE
		selected_stat[7]=boot.pval
	}
	if ( (sum(indep.measures=="hsic")>0) | (sum(indep.measures=="all")>0) ) {
		selected_stat[8]=TRUE
		selected_stat[9]=TRUE
		selected_stat[10]=TRUE
		selected_stat[11]=boot.pval
	}

	#bootstrapped statistics:
	boot_res=matrix(nrow=NB_MAX_STATISTICS,ncol=B)
	rownames(boot_res)=c("rs.x.stat", "rs.y.stat", "diff.rs", "dcor.x.stat", "dcor.y.stat", "diff.dcor.stat", "diff.dcor.pvalue", "hsic.x.stat", "hsic.y.stat", "diff.hsic.stat", "diff.hsic.gamma_pval")
	boot_res=boot_res[selected_stat,]

	#sample selection:
	bootstrap.select=matrix( sample(x=1:n, size=(B*n),replace=TRUE) , ncol=B )

	#computation and storage (for each sample) of the bootrapped estimated statistics
	if (!boot.pval){
		for (i in 1:B){
			data.boot=data[bootstrap.select[,i],]
			boot_sample=causation_stat(data=data.boot,assumed_cause_indicators=assumed_cause_indicators, assumed_consequence_indicators=assumed_consequence_indicators, indep.measures=indep.measures, dcor.pvalues=FALSE, HSIC.pvalues=c("none"), dcor.R=dcor.R, HSIC.B=HSIC.B)
			boot_res[,i]=c(boot_sample$spearman[1,"rs.x"],boot_sample$spearman[1,"rs.y"],boot_sample$spearman[1,"diff.rs"],boot_sample$dcor["statistic","dcor.x"],boot_sample$dcor["statistic","dcor.y"],boot_sample$dcor["statistic","diff.dcor"],boot_sample$hsic["statistic","hsic.x"],boot_sample$hsic["statistic","hsic.y"],boot_sample$hsic["statistic","diff.hsic"])
		}
	}
	else {
		for (i in 1:B){
			data.boot=data[bootstrap.select[,i],]
			boot_sample=causation_stat(data=data.boot,assumed_cause_indicators=assumed_cause_indicators, assumed_consequence_indicators=assumed_consequence_indicators, indep.measures=indep.measures, dcor.pvalues=TRUE, HSIC.pvalues=c("gamma"), dcor.R=dcor.R, HSIC.B=HSIC.B)
			boot_res[,i]=c(boot_sample$spearman[1,"rs.x"],boot_sample$spearman[1,"rs.y"],boot_sample$spearman[1,"diff.rs"],boot_sample$dcor["statistic","dcor.x"],boot_sample$dcor["statistic","dcor.y"],boot_sample$dcor["statistic","diff.dcor"],boot_sample$dcor["p.value","diff.dcor"],boot_sample$hsic["statistic","hsic.x"],boot_sample$hsic["statistic","hsic.y"],boot_sample$hsic["statistic","diff.hsic"],boot_sample$hsic["p.value - gamma","diff.hsic"])
		}
	}

	nb_problematic_bootstrapped_samples=sum((apply(is.na(boot_res),2,sum))!=0)
	info=(paste("There is", nb_problematic_bootstrapped_samples, "bootstrapped samples where some statistics could not be computed in boot.results. Please check it before interpreting boot.CI outputs !" ))	

	boot.CI=t( apply(boot_res, 1, quantile, probs=CI, na.rm=TRUE))

	output=c(output, list("boot.results"=boot_res),list("boot.info"=info), list("boot.CI"=boot.CI))

	return(output)
}

 
############################################################################################
############################################################################################

causal.interpret=function(outputs,assumed_cause_indicators, assumed_consequence_indicators,round.digits=3){

	#The goal of this function is to help interpreting the results from outputs given by the LC function.
	# "assumed_cause_indicators" is a vector of characters containing the names of the indicators measuring the assumed cause
	# "assumed_consequence_indicators" is a vector of characters containing the names of the indicators measuring the assumed consequence
	# when "assumed_cause_indicators" contains only one variable name, it is assumed this indicator has no error of measurement (which means the observed variable and the latent concept are the same). The same is also true for the "assumed_consequence_indicators" vector.
	# when either the "assumed_cause_indicators" or "assumed_consequence_indicators" (or both) is only one observed variable, and when some relevant additional tests are significants, more information will be returned about possible existing latent confounders.

	# round.digit is the number of decimals to show for numbers in the outputs. It is not used for statistical decisions.

	#The return "output.interpret" is a list of available methods containing the results and the associated interpretation for each test.

	#Warning: the interpretations given by this function can be considered only if some hypotheses are met : e.g.: linearity of the relationships, non normality of the disturbance, ...
	#Furthermore, to interpret tests speaking about confounder(s), it is first assumed the conclusion related to the CI of the bootstrapped diff is correct (when 0 is not included in the CI interval related to the diff...).
	#It is also assumed there is no cyclic patterns in the true model and no bi-directional relationships. 

	if (2!=dim(outputs$boot.CI)[2]) {
		stop("CI intervals must be defined by exactly 2 bounds in outputs !")
	}

	CI_perc=as.numeric(gsub("%","", colnames(outputs$boot.CI)))
	alpha=(100-abs(CI_perc[2]-CI_perc[1]))/100

	output.interpret=list()

	#strings to prints (if needed) in output.interpret: 
	print_xi_arrow_eta=paste("(",paste(assumed_cause_indicators, collapse=", "),")", "-->", "(",paste(assumed_consequence_indicators, collapse=", "),")","is preferred")
	print_eta_arrow_xi=paste("(",paste(assumed_consequence_indicators, collapse=", "),")", "-->", "(",paste(assumed_cause_indicators, collapse=", "),")","is preferred")
	print_confounder=paste("And, a latent confounder U must be:","(",paste(assumed_cause_indicators, collapse=", "),")","<-- U -->","(",paste(assumed_consequence_indicators, collapse=", "),")")
	print_no_ccl="Based on this CI, data do not allow to show a causal direction"
	print_eta_arrow_xi_or_conf=paste("The hypothesis: <<","(",paste(assumed_cause_indicators, collapse=", "),")", "-->", "(",paste(assumed_consequence_indicators, collapse=", "),")","WITH no confounder >> might be rejected")#paste("(",paste(assumed_consequence_indicators, collapse=", "),")", "-->", "(",paste(assumed_cause_indicators, collapse=", "),")","is possible","AND/OR, a latent confounder U can be like:","(",paste(assumed_cause_indicators, collapse=", "),")","<-- U -->","(",paste(assumed_consequence_indicators, collapse=", "),")") 
	print_xi_arrow_eta_or_conf=paste("The hypothesis: <<","(",paste(assumed_consequence_indicators, collapse=", "),")", "-->", "(",paste(assumed_cause_indicators, collapse=", "),")","WITH no confounder >> might be rejected")#paste("(",paste(assumed_cause_indicators, collapse=", "),")", "-->", "(",paste(assumed_consequence_indicators, collapse=", "),")","is possible","AND/OR, a latent confounder U can be like:","(",paste(assumed_cause_indicators, collapse=", "),")","<-- U -->","(",paste(assumed_consequence_indicators, collapse=", "),")")

	if (sum(names(outputs$original_sample)=="spearman")){

		spearman_mat=matrix(nrow=1,ncol=2)
		colnames(spearman_mat)=c("results","interpretation")
		rownames(spearman_mat)="bootCI on diff"

		spearman_mat["bootCI on diff","results"]=paste(round(outputs$boot.CI["diff.rs",],round.digits), collapse=" ; ")

		# if CI diff.rs >0 : 	preference for xi -> eta
		if (outputs$boot.CI["diff.rs",1]>0 & outputs$boot.CI["diff.rs",2]>0) {
			spearman_mat["bootCI on diff","interpretation"]=print_xi_arrow_eta

			#Additionally, if the preferred cause is an obs_var and if CI rs.y sig. : add xi <-U-> eta
			if ((length(assumed_cause_indicators)==1) & ((outputs$boot.CI["rs.y.stat",1]*outputs$boot.CI["rs.y.stat",2])>0)) {
				spearman_mat=rbind(spearman_mat, matrix(c(paste(round(outputs$boot.CI["rs.y.stat",],round.digits), collapse=" ; "), print_confounder),nrow=1,ncol=2) )
				rownames(spearman_mat)[dim(spearman_mat)[1]]="rs.y.stat"
			}
		}

		else {
			# if CI diff.rs <0 : 	preference for eta -> xi
			if (outputs$boot.CI["diff.rs",1]<0 & outputs$boot.CI["diff.rs",2]<0) {
				spearman_mat["bootCI on diff","interpretation"]=print_eta_arrow_xi
	
				#Additionally, if the preferred cause is an obs_var and if CI rs.x sig. : add xi <-U-> eta
				if ((length(assumed_consequence_indicators)==1) & ((outputs$boot.CI["rs.x.stat",1]*outputs$boot.CI["rs.x.stat",2])>0)) {
					spearman_mat=rbind(spearman_mat, matrix(c(paste(round(outputs$boot.CI["rs.x.stat",],round.digits), collapse=" ; "), print_confounder),nrow=1,ncol=2) )
					rownames(spearman_mat)[dim(spearman_mat)[1]]="rs.x.stat"
				}
			}
			
			else {
				# if CI diff.rs not sig:	"data do not allow to show a causal direction"
				spearman_mat["bootCI on diff","interpretation"]=print_no_ccl
	
				#Additionally, if CI rs.x sig and CI rs.y sig. : 
				if ( ((outputs$boot.CI["rs.x.stat",1]*outputs$boot.CI["rs.x.stat",2])>0) & ((outputs$boot.CI["rs.y.stat",1]*outputs$boot.CI["rs.y.stat",2])>0) & ( (length(assumed_cause_indicators)==1) | (length(assumed_consequence_indicators)==1) ) ) {
					rs.x_rs.y_mat=matrix(nrow=3,ncol=2)
					rownames(rs.x_rs.y_mat)=c("rs.x.stat","rs.y.stat","rs.stat")
					rs.x_rs.y_mat["rs.x.stat",]=c(paste(round(outputs$boot.CI["rs.x.stat",],round.digits), collapse=" ; ") ,"")
					rs.x_rs.y_mat["rs.y.stat",]=c(paste(round(outputs$boot.CI["rs.y.stat",],round.digits), collapse=" ; ") ,"")	

					# if xi and eta are 2 obs_var : add xi <-U-> eta
					if ( (length(assumed_cause_indicators)==1) & (length(assumed_consequence_indicators)==1) ) {
						rs.x_rs.y_mat["rs.stat",]=c("",print_confounder)
						spearman_mat=rbind(spearman_mat, rs.x_rs.y_mat)
					}
					else {
						#if xi is obs_var and eta is latent:
						if ( length(assumed_cause_indicators)==1 ) {# & (length(assumed_consequence_indicators)>1)
							rs.x_rs.y_mat["rs.stat",]=c("",print_eta_arrow_xi_or_conf)
							spearman_mat=rbind(spearman_mat, rs.x_rs.y_mat)
						}

						#if eta is obs_var and xi is latent:
						if ( length(assumed_consequence_indicators)==1 ) {# & (length(assumed_cause_indicators)>1)
							rs.x_rs.y_mat["rs.stat",]=c("", print_xi_arrow_eta_or_conf)
							spearman_mat=rbind(spearman_mat, rs.x_rs.y_mat)
						}
					}		

					
				}
			}
		} #end else

		output.interpret=c(output.interpret, list("spearman"=spearman_mat))
	}
	if (sum(names(outputs$original_sample)=="dcor")){
		output.interpret=c(output.interpret, list("dcor"=list()))
		
		dcor_stat_mat=matrix(nrow=1,ncol=2)
		colnames(dcor_stat_mat)=c("results","interpretation")
		rownames(dcor_stat_mat)="bootCI on diff"

		dcor_stat_mat["bootCI on diff","results"]=paste(round(outputs$boot.CI["diff.dcor.stat",],round.digits), collapse=" ; ")

		# if CI diff.dcor.stat >0 : 	preference for xi -> eta
		if (outputs$boot.CI["diff.dcor.stat",1]>0 & outputs$boot.CI["diff.dcor.stat",2]>0) {
			dcor_stat_mat["bootCI on diff","interpretation"]=print_xi_arrow_eta

			#Additionally, if the preferred cause is an obs_var and if p.value dcor.y sig. : add xi <-U-> eta
			if ((length(assumed_cause_indicators)==1) & ((outputs$original_sample$dcor["p.value","dcor.y"])<alpha)) {
				dcor_stat_mat=rbind(dcor_stat_mat, matrix(c(round(outputs$original_sample$dcor["p.value","dcor.y"],round.digits), print_confounder),nrow=1,ncol=2) )
				rownames(dcor_stat_mat)[dim(dcor_stat_mat)[1]]="p.value_dcor.y"
			}
		}

		else {
			# if CI diff.dcor.stat <0 : 	preference for eta -> xi
			if (outputs$boot.CI["diff.dcor.stat",1]<0 & outputs$boot.CI["diff.dcor.stat",2]<0) {
				dcor_stat_mat["bootCI on diff","interpretation"]=print_eta_arrow_xi
	
				#Additionally, if the preferred cause is an obs_var and if p.value dcor.x sig. : add xi <-U-> eta
				if ((length(assumed_consequence_indicators)==1) & ((outputs$original_sample$dcor["p.value","dcor.x"])<alpha)) {
					dcor_stat_mat=rbind(dcor_stat_mat, matrix(c(round(outputs$original_sample$dcor["p.value","dcor.x"],round.digits), print_confounder),nrow=1,ncol=2) )
					rownames(dcor_stat_mat)[dim(dcor_stat_mat)[1]]="p.value_dcor.x"
				}
			}
			
			else {
				# if CI diff.dcor.stat not sig:	"data do not allow to show a causal direction"
				dcor_stat_mat["bootCI on diff","interpretation"]=print_no_ccl
	
				#Additionally, if p.value dcor.x sig and p.value dcor.y sig. : 
				if ( ((outputs$original_sample$dcor["p.value","dcor.x"])<alpha) & ((outputs$original_sample$dcor["p.value","dcor.y"])<alpha) & ( (length(assumed_cause_indicators)==1) | (length(assumed_consequence_indicators)==1) ) ) {
					dcor.x_dcor.y_mat=matrix(nrow=3,ncol=2)
					rownames(dcor.x_dcor.y_mat)=c("p.value_dcor.x","p.value_dcor.y","p.values_dcor")
					dcor.x_dcor.y_mat["p.value_dcor.x",]=c(round(outputs$original_sample$dcor["p.value","dcor.x"],round.digits),"")
					dcor.x_dcor.y_mat["p.value_dcor.y",]=c(round(outputs$original_sample$dcor["p.value","dcor.y"],round.digits),"")	

					# if xi and eta are 2 obs_var : add xi <-U-> eta
					if ( (length(assumed_cause_indicators)==1) & (length(assumed_consequence_indicators)==1) ) {
						dcor.x_dcor.y_mat["p.values_dcor",]=c("",print_confounder)
						dcor_stat_mat=rbind(dcor_stat_mat, dcor.x_dcor.y_mat)
					}
					else {
						#if xi is obs_var and eta is latent:
						if ( length(assumed_cause_indicators)==1 ) {# & (length(assumed_consequence_indicators)>1)
							dcor.x_dcor.y_mat["p.values_dcor",]=c("",print_eta_arrow_xi_or_conf)
							dcor_stat_mat=rbind(dcor_stat_mat, dcor.x_dcor.y_mat)
						}

						#if eta is obs_var and xi is latent:
						if ( length(assumed_consequence_indicators)==1 ) {# & (length(assumed_cause_indicators)>1)
							dcor.x_dcor.y_mat["p.values_dcor",]=c("", print_xi_arrow_eta_or_conf)
							dcor_stat_mat=rbind(dcor_stat_mat, dcor.x_dcor.y_mat)
						}
					}		

					
				}
			}
		} #end else

		output.interpret$dcor=c(output.interpret$dcor, list("dcor_stat"=dcor_stat_mat))

		if (as.logical(sum(rownames(outputs$boot.CI)=="diff.dcor.pvalue"))) {
		
			dcor_pval_mat=matrix(nrow=1,ncol=2)
			colnames(dcor_pval_mat)=c("results","interpretation")
			rownames(dcor_pval_mat)="bootCI on diff"
	
			dcor_pval_mat["bootCI on diff","results"]=paste(round(outputs$boot.CI["diff.dcor.pvalue",],round.digits), collapse=" ; ")
	
			# if CI diff.dcor.pvalue >0 : 	preference for eta -> xi
			if (outputs$boot.CI["diff.dcor.pvalue",1]>0 & outputs$boot.CI["diff.dcor.pvalue",2]>0) {
				dcor_pval_mat["bootCI on diff","interpretation"]=print_eta_arrow_xi
	
				#Additionally, if the preferred cause is an obs_var and if p.value dcor.x sig. : add xi <-U-> eta
				if ((length(assumed_consequence_indicators)==1) & ((outputs$original_sample$dcor["p.value","dcor.x"])<alpha)) {
					dcor_pval_mat=rbind(dcor_pval_mat, matrix(c(round(outputs$original_sample$dcor["p.value","dcor.x"],round.digits), print_confounder),nrow=1,ncol=2) )
					rownames(dcor_pval_mat)[dim(dcor_pval_mat)[1]]="p.value_dcor.x"
				}
			}	
	
			else {
				# if CI diff.dcor.pvalue <0 : 	preference for xi -> eta
				if (outputs$boot.CI["diff.dcor.pvalue",1]<0 & outputs$boot.CI["diff.dcor.pvalue",2]<0) {
					dcor_pval_mat["bootCI on diff","interpretation"]=print_xi_arrow_eta
		
					#Additionally, if the preferred cause is an obs_var and if p.value dcor.y sig. : add xi <-U-> eta
					if ((length(assumed_cause_indicators)==1) & ((outputs$original_sample$dcor["p.value","dcor.y"])<alpha)) {
						dcor_pval_mat=rbind(dcor_pval_mat, matrix(c(round(outputs$original_sample$dcor["p.value","dcor.y"],round.digits), print_confounder),nrow=1,ncol=2) )
						rownames(dcor_pval_mat)[dim(dcor_pval_mat)[1]]="p.value_dcor.y"
					}
				}
				
				else {
					# if CI diff.dcor.pvalue not sig:	"data do not allow to show a causal direction"
					dcor_pval_mat["bootCI on diff","interpretation"]=print_no_ccl
		
					#Additionally, if p.value dcor.x sig and p.value dcor.y sig. : 
					if ( ((outputs$original_sample$dcor["p.value","dcor.x"])<alpha) & ((outputs$original_sample$dcor["p.value","dcor.y"])<alpha) & ( (length(assumed_cause_indicators)==1) | (length(assumed_consequence_indicators)==1) ) ) {
						dcor.x_dcor.y_mat=matrix(nrow=3,ncol=2)
						rownames(dcor.x_dcor.y_mat)=c("p.value_dcor.x","p.value_dcor.y","p.values_dcor")
						dcor.x_dcor.y_mat["p.value_dcor.x",]=c(round(outputs$original_sample$dcor["p.value","dcor.x"],round.digits),"")
						dcor.x_dcor.y_mat["p.value_dcor.y",]=c(round(outputs$original_sample$dcor["p.value","dcor.y"],round.digits),"")	
	
						# if xi and eta are 2 obs_var : add xi <-U-> eta
						if ( (length(assumed_cause_indicators)==1) & (length(assumed_consequence_indicators)==1) ) {
							dcor.x_dcor.y_mat["p.values_dcor",]=c("",print_confounder)
							dcor_pval_mat=rbind(dcor_pval_mat, dcor.x_dcor.y_mat)
						}
						else {
							#if xi is obs_var and eta is latent:
							if ( length(assumed_cause_indicators)==1 ) {# & (length(assumed_consequence_indicators)>1)
								dcor.x_dcor.y_mat["p.values_dcor",]=c("",print_eta_arrow_xi_or_conf)
								dcor_pval_mat=rbind(dcor_pval_mat, dcor.x_dcor.y_mat)
							}
	
							#if eta is obs_var and xi is latent:
							if ( length(assumed_consequence_indicators)==1 ) {# & (length(assumed_cause_indicators)>1)
								dcor.x_dcor.y_mat["p.values_dcor",]=c("", print_xi_arrow_eta_or_conf)
								dcor_pval_mat=rbind(dcor_pval_mat, dcor.x_dcor.y_mat)
							}
						}		
	
						
					}
				}
			} #end else

			output.interpret$dcor=c(output.interpret$dcor, list("dcor_p.value"=dcor_pval_mat))	
		}

	}
	if (sum(names(outputs$original_sample)=="hsic")){
		output.interpret=c(output.interpret, list("hsic"=list()))

		
		hsic_stat_mat=matrix(nrow=1,ncol=2)
		colnames(hsic_stat_mat)=c("results","interpretation")
		rownames(hsic_stat_mat)="bootCI on diff"

		hsic_stat_mat["bootCI on diff","results"]=paste(round(outputs$boot.CI["diff.hsic.stat",],round.digits), collapse=" ; ")

		# if CI diff.hsic.stat >0 : 	preference for xi -> eta
		if (outputs$boot.CI["diff.hsic.stat",1]>0 & outputs$boot.CI["diff.hsic.stat",2]>0) {
			hsic_stat_mat["bootCI on diff","interpretation"]=print_xi_arrow_eta

			#Additionally, if the preferred cause is an obs_var and if p.value - bootstrap hsic.y sig. : add xi <-U-> eta
			if ((length(assumed_cause_indicators)==1) & ((outputs$original_sample$hsic["p.value - bootstrap","hsic.y"])<alpha)) {
				hsic_stat_mat=rbind(hsic_stat_mat, matrix(c(round(outputs$original_sample$hsic["p.value - bootstrap","hsic.y"],round.digits), print_confounder),nrow=1,ncol=2) )
				rownames(hsic_stat_mat)[dim(hsic_stat_mat)[1]]="p.value - bootstrap_hsic.y"
			}
		}

		else {
			# if CI diff.hsic.stat <0 : 	preference for eta -> xi
			if (outputs$boot.CI["diff.hsic.stat",1]<0 & outputs$boot.CI["diff.hsic.stat",2]<0) {
				hsic_stat_mat["bootCI on diff","interpretation"]=print_eta_arrow_xi
	
				#Additionally, if the preferred cause is an obs_var and if p.value - bootstrap hsic.x sig. : add xi <-U-> eta
				if ((length(assumed_consequence_indicators)==1) & ((outputs$original_sample$hsic["p.value - bootstrap","hsic.x"])<alpha)) {
					hsic_stat_mat=rbind(hsic_stat_mat, matrix(c(round(outputs$original_sample$hsic["p.value - bootstrap","hsic.x"],round.digits), print_confounder),nrow=1,ncol=2) )
					rownames(hsic_stat_mat)[dim(hsic_stat_mat)[1]]="p.value - bootstrap_hsic.x"
				}
			}
			
			else {
				# if CI diff.hsic.stat not sig:	"data do not allow to show a causal direction"
				hsic_stat_mat["bootCI on diff","interpretation"]=print_no_ccl
	
				#Additionally, if p.value - bootstrap hsic.x sig and p.value - bootstrap hsic.y sig. : 
				if ( ((outputs$original_sample$hsic["p.value - bootstrap","hsic.x"])<alpha) & ((outputs$original_sample$hsic["p.value - bootstrap","hsic.y"])<alpha) & ( (length(assumed_cause_indicators)==1) | (length(assumed_consequence_indicators)==1) ) ) {
					hsic.x_hsic.y_mat=matrix(nrow=3,ncol=2)
					rownames(hsic.x_hsic.y_mat)=c("p.value - bootstrap_hsic.x","p.value - bootstrap_hsic.y","p.value - bootstrap_hsic")
					hsic.x_hsic.y_mat["p.value - bootstrap_hsic.x",]=c(round(outputs$original_sample$hsic["p.value - bootstrap","hsic.x"],round.digits),"")
					hsic.x_hsic.y_mat["p.value - bootstrap_hsic.y",]=c(round(outputs$original_sample$hsic["p.value - bootstrap","hsic.y"],round.digits),"")	

					# if xi and eta are 2 obs_var : add xi <-U-> eta
					if ( (length(assumed_cause_indicators)==1) & (length(assumed_consequence_indicators)==1) ) {
						hsic.x_hsic.y_mat["p.value - bootstrap_hsic",]=c("",print_confounder)
						hsic_stat_mat=rbind(hsic_stat_mat, hsic.x_hsic.y_mat)
					}
					else {
						#if xi is obs_var and eta is latent:
						if ( length(assumed_cause_indicators)==1 ) {# & (length(assumed_consequence_indicators)>1)
							hsic.x_hsic.y_mat["p.value - bootstrap_hsic",]=c("",print_eta_arrow_xi_or_conf)
							hsic_stat_mat=rbind(hsic_stat_mat, hsic.x_hsic.y_mat)
						}

						#if eta is obs_var and xi is latent:
						if ( length(assumed_consequence_indicators)==1 ) {# & (length(assumed_cause_indicators)>1)
							hsic.x_hsic.y_mat["p.value - bootstrap_hsic",]=c("", print_xi_arrow_eta_or_conf)
							hsic_stat_mat=rbind(hsic_stat_mat, hsic.x_hsic.y_mat)
						}
					}		

					
				}
			}
		} #end else

		output.interpret$hsic=c(output.interpret$hsic, list("hsic_stat"=hsic_stat_mat))

		if (as.logical(sum(rownames(outputs$boot.CI)=="diff.hsic.gamma_pval"))) {
		
			hsic_pval_gamma_mat=matrix(nrow=1,ncol=2)
			colnames(hsic_pval_gamma_mat)=c("results","interpretation")
			rownames(hsic_pval_gamma_mat)="bootCI on diff"
	
			hsic_pval_gamma_mat["bootCI on diff","results"]=paste(round(outputs$boot.CI["diff.hsic.gamma_pval",],round.digits), collapse=" ; ")
	
			# if CI diff.hsic.gamma_pval >0 : 	preference for eta -> xi
			if (outputs$boot.CI["diff.hsic.gamma_pval",1]>0 & outputs$boot.CI["diff.hsic.gamma_pval",2]>0) {
				hsic_pval_gamma_mat["bootCI on diff","interpretation"]=print_eta_arrow_xi
	
				#Additionally, if the preferred cause is an obs_var and if p.value - gamma hsic.x sig. : add xi <-U-> eta
				if ((length(assumed_consequence_indicators)==1) & ((outputs$original_sample$hsic["p.value - gamma","hsic.x"])<alpha)) {
					hsic_pval_gamma_mat=rbind(hsic_pval_gamma_mat, matrix(c(round(outputs$original_sample$hsic["p.value - gamma","hsic.x"],round.digits), print_confounder),nrow=1,ncol=2) )
					rownames(hsic_pval_gamma_mat)[dim(hsic_pval_gamma_mat)[1]]="p.value - gamma_hsic.x"
				}
			}	
	
			else {
				# if CI diff.hsic.gamma_pval <0 : 	preference for xi -> eta
				if (outputs$boot.CI["diff.hsic.gamma_pval",1]<0 & outputs$boot.CI["diff.hsic.gamma_pval",2]<0) {
					hsic_pval_gamma_mat["bootCI on diff","interpretation"]=print_xi_arrow_eta
		
					#Additionally, if the preferred cause is an obs_var and if p.value - gamma hsic.y sig. : add xi <-U-> eta
					if ((length(assumed_cause_indicators)==1) & ((outputs$original_sample$hsic["p.value - gamma","hsic.y"])<alpha)) {
						hsic_pval_gamma_mat=rbind(hsic_pval_gamma_mat, matrix(c(round(outputs$original_sample$hsic["p.value - gamma","hsic.y"],round.digits), print_confounder),nrow=1,ncol=2) )
						rownames(hsic_pval_gamma_mat)[dim(hsic_pval_gamma_mat)[1]]="p.value - gamma_hsic.y"
					}
				}
				
				else {
					# if CI diff.hsic.gamma_pval not sig:	"data do not allow to show a causal direction"
					hsic_pval_gamma_mat["bootCI on diff","interpretation"]=print_no_ccl
		
					#Additionally, if p.value - gamma hsic.x sig and p.value - gamma hsic.y sig. : 
					if ( ((outputs$original_sample$hsic["p.value - gamma","hsic.x"])<alpha) & ((outputs$original_sample$hsic["p.value - gamma","hsic.y"])<alpha) & ( (length(assumed_cause_indicators)==1) | (length(assumed_consequence_indicators)==1) ) ) {
						hsic.x_hsic.y_mat=matrix(nrow=3,ncol=2)
						rownames(hsic.x_hsic.y_mat)=c("p.value - gamma_hsic.x","p.value - gamma_hsic.y","p.value - gamma_hsic")
						hsic.x_hsic.y_mat["p.value - gamma_hsic.x",]=c(round(outputs$original_sample$hsic["p.value - gamma","hsic.x"],round.digits),"")
						hsic.x_hsic.y_mat["p.value - gamma_hsic.y",]=c(round(outputs$original_sample$hsic["p.value - gamma","hsic.y"],round.digits),"")	
	
						# if xi and eta are 2 obs_var : add xi <-U-> eta
						if ( (length(assumed_cause_indicators)==1) & (length(assumed_consequence_indicators)==1) ) {
							hsic.x_hsic.y_mat["p.value - gamma_hsic",]=c("",print_confounder)
							hsic_pval_gamma_mat=rbind(hsic_pval_gamma_mat, hsic.x_hsic.y_mat)
						}
						else {
							#if xi is obs_var and eta is latent:
							if ( length(assumed_cause_indicators)==1 ) {# & (length(assumed_consequence_indicators)>1)
								hsic.x_hsic.y_mat["p.value - gamma_hsic",]=c("",print_eta_arrow_xi_or_conf)
								hsic_pval_gamma_mat=rbind(hsic_pval_gamma_mat, hsic.x_hsic.y_mat)
							}
	
							#if eta is obs_var and xi is latent:
							if ( length(assumed_consequence_indicators)==1 ) {# & (length(assumed_cause_indicators)>1)
								hsic.x_hsic.y_mat["p.value - gamma_hsic",]=c("", print_xi_arrow_eta_or_conf)
								hsic_pval_gamma_mat=rbind(hsic_pval_gamma_mat, hsic.x_hsic.y_mat)
							}
						}		
	
						
					}
				}
			} #end else

			output.interpret$hsic=c(output.interpret$hsic, list("hsic_p.value_gamma"=hsic_pval_gamma_mat))	
		}


		if (as.logical(sum(rownames(outputs$boot.CI)=="diff.hsic.gamma_pval"))) {
		
			hsic_pval_boot_mat=matrix(nrow=1,ncol=2)
			colnames(hsic_pval_boot_mat)=c("results","interpretation")
			rownames(hsic_pval_boot_mat)="bootCI on diff"
	
			hsic_pval_boot_mat["bootCI on diff","results"]=paste(round(outputs$boot.CI["diff.hsic.gamma_pval",],round.digits), collapse=" ; ")
	
			# if CI diff.hsic.gamma_pval >0 : 	preference for eta -> xi
			if (outputs$boot.CI["diff.hsic.gamma_pval",1]>0 & outputs$boot.CI["diff.hsic.gamma_pval",2]>0) {
				hsic_pval_boot_mat["bootCI on diff","interpretation"]=print_eta_arrow_xi
	
				#Additionally, if the preferred cause is an obs_var and if p.value - bootstrap hsic.x sig. : add xi <-U-> eta
				if ((length(assumed_consequence_indicators)==1) & ((outputs$original_sample$hsic["p.value - bootstrap","hsic.x"])<alpha)) {
					hsic_pval_boot_mat=rbind(hsic_pval_boot_mat, matrix(c(round(outputs$original_sample$hsic["p.value - bootstrap","hsic.x"],round.digits), print_confounder),nrow=1,ncol=2) )
					rownames(hsic_pval_boot_mat)[dim(hsic_pval_boot_mat)[1]]="p.value - bootstrap_hsic.x"
				}
			}	
	
			else {
				# if CI diff.hsic.gamma_pval <0 : 	preference for xi -> eta
				if (outputs$boot.CI["diff.hsic.gamma_pval",1]<0 & outputs$boot.CI["diff.hsic.gamma_pval",2]<0) {
					hsic_pval_boot_mat["bootCI on diff","interpretation"]=print_xi_arrow_eta
		
					#Additionally, if the preferred cause is an obs_var and if p.value - bootstrap hsic.y sig. : add xi <-U-> eta
					if ((length(assumed_cause_indicators)==1) & ((outputs$original_sample$hsic["p.value - bootstrap","hsic.y"])<alpha)) {
						hsic_pval_boot_mat=rbind(hsic_pval_boot_mat, matrix(c(round(outputs$original_sample$hsic["p.value - bootstrap","hsic.y"],round.digits), print_confounder),nrow=1,ncol=2) )
						rownames(hsic_pval_boot_mat)[dim(hsic_pval_boot_mat)[1]]="p.value - bootstrap_hsic.y"
					}
				}
				
				else {
					# if CI diff.hsic.gamma_pval not sig:	"data do not allow to show a causal direction"
					hsic_pval_boot_mat["bootCI on diff","interpretation"]=print_no_ccl
		
					#Additionally, if p.value - bootstrap hsic.x sig and p.value - bootstrap hsic.y sig. : 
					if ( ((outputs$original_sample$hsic["p.value - bootstrap","hsic.x"])<alpha) & ((outputs$original_sample$hsic["p.value - bootstrap","hsic.y"])<alpha) & ( (length(assumed_cause_indicators)==1) | (length(assumed_consequence_indicators)==1) ) ) {
						hsic.x_hsic.y_mat=matrix(nrow=3,ncol=2)
						rownames(hsic.x_hsic.y_mat)=c("p.value - bootstrap_hsic.x","p.value - bootstrap_hsic.y","p.value - bootstrap_hsic")
						hsic.x_hsic.y_mat["p.value - bootstrap_hsic.x",]=c(round(outputs$original_sample$hsic["p.value - bootstrap","hsic.x"],round.digits),"")
						hsic.x_hsic.y_mat["p.value - bootstrap_hsic.y",]=c(round(outputs$original_sample$hsic["p.value - bootstrap","hsic.y"],round.digits),"")	
	
						# if xi and eta are 2 obs_var : add xi <-U-> eta
						if ( (length(assumed_cause_indicators)==1) & (length(assumed_consequence_indicators)==1) ) {
							hsic.x_hsic.y_mat["p.value - bootstrap_hsic",]=c("",print_confounder)
							hsic_pval_boot_mat=rbind(hsic_pval_boot_mat, hsic.x_hsic.y_mat)
						}
						else {
							#if xi is obs_var and eta is latent:
							if ( length(assumed_cause_indicators)==1 ) {# & (length(assumed_consequence_indicators)>1)
								hsic.x_hsic.y_mat["p.value - bootstrap_hsic",]=c("",print_eta_arrow_xi_or_conf)
								hsic_pval_boot_mat=rbind(hsic_pval_boot_mat, hsic.x_hsic.y_mat)
							}
	
							#if eta is obs_var and xi is latent:
							if ( length(assumed_consequence_indicators)==1 ) {# & (length(assumed_cause_indicators)>1)
								hsic.x_hsic.y_mat["p.value - bootstrap_hsic",]=c("", print_xi_arrow_eta_or_conf)
								hsic_pval_boot_mat=rbind(hsic_pval_boot_mat, hsic.x_hsic.y_mat)
							}
						}		
	
						
					}
				}
			} #end else

			output.interpret$hsic=c(output.interpret$hsic, list("hsic_p.value_boot"=hsic_pval_boot_mat))	
		}

	}

	return(output.interpret)
}


rbvnorm=function(n,rho){
	x=rnorm(n)
	z=rnorm(n)
	y=rho*x+sqrt(1-rho^2)*z
	return(cbind(x,y))
}

glink=function(to_transform, based_on){
	#correction made on March 26th 2022 : it concerns cumfreq computation
	result=rep(NA,length(to_transform))
	for (i in 1:length(to_transform)){
		cumfreq=sum(to_transform[i]>=to_transform)/length(to_transform)
		result[i]=quantile(based_on, probs=cumfreq)
	}
	return(result)
}

HSIC_HT=function(data, assumed_cause_indicators, assumed_consequence_indicators, CI=c(0.025,0.975), B=1000, seed=1234)
{

	#This function helps to compute a percentile bootstrap confidence interval for the difference of two p.values (each one corresponding to a possible causal direction). Each p.value comes from a gamma approximation and informs about an independence hypothesis (using HSIC).  
	#Building samples without the concerned causation ("H0 samples"), this function also gives an estimated proportion of differences (assuming H0: there is no causation) below or equal the observed difference in the original dataset.
	#data is a data.frame
	#"assumed_cause_indicators" is a vector of characters giving the names of the indicators for the assumed (latent) cause.
	#"assumed_consequence_indicators" is a vector of characters giving the names of the indicators for the assumed (latent) consequence.
	#CI is a vector (length equal to 2) giving the bounds for the percentiles confidence intervals to compute by bootstrap. Default is c(0.025,0.975).
	#B is the number of bootstrap samples to draw before computing the additional percentiles bootstrap confidence intervals.
	#seed is a value used to initialize the randoms generators before computing first statistics on the original dataset and before generating the bootstrapped samples (and computing the related statistics).

	#The return value "output" is a list giving:
	#	[1] statistics from the original sample (uses "causation_stat()" function); 
	#	[2] differences of directional p.values computed with "causation_stat()" for each bootstrap sample + associated differences for H0 samples  
	#	[3] information about the number of bootstrapped sample where some statistics could not be computed. 
	#	[4] the asked additional bootstrapped CI.  
	#	[5] the estimated proportion of datasets (assuming H0) giving a difference score below or equal to the observed difference (of HSIC gamma's directional probabilities) in the original dataset.

	if (is.numeric(seed))	set.seed(seed)
	
	n=dim(data)[1]

	output=list()

	#Analysis of the original sample:
	original_sample=causation_stat(data=data, assumed_cause_indicators=assumed_cause_indicators, assumed_consequence_indicators=assumed_consequence_indicators, indep.measures="hsic", dcor.pvalues=FALSE, HSIC.pvalues=c("gamma"))
	diff_p.value_gamma=original_sample$hsic["p.value - gamma","diff.hsic"]
	output=c(output, list("original_sample"=original_sample))

	if (is.numeric(seed))	set.seed(seed)

	#bootstrapped statistics:
	boot_res=matrix(nrow=3,ncol=B)
	rownames(boot_res)=c("diff.hsic.gamma_pval", "diff.hsic.gamma_pval_HT_bvnorm", "diff.hsic.gamma_pval_HT_keep_marg")

	#sample selection:
	bootstrap.select=matrix( sample(x=1:n, size=(B*n),replace=TRUE) , ncol=B )


	#computation and storage (for each sample) of the bootrapped estimated statistics

	for (i in 1:B){
		data.boot=data[bootstrap.select[,i],]
		boot_sample=causation_stat(data=data.boot,assumed_cause_indicators=assumed_cause_indicators, assumed_consequence_indicators=assumed_consequence_indicators, indep.measures="hsic", dcor.pvalues=FALSE, HSIC.pvalues=c("gamma"))
		boot_res[1,i]=boot_sample$hsic["p.value - gamma","diff.hsic"]
		
		#computation of H0 samples differences:

		resXi=PCA(data.boot[,assumed_cause_indicators], graph=FALSE)
		Fx=resXi$ind$coord[,1]
		resEta=PCA(data.boot[,assumed_consequence_indicators], graph=FALSE)
		Fy=resEta$ind$coord[,1]

		rp=cor(Fx,Fy)


		z=rbvnorm(n,rp)

		Fx_ctrl_norm=z[,1]
		Fy_ctrl_norm=z[,2]
		resx_ctrl_norm=residuals(lm(Fx_ctrl_norm~Fy_ctrl_norm))
		resy_ctrl_norm=residuals(lm(Fy_ctrl_norm~Fx_ctrl_norm))

		save_seed=.Random.seed

		hsic.gamma.y_ctrl_norm = dhsic.test(Fx_ctrl_norm, resy_ctrl_norm, method = "gamma", kernel = "gaussian")	 
		hsic.gamma.x_ctrl_norm = dhsic.test(Fy_ctrl_norm, resx_ctrl_norm, method = "gamma", kernel = "gaussian")

		.Random.seed<<-save_seed #assignation to a global variable ".Random.seed"

		boot_res[2,i]=hsic.gamma.x_ctrl_norm$p.value-hsic.gamma.y_ctrl_norm$p.value #computing a difference score between both p-values (a negative difference is assumed to mean Xi->Eta) 


		Fx_ctrl=glink(z[,1],Fx)
		Fy_ctrl=glink(z[,2],Fy)
			
		resx_ctrl=residuals(lm(Fx_ctrl~Fy_ctrl))
		resy_ctrl=residuals(lm(Fy_ctrl~Fx_ctrl))

		save_seed=.Random.seed

		hsic.gamma.y_ctrl = dhsic.test(Fx_ctrl, resy_ctrl, method = "gamma", kernel = "gaussian")	 
		hsic.gamma.x_ctrl = dhsic.test(Fy_ctrl, resx_ctrl, method = "gamma", kernel = "gaussian")

		.Random.seed<<-save_seed #assignation to a global variable ".Random.seed"

		boot_res[3,i]=hsic.gamma.x_ctrl$p.value-hsic.gamma.y_ctrl$p.value #computing a difference score between both p-values (a negative difference is assumed to mean Xi->Eta) 

	}

	nb_problematic_bootstrapped_samples=sum((apply(is.na(boot_res),2,sum))!=0)
	info=(paste("There is", nb_problematic_bootstrapped_samples, "bootstrapped samples where some statistics could not be computed in boot.results. Please check it before interpreting boot.CI outputs !" ))	

	boot.CI=t( apply(boot_res, 1, quantile, probs=CI, na.rm=TRUE))

	estimated_prop_bvnorm=sum(diff_p.value_gamma>=boot_res[2,])/length(boot_res[2,])
	estimated_prop=sum(diff_p.value_gamma>=boot_res[3,])/length(boot_res[3,])

	output=c(output, list("boot.results"=boot_res),list("boot.info"=info), list("boot.CI"=boot.CI), list("estimated_proportion_bvnorm"=estimated_prop_bvnorm), list("estimated_proportion_kmarg"=estimated_prop))

	return(output)

}


############################################################################################
############################################################################################

MCsim=function(model, n_MC_samples, sample_size, obs_var="all", seed=4321, assumed_cause_indicators, assumed_consequence_indicators, indep.measures="all", CI=c(0.025,0.975), boot.pval=TRUE, dcor.R=500, HSIC.B=500, B=1000){

	#Function for testing on samples (generated by a MonteCarlo approach) the causal inference algorithm implemented in the function "LC()". It uses "causal_sim()", "LC()" and "causal.interpret()"
	#Warning: in comparaison with the first version of this fonction, some "print()" have been added.   

	# model is a string like:'
	#	#Exogenous variables (come first)
	#	U ~ unif
	#	#Latent Endogeneous variables
	#	X ~ 0.2*U + chisq
	#	Y ~ 0.3*U + 0.5*X + unif
	#
	#	#Observed indicators
	#	
	#	x1 ~ 0.7*X + norm
	#	x2 ~ 0.8*X + norm
	#	x3 ~ 0.9*X + norm
	#	
	#	y1 ~ 0.7*Y + norm
	#	y2 ~ 0.8*Y + norm
	#	y3 ~ 0.9*Y + norm	
	#'

	#n_MC_samples: number of MonteCarlo samples to generate
	#sample_size : the number of statistical individuals in each dataset
	#obs_var: indicates for which variables to keep records in the generated dataframes
	#seed: a value to initialize seeds in "causal_sim()" and "LC()".	

	#"assumed_cause_indicators" is a vector of characters giving the names of the indicators for the assumed (latent) cause.
	#"assumed_consequence_indicators" is a vector of characters giving the names of the indicators for the assumed (latent) consequence. 
	# when "assumed_cause_indicators" contains only one variable name, it is assumed this indicator has no error of measurement (which means the observed variable and the latent concept are the same). The same is also true for the "assumed_consequence_indicators" vector.
	# when either the "assumed_cause_indicators" or "assumed_consequence_indicators" (or both) is only one observed variable, and when some relevant additional tests are significants, more information will be returned about possible existing latent confounders.

	#indep.measures: which measures should be used to assess independency in LC()? indep.measures is a vector which can include "spearman", "dcor", "hsic".
	#CI= a vector of 2 numerical values giving the quantiles for the lower and upper bound of the confidence intervals. (CI is also used to compute the alpha to interpret p.values in "causal.interpret")
	#boot.pval: logical, do we need additional bootstrap CI for the difference of dcor's p.values and for the difference of hsic's gamma_p.value.
	#dcor.R: in "LC()", the number of replicates used for estimating each dcor.pvalue
	#HSIC.B: in "LC()", the number of resamples used to compute each bootstrap's HSIC p-value
	#B: in "LC()", the number of resamples used to compute the additional percentiles bootstrap CI.

	#a list is returned where :
		#"datasets" is a list of generated datasets by MonteCarlo (causal_sim())
		#"results" contains the various interpretations given for each method for each generated dataset.
		#"seed_info"= the seed used initially
		#"assumed_directionality": the directionality first specified by the user (which can be different from the conclusions of the algorithm).
		#...


	methods=c()
	if ( (sum(indep.measures=="spearman")>0) | (sum(indep.measures=="all")>0) ) {
		methods=c(methods,"spearman")
	}
	if ( (sum(indep.measures=="dcor")>0) | (sum(indep.measures=="all")>0) ) {

		methods=c(methods,"dcor$dcor_stat")

		if ( boot.pval ) {
			methods=c(methods, "dcor$dcor_p.value")
		}
	}
	if ( (sum(indep.measures=="hsic")>0) | (sum(indep.measures=="all")>0) ) {

		methods=c(methods,"hsic$hsic_stat")

		if ( boot.pval ) {
			methods=c(methods, "hsic$hsic_p.value_gamma", "hsic$hsic_p.value_boot")
		}
	}

	

	datasets_list=causal_sim(model=model, n=sample_size, n_samples=n_MC_samples, obs_var=obs_var, seed=seed)

	print("model")
	print(model)
	print("datasets")
	print(datasets_list)

	results_mat=matrix(nrow=length(methods),ncol=n_MC_samples)
	rownames(results_mat)=c(methods)

	for (sampl_index in 1:n_MC_samples){
		res_tmp=LC( (datasets_list[[sampl_index]])[,c(assumed_cause_indicators, assumed_consequence_indicators)] , assumed_cause_indicators=assumed_cause_indicators, assumed_consequence_indicators=assumed_consequence_indicators, indep.measures=indep.measures , CI=CI, boot.pval=boot.pval, dcor.R=dcor.R, HSIC.B=HSIC.B, B=B, seed=seed)
		res_interpr_tmp=causal.interpret(res_tmp, assumed_cause_indicators=assumed_cause_indicators, assumed_consequence_indicators=assumed_consequence_indicators)

		for (method_index in 1:length(methods)){
			method_tmp=eval(parse( text=paste("res_interpr_tmp$",methods[method_index],sep="") ))
			results_mat[method_index,sampl_index]=paste(method_tmp[,"interpretation"], collapse=" ")
			
		}
		print("original dataset")
		print(res_tmp[[1]])
		print("results")
		print(results_mat[,sampl_index])
	}

	print("seed_info")
	print(seed)
	print("assumed_directionality")
	print(paste("(",paste(assumed_cause_indicators, collapse=", "),")", "-->", "(",paste(assumed_consequence_indicators, collapse=", "),")","was first assumed"))

	print("CI")
	print(CI)
	print("dcor.R")
	print(dcor.R)
	print("HSIC.B")
	print(HSIC.B)
	print("number_of_bootstrap_B")
	print(B)

	return(list("model"=model,"sample_size"=sample_size, "datasets"=datasets_list,"results"=results_mat, "seed_info"=seed, "assumed_directionality"=paste("(",paste(assumed_cause_indicators, collapse=", "),")", "-->", "(",paste(assumed_consequence_indicators, collapse=", "),")","was first assumed"), "CI"=CI, "dcor.R"=dcor.R, "HSIC.B"=HSIC.B, "number_of_bootstrap_B"=B ))

}

MCsimHT=function(model, n_MC_samples, sample_size, obs_var="all", seed=4321, assumed_cause_indicators, assumed_consequence_indicators, CI=c(0.025,0.975), B=1000){

	#Function for testing on samples (generated by a MonteCarlo approach) the causal inference algorithm implemented in the function "HSIC_HT()". It uses "causal_sim()", "HSIC_HT()"

	# model is a string like:'
	#	#Exogenous variables (come first)
	#	U ~ unif
	#	#Latent Endogeneous variables
	#	X ~ 0.2*U + chisq
	#	Y ~ 0.3*U + 0.5*X + unif
	#
	#	#Observed indicators
	#	
	#	x1 ~ 0.7*X + norm
	#	x2 ~ 0.8*X + norm
	#	x3 ~ 0.9*X + norm
	#	
	#	y1 ~ 0.7*Y + norm
	#	y2 ~ 0.8*Y + norm
	#	y3 ~ 0.9*Y + norm	
	#'

	#n_MC_samples: number of MonteCarlo samples to generate
	#sample_size : the number of statistical individuals in each dataset
	#obs_var: indicates for which variables to keep records in the generated dataframes
	#seed: a value to initialize seeds in "causal_sim()" and "HSIC_HT()".	

	#"assumed_cause_indicators" is a vector of characters giving the names of the indicators for the assumed (latent) cause.
	#"assumed_consequence_indicators" is a vector of characters giving the names of the indicators for the assumed (latent) consequence. 
	# when "assumed_cause_indicators" contains only one variable name, it is assumed this indicator has no error of measurement (which means the observed variable and the latent concept are the same). The same is also true for the "assumed_consequence_indicators" vector.

	#CI= a vector of 2 numerical values giving the quantiles for the lower and upper bound of the confidence intervals.
	#B: in "HSIC_HT()", the number of resamples used to compute the additional percentiles bootstrap CI.

	#a list is returned where :
		#"datasets" is a list of generated datasets by MonteCarlo (causal_sim())
		#"results" contains the various statistics given by each method for each generated dataset.
		#"seed_info"= the seed used initially
		#"assumed_directionality": the directionality first specified by the user (which can be different from the conclusions of the algorithm).
		#...

	datasets_list=causal_sim(model=model, n=sample_size, n_samples=n_MC_samples, obs_var=obs_var, seed=seed)

	print("model")
	print(model)
	#print("datasets")
	#print(datasets_list)

	results_mat=matrix(nrow=n_MC_samples,ncol=14)
	colnames(results_mat)=c("orig_sample_hsic.x_stat","orig_sample_hsic.y_stat","orig_sample_diff_hsic_stat",
	"orig_sample_hsic.x_pval_gamma","orig_sample_hsic.y_pval_gamma","orig_sample_diff_hsic_pval_gamma",
	"bootCIdiff.hsic.gamma_pval_lower_bound", "bootCIdiff.hsic.gamma_pval_upper_bound",
	"bootCIdiff.hsic.gamma_pval_lower_bound_for_HT_bvnorm", "bootCIdiff.hsic.gamma_pval_upper_bound_for_HT_bvnorm",
	"estimated_proportion_of_BvNormSamples_be_obs_diff",
	"bootCIdiff.hsic.gamma_pval_lower_bound_for_HT_kmarg", "bootCIdiff.hsic.gamma_pval_upper_bound_for_HT_kmarg",
	"estimated_proportion_of_KMargSamples_be_obs_diff")

	cat(colnames(results_mat))
	cat("\n")

	for (sampl_index in 1:n_MC_samples){
		res_tmp=HSIC_HT( (datasets_list[[sampl_index]])[,c(assumed_cause_indicators, assumed_consequence_indicators)] , assumed_cause_indicators=assumed_cause_indicators, assumed_consequence_indicators=assumed_consequence_indicators, CI=CI, B=B, seed=seed)
		results_mat[sampl_index,c("orig_sample_hsic.x_stat", "orig_sample_hsic.y_stat", "orig_sample_diff_hsic_stat")]=res_tmp$original_sample$hsic["statistic",c("hsic.x","hsic.y","diff.hsic")]
		results_mat[sampl_index,c("orig_sample_hsic.x_pval_gamma", "orig_sample_hsic.y_pval_gamma", "orig_sample_diff_hsic_pval_gamma")]=res_tmp$original_sample$hsic["p.value - gamma",c("hsic.x","hsic.y","diff.hsic")]
		results_mat[sampl_index,c("bootCIdiff.hsic.gamma_pval_lower_bound", "bootCIdiff.hsic.gamma_pval_upper_bound")]=res_tmp$boot.CI["diff.hsic.gamma_pval",]
		results_mat[sampl_index,c("bootCIdiff.hsic.gamma_pval_lower_bound_for_HT_bvnorm", "bootCIdiff.hsic.gamma_pval_upper_bound_for_HT_bvnorm")]=res_tmp$boot.CI["diff.hsic.gamma_pval_HT_bvnorm",]
		results_mat[sampl_index,"estimated_proportion_of_BvNormSamples_be_obs_diff"]=res_tmp$estimated_proportion_bvnorm
		results_mat[sampl_index,c("bootCIdiff.hsic.gamma_pval_lower_bound_for_HT_kmarg", "bootCIdiff.hsic.gamma_pval_upper_bound_for_HT_kmarg")]=res_tmp$boot.CI["diff.hsic.gamma_pval_HT_keep_marg",]
		results_mat[sampl_index,"estimated_proportion_of_KMargSamples_be_obs_diff"]=res_tmp$estimated_proportion_kmarg

		cat(results_mat[sampl_index,])
		cat("\n")
	}

	print("seed_info")
	print(seed)
	print("assumed_directionality")
	print(paste("(",paste(assumed_cause_indicators, collapse=", "),")", "-->", "(",paste(assumed_consequence_indicators, collapse=", "),")","was first assumed"))

	print("CI")
	print(CI)
	print("number_of_bootstrap_B")
	print(B)

	return(list("model"=model,"sample_size"=sample_size, "datasets"=datasets_list,"results"=results_mat, "seed_info"=seed, "assumed_directionality"=paste("(",paste(assumed_cause_indicators, collapse=", "),")", "-->", "(",paste(assumed_consequence_indicators, collapse=", "),")","was first assumed"), "CI"=CI, "number_of_bootstrap_B"=B))

}

############################################################################################
############################################################################################

model='
	#A latent confounder U and no causal effect between X and Y 
	#Exogenous variables (come first)
	U ~ chisq
	#Latent Endogeneous variables
	X ~ 0.7*U + chisq
	Y ~ 0.7*U + norm # test
	#Observed indicators
	
	x1 ~ 0.7*X + norm
	x2 ~ 0.8*X + norm
	x3 ~ 0.9*X + norm
	
	y1 ~ 0.7*Y + norm
	y2 ~ 0.8*Y + norm
	y3 ~ 0.9*Y + norm
	
'


resu=MCsim(model=model, n_MC_samples=1000, sample_size=500, obs_var=c("x1","x2","x3","y1","y2","y3"), seed=4321, assumed_cause_indicators=c("x1","x2","x3"), assumed_consequence_indicators=c("y1","y2","y3"), indep.measures="all", CI=c(0.025,0.975), boot.pval=TRUE, dcor.R=500, HSIC.B=500, B=1000)
resu

resu2=MCsimHT(model=model, n_MC_samples=1000, sample_size=500, obs_var=c("x1","x2","x3","y1","y2","y3"), seed=4321, assumed_cause_indicators=c("x1","x2","x3"), assumed_consequence_indicators=c("y1","y2","y3"), CI=c(0.025,0.975), B=1000)
resu2[c(1:2,4:8)]