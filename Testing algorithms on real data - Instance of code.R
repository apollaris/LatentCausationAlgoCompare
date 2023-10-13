#Main author : A. Pollaris
#Date : October 13, 2023
#Concern: implementations in the R language to generate data for noisy measures of latent variables and perform causal inference. This code can help providing results to compare 3 causal inference algorithms:

#1. The factor scores computation followed by the Direction Dependence Analysis independence component (FSC+DDA)
#2. The Latent Causation Algorithm (LC)
#3. The Hypothesis Test Algorithm (HT)

This code is adapted to run on some real data from the Tuebingen's website (i.e. https://webdav.tuebingen.mpg.de/cause-effect/). Noisy data can then be generated from real data.

#In this instance of code, causal inference algorithms for latent pairs of variables are implemented. Noisy data are then introduced as an input for causal inference.

#For every statistical analysis, the author did his best to make sure this code is working properly. However, if some problems or inconsistencies are found using the code below,
#please inform the author, A. Pollaris (Arnaud.Pollaris@ulb.ac.be), to help him to improve this work.  

#By using this code, you accept to take all responsibility for any negative consequences that might result. This code is provided "as is", without warranty of any kind,
#express or implied. In no event shall the author or copyright holders be liable for any claim, damages or other liability.

#N.B.: 
#1. Some parts of the code below come from the DDA project (available online : https://www.ddaproject.com/ ; see also Wolfgang Wiedermann, Xintong Li and Alexander von Eye)
#   and were not written by A. Pollaris.
#2. Some functions below were already available in https://github.com/apollaris/LatentCausation
#3. Note that some parts of the current implementation are still in a beta version.

###############################################################################################################################################################################
###############################################################################################################################################################################

#printing options:
options(width=10000)

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

sim_measurement=function(model, latent_dataset, seed=4321){
	#	#Observed indicators
	#	
	#	x1 ~ 0.7*X + norm
	#	x2 ~ 0.8*X + norm
	#	x3 ~ 0.9*X + norm
	#	
	#	y1 ~ 0.7*Y + norm
	#	y2 ~ 0.8*Y + norm
	#	y3 ~ 0.9*Y + norm	
	#
	#	#where X and Y are latent variables
	#'

	# latent dataset: a data.frame containing every (and only) the latent variables of the considered model.
	# seed: to indicate a seed at the beginning of the simulations

	# Returns a dataframe 

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

	n=dim(latent_dataset)[1]
	VL_names=colnames(latent_dataset)

	#standardization of the VL:
	for (Vcol in VL_names){
		latent_dataset[,Vcol]=(latent_dataset[,Vcol]-mean(latent_dataset[,Vcol]))/sd(latent_dataset[,Vcol])
	}

	data_mat=cbind(as.matrix(latent_dataset), matrix(nrow=n,ncol=dim(resid_shapes)[2]))
	colnames(data_mat)=c(VL_names,colnames(resid_shapes))
	


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
	return(as.data.frame(data_mat[,colnames(resid_shapes)]))
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

HSIC_HT_v2=function(data, assumed_cause_indicators, assumed_consequence_indicators, CI=c(0.025,0.975), B=1000, seed=1234)
{

	#Update HSIC_HT -> HSIC_HT_v2 : remplacement of the line : "boot_res[1,i]=boot_sample$hsic["p.value - gamma","diff.hsic"]" by "boot_res[1,i]=round(boot_sample$hsic["p.value - gamma","hsic.x"],3)-round(boot_sample$hsic["p.value - gamma","hsic.y"],3)" (for rounded p.values) 

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
		boot_res[1,i]=round(boot_sample$hsic["p.value - gamma","hsic.x"],3)-round(boot_sample$hsic["p.value - gamma","hsic.y"],3)
		
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

tuebingen_data=function(metadta, dataset_vec, folder_datasets, measurement_model, CI=c(0.025,0.975), dcor.R=500, HSIC.B=500, B=1000, seed=4321){
	
	if (is.numeric(seed))	set.seed(seed)
	
	#Matrix of causal results:
	coln=c("weight", "sample_size",
	"orig_sample_rs.x", "orig_sample_rs.y", "orig_sample_diff.rs", "bootCI_rs_low", "bootCI_rs_high", 
	"orig_sample_dcor.x_stat","orig_sample_dcor.y_stat","orig_sample_diff_dcor_stat", "bootCI_dcor_stat_low", "bootCI_dcor_stat_high",
	"orig_sample_dcor.x_pval", "orig_sample_dcor.y_pval", "orig_sample_diff_dcor_pval", "bootCI_dcor_pval_low", "bootCI_dcor_pval_high",
	"orig_sample_hsic.x_stat","orig_sample_hsic.y_stat","orig_sample_diff_hsic_stat", "bootCI_hsic_stat_low", "bootCI_hsic_stat_high",
	"orig_sample_hsic.x_pval_gamma", "orig_sample_hsic.y_pval_gamma", "orig_sample_diff_hsic_pval_gamma", "bootCI_hsic_pval_gamma_low", "bootCI_hsic_pval_gamma_high",
	"orig_sample_hsic.x_pval_boot", "orig_sample_hsic.y_pval_boot", "orig_sample_diff_hsic_pval_boot",
	"bootCIdiff.hsic.gamma_pval_r3_lower_bound", "bootCIdiff.hsic.gamma_pval_r3_upper_bound",
	"bootCIdiff.hsic.gamma_pval_lower_bound_for_HT_bvnorm", "bootCIdiff.hsic.gamma_pval_upper_bound_for_HT_bvnorm",
	"estimated_proportion_of_BvNormSamples_be_obs_diff",
	"bootCIdiff.hsic.gamma_pval_lower_bound_for_HT_kmarg", "bootCIdiff.hsic.gamma_pval_upper_bound_for_HT_kmarg",
	"estimated_proportion_of_KMargSamples_be_obs_diff"
	)
	results_mat=matrix(ncol=length(coln),nrow=length(dataset_vec))
	colnames(results_mat)=coln
	rownames(results_mat)=dataset_vec #rownames(results_mat) is a vector of char

	#building of a list of datasets:

	datasets_list=list()

	for (i in dataset_vec) {
	
		dtapath=paste(folder_datasets,"pair0","0"[i<100],"0"[i<10],i,".txt", sep="")
		dta=read.table(dtapath)

		mtadta=metadta[metadta$V1==i,]

		if (length(mtadta$V2:mtadta$V3)!=1 | length(mtadta$V4:mtadta$V5)!=1) stop("'Latent' cause or consequence can not be multivariate here.")

		dta=data.frame(dta[,mtadta$V2],dta[,mtadta$V4]) #also permute the columns if the consequence is before the cause
		colnames(dta)=c("Vca","Vcsq")
	
		if (is.numeric(seed))	set.seed(seed)   
   
    		if (1000<dim(dta)[1]) {
        		dta=dta[sample(1:(dim(dta)[1]), size=1000, replace = FALSE),]
    		}

		if (is.numeric(seed))	set.seed(seed)

		datasets_list[[i]]=sim_measurement(measurement_model, dta, seed=seed)

	}

	#treatement of each dataset (dataset_vec):

	for (i in dataset_vec) {

		mtadta=metadta[metadta$V1==i,]
	
		results_mat[as.character(i),"weight"]=(mtadta$V6)[1]
		results_mat[as.character(i),"sample_size"]=dim(datasets_list[[i]])[1]

		assumed_cause_indicators=colnames(datasets_list[[i]])[grep("Vca", colnames(datasets_list[[i]]))]
		assumed_consequence_indicators=colnames(datasets_list[[i]])[grep("Vcsq", colnames(datasets_list[[i]]))]
		
		res_tmp=LC( (datasets_list[[i]])[,c(assumed_cause_indicators, assumed_consequence_indicators)] , assumed_cause_indicators=assumed_cause_indicators, assumed_consequence_indicators=assumed_consequence_indicators, indep.measures="all", CI=CI, boot.pval=TRUE, dcor.R=dcor.R, HSIC.B=HSIC.B, B=B, seed=seed)
	
		print(res_tmp[c(1,3,4)])
		print("\n")

		results_mat[as.character(i),"orig_sample_rs.x"]=res_tmp$original_sample$spearman[,"rs.x"]
		results_mat[as.character(i),"orig_sample_rs.y"]=res_tmp$original_sample$spearman[,"rs.y"]
		results_mat[as.character(i),"orig_sample_diff.rs"]=res_tmp$original_sample$spearman[,"diff.rs"]
		results_mat[as.character(i),"bootCI_rs_low"]=res_tmp$boot.CI["diff.rs",1]
		results_mat[as.character(i),"bootCI_rs_high"]=res_tmp$boot.CI["diff.rs",2]

		results_mat[as.character(i),"orig_sample_dcor.x_stat"]=res_tmp$original_sample$dcor["statistic","dcor.x"]
		results_mat[as.character(i),"orig_sample_dcor.y_stat"]=res_tmp$original_sample$dcor["statistic","dcor.y"]
		results_mat[as.character(i),"orig_sample_diff_dcor_stat"]=res_tmp$original_sample$dcor["statistic","diff.dcor"]
		results_mat[as.character(i),"bootCI_dcor_stat_low"]=res_tmp$boot.CI["diff.dcor.stat",1]
		results_mat[as.character(i),"bootCI_dcor_stat_high"]=res_tmp$boot.CI["diff.dcor.stat",2]

		results_mat[as.character(i),"orig_sample_dcor.x_pval"]=res_tmp$original_sample$dcor["p.value","dcor.x"]
		results_mat[as.character(i),"orig_sample_dcor.y_pval"]=res_tmp$original_sample$dcor["p.value","dcor.y"]
		results_mat[as.character(i),"orig_sample_diff_dcor_pval"]=res_tmp$original_sample$dcor["p.value","diff.dcor"]
		results_mat[as.character(i),"bootCI_dcor_pval_low"]=res_tmp$boot.CI["diff.dcor.pvalue",1]
		results_mat[as.character(i),"bootCI_dcor_pval_high"]=res_tmp$boot.CI["diff.dcor.pvalue",2]

		results_mat[as.character(i),"orig_sample_hsic.x_stat"]=res_tmp$original_sample$hsic["statistic","hsic.x"]
		results_mat[as.character(i),"orig_sample_hsic.y_stat"]=res_tmp$original_sample$hsic["statistic","hsic.y"]
		results_mat[as.character(i),"orig_sample_diff_hsic_stat"]=res_tmp$original_sample$hsic["statistic","diff.hsic"]
		results_mat[as.character(i),"bootCI_hsic_stat_low"]=res_tmp$boot.CI["diff.hsic.stat",1]
		results_mat[as.character(i),"bootCI_hsic_stat_high"]=res_tmp$boot.CI["diff.hsic.stat",2]

		results_mat[as.character(i),"orig_sample_hsic.x_pval_gamma"]=res_tmp$original_sample$hsic["p.value - gamma","hsic.x"]
		results_mat[as.character(i),"orig_sample_hsic.y_pval_gamma"]=res_tmp$original_sample$hsic["p.value - gamma","hsic.y"]
		results_mat[as.character(i),"orig_sample_diff_hsic_pval_gamma"]=res_tmp$original_sample$hsic["p.value - gamma","diff.hsic"]
		results_mat[as.character(i),"bootCI_hsic_pval_gamma_low"]=res_tmp$boot.CI["diff.hsic.gamma_pval",1]
		results_mat[as.character(i),"bootCI_hsic_pval_gamma_high"]=res_tmp$boot.CI["diff.hsic.gamma_pval",2]
	
		results_mat[as.character(i),"orig_sample_hsic.x_pval_boot"]=res_tmp$original_sample$hsic["p.value - bootstrap","hsic.x"]
		results_mat[as.character(i),"orig_sample_hsic.y_pval_boot"]=res_tmp$original_sample$hsic["p.value - bootstrap","hsic.y"]
		results_mat[as.character(i),"orig_sample_diff_hsic_pval_boot"]=res_tmp$original_sample$hsic["p.value - bootstrap","diff.hsic"]

		####Results of using HSIC_HT_v2() to generate "H0" datasets
		res_tmp=NULL
		res_tmp=HSIC_HT_v2( (datasets_list[[i]])[,c(assumed_cause_indicators, assumed_consequence_indicators)] , assumed_cause_indicators=assumed_cause_indicators, assumed_consequence_indicators=assumed_consequence_indicators, CI=CI, B=B, seed=seed)

		print(res_tmp[c(1,3,4,5,6)])
		print("\n")

		results_mat[as.character(i),"bootCIdiff.hsic.gamma_pval_r3_lower_bound"]=res_tmp$boot.CI["diff.hsic.gamma_pval",1]
		results_mat[as.character(i),"bootCIdiff.hsic.gamma_pval_r3_upper_bound"]=res_tmp$boot.CI["diff.hsic.gamma_pval",2]
		results_mat[as.character(i),"bootCIdiff.hsic.gamma_pval_lower_bound_for_HT_bvnorm"]=res_tmp$boot.CI["diff.hsic.gamma_pval_HT_bvnorm",1]
		results_mat[as.character(i),"bootCIdiff.hsic.gamma_pval_upper_bound_for_HT_bvnorm"]=res_tmp$boot.CI["diff.hsic.gamma_pval_HT_bvnorm",2]
		results_mat[as.character(i),"estimated_proportion_of_BvNormSamples_be_obs_diff"]=res_tmp$estimated_proportion_bvnorm
		results_mat[as.character(i),"bootCIdiff.hsic.gamma_pval_lower_bound_for_HT_kmarg"]=res_tmp$boot.CI["diff.hsic.gamma_pval_HT_keep_marg",1]
		results_mat[as.character(i),"bootCIdiff.hsic.gamma_pval_upper_bound_for_HT_kmarg"]=res_tmp$boot.CI["diff.hsic.gamma_pval_HT_keep_marg",2]
		results_mat[as.character(i),"estimated_proportion_of_KMargSamples_be_obs_diff"]=res_tmp$estimated_proportion_kmarg
		####
		
   
    		print(paste("End of display for results from sampleset number ", as.character(i), sep=""))
		
	}

	return(list("folder_datasets"=folder_datasets,"used_datasets"=dataset_vec,"measurement_model"=measurement_model, "sim_datasets"=datasets_list,"results"=results_mat, "seed_info"=seed, "CI"=CI, "number_of_bootstrap_B"=B, "dcor.R"=dcor.R, "HSIC.B"=HSIC.B))

}

############################################################################################
############################################################################################



mes_model=
'#Observed indicators

	Vca1 ~ 0.96*Vca + norm
	Vca2 ~ 0.97*Vca + norm
	Vca3 ~ 0.98*Vca + norm

	Vcsq1 ~ 0.96*Vcsq + norm
	Vcsq2 ~ 0.97*Vcsq + norm
	Vcsq3 ~ 0.98*Vcsq + norm

	#where Vca and Vcsq are latent variables
'



metadta=read.table("/home/users/a/p/apollari/LC_analyses/tuebingen/pairmeta.txt")
folder_datasets="/home/users/a/p/apollari/LC_analyses/tuebingen/pairs/"

dataset_vec=1:(dim(metadta)[1])
dataset_vec=dataset_vec[-c(52,53,54,55,71,105,5:11,49:51)] #exclusion of multivariate datasets

resul=tuebingen_data(metadta,dataset_vec,folder_datasets, mes_model, CI=c(0.025,0.975), dcor.R=500, HSIC.B=500, B=1000, seed=4321)
resul
