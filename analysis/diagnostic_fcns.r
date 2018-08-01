#functions providng various diagnostics for evaluating model assumptions and stability and
#for helping with setting an appropriate random slopes structure for (G)LMMs
#written by Roger Mundry, last modified early 2016
ranef.diagn.plot<-function(model.res, QQ=F){
  old.par = par(no.readonly = TRUE)
	n.plots=sum(unlist(lapply(ranef(model.res), length)))
	x=ifelse(n.plots%%2==1,n.plots+1,n.plots)
	xmat=outer(1:x, 1:(1+x/2), "*")-n.plots
	colnames(xmat)=1:(1+x/2)
	rownames(xmat)=1:x
	xmat=as.data.frame(as.table(xmat))
	xmat=subset(xmat, as.numeric(as.character(xmat$Var1))>=as.numeric(as.character(xmat$Var2)) & xmat$Freq>=0)
	sum.diff=as.numeric(as.character(xmat$Var1))-as.numeric(as.character(xmat$Var2))+xmat$Freq
	xmat=xmat[sum.diff==min(sum.diff),]
	xmat=xmat[which.min(xmat$Freq),]
	par(mfrow=c(xmat$Var2, xmat$Var1))
	par(mar=c(rep(2, 3), 1))
	par(mgp=c(1, 0.5, 0))
	for(i in 1:length(ranef(model.res))){
		to.plot=ranef(model.res)[[i]]
		for(k in 1:ncol(to.plot)){
			if(QQ){
				qqnorm(to.plot[,k], main="", tcl=-0.25, xlab="", ylab="", pch=19, col=grey(0.5, alpha=0.75))
				qqline(to.plot[,k])
			}else{
				hist(to.plot[,k], main="", tcl=-0.25, xlab="", ylab="")
			}
			mtext(text=paste(c(names(ranef(model.res)[i]), colnames(to.plot)[k]), collapse=", "), side=3)
		}
	}
	par(old.par)
}

how.many.uniques.aao<-function(fe, re, data){
	print("you should consider using the function fe.re.tab instead")
	xx=setdiff(c(fe, re), names(data))
	if(length(xx)>0){
		stop(paste(c("Error: ", paste(xx,  collapse=", "), " are not in data"), collapse=""))
	}
	data=droplevels(as.data.frame(na.omit(data[, c(fe, re)])))
	to.do=data.frame(expand.grid(re, fe))
	names(to.do)=c("re", "fe")
	to.do$re=as.character(to.do$re)
	to.do$fe=as.character(to.do$fe)
	#browser()
	res.detailed=lapply(1:nrow(to.do), function(xrow){
		table(data[,to.do$re[xrow]], data[,to.do$fe[xrow]])
	})
	res.summary=lapply(res.detailed, function(xtab){
		table(apply(xtab>0, 1, sum))
	})
	xnames=paste(to.do$fe, to.do$re, sep="_within_")
	names(res.detailed)=xnames
	names(res.summary)=xnames
	return(list(detailed=res.detailed, summary=res.summary))
}

diagnostics.plot<-function(mod.res, col=grey(level=0.25, alpha=0.5)){
  old.par = par(no.readonly = TRUE)
  par(mfrow=c(2, 2))
  par(mar=c(3, 3, 1, 0.5))
  hist(residuals(mod.res), probability=T, xlab="", ylab="", main="")
  mtext(text="histogram of residuals", side=3, line=0)
  x=seq(min(residuals(mod.res)), max(residuals(mod.res)), length.out=100)
  lines(x, dnorm(x, mean=0, sd=sd(residuals(mod.res))))
  qqnorm(residuals(mod.res), main="", pch=19)
  qqline(residuals(mod.res))
  mtext(text="qq-plot of residuals", side=3, line=0)
  plot(fitted(mod.res), residuals(mod.res), pch=19, col=col)
  abline(h=0, lty=2)
  mtext(text="residuals against fitted values", side=3, line=0)
  par(old.par)
}

lev.thresh<-function(model.res){
	k=length(coefficients(model.res))
	n=length(residuals(model.res))
 return(2*(k+1)/n)
}

overdisp.test<-function(x){
  pr=residuals(x, type ="pearson")
  sum.dp=sum(pr^2)
  if(class(x)[[1]]=="lmerMod"){
    xdf=length(residuals(x))-length(fixef(x))
  }else if(class(x)[[1]]=="glmerMod"){
    xdf=length(residuals(x))-length(fixef(x))
  }else{
    xdf=length(residuals(x))-length(x$coefficients)
  }
  return(data.frame(chisq=sum.dp, df=xdf, P=1-pchisq(sum.dp, xdf), dispersion.parameter=sum.dp/xdf))
}

how.many.uniques<-function(xfac, xcov){
	print("you should consider using the function fe.re.tab instead")
  ii.data=data.frame(xfac=factor(xfac), xcov=xcov)
  ii.data=data.frame(na.omit(ii.data))
  ires=tapply(ii.data$xcov, ii.data$xfac, function(ii){
    length(unique(ii))
  })
  return(ires)
}

fe.re.tab<-function(fe.model, re, other.vars=NULL, data, treat.covs.as.factors=F){
	#function helping in determinining which random slopes are needed
	#last updated: 2015, Nov 25
	#input/arguments:
		#data: a dataframe with all relevant variables (including the response)
		#fe.model: character; the model wrt the fixed effect (including the response); e.g., "r~f1*c*f2"
		#re: character; either a vector with the names of the random effects (e.g., c("re1", "re2")) or a random intercepts expression (e.g., "(1|re1)+(1|re2)")
		#other.vars: character, optional; a vector with the names of variables which are to be kept in the data considered and returned
		#treat.covs.as.factors: logical, when set to TRUE covariates will be treated like factors (see value/summary for details)
	#value: list with the following entries:
		#detailed: list with cross-tabulations ffor each combination of (main) fixed and random effect 
		#summary: list tables...
			#telling for each combination of (main) fixed and random effect...
				#the number of levels of the random effect with a given number of unique values of the fixed effect (in case of a covariate)
				#the number of levels of the random effect with a given number of levels of the fixed effect for which at least two cases do exist (in case of a factor)
			#telling for each combination of interaction and random effect...
				#the combination of the above two informations, i.e., the number of individuals with a given number of unique values of the covariate
					#and a given number of factor levels for which at least two cases exist
		#data: data frame containing all relevant variables (i.e., response, fixed and random effects as well as those indicated in other.vars (e.g., offset terms)
			#also includes columns for dummy variables coding the levels (except the reference level) of all factors
		#pot.terms: length one vector comprising the model wrt the random slopes (for all combinations of fixed and random effects, i.e., most likely some will need to be omitted)
			#note that this comprises only the random slopes but not the correlations between random slopes and intercepts not the random intercept itself
		#pot.terms.with.corr: length one vector comprising the model wrt the random slopes (for all combinations of fixed and random effects, i.e., most likely some will need to be omitted)
			#note that this comprises the random slopes and intercepts  and also all correlations among them
	if(sum(grepl(x=re, pattern="", fixed=T))>0 & length(re)==1){#if random effects are handed over as formula
		re=gsub(x=re, pattern="(1|", replacement="", fixed=T)
		re=gsub(x=re, pattern=")", replacement="", fixed=T)
		re=gsub(x=re, pattern=" ", replacement="", fixed=T)
		re=unlist(strsplit(re, split="+", fixed=T))
	}
	fe.model=gsub(x=fe.model, pattern=" ", replacement="", fixed=T)#remove spaces
	model.terms=attr(terms(as.formula(fe.model)), "term.labels")#get individual terms from fixed effects model
	fe.me=model.terms[!grepl(x=model.terms, pattern=":", fixed=T)]#remove interactions
	fe.me=fe.me[!grepl(x=fe.me, pattern="^", fixed=T)]#remove squares terms
	resp=unlist(strsplit(fe.model, split="~", fixed=T))[1]#determine response
	if(substr(resp, start=1, stop=6)=="cbind("){
		resp=gsub(x=resp, pattern="cbind(", replacement="", fixed=T)
		resp=gsub(x=resp, pattern=")", replacement="", fixed=T)
		resp=gsub(x=resp, pattern=" ", replacement="", fixed=T)
		resp=unlist(strsplit(resp, split=",", fixed=T))
	}
	xx=setdiff(c(resp, fe.me, re, other.vars), names(data))
	if(length(xx)>0){
		stop(paste(c("error: preditor(s) missing in the data is/are ", paste(xx, collapse=", ")), collapse=""))
	}
	data=droplevels(as.data.frame(na.omit(data[, c(resp, fe.me, re, other.vars)])))#keep complete data wrt all relevant variables
	model.terms=model.terms[!grepl(x=model.terms, pattern="^", fixed=T)]#remove squares terms
	modes=rep(NA, length(model.terms))#initialize vector storing whether PVs are factors or not
	effect=c(rep("main", length(fe.me)), rep("int", length(model.terms)-length(fe.me)))#create vector telling for each model term whether it is a main effect of not
	for(i in 1:ncol(data)){#for all columns in data
		if(is.character(data[, i])){data[, i]=as.factor(data[, i])}#turn character column into factor
		modes[i]=class(data[, i])#and determine its class
	}
	names(modes)=names(data)#name 'm
	to.do=data.frame(expand.grid(re=re, fe=fe.me))#create data frame with one column for each combination of fixed main and random effect
	to.do$re=as.character(to.do$re)#reformat to character (for later addressing by name)
	to.do$fe=as.character(to.do$fe)#reformat to character (for later addressing by name)
	res.detailed=lapply(1:nrow(to.do), function(xrow){#create detailed results by lapply-ing over the rows of to.do
			table(data[,to.do$re[xrow]], data[,to.do$fe[xrow]])#cross tabulate the respective fixed and random effect
	})
	names(res.detailed)=paste(to.do$fe, to.do$re, sep="_within_")#name it
	res.summary=lapply(1:nrow(to.do), function(xrow){#begin with creating summary results by lapply-ing over the rows of to.do
		if(modes[to.do[xrow, "fe"]]!="factor" & !treat.covs.as.factors){#if fixed effect is not a factor
			table(apply(res.detailed[[xrow]]>0, 1, sum))#determine number of levels of the random effect per number of unique cases of the fixed effect
		}else{#if fixed effect is a factor
			table(apply(res.detailed[[xrow]]>1, 1, sum))#determine number of levels of the random effect per number of unique cases of the fixed effect with at least two observations
		}
	})
	names(res.summary)=paste(to.do$fe, to.do$re, sep="_within_")#name it
	#append 'factor' or 'covariate' to the names:
	names(res.summary)[modes[to.do$fe]=="factor"]=paste(names(res.summary)[modes[to.do$fe]=="factor"], "(factor)", sep=" ")
	names(res.summary)[!modes[to.do$fe]=="factor"]=paste(names(res.summary)[!modes[to.do$fe]=="factor"], "(covariate)", sep=" ")
	to.do=data.frame(expand.grid(re=re, int=setdiff(model.terms, fe.me)))#create data frame with one row for each combination of interaction and random effect
	if(nrow(to.do)>0){
		to.do[, "int"]=as.character(to.do[, "int"])#reformat to character (for later addressing by name)
		to.do[, "re"]=as.character(to.do[, "re"])#reformat to character (for later addressing by name)
		to.add=lapply(1:nrow(to.do), function(xrow){#treat combinations of fixed and random effect by lapply-ing over the rows of to.do
			iterms=unlist(strsplit(to.do[xrow, "int"], split=":", fixed=T))#determine main effects involved in the interaction...
			imodes=modes[iterms]#... and their classes
			all.tabs=lapply(1:length(iterms), function(yrow){#for each term
				#if its not a factor, determine determine the number of unique cases of the fixed effect per level of the random effect
				#otherwise, determine per level of thee random effect the number of levels of the fixed effect with at least two observations
				apply(res.detailed[[paste(c(iterms[yrow], to.do[xrow, "re"]), collapse="_within_")]]>ifelse(imodes[iterms[yrow]]=="factor" | treat.covs.as.factors, 1, 0), 1, sum)
			})
			all.tabs=matrix(unlist(all.tabs), ncol=length(all.tabs), byrow=F)#reformat to matrix (one row per level of the random effect, one column per term)
			all.tabs=aggregate(1:nrow(all.tabs), c(data.frame(all.tabs)), length)#summarize it ((wrt the number of levels of the random effec per combination of values)
			colnames(all.tabs)=c(iterms, paste(c("n", to.do[xrow, "re"]), collapse="."))#name it
			return(all.tabs)
		})
		names(to.add)=gsub(x=paste(to.do$int, to.do$re, sep="_within_"), pattern=":", replacement="_", fixed=T)#name it
		res.summary=c(res.summary, to.add)#and append to.add to res.summary
	}
	#add columns with the dummy coded factor levels to data:
	to.code=fe.me[modes[fe.me]=="factor"]#determine factors among the fixed effects:
	if(length(to.code)>0){
		coded=lapply(to.code, function(xc){#for each factor
			lapply(levels(data[, xc])[-1], function(xl){#for all levels except the reference level
				as.numeric(data[, xc]==xl)#code it
			})
		})
		coded=matrix(unlist(coded), nrow=nrow(data), byrow=F)#reformat to matrix
		xnames=unlist(lapply(to.code, function(xc){#determine column names to be given to the matrix...
			paste(xc, levels(data[, xc])[-1], sep=".")
		}))
		colnames(coded)=xnames#... and use 'm
		data=data.frame(data, coded)#and append 'm to data
	}else{
		xnames=""
	}
	#now the model expression wrt random slopes; first no correllations between random slopes and intercepts:
	pot.terms=c(xnames, fe.me[modes[fe.me]!="factor"])#create vector with names of dummy variables and names of fixed main effects not being factors
	pot.terms=outer(pot.terms, re, Vectorize(function(x, y){#create matrix with separate model term for each combination of fixed and random effect
		paste(c("(0+", x, "|", y, ")"), collapse="")
	}))
	pot.terms=paste(c(t(pot.terms)), collapse="+")#put them all in a single entry
	#now the random slopes part of the model including random intercepts and slopes and their correlation:
	pot.terms.with.corr=paste(c(xnames, fe.me[modes[fe.me]!="factor"]), collapse="+")#get all fixed effects terms together...
	pot.terms.with.corr=paste(paste("(1+", pot.terms.with.corr, "|", re, ")", sep=""), collapse="+")#... and paste random effects, brackets and all that 
		#(and everything in a single entry)
	return(list(detailed=res.detailed, summary=res.summary, data=data, pot.terms=pot.terms, pot.terms.with.corr=pot.terms.with.corr))
}

write.fe.re.summary2file<-function(x, file){
	append=F
	for(i in 1:length(x)){
		write.table(x=names(x[i]), file=file, row.names=F, col.names=F, sep="\t", quote=F, append=append)
		append=T
		xx=x[[i]]
		if(length(dim(xx))==1){
			xx=matrix(xx, ncol=1)
			row.names(xx)=names(x[[i]])
			xx=c.tab(x=xx, add.hash=F)
			xx=matrix(xx, ncol=1)
		}else{
			xx=c.tab(x=xx, add.hash=F)
			xx=matrix(xx, ncol=1)
		}
		write.table(x=xx, file=file, row.names=F, col.names=F, sep="\t", quote=F, append=append)
		write.table(x="", file=file, row.names=F, col.names=F, sep="\t", quote=F, append=append)
	}
}

m.stab.plot<-function(est, lower=NULL, upper=NULL, xnames=NULL, col="black", center.at.null=F){
	#version from dec. 1 2015
	#function to plot model stability or CIs (not really supposed to be nice, but to give a quick overview/rapid diagnostic)
	#input:
		#either a three columns data frame or matrix (with rownames) handed over to argument est (first column needs to comprise the original estimate)
		#or several vectors handed over as follows:
			#est: numeric; estimated coefficients of the model
			#lower: numeric; lower limits of the estimates (either from model stability of from bootstrap)
			#upper: numeric; upper limits of the estimates (either from model stability of from bootstrap)
			#xnames: character; names to be depicted besides the error bars
		#col: character; name of the color with which data should be depicted
	if(ncol(est)==3){
		lower=est[, 2]
		upper=est[, 3]
		xnames=rownames(est)
		est=est[, 1]
		x.at=est
	}
	par(mar=c(3, 0.5, 0.5, 0.5), mgp=c(1, 0.4, 0), tcl=-0.2)
	plot(x=est, y=1:length(est), pch=18, xlab="estimate", ylab="", yaxt="n", xlim=range(c(lower, upper)), type="n", ylim=c(1, length(est)+1))
	abline(v=0, lty=3)
	if(center.at.null){
		x.at=rep(0, length(est))
		text(labels=xnames, x=x.at, y=(1:length(est))+0.3, cex=0.8)
	}else{
		text(labels=xnames, x=x.at, y=(1:length(est))+0.3, cex=0.8, pos=c(2, 4)[1+as.numeric(est<0)])
	}
	points(x=est, y=1:length(est), pch=18, col=col)
	segments(x0=lower, x1=upper, y0=1:length(est), y1=1:length(est), col=col)
}
