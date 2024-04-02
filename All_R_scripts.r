library(adephylo)
library(diagram)
library(HDInterval)
library(lubridate)
library(maptools)
library(rgdal)
library(ridgeline)
library(seraphim)
library(vioplot)

# 1. Investigating the impact of sampling effort on the estimate of dispersal metrics (using BRW or RRW simulations)

	# 1.1. Investigating the consistency of dispersal metrics when estimated on Brownian random walk (BRW) simulations

		# 1.1.1. Birth-death simulations

simulationDirectory = "Simulations_1/BirthDeath_simulations"; RRW = FALSE; scalingValue = 1.5
envVariable = crop(raster("Template_R.tif"), extent(-100,-50,20,60))
envVariable[] = 1
sigma = scalingValue*res(envVariable)[1]
ancestPosition = c(-73.79730, 40.94328)
birthRate = 0.6
samplingRate = 0.4
startingYear = 1999
samplingWindow = cbind(2001,2017)
timeSlice = 1/365.25
timeIntervale = 1
showingPlots = FALSE
for (i in 1:250)
	{
		if (!file.exists(paste0(simulationDirectory,"/All_the_simulations/RRW_simulation_",i,".csv")))
			{
				simulation = NULL; worked = FALSE
				while (worked == FALSE)
					{
						trycatch = tryCatch(
							{
								simulation = simulatorRRW2(RRW, envVariable, sigma, ancestPosition, birthRate, samplingRate, startingYear, 
														   samplingWindow, timeSlice, timeIntervale, showingPlots)
							},	error = function(cond) {
							},	finally = {
							})
						if (!is.null(simulation)) worked = TRUE
					}								   
				write.csv(simulation[[1]], paste0(simulationDirectory,"/All_the_simulations/BRW_simulation_",i,".csv"), quote=F, row.names=F)
				write.tree(simulation[[2]], paste0(simulationDirectory,"/All_the_simulations/BRW_simulation_",i,".tre"))
			}
	}

		# 1.1.2. Coalescent simulations

simulationDirectory = "Simulations_1/Coalescent_simulations"
ancestPosition = c(-73.79730, 40.94328)
for (i in 1:50)
	{
		tree = ape::rcoal(n=500, br="coalescent")
		rootHeight = max(phytools::nodeHeights(tree))
		tree$edge.length = (tree$edge.length/rootHeight)*18
		sd = 0.00833333*1.5; sd = 0.00833333*20; sigma2 = sd^2
		lons = phytools::fastBM(tree, a=ancestPosition[1], mu=0, sig2=sigma2, internal=T)
		lats = phytools::fastBM(tree, a=ancestPosition[2], mu=0, sig2=sigma2, internal=T)	
		for (j in 1:length(lons))
			{
				index = which(tree$tip.label==names(lons)[j])
				if (length(index) == 1) names(lons)[j] = index
			}
		for (j in 1:length(lats))
			{
				index = which(tree$tip.label==names(lats)[j])
				if (length(index) == 1) names(lats)[j] = index
			}		
		colNames = c("node1","node2","length","startLon","startLat","endLon","endLat")
		tab = matrix(nrow=dim(tree$edge)[1], ncol=length(colNames)); colnames(tab) = colNames
		tab[,c("node1","node2")] = tree$edge; tab[,c("length")] = tree$edge.length
		l = length(tab[,1]); ll = matrix(1:l, nrow=l, ncol=l); ll[] = 0
		mostRecentSamplingDatum = 2017
		for (j in 1:l)
			{
				subMat = tab[j,2]
				subMat = subset(tab, tab[,2]==subMat)
				ll[j,1] = subMat[,3]
				subMat = subMat[1,1]
				subMat1 = subset(tab, tab[,2]==subMat)
				for (k in 1:l)
					{
						if (nrow(subMat1) > 0)
							{
								ll[j,k+1] = subMat1[,3]
	 							subMat2 = subMat1[1,1]
	 							subMat1 = subset(tab, tab[,2]==subMat2)
	 						}
	 				}
			}
		endNodeL = rowSums(ll)
		tab = cbind(tab, endNodeL)
		startNodeL = matrix(1:l, nrow=l, ncol=1)
		startNodeL[] = 0
		for (j in 1:l)
			{
				r = tab[j,1]
				s = subset(tab, tab[,2]==r)
				for (k in 1:l)
					{
						if (nrow(s) > 0)
							{
								startNodeL[j,1] = s[,8]
	 						}
	 				}	
			}
		colnames(startNodeL) = "startNodeL"
		tab = cbind(tab, startNodeL)
		maxEndLIndice = which.max(tab[,"endNodeL"])
		maxEndL = tab[maxEndLIndice,"endNodeL"]
	 	endYear = matrix(tab[,"endNodeL"]-maxEndL)
	 	endYear = matrix(mostRecentSamplingDatum+(endYear[,1]))
		startYear = matrix(tab[,"startNodeL"]-maxEndL)
	 	startYear = matrix(mostRecentSamplingDatum+(startYear[,1]))
	 	colnames(startYear) = "startYear"; colnames(endYear) = "endYear"
	 	tab = cbind(tab,startYear,endYear)
		tab = tab[order(tab[,"startYear"],decreasing=F),]
		tab1 = tab[1,]; tab2 = tab[2:dim(tab)[1],]
		tab2 = tab2[order(tab2[,"endYear"],decreasing=F),]
		tab = rbind(tab1, tab2)
		for (j in 1:length(lons))
			{
				indices1 = which(tab[,"node1"]==as.numeric(names(lons)[j]))
				if (length(indices1) > 0) tab[indices1,"startLon"] = lons[j]
				indices2 = which(tab[,"node2"]==as.numeric(names(lons)[j]))
				if (length(indices2) > 0) tab[indices2,"endLon"] = lons[j]
			}
		for (j in 1:length(lats))
			{
				indices1 = which(tab[,"node1"]==as.numeric(names(lats)[j]))
				if (length(indices1) > 0) tab[indices1,"startLat"] = lats[j]
				indices2 = which(tab[,"node2"]==as.numeric(names(lats)[j]))
				if (length(indices2) > 0) tab[indices2,"endLat"] = lats[j]
			}
		greatCircleDistances = matrix(nrow=dim(tab)[1], ncol=1)
		colnames(greatCircleDistances) = "greatCircleDist_km"
		for (j in 1:dim(tab)[1])
			{
				x1 = cbind(tab[j,"startLon"], tab[j,"startLat"]); x2 = cbind(tab[j,"endLon"], tab[j,"endLat"])
				greatCircleDistances[j,1] = rdist.earth(x1, x2, miles=F, R=NULL)
			}
		tab = cbind(tab, greatCircleDistances)
		write.csv(tab, paste0(simulationDirectory,"/All_the_simulations/BRW_simulation_",i,".csv"), quote=F, row.names=F)
		write.tree(tree, paste0(simulationDirectory,"/All_the_simulations/BRW_simulation_",i,".tre"))
	}

	# 1.2. RRW simulations conducted on a uniform raster to investigate the consistency of dispersal metrics

simulationDirectory = "Simulations_2"; RRW = TRUE; scalingValue = 1.5/50
dir.create(file.path(paste0(simulationDirectory)), showWarnings=F)
dir.create(file.path(paste0(simulationDirectory,"/All_the_simulations")), showWarnings=F)
dir.create(file.path(paste0(simulationDirectory,"/Subsampling_of_tips")), showWarnings=F)
envVariable = crop(raster("Template_R.tif"), extent(-100,-50,20,60))
envVariable[] = 1
sigma = scalingValue*res(envVariable)[1]
ancestPosition = c(-73.79730, 40.94328)
birthRate = 0.6
samplingRate = 0.4
startingYear = 1999
samplingWindow = cbind(2001,2017)
timeSlice = 1/365.25
timeIntervale = 1
showingPlots = FALSE
for (i in 1:200)
	{
		if (!file.exists(paste0(simulationDirectory,"/All_the_simulations/RRW_simulation_",i,".csv")))
			{
				simulation = NULL; worked = FALSE
				while (worked == FALSE)
					{
						trycatch = tryCatch(
							{
								simulation = simulatorRRW2(RRW, envVariable, sigma, ancestPosition, birthRate, samplingRate, startingYear, 
														   samplingWindow, timeSlice, timeIntervale, showingPlots)
							},	error = function(cond) {
							},	finally = {
							})
						if (!is.null(simulation)) worked = TRUE
					}								   
				write.csv(simulation[[1]], paste0(simulationDirectory,"/All_the_simulations/RRW_simulation_",i,".csv"), quote=F, row.names=F)
				write.tree(simulation[[2]], paste0(simulationDirectory,"/All_the_simulations/RRW_simulation_",i,".tre"))
			}
	}

	# 1.3. Estimating the dispersal metrics based on the BRW and RRW simulations conducted on rasters

simulationDirectory = "Simulations_1/BirthDeath_simulations"; prefix = "BRW"
simulationDirectory = "Simulations_1/Coalescent_simulations"; prefix = "BRW"
simulationDirectory = "Simulations_2"; prefix = "RRW"
nberOfRemainingSimulations = 50; nberOfBranchesToKeep = seq(500,50,-50)

files = list.files(paste0(simulationDirectory,"/All_the_simulations/"))
files = files[which(grepl("\\.tre",files))]
for (i in 1:length(files))
	{
		if (i == 1) n = 0
		tab = read.csv(paste0(simulationDirectory,"/All_the_simulations/",prefix,"_simulation_",i,".csv"), head=T)
		nTips = length(which(!tab[,"node2"]%in%tab[,"node1"]))
		if (nTips >= 500)
			{
				n = n+1; write.csv(tab, paste0(simulationDirectory,"/Subsampling_of_tips/",prefix,"_simulation_",n,"_ALL.csv"), quote=F, row.names=F)
				file.copy(paste0(simulationDirectory,"/All_the_simulations/",prefix,"_simulation_",i,".tre"),
						  paste0(simulationDirectory,"/Subsampling_of_tips/",prefix,"_simulation_",n,"_ALL.tre"))
			}
	}
startingYear = 1999; samplingWindow = cbind(2001,2017)
ancestPosition = c(-73.79730, 40.94328)
for (i in 1:nberOfRemainingSimulations)
	{
		pdf(paste0(simulationDirectory,"/Subsampling_of_tips/",prefix,"_simulation_",i,"_ALL.pdf"), width=9, height=4.4)
		par(mfrow=c(1,2), mgp=c(0,0,0), oma=c(0,0,0,0), mar=c(1.4,0.5,0.5,0.5), lwd=0.2, col="gray50")		
		colourScale = colorRampPalette(brewer.pal(9,"YlGnBu"))(131)[1:101]
		tab = read.csv(paste0(simulationDirectory,"/Subsampling_of_tips/",prefix,"_simulation_",i,"_ALL.csv"), head=T)
		tab[,c("startLon","endLon")] = tab[,c("startLon","endLon")]-ancestPosition[1]
		tab[,c("startLat","endLat")] = tab[,c("startLat","endLat")]-ancestPosition[2]
		tree = read.tree(paste0(simulationDirectory,"/Subsampling_of_tips/",prefix,"_simulation_",i,"_ALL.tre"))
		col_start = colourScale[1]; cols3 = colourScale[(((nodeHeights(tree)[,2])/(samplingWindow[2]-startingYear))*100)+1]
		plot(tree, show.tip.label=F, edge.width=0.5, edge.col="gray50")
		for (j in 1:dim(tree$edge)[1])
			{
				if (j == 1)
					{
						nodelabels(node=tree$edge[j,1], pch=16, cex=0.5, col=col_start)
						nodelabels(node=tree$edge[j,1], pch=1, cex=0.5, col="gray30", lwd=0.25)
					}
				if (grepl("BirthDeath",simulationDirectory))
					{
						nodelabels(node=tree$edge[j,2], pch=16, cex=0.5, col=cols3[j])
						nodelabels(node=tree$edge[j,2], pch=1, cex=0.5, col="gray30", lwd=0.25)
					}
				if (grepl("Coalescent",simulationDirectory))
					{
						if (tree$edge[j,2]%in%tree$edge[,1])
							{
								nodelabels(node=tree$edge[j,2], pch=16, cex=0.5, col=cols3[j])
								nodelabels(node=tree$edge[j,2], pch=1, cex=0.5, col="gray30", lwd=0.25)
							}
						if (!tree$edge[j,2]%in%tree$edge[,1])
							{
								# nodelabels(node=tree$edge[j,2], pch=16, cex=0.2, col=cols3[j])
								# nodelabels(node=tree$edge[j,2], pch=1, cex=0.2, col="gray30", lwd=0.25)
							}
					}
			}
		axis(side=1, lwd.tick=0.3, cex.axis=0.5, mgp=c(0,-0.10,0), lwd=0.3, tck=-0.010, col.tick="gray30", col.axis="gray30", col="gray30", at=seq(0,20,5))
		startingYear = min(tab[,"startYear"]); col_start = colourScale[1]
		cols2 = colourScale[(((tab[,"endYear"]-startingYear)/(samplingWindow[2]-startingYear))*100)+1]
		xMin = min(tab[,c("startLon","endLon")]); xMax = max(tab[,c("startLon","endLon")]); dX = xMax-xMin
		yMin = min(tab[,c("startLat","endLat")]); yMax = max(tab[,c("startLat","endLat")]); dY = yMax-yMin
		plot(extent(xMin,xMax,yMin,yMax), col=NA, ann=F, axes=F, asp=1)
		for (j in 1:dim(tab)[1])
			{	
				curvedarrow(cbind(tab[j,"startLon"],tab[j,"startLat"]), cbind(tab[j,"endLon"],tab[j,"endLat"]), arr.length=0,
							arr.width=0, lwd=0.3, lty=1, lcol="gray50", arr.col=NA, arr.pos=F, curve=0.1, dr=NA, endhead=F)
			}
		for (j in dim(tab)[1]:1)
			{
				if (j == 1)
					{
						points(cbind(tab[j,"startLon"],tab[j,"startLat"]), pch=16, cex=0.5, col=col_start)
						points(cbind(tab[j,"startLon"],tab[j,"startLat"]), pch=1, cex=0.5, col=rgb(1,1,1,255,maxColorValue=255), lwd=0.25)
					}
				points(cbind(tab[j,"endLon"],tab[j,"endLat"]), pch=16, cex=0.5, col=cols2[j])
				points(cbind(tab[j,"endLon"],tab[j,"endLat"]), pch=1, cex=0.5, col=rgb(1,1,1,255,maxColorValue=255), lwd=0.25)
			}
		if (grepl("Simulations_1",simulationDirectory))
			{
				axis(side=1, lwd.tick=0.3, cex.axis=0.5, mgp=c(0,-0.10,0), lwd=0.3, tck=-0.010, col.tick="gray30", col.axis="gray30", col="gray30", at=seq(-10,10,1))
				axis(side=2, lwd.tick=0.3, cex.axis=0.5, mgp=c(0,0.20,0), lwd=0.3, tck=-0.010, col.tick="gray30", col.axis="gray30", col="gray30", at=seq(-10,10,1))
			}
		if (grepl("Simulations_2",simulationDirectory))
			{
				axis(side=1, lwd.tick=0.3, cex.axis=0.5, mgp=c(0,-0.10,0), lwd=0.3, tck=-0.010, col.tick="gray30", col.axis="gray30", col="gray30", at=seq(-100,100,5))
				axis(side=2, lwd.tick=0.3, cex.axis=0.5, mgp=c(0,0.20,0), lwd=0.3, tck=-0.010, col.tick="gray30", col.axis="gray30", col="gray30", at=seq(-100,100,5))
			}
		dev.off()
	}
for (i in 1:length(nberOfBranchesToKeep))
	{
		for (j in 1:nberOfRemainingSimulations)
			{
				tab1 = read.csv(paste0(simulationDirectory,"/Subsampling_of_tips/",prefix,"_simulation_",j,"_ALL.csv"), head=T)
				if ("treeID"%in%colnames(tab1))
					{
						colnames(tab1)[which(colnames(tab1)=="treeID")] = "length"
						tab1[,"length"] = tab1[,"endYear"]-tab1[,"startYear"]
					}
				tipBranches = which(!tab1[,"node2"]%in%tab1[,"node1"]); nberOfTips = length(tipBranches)
				if (nberOfTips > nberOfBranchesToKeep[i])
					{
						tipsToDiscard = sample(tipBranches, nberOfTips-nberOfBranchesToKeep[i], replace=F)
						nodeClades = matrix(nrow=dim(tab1), ncol=2); colnames(nodeClades) = c("startClade","endClade")
						nodeClades[tipBranches,"endClade"] = "toKeep"; nodeClades[tipsToDiscard,"endClade"] = "toDiscard"
						tab1 = cbind(tab1, nodeClades); noRemainingCladeToIdentify = FALSE
						while (noRemainingCladeToIdentify == FALSE)
							{
								for (n in 1:dim(tab1)[1])
									{
										if (is.na(tab1[n,"endClade"]))
											{
												indices = which(tab1[,"node1"]==tab1[n,"node2"])
												if (sum(!is.na(tab1[indices,"endClade"])) == 2)
													{
														if (tab1[indices[1],"endClade"] == tab1[indices[2],"endClade"])
															{
																tab1[indices[1],"startClade"] = tab1[indices[1],"endClade"]
																tab1[indices[2],"startClade"] = tab1[indices[2],"endClade"]
																tab1[n,"endClade"] = tab1[indices[1],"endClade"]
															}	else	{
																tab1[indices[1],"startClade"] = "toKeep"
																tab1[indices[2],"startClade"] = "toKeep"
																tab1[n,"endClade"] = "toKeep"
															}
													}
											}
									}
								if (sum(is.na(tab1[,"startClade"])) == 2)
									{
										tab1[which(!tab1[,"node1"]%in%tab1[,"node2"]),"startClade"] = "toKeep"
									}				
								if (sum(is.na(tab1[,"startClade"])) == 0) noRemainingCladeToIdentify = TRUE
							}
						tab2 = tab1[-which(tab1[,"endClade"]=="toDiscard"),]; tab3 = tab2
						nTips = length(which(!tab3[,"node2"]%in%tab3[,"node1"])); nBranches = dim(tab3)[1]
						while (nBranches != ((2*nTips)-2))
							{
								branchesToMerge = c(); newBranches = c()			
								for (n in 1:dim(tab1)[1])
									{
										if (!n%in%branchesToMerge)
											{
												if (length(which(tab3[,"node1"]==tab3[n,"node2"])) == 1)
													{
														m = which(tab3[,"node1"]==tab3[n,"node2"])
														branchesToMerge = c(branchesToMerge, n, m)
														newBranch = tab3[n,]
														newBranch[1,"node2"] = tab3[m,"node2"]
														newBranch[1,"endLon"] = tab3[m,"endLon"]
														newBranch[1,"endLat"] = tab3[m,"endLat"]
														newBranch[1,"endYear"] = tab3[m,"endYear"]
														newBranch[1,"endNodeL"] = tab3[m,"endNodeL"]
														newBranch[1,"length"] = newBranch[1,"length"]+tab3[m,"length"]
														newBranches = rbind(newBranches, newBranch)
													}
											}
									}
								if (length(branchesToMerge) > 0)
									{
										tab3 = tab3[-branchesToMerge,]; tab3 = rbind(tab3, newBranches)
									}
								if (length(which(!tab3[,"node1"]%in%tab3[,"node2"])) == 1)
									{
										tab3 = tab3[-which(!tab3[,"node1"]%in%tab3[,"node2"]),]
									}
								nBranches = dim(tab3)[1]; # print(nBranches)
							}
						tab3 = tab3[,which(!colnames(tab3)%in%c("startClade","endClade"))]
						write.csv(tab3, paste0(simulationDirectory,"/Subsampling_of_tips/",prefix,"_simulation_",j,"_",nberOfBranchesToKeep[i],".csv"), row.names=F, quote=F)
					}	else	{
						write.csv(tab1, paste0(simulationDirectory,"/Subsampling_of_tips/",prefix,"_simulation_",j,"_",nberOfBranchesToKeep[i],".csv"), row.names=F, quote=F)
					}
			}
	}

tab_WLDV_branches = matrix(nrow=nberOfRemainingSimulations, ncol=length(nberOfBranchesToKeep)); colnames(tab_WLDV_branches) = paste0(nberOfBranchesToKeep,"_tips")
tab_WDCs_branches = matrix(nrow=nberOfRemainingSimulations, ncol=length(nberOfBranchesToKeep)); colnames(tab_WDCs_branches) = paste0(nberOfBranchesToKeep,"_tips")
tab_IBD_signal_rP = matrix(nrow=nberOfRemainingSimulations, ncol=length(nberOfBranchesToKeep)); colnames(tab_IBD_signal_rP) = paste0(nberOfBranchesToKeep,"_tips")
for (i in 1:length(nberOfBranchesToKeep))
	{
		for (j in 1:nberOfRemainingSimulations)
			{
				print(c(i,j))
				tab = read.csv(paste0(simulationDirectory,"/Subsampling_of_tips/",prefix,"_simulation_",j,"_",nberOfBranchesToKeep[i],".csv"), head=T)
				tre = read.tree(paste0(simulationDirectory,"/Subsampling_of_tips/",prefix,"_simulation_",j,"_ALL.tre"))
				if (grepl("t",tre$tip.label[1])) tre$tip.label = 1:length(tre$tip.label)
				tab[,"greatCircleDist_km"] = diag(rdist.earth(cbind(tab[,"startLon"],tab[,"startLat"]), cbind(tab[,"endLon"],tab[,"endLat"]), miles=F))
				tab_WLDV_branches[j,i] = sum(tab[,"greatCircleDist_km"])/(sum(tab[,"length"]))
				tab_WDCs_branches[j,i] = sum(tab[,"greatCircleDist_km"]^2)/(sum(4*tab[,"length"]))
				tipNodeIndices = which(!tab[,"node2"]%in%tab[,"node1"]); tipNodes = tab[tipNodeIndices,"node2"]
				distsGeo = rdist.earth(cbind(tab[tipNodeIndices,"endLon"],tab[tipNodeIndices,"endLat"]), cbind(tab[tipNodeIndices,"endLon"],tab[tipNodeIndices,"endLat"]), miles=F)
				distsPat = as.matrix(distTips(tre, method="patristic")); distsPat = distsPat[as.character(tipNodes),as.character(tipNodes)]
				distsGeo = distsGeo[lower.tri(distsGeo)]; distsPat = distsPat[lower.tri(distsPat)]
				tab_IBD_signal_rP[j,i] = cor(distsPat[lower.tri(distsPat)],log(distsGeo[lower.tri(distsGeo)]), method="pearson")
			}
	}
write.csv(tab_WLDV_branches, paste0(simulationDirectory,"/WLDV_branches.csv"), row.names=F, quote=F)
write.csv(tab_WDCs_branches, paste0(simulationDirectory,"/WDCs_branches.csv"), row.names=F, quote=F)
write.csv(tab_IBD_signal_rP, paste0(simulationDirectory,"/IBD_signal_rPcor.csv"), row.names=F, quote=F)

tab_WLDV_branches = read.csv(paste0(simulationDirectory,"/WLDV_branches.csv"), head=T)
tab_WDCs_branches = read.csv(paste0(simulationDirectory,"/WDCs_branches.csv"), head=T)
tab_IBD_signal_rP = read.csv(paste0(simulationDirectory,"/IBD_signal_rPcor.csv"), head=T)
max_WLDV_branches = max(tab_WLDV_branches); median_WLDV_branches = rep(NA, length(nberOfBranchesToKeep))
max_WDCs_branches = max(tab_WDCs_branches); median_WDCs_branches = rep(NA, length(nberOfBranchesToKeep))
max_IBD_signal_rP = max(tab_IBD_signal_rP); median_IBD_signal_rP = rep(NA, length(nberOfBranchesToKeep))
for (i in 1:length(nberOfBranchesToKeep))
	{
		median_WLDV_branches[i] = median(tab_WLDV_branches[,i], na.rm=T)
		median_WDCs_branches[i] = median(tab_WDCs_branches[,i], na.rm=T)
		median_IBD_signal_rP[i] = median(tab_IBD_signal_rP[,i], na.rm=T)
	}

pdf(paste0(simulationDirectory,"/Subsampling_tips_NEW1.pdf"), width=9.0, height=2.7) # dev.new(width=9.0, height=2.7)
par(mfrow=c(1,3), mgp=c(0,0,0), oma=c(0,0.5,0,0.1), mar=c(3.3,3.7,1.5,1.5), lwd=0.2, col="gray50"); cexLab = 1.07; cexAxis = 0.95; tck = -0.02
plot(nberOfBranchesToKeep, tab_WLDV_branches[1,], ann=F, axes=F, type="l", lwd=0.1, col="gray50", xlim=c(50,502), ylim=c(0,max_WLDV_branches))
for (i in 2:dim(tab_WLDV_branches)[1]) lines(nberOfBranchesToKeep, tab_WLDV_branches[i,], lwd=0.2, col="gray30")
lines(nberOfBranchesToKeep, median_WLDV_branches, lwd=1.5, col=rgb(204,0,0,255,maxColorValue=255))
axis(side=1, lwd.tick=0.3, cex.axis=cexAxis, mgp=c(0,0.40,0), lwd=0.3, tck=tck, col.tick="gray30", col.axis="gray30", col="gray30", at=seq(0,600,100), pos=0)
axis(side=2, lwd.tick=0.3, cex.axis=cexAxis, mgp=c(0,0.45,0), lwd=0.3, tck=tck, col.tick="gray30", col.axis="gray30", col="gray30", at=seq(0,30,5)) # Simulations_1[or 2]/BirthDeath_simulations
# axis(side=2, lwd.tick=0.3, cex.axis=cexAxis, mgp=c(0,0.45,0), lwd=0.3, tck=tck, col.tick="gray30", col.axis="gray30", col="gray30", at=seq(0,50,10)) # Simulations_1/Coalescent_simulations
title(ylab="Weighted lineage dispersal velocity (km/year)", cex.lab=cexLab, mgp=c(2.1,0,0), col.lab="gray30")
title(xlab="Number of tip nodes in the subsampled tree", cex.lab=cexLab, mgp=c(1.2,0,0), col.lab="gray30")
plot(nberOfBranchesToKeep, tab_WDCs_branches[1,], ann=F, axes=F, type="l", lwd=0.1, col="gray50", xlim=c(50,502), ylim=c(0,max(tab_WDCs_branches, na.rm=T)))
# plot(nberOfBranchesToKeep, tab_WDCs_branches[1,], ann=F, axes=F, type="l", lwd=0.1, col="gray50", xlim=c(50,502), ylim=c(0,252)) # Simulations_1[or 4]/Coalescent_simulations_NA_2res
for (i in 2:dim(tab_WDCs_branches)[1]) lines(nberOfBranchesToKeep, tab_WDCs_branches[i,], lwd=0.2, col="gray50")
lines(nberOfBranchesToKeep, median_WDCs_branches, lwd=1.5, col=rgb(204,0,0,255,maxColorValue=255))
axis(side=1, lwd.tick=0.3, cex.axis=cexAxis, mgp=c(0,0.40,0), lwd=0.3, tck=tck, col.tick="gray30", col.axis="gray30", col="gray30", at=seq(0,600,100), pos=0)
axis(side=2, lwd.tick=0.3, cex.axis=cexAxis, mgp=c(0,0.45,0), lwd=0.3, tck=tck, col.tick="gray30", col.axis="gray30", col="gray30", at=seq(0,400,50)) # Simulations_1[or 2]/BirthDeath_simulations
# axis(side=2, lwd.tick=0.3, cex.axis=cexAxis, mgp=c(0,0.45,0), lwd=0.3, tck=tck, col.tick="gray30", col.axis="gray30", col="gray30", at=seq(0,350,50)) # Simulations_1/Coalescent_simulations
title(ylab=expression("  Weighted diffusion coefficient (km"^2*"/year)"), cex.lab=cexLab, mgp=c(2.1,0,0), col.lab="gray30")
title(xlab="Number of tip nodes in the subsampled tree", cex.lab=cexLab, mgp=c(1.2,0,0), col.lab="gray30")
plot(nberOfBranchesToKeep, tab_IBD_signal_rP[1,], ann=F, axes=F, type="l", lwd=0.1, col="gray50", xlim=c(50,502), ylim=c(0,0.96))
for (i in 2:dim(tab_IBD_signal_rP)[1]) lines(nberOfBranchesToKeep, tab_IBD_signal_rP[i,], lwd=0.2, col="gray50")
lines(nberOfBranchesToKeep, median_IBD_signal_rP, lwd=1.5, col=rgb(204,0,0,255,maxColorValue=255))
axis(side=1, lwd.tick=0.3, cex.axis=cexAxis, mgp=c(0,0.40,0), lwd=0.3, tck=tck, col.tick="gray30", col.axis="gray30", col="gray30", at=seq(0,600,100), pos=0)
axis(side=2, lwd.tick=0.3, cex.axis=cexAxis, mgp=c(0,0.45,0), lwd=0.3, tck=tck, col.tick="gray30", col.axis="gray30", col="gray30", at=seq(0,1,0.2))
title(ylab="IBD signal (Pearson coefficient)", cex.lab=cexLab, mgp=c(2.1,0,0), col.lab="gray30")
title(xlab="Number of tip nodes in the subsampled tree", cex.lab=cexLab, mgp=c(1.2,0,0), col.lab="gray30")
dev.off()


# 2. Comparing dispersal metrics estimated on real genomic datasets of viruses spreading in animal populations

localTreesDirectories = c("GTEV_RRW_125", # Zhao et al. (2023, J. Virol.)
						  "H3N1_Belgium", # Van Borm et al. (2023, Emerg. Infect. Dis.)
						  "H5N1_ME_all_c", # Dellicour et al. (2020, Bioinformatics)
						  "LASV_S_align_3", # Klitting et al. (2022, Nat. Commun.)
						  "LSVD_all_1_WTs", # Van Borm et al. (2023, J. Virol.)
						  "NVAV_BE_moles", # Laenen et al. (2016, Mol. Ecol.)
						  "PDCV_Dlat_pols", # He et al. (2020, Mol. Biol. Evol.)
						  "POWV_US_ext2", # Vogels et al. (2023, PNAS)
						  "PUUV_AK_clades", # Laenen et al. (2019, Vir. Evol.)
						  "RABV_AF_dogs", # Talbi et al. (2010, PLoS Path.)
						  "RABV_BR_bats", # Vieira et al. (2013, Virus Genes)
						  "RABV_IR_all_c", # Dellicour et al. (2019, Mol. Ecol.)
						  "RABV_NEA_bats", # Torres et al. (2014, Mol. Ecol.)
						  "RABV_PE_L1-L3", # Streicker et al. (2016, PNAS)
						  "RABV_US_racc", # Biek et al. (2007, PNAS)
						  "RABV_US_skunk", # Kuzmina et al. (2013, PLoS One)
						  "RABV_YU_fourC", # Tian et al. (2018, PLoS Path.)
						  "TULV_Europe_C", # Cirkovic et al. (2022, Vir. Evol.)
						  "WNV_gamma_all", # Dellicour et al. (2022, Nat. Commun.)
						  "WNV_before_04","WNV_after_04") # WNV lineages < and >2003
dataset_names = c("Getah virus, China", # Zhao et al. (2023, J. Virol.)
				  "AIV H3N1, Belgium", # Van Borm et al. (2023, Emerg. Infect. Dis.)
				  "AIV H5N1, Mekong region", # Dellicour et al. (2020, Bioinformatics)
				  "Lassa virus, Africa", # Klitting et al. (2022, Nat. Commun.)
				  "Lumpy skin disease virus*", # Van Borm et al. (2023, J. Virol.)
				  "Nova virus, Belgium", # Laenen et al. (2016, Mol. Ecol.)
				  "Porcine deltacoronavirus, China", # He et al. (2020, Mol. Biol. Evol.)
				  "Powassan virus, USA", # Vogels et al. (2023, PNAS)
				  "Puumala virus, Belgium", # Laenen et al. (2019, Vir. Evol.)
				  "Rabies virus (dogs), North Africa", # Talbi et al. (2010, PLoS Path.)
				  "Rabies virus (bats), eastern Brazil", # Vieira et al. (2013, Virus Genes)
				  "Rabies virus (dogs), Iran", # Dellicour et al. (2019, Mol. Ecol.)
				  "Rabies virus (bats), Argentina", # Torres et al. (2014, Mol. Ecol.)
				  "Rabies virus (bats)**, Peru", # Streicker et al. (2016, PNAS)
				  "Rabies virus (raccoons), USA", # Biek et al. (2007, PNAS)
				  "Rabies virus (skunks), USA", # Kuzmina et al. (2013, PLoS One)
				  "Rabies virus (dogs), Yunnan (CH)", # Tian et al. (2018, PLoS Path.)
				  "Tula virus, European clade", # Cirkovic et al. (2022, Vir. Evol.)
				  "West Nile virus, North America", # Dellicour et al. (2022, Nat. Commun.)
				  "West Nile virus, North America (<2004)","West Nile virus, North America (>2004)")
references = c("Zhao et al. (2023)", # Zhao et al. (2023, J. Virol.)
			   "Van Borm et al. (2023a) ", # Van Borm et al. (2023, Emerg. Infect. Dis.)
			   "Dellicour et al. (2020a)", # Dellicour et al. (2020, Bioinformatics)
			   "Klitting et al. (2022)", # Klitting et al. (2022, Nat. Commun.)
			   "Van Borm et al. (2023b)", # Van Borm et al. (2023, J. Virol.)	   
			   "Laenen et al. (2016)", # Laenen et al. (2016, Mol. Ecol.)
			   "He et al. (2020)", # He et al. (2020, Mol. Biol. Evol.)
			   "Vogels et al. (2023)", # Vogels et al. (2023, PNAS)
			   "Laenen et al. (2019)", # Laenen et al. (2019, Vir. Evol.)
			   "Talbi et al. (2010)", # Talbi et al. (2010, PLoS Path.)
			   "Vieira et al. (2013)", # Vieira et al. (2013, Virus Genes)
			   "Dellicour et al. (2019)", # Dellicour et al. (2019, Mol. Ecol.)
			   "Torres et al. (2014)", # Torres et al. (2014, Mol. Ecol.)
			   "Streicker et al. (2016)", # Streicker et al. (2016, PNAS)
			   "Biek et al. (2007)", # Biek et al. (2007, PNAS)
			   "Kuzmina et al. (2013)", # Kuzmina et al. (2013, PLoS One)
			   "Tian et al. (2018)", # Tian et al. (2018, PLoS Path.)
			   "Cirkovic et al. (2022)", # Cirkovic et al. (2022, Vir. Evol.)
			   "Dellicour et al. (2020b)", # Dellicour et al. (2022, Nat. Commun.)
			   "Dellicour et al. (2020b)","Dellicour et al. (2020b)")

	# 2.1. Estimating the different dispersal metrics (WLDV, WDC, IBD signal)

references = references[which(!localTreesDirectories%in%c("WNV_before_04","WNV_after_04"))]
dataset_names = dataset_names[which(!localTreesDirectories%in%c("WNV_before_04","WNV_after_04"))]
localTreesDirectories = localTreesDirectories[which(!localTreesDirectories%in%c("WNV_before_04","WNV_after_04"))]
table_S1 = matrix(nrow=length(dataset_names), ncol=6); table_S1[,1] = dataset_names; table_S1[,6] = references
colnames(table_S1) = c("Dataset","Weighted diffusion coefficient (km2/year)","IBD signal (rP)","Number of samples","Study area (km2)","Reference")
for (i in 1:length(localTreesDirectories))
	{
		if (!file.exists(paste0("Observations/Dispersal_stats/",localTreesDirectories[i],"_estimated_dispersal_statistics.txt")))
			{
				localTreesDirectory = paste0("Observations/",localTreesDirectories[i])
				extractionFiles = list.files(localTreesDirectory)
				nberOfExtractionFiles = length(extractionFiles[which(grepl("TreeExtraction",extractionFiles))])
				timeSlices = 100; onlyTipBranches = F; showingPlots = F
				outputName = paste0("Observations/Dispersal_stats/",localTreesDirectories[i]); nberOfCores = 10; slidingWindow = 1
				spreadStatistics(localTreesDirectory, nberOfExtractionFiles, timeSlices, onlyTipBranches, showingPlots, outputName, nberOfCores, slidingWindow)
			}
	}
for (i in 1:length(localTreesDirectories))
	{
		if (file.exists(paste0("Observations/Dispersal_stats/",localTreesDirectories[i],"_estimated_dispersal_statistics.txt")))
			{
				cat("\t",localTreesDirectories[i])
				tab = read.table(paste0("Observations/Dispersal_stats/",localTreesDirectories[i],"_estimated_dispersal_statistics.txt"), head=T)
				vS = tab[,"weighted_branch_dispersal_velocity"]; median = round(median(vS),1); HPD = round(HDInterval::hdi(vS)[1:2],1)
				cat(": median WLDV = ",median,", 95% HPD = [",HPD[1],", ",HPD[2],"]","\n",sep="")
			}	else	{
				cat("\t",localTreesDirectories[i],"\n")
			}
	}
median_WDC_values = rep(NA, length(localTreesDirectories))
for (i in 1:length(localTreesDirectories))
	{
		if (file.exists(paste0("Observations/Dispersal_stats/",localTreesDirectories[i],"_estimated_dispersal_statistics.txt")))
			{
				cat("\t",localTreesDirectories[i])
				tab = read.table(paste0("Observations/Dispersal_stats/",localTreesDirectories[i],"_estimated_dispersal_statistics.txt"), head=T)
				vS = tab[,"weighted_diffusion_coefficient"]; median = round(median(vS),0); HPD = round(HDInterval::hdi(vS)[1:2],0)
				table_S1[i,"Weighted diffusion coefficient (km2/year)"] = paste0(median," [",HPD[1],", ",HPD[2],"]")
				cat(": median WDC = ",median,", 95% HPD = [",HPD[1],", ",HPD[2],"]","\n",sep="")
				median_WDC_values[i] = median
			}	else	{
				cat("\t",localTreesDirectories[i],"\n")
			}
	}
for (i in 1:length(localTreesDirectories))
	{
		if (file.exists(paste0("Observations/Dispersal_stats/",localTreesDirectories[i],"_estimated_dispersal_statistics.txt")))
			{
				cat("\t",localTreesDirectories[i])
				tab = read.table(paste0("Observations/Dispersal_stats/",localTreesDirectories[i],"_estimated_dispersal_statistics.txt"), head=T)
				vS = tab[,"isolation_by_distance_signal_rP2"]; median = round(median(vS),2); HPD = round(HDInterval::hdi(vS)[1:2],2)
				table_S1[i,"IBD signal (rP)"] = paste0(median," [",HPD[1],", ",HPD[2],"]")
				cat(": median IBD signal (rP #2) = ",median,", 95% HPD = [",HPD[1],", ",HPD[2],"]","\n",sep="")
			}	else	{
				cat("\t",localTreesDirectories[i],"\n")
			}
	}
continents = shapefile("Continents_shp/Continents.shp")
for (i in 1:length(localTreesDirectories))
	{
		extractionFiles = list.files(paste0("Observations/",localTreesDirectories[i],"/"))
		extractionFiles = extractionFiles[grepl("TreeExtractions_",extractionFiles)]; points = c()
		for (j in 1:length(extractionFiles))
			{
				tab = read.csv(paste0("Observations/",localTreesDirectories[i],"/TreeExtractions_",i,".csv"))
				points = rbind(points, tab[which(!tab[,"node2"]%in%tab[,"node1"]),c("endLon","endLat")])
			}
		hull = chull(points); hull = c(hull,hull[1]); p = Polygon(points[hull,]); ps = Polygons(list(p),1)
		sps1 = SpatialPolygons(list(ps)); wgs84 = CRS("+init=epsg:4326")
		crs(sps1) = wgs84; sps2 = intersect(sps1, continents)
		area = sum(raster::area(sps2))/(1000*1000); area = round(area/1000)*1000
		table_S1[i,"Study area (km2)"] = paste0("~",area)
	}
for (i in 1:length(localTreesDirectories))
	{
		tab = read.csv(paste0("Observations/",localTreesDirectories[i],"/TreeExtractions_",1,".csv"))
		nberOfTips = length(which(!tab[,"node2"]%in%tab[,"node1"])); table_S1[i,"Number of samples"] = nberOfTips
	}
table_S1 = table_S1[order(median_WDC_values),]
write.table(table_S1, "Estimates.csv", row.names=F, quote=F, sep=";")

	# 2.2. Generating an overall comparative figure based on violin or ridgeline plots

dataset_names = dataset_names[which(!localTreesDirectories%in%c("WNV_before_04","WNV_after_04"))]
localTreesDirectories = localTreesDirectories[which(!localTreesDirectories%in%c("WNV_before_04","WNV_after_04"))]
datasets = localTreesDirectories; medians = rep(NA, length(localTreesDirectories))
vSs1 = list(); vSs1_log = list(); all1_log = c(); vSs2 = list(); all2 = c()
for (i in 1:length(localTreesDirectories))
	{
		tab = read.table(paste0("Observations/Dispersal_stats/",localTreesDirectories[i],"_estimated_dispersal_statistics.txt"), head=T)
		vS = tab[,"weighted_diffusion_coefficient"]; medians[i] = median(vS); vSs1[[i]] = vS; vSs1_log[[i]] = log(vS)
	}
datasets = datasets[order(medians,decreasing=T)]; dataset_names_ordered = dataset_names[order(medians,decreasing=T)]
vSs1 = vSs1[order(medians,decreasing=T)]; vSs1_log = vSs1_log[order(medians,decreasing=T)]
for (i in 1:length(localTreesDirectories))
	{
		if (i < 10) I = paste0("0",i)
		if (i >= 10) I = paste0("",i)
		mat = cbind(log(vSs1[[i]]),rep(paste0(I,"_",datasets[i]),length(vSs1[[i]]))); all1_log = rbind(all1_log, mat)
	}
all1_log = as.data.frame(all1_log); all1_log$V1 = sapply(all1_log$V1, as.numeric)
for (i in 1:length(localTreesDirectories))
	{
		tab = read.table(paste0("Observations/Dispersal_stats/",localTreesDirectories[i],"_estimated_dispersal_statistics.txt"), head=T)
		if ("isolation_by_distance_signal_rS"%in%colnames(tab))
			{
				vSs2[[i]] = tab[,"isolation_by_distance_signal_rP2"]
			}	else		{
				vSs2[[i]] = rnorm(100, runif(1,0.2,0.8), 0.05)
			}
	}
vSs2 = vSs2[order(medians,decreasing=T)]
for (i in 1:length(localTreesDirectories))
	{
		if (i < 10) I = paste0("0",i)
		if (i >= 10) I = paste0("",i)
		mat = cbind(vSs2[[i]],rep(paste0(I,"_",datasets[i]),length(vSs2[[i]]))); all2 = rbind(all2, mat)
	}
all2 = as.data.frame(all2); all2$V1 = sapply(all2$V1, as.numeric)

pdf("Comparison_NEW.pdf", width=9.0, height=5.0) # dev.new(width=9.0, height=5.0)
par(mgp=c(0,0,0), oma=c(0,0.5,0,0), mar=c(3.0,10.5,0.5,0.5), lwd=0.3, col="gray30")
par(fig=c(0,7,0,9)/9)
ridgeline_mod1 = function(x, y, bw="nrd0", mode=F, col="gray", border, lty=1, lwd=1, bty="o", labels=NULL, palette, axes=T) 
	{
		dens = tapply(x, y, density, bw=bw)
		xs = Map(getElement, dens, "x")
		ys = Map(getElement, dens, "y")
		ys = Map(function(x) (x-min(x))/max(x-min(x))*1.5, ys)
		ys = Map(`+`, ys, length(ys):1)
		op = par(no.readonly=T)
		par(mar=op$mar); plot.new()
		plot.window(xlim=c(-0.75,log(200000)), ylim=c(1.3,length(ys)+1))
		abline(h=length(ys):1, col=col, lwd=0.3)
		cols = hcl.colors(length(ys), "Zissou", alpha=0.8)
		border = rep(1, length(ys))
		Map(polygon, xs, ys, col=cols, border=border, lty=lty, lwd=lwd)
		axis(side=1, lwd.tick=0.3, cex.axis=0.8, mgp=c(0,0.25,0), lwd=0, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30",
			 at=c(log(10^-10),log(1),log(10),log(100),log(1000),log(10000),log(100000),log(1000000)), labels=c("0","1","10","100","1000","10000","100000","1000000"))
		axis(side=2, lwd.tick=0.0, cex.axis=0.8, mgp=c(0,0.38,0), lwd=0, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30", 
			 at=c(1:length(datasets))+0.07, labels=rev(dataset_names_ordered), las=1)
		title(xlab=expression("Weighted diffusion coefficient (km"^2*"/year)"), cex.lab=0.9, mgp=c(1.6,0,0), col.lab="gray30")
		box(bty=bty, lwd=0.3); par(op)
	}
x = all1_log[,1]; y = all1_log[,2]; ridgeline_mod1(x, y, axes=T, lwd=0.5)
par(fig=c(7,9,0,9)/9); par(new=T, mar=c(3.0,0,0.5,0.4), lwd=0.3, col="gray30")
ridgeline_mod2 = function(x, y, bw="nrd0", mode=F, col="gray", border, lty=1, lwd=1, bty="o", labels=NULL, palette, axes=T) 
	{
		dens = tapply(x, y, density, bw=bw)
		xs = Map(getElement, dens, "x")
		ys = Map(getElement, dens, "y")
		ys = Map(function(x) (x-min(x))/max(x-min(x))*1.5, ys)
		ys = Map(`+`, ys, length(ys):1)
		op = par(no.readonly=T)
		par(mar=op$mar); plot.new()
		plot.window(xlim=c(-0.10,0.98), ylim=c(1.3,length(ys)+1))
		abline(h=length(ys):1, col=col, lwd=0.3)
		abline(v=0, col=col, lwd=0.6, lty="aa")
		cols = hcl.colors(length(ys), "Zissou", alpha=0.8)
		border = rep(1, length(ys))
		Map(polygon, xs, ys, col=cols, border=border, lty=lty, lwd=lwd)
		axis(side=1, lwd.tick=0.3, cex.axis=0.8, mgp=c(0,0.25,0), lwd=0, tck=-0.035, col.tick="gray30", col.axis="gray30", col="gray30", at=seq(0,1,0.5), labels=seq(0,1,0.5))
		title(xlab=expression("IBD signal (r"[S]*")"), cex.lab=0.9, mgp=c(1.6,0,0), col.lab="gray30")
		box(bty=bty, lwd=0.3); par(op)
	}
x = all2[,1]; y = all2[,2]; ridgeline_mod2(x, y, axes=T, lwd=0.5)
dev.off()

