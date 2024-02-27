library(adephylo)
library(diagram)
library(HDInterval)
library(lubridate)
library(maptools)
library(rgdal)
library(seraphim)

source("R_functions/simulatorBRW.r")
source("R_functions/simulatorRRW2.r")

# 1. Investigating the impact of sampling effort on the estimate of dispersal metrics (using BRW or RRW simulations)

mostRecentSamplingDatum = 2016.6475
nberOfExtractionFiles = 100
mcc_tab = read.csv("WNV_MCC.csv", head=T)

	# 1.1. Investigating the consistency of dispersal metrics when estimated on Brownian random walk (BRW) simulations

simulationDirectory = "Simulations_1"
dir.create(file.path(paste0(simulationDirectory)), showWarnings=F)
dir.create(file.path(paste0(simulationDirectory,"/All_200_simulations")), showWarnings=F)
dir.create(file.path(paste0(simulationDirectory,"/Subsampling_of_tips")), showWarnings=F)
envVariable = raster("Template_R.tif"); envVariable[] = 1
d = res(envVariable)[1]*2
mcc_tab = read.csv("WNV_MCC.csv", head=T)
ancestID = which(!mcc_tab[,"node1"]%in%mcc_tab[,"node2"])[1]
ancestPosition = c(mcc_tab[ancestID,"startLon"], mcc_tab[ancestID,"startLat"])
birthRate = 0.6
samplingRate = 0.4
startingYear = 1999
samplingWindow = cbind(2001,2017)
timeSlice = 1/365.25
timeIntervale = 1
showingPlots = FALSE
extractionOfValuesOnMatrix = FALSE
for (i in 1:200)
	{
		if (!file.exists(paste0(simulationDirectory,"/All_200_simulations/BRW_simulation_",i,".csv")))
			{
				simulation = NULL; worked = FALSE
				while (worked == FALSE)
					{
						trycatch = tryCatch(
							{
								simulation = simulatorBRW(envVariable, d, ancestPosition, birthRate, samplingRate, startingYear, samplingWindow,
														  timeSlice, timeIntervale, showingPlots, extractionOfValuesOnMatrix)
							},	error = function(cond) {
							},	finally = {
							})
						if (!is.null(simulation)) worked = TRUE
					}
				write.csv(simulation[[1]], paste0(simulationDirectory,"/All_200_simulations/BRW_simulation_",i,".csv"), quote=F, row.names=F)
				write.tree(simulation[[2]], paste0(simulationDirectory,"/All_200_simulations/BRW_simulation_",i,".tre"))
			}
	}

	# 1.2. RRW simulations conducted on a uniform raster to investigate the consistency of dispersal metrics

simulationDirectory = "Simulations_2"; scalingValue = 1.5
dir.create(file.path(paste0(simulationDirectory)), showWarnings=F)
dir.create(file.path(paste0(simulationDirectory,"/All_200_simulations")), showWarnings=F)
dir.create(file.path(paste0(simulationDirectory,"/Subsampling_of_tips")), showWarnings=F)
envVariable = crop(raster("Template_R.tif"), extent(-100,-50,20,60))
envVariable[] = 1; resistance = FALSE
mcc_tab = read.csv("WNV_MCC.csv", head=T)
ancestID = which(!mcc_tab[,"node1"]%in%mcc_tab[,"node2"])[1]
ancestPosition = c(mcc_tab[ancestID,"startLon"], mcc_tab[ancestID,"startLat"])
birthRate = 0.6
samplingRate = 0.4
startingYear = 1999
samplingWindow = cbind(2001,2017)
timeSlice = 1/365.25
timeIntervale = 1
showingPlots = FALSE
extractionOfValuesOnMatrix = FALSE
for (i in 1:200)
	{
		if (!file.exists(paste0(simulationDirectory,"/All_200_simulations/RRW_simulation_",i,".csv")))
			{
				simulation = NULL; worked = FALSE
				while (worked == FALSE)
					{
						trycatch = tryCatch(
							{
								simulation = simulatorRRW2(envVariable, resistance, scalingValue, ancestPosition, birthRate, samplingRate, startingYear, 
														   samplingWindow, timeSlice, timeIntervale, showingPlots, extractionOfValuesOnMatrix)
							},	error = function(cond) {
							},	finally = {
							})
						if (!is.null(simulation)) worked = TRUE
					}								   
				write.csv(simulation[[1]], paste0(simulationDirectory,"/All_200_simulations/RRW_simulation_",i,".csv"), quote=F, row.names=F)
				write.tree(simulation[[2]], paste0(simulationDirectory,"/All_200_simulations/RRW_simulation_",i,".tre"))
			}
	}

	# 1.3. Estimating the dispersal metrics based on the BRW and RRW simulations conducted on rasters

simulationDirectory = "Simulations_1"; prefix = "BRW"
simulationDirectory = "Simulations_2"; prefix = "RRW"

nberOfRemainingSimulations = 50; nberOfBranchesToKeep = seq(500,50,-50)
for (i in 1:200)
	{
		if (i == 1) n = 0
		tab = read.csv(paste0(simulationDirectory,"/All_200_simulations/",prefix,"_simulation_",i,".csv"), head=T)
		nTips = length(which(!tab[,"node2"]%in%tab[,"node1"]))
		if (nTips >= 500)
			{
				n = n+1; write.csv(tab, paste0(simulationDirectory,"/Subsampling_of_tips/",prefix,"_simulation_",n,"_ALL.csv"), quote=F, row.names=F)
				file.copy(paste0(simulationDirectory,"/All_200_simulations/",prefix,"_simulation_",i,".tre"),
						  paste0(simulationDirectory,"/Subsampling_of_tips/",prefix,"_simulation_",n,"_ALL.tre"))
			}
	}
startingYear = 1999; samplingWindow = cbind(2001,2017)
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
				nodelabels(node=tree$edge[j,2], pch=16, cex=0.5, col=cols3[j])
				nodelabels(node=tree$edge[j,2], pch=1, cex=0.5, col="gray30", lwd=0.25)
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
		axis(side=1, lwd.tick=0.3, cex.axis=0.5, mgp=c(0,-0.10,0), lwd=0.3, tck=-0.010, col.tick="gray30", col.axis="gray30", col="gray30", at=seq(-10,10,1))
		axis(side=2, lwd.tick=0.3, cex.axis=0.5, mgp=c(0,0.20,0), lwd=0.3, tck=-0.010, col.tick="gray30", col.axis="gray30", col="gray30", at=seq(-10,10,1))
		dev.off()
	}
for (i in 1:length(nberOfBranchesToKeep))
	{
		for (j in 1:nberOfRemainingSimulations)
			{
				tab1 = read.csv(paste0(simulationDirectory,"/Subsampling_of_tips/",prefix,"_simulation_",j,"_ALL.csv"), head=T)
				colnames(tab1)[which(colnames(tab)=="treeID")] = "length"; tab1[,"length"] = tab1[,"endYear"]-tab1[,"startYear"]
				tipBranches = which(!tab1[,"node2"]%in%tab1[,"node1"]); nberOfTips = length(tipBranches)
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
			}
	}
tab_WLDV = matrix(nrow=nberOfRemainingSimulations, ncol=length(nberOfBranchesToKeep)); colnames(tab_WLDV) = paste0(nberOfBranchesToKeep,"_tips")
tab_WDCs = matrix(nrow=nberOfRemainingSimulations, ncol=length(nberOfBranchesToKeep)); colnames(tab_WDCs) = paste0(nberOfBranchesToKeep,"_tips")
tab_IBDs = matrix(nrow=nberOfRemainingSimulations, ncol=length(nberOfBranchesToKeep)); colnames(tab_IBDs) = paste0(nberOfBranchesToKeep,"_tips")
for (i in 1:length(nberOfBranchesToKeep))
	{
		for (j in 1:nberOfRemainingSimulations)
			{
				print(c(i,j))
				tab = read.csv(paste0(simulationDirectory,"/Subsampling_of_tips/",prefix,"_simulation_",j,"_",nberOfBranchesToKeep[i],".csv"), head=T)
				tre = read.tree(paste0(simulationDirectory,"/Subsampling_of_tips/",prefix,"_simulation_",j,"_ALL.tre"))
				tab[,"greatCircleDist_km"] = diag(rdist.earth(cbind(tab[,"startLon"],tab[,"startLat"]), cbind(tab[,"endLon"],tab[,"endLat"]), miles=F))
				tab_WLDV[j,i] = sum(tab[,"greatCircleDist_km"])/(sum(tab[,"length"]))
				tab_WDCs[j,i] = sum(tab[,"greatCircleDist_km"]^2)/(sum(4*tab[,"length"]))
				tipNodeIndices = which(!tab[,"node2"]%in%tab[,"node1"]); tipNodes = tab[tipNodeIndices,"node2"]
				distsGeo = rdist.earth(cbind(tab[tipNodeIndices,"endLon"],tab[tipNodeIndices,"endLat"]), cbind(tab[tipNodeIndices,"endLon"],tab[tipNodeIndices,"endLat"]), miles=F)
				distsPat = as.matrix(distTips(tre, method="patristic")); distsPat = distsPat[as.character(tipNodes),as.character(tipNodes)]
				distsGeo = distsGeo[lower.tri(distsGeo)]; distsPat = distsPat[lower.tri(distsPat)]
				tab_IBDs[j,i] = cor(distsPat[lower.tri(distsPat)],distsGeo[lower.tri(distsGeo)], method="spearman")
			}
	}
write.csv(tab_WLDV, paste0(simulationDirectory,"/WLDV_estimates.csv"), row.names=F, quote=F)
write.csv(tab_WDCs, paste0(simulationDirectory,"/WDCs_estimates.csv"), row.names=F, quote=F)
write.csv(tab_IBDs, paste0(simulationDirectory,"/IBDs_estimates.csv"), row.names=F, quote=F)
pdf(paste0(simulationDirectory,"/Subsampling_tips_N.pdf"), width=8, height=3.5) # dev.new(width=8, height=3.5)
par(mfrow=c(1,2), mgp=c(0,0,0), oma=c(0,1,0,0), mar=c(2.5,2.2,1.3,1.3), lwd=0.2, col="gray50")
plot(nberOfBranchesToKeep, tab_WLDV[1,], ann=F, axes=F, type="l", lwd=0.1, col="gray50", ylim=c(0,max(tab_WLDV, na.rm=T)))
for (i in 2:dim(tab_WLDV)[1]) lines(nberOfBranchesToKeep, tab_WLDV[i,], lwd=0.2, col="gray30")
lines(nberOfBranchesToKeep, median_WLDV, lwd=1.5, col=rgb(204,0,0,255,maxColorValue=255))
axis(side=1, lwd.tick=0.3, cex.axis=0.6, mgp=c(0,0.00,0), lwd=0.3, tck=-0.02, col.tick="gray30", col.axis="gray30", col="gray30", at=seq(0,600,100), pos=0)
axis(side=2, lwd.tick=0.3, cex.axis=0.6, mgp=c(0,0.25,0), lwd=0.3, tck=-0.02, col.tick="gray30", col.axis="gray30", col="gray30", at=seq(0,30,5))
title(ylab="Weighted lineage dispersal velocity (km/year)", cex.lab=0.7, mgp=c(1.4,0,0), col.lab="gray30")
title(xlab="Number of tip nodes in the subsampled tree", cex.lab=0.7, mgp=c(0.5,0,0), col.lab="gray30")
plot(nberOfBranchesToKeep, tab_WDCs[1,], ann=F, axes=F, type="l", lwd=0.1, col="gray50", ylim=c(0,max(tab_WDCs, na.rm=T)))
for (i in 2:dim(tab_WDCs)[1]) lines(nberOfBranchesToKeep, tab_WDCs[i,], lwd=0.2, col="gray50")
lines(nberOfBranchesToKeep, median_WDCs, lwd=1.5, col=rgb(204,0,0,255,maxColorValue=255))
axis(side=1, lwd.tick=0.3, cex.axis=0.6, mgp=c(0,0.00,0), lwd=0.3, tck=-0.02, col.tick="gray30", col.axis="gray30", col="gray30", at=seq(0,600,100), pos=0)
axis(side=2, lwd.tick=0.3, cex.axis=0.6, mgp=c(0,0.25,0), lwd=0.3, tck=-0.02, col.tick="gray30", col.axis="gray30", col="gray30", at=seq(0,350,50))
title(ylab="Weighted diffusion coefficient (km2/year)", cex.lab=0.7, mgp=c(1.4,0,0), col.lab="gray30")
title(xlab="Number of tip nodes in the subsampled tree", cex.lab=0.7, mgp=c(0.5,0,0), col.lab="gray30")
dev.off()

tab_WLDV = read.csv(paste0(simulationDirectory,"/WLDV_estimates.csv"), head=T)
tab_WDCs = read.csv(paste0(simulationDirectory,"/WDCs_estimates.csv"), head=T)
tab_IBDs = read.csv(paste0(simulationDirectory,"/IBDs_estimates.csv"), head=T)
median_WLDV = rep(NA, length(nberOfBranchesToKeep)); median_WDCs = rep(NA, length(nberOfBranchesToKeep)); median_IBDs = rep(NA, length(nberOfBranchesToKeep))
for (i in 1:length(nberOfBranchesToKeep))
	{
		median_WLDV[i] = median(tab_WLDV[,i], na.rm=T); median_WDCs[i] = median(tab_WDCs[,i], na.rm=T); median_IBDs[i] = median(tab_IBDs[,i], na.rm=T)
	}
pdf(paste0(simulationDirectory,"/Subsampling_tips.pdf"), width=9.0, height=2.7) # dev.new(width=9.0, height=2.7)
par(mfrow=c(1,3), mgp=c(0,0,0), oma=c(0,0.5,0,0), mar=c(2.7,3.0,1.3,1.3), lwd=0.2, col="gray50")
plot(nberOfBranchesToKeep, tab_WLDV[1,], ann=F, axes=F, type="l", lwd=0.1, col="gray50", xlim=c(50,495), ylim=c(0,max(tab_WLDV, na.rm=T)))
for (i in 2:dim(tab_WLDV)[1]) lines(nberOfBranchesToKeep, tab_WLDV[i,], lwd=0.2, col="gray30")
lines(nberOfBranchesToKeep, median_WLDV, lwd=1.5, col=rgb(204,0,0,255,maxColorValue=255))
axis(side=1, lwd.tick=0.3, cex.axis=0.7, mgp=c(0,0.10,0), lwd=0.3, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30", at=seq(0,600,100), pos=0)
axis(side=2, lwd.tick=0.3, cex.axis=0.7, mgp=c(0,0.28,0), lwd=0.3, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30", at=seq(0,30,5))
title(ylab="Weighted lineage dispersal velocity (km/year)", cex.lab=0.9, mgp=c(1.6,0,0), col.lab="gray30")
title(xlab="Number of tip nodes in the subsampled tree", cex.lab=0.9, mgp=c(0.7,0,0), col.lab="gray30")
plot(nberOfBranchesToKeep, tab_WDCs[1,], ann=F, axes=F, type="l", lwd=0.1, col="gray50", xlim=c(50,495), ylim=c(0,max(tab_WDCs, na.rm=T)))
for (i in 2:dim(tab_WDCs)[1]) lines(nberOfBranchesToKeep, tab_WDCs[i,], lwd=0.2, col="gray50")
lines(nberOfBranchesToKeep, median_WDCs, lwd=1.5, col=rgb(204,0,0,255,maxColorValue=255))
axis(side=1, lwd.tick=0.3, cex.axis=0.7, mgp=c(0,0.10,0), lwd=0.3, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30", at=seq(0,600,100), pos=0)
axis(side=2, lwd.tick=0.3, cex.axis=0.7, mgp=c(0,0.28,0), lwd=0.3, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30", at=seq(0,350,50))
title(ylab="Weighted diffusion coefficient (km2/year)", cex.lab=0.9, mgp=c(1.6,0,0), col.lab="gray30")
title(xlab="Number of tip nodes in the subsampled tree", cex.lab=0.9, mgp=c(0.7,0,0), col.lab="gray30")
plot(nberOfBranchesToKeep, tab_IBDs[1,], ann=F, axes=F, type="l", lwd=0.1, col="gray50", xlim=c(50,495), ylim=c(0,0.96))
for (i in 2:dim(tab_IBDs)[1]) lines(nberOfBranchesToKeep, tab_IBDs[i,], lwd=0.2, col="gray50")
lines(nberOfBranchesToKeep, median_IBDs, lwd=1.5, col=rgb(204,0,0,255,maxColorValue=255))
axis(side=1, lwd.tick=0.3, cex.axis=0.7, mgp=c(0,0.10,0), lwd=0.3, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30", at=seq(0,600,100), pos=0)
axis(side=2, lwd.tick=0.3, cex.axis=0.7, mgp=c(0,0.28,0), lwd=0.3, tck=-0.015, col.tick="gray30", col.axis="gray30", col="gray30", at=seq(0,1,0.1))
title(ylab="IBD rS (patristic vs geographic distance)", cex.lab=0.9, mgp=c(1.6,0,0), col.lab="gray30")
title(xlab="Number of tip nodes in the subsampled tree", cex.lab=0.9, mgp=c(0.7,0,0), col.lab="gray30")
dev.off()

# 2. Comparing dispersal metrics estimated on real genomic datasets of viruses spreading in animal populations

localTreeDirectories = c("GTEV_RRW_125", # Zhao et al. (2023, J. Virol.)
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
						 "RABV_US_racc", # Biek et al. (2017, PNAS)
						 "RABV_US_skunk", # Kuzmina et al. (2013, PLoS One)
						 "RABV_YU_fourC", # Tian et al. (2018, PLoS Path.)
						 "TULV_central_E", # Cirkovic et al. (2022, Vir. Evol.)
						 "WNV_gamma_all") # Dellicour et al. (2022, Nat. Commun.)

for (i in 1:length(localTreeDirectories))
	{
		localTreesDirectory = paste0("Observations/",localTreeDirectories[i])
		extractionFiles = list.files(localTreeDirectories[i])
		nberOfExtractionFiles = length(extractionFiles[which(grepl("TreeExtraction",extractionFiles))])
		timeSlices = 100; onlyTipBranches = F; showingPlots = F
		outputName = paste0("Dispersal_stats/",localTreesDirectory); nberOfCores = 5; slidingWindow = 1
		spreadStatistics(localTreesDirectory, nberOfExtractionFiles, timeSlices, onlyTipBranches, showingPlots, outputName, nberOfCores, slidingWindow)
	}
for (i in 1:length(localTreeDirectories))
	{
		tab = read.csv(paste0("Observations/",localTreeDirectories[i],"/TreeExtractions_",i,".csv"))
		tips = length(which(!tab[,"node2"]%in%tab[,"node1"]))
		cat("\t",localTreeDirectories[i],": number of samples = ",tips,"\n",sep="")
	}
for (i in 1:length(localTreeDirectories))
	{
		if (file.exists(paste0("Dispersal_stats/",localTreeDirectories[i],"_estimated_dispersal_statistics.txt")))
			{
				cat("\t",localTreeDirectories[i])
				tab = read.table(paste0("Dispersal_stats/",localTreeDirectories[i],"_estimated_dispersal_statistics.txt"), head=T)
				vS = tab[,"weighted_branch_dispersal_velocity"]; median = round(median(vS),1); HPD = round(HDInterval::hdi(vS)[1:2],1)
				cat(": median WLDV = ",median,", 95% HPD = [",HPD[1],", ",HPD[2],"]","\n",sep="")
			}	else	{
				cat("\t",localTreeDirectories[i],"\n")
			}
	}
for (i in 1:length(localTreeDirectories))
	{
		if (file.exists(paste0("Dispersal_stats/",localTreeDirectories[i],"_estimated_dispersal_statistics.txt")))
			{
				cat("\t",localTreeDirectories[i])
				tab = read.table(paste0("Dispersal_stats/",localTreeDirectories[i],"_estimated_dispersal_statistics.txt"), head=T)
				vS = tab[,"isolation_by_distance_signal_rS"]; median = round(median(vS),1); HPD = round(HDInterval::hdi(vS)[1:2],1)
				cat(": median IBD signal (rS) = ",median,", 95% HPD = [",HPD[1],", ",HPD[2],"]","\n",sep="")
			}	else	{
				cat("\t",localTreeDirectories[i],"\n")
			}
	}

