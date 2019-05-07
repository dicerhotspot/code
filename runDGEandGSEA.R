library("limma");
library("edgeR");

# rnaseq, studies
load("data/ALL-TP.RNAseq_RSEM_log2-20140115.RData");
# domMuts, allMuts, sequencedCases, excludedCases
load("data/TP.RData");
# miRfamilyGeneSets
load("data/mirFamiliesAll.RData");
# Hs.c2
load("data/human_c2_v4.rdata");
reactome <- Hs.c2[ grep("REACTOME", names(Hs.c2)) ];

# this is for reproducibility (today's date)
set.seed(20040403);

rsamples <- colnames(rnaseq);
#inferredSamples <- intersect(rsamples, names(which(m53o < max(m53mut))));
inferredSamples <- intersect(rsamples, names(which(m53o < 0)));
badSamples <- intersect(domMuts, names(which(m53mut > 0)));
#samplesOfInterest <- intersect(c(setdiff(rsamples, c(setdiff(allMuts, domMuts), excludedCases)), inferredSamples), sequencedCases); 
samplesOfInterest <- intersect(rsamples, sequencedCases);
mutSamples <- setdiff(intersect(c(setdiff(as.vector(domMuts), excludedCases), inferredSamples), samplesOfInterest), badSamples);
otherMutSamples <- setdiff( intersect(allMuts, samplesOfInterest), mutSamples );

rna <- rnaseq[, samplesOfInterest];
src <- studies[samplesOfInterest];
rm(rnaseq);

mutDist <- table(src[mutSamples]);
mut2Dist <- table(src[otherMutSamples]);
relevantStudies <- intersect( names(mutDist[ which(mutDist > 0 )]), names(mut2Dist[ which(mut2Dist > 0 )]) );

rrna <- rna[, which(src %in% relevantStudies)];
# Assume all NAs as 0 counts
rrna <- 2^rrna;
rrna[ is.na(rrna) ] <- 0;
rsrc <- as.factor(as.matrix(src[which(src %in% relevantStudies)]));
names(rsrc) <- colnames(rrna);
rstatus <- rep("wt", length(rsrc));
names(rstatus) <- names(rsrc);
rstatus[names(rstatus) %in% mutSamples] <- "mut";
rstatus[names(rstatus) %in% otherMutSamples] <- "other";
rstatus <- as.factor(rstatus);
rm(rna);


runDGE <- function(sumStats, srna, sdesign, sstatus, study, sstudy, miRfamilyGeneSets) {
	cat("running", study, "\n");

	#sdesign <- model.matrix(~0+sstatus);

	y <- DGEList(srna);
	y <- calcNormFactors(y);
	v <- voom(y, sdesign);
	fit <- lmFit(v, sdesign);
	cmatrix <- makeContrasts("sstatusmut-sstatuswt", "sstatusmut-sstatusother", "sstatusother-sstatuswt", levels=sdesign);
	fit <- contrasts.fit(fit, cmatrix);
	fit <- eBayes(fit);

	for( contrast in colnames(cmatrix) ) {
		groupAB <- unlist( strsplit(gsub("sstatus", "", contrast), "-") );
		groupA <- groupAB[1];
		groupB <- groupAB[2]; 
		cntA <- length( which(sstatus == groupA) );
		cntB <- length( which(sstatus == groupB) );
		upDownStats <- decideTests(fit, adjust.method="bonferroni")[, contrast];
	
		noOfUp <- sum(upDownStats > 0);
		noOfDown <- sum(upDownStats < 0);
		noOfSig <- noOfUp + noOfDown;
		
		sumStats <- rbind(sumStats, c(study, sstudy, contrast, cntA, cntB, noOfUp, noOfDown)); 

		stab <- topTable(fit, coef=contrast, adjust.method="bonferroni", sort.by="p", n=nrow(srna));
		write.table(stab, file=sprintf("data/DGE-%s-%s-all.tsv", sstudy, contrast), sep="\t", quote=F);
		sigtab <- stab[1:noOfSig, ];
		write.table(sigtab, file=sprintf("data/DGE-%s-%s-significant.tsv", sstudy, contrast), sep="\t", quote=F);

		gsymbols <- gsub("\\|.*$", "", rownames(v));
		geneids <- gsub("^.*\\|", "", rownames(v));
		geneIdx <- symbols2indices(miRfamilyGeneSets, gsymbols);
		geneIdx2 <- symbols2indices(reactome, geneids);
		allGeneIdx <- c(geneIdx, geneIdx2);
		sromer <- romer(allGeneIdx, v, sdesign, contrast=cmatrix[, contrast], nrot=10000);

		sromertab <- topRomer(sromer, alt="up", n=length(miRfamilyGeneSets)+length(reactome));
		sdownpvals <- sromertab[, "Down"];
		sdownpvals.adj <- p.adjust(sdownpvals, method="fdr");
		sromertab <- cbind(sromertab, sdownpvals.adj);
		colnames(sromertab)[ncol(sromertab)] <- "FDR.Down";	
		suppvals <- sromertab[, "Up"];
		suppvals.adj <- p.adjust(suppvals, method="fdr");
		sromertab <- cbind(sromertab, suppvals.adj);
		colnames(sromertab)[ncol(sromertab)] <- "FDR.Up";	
		write.table(sromertab, file=sprintf("data/mirnaGSEA-%s-%s.tsv", sstudy, contrast, row.names=F), sep="\t", quote=F);

		pdf(sprintf("plots/DGE-%s-%s-GSEA-genes.pdf", sstudy, contrast), useDingbats=F);
		#maxDiff <- max(stab[, "logFC"]);
		#minDiff <- min(stab[, "logFC"]);
		maxDiff <- max( abs(stab[, "logFC"]) );
		minDiff <- -1 * maxDiff;

		plot(stab[, "logFC"], rep(0, nrow(stab)), xlim=c(minDiff, maxDiff), main="All Genes", ylim=c(-1, 1), pch="|");
		abline(v=0, col="gray");
		abline(h=0, col="gray");

		significantSets <- names(which(sdownpvals.adj < .25 | suppvals.adj < .25));
		if( length(significantSets) > 0 ) {
			for(sigset in significantSets) {
				setGenes <- rownames(v)[allGeneIdx[[sigset]]];
				
				plot(stab[setGenes, "logFC"], rep(0, length(setGenes)), xlim=c(minDiff, maxDiff), main=sigset, ylim=c(-1, 1), pch="|");
				abline(v=0, col="gray");
				abline(h=0, col="gray");
			}
		}
		tmp <- dev.off();

		save(srna, sstatus, sdesign, sstudy, y, v, fit, stab, sigtab, gsymbols, geneIdx, sromer, sromertab, significantSets, allGeneIdx, contrast, file=sprintf("data/rnaSeq-%s-%s.RData", sstudy, contrast));
	}

	return(sumStats);
}

sumStats <- NULL;

#for(study in relevantStudies[length(relevantStudies)]) {
for(study in relevantStudies) {
	sidx <- which(rsrc == study);
	srna <- rrna[, sidx];
	sstatus <- as.factor(rstatus[sidx]);
	sstudy <- gsub("\\..*$", "", study);
	sdesign <- model.matrix(~0+sstatus);

	sumStats <- runDGE(sumStats, srna, sdesign, sstatus, study, sstudy, miRfamilyGeneSets);
}

sstatus <- rstatus;
sdesign <- model.matrix(~0+sstatus);
sumStats <- runDGE(sumStats, rrna, sdesign, sstatus, "ALL", "ALL", miRfamilyGeneSets);

colnames(sumStats) <- c("OrgStudy", "Study", "Contrast", "GroupA", "GroupB", "NoOfUpGenes", "NoOfDownGenes");

write.table(sumStats, file="data/DGE_summary.tsv", sep="\t", quote=F);
save.image(file="data/rnaseq.RData");
