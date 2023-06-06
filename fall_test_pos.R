library(xcms)
library('magrittr')
library(CAMERA)
library(SummarizedExperiment)
library(stats)
library(doParallel)
library("kableExtra")
library("plotly")
library("dplyr")
library("ggplot2")
library("devtools")

pos_path <- c('C:/Users/shiche/Box/EcoChem/Student Projects/Cheng/Data/Fall Creek/Centroid_positive')
pos_files <- list.files(pos_path, pattern = "*.mzXML", recursive = FALSE, full.names = TRUE)
group_names <- sub(".*w_", "", sub("-1.*", "", pos_files))
group_names <- replace(group_names, which(group_names=="170423"), "170412")
pd <- data.frame(sample_name = pos_files, sample_group = group_names)
raw_data <- readMSData(pos_files, pdata = AnnotatedDataFrame(pd), mode = "onDisk")

qual_col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

group_colors <- paste0(c(col_vector, RColorBrewer::brewer.pal(5, "Set1")), "60")
names(group_colors) <- unique(sub(".*w_", "", sub("-1.*", "", pos_files)))
sample_colors <- group_colors[raw_data$sample_group]

tc <- chromatogram(raw_data, aggregationFun = 'sum')
png(file= "tic_positive.png", width = 1560, height = 960)
plot(tc, col = sample_colors)
dev.off()
tcb <- split(tic(raw_data), f = fromFile(raw_data))
png(file= "tic_performance.png", width = 1260, height = 960)
boxplot(tcb, col = sample_colors,ylab = "intensity", main = "Total ion current")
dev.off()

cwp_pos <- CentWaveParam(peakwidth  = c(4.44, 43.15), #Min/Max Peak Width
                         ppm             = 5, #Maximum tolerated m/z deviation
                         noise           = 10000, #Minimum intensity of true peaks
                         snthresh        = 10, #Signal to noise ratio cutoff
                         prefilter       = c(3, 20000)) #Prefilter for ROI, c(peaks, intensity)
xsn <- findChromPeaks(raw_data, param = cwp_pos)
png(file= "xic_positive.png", width = 1260, height = 960)
boxplot(split((chromPeaks(xsn)[,'rtmax'] - chromPeaks(xsn)[,'rtmin']), f = chromPeaks(xsn)[,'sample']),
        col = sample_colors, ylab = 'Peak Width', xlab = 'File', main = 'Distibution of Peak Width')
dev.off()


mzrt_dt <- read.csv('C:/Users/shiche/Documents/ms_rt1.csv', row.names = 1)
for ( n in rownames(mzrt_dt)){
  mzr_1 <- mzrt_dt[n,2] + c(-0.005, 0.005)
  rtr_1 <- mzrt_dt[n,1] + c(-60, 60)
  chr_raw <- chromatogram(raw_data, mz = mzr_1, rt = rtr_1)
  png(file= paste(n, "_raw_01.png",sep="_"), width = 960, height = 960)
  plot(chr_raw, col=sample_colors)
  dev.off()
}

mpp <- MergeNeighboringPeaksParam(expandRt = 2, ppm = 2)
xdata_pp <- refineChromPeaks(xsn, mpp)
xsnc_1 <- refineChromPeaks(xdata_pp, param = CleanPeaksParam(maxPeakwidth = 43.15))
png(file= "refine_positive.png", width = 1260, height = 960)
boxplot(split((chromPeaks(xsnc_1)[,'rtmax'] - chromPeaks(xsnc_1)[,'rtmin']), f = chromPeaks(xsnc_1)[,'sample']),
        col = sample_colors, ylab = 'Peak Width', xlab = 'File', main = 'Distibution of Peak Width')
dev.off()

pdp_pos <- PeakDensityParam(sampleGroups = xsnc_1$sample_group,
                            bw      = 10, #Bandwidth of the smoothing kernel
                            binSize   = 0.005, #Binsize to slice up m/z dimension into
                            minFraction = 0.5, #Minimum fraction of samples a group must be present in
                            minSamples = 2, #Minimum number of samples a group must be present in 
                            maxFeatures     = 20) #maximum number of features per group

# evaluate the internal standard peaks with the grouping settings
#x_is <- matrix(nrow = 31, ncol = 306)
ms_is <- matrix(c(rep(mzrt_dt[,1],306)),nrow = nrow(mzrt_dt), ncol = 306)
rownames(ms_is) = rownames(mzrt_dt)
for ( n in rownames(mzrt_dt)){
  mzr_1 <- mzrt_dt[n,2] + c(-0.005, 0.005)
  rtr_1 <- mzrt_dt[n,1] + c(-40, 40)
  chr_raw <- chromatogram(xsnc_1, mz = mzr_1, rt = rtr_1)
  xic_peaks <- chromPeaks(chr_raw)
  if (length(which(duplicated(xic_peaks[,'sample'])==TRUE))==0){
    print(paste0('No duplicated sample in ',n))
    ms_is[n,xic_peaks[,'sample']] <- xic_peaks[,'rt']
  }
  else {
    dup_samples <- which(duplicated(xic_peaks[,'sample'])==TRUE)
    ms_is[n,xic_peaks[-dup_samples,'sample']] <- xic_peaks[-dup_samples,'rt']
    print(paste0('Duplicated sample in ',n))
    print(xic_peaks[dup_samples,'sample'])
    for (du_nm in xic_peaks[dup_samples,'sample']){
      du_peaks <- xic_peaks[which(xic_peaks[,'sample']==du_nm),]
      uni_peaks <- du_peaks[which(du_peaks[,'into']==max(du_peaks[,'into'])),]
      ms_is[n,du_nm] <- uni_peaks['rt']
    }
  }
  #png(file= paste(n, "_group.png",sep="_"), width = 960, height = 960)
  #plotChromPeakDensity(chr_raw,param = pdp_neg, col = sample_colors, peakBg = sample_colors[chromPeaks(chr_raw)[, "sample"]],
  #                    peakCol = sample_colors[chromPeaks(chr_raw)[, "sample"]],peakPch = 16,xlab = "retention time")
  #dev.off()
}

# identify non-detects
x_nul <- ms_is - matrix(c(rep(mzrt_dt[,1],306)),nrow = 31, ncol = 306)
colnames(as.data.frame(x_nul))[colSums(x_nul==0)>15]
# non-detects of more than 15 IS
group_names[colSums(x_nul==0)>15]
#non detects in more than 15 samples
#which(rowSums(x_nul==0)>15)

xic_peaks <- chromPeaks(chr_raw)
unq_xic <- xic_peaks[!rownames(xic_peaks) %in% rownames(subset(xic_peaks ,duplicated(sample))), ]
rt_plot <- ggplot(data.frame(chromPeaks(chr_raw)), aes(x = sample, y = rt, color=column))+
  geom_point()
mz_plot <- ggplot(data.frame(chromPeaks(chr_raw)), aes(x = sample, y = mz, color=column))+
  geom_point()

# plot raw rata signals to check their distribution
filterFile(xsn, file =  c(3,4,7,8,11,12,64,65,66,274,275)) %>%
  filterRt(rt = rtr_1) %>%
  filterMz(mz = mzr_1) %>%
  plot(type = "XIC")

xgroup2 <- groupChromPeaks(xsnc_1, param = pdp_pos)
res_neg2 <- quantify(xgroup2, method = "maxint", value = 'index')
quality <- featureSummary(xgroup2)

rofeature_chroms1 <- featureChromatograms(xgroup2, features = 7254:7255)
plot(rofeature_chroms1, col = sample_colors, peakBg = sample_colors[chromPeaks(rofeature_chroms1)[, "sample"]])

## peak group alignment
pgp <- PeakGroupsParam(minFraction = 0.9,
                       extraPeaks = 3,
                       smooth = "loess",
                       span=0.6,
                       family = "gaussian",
                       peakGroupsMatrix = ms_is)
xsn_a <- adjustRtime(xgroup2, param = pgp)
png(file= "adjust_peakgroup.png", width = 2160, height = 2160)
plotAdjustedRtime(xsn_a, col = sample_colors, adjusted = FALSE) 
dev.off()

xgroup3 <- groupChromPeaks(xsn_a, param = pdp_pos)
peaks <- featureDefinitions(xgroup3)
hist(peaks$mzmax-peaks$mzmin, 100, xlab = "Da", main="Histogram of m/z deviation between matched peaks")

peaks_aligned <- featureValues(xgroup3, value = "into", method = "maxint", intensity = "into", missing = 0)
kable(peaks_aligned[1:5,], format = "html") %>%
  kable_styling(font_size = 10)

xgroup3 <- filterMsLevel(xgroup3, msLevel = 1L) #convert to MS level1, required for analyze.xcms.group, CAMERA, and warpgroup
saveRDS(xgroup3, 'Fall_pos_XCMS_result_031923.rds')
devtools::source_url("https://gitlab.com/R_packages/chemhelper/raw/master/R/xcms_helpers.R")
analyze.xcms.group(as(xgroup3, "xcmsSet"), mz = mzrt_dt[2,2], rt = mzrt_dt[2,1], rt_tol_sample=60,
                   mz_tol_sample=0.008, rt_tol_group=60, mz_tol_group=0.008)

xgroup3 <- filterMsLevel(xgroup3, msLevel = 1L)
xgroup3 <- filterRt(xgroup3, rt = c(360, 1260))

## warpgroup
#cl = makeCluster(detectCores() - 1)
#registerDoParallel(cl)
#xr.l = llply(xgroup3@phenoData@data[,1], xcmsRaw, profstep=0)
#xs.warpgroup = group.warpgroup(as(xgroup3, "xcmsSet"), xr.l = xr.l, rt.max.drift = 20, ppm.max.drift = 3, rt.aligned.lim = 5)
## filter function problem
#CAMERA
xsa <- xsAnnotate(as(xgroup3, "xcmsSet"), sample=NA, nSlaves = 4, polarity = 'positive')
xsaF <- groupFWHM(xsa, perfwhm =0.1, intval = "into", sigma = 6) #perfwhm:percentage of the width of the FWHM,Same time-point is defined about the Rt_med +/- FWHM * perfwhm
## FWHM = 2.35 * SD = 2.35 * (rt.max-rt.min)/sigma
xsaFI <- findIsotopes(xsaF, ppm = 5, mzabs= 0.002,  intval = "into")
xsaC <-  groupCorr(xsaFI, calcIso = TRUE, calcCiS = TRUE, calcCaS = TRUE, cor_eic_th=0.85,
                   cor_exp_th=0.85, pval= 0.000001, graphMethod="lpc",intval="into")


#+devtools::install_github("stanstrup/commonMZ")
#library(commonMZ)
#rules <- MZ_CAMERA(mode = "pos", warn_clash = TRUE, clash_ppm = 5)
#rules <- as.data.frame(rules)
rule_file <- system.file('rules/primary_adducts_pos.csv', package = "CAMERA")
rules_prim <- read.csv(rule_file)
xsaFA <- findAdducts(xsaC, ppm=2.5, mzabs=0.001, multiplier=4, polarity="positive", rules=rules_prim)
#xsaFA_def <- findAdducts(xsaFI, ppm=5, mzabs=0.002, multiplier=4, polarity="positive")

rules_count <- 
  xsaFA@annoID[, "ruleID"] %>% 
  factor(levels=1:nrow(rules_prim)) %>% 
  table %>% as.matrix %>% as.data.frame %>% 
  bind_cols(rules_prim[,'name',drop=FALSE],.) %>% 
  setNames(c("Rule","Count")) %>% 
  arrange(desc(Count)) %>% 
  as.data.frame
pie(as.numeric(rules_count[,"Count"]),labels=rules_count[,"Rule"])

saveRDS(xsaFA, 'Fall_pos_CAMERA_result_031923.rds')
fall_camera <- readRDS("C:/Users/shiche/Documents/CAMERA_result.rds")
peaklist <- getPeaklist(xsaFA)
write.csv(peaklist, 'Fall_pos_CAMERA_peaktable.csv')
#peaklist_df <- getPeaklist(xsaFA_def)

reducedpklist <- getReducedPeaklist(xsaFA, method = "sum", cleanup= TRUE)
peaklist %>% select("mz", "rt", "isotopes", "adduct", "pcgroup") %>% 
  kable(format="html", padding=0) %>% 
  kable_styling(font_size = 10)

peaklist %>% filter(abs(mz-508.3334)<0.001) %>% 
  pull(pcgroup) %>% 
  {filter(peaklist, pcgroup %in% .)} %>% 
  select("mz", "rt", "isotopes", "adduct", "pcgroup") %>% 
  arrange(pcgroup, mz) %>% 
  kable(format="html", padding=0) %>% 
  kable_styling(font_size = 15)

#plotPsSpectrum(xsaFA,1,maxlabel=10)
#plotEICs(xsaFA,2)

# find best IS for normalization

#peaklist.norm <- peaklist[,87:392]

final_is <- matrix(nrow = nrow(mzrt_dt), ncol = 306)
rownames(final_is) = rownames(mzrt_dt)
is_detects <- c()
for ( n in c(1:nrow(mzrt_dt))){
  detects <- which(with(peaklist,
                        mz<= mzrt_dt[n,2]+0.005 & mz >= mzrt_dt[n,2]-0.005 & 
                          rt <= mzrt_dt[n,1]+30 & rt >= mzrt_dt[n,1]-30))
  if (length(detects)>1){
    print(mzrt_dt[n,])
    print(length(detects))
    max_detect <- which.max(peaklist[detects,'npeaks'])
    detects <- detects[max_detect]
  }
  is_detects <- c(is_detects, detects)
  final_is[n, ] <- as.numeric(peaklist[detects,87:392])
}

scale_is <- rowMeans(final_is, na.rm = TRUE)/peaklist[is_detects,87:392]
scale_is <- replace(scale_is, is.na(scale_is), 1)
peaklist.norm <- peaklist[,87:392]*colMeans(scale_is)
merged.names <- unique(xgroup3$sample_group)

# merged.peaktable <- matrix(nrow = nrow(reducedpklist), ncol = length(merged.names))
merged.peaktable.full <- matrix(nrow = nrow(peaklist.norm), ncol = length(merged.names))
colnames(merged.peaktable.full) <- merged.names
peaklist.norm <- replace(peaklist.norm, is.na(peaklist.norm), 0)
for (i in c(1:length(merged.names))) {
  grp <- merged.names[i]
  grp_peaks <- peaklist.norm[,group_names==grp]
  if (ncol(grp_peaks)==4) {
    blank <- grp_peaks[,1]
    min_sample <- apply(grp_peaks[,2:4], 1, FUN = min)
    real <- blank < min_sample/20 & min_sample > 10^5
    merged.peaktable.full[real,grp] <- rowMeans(grp_peaks[,2:4][real,])
  }
  else if (ncol(grp_peaks)==3){
    print(paste(grp,"with 3 samples"))
    print(pos_files[group_names==merged.names[i-1]][1])
    print(pos_files[group_names==merged.names[i+1]][1])
    blank <- (peaklist.norm[,group_names==merged.names[i-1]][,1] + 
                peaklist.norm[,group_names==merged.names[i+1]][,1])/2
    min_sample <- apply(grp_peaks[,1:3], 1, FUN = min)
    real <- blank < min_sample/20 & min_sample > 10^5
    merged.peaktable.full[real,grp] <- rowMeans(grp_peaks[,1:3][real,])
  }
  else {
    print(paste(grp,"with less than 3 samples"))
    print(pos_files[group_names==merged.names[i-1]][1])
    print(pos_files[group_names==merged.names[i+1]][1])
    blank <- (peaklist.norm[,group_names==merged.names[i-1]][,1] + 
                peaklist.norm[,group_names==merged.names[i+1]][,1])/2
    min_sample <- apply(grp_peaks[,1:2], 1, FUN = min)
    real <- blank < min_sample/20 & min_sample > 10^5
    merged.peaktable.full[real,grp] <- rowMeans(grp_peaks[,1:2][real,])
  }
}
row.names(merged.peaktable.full) <- row.names(peaklist)
#merged_raw_peaktable <- merged.peaktable
merged_raw_peaktable.full <- merged.peaktable.full
# remove empty rows (absent in all samples)
merged.peaktable.full <- merged.peaktable.full[rowSums(is.na(merged.peaktable.full)) != ncol(merged.peaktable.full), ]
#merged.peaktable.full <- merged.peaktable.full[rowSums(merged.peaktable.full) != 0, ]
merged.features.def <- peaklist[row.names(merged.peaktable.full),c(c(1:7),c(393:395))]
write.csv(merged.peaktable.full, 'fall_pos_merged_peak_list_031923.csv')
write.csv(merged.features.def, 'fall_pos_merged_peak_definition_031923.csv')


chem.count <- colSums(is.na(merged.peaktable.full)==FALSE)
chems <- data.frame(time = merged.names, count = chem.count)
tp <- ggplot(chems, aes(x = time, y = count )) +
  geom_point()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

plotly::ggplotly(tp)

merged.peaktable <- replace(merged.peaktable, is.na(merged.peaktable), 0)
merged.peaktable <- log10(merged.peaktable+1)

# present in at least 4 samples (maximum 74 zeros)
merged.peaktable.4 <- over1.peaktable[which(rowSums(over1.peaktable==0)<75),]

peaks.trans <- (merged.peaktable - apply(merged.peaktable, 1, mean))/apply(merged.peaktable, 1, sd)
ft_df <- mvfft(t(peaks.trans))
#mag <- 1/2*sqrt(Re(ft_df)^2 + Im(ft_df)^2)
#pha <- atan(Im(ft_df)/Re(ft_df))
#mon_f <- seq(ncol(merged.peaktable))
#plot(mon_f, mag[,836],type='b', col='forestgreen', lwd=3)
#plot(mon_f, merged.peaktable[836,],type='b', col='orange', lwd=3)
order(rowSums(merged.peaktable>0))[nrow(merged.peaktable)] # the most abundant chemical
