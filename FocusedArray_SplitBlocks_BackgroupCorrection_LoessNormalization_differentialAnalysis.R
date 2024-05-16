#split blocks according to the block_ind_mappint file
#block_ind_mappint is the experiment design of Focused array indicates which block corresponding to which sample
setwd('D:/LC autoantibodies/LC_diagnosis autoantibodies/LUADlowdensity')

mapping = read.csv('./splitgpr/block_ind_mappint_LUAD.csv', row.names=1, check.names=F)
gprdir<-'./rawgpr/IgG'
outdir<-'./splitgpr/IgG'
gpr_files = list.files(gprdir, patter='.gpr',full.names=T)

for (k in 1:length(gpr_files)) {
  gpr_file = gpr_files[k]	
  simple_file_name = unlist(sapply(colnames(mapping), function(x)grep(paste0(x,'-'), gpr_file, value=T))) ###以芯片号对应到相应的原始gpr文件名
  ## find the header line
  head_row = grep('^\\"Block', scan(gpr_file, n=50, what='character', sep='\n'), value=F)
  gpr_data = read.table(gpr_file, header=T, sep='\t', check.names=F, as.is=T, skip=head_row-1)
  gpr_split_data = split(gpr_data, f=gpr_data$Block)
  
  for (j in 1:length(gpr_split_data)){
    blockName = names(gpr_split_data)[j]	
    blockData = gpr_split_data[[j]]
    write.table(blockData, file=file.path(outdir, paste0(mapping[blockName,names(simple_file_name)], '.gpr')), sep='\t', quote=F, row.names=F)
  }
  
}


###Loess Normalization
##IgG的数据
setwd('D:/LC autoantibodies/LC_diagnosis autoantibodies/LUADlowdensity')

inputDir = './splitgpr/IgG'
targetFile = './targetFile.txt'

## Notice: Before running the script, please check your chip version!
## If such protein ID as "JHU00055.B1C2R1" exists in the "ID" column of your gpr file, set "sameID_for_duplicated_protein" to FALSE. Otherwise, set it to TRUE!
sameID_for_duplicated_protein = F
IDsuffix_for_duplicated_protein = 'R'

## Notice: If "removeFlaggedSpot" is set to TRUE, any spots with -100 flags will be removed before analysis
removeFlaggedSpot = TRUE

####################  In general, no need to change the code below #############
mandatory.Columns = c("Block", "Column", "Row", "Name", "ID", "Flags")
density.Columns = c('F635 Median', 'B635 Median')
aggregation.method = 'mean'
returned = file.path(inputDir, 'rawList_ID.RData')

#############################################################################

inputFiles = list.files(inputDir, recursive=T, pattern='gpr$', full.names=T)
targetData = read.table(targetFile, header=T, sep='\t', as.is=T)

## check headers of the target file
if (!all(c("FileName", "SampleID") %in% colnames(targetData))) stop(sprintf("FileName and SampleID must be in the header of [%s]", targetFile))


E.list = list()
Eb.list = list()

for (i in 1:length(inputFiles)){
  
  fileName = inputFiles[i]
  
  message(sprintf('Processing [%s]', fileName))
  
  skipRow = grep('(\\"Block|Block)', scan(fileName, n=40, what='character', sep='\n'), value=F) - 1
  
  if (skipRow > 0) {
    Data = read.table(fileName, header=T, as.is=T, check.names=F, fill=T, skip=skipRow, sep='\t')
  } else {
    Data = read.table(fileName, header=T, as.is=T, check.names=F, fill=T, sep='\t')
  }
  
  
  ## "Control" "buffer"  "empty"   "ND"
  Data_clean = Data[!(Data$ID %in% c('empty', 'buffer', 'Control', 'ND') | Data$Name %in% c('ND')), ]
  Data_use = Data_clean[c(mandatory.Columns, density.Columns)]
  
  ## flags -100
  index_flags = which(Data_use$Flags == -100)
  
  if (removeFlaggedSpot) {
    Data_use[index_flags, density.Columns[1]] = NA
    Data_use[index_flags, density.Columns[2]] = NA   
  }
  
  if (!sameID_for_duplicated_protein & IDsuffix_for_duplicated_protein == 'R') {
    Data_use_backup = Data_use
    Data_use$ID= sub('(JHU.+)R\\d+', '\\1', Data_use$ID)
    Data_use$ID= sub('(Auto-antigen.+)R\\d+', '\\1', Data_use$ID)
  }
  
  Symbol = Data_use[!duplicated(Data_use[,'ID']), mandatory.Columns]
  
  
  if (aggregation.method == 'mean') {
    #E = tapply(Data_use[,density.Columns[1]], INDEX=Data_use[,'Name'], FUN=function(x)mean(x, na.rm=T))
    #Eb = tapply(Data_use[,density.Columns[2]], INDEX=Data_use[,'Name'], FUN=function(x)mean(x, na.rm=T))
    E = tapply(Data_use[,density.Columns[1]], INDEX=Data_use[,'ID'], FUN=function(x)mean(x, na.rm=T))
    Eb = tapply(Data_use[,density.Columns[2]], INDEX=Data_use[,'ID'], FUN=function(x)mean(x, na.rm=T))
    
    E.list[[i]] = E
    Eb.list[[i]] = Eb
  }
}

E.matrix = do.call(cbind, E.list)
colnames(E.matrix) = targetData$SampleID[match(basename(inputFiles), basename(targetData$FileName))]
Eb.matrix = do.call(cbind, Eb.list)
colnames(Eb.matrix) = targetData$SampleID[match(basename(inputFiles), basename(targetData$FileName))]

raw.List = list(E=E.matrix, Eb=Eb.matrix, Symbol=Symbol)
save(raw.List, file=returned)


#####normalization
#1. impute missing values in raw data due to spots with -100 flags.

library(limma)
library(pheatmap)
library(impute)

########################################## Parameter Sets ###################################
## 1. normalize
setwd('D:/LC autoantibodies/LC_diagnosis autoantibodies/LUADlowdensity')
load('./splitgpr/IgG/rawList_ID.RData') 
normalizationDir = 'normalization.IgG'
normalizationFile = 'normalizedData_fastCyclicloessIgG.RData'

## 2. differential analysis
targetFile = './targetFile.txt'
differentialDir = 'differentialAnalysis.IgG'
p_threshold = 0.05
FC_threshold = c(1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.8, 1.9, 2)
ctrLabels = c('BLD','NHC','Control','BLD','BLD','BLD','NHC','NHC','NHC','Control','Control','Control','NHC')
caseLabels = c('EarlyLUAD','EarlyLUAD','EarlyLUAD','EarlyLUAD.Stage0','EarlyLUAD.StageIA','EarlyLUAD.StageIB','EarlyLUAD.Stage0','EarlyLUAD.StageIA','EarlyLUAD.StageIB','EarlyLUAD.Stage0','EarlyLUAD.StageIA','EarlyLUAD.StageIB','BLD')
columnLabels.matched = c('Group1','Group1','Group2','Group3','Group3','Group3','Group3','Group3','Group3','Group4','Group4','Group4','Group1')

############################################################################################
## 0. check parameters
if (length(ctrLabels) != length(caseLabels)) stop('ctrLabels and caseLabels must have the same length!')

## 1. normalize
dir.create(normalizationDir)
normalizedData = file.path(normalizationDir, normalizationFile)

E.bgCorrected = backgroundCorrect.matrix(E=raw.List$E, Eb=raw.List$Eb, method="normexp", normexp.method="saddle", offset=0, printer=NULL, verbose=TRUE)
E.bgCorrected.logScale = log2(E.bgCorrected)
E.normalized = normalizeBetweenArrays(E.bgCorrected.logScale, method='cyclicloess',  cyclic.method="fast")###Loess归一化
save(E.normalized, file=normalizedData)

E.normalized.addPos= cbind(raw.List$Symbol[match(rownames(E.normalized), raw.List$Symbol[,'ID']), ], E.normalized)
write.csv(E.normalized.addPos, sub('.RData', '.csv', normalizedData), row.names=F)

## 2. differential analysis
dir.create(differentialDir)

targetData = read.table(targetFile, header=T, sep='\t', as.is=T)


for (k in seq_along(caseLabels)) {
  
  caseLabel = caseLabels[k]
  ctrLabel = ctrLabels[k]
  columnLabel = columnLabels.matched[k]
  ### ebayes test in limma
  message(sprintf('Performing the differential analysis between [%s] and [%s] via limma', caseLabel, ctrLabel)) 
  
  id.Poor_early = targetData$SampleID[which(targetData[,columnLabel] == ctrLabel)]
  id.Good_early = targetData$SampleID[which(targetData[,columnLabel] == caseLabel)]
  E.use = E.normalized[,c(id.Poor_early, id.Good_early)]
  
  ### impute missing values
  E.use.bk = E.use
  if (any(is.nan(E.use))) E.use = impute.knn(E.use)$data
  
  if (columnLabel == 'Efficacy') {
    group.early = factor(rep(c('Poor', 'Good'), c(length(id.Poor_early), length(id.Good_early))), levels=c('Poor', 'Good'))
    
    design.matrix = model.matrix( ~ 0 + group.early)
    colnames(design.matrix) = c('Poor', 'Good')
    fit_1 = lmFit(E.use, design.matrix)
    contrasts.matrix = makeContrasts(contrasts=c('Good-Poor'), levels=design.matrix)
    fit_2 = contrasts.fit(fit_1, contrasts.matrix)
    fit_2 = eBayes(fit_2) 
    
    result_early = topTable(fit_2, adjust='BH', number=Inf, coef='Good-Poor')
    
    ## HC heatmap
    index_sig = which(result_early$P.Value <= p_threshold)
    
    for (z in seq_along(FC_threshold)) {
      
      thres = FC_threshold[z]
      index_sig_for_heatmap = which(result_early$P.Value <= p_threshold & abs(result_early$logFC) > log2(thres) )
      if (length(index_sig_for_heatmap) < 3) next
      de = rownames(result_early)[index_sig_for_heatmap]
      de.E.use = E.use[de, ]
      annotation_col = data.frame(Response=group.early)
      rownames(annotation_col) = c(id.Poor_early, id.Good_early)
      
      pdf(file.path(differentialDir, sprintf('HC_heatmap_of_differential_autoAb_%svs%s_FC%s.pdf', caseLabel, ctrLabel, as.character(thres))))
      pheatmap(de.E.use, scale='row', show_colnames=T, show_rownames=F, clustering_method='ward.D2', clustering_distance_rows='correlation', clustering_distance_cols='euclidean', annotation_col=annotation_col['Response'])
      dev.off() 
      
    }
    
  } else {
    
    group.early = factor(rep(c(ctrLabel, caseLabel), c(length(id.Poor_early), length(id.Good_early))), levels=c(ctrLabel, caseLabel))
    
    design.matrix = model.matrix( ~ 0 + group.early)
    colnames(design.matrix) = c(ctrLabel, caseLabel)
    fit_1 = lmFit(E.use, design.matrix)
    contrasts.matrix = makeContrasts(contrasts=c(sprintf('%s-%s', caseLabel, ctrLabel)), levels=design.matrix)
    fit_2 = contrasts.fit(fit_1, contrasts.matrix)
    fit_2 = eBayes(fit_2)
    
    result_early = topTable(fit_2, adjust='BH', number=Inf, coef=sprintf('%s-%s', caseLabel, ctrLabel))
    
    ## HC heatmap
    index_sig = which(result_early$P.Value <= p_threshold)
    
    for (z in seq_along(FC_threshold)) {
      
      thres = FC_threshold[z]
      message(sprintf('Drawing the heatmap based on differential sites meeting %s fold change between [%s] and [%s] via limma', thres, caseLabel, ctrLabel))
      
      index_sig_for_heatmap = which(result_early$P.Value <= p_threshold & abs(result_early$logFC) > log2(thres) )
      if (length(index_sig_for_heatmap) < 3) next
      de = rownames(result_early)[index_sig_for_heatmap]
      de.E.use = E.use[de, ]
      annotation_col = data.frame(Category=group.early)
      rownames(annotation_col) = c(id.Poor_early, id.Good_early)
      
      pdf(file.path(differentialDir, sprintf('HC_heatmap_of_differential_autoAb_%svs%s_FC%s.pdf', caseLabel, ctrLabel, as.character(thres))))
      pheatmap(de.E.use, scale='row', show_colnames=T, show_rownames=F, clustering_method='ward.D2', clustering_distance_rows='correlation', clustering_distance_cols='euclidean', annotation_col=annotation_col['Category'])
      dev.off()
      
    }
    
  }
  
  #### merge results
  #write.csv(result_early, file.path(differentialDir, sprintf('differentialTest_%svs%s.csv', caseLabel, ctrLabel)), quote=F)
  #write.csv(result_early[index_sig, ], file.path(differentialDir, sprintf('significantDifferentialTest_%svs%s.csv', caseLabel, ctrLabel)), quote=F)
  featureName = E.normalized.addPos$Name[match(rownames(result_early), E.normalized.addPos$ID)]
  result_early_addName = cbind(featureName, result_early)
  write.csv(result_early_addName, file.path(differentialDir, sprintf('differentialTest_%svs%s.csv', caseLabel, ctrLabel)), quote=F)
  write.csv(result_early_addName[index_sig, ], file.path(differentialDir, sprintf('significantDifferentialTest_%svs%s.csv', caseLabel, ctrLabel)), quote=F)
  
  
}



###saving files###
library(dplyr)
library(pROC)
if (any(is.nan(E.normalized))) E.normknn= impute.knn(E.normalized)$data
returned1 = file.path('./normalization.IgG', 'E.normknn.RData')
save(E.normknn,file = returned1)
load('./normalization.IgG/E.normknn.RData')
normalizedknn.dat<-E.normknn%>%t()%>%as.data.frame()%>%mutate(Sample=row.names(.))%>%select(513,1:512)
write.csv(normalizedknn.dat,'./normalization.IgG/normalization data of IgG after impute.csv',row.names = F)





