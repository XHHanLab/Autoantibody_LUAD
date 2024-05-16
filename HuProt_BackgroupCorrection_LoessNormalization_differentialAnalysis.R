#####Loess Normalization
#preprocess

#### 1.check the header of target file. "FileName" and "SampleID" must be in the target file.
#### target file is the metadata of your samples include the path of gprfile (derived from GenePix Pro) and group data
#### 2.remove the restrictions on the option "inputDir" by adding the "basename" function
#### 3.set splots with -100 flags to zero
#### 4. disentangle 12 proteins in more than 1 blocks
#### 5. disentangle additional 1 proteins in more than 1 blocks, because a new chip version from shasha
#### 6. deal with proteins starting with "Auto-antigen"
#### 7. remove name=ND
###############################   options  ######################################
setwd('K:/LC_diagnosis autoantibodies/112samples/Loess/originalData/highdensity')
inputDir = './gpr/IgM'
targetFile = './targetFile.IgM.txt'

## Notice: Before running the script, please check your chip version!
## If such protein ID as "JHU00055.B1C2R1" exists in the "ID" column of your gpr file, set "sameID_for_duplicated_protein" to FALSE. Otherwise, set it to TRUE!
sameID_for_duplicated_protein = F
IDsuffix_for_duplicated_protein = 'R'

## Notice: If "removeFlaggedSpot" is set to TRUE, any spots with -100 flags will be removed before analysis
removeFlaggedSpot = TRUE

####################  In general, no need to change the code below #############
mandatory.Columns = c("Block", "Column", "Row", "Name", "ID", "Flags")
density.Columns = c('F532 Median', 'B532 Median')
aggregation.method = 'mean'
returned = file.path(inputDir, 'rawList_ID.IgM.RData')

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
  
  ## disentangle 13 proteins and auto proteins
  Data_use$ID[grep('JHU16058.B5C17R19|JHU16058.B5C17R20', Data_use_backup$ID, perl=T)] = paste0(Data_use$ID[grep('JHU16058.B5C17R19|JHU16058.B5C17R20', Data_use_backup$ID, perl=T)], '_R1920')
  Data_use$ID[grep('JHU16058.B5C17R21|JHU16058.B5C17R22', Data_use_backup$ID, perl=T)] = paste0(Data_use$ID[grep('JHU16058.B5C17R21|JHU16058.B5C17R22', Data_use_backup$ID, perl=T)], '_R2122')
  
  Data_use$ID[grep('JHU17260.B16C10R31|JHU17260.B16C10R32', Data_use_backup$ID, perl=T)] = paste0(Data_use$ID[grep('JHU17260.B16C10R31|JHU17260.B16C10R32', Data_use_backup$ID, perl=T)], '_R3132')
  Data_use$ID[grep('JHU17260.B16C10R35|JHU17260.B16C10R36', Data_use_backup$ID, perl=T)] = paste0(Data_use$ID[grep('JHU17260.B16C10R35|JHU17260.B16C10R36', Data_use_backup$ID, perl=T)], '_R3536')
  
  Data_use$ID[grep('JHU17633.B16C22R55|JHU17633.B16C22R56', Data_use_backup$ID, perl=T)] = paste0(Data_use$ID[grep('JHU17633.B16C22R55|JHU17633.B16C22R56', Data_use_backup$ID, perl=T)], '_R5556')
  Data_use$ID[grep('JHU17633.B16C22R57|JHU17633.B16C22R58', Data_use_backup$ID, perl=T)] = paste0(Data_use$ID[grep('JHU17633.B16C22R57|JHU17633.B16C22R58', Data_use_backup$ID, perl=T)], '_R5758')
  
  Data_use$ID[grep('JHU16963.B13C31R25|JHU16963.B13C31R26', Data_use_backup$ID, perl=T)] = paste0(Data_use$ID[grep('JHU16963.B13C31R25|JHU16963.B13C31R26', Data_use_backup$ID, perl=T)], '_R2526')
  Data_use$ID[grep('JHU16963.B13C31R27|JHU16963.B13C31R28', Data_use_backup$ID, perl=T)] = paste0(Data_use$ID[grep('JHU16963.B13C31R27|JHU16963.B13C31R28', Data_use_backup$ID, perl=T)], '_R2728')
  
  Data_use$ID[grep('JHU17309.B16C5R37|JHU17309.B16C5R38', Data_use_backup$ID, perl=T)] = paste0(Data_use$ID[grep('JHU17309.B16C5R37|JHU17309.B16C5R38', Data_use_backup$ID, perl=T)], '_R3738')
  Data_use$ID[grep('JHU17309.B16C5R39|JHU17309.B16C5R40', Data_use_backup$ID, perl=T)] = paste0(Data_use$ID[grep('JHU17309.B16C5R39|JHU17309.B16C5R40', Data_use_backup$ID, perl=T)], '_R3940')
  
  Data_use$ID[grep('JHU19752.B14C19R1|JHU19752.B14C19R2', Data_use_backup$ID, perl=T)] = paste0(Data_use$ID[grep('JHU19752.B14C19R1|JHU19752.B14C19R2', Data_use_backup$ID, perl=T)], '_R12')
  Data_use$ID[grep('JHU19752.B14C19R3|JHU19752.B14C19R4', Data_use_backup$ID, perl=T)] = paste0(Data_use$ID[grep('JHU19752.B14C19R3|JHU19752.B14C19R4', Data_use_backup$ID, perl=T)], '_R34')
  
  Data_use$ID[grep('JHU17272.B16C21R33|JHU17272.B16C21R34', Data_use_backup$ID, perl=T)] = paste0(Data_use$ID[grep('JHU17272.B16C21R33|JHU17272.B16C21R34', Data_use_backup$ID, perl=T)], '_R3334')
  Data_use$ID[grep('JHU17272.B16C21R35|JHU17272.B16C21R36', Data_use_backup$ID, perl=T)] = paste0(Data_use$ID[grep('JHU17272.B16C21R35|JHU17272.B16C21R36', Data_use_backup$ID, perl=T)], '_R3536')
  
  Data_use$ID[grep('JHU14199.B9C2R47|JHU14199.B9C2R48', Data_use_backup$ID, perl=T)] = paste0(Data_use$ID[grep('JHU14199.B9C2R47|JHU14199.B9C2R48', Data_use_backup$ID, perl=T)], '_R4748')
  Data_use$ID[grep('JHU14199.B9C2R53|JHU14199.B9C2R54', Data_use_backup$ID, perl=T)] = paste0(Data_use$ID[grep('JHU14199.B9C2R53|JHU14199.B9C2R54', Data_use_backup$ID, perl=T)], '_R5354')
  
  Data_use$ID[grep('JHU14238.B9C3R43|JHU14238.B9C3R44', Data_use_backup$ID, perl=T)] = paste0(Data_use$ID[grep('JHU14238.B9C3R43|JHU14238.B9C3R44', Data_use_backup$ID, perl=T)], '_R4344')
  Data_use$ID[grep('JHU14238.B9C3R49|JHU14238.B9C3R50', Data_use_backup$ID, perl=T)] = paste0(Data_use$ID[grep('JHU14238.B9C3R49|JHU14238.B9C3R50', Data_use_backup$ID, perl=T)], '_R4950')
  
  Data_use$ID[grep('JHU30455.B17C15R21|JHU30455.B17C15R22', Data_use_backup$ID, perl=T)] = paste0(Data_use$ID[grep('JHU30455.B17C15R21|JHU30455.B17C15R22', Data_use_backup$ID, perl=T)], '_R2122')
  Data_use$ID[grep('JHU30455.B17C15R39|JHU30455.B17C15R40', Data_use_backup$ID, perl=T)] = paste0(Data_use$ID[grep('JHU30455.B17C15R39|JHU30455.B17C15R40', Data_use_backup$ID, perl=T)], '_R3940')
  
  Data_use$ID[grep('JHU18327.B13C22R47|JHU18327.B13C22R48', Data_use_backup$ID, perl=T)] = paste0(Data_use$ID[grep('JHU18327.B13C22R47|JHU18327.B13C22R48', Data_use_backup$ID, perl=T)], '_R4748')
  Data_use$ID[grep('JHU18327.B13C22R53|JHU18327.B13C22R54', Data_use_backup$ID, perl=T)] = paste0(Data_use$ID[grep('JHU18327.B13C22R53|JHU18327.B13C22R54', Data_use_backup$ID, perl=T)], '_R5354')
  
  Data_use$ID[grep('JHU14257.B11C15R47|JHU14257.B11C15R48', Data_use_backup$ID, perl=T)] = paste0(Data_use$ID[grep('JHU14257.B11C15R47|JHU14257.B11C15R48', Data_use_backup$ID, perl=T)], '_R4748')
  Data_use$ID[grep('JHU14257.B11C15R53|JHU14257.B11C15R54', Data_use_backup$ID, perl=T)] = paste0(Data_use$ID[grep('JHU14257.B11C15R53|JHU14257.B11C15R54', Data_use_backup$ID, perl=T)], '_R5354')
  
  
  Data_use$ID[grep('JHU15954.B9C6R67|JHU15954.B9C6R68', Data_use_backup$ID, perl=T)] = paste0(Data_use$ID[grep('JHU15954.B9C6R67|JHU15954.B9C6R68', Data_use_backup$ID, perl=T)], '_R6768')
  Data_use$ID[grep('JHU15954.B9C6R87|JHU15954.B9C6R88', Data_use_backup$ID, perl=T)] = paste0(Data_use$ID[grep('JHU15954.B9C6R87|JHU15954.B9C6R88', Data_use_backup$ID, perl=T)], '_R8788')
  
  Data_use$ID[grep('Auto-antigen.B20C16R37|Auto-antigen.B20C16R38', Data_use_backup$ID, perl=T)] = paste0(Data_use$ID[grep('Auto-antigen.B20C16R37|Auto-antigen.B20C16R38', Data_use_backup$ID, perl=T)], '_R3738')
  Data_use$ID[grep('Auto-antigen.B20C16R39|Auto-antigen.B20C16R40', Data_use_backup$ID, perl=T)] = paste0(Data_use$ID[grep('Auto-antigen.B20C16R39|Auto-antigen.B20C16R40', Data_use_backup$ID, perl=T)], '_R3940')
  Data_use$ID[grep('Auto-antigen.B20C16R41|Auto-antigen.B20C16R42', Data_use_backup$ID, perl=T)] = paste0(Data_use$ID[grep('Auto-antigen.B20C16R41|Auto-antigen.B20C16R42', Data_use_backup$ID, perl=T)], '_R4142')
  
  Data_use$ID[grep('Auto-antigen.B20C24R37|Auto-antigen.B20C24R38', Data_use_backup$ID, perl=T)] = paste0(Data_use$ID[grep('Auto-antigen.B20C24R37|Auto-antigen.B20C24R38', Data_use_backup$ID, perl=T)], '_R3738')
  Data_use$ID[grep('Auto-antigen.B20C24R39|Auto-antigen.B20C24R40', Data_use_backup$ID, perl=T)] = paste0(Data_use$ID[grep('Auto-antigen.B20C24R39|Auto-antigen.B20C24R40', Data_use_backup$ID, perl=T)], '_R3940')
  
  Data_use$ID[grep('Auto-antigen.B20C4R37|Auto-antigen.B20C4R38', Data_use_backup$ID, perl=T)] = paste0(Data_use$ID[grep('Auto-antigen.B20C4R37|Auto-antigen.B20C4R38', Data_use_backup$ID, perl=T)], '_R3738')
  Data_use$ID[grep('Auto-antigen.B20C4R39|Auto-antigen.B20C4R40', Data_use_backup$ID, perl=T)] = paste0(Data_use$ID[grep('Auto-antigen.B20C4R39|Auto-antigen.B20C4R40', Data_use_backup$ID, perl=T)], '_R3940')
  Data_use$ID[grep('Auto-antigen.B20C4R41|Auto-antigen.B20C4R42', Data_use_backup$ID, perl=T)] = paste0(Data_use$ID[grep('Auto-antigen.B20C4R41|Auto-antigen.B20C4R42', Data_use_backup$ID, perl=T)], '_R4142')
  
  Data_use$ID[grep('Auto-antigen.B20C6R37|Auto-antigen.B20C6R38', Data_use_backup$ID, perl=T)] = paste0(Data_use$ID[grep('Auto-antigen.B20C6R37|Auto-antigen.B20C6R38', Data_use_backup$ID, perl=T)], '_R3738')
  Data_use$ID[grep('Auto-antigen.B20C6R39|Auto-antigen.B20C6R40', Data_use_backup$ID, perl=T)] = paste0(Data_use$ID[grep('Auto-antigen.B20C6R39|Auto-antigen.B20C6R40', Data_use_backup$ID, perl=T)], '_R3940')
  Data_use$ID[grep('Auto-antigen.B20C6R41|Auto-antigen.B20C6R42', Data_use_backup$ID, perl=T)] = paste0(Data_use$ID[grep('Auto-antigen.B20C6R41|Auto-antigen.B20C6R42', Data_use_backup$ID, perl=T)], '_R4142')
  
  
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
setwd('K:/LC_diagnosis autoantibodies/112samples/Loess/originalData/highdensity')
load('./gpr/IgM/rawList_ID.IgM.RData') ###һ??ע???????ļ??????ļ???
normalizationDir = 'normalization.IgM'
normalizationFile = 'normalizedData_fastCyclicloess.IgM.RData'

## 2. differential analysis
targetFile = './targetFile.IgM.txt'
differentialDir = 'differentialAnalysis.IgM'
p_threshold = 0.05
FC_threshold = c(1.5, 1.6, 1.7, 1.8, 1.9, 2)
ctrLabels = c('Benign','Normal','Control','Benign','Benign','Benign','Normal','Normal','Normal','Control','Control','Control','Normal')
caseLabels = c('Early_LC','Early_LC','Early_LC','Early_LC.0','Early_LC.I','Early_LC.II','Early_LC.0','Early_LC.I','Early_LC.II','Early_LC.0','Early_LC.I','Early_LC.II','Benign')
columnLabels.matched = c('Group1','Group1','Group2','Group3','Group3','Group3','Group3','Group3','Group3','Group4','Group4','Group4','Group1')

############################################################################################
## 0. check parameters
if (length(ctrLabels) != length(caseLabels)) stop('ctrLabels and caseLabels must have the same length!')

## 1. normalize
dir.create(normalizationDir)
normalizedData = file.path(normalizationDir, normalizationFile)

E.bgCorrected = backgroundCorrect.matrix(E=raw.List$E, Eb=raw.List$Eb, method="normexp", normexp.method="saddle", offset=0, printer=NULL, verbose=TRUE)
E.bgCorrected.logScale = log2(E.bgCorrected)
E.normalized = normalizeBetweenArrays(E.bgCorrected.logScale, method='cyclicloess',  cyclic.method="fast")
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
      dev.off() #????ͼƬ?ļ?
      
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
      
      pdf(file.path(differentialDir, sprintf('HC_heatmap_of_differential_autoAb_%svs%s_FC%s.IgM.pdf', caseLabel, ctrLabel, as.character(thres))))
      pheatmap(de.E.use, scale='row', show_colnames=T, show_rownames=F, clustering_method='ward.D2', clustering_distance_rows='correlation', clustering_distance_cols='euclidean', annotation_col=annotation_col['Category'])
      dev.off()
      
    }
    
  }
  
  #### merge results
  #write.csv(result_early, file.path(differentialDir, sprintf('differentialTest_%svs%s.csv', caseLabel, ctrLabel)), quote=F)
  #write.csv(result_early[index_sig, ], file.path(differentialDir, sprintf('significantDifferentialTest_%svs%s.csv', caseLabel, ctrLabel)), quote=F)
  featureName = E.normalized.addPos$Name[match(rownames(result_early), E.normalized.addPos$ID)]
  result_early_addName = cbind(featureName, result_early)
  write.csv(result_early_addName, file.path(differentialDir, sprintf('differentialTest_%svs%s.IgM.csv', caseLabel, ctrLabel)), quote=F)
  write.csv(result_early_addName[index_sig, ], file.path(differentialDir, sprintf('significantDifferentialTest_%svs%s.IgM.csv', caseLabel, ctrLabel)), quote=F)
  
  
}

