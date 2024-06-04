library(data.table)

sampColors = fread('~/Phd/Projects/2021-snap_hap/sample_info/colour_info/samples+colour_merged.csv', header=TRUE)
cols2keep = c('PlantID', 'ColorClass_1','Comments')
sampColors = sampColors[,..cols2keep]
str(sampColors)

AmSamples = fread('~/Phd/Projects/2021-snap_hap/sample_info/samples_oldUpdates/samples_Amajus_SnapHap_LastUpdate-2023-10.txt')
allSamples = fread('~/Phd/Projects/2021-snap_hap/sample_info/samples_oldUpdates/samples_ALL_SnapHap_LastUpdate-2023-10.txt')
str(AmSamples)
str(allSamples)

merged_AmSamples = merge(AmSamples[,-'Comments'], sampColors, by='PlantID', all = TRUE)
str(merged_AmSamples)
# fwrite(merged_AmSamples, '~/Phd/Projects/2021-snap_hap/sample_info/samples_oldUpdates/samples_Amajus_SnapHap_LastUpdate-2024-05.txt', 
#        col.names = T, row.names = F, quote = F, sep = '\t')

merged_allSamples = merge(allSamples[,-'Comments'], sampColors, by='PlantID', all=TRUE)
str(merged_allSamples)
# fwrite(merged_allSamples, '~/Phd/Projects/2021-snap_hap/sample_info/samples_oldUpdates/samples_ALL_SnapHap_LastUpdate-2024-05.txt', 
#        col.names = T, row.names = F, quote = F, sep = '\t')

n96 = fread('./Phd/Projects/2021-snap_hap/sample_info/samples_n96_2021/samples_n96_colorInfo.csv', header=T)
colnames(n96)[1] = 'PlantID'
str(n96)

merged_AmSamples_n96 = merge(merged_AmSamples, n96[,c(1,5)], by='PlantID', all.x = TRUE)
str(merged_AmSamples_n96)
fwrite(merged_AmSamples_n96, '~/Phd/Projects/2021-snap_hap/sample_info/samples_Amajus_SnapHapColor_LastUpdate-2024-05.txt', 
       col.names = T, row.names = F, quote = F, sep = '\t')

merged_allSamples_n96 = merge(merged_allSamples, n96[,c(1,5)], by='PlantID', all.x = TRUE)
str(merged_allSamples_n96)
fwrite(merged_allSamples_n96, '~/Phd/Projects/2021-snap_hap/sample_info/samples_ALL_SnapHapColor_LastUpdate-2024-05.txt', 
       col.names = T, row.names = F, quote = F, sep = '\t')
