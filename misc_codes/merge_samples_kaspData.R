library(data.table)

## Read 2021 KASP data
kasp = fread('/Users/apal/Phd/Projects/2021-snap_hap/Amajus_KASP/KaspGenotyping-2021.csv', header=T)
str(kasp)

### Markers to keep ###

## Ros: ros_assembly_543443; s678-C8923637_543443; 
## El: ros_assembly_715001; s117reverse_156993
## SULF: s91_78256; s91_78256
## CREMOSA: s1187_290152; s1187_290152
## RUBIA: s261_720757 ## NOT FOUND
markers = c('PlantID',
            's678-C8923637_543443', #ROS
            's117reverse_156993', #EL
            's91_78256', #SULF
            's1187_290152' #CREMOSA
            # 's261_720757' #RUBIA
            )
kasp = kasp[, ..markers]
colnames(kasp) = c('PlantID_UPPER', 'Ros__s678-C8923637_543443', 'EL__s117reverse_156993', 'SULF__s91_78256', 'CRE__s1187_290152')

## Read samples
samples = fread('./Phd/Projects/2021-snap_hap/sample_info/samples_ALL_SnapHapColor_LastUpdate-2024-05.txt', header=T)
str(samples)

dat = merge(samples, kasp, by='PlantID_UPPER', all.x = TRUE)
View(dat)

fwrite(dat, './Phd/Projects/2021-snap_hap/sample_info/test.txt', quote = F, sep = '\t', row.names = F, col.names = T)
