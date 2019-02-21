#Command 1 (execute this first)
#30 and 42 are latitudes I want, and 280 to 292 is the longitude range (equiv. to 280-360 and 292-360 degrees)
#awk '{if(NF>7 && $3>=30){c=1} if(NF>7 && $3>=42){c=0} if(c>0){print $0}}' multi_1.partition.at_4m.201002 > partition_ww3_at_201002_ds1
#Command  2
#awk '{if(NF>7 && $4>=280){c=1} if(NF>7 && $4>=292){c=0} if(c>0){print $0}}' partition_ww3_at_201002_ds1 > partition_ww3_at_201002_2_ds2

#awk '{if(NF>7 && $3>=34.25){c=1} if(NF>7 && $3>=34.30){c=0} if(c>0){print $0}}' multi_1.partition.wc_4m.200601 > partition_ww3_wc_200601_ds1

#awk '{if(NF>7 && $4>=239.28){c=1} if(NF>7 && $4>=239.35){c=0} if(c>0){print $0}}' partition_ww3_wc_200601_ds1 > partition_ww3_wc_200601_ds2

#awk '{if(NF>7 && $3>=32.5){c=1} if(NF>7 && $3>=34.5){c=0} if(c>0){print $0}}' multi_1.partition.wc_4m.200606 > partition_ww3_wc_200606_scb1

awk '{if(NF>7 && $4>=238.933){c=1} if(NF>7 && $4>=243){c=0} if(c>0){print $0}}' multi_1.partition.wc_4m.200606 > partition_ww3_wc_200606_scb11
