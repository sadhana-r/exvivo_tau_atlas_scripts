#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

thk_list=read.table(args[1]);
id=args[2];

#Get rid of outlier histo label values ( < 1)
thk_list <- thk_list[thk_list[,2] >= 1,];

#Loop through all the unique label values
for (label in unique(thk_list[,2])){

  thk_select = thk_list[thk_list[,2] == label,1];
  if(ncol(thk_list) == 3){
  	area_select = thk_list[thk_list[,2] == label,3];
  	area_sum = sum(area_select,na.rm = TRUE);
	num_nan = sum(is.na(thk_select));
  	num_thk = sum(thk_select*area_select, na.rm = TRUE);
	if(num_nan/length(thk_select) > 0.8){
                mean_thk = NA;
                median_thk = NA;
        }
        else{
  		mean_thk = num_thk*2/area_sum;
  		median_thk = median(thk_select,na.rm = TRUE)*2;
	}
  	cat(id,label, mean_thk, median_thk,"\n");
  }
  else{
	num_nan = sum(is.na(thk_select));
	num_thk = sum(thk_select, na.rm = TRUE);
	if(num_nan/length(thk_select) > 0.8){
		mean_thk = NA;
                median_thk = NA;	
	}
	else{
        	mean_thk = num_thk*2/length(na.omit(thk_select));
        	median_thk = median(thk_select,na.rm = TRUE)*2;
	}
	cat(id,label,mean_thk,median_thk,"\n");
   }
}
  
