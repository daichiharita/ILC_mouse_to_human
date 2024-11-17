# SLDSC_Implementation
## このtutorialの目標
- bedファイル(/home/k1_taka/reference/LDSCORE/1000G_Phase3_cell_type_groups)からldscoreを作成
- 作成したldscoreを使用してS-LDSCを実施する
- barplotで図示する

## -log10pのbarplotの作成
##ここからはR

library(ggplot2) 

set.seed(1)
library(fs)
library(data.table)
library(dplyr)
library(tidyverse)

## X: CC毎の比較, y: -log10p, facet: sumstats (noexp)
dir_path="~/LDSC/S_LDSC_tutorial/result/20240902_tutorial/sldsc"
filelist=dir(dir_path)

#modified LDSC (全mark含む)
enrichment_allsumstats_df=data.frame(matrix(rep(NA, 11), ncol=11))[0,]
## このtutorialの目標
- bedファイル(/home/k1_taka/reference/LDSCORE/1000G_Phase3_cell_type_groups)からldscoreを作成
- 作成したldscoreを使用してS-LDSCを実施する
- barplotで図示する
for (file_n in filelist){
    sumstats_list=c("Lupus_langefeld.results", "PASS_Intelligence_SavageJansen2018.results",  "PASS_Height1.results", "PASS_CD_deLange2017.results","PASS_Schizophrenia_Pardinas2018.results",  "PASS_BMI1.results", "RA_ishigaki.results", "PASS_SleepDuration_Dashti2019.results","PASS_Type_2_Diabetes.results")
    for (sumstats in sumstats_list){
        df=read.table(paste0(dir_path, "/", file_n, "/", sumstats), header=1)
        annotation_df=df[1, 2:10]
        
        new_df=data.frame(sumstatsname=sumstats,celltype=file_n, annotation_df)
        enrichment_allsumstats_df=rbind(enrichment_allsumstats_df, new_df)
    }
}

enrichment_allsumstats_df$minuslog10p=-log10(pnorm(-enrichment_allsumstats_df$Coefficient_z.score,0,1))
out_f="~/LDSC/S_LDSC_tutorial/result/20240902_tutorial/plot"
dir.create(out_f)
metadata=read.table("/home/k1_taka/reference/LDSCORE/1000G_Phase3_cell_type_groups/names", header=1)
head(metadata)
colnames(metadata)=c("celltype", "cell_type")
metadata$celltype=as.character(metadata$celltype)
metadata
merge_df=left_join(enrichment_allsumstats_df, metadata, by="celltype")
#metadataのcelltype列は元はint型だったのでcharacter型に変換する
head(merge_df)
merge_df2=merge_df %>% separate(cell_type, c("cell_type2", ".bed"), sep="\\.")

merge_df2$sumstatsname=factor(merge_df2$sumstatsname, levels=sumstats_list)
ggplot(merge_df2, aes(x=minuslog10p,y=cell_type2,fill=cell_type2))+ geom_bar(stat = "identity")+facet_wrap(~sumstatsname)
ggsave(paste0(out_f, "/barplot_minuslog10p_all_sumstats.pdf"), width=10, height=10)
