library("data.table")
library(ggplot2)
library(dplyr)
library(ggridges)
library(epitools)
library(forcats)


infile="/Volumes/groups/smontgom/shared/UDN/PreprocessingHG38/Splicing/LeafCutterMD/Constraint/combined_splice_constraint.tsv.gz"

data=fread(infile)
print(nrow(data))
colnames(data)<-c("ensg","sample","z","enst","genic_score",
                  "udn_id","hgnc","maf","cadd","regional_score","snp_or_indel","status","tissue")
data<-data %>% dplyr::filter(ensg!=".")
print(nrow(data))
constraint_genic<-apply(data,1,function(x){
  this_genic<-as.numeric(x[which(colnames(data)=="genic_score")])
  if(this_genic>=0.9 & !is.na(this_genic)){"ge0.9"}
  else if(this_genic>=0.75 & this_genic<0.9 & !is.na(this_genic)){"0.75_to_0.9"}
  else if(this_genic>0.25 & this_genic<0.75  & !is.na(this_genic)){"0.25_to_0.75"}
  else if(this_genic<=0.25 & this_genic>0.1  & !is.na(this_genic)){"0.1_to_0.25"}
  else if(this_genic<=0.1 & !is.na(this_genic)){"le0.1"}
  else{"NA"}
})
contraint_regional<-apply(data,1,function(x){
  this_regional<-as.numeric(x[which(colnames(data)=="regional_score")])
  # if(this_regional>0.8 & !is.na(this_regional)){"high"}else if(this_regional<0.2 & !is.na(this_regional)){"low"}else{"NA"}
  if(this_regional>=0.9 & !is.na(this_regional)){"ge0.9"}
  else if(this_regional>=0.75 & this_regional<0.9 & !is.na(this_regional)){"0.75_to_0.9"}
  else if(this_regional>0.25 & this_regional<0.75  & !is.na(this_regional)){"0.25_to_0.75"}
  else if(this_regional<=0.25 & this_regional>0.1  & !is.na(this_regional)){"0.1_to_0.25"}
  else if(this_regional<=0.1 & !is.na(this_regional)){"le0.1"}
  else{"NA"}
  
})
data$binned_genic<- constraint_genic
data$binned_regional<- contraint_regional

data<-data %>% group_by(ensg) %>% mutate(maxZ=max(z))
data$maf<-as.numeric(data$maf)

zmin=2
data_AtLeastOneOutlier<-data%>% dplyr::filter(maxZ>zmin)
data_outliers<-data%>% filter(z>zmin) #data_AtLeastOneOutlier %>% filter(z>zmin) # & tissue=="Fibroblast")
data_nooutliers<-data%>% filter(z<=zmin) #data_AtLeastOneOutlier %>% filter(z<=zmin)# & tissue=="Fibroblast")
# data_outliers<-data_AtLeastOneOutlier %>% filter(z>zmin) # & tissue=="Fibroblast") 
# data_nooutliers<-data_AtLeastOneOutlier %>% filter(z<=zmin)# & tissue=="Fibroblast")

bin_list<-c("ge0.9","0.75_to_0.9","0.25_to_0.75","0.1_to_0.25","le0.1")
constraint_types<-c("binned_genic","binned_regional")
statuses<-c("Case","Control","All")
filter_constraint<-c(TRUE,FALSE)
all_risks<-data.frame()
#general rr
for(this_constraint_type in constraint_types){
  for(this_group in statuses){
    for(filter_type in filter_constraint){
      if(this_group=="All"){
        group_outliers=data_outliers #%>% dplyr::filter(!!sym(this_constraint_type) !="NA")
        group_nonoutliers=data_nooutliers #%>% dplyr::filter(!!sym(this_constraint_type) !="NA")
      }else{
        group_outliers=data_outliers%>% dplyr::filter(status== this_group)  #%>% dplyr::filter(!!sym(this_constraint_type) !="NA")
        group_nonoutliers=data_nooutliers%>% dplyr::filter(status == this_group)  #%>% dplyr::filter(!!sym(this_constraint_type) !="NA")
      }
      if(filter_type==T){
        group_outliers =group_outliers%>% dplyr::filter(!!sym(this_constraint_type) !="NA")
        group_nonoutliers =group_nonoutliers%>% dplyr::filter(!!sym(this_constraint_type) !="NA")
      }
      for(this_bin in bin_list){
        
        exp_nn = nrow(group_nonoutliers %>% dplyr::filter(is.na(maf)) )#common
        #!!sym(this_constraint_type) forces dplyr to view this as a column correctly
        exp_ny = nrow(group_nonoutliers%>% dplyr::filter(!!sym(this_constraint_type) == this_bin)%>% dplyr::filter(!is.na(maf))) #   #rare
        exp_yn = nrow(group_outliers %>% dplyr::filter(is.na(maf))) #common
        exp_yy = nrow(group_outliers %>% dplyr::filter(!!sym(this_constraint_type) == this_bin) %>% dplyr::filter( !is.na(maf))) #%>%  #rare
        exptable = rbind(c(exp_nn,exp_ny),c(exp_yn,exp_yy))
        err = epitab(exptable, method = 'riskratio')
        this_risk<-data.frame(Risk = err$tab[2,5],
                              Lower = err$tab[2,6],
                              Upper = err$tab[2,7],
                              Pval = err$tab[2,8],
                              constraint_type=this_constraint_type,
                              constraint_level=this_bin,
                              filtered_constraint=filter_type,
                              status=this_group
        )
        all_risks<-rbind(all_risks,this_risk)
      }
      
      ### any score
      exp_nn = nrow(group_nonoutliers %>% dplyr::filter( is.na(maf) ) %>% dplyr::filter(!!sym(this_constraint_type) !="NA"))#common
      exp_ny = nrow(group_nonoutliers  %>% dplyr::filter(maf != "NA")%>% dplyr::filter(!!sym(this_constraint_type) !="NA")) # %>% dplyr::filter(binned_constraint == "low")  #rare
      exp_yn = nrow(group_outliers%>% dplyr::filter(is.na(maf))%>% dplyr::filter(!!sym(this_constraint_type) !="NA")) #common
      exp_yy = nrow(group_outliers  %>% dplyr::filter( maf != "NA")%>% dplyr::filter(!!sym(this_constraint_type) !="NA")) #%>% dplyr::filter(binned_constraint == "low") #rare
      
      exptable = rbind(c(exp_nn,exp_ny),c(exp_yn,exp_yy))
      err = epitab(exptable, method = 'riskratio')
      this_risk<-data.frame(Risk = err$tab[2,5],
                            Lower = err$tab[2,6],
                            Upper = err$tab[2,7],
                            Pval = err$tab[2,8],
                            constraint_type=this_constraint_type,
                            constraint_level="AnyConstraintScore",
                            filtered_constraint=filter_type,
                            status=this_group)
      all_risks<-rbind(all_risks,this_risk)
      ### all rv
      exp_nn = nrow(group_nonoutliers %>% dplyr::filter( is.na(maf) ) )#common
      exp_ny = nrow(group_nonoutliers%>% dplyr::filter(maf != "NA")) # %>% dplyr::filter(binned_constraint == "low")  #rare
      exp_yn = nrow(group_outliers %>% dplyr::filter(is.na(maf))) #common
      exp_yy = nrow(group_outliers  %>% dplyr::filter( maf != "NA")) #%>% dplyr::filter(binned_constraint == "low") #rare
      
      exptable = rbind(c(exp_nn,exp_ny),c(exp_yn,exp_yy))
      err = epitab(exptable, method = 'riskratio')
      this_risk<-data.frame(Risk = err$tab[2,5],
                            Lower = err$tab[2,6],
                            Upper = err$tab[2,7],
                            Pval = err$tab[2,8],
                            constraint_type=this_constraint_type,
                            constraint_level="NoFilter",
                            filtered_constraint=filter_type,
                            status=this_group)
      all_risks<-rbind(all_risks,this_risk)
    }
  }
}

##cadd
bin_list<-c("0","1","5","10","15","25")
for(this_constraint_type in "cadd"){
  for(this_group in statuses){
    if(this_group=="All"){
      group_outliers=data_outliers
      group_nonoutliers=data_nooutliers
    }else{
      group_outliers=data_outliers%>% dplyr::filter(status== this_group) 
      group_nonoutliers=data_nooutliers%>% dplyr::filter(status == this_group) 
    }
    for(this_bin in bin_list){
      
      exp_nn = nrow(group_nonoutliers %>% dplyr::filter( is.na(maf) ) )#common
      #!!sym(this_constraint_type) forces dplyr to view this as a column correctly
      exp_ny = nrow(group_nonoutliers%>% dplyr::filter(!!sym(this_constraint_type) >= this_bin)%>% dplyr::filter(maf != "NA")) #   #rare
      exp_yn = nrow(group_outliers %>% dplyr::filter(is.na(maf))) #common
      exp_yy = nrow(group_outliers %>% dplyr::filter(!!sym(this_constraint_type) >= this_bin) %>% dplyr::filter( maf != "NA")) #%>%  #rare
      exptable = rbind(c(exp_nn,exp_ny),c(exp_yn,exp_yy))
      err = epitab(exptable, method = 'riskratio')
      this_risk<-tryCatch({
        err = epitab(exptable, method = 'riskratio')
        data.frame(Risk = err$tab[2,5],
                              Lower = err$tab[2,6],
                              Upper = err$tab[2,7],
                              Pval = err$tab[2,8],
                              constraint_type=this_constraint_type,
                              constraint_level=this_bin,
                              status=this_group)
          },warning=function(w) {
                    print("WARNING")
                 data.frame(Risk = "NA",
                                Lower = "NA",
                                Upper = "NA",
                                Pval = "NA",
                                constraint_type=this_constraint_type,
                                constraint_level=this_bin,
                                status=this_group)
          #print(this_risk)
        })

      all_risks<-rbind(all_risks,this_risk)
    }
    
    ### any score
    exp_nn = nrow(group_nonoutliers %>% dplyr::filter(!is.na(!!sym(this_constraint_type) ) )%>% dplyr::filter( is.na(maf) ) )#common
    exp_ny = nrow(group_nonoutliers %>% dplyr::filter(!is.na(!!sym(this_constraint_type)  ) ) %>% dplyr::filter(maf != "NA")) # %>% dplyr::filter(binned_constraint == "low")  #rare
    exp_yn = nrow(group_outliers %>% dplyr::filter(!is.na(!!sym(this_constraint_type)  ) ) %>% dplyr::filter(is.na(maf))) #common
    exp_yy = nrow(group_outliers  %>% dplyr::filter(!is.na(!!sym(this_constraint_type)  ) ) %>% dplyr::filter( maf != "NA")) #%>% dplyr::filter(binned_constraint == "low") #rare
    
    exptable = rbind(c(exp_nn,exp_ny),c(exp_yn,exp_yy))
    this_risk<-tryCatch({
      err = epitab(exptable, method = 'riskratio')
     data.frame(Risk = err$tab[2,5],
                            Lower = err$tab[2,6],
                            Upper = err$tab[2,7],
                            Pval = err$tab[2,8],
                            constraint_type=this_constraint_type,
                            constraint_level="AnyScore",
                            status=this_group) 
      },warning=function(w) {
        data.frame(Risk = "NA",
                   Lower = "NA",
                   Upper = "NA",
                   Pval = "NA",
                   constraint_type=this_constraint_type,
                   constraint_level="AnyScore",
                   status=this_group) 
      }
     )
    all_risks<-rbind(all_risks,this_risk)
    # ### all rv
    # exp_nn = nrow(group_nonoutliers %>% dplyr::filter( is.na(maf) ) )#common
    # exp_ny = nrow(group_nonoutliers%>% dplyr::filter(maf != "NA")) # %>% dplyr::filter(binned_constraint == "low")  #rare
    # exp_yn = nrow(group_outliers %>% dplyr::filter(is.na(maf))) #common
    # exp_yy = nrow(group_outliers  %>% dplyr::filter( maf != "NA")) #%>% dplyr::filter(binned_constraint == "low") #rare
    # 
    # exptable = rbind(c(exp_nn,exp_ny),c(exp_yn,exp_yy))
    # err = epitab(exptable, method = 'riskratio')
    # this_risk<-data.frame(Risk = err$tab[2,5],
    #                       Lower = err$tab[2,6],
    #                       Upper = err$tab[2,7],
    #                       Pval = err$tab[2,8],
    #                       constraint_type=this_constraint_type,
    #                       constraint_level="NoFilter",
    #                       status=this_group)
    # all_risks<-rbind(all_risks,this_risk)
  }
}

#oulier only
group_outliers<-data_outliers %>% dplyr::filter(status=="Case")
group_nonoutliers<-data_outliers%>% dplyr::filter(status=="Control")
all_risks<-data.frame()
for(this_constraint_type in constraint_types){
    for(this_bin in bin_list){
      
      exp_nn = nrow(group_nonoutliers%>%dplyr::filter( is.na(maf) ) )#common
      #!!sym(this_constraint_type) forces dplyr to view this as a column correctly
      exp_ny = nrow(group_nonoutliers%>% dplyr::filter(!!sym(this_constraint_type) == this_bin)%>% dplyr::filter(maf != "NA")) #   #rare
      exp_yn = nrow(group_outliers %>% dplyr::filter(is.na(maf))) #common
      exp_yy = nrow(group_outliers %>% dplyr::filter(!!sym(this_constraint_type) == this_bin) %>% dplyr::filter( maf != "NA")) #%>%  #rare
      exptable = rbind(c(exp_nn,exp_ny),c(exp_yn,exp_yy))
      tryCatch({
        err = epitab(exptable, method = 'riskratio')
        this_risk<-data.frame(Risk = err$tab[2,5],
                              Lower = err$tab[2,6],
                              Upper = err$tab[2,7],
                              Pval = err$tab[2,8],
                              constraint_type=this_constraint_type,
                              constraint_level=this_bin,
                              status="CaseVsControl"
        )
        all_risks<-rbind(all_risks,this_risk)}, 
      warning=function(w) {
        print("WARNING")
        this_risk<-data.frame(Risk = "NA",
                              Lower = "NA",
                              Upper = "NA",
                              Pval = "NA",
                              constraint_type=this_constraint_type,
                              constraint_level=this_bin,
                              status="CaseVsControl")
      })
    }
  
  ### any score
  exp_nn = nrow(group_nonoutliers %>% dplyr::filter(!!sym(this_constraint_type) !="NA")%>% dplyr::filter( is.na(maf) ) )#common
  exp_ny = nrow(group_nonoutliers %>% dplyr::filter(!!sym(this_constraint_type) !="NA") %>% dplyr::filter(maf != "NA")) # %>% dplyr::filter(binned_constraint == "low")  #rare
  exp_yn = nrow(group_outliers %>% dplyr::filter(!!sym(this_constraint_type) !="NA") %>% dplyr::filter(is.na(maf))) #common
  exp_yy = nrow(group_outliers  %>% dplyr::filter(!!sym(this_constraint_type) !="NA") %>% dplyr::filter( maf != "NA")) #%>% dplyr::filter(binned_constraint == "low") #rare
  
  exptable = rbind(c(exp_nn,exp_ny),c(exp_yn,exp_yy))
  err = epitab(exptable, method = 'riskratio')
  this_risk<-data.frame(Risk = err$tab[2,5],
                        Lower = err$tab[2,6],
                        Upper = err$tab[2,7],
                        Pval = err$tab[2,8],
                        constraint_type=this_constraint_type,
                        constraint_level="AnyConstraintScore",
                        status="CaseVsControl")
  all_risks<-rbind(all_risks,this_risk)
  ### all rv
  exp_nn = nrow(group_nonoutliers %>% dplyr::filter( is.na(maf) ) )#common
  exp_ny = nrow(group_nonoutliers%>% dplyr::filter(maf != "NA")) # %>% dplyr::filter(binned_constraint == "low")  #rare
  exp_yn = nrow(group_outliers %>% dplyr::filter(is.na(maf))) #common
  exp_yy = nrow(group_outliers  %>% dplyr::filter( maf != "NA")) #%>% dplyr::filter(binned_constraint == "low") #rare
  
  exptable = rbind(c(exp_nn,exp_ny),c(exp_yn,exp_yy))
  err = epitab(exptable, method = 'riskratio')
  this_risk<-data.frame(Risk = err$tab[2,5],
                        Lower = err$tab[2,6],
                        Upper = err$tab[2,7],
                        Pval = err$tab[2,8],
                        constraint_type=this_constraint_type,
                        constraint_level="NoFilter",
                        status="CaseVsControl")
  all_risks<-rbind(all_risks,this_risk)
}



# to_plot<-all_risks %>% mutate(constraint_level=fct_relevel(constraint_level,"le0.1","0.1_to_0.25","0.25_to_0.75","0.75_to_0.9","ge0.9","AnyConstraintScore","NoFilter"))%>% dplyr::filter(status!="Control")
to_plot<-all_risks %>% dplyr::filter(constraint_type=="binned_regional" ) %>%
  mutate(constraint_level=fct_relevel(constraint_level,"le0.1","0.1_to_0.25","0.25_to_0.75","0.75_to_0.9","ge0.9","AnyConstraintScore","NoFilter"))
  #mutate(constraint_level=fct_relevel(constraint_level,"0","1","5","10","15","25"))
  #  mutate(constraint_level=fct_relevel(constraint_level,"le0.1","0.1_to_0.25","0.25_to_0.75","0.75_to_0.9","ge0.9","AnyConstraintScore","NoFilter"))
#%>% dplyr::filter(status!="Control")

to_plot$Risk<-as.numeric(to_plot$Risk)
to_plot$Lower<-as.numeric(to_plot$Lower)
to_plot$Upper<-as.numeric(to_plot$Upper)

#"ge0.9","0.75_to_0.9","0.25_to_0.75","0.1_to_0.25","le0.1")
ggplot(to_plot,aes(x=status,y=Risk,fill=constraint_level,group=constraint_level))+ 
  ggtitle(paste0("RR, z=",zmin))+
  # ggtitle("RR of having nearby RV with x constraint score given CASE/CONTROL status")+
  facet_wrap(~filtered_constraint,ncol=2)+
  #scale_fill_manual(values=c("#c2b2b5","#f7b0bc","#ed7e93","#f05471","#eb2d50","#913c3c","#e8133a"))+
   scale_fill_manual(values=c("#1f499c","#7594d1","#75d1c3","#75d180","#1e872a","#155940","#913c3c"))+
  theme(axis.text.x = element_text(angle = 45,  hjust=1))+
  #scale_x_discrete(expand = expansion(add = 1))+
  theme_linedraw(base_size=15)+#ylim(c(0.9,3.5))+
  xlab("Tissue Type")+
  ylab("Relative Risk")+
  geom_crossbar(aes(x=status,y=Risk,ymin=Lower,ymax=Upper),position="dodge") +
  geom_hline(yintercept=1,color="red")

ggplot(data_outliers,aes(x=cadd,y=genic_score,alpha=0.1))+
  geom_point()+ggtitle("Z vs. genic score")
ggplot(data_outliers,aes(x=sample,y=regional_score))+geom_density_ridges()

regional_corr<-lm(data=data_outliers,z~regional_score)
genic_corr<-lm(data=data_outliers,z~genic_score)



data_outliers$genic_rank<-frankv(x=-data_outliers$genic_score,na.last="keep")/length(na.omit(data_outliers$genic_score))
data_outliers$regional_rank<-frankv(x=-data_outliers$regional_score,na.last="keep")/length(na.omit(data_outliers$regional_score))

 # data_outliers$regional_rank<-as.numeric(rank(-data_outliers$regional_score)/nrow(data_outliers))

data_outliers_case<-data_outliers%>%filter(status=="Case")#%>%filter(z>3)
data_outliers_control<-data_outliers%>%filter(status=="Control")# %>%filter(z>3)
qqplot(y=data_outliers_case$genic_rank,
       x=data_outliers_control$genic_rank,
       ylab='Case outlier genes',
       xlab = 'Control outlier genes') 
title("genic_rank (Z=2) ")
abline(0, 1, col = 'red')


solved=fread("/Volumes/groups/smontgom/shared/UDN/ReferenceFiles/solved_cases.txt")
all_solved_data<-data.frame()
for(rd in solved$RD){
  print(rd)
  solved_ensg<-solved %>% filter(RDID==rd) %>% pull(ENSG)
  solved_data<-data %>%filter(sample==rd)
  solved_data$genic_rank<-frankv(x=-solved_data$genic_score)
  solved_data$regional_rank<-frankv(x=-solved_data$regional_score)
  solved_data$cadd_rank<-frankv(x=-solved_data$cadd)
  ##outlier only
  solved_data$genic_rank_outlierONLY<-frankv(x=-(solved_data %>% 
                                                   mutate(regional_outlier=ifelse(z>2,genic_score,NA))
                                                 %>%pull(regional_outlier)))
  solved_data$regional_rank_outlierONLY<-frankv(x=-(solved_data %>% 
                                                      mutate(regional_outlier=ifelse(z>2,regional_score,NA))
                                                    %>%pull(regional_outlier)))
  solved_data$cadd_rank_outlierONLY<-frankv(x=-(solved_data %>% 
                                                  mutate(regional_outlier=ifelse(z>2,cadd,NA))
                                                %>%pull(regional_outlier)))
  ##now pull solved ensg
  solved_data<-solved_data %>%filter(ensg==solved_ensg)
  solved_data$SplicingOutlier<-solved %>% filter(RDID==rd)  %>% dplyr::pull(Splicing_Outlier)
  solved_data$eOutlier<-solved %>% filter(RDID==rd)  %>% dplyr::pull(Expression_Outlier)
  solved_data$aseOutlier<-solved %>% filter(RDID==rd)  %>% dplyr::pull(ASE_Outlier)
  all_solved_data<-rbind(all_solved_data,solved_data)
  }

ggplot(data=all_solved_data,aes(x=z,y=cadd,color=SplicingOutlier,alpha=0.9))+
  ggtitle("Rank of Score vs. Z score")+
  scale_color_manual(values=c("blue","orange"))+
  geom_point()+theme_linedraw(base_size=15)

###solved ranking





