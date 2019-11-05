library(dplyr)
library(ggplot2)




###parameters

cumfreqcut=0.1
inversecumfreqcut=1-cumfreqcut
minval=5
scaleall=1
filtercumfreq=1
filtercumfreq2=1 
medianscale=1




ShortenGuideId<-function(fcounts){
    fcounts$guide_id=gsub("_othernuclear$","",fcounts$guide_id,perl=TRUE)
    fcounts$guide_id=gsub("_TransFac$","",fcounts$guide_id,perl=TRUE)
    fcounts$guide_id=gsub("_Epigenetic_Reg$","",fcounts$guide_id,perl=TRUE)
    fcounts$guide_id=gsub("_dnaRepair$","",fcounts$guide_id,perl=TRUE)
    fcounts$guide_id=gsub("_Handpicked$","",fcounts$guide_id,perl=TRUE)
    fcounts$guide_id=gsub("_OT_predicted$","",fcounts$guide_id,perl=TRUE)
    return(fcounts)
}


CountFilter<-function(fcounts,fminval,fcumfreqcut,expname,flibsize,fscalevar,filtercumfreq1){
    
    fscalefactor<-flibsize$scalefactor[flibsize$expid == expname]

    ## filter unscaled count-table at fixed abundance cutoff
    fcounts<-fcounts[fcounts$guide_umi_count>0,]
    fcounts<-fcounts[fcounts$guide_umi_count>=fminval,]

    ## filter unscaled count-table at fixed abundance cutoff
    if(filtercumfreq1){
        
        inversefcumfreqcut=1-fcumfreqcut
        fcounts<-fcounts[order(fcounts$guide_id, -fcounts$guide_umi_count),]
         fcounts <- fcounts %>%                                   
            group_by(guide_id) %>%                       
                mutate(guide_umi_count_freq = guide_umi_count/sum(guide_umi_count, na.rm=TRUE)) # %>%     
        fcounts <- fcounts %>%                                   
            group_by(guide_id) %>%                       
                mutate(cumulative_percentage = cumsum(guide_umi_count)/sum(guide_umi_count, na.rm=TRUE)) # %>%     
        tol = 1e-15 # tolerance for float precision
        fcounts <- fcounts %>%                                   
            group_by(guide_id) %>%                       
                mutate(cumulative_percentage_cutoff = max(inversefcumfreqcut,(min(cumulative_percentage)+tol)))
        fcounts <- fcounts %>%                                   
            group_by(guide_id) %>%                       
                mutate(consider4mincount_cutoff = ifelse(cumulative_percentage<=cumulative_percentage_cutoff,guide_umi_count,0))
        fcounts <- fcounts %>%                                   
            group_by(guide_id) %>%                       
                mutate(cumfreq_count_cutoff = min(consider4mincount_cutoff[consider4mincount_cutoff>0]))

        ## filter
        fcounts<-fcounts[fcounts$guide_umi_count>=fcounts$cumfreq_count_cutoff | fcounts$guide_umi_count_freq>=fcumfreqcut ,]    }
   ## normalize by library size
    if(fscalevar){
        ##scale using scaling factor based on all experiments to keep counts close to original values
        fcounts$guide_umi_count<-fcounts$guide_umi_count/fscalefactor
    }
    else{
        ##scale per sample - CPM like
        fcounts$guide_umi_count<-1000000*fcounts$guide_umi_count/sum(fcounts$guide_umi_count)
    }
   
   
   fcounts$gene = as.character(lapply(strsplit(as.character(fcounts$guide_id), split="_"), "[", 1))
   

   return(fcounts)
}




CumFilterMerged<-function(fcounts){

    cumfreqcutmin<-0.01
    inversecumfreqcutmin=1-cumfreqcutmin

    colnames.fcounts<-colnames(fcounts)
    
    fcounts<-fcounts[order(fcounts$guide_id.condition, -fcounts$guide_umi_count),]
    
    fcounts <- fcounts %>%                                   
        group_by(guide_id.condition) %>%                       
            mutate(fcounts_cumulative_percentage = cumsum(guide_umi_count)/sum(guide_umi_count, na.rm=TRUE)) # %>%     

    tol = 1e-15 # tolerance for float precision
    fcounts <- fcounts %>%                                   
        group_by(guide_id.condition) %>%                       
            mutate(fcounts_cumulative_percentage_cutoff = max(inversecumfreqcutmin,(min(fcounts_cumulative_percentage)+tol)))
    
    fcounts <- fcounts %>%                                   
        group_by(guide_id.condition) %>%                       
            mutate(fcountsconsider4mincount_cutoff = ifelse(fcounts_cumulative_percentage<=fcounts_cumulative_percentage_cutoff,guide_umi_count,0))
    
    fcounts <- fcounts %>%                                   
        group_by(guide_id.condition) %>%                       
            mutate(fcounts_cumfreq_count_cutoff = min(fcountsconsider4mincount_cutoff[fcountsconsider4mincount_cutoff>0]))
    
    fcounts <- fcounts %>%                                   
        group_by(guide_id.condition) %>%                       
            mutate(fcounts_guide_umi_count_freq = guide_umi_count/sum(guide_umi_count)) # %>%     

    ## do filtering
    fcounts<-fcounts[fcounts$guide_umi_count>=fcounts$fcounts_cumfreq_count_cutoff | fcounts$fcounts_guide_umi_count_freq>=cumfreqcutmin,]
    
    ## keep only inital columns
    fcounts<-fcounts[,colnames.fcounts]

    return(fcounts)
}

ParameterToFileName<-function(fminval,ffiltercumfreq,ffiltercumfreq2,fscaleall,fmedianscale){
    fextension<-paste0("minval",fminval,ifelse(ffiltercumfreq, paste0("_CumFreqFilterSamp",cumfreqcut), ""),ifelse(fmedianscale,"_MEDscale", ""),ifelse(ffiltercumfreq2,"_CumFreqFilterMerged", ""),ifelse(fscaleall,"_scaleALL", ""))
    return(fextension)
}






















## make output filename to include parameter settings

outfilenamevar<-ParameterToFileName(minval,filtercumfreq,filtercumfreq2,scaleall,medianscale)

### data input

results.dir<-Sys.getenv(c("RESULTS_DIR"))

targets.file<-Sys.getenv(c("TARGETS_FILE"))
targets <- read.table(targets.file,sep="",header=FALSE)

libsize.file<-paste0(results.dir,"/libsize.txt")
libsize<-read.table(libsize.file,header=TRUE)

ids2plot.file<-Sys.getenv(c("WANTEDGUIDES_FILE"))
ids2plot.ids=read.table(ids2plot.file)
ids2plot.ids<-as.character(levels(ids2plot.ids[[1]]))




exp.list <- list()







                                        #import counts begin
for (i in 1:nrow(targets))
    {
        name.expid<-paste0(targets[i,2],"_",targets[i,3])
        name.exp<-targets[i,2]
        name.condition<-ifelse(grepl("^E",name.exp),"exp","ctr")
        
        count.filename=paste0(results.dir,"/",name.expid,"_guide_UMI_counts.txt")
        
        print(paste0("reading ",count.filename))
        stopifnot(file.exists(count.filename))
        
        counts<-read.table(count.filename,header=FALSE,sep=" ")
        colnames(counts)<-c("guide_id","umi_seq","guide_umi_count")

        counts$guide_umi_count<-as.numeric(counts$guide_umi_count)
        
        ##flag guideids that are to be plotted
        counts$wanted<-ifelse(counts$guide_id %in% ids2plot.ids,1,0)
        
        ##filter: keep counts of 5 and above; for each guide keep count values up to the cumulative percentage of 90%
        counts<-CountFilter(counts,minval,cumfreqcut,name.expid,libsize,scaleall,filtercumfreq)

        ## set meta info to allow summarization of guide_umi_counts at various levels
        counts$expid<-name.expid #e.g. E1_A
        counts$experiment<-name.exp #e.g. E1
        counts$condition<-name.condition #e.g. exp
                
        ## median scale remaining guide_umi_counts to allow visualization of guide_umi_counts relative to the median of the experiment
        counts.median<-median(counts$guide_umi_count,na.rm=TRUE)
        counts$guide_umi_count.medscaled<-as.numeric(as.character(counts$guide_umi_count/counts.median))


        ## keep only columns to be used on
        counts<-counts[,c("guide_id","expid","experiment","condition","guide_umi_count","guide_umi_count.medscaled","gene","wanted")]

        counts$guide_id<-as.character(counts$guide_id)

        ## store count table in list for all exp counts
        exp.list[[name.expid]]<-counts
        
    }
                                        #import counts end



## merge all exp counts into one table
counts.allexp<-do.call(rbind, exp.list)

## mild filter for merged data- remove outlier by keeping counts that show cumulative frequency of up to 99%
counts.allexp$guide_id.condition<-paste0(counts.allexp$guide_id,counts.allexp$condition)
counts.allexp<-CumFilterMerged(counts.allexp)

if (medianscale){
    counts.allexp$guide_umi_count<-counts.allexp$guide_umi_count.medscaled
}



## keep only data from exp and not control 
counts.allexp<-counts.allexp[ counts.allexp$condition == "exp",]

## for plotting set guide_id for not selected ids to OTHER
counts.allexp$guide_id[counts.allexp$wanted==0]<-"OTHER"



## determine order for plotting
counts.allexp <- counts.allexp %>%                                 
    group_by(gene) %>%                       
    mutate(median_perg = median(guide_umi_count,na.rm=TRUE)) # %>%     

counts.allexp <- counts.allexp %>%                                 
    group_by(guide_id) %>%                       
    mutate(median_perseqid = median(guide_umi_count,na.rm=TRUE)) # %>%       

counts.allexp<-counts.allexp[order(counts.allexp$median_perg,counts.allexp$median_perseqid),]
geneorder<-unique(counts.allexp$guide_id)

counts.allexp$guide_id <- factor(counts.allexp$guide_id, levels = geneorder)




## plot counts for selected targets
outfile<-paste0(results.dir,"/",outfilenamevar,"_genesort.pdf")

p<-ggplot(counts.allexp, aes(x=guide_id, y=log2(guide_umi_count))) + 
    geom_boxplot(width=0.5,outlier.colour = "black",outlier.size = 0.3,fill="#56B4E9", color="black")  +  
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text.x = element_text(face="bold", color="#993333", size=10, angle=90))   +  xlab("guide")

ggsave(outfile,p,width=20,height=5,limitsize = FALSE)


print(paste0("done plotting to ",outfile))
