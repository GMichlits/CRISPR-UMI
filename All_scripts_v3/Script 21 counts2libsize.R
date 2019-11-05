library(dplyr)


results.dir<-Sys.getenv(c("RESULTS_DIR"))
targets.file<-Sys.getenv(c("TARGETS_FILE"))

libsize.file<-paste0(results.dir,"/libsize.txt")


targets <- read.table(targets.file,sep="",header=FALSE)



exp.list <- list()

for (i in 1:nrow(targets))
            {
                expid<-paste0(targets[i,2],"_",targets[i,3])
                count.filename=paste0(results.dir,"/",expid,"_guide_UMI_counts.txt")
                print(paste0("reading ",count.filename))
                stopifnot(file.exists(count.filename))
                counts<-read.table(count.filename,header=FALSE,sep=" ")
                colnames(counts)<-c("guide_id","umi_seq","guide_umi_count")
                counts<-counts[counts$guide_umi_count>0,]
                counts$expid<-expid
                counts$guide_id<-as.character(counts$guide_id)
                exp.list[[expid]]<-counts

            }

counts.allexp<-do.call(rbind, exp.list)

libsize <- counts.allexp %>%
    group_by(expid) %>%
    summarize(
        total_count = n(),
        guide_umi_sum = sum(guide_umi_count, na.rm=TRUE)
        )

libsize$max=max(libsize$guide_umi_sum)
libsize$scalefactor=format(round(libsize$guide_umi_sum/libsize$max, digits =8), nsmall = 8)
write.table(libsize,libsize.file,row.names=FALSE, quote=FALSE,sep="\t")

print(paste0("done writing ",libsize.file))
