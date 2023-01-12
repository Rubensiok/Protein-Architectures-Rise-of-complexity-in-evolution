library(mixOmics)
library('sgPLS')
library(stringr)
data_tdidf <- read.csv('/home/ruben/SINCRONIZACION_CABD/bigrams_n_dataset/corrected/refract/uniarch_and_bigrams_tfidf.csv',row.names = 1, header = TRUE)

data1 <- data_tdidf

data1<-data1[data1$superkingdom != 'Eukaryota',]
data1<-data1[data1$superkingdom == 'Bacteria',]

myvars <- names(data1) %in% c("superkingdom", "phylum", "classes")
newdata <- data1[!myvars]
newdata = newdata[,colSums(newdata) > 0]

grac = c('Proteobacteria','Spirochaetota','Bacteroidota','Acidobacteriota',
         'Planctomycetota','Verrucomicrobiota','Chlamydiota','Bdellovibrionota','Campylobacterota','Desulfobacterota')
terra = c('Actinobacteriota','Firmicutes','Cyanobacteria','Chloroflexota',
          'Aquificota','Deinococcota','Synergistota','Fusobacteriota','Thermotogota')

i = 'Eukaryota'
Z<-data1$superkingdom
Z[(Z!=i)] = 'Prokaryotes'
Y <- as.factor(Z)

#########
MyResult.plsda2 <- splsda(newdata,Y, ncomp=3)#,near.zero.var=TRUE)
MyPerf.plsda <- perf(MyResult.plsda2, validation = "Mfold", folds = 5, 
                     progressBar = TRUE, nrepeat = 50) # we suggest nrepeat = 50
par(mfrow=c(1,1))
plot(MyPerf.plsda, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")
MyPerf.plsda$choice.ncomp
list.keepX <- c(5:10,  seq(10, 100, 10))

set.seed(30) # for reproducbility in this vignette, otherwise increase nrepeat
tune.splsda.srbct <- tune.splsda(newdata, Y, ncomp =4, # we suggest to push ncomp a bit more, e.g. 4
                                 validation = 'Mfold',
                                 folds = 5, 
                                 dist = 'max.dist', 
                                 progressBar = TRUE,
                                 measure = "overall", 
                                 test.keepX = list.keepX,
                                 nrepeat = 50,
                                 cpus=5)
# cpus=2)   # we suggest nrepeat = 50

error <- tune.splsda.srbct$error.rate
ncomp <- tune.splsda.srbct$choice.ncomp$ncomp # optimal number of components based on t-tests on the error rate
ncomp

select.keepX <- tune.splsda.srbct$choice.keepX[1:ncomp]  # optimal number of variables to select
# select.keepX<-c(100,20,100)
select.keepX

plot(tune.splsda.srbct, col = color.jet(4),title = 'mahalanobis 5 ncomp')

set.seed(30) # for reproducbility in this vignette, otherwise increase nrepeat
tune.splsda.srbct <- tune.splsda(newdata, Y, ncomp=4, # we suggest to push ncomp a bit more, e.g. 4
                                 validation = 'Mfold',
                                 folds = 5, 
                                 dist = 'centroids.dist', 
                                 progressBar = TRUE,
                                 measure = "overall", 
                                 test.keepX = list.keepX,
                                 nrepeat = 50,
                                 cpus=4)

error <- tune.splsda.srbct$error.rate
ncomp <- tune.splsda.srbct$choice.ncomp$ncomp # optimal number of components based on t-tests on the error rate
ncomp

select.keepX <- tune.splsda.srbct$choice.keepX[1:ncomp] 

plot(tune.splsda.srbct, col = color.jet(4),title = 'mahalanobis 5 ncomp')

MyResult.splsda.final <- splsda(newdata, Y, ncomp = ncomp, keepX = select.keepX)#,near.zero.var=TRUE)

par(mfrow=c(1,ncomp+1))    # set the plotting area into a 1*2 array
plotLoadings(MyResult.splsda.final, contrib = 'max', method = 'mean',size.name=0.65,comp=1,plot=TRUE,ndisplay=select.keepX['comp1'],legend = FALSE)
plotLoadings(MyResult.splsda.final, contrib = 'max', method = 'mean',size.name=0.65,comp=2,plot=TRUE,ndisplay=select.keepX['comp2'],legend = FALSE)
plotLoadings(MyResult.splsda.final, contrib = 'max', method = 'mean',size.name=0.65,comp=3,plot=TRUE,ndisplay=select.keepX['comp3'],legend = TRUE)

for (x in 1:2){
  pl = plotLoadings(MyResult.splsda.final, contrib = 'max', method = 'mean',size.name=0.65,comp=x,ndisplay=select.keepX[x],plot=FALSE)
  write.csv(pl,paste("/home/ruben/SINCRONIZACION_CABD/bigrams_n_dataset/corrected/refract/texto/nuevos_plots/sPLS-DA binary_koonin/50 repeats/uniarch_and_bigrams_TDIDF_Cyano_vs_Terra_", x, ".csv",sep = ''), row.names = TRUE)
}


plotIndiv(MyResult.splsda.final, ind.names = FALSE, legend=TRUE,
          ellipse = TRUE, title="sPLS-DA", legend.position='bottom',legend.title = '')