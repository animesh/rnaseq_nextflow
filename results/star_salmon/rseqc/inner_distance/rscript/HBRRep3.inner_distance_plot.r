out_file = 'HBRRep3'
pdf('HBRRep3.inner_distance_plot.pdf')
fragsize=rep(c(-247.5,-242.5,-237.5,-232.5,-227.5,-222.5,-217.5,-212.5,-207.5,-202.5,-197.5,-192.5,-187.5,-182.5,-177.5,-172.5,-167.5,-162.5,-157.5,-152.5,-147.5,-142.5,-137.5,-132.5,-127.5,-122.5,-117.5,-112.5,-107.5,-102.5,-97.5,-92.5,-87.5,-82.5,-77.5,-72.5,-67.5,-62.5,-57.5,-52.5,-47.5,-42.5,-37.5,-32.5,-27.5,-22.5,-17.5,-12.5,-7.5,-2.5,2.5,7.5,12.5,17.5,22.5,27.5,32.5,37.5,42.5,47.5,52.5,57.5,62.5,67.5,72.5,77.5,82.5,87.5,92.5,97.5,102.5,107.5,112.5,117.5,122.5,127.5,132.5,137.5,142.5,147.5,152.5,157.5,162.5,167.5,172.5,177.5,182.5,187.5,192.5,197.5,202.5,207.5,212.5,217.5,222.5,227.5,232.5,237.5,242.5,247.5),times=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,123,2083,2467,3228,3770,4310,4733,4966,5125,5076,5022,4955,4656,4571,4249,3992,3685,3511,3341,2978,2909,2555,2334,2202,2026,1816,1719,1565,1415,1272,1148,1081,942,866,781,693,687,618,579,511,455,429,381,355,349,305,238,236,213,213,161,170,142,144,113,87,83,69,76,73,58,17,12,6,6,5,5,3,7,3,5))
frag_sd = sd(fragsize)
frag_mean = mean(fragsize)
frag_median = median(fragsize)
write(x=c("Name","Mean","Median","sd"), sep="	", file=stdout(),ncolumns=4)
write(c(out_file,frag_mean,frag_median,frag_sd),sep="	", file=stdout(),ncolumns=4)
hist(fragsize,probability=T,breaks=100,xlab="mRNA insert size (bp)",main=paste(c("Mean=",frag_mean,";","SD=",frag_sd),collapse=""),border="blue")
lines(density(fragsize,bw=10),col='red')
dev.off()
