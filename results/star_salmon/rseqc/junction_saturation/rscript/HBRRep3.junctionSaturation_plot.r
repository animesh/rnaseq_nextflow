pdf('HBRRep3.junctionSaturation_plot.pdf')
x=c(5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100)
y=c(835,1295,1591,1813,1974,2108,2231,2344,2442,2529,2596,2656,2709,2775,2815,2850,2891,2922,2956,2980)
z=c(844,1312,1624,1854,2033,2172,2306,2429,2535,2629,2709,2780,2848,2923,2966,3011,3056,3097,3136,3169)
w=c(9,17,33,41,59,64,75,85,93,100,113,124,139,148,151,161,165,175,180,189)
m=max(2,3,0)
n=min(0,0,0)
plot(x,z/1000,xlab='percent of total reads',ylab='Number of splicing junctions (x1000)',type='o',col='blue',ylim=c(n,m))
points(x,y/1000,type='o',col='red')
points(x,w/1000,type='o',col='green')
legend(5,3, legend=c("All junctions","known junctions", "novel junctions"),col=c("blue","red","green"),lwd=1,pch=1)
dev.off()
