pdf('UHRRep3.junctionSaturation_plot.pdf')
x=c(5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100)
y=c(974,1474,1833,2077,2267,2403,2523,2623,2717,2780,2838,2905,2955,2994,3040,3077,3110,3140,3163,3185)
z=c(991,1506,1890,2151,2350,2499,2626,2730,2836,2913,2980,3053,3113,3164,3220,3263,3308,3346,3382,3412)
w=c(17,32,57,74,83,96,103,107,119,133,142,148,158,170,180,186,198,206,219,227)
m=max(3,3,0)
n=min(0,0,0)
plot(x,z/1000,xlab='percent of total reads',ylab='Number of splicing junctions (x1000)',type='o',col='blue',ylim=c(n,m))
points(x,y/1000,type='o',col='red')
points(x,w/1000,type='o',col='green')
legend(5,3, legend=c("All junctions","known junctions", "novel junctions"),col=c("blue","red","green"),lwd=1,pch=1)
dev.off()
