args = commandArgs(TRUE)

require(ggplot2)
require(RColorBrewer)
library(gridExtra)


dat = read.table(args[1],sep="\t",h=T)
outfile = args[2]

# Cohort assignment
dat$cohort = "CAM"
dat[grepl("TUM_",dat$IID) == T,"cohort"] = "TUM"

# Scaling
for (j in 2:11){
	dat[,j] = scale(dat[,j])
}

c_values = c("darkblue","red")

# Plo
leg_rows=1
pl1 <- ggplot(data=dat, aes_string(x="C1", y="C2", colour="cohort")) +geom_point(size=I(1)) +theme_bw() +theme(legend.position="bottom",panel.grid.major=element_line(colour="grey60"),panel.grid.minor=element_blank()) +scale_colour_manual(values=c_values) +scale_x_continuous(breaks=seq(-100,100,1)) +scale_y_continuous(breaks=seq(-100,100,1)) +xlab("PC 1") +ylab("PC 2") +guides(col=guide_legend(nrow=leg_rows,override.aes=list(size=I(2))))
pl2 <- ggplot(data=dat, aes_string(x="C3", y="C4", colour="cohort")) +geom_point(size=I(1)) +theme_bw() +theme(legend.position="bottom",panel.grid.major=element_line(colour="grey60"),panel.grid.minor=element_blank()) +scale_colour_manual(values=c_values) +scale_x_continuous(breaks=seq(-100,100,1)) +scale_y_continuous(breaks=seq(-100,100,1)) +xlab("PC 3") +ylab("PC 4") 
pl3 <- ggplot(data=dat, aes_string(x="C5", y="C6", colour="cohort")) +geom_point(size=I(1)) +theme_bw() +theme(legend.position="bottom",panel.grid.major=element_line(colour="grey60"),panel.grid.minor=element_blank()) +scale_colour_manual(values=c_values) +scale_x_continuous(breaks=seq(-100,100,1)) +scale_y_continuous(breaks=seq(-100,100,1)) +xlab("PC 5") +ylab("PC 6") 
pl4 <- ggplot(data=dat, aes_string(x="C7", y="C8", colour="cohort")) +geom_point(size=I(1)) +theme_bw() +theme(legend.position="bottom",panel.grid.major=element_line(colour="grey60"),panel.grid.minor=element_blank()) +scale_colour_manual(values=c_values) +scale_x_continuous(breaks=seq(-100,100,1)) +scale_y_continuous(breaks=seq(-100,100,1)) +xlab("PC 7") +ylab("PC 8") 

png(outfile)
grid.arrange(pl1,pl2,pl3,pl4,nrow=2)
dev.off()