# Principal Components Analysis

# install adegenet if you haven't already
#install.packages('adegenet')
install.packages('adegenet', dep=TRUE)

# load adegenet package
library(adegenet)
?adegenet

# read in data as structure file (.str)
# adjust file name, n.ind, n.loci; n.ind = number of individuals; n.loci = number of unlinked SNPs
data <- read.structure("Atytest.stru", n.ind=85, n.loc=5421, onerowperind=FALSE, col.lab=1, col.pop=0, col.others=NULL, row.marknames=NULL, NA.char="-9", pop=NULL, ask=FALSE, quiet=FALSE)
data
View(data)
#

dat = readLines("Groups_Structure2.csv")
data.frame(dat)
dapc2 <- dapc(data,pop = dat)
scatter(dapc2, scree.da = FALSE, bg = "white", pch = 20, cell = 0, cstar = 0, col = myCol, 
        solid = .4, cex = 2, clab = 0, leg= FALSE, txt.leg =paste ("Cluster",1:6))


# put data values into a matrix
#data_scaled <- scaleGen(data, center=FALSE, scale=FALSE, NA.method=c("zero"))

test_data <- find.clusters(data, max.n.clust = 20)
test_data
#test_group <- find.clusters(data_scaled, max.n.clust=20) 
names(test_data)
head(test_data$Kstat, 4)
head(test_data$stat, 4)
head(test_data$size, 4)
head(test_data$grp, 4)

scatter(test_data)


test_data
table<-table(pop(test_data) , test_data$grp)



dapc1 <- dapc(data, test_data$grp)
#dapc1 <- dapc(data, nmax = 12)

dapc1
scatter(dapc1)
scatter(dapc1,posi.da="bottomright", bg="white", pch = 17:22)

myCol <- c("red", "darkblue", "darkgreen", "purple", "orange", "skyblue")
scatter(dapc1, bg="white", pch = 17:22, cstar = 0, col=myCol, 
        scree.pca=TRUE, psoi.pca="bottomleft" )

scatter(dapc1, scree.da = FALSE, bg = "white", pch = 20, cell = 0, cstar = 0, col = myCol, 
        solid = .4, cex = 2, clab = 0, leg= TRUE, txt.leg =paste ("Cluster",1:5))

scatter(dapc2, scree.da = FALSE, bg = "white", pch = 20, cell = 0, cstar = 0, col = myCol, 
        solid = .4, cex = 2, clab = 0, leg= FALSE, txt.leg =paste ("Cluster",1:6))




