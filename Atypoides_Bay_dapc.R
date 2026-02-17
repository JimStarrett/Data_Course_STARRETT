# Principal Components Analysis
#analysis of multivariate data, complex data with many dimensions reduced to few dimensions

#clone data and scripts from github
#git clone https://github.com/JimStarrett/Data_Course_STARRETT/DAPC_test.git


#Use adegenet package
#adegenet package: investigation of population structure from DNA sequence data
#follow DAPC tutorial (Jombart and Collins, June 5, 2022)

#Data type to analyze
#SNP data, Single Nucleotide Polymorphism, in Structure format
#Structure format: alleles represented by numbers, rows=individuals, columns=DNA marker 


# install adegenet if you haven't already
#install.packages('adegenet')
#install.packages('adegenet', dep=TRUE)

# load adegenet package
library(adegenet)
?adegenet

######### Read in data 

# read in data as structure file (.str)
#data is unlinked SNPs (1 SNP per locus) from turret spiders around SF Bay Area
#adjust file name, n.ind, n.loci; n.ind = number of individuals; n.loci = number of unlinked SNPs
Bay_stru_data <- read.structure("Atytest.stru", n.ind=85, n.loc=5421, onerowperind=FALSE, col.lab=1, col.pop=0, col.others=NULL, row.marknames=NULL, NA.char="-9", pop=NULL, ask=FALSE, quiet=FALSE)
Bay_stru_data               #genind data
#View(Bay_stru_data)
#goal is to assign individuals to populations (or clusters or groups) based on genetic variation
#optimizing variation between groups and minimizing variation within groups


########### Find hypothesized number of groups individuals form

#Use k-means clustering analysis to find groups that maximize variation between groups
Aty_test_data <- find.clusters(Bay_stru_data, max.n.clust = 20)      #k-means clustering
          #prompt in console: choose number of principal components to retain, choose 100
          #data transformed to principal components
          #principal components: uncorrelated variables that explain maximum possible variance
          #Bayesian Information Criterion: goodness-of-fit test statistic for k values, larger k values penalized for over-fitting
          #prompt in console: choose number of clusters, find 'bend' in data, choose 4
          #individuals are assigned to one of four groups (Levels 1,2,3,4), based on genetic association


Aty_test_data                  
     #Kstat: shows Bayesian incormation criteria for different cluster numbers
     #stat: BIC value for K=4
     #grp: the group each individual is assigned to 
names(Aty_test_data)
head(Aty_test_data$Kstat, 10)   #BIC values for first ten cluster values
head(Aty_test_data$stat)        #BIC valuse for chosen cluster value, which is 4 in this case
head(Aty_test_data$size)      #number of individuals assigned to each group
head(Aty_test_data$grp, 10)    #group association for first ten individuals


##########Run DAPC: discriminant analysis of principal components

#data transformed using PCA, then discriminant analysis on retained principal components
#test based on associations of individuals assigned to groups from find_clusters
?dapc
Aty_dapc_test <- dapc(Bay_stru_data, Aty_test_data$grp) 
#prompt in console: choose number of principal components to retain
            #choosing too few PCs will not provide enough info to test group association  
            #choosing many PCs will overfit data, or assign with over-confidence
            #check sensitivity to number of PCs retained

#prompt in console: choose number of discriminant functions (linear discriminants) to retain
          #discriminant functions: synthetic variables, show differences between clusters while minimizing variation within clusters
          #first two dimensions cover most of the information, choose two
          #Eigenvector: characteristic vector of data, direction is linear 
          #Eigenvalues: factor by which eigenvector is transformed linearily
          #Eigenvalues correspond to ratio of variance between groups over variance within groups for each discriminant function

Aty_dapc_test     #view statistical results of Discriminant Analysis of Principal Componenets


########## Visualize result of DAPC with scatterplots

# each point indicates placement of individual in statistical space relative to others

myCol <- c("red", "darkblue", "darkgreen", "purple", "orange", "skyblue")   #colors for scatterplot

scatter(Aty_dapc_test)
#scatter(Aty_dapc_test, posi.da="bottomright", bg="white", pch = 17:22)

scatter(Aty_dapc_test, posi.da="bottomright", bg = "white", pch = 17:22, 
        cstar = 0, col = myCol, scree.pca = TRUE, posi.pca="bottomleft")

scatter(Aty_dapc_test, scree.da = FALSE, bg = "white", pch = 20, cell = 0, cstar = 0, col = myCol, 
        solid = .4, cex = 2, clab = 0, leg = TRUE, txt.leg = paste ("Cluster",1:4))


####### Rerun DAPC, check sensitivity of using different number of Principal Components

#Rerun DAPC test with many PCs
Aty_dapc_test <- dapc(Bay_stru_data, Aty_test_data$grp) 
#prompt in console: choose number of principal components to retain
      #choose a high number of PCs, >40, document the number you picked: 
      #choose 2 for linear functions
      #generate scatterplots using scatter commands above
      #How has the distribution of clusters changed?

#Rerun DAPC test with few PCs
Aty_dapc_test <- dapc(Bay_stru_data, Aty_test_data$grp) 
#prompt in console: choose number of principal components to retain
        #choose a low number of PCs, <10, document the number you picked: 
        #choose 2 for linear functions
        #generate scatterplots using scatter commands above
        #How has the distribution of clusters changed?


#Export one of your scatterplots
#Save as image, jpeg
#name the file YOURNAME_##PCsYouPicked.jpg

#email me your jpg file (5 points)





