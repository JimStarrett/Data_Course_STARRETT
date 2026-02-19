#################### Principal Components Analysis ############################
#analysis of multivariate data, complex data with many dimensions reduced to few dimensions


###################################################################################### 
####### Investigation of population structure in the California Turret Spider ########

#goal is to assign individuals to populations (or clusters or groups) based on genetic variation
#optimizing variation between groups and minimizing variation within groups


############################################################
############ Download data files ###########################

#https://github.com/JimStarrett/Data_Course_STARRETT/DAPC_test.git
#download files:
#in terminal, clone data and scripts from github
#git clone https://github.com/JimStarrett/Data_Course_STARRETT/DAPC_test.git


############################################################
############   Use adegenet package   ######################
#adegenet package: investigation of population structure from DNA sequence data
#follow DAPC tutorial (Jombart and Collins, June 5, 2022)

# install adegenet, if you haven't already
install.packages('adegenet', dep=TRUE)

#install packages for mapping and graphing, if needed
install.packages('ggplot2')
install.packages('sf')
install.packages('rnaturalearth')
install.packages('rnaturalearthdata')
install.packages('maps')


###########################################################
########### Background on data being analyzed #############

#85 individuals of Atypoides riversi, 'California turret spider'
#individuals sampled from 53 localities near San Francisco, California

#SNP data, Single Nucleotide Polymorphism, in Structure format
#SNP data from 5421 loci generated with Restriction-site Associated DNA sequencing (3RAD)
#Structure format: alleles represented by numbers, rows=individuals, columns=DNA marker 


##############################################
####### load adegenet package ################
library(adegenet)
?adegenet       #access manual for adegenet package

#############################################
######### Read in SNP data ################## 

# read in data as structure file (.str)
#data is unlinked SNPs (1 SNP per locus)
#adjust file name, n.ind, n.loci; n.ind = number of individuals; n.loci = number of unlinked SNPs
Bay_stru_data <- read.structure("Atytest.stru", n.ind=85, n.loc=5421, onerowperind=FALSE, 
                                col.lab=1, col.pop=0, col.others=NULL, row.marknames=NULL, 
                                NA.char="-9", pop=NULL, ask=FALSE, quiet=FALSE)
Bay_stru_data       #data converted to genind format
#View(Bay_stru_data)



?find.clusters
####################################################################
########### Find hypothesized number of groups #####################

#Use k-means clustering analysis to find groups that maximize variation between groups
# data is transformed using PCA and reduced to principal components
#k-means clustering: runs successive k-means
#?find.clusters
Aty_test_data <- find.clusters(Bay_stru_data, max.n.clust = 20)      #k-means clustering
          #prompt in console: choose number of principal components to retain, choose 100
          #data transformed to principal components
          #principal components: uncorrelated variables that explain maximum possible variance
          #Bayesian Information Criterion: goodness-of-fit test statistic for k values, larger k values penalized for over-fitting
          #prompt in console: choose number of clusters, find 'bend' in data, choose 4
          #individuals are assigned to one of four groups (Levels 1,2,3,4), based on genetic association

Aty_test_data  #view statistical results of k-means clustering analysis                 
     #Kstat: shows Bayesian incormation criteria for different cluster numbers
     #stat: BIC value for K=4
     #grp: the group each individual is assigned to 

#Brief view of results 
#names(Aty_test_data)
#head(Aty_test_data$Kstat, 10)   #BIC values for first ten cluster values
#head(Aty_test_data$stat)        #BIC valuse for chosen cluster value, which is 4 in this case
#head(Aty_test_data$size)        #number of individuals assigned to each group
#head(Aty_test_data$grp, 10)     #group association for first ten individuals


######################################################################
############### Map the individuals with group assignments ###########

#load libraries for mapping
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(maps)

cluster_assign <- data.frame(Aty_test_data$grp)     #cluster asignment from k-means clustering

occurence_data <- read.delim("Aty_Bay_coord_v2.txt", header = FALSE, sep = "\t")  #read long lat coordinates 

clust_occur_data <- cbind(cluster_assign, occurence_data)  #bind cluster assignment with occurrence data 
#View(clust_occur_data)

colnames(clust_occur_data) <- c("Group", "longitude", "latitude")    #label columns
clust_occur_frame <- data.frame(clust_occur_data)               #make dataframe
View(clust_occur_frame)

# Load simple background map
world <- ne_countries(scale = "medium", returnclass = "sf")

#plot simple map showing individual localities with group assignment
ggplot() +                              
  borders("world", colour="gray85", fill="gray80") +         #map background aesthetic
  geom_point(data=clust_occur_frame, aes(x=longitude, y=latitude, color=Group, shape=Group), size=2) +  #localities aesthetic
  theme_minimal() +
  labs(title="Distribution of Atypoides populations", x="Longitude", y="Latitude") +
  coord_sf(xlim = c(-123.75, -121.1), ylim = c(36.75, 39.25), expand = FALSE)      #focus map to relevant region


#############################################################################
########## Run DAPC: discriminant analysis of principal components ##########

#Discriminant analysis:

#data transformed using PCA, then discriminant analysis on retained principal components
#test based on associations of individuals assigned to groups from k-means clustering (find_clusters)
?dapc     #access dapc manual if needed

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

Aty_dapc_test     #view statistical results of Discriminant Analysis of Principal Components


####################################################################
########## Visualize result of DAPC with scatterplots ##############

# each point indicates placement of individual in statistical space relative to others

myCol <- c("red", "darkblue", "darkgreen", "purple", "orange", "skyblue")   #colors for scatterplot

#scatter(Aty_dapc_test)
scatter(Aty_dapc_test, posi.da="bottomright", bg="white", pch = 17:22)

scatter(Aty_dapc_test, posi.da="bottomright", bg = "white", pch = 17:22, 
        cstar = 0, col = myCol, scree.pca = TRUE, posi.pca="bottomleft")

scatter(Aty_dapc_test, scree.da = FALSE, bg = "white", pch = 20, cell = 0, cstar = 0, col = myCol, 
        solid = .4, cex = 2, clab = 0, leg = TRUE, txt.leg = paste ("Cluster",1:5))


###################################################################################################
####### Rerun DAPC, check sensitivity of using different number of Principal Components ###########

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


######################################################
############ Thought Questions #######################

#1) How does the different levels or retained PCs impact your understanding of population structure?
# Think about the implications of uncertainty

#2) In what ways could you apply the results you obtained?
# Think about population connectivity, conservation in an urban landscape, past vs current genetic barriers.

#3) What other kinds of datasets could you apply these methods to?
# Think about topics related to genetics, ecology, medicine, conservation, etc.





#submit to canvas (5 points)





