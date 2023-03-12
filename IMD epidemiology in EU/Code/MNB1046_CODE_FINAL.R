###############################################################################
########## R script for generating the figures of the publication    ##########
##########          Regional Epidemiology Review (MNB1046)           ##########
##########    Author: Charles NUTTENS, charles.nuttens@pfizer.com    ##########
##########                  Version 0.4 - 21/07/2019                 ##########
##########           Confidential - For internal use only            ##########
###############################################################################

# WORKING DIRECTORY, LIBRARY, FUNCTIONS AND COLOR PALETTES
    ## Set working directory

setwd("C:/Users/nuttec/OneDrive - Pfizer/14 - R Epidemiology/IMD/Publication (MNB1046)")

# Loading -----------------------------------------------------------------


    ## Load library
library(tidyverse) #include ggplot2, dplyr and tidyr
library(beeswarm) #organized stripchart (alternative to stripchat())
library(RColorBrewer) #creation of the color palette for heatmap
library(gplots) #generation of heatmap
library(ggridges)

    ## Define functions
fun.stats <- function(x) {
    c(min = min(x, na.rm = TRUE), mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE), median = median(x, na.rm = TRUE), max = max(x, na.rm = TRUE))
}

    ## Define color palettes
color.sg <- c("black",       #all > black
              "dodgerblue2", #B > blue
              "green4",      #C > green     /!\ Alternative is "#E31A1C" (red)
              "#FF7F00",     #W > orange
              "#6A3D9A",     #W > purple
              "grey"         #Other > grey
              )

color.age <- c("#E31A1C",    #<01    > red
               "#FF7F00",    #01-04  > orange
               "gold1",      #05-14  > yellow
               "green4",     #15-24  > green
               "dodgerblue2",#25-49  > blue
               "#6A3D9A",    #50-64  > purple
               "black"       #65+    > black
               )

color.country <- rep(c("black",       
                       "dodgerblue2", 
                       "green4",      
                       "#FF7F00",    
                       "#6A3D9A",     
                       "grey",
                       "#E31A1C",
                       "gold1"
                       ),
                 5)

# DATABASE LOADING, CLEANING AND SUBSETTING
 ## Database loading 
df.full <- read.csv2("MNB1046_DATABASE_FINAL.csv", #load database
                     header = TRUE,
                     na.strings = "NA",
                     sep = ";",
                     dec = ".",
                     colClasses = c("Continent" = "character",
                                    "Region" = "character",
                                    "Country" = "character",
                                    "ISO3" = "character",
                                    "Year" = "integer",
                                    "Age" = "character",
                                    "Serogroup" = "character",
                                    "Outcome" = "character",
                                    "Variable" = "character",
                                    "By" = "character",
                                    "Unit" = "character",
                                    "Value" = "numeric",
                                    "Source" = "character",
                                    "Source.details" = "character",
                                    "Filename" = "character",
                                    "Hyperlink" = "character",
                                    "Confirmed" = "character"
                     ))

    ## Database Cleaning
df <- df.full %>%
        select(3:13, 17) %>% #remove unnecessary columns (i.e. continent, region, sources, files)
        filter(!Country %in% c("Malta", "Cyprus", "Luxembourg", "Iceland", "Russia", "Australia", "New Zealand", "South Korea", "Switzerland", "England"))%>% #remove non-EU countries and countries excluded from analysis (i.e. < 1 million population)
        filter(Outcome == "all") #select 'all' as outcome

rm(df.full) #remove 'df.full' from the working space

    ## Database subseting
df.ecdc <- filter(df, Source == "ECDC")
df.population <- filter(df, Variable == "population")

# GENERATING FIGURES

# Table 1 -----------------------------------------------------------------

    ## TABLE 1 - General statistics
        ### Selection of variables
year <- c(2008, 2017)

        ### Generation of incidence statistics at country level
df.inc <- filter(
    df.ecdc,
    Age == "all",
    Variable == "incidence",
    By  %in% c("serogroup", "overall"),
    Year %in% year
)

split.inc <- split(df.inc, df.inc$Year) #create a list of values split by year
spread.inc <- lapply(split.inc, function(i) {spread(i[,c(1,5,10)], key = Serogroup, value = Value)}) #spread the values by serogroup into a dataframe for each year of the list
spread.inc.other <- lapply(spread.inc, function(i) {i$Other <- apply(i[-1], 1, function(j) {j[1]-sum(j[2:5])}); return(i)}) #calculate the incidence for 'other" and add the values to the dataframes
stat.inc <- lapply(spread.inc.other, function (i) {apply(i[-1], 2, fun.stats)}) #calculate the statistics per year

        ### Generation of proportion statistics at country level
df.prop <- filter(
    df.ecdc,
    Age == "all",
    Variable == "distribution",
    Serogroup %in% c("B", "C", "W", "Y"),
    By  == "serogroup_in_overall",
    Year %in% year
)

split.prop <- split(df.prop, df.prop$Year) #create a list of values split by year
spread.prop <- lapply(split.prop, function(i) {spread(i[,c(1,5,10)], key = Serogroup, value = Value)}) #spread the values by serogroup into a dataframe for each year of the list
spread.prop.other <- lapply(spread.prop, function(i) {i$Other <- apply(i[-1], 1, function(j) {100-sum(j[1:4])}); return(i)}) #calculate the incidence for 'other" and add the values to the dataframes
stat.prop <- lapply(spread.prop.other, function (i) {apply(i[-1], 2, fun.stats)}) #calculate the statistics per year

        ### Generation of incidence and proportion data at EU/EEA level
#### Calulate the population in EU/EEA
df.pop <- filter(
    df.population,
    Age == "all",
    Variable == "population",
    By  == "overall",
    Year %in% year
)

split.pop <- split(df.pop, df.pop$Year) #create a list of values split by year
select.pop <- lapply(split.pop, function(i) {select(i, 1,10)}) #select the needed columns

#### Calulate the number of cases in EU/EEA
df.cases <- filter(
    df.ecdc,
    Age == "all",
    Variable == "cases",
    By  %in% c("serogroup", "overall"),
    Year %in% year
)

split.cases <- split(df.cases, df.cases$Year) #create a list of values split by year
spread.cases <- lapply(split.cases, function(i) {spread(i[,c(1,5,10)], key = Serogroup, value = Value)}) #spread the values by serogroup into a dataframe for each year of the list
spread.cases.other <- lapply(spread.cases, function(i) {i$Other <- apply(i[-1], 1, function(j) {j[1]-sum(j[2:5])}); return(i)}) #calculate the incidence for 'other" and add the values to the dataframes

#### Calulate the incidence by taking in consideration the NA value in the cases dataframe
select.pop.ordered <- mapply(function(i,j) {left_join(i[1], j, by = "Country")}, spread.cases.other, select.pop, SIMPLIFY = F) #align the population data with the order use in the case data frame
EU.inc <- mapply(function(i,j) {apply(i[-1],
                                      2,
                                      function(k) {(sum(k, na.rm = T)/sum(j[,2][!is.na(k)], na.rm = T))*100000}
                                      )},
                 spread.cases.other,
                 select.pop.ordered,
                 SIMPLIFY = F
                 )

#### Calulate the proportion by removing the countries without serogroups data
spread.cases.other.na <- lapply(spread.cases.other, function(i) {drop_na(i)})
sum.cases <- lapply(spread.cases.other.na, function(i) {apply(i[-1], 2, sum)})
EU.prop <- lapply(sum.cases, function(i) {(i/i[1])*100})

        ### Printing data
lapply(stat.inc, function(i) {t(round(i[c(1,4,5),],2))})
lapply(EU.inc, function(i) {round(i,8)})
lapply(EU.prop, function(i) {round(i,4)})

lapply(EU.inc, function(i) {round(i,4)})$`2017`/lapply(EU.inc, function(i) {round(i,4)})$`2008` #Fold Change

#       ###IQR + rIQR
# stat.inc$`2008`[7,]-stat.inc$`2008`[6,]
# (stat.inc$`2008`[7,]-stat.inc$`2008`[6,])/stat.inc$`2008`[5,]
# stat.inc$`2017`[7,]-stat.inc$`2017`[6,]
# (stat.inc$`2017`[7,]-stat.inc$`2017`[6,])/stat.inc$`2017`[5,]

# Supplemental Figure 1 ---------------------------------------------------

    ## SUPPLEMENTAL FIGURE 1 - Graphical representation of the dispersion of country incidence by serogroup and years
        ### selection of variables
year <- c(2008, 2017) #/!\ should be the same as TABLE 1 as the EU/EEA incidence value are comming from the calculation perfomed for TABLE 1

        ### Data sorting
#### Generate the dataframe
df.inc <- filter(df.ecdc,
                Age == "all",
                Variable == "incidence",
                By == "serogroup",
                Year %in% year
)

#### Generate the dataframe restrected to the ouliers
split.outliers <- split(df.inc, df.inc[,c("Year", "Serogroup")])
df.outliers <- lapply(split.outliers, 
                      function(i) {
                          q <- quantile(i$Value, na.rm = TRUE)
                          filter(i, Value > q[4]+2*(q[4]-q[2]))
})

        ### Setting the graphical areas
par(oma= c(1,3,2,1)) #outer margin area

layout(matrix(c(1,2), nrow = 2, byrow = TRUE), #top and bottom panel with a 1/3 ratio
       heights = c(1,3)
       )
#layout.show(n=2)

        ### Plotting the top panel
par(mar= c(0,2,1,1), bty = "7", cex = 0.75, lwd = 1) #setting the graphic margins

#### Generating the box plot
boxplot(Value~as.factor(Year)*as.factor(Serogroup),
        data=df.inc,
        ylim = c(1.5,4), #to modify according the the data
        border = "gray50",
        at = c(1,2,4,5,7,8,10,11),
        outline = FALSE,
        xaxt='n'
)

#### Generating the points layover
beeswarm(Value~as.factor(Year)*as.factor(Serogroup),data=df.inc,
           add = TRUE,
           pch = 20,
           col = "black",
           method = "swarm",
           spacing = 0.85, #0.85 for PDF generation
           at = c(1,2,4,5,7,8,10,11)
         )
#### Print the name of the ouliers
text(1, df.outliers$`2008.B`$Value, df.outliers$`2008.B`$ISO3, pos = 4, cex = 0.75)

abline(v=c(3,6,9), col = "black", lty = 4)
mtext(c("serogroup B", "serogroup C", "serogroup W", "serogroup Y"), side=3, line=0.7, at=c(1.5, 4.5, 7.5, 10.5))
legend("topright", col = c("black", "red"), pch = c(20,18), legend = c("individual country", "Europe"), bty = "n")

        ### Plotting the bottom panel
par(mar= c(3,2,1,1), bty = "u", cex = 0.75, lwd = 1) #setting the graphic margins

#### Generating the box plot
boxplot(Value~as.factor(Year)*as.factor(Serogroup),
        data=df.inc,
        ylim = c(0,1.2), #to modify according the the data
        border = "gray50",
        at = c(1,2,4,5,7,8,10,11),
        names = rep(c(year[1], year[2]), 4),
        outline = FALSE
)

#### Generating the points layover
beeswarm(Value~as.factor(Year)*as.factor(Serogroup),
           data=df.inc,
           add = TRUE,
           pch = 20,
           col = "black",
           method = "swarm",
           spacing = 0.85, #0.85 for PDF generation
           at = c(1,2,4,5,7,8,10,11)
         )

#### Add EU/EEA point
points(c(1,2,4,5,7,8,10,11), t(sapply(EU.inc, function(i) i[2:5])), col = "red", pch = 18, cex = 1) #use the data generated for TABLE 1

#### Adding the name of the outliers
text(1, df.outliers$`2008.B`$Value, df.outliers$`2008.B`$ISO3, pos = 4, cex = 0.75)
text(2, df.outliers$`2017.B`$Value, df.outliers$`2017.B`$ISO3, pos = 4, cex = 0.75)
text(4, df.outliers$`2008.C`$Value, df.outliers$`2008.C`$ISO3, pos = 4, cex = 0.75)
text(5, df.outliers$`2017.C`$Value, df.outliers$`2017.C`$ISO3, pos = 4, cex = 0.75)
text(7, df.outliers$`2008.W`$Value, df.outliers$`2008.W`$ISO3, pos = 4, cex = 0.75)
text(8, df.outliers$`2017.W`$Value, df.outliers$`2017.W`$ISO3, pos = 4, cex = 0.75)
text(10, df.outliers$`2008.Y`$Value, df.outliers$`2008.Y`$ISO3, pos = 4, cex = 0.75)
text(11, df.outliers$`2017.Y`$Value, df.outliers$`2017.Y`$ISO3, pos = 4, cex = 0.75)

abline(v=c(3,6,9), col = "black", lty = 4)

        ### Add legend and graphical items
mtext(side=2, line=1, cex=1, adj = 0.5, outer=T, "incidence (cases/100,000/year)")

        ### Formating the mulpitle frame
pdf("Supplemental Figure 1.pdf", width=8.3, height=6.5, paper = "a4", title = "Rplot")

dev.off()


# Figure 1 ----------------------------------------------------------------

    ## FIGURE 1 - Incidence by age group for serogroup B, C, W and Y
        ### Selection of variables
year <- 2017
sg <- "B" #to modify according to the serogroup

        ### Data selection
#### Main dataframe
df.inc <- filter(df.ecdc,
                Age != "50+",
                Variable == "incidence",
                By == "age*serogroup",
                Serogroup == sg,
                Year == year
                )

#### Calulate the number of cases by age group in EU/EEA using the total number of case and the proportion by age group for each individual country
df.cases.all <- filter(df.ecdc,
       Age != "50+",
       Variable == "cases",
       By == "serogroup",
       Serogroup == sg,
       Year == year
       )

df.age.prop <- filter(df.ecdc,
                Age != "50+",
                Variable == "distribution",
                By == "age_in_serogroup",
                Serogroup == sg,
                Year == year
                )

spread.df.age.prop <- spread(df.age.prop[,c(1,4,10)], key = "Age", value = Value)
df.age.cases <- left_join(spread.df.age.prop, df.cases.all[c(1,10)], by = "Country")
df.age.cases[2:8] <- (df.age.cases[,2:8]*df.age.cases[,9])/100 #calculate the number of cases per age group (this data is not provided by ECDC)
df.sum.cases <- apply(df.age.cases[2:8], 2, function(i) sum(i, na.rm = TRUE)) #calculate the sum of the number of cases

#### Calulate the population by age group in EU/EEA taking in consideration the NA value in the cases dataframe
df.pop <- filter(
    df.population,
    Age %in% c("<01", "01-04", "05-14", "15-24", "25-49", "50-64", "65+"),
    By  == "age",
    Year == year
)
spread.pop <- spread(df.pop[,c(1,4,10)], key = "Age", value = Value)
spread.pop.ordered <- left_join(df.age.cases[1], spread.pop, by = "Country") #align the population data with the order use in the case data frame
spread.pop.ordered[-1] <- apply(spread.pop.ordered[-1], 2, function(i) {i[is.na(df.age.cases[9])] <- NA; return(as.numeric(i))}) #remove the population values for countries without IMD cases data
df.sum.pop <- apply(spread.pop.ordered[-1], 2, function(i) sum(i, na.rm = TRUE)) #calculate the sum of the population

#### Calculate the incidence by age group
EU.inc <- (df.sum.cases/df.sum.pop)*100000

### Plot the name of outliers
split.outliers <- split(df.inc, df.inc[,c("Age")])
q.sg <- quantile(df.inc$Value, na.rm = T)
offset <- (q.sg[5]-q.sg[1])*0.003 #use a threshold for printing the name of outliers
df.outliers <- lapply(split.outliers,
                      function(i) {
                          q <- quantile(i$Value, na.rm = TRUE)
                          filter(i, Value > (q[4]+1.5*(q[4]-q[2]))+offset)
})

        ### Setting the graphical areas
par(oma= c(1,3,2,1)) #outer margin area

layout(matrix(c(1,2), nrow = 2, byrow = TRUE), #top and bottom panel with a 1/3 ratio
       heights = c(1,4)
       )
#layout.show(n=2)

#layout(matrix(c(1), nrow = 1, byrow = TRUE)) #Unique panel for serogroup B
#layout.show(n=1)

        ### Plotting the top panel
par(mar= c(0,2,1,1), bty = "7", cex = 0.65, lwd = 0.75) #setting the graphic margins
        
#### Generating the box plot
boxplot(Value~as.factor(Age),
        data=df.inc,
        outline = FALSE,
        ylim=c(1,2),
        border = "gray50",
        xaxt='n',
        yaxt='n'
        )

axis(2, lwd = 0.75)

#### Generating the points layover
beeswarm(Value~as.factor(Age),
         data=df.inc,
         add = TRUE,
         pch = 20,
         col = "black",
         method = "swarm",
         spacing = 0.60,
         cex = 0.65
         )

#### Print the name of the ouliers
text(1, df.outliers$`<01`$Value, df.outliers$`<01`$ISO3, pos = 4, cex = 0.50)
text(2, df.outliers$`01-04`$Value, df.outliers$`01-04`$ISO3, pos = 4, cex = 0.50)

### Add legend
legend("topright", col = c("black", "red"), pch = c(20,18), legend = c("individual country", "Europe"), bty = "n")


        ### Plotting the bottom panel
par(mar= c(2,2,1,1), bty = "u", cex = 0.65, lwd = 0.75) #setting the graphic margins
#par(mar= c(2,2,1,1), bty = "o", cex = 0.65, lwd = 0.75) #for serogoup B
        
#### Generating the box plot
boxplot(Value~as.factor(Age),
        data=df.inc,
        outline = FALSE,
        ylim=c(0,1),
        border = "gray50",
        xaxt='n',
        yaxt='n'
        )

axis(1, lwd = 0.75, at = seq(1,7) , labels = c("<1", "1-4", "5-14", "15-24", "25-49",  "50-64", "65+"))
#axis(2, lwd = 0.75, at = seq(0,16,2)) For B only
axis(2, lwd = 0.75)

#### Generating the points layover
beeswarm(Value~as.factor(Age),
         data=df.inc,
         add = TRUE,
         pch = 20,
         col = "black",
         method = "swarm",
         spacing = 0.60,
         cex = 0.65
         )

#### Print the EU/EAA points and the number of cases by age group
points(1:7, EU.inc, col = "red", pch = 18, cex = 1)
text(1:7, 0, pos = 1, offset = 0.25, labels = paste0("n=", round(df.sum.cases)), col = "red", cex = 0.60)
#text(1:7, 0, pos = 1, offset = 0.4, labels = paste0("n=", round(df.sum.cases)), col = "red", cex = 0.65) #B only
        ### Plot the name of outliers
text(1, df.outliers$`<01`$Value, df.outliers$`<01`$ISO3, pos = 4, cex = 0.50)
text(2, df.outliers$`01-04`$Value, df.outliers$`01-04`$ISO3, pos = 4, cex = 0.50)
text(3, df.outliers$`05-14`$Value, df.outliers$`05-14`$ISO3, pos = 4, cex = 0.50)
text(4, df.outliers$`15-24`$Value, df.outliers$`15-24`$ISO3, pos = 4, cex = 0.50)
text(5, df.outliers$`25-49`$Value, df.outliers$`25-49`$ISO3, pos = 4, cex = 0.50)
text(6, df.outliers$`50-64`$Value, df.outliers$`50-64`$ISO3, pos = 4, cex = 0.50)
text(7, df.outliers$`65+`$Value, df.outliers$`65+`$ISO3, pos = 4, cex = 0.50)

        ### Add legend and graphical items
mtext(side=2, line=1, cex=0.75, adj = 0.5, outer=T, "incidence (cases/100,000)")
mtext(side=3, line=0, cex=0.75, adj = 0.5, outer=T, paste0("serogroup ", sg))

        ### Formating the mulpitle frame
pdf(paste0("Figure 1_", sg, ".pdf"), width=4, height=4, paper = "a4", title = "Rplot")

dev.off()


# Figure 2 ----------------------------------------------------------------

    ## FIGURE 2 - Evolution of incidence and number of cases in EU/EEA per serogroups
### Setting the graphical areas
par(mfrow=c(2,2),
    oma = c(1, 1, 1, 1),
    mar=c(3, 4, 3, 4)+0.1,
    cex = 0.60,
    lwd = 0.8
    )

### Selection of variables
sg <- "Y"

### Generation of dataframe
df.inc <- filter(df.ecdc,
                Age == "all",
                Variable == "incidence",
                By == "serogroup",
                Serogroup == sg,
                Year >= 2008
)

df.cases <- filter(df.ecdc,
                 Age == "all",
                 Variable == "cases",
                 By == "serogroup",
                 Serogroup == sg,
                 Year >= 2008
)

df.pop <- filter(df.population,
                  Age == "all",
                  By  == "overall",
                  Year >= 2008
)

### Generating the box plot
boxplot(Value ~ as.factor(Year),
        data=df.inc,
        outline = FALSE, 
        ylim=c(0,max(df.inc$Value, na.rm=TRUE)),
        ylab="incidence (cases/100,000/year)",
        border = "gray50",
        lwd = 0.6,
        cex = 0.60,
        xaxt='n',
        yaxt='n'
        )
axis(1, lwd = 0.6, at = seq(1,10) , labels = seq(2008,2017))
axis(2, lwd = 0.6)


### Generating the points layover (circle)
pos = 0
for (i in unique(df.inc$Year)) {
    pos <- pos+1
    beeswarm(filter(df.inc, Year == i)$Value,
               add = TRUE,
               pch = 20,
               at = pos,
               col = "black",
               method = "swarm",
               spacing = 0.85,
               cex = 0.60,
               lwd = 0.8
    )
}

### Plotting the name of outliers
c <- 0
for (i in unique(df.inc$Year)) {
    yr <- i
    c <- c+1
    q <- quantile(filter(df.inc, Year == i)$Value, na.rm = TRUE)
    outliers <- filter(df.inc, Year == i, Value > q[4]+1.5*(q[4]-q[2]))
    if (nrow(outliers) == 0) {
    } else {
        text(c, outliers$Value, outliers$ISO3, pos = 4, cex = 0.60, offset = 0.5)
    }
}

### Add legend and graphical items
mtext(paste0("serogroup ", sg), side = 3, line = 1, cex = 0.80)
legend("topleft",
       col = c("black", "red"),
       pch = c(20, 18),
       legend = c("individual county", "Europe"),
       x.intersp = 2,
       bty = "n")


      ## Adding the total number of cases in EU
          ### Calculate the number of cases in EU/EEA
split.cases <- split(df.cases, df.cases$Year)
EU.cases <- sapply(split.cases, function(i) sum(i$Value, na.rm = TRUE))

        #### Calulate the population in EU/EEA taking in consideration the NA value in the cases dataframe
spread.df.pop <- spread(df.pop[,c(1,3,10)], key = "Year", value = Value)
spread.df.cases <- spread(df.cases[,c(1,3,10)], key = "Year", value = Value)
spread.pop.ordered <- left_join(spread.df.cases[1], spread.df.pop, by = "Country") #align the population data with the order use in the case data frame
spread.pop.ordered[is.na(spread.df.cases)] <- NA #remove the population values for countries without IMD cases data
df.pop.sum <- apply(spread.pop.ordered[-1], 2, function(j) {sum(j, na.rm = T)}) #calculate the sum of the population

EU.incidence <- (EU.cases/df.pop.sum)*100000

par(new = TRUE, mar=c(3.1, 5.95,3.1, 5.95))
plot(EU.cases, type = "b", pch=18, xaxt = "n", yaxt = "n", ylab = "", xlab = "", col = "red", lty = 1, bty="n", ylim = c(0, max(EU.cases)), cex = 0.60)

par(mar=c(3, 4, 3, 4)+0.1)
axis(side = 4, col.axis="red", lwd = 0.6, cex = 0.60)
mtext("number of cases in Europe (cases/year)", side = 4, line = 3, col = "red", cex = 0.60)
        

        ### Formating the mulpitle frame
pdf("Figure 2.pdf", width=11.7, height=8.3, paper = "a4r", title = "Rplot")

dev.off()


# Figure 4 ----------------------------------------------------------------

        ## FIGURE 4 - Heatmap
    ### Setting the graphical areas
par(mfrow=c(2,2),
    mar=c(2.1, 4.1, 2.1, 1),
    oma = c(4, 2, 2, 2),
    cex = 0.65,
    lwd = 0.55
)

    ### Selection of variables
sg <- "Y"
year.ini <- 2015
    
    ### Data selection
df.end <- filter(df.ecdc,
                     Variable == "incidence",
                     By %in% c("serogroup", "age*serogroup"),
                     Serogroup == sg,
                     Year == 2017
)

spread.df.end <- spread(df.end[c(1,4,10)], key = Age, value = Value)

df.ini <- filter(df.ecdc,
                     Variable == "incidence",
                     By %in% c("serogroup", "age*serogroup"),
                     Serogroup == sg,
                     Year == year.ini
)

spread.df.ini <- spread(df.ini[c(1,4,10)], key = Age, value = Value)


evo <- spread.df.end[-1]/spread.df.ini[-1]
rownames(evo) <- spread.df.end[,1]
evo[evo == "NaN"] <- 1
evo.adj <- evo

evo.adj[is.na(evo) | evo >= 1] <- evo[is.na(evo) | evo >= 1]-1
evo.adj[is.na(evo) | evo < 1] <- -1/evo[is.na(evo) | evo < 1] +1
evo.adj[evo.adj == "Inf"] <- 22
evo.adj[evo.adj == "-Inf"] <- -22
evo.adj <- evo.adj[,c(8, 1:7)]
evo.adj.order <- evo.adj[order(evo.adj[,1],decreasing=T),]



#### Calulate the number of cases by age group in EU/EEA using the total number of case and the proportion by age group for each individual country
df.cases.all <- filter(df.ecdc,
       Age == "all",
       Variable == "cases",
       By == "serogroup",
       Serogroup == sg,
       Year %in% c(2015, 2017)
       )

df.age.prop <- filter(df.ecdc,
                Age %in% c("<01", "01-04", "05-14", "15-24", "25-49", "50-64", "65+"),
                Variable == "distribution",
                By == "age_in_serogroup",
                Serogroup == sg,
                Year %in% c(2015, 2017)
                )

spread.df.cases.all <- spread(df.cases.all[,c(1,3,10)], key = "Year", value = Value)
spread.df.cases.all.sub <- round(spread.df.cases.all[,3] - spread.df.cases.all[,2], 0)
    
split.df.age.prop <- split(df.age.prop, df.age.prop$Age)
spread.df.age.prop <- lapply(split.df.age.prop, function(i) {spread(i[,c(1,3,10)], key = "Year", value = Value)})
df.age.cases <- lapply(spread.df.age.prop, function(i){i[-1] <- (i[-1]*spread.df.cases.all[-1])/100; return(i)})
df.age.cases.sub <- as.data.frame(sapply(df.age.cases, function(i) {round(i[,3] - i[,2], 0)}))
row.names(df.age.cases.sub) <- df.age.cases$`<01`$Country
df.age.cases.sub$all <- spread.df.cases.all.sub
df.age.cases.sub <- df.age.cases.sub[,c(8, 1:7)]
df.age.cases.sub.order <- df.age.cases.sub[order(evo.adj[,1],decreasing=T),]

#df.age.cases.sub.order.over <- df.age.cases.sub.order
#df.age.cases.sub.order.over[evo.adj.order<=3 & evo.adj.order>= -3] <- NA 

    ### Heatmap setting
#Breaks=c(seq(-5,-1,0.5), seq(1,5,0.5))
breaks=seq(-3,3,0.35)
colors=rev(brewer.pal(11,"RdBu"))
colors=colorRampPalette(colors)(17)

par(cex.main=0.8)
heatmap.2(as.matrix(evo.adj.order),
          dendrogram="none",
          Rowv=FALSE,
          Colv=FALSE, 
          col = colors,
          breaks = breaks,
          scale="none",
          key=TRUE,
          density.info="none", 
          trace="none",
          cexRow=0.8,
          cexCol=0.8,
          symm=F,
          symkey=T,
          symbreaks=T,
          sepcolor="white",
          sepwidth=c(0.05,0.05),
          colsep=1:ncol(evo),
          rowsep=1:nrow(evo),
          na.color = "grey",
          cellnote = df.age.cases.sub.order,
          notecol = "white",
          notecex = 0.5,
          #keysize = 0.5,
          margins=c(3,8),
          srtCol = 45,
          adjCol = c(1,0),
          adjRow = c(0,0.5),
          offsetRow = 0,
          offsetCol = 0,
          lmat = rbind(4:3,2:1),
          lwid = c(1.5,4),
          lhei = c(1.4,4),
          main=paste0("Serogroup ", sg),
          labCol = c("all", "<1", "1-4", "5-14", "15-24", "25-49",  "50-64", "65+")
)

### Formating the mulpitle frame
pdf(paste0("Figure 4 (serogroup ", sg, ").pdf"), width=5, height=6, paper = "a4", title = "Rplot")

dev.off()



heatmap.2(as.matrix(evo.adj.order),
          dendrogram="none",
          Rowv=FALSE,
          Colv=FALSE, 
          col = colors,
          breaks = breaks,
          scale="none",
          key=TRUE,
          density.info="none",
          trace="none",
          cexRow=0.8,
          cexCol=0.8,
          symm=F,
          symkey=T,
          lmat = rbind(4:3,2:1),
          lwid = c(4,4),
          lhei = c(1.5,4)
          )
### Formating the mulpitle frame
pdf(paste0("Figure 4_Key.pdf"), width=4, height=5, paper = "a4", title = "Rplot")

dev.off()





# Figure 3 ----------------------------------------------------------------

    ## FIGURE 3 - Evolution of incidence per age groups
        ### Generation of dataframe for EU/EEA

par(mfrow=c(3,2),
    oma = c(1, 1, 1, 1),
    mar=c(3, 4, 2, 2)+0.1,
    cex = 0.60,
    lwd = 0.8
    )

sg = "Y"

#### Calulate the number of cases by age group in EU/EEA using the total number of case and the proportion by age group for each individual country
df.cases.all <- filter(df.ecdc,
       Age == "all",
       Variable == "cases",
       By == "serogroup",
       Serogroup == sg,
       Year >= 2008
       )

df.age.prop <- filter(df.ecdc,
                Age %in% c("<01", "01-04", "05-14", "15-24", "25-49", "50-64", "65+"),
                Variable == "distribution",
                By == "age_in_serogroup",
                Serogroup == sg,
                Year >= 2008
                )

spread.df.cases.all <- spread(df.cases.all[,c(1,3,10)], key = "Year", value = Value)

split.df.age.prop <- split(df.age.prop, df.age.prop$Age)
spread.df.age.prop <- lapply(split.df.age.prop, function(i) {spread(i[,c(1,3,10)], key = "Year", value = Value)})
df.age.cases <- lapply(spread.df.age.prop, function(i){i[-1] <- (i[-1]*spread.df.cases.all[-1])/100; return(i)})
df.age.sum <- t(sapply(df.age.cases, function(i) {apply(i[-1], 2, function(j) {sum(j, na.rm = T)})}))


#### Calulate the population by age group in EU/EEA taking in consideration the NA value in the cases dataframe
df.pop <- filter(
    df.population,
    Age %in% c("<01", "01-04", "05-14", "15-24", "25-49", "50-64", "65+"),
    By  == "age",
    Year >= 2008
)

split.df.pop <- split(df.pop, df.pop$Age)
spread.df.pop <- lapply(split.df.pop, function(i) {spread(i[,c(1,3,10)], key = "Year", value = Value)})
spread.pop.ordered <- lapply(spread.df.pop,function(i) {left_join(df.age.cases$`<01`[1], i, by = "Country")}) #align the population data with the order use in the case data frame
spread.pop.ordered.NA <- lapply(spread.pop.ordered, function(i) {i[is.na(df.age.cases$`<01`)] <- NA; return(i)})#remove the population values for countries without IMD cases data
df.pop.sum <- t(sapply(spread.pop.ordered.NA, function(i) {apply(i[-1], 2, function(j) {sum(j, na.rm = T)})})) #calculate the sum of the population

#### Calculate the incidence by age group
EU.inc <- (df.age.sum/df.pop.sum)*100000

#### Calculate the evolution
EU.evo <- EU.inc/EU.inc[,1]
EU.evo[EU.evo<1]<- 1-1/EU.evo[EU.evo<1]
EU.evo[EU.evo>=1]<- EU.evo[EU.evo>=1]-1

plot(colnames(EU.evo),
     EU.evo[1,],
     type = "n",
     ylab = "fold change compared to 2008",
     ylim = c(ifelse(min(EU.evo)<=-1,floor(min(EU.evo)),-1), ifelse(max(EU.evo)>=1,ceiling(max(EU.evo)),1)),
     yaxt='n',
     xaxt='n',
     xlab = ""
     )

seq <- seq(floor(min(EU.evo)),ceiling(max(EU.evo)))

axis(side = 2,
     at=seq,
     labels=c(paste0("1/", -(seq[seq<0]-1)), 1, (seq[seq>0])+1),
     lwd = 0.8
     )

axis(side = 1, at = seq(2008, 2017), labels = seq(2008, 2017), lwd = 0.8)

abline(h=0, col= "grey", lty=4)

for (i in seq(nrow(EU.evo))) {
    lines(colnames(EU.evo), EU.evo[i,], type = "p", col = color.age[i], pch = 16, lwd = 0.65, cex = 0.65)
    lines(lowess(colnames(EU.evo), EU.evo[i,], f = 0.75), col = color.age[i], type = "b", lwd = 0.65, cex = 0.65)
}

legend(ifelse(sg %in% c("B", "C"),"bottomleft", "topleft"), legend = c("<1", "1-4", "5-14", "15-24", "25-49",  "50-64", "65+"), col=color.age, pch=16, bty = "n", x.intersp = 0.5, y.intersp = 0.80, cex = 0.90)
mtext(paste0("serogroup ", sg), side = 3, line = 0.5, cex = 0.75)

### Formating the mulpitle frame
pdf("Figure 3.pdf", width=8.27, height=10, paper = "a4", title = "Rplot")

dev.off()


# Figure 5 ----------------------------------------------------------------

## FIGURE 5 - Evolution of incidence per age groups for all serogroups
### Generation of dataframe for EU/EEA

par(mfrow=c(1,1),
    oma = c(1, 1, 1, 1),
    mar=c(3, 4, 2, 2)+0.1,
    cex = 0.65,
    lwd = 0.80
)

#### Calulate the number of cases by age group in EU/EEA using the total number of case and the proportion by age group for each individual country
df.cases.all <- filter(df.ecdc,
                       Age == "all",
                       Variable == "cases",
                       By == "overall",
                       Year >= 2008
)

df.age.prop <- filter(df.ecdc,
                      Age %in% c("<01", "01-04", "05-14", "15-24", "25-49", "50-64", "65+"),
                      Variable == "distribution",
                      By == "age_in_overall",
                      Year >= 2008
)

spread.df.cases.all <- spread(df.cases.all[,c(1,3,10)], key = "Year", value = Value)

split.df.age.prop <- split(df.age.prop, df.age.prop$Age)
spread.df.age.prop <- lapply(split.df.age.prop, function(i) {spread(i[,c(1,3,10)], key = "Year", value = Value)})
df.age.cases <- lapply(spread.df.age.prop, function(i){i[-1] <- (i[-1]*spread.df.cases.all[-1])/100; return(i)})
df.age.sum <- t(sapply(df.age.cases, function(i) {apply(i[-1], 2, function(j) {sum(j, na.rm = T)})}))


#### Calulate the population by age group in EU/EEA taking in consideration the NA value in the cases dataframe
df.pop <- filter(
  df.population,
  Age %in% c("<01", "01-04", "05-14", "15-24", "25-49", "50-64", "65+"),
  By  == "age",
  Year >= 2008
)

split.df.pop <- split(df.pop, df.pop$Age)
spread.df.pop <- lapply(split.df.pop, function(i) {spread(i[,c(1,3,10)], key = "Year", value = Value)})
spread.pop.ordered <- lapply(spread.df.pop,function(i) {left_join(df.age.cases$`<01`[1], i, by = "Country")}) #align the population data with the order use in the case data frame
spread.pop.ordered.NA <- lapply(spread.pop.ordered, function(i) {i[is.na(df.age.cases$`<01`)] <- NA; return(i)})#remove the population values for countries without IMD cases data
df.pop.sum <- t(sapply(spread.pop.ordered.NA, function(i) {apply(i[-1], 2, function(j) {sum(j, na.rm = T)})})) #calculate the sum of the population

#### Calculate the incidence by age group
EU.inc <- (df.age.sum/df.pop.sum)*100000

#### Calculate the evolution
EU.evo <- EU.inc/EU.inc[,1]
EU.evo[EU.evo<1]<- 1-1/EU.evo[EU.evo<1]
EU.evo[EU.evo>=1]<- EU.evo[EU.evo>=1]-1

plot(colnames(EU.evo),
     EU.evo[1,],
     type = "n",
     ylab = "fold change compared to 2008",
     ylim = c(ifelse(min(EU.evo)<=-1,floor(min(EU.evo)),-1), ifelse(max(EU.evo)>=1,ceiling(max(EU.evo)),1)),
     yaxt='n',
     xaxt='n',
     xlab = ""
)

seq <- seq(floor(min(EU.evo)),ceiling(max(EU.evo)))

axis(side = 2,
     at=seq,
     labels=c(paste0("1/", -(seq[seq<0]-1)), 1, (seq[seq>0])+1),
     lwd = 0.8
)

axis(side = 1, at = seq(2008, 2017), labels = seq(2008, 2017), lwd = 0.8)

abline(h=0, col= "grey", lty=4)

for (i in seq(nrow(EU.evo))) {
  lines(colnames(EU.evo), EU.evo[i,], type = "p", col = color.age[i], pch = 16, lwd = 0.65, cex = 0.65)
  lines(lowess(colnames(EU.evo), EU.evo[i,], f = 0.75), col = color.age[i], type = "b", lwd = 0.65, cex = 0.65)
}

legend("bottomleft", legend = c("<1", "1-4", "5-14", "15-24", "25-49",  "50-64", "65+"), col=color.age, pch=16, bty = "n", x.intersp = 0.5, y.intersp = 0.80, cex = 0.90)
mtext(paste0("IMD (ie, all serogroups and non-groupable)"), side = 3, line = 0.5, cex = 0.85)

### Formating the mulpitle frame
pdf("Figure 5.pdf", width=5.36, height=4.55, paper = "a4", title = "Rplot")

dev.off()


### Alternative plotting (Incidence)

### Setting the graphical areas
par(oma= c(1,3,2,1)) #outer margin area

layout(matrix(c(1,2), nrow = 2, byrow = TRUE), #top and bottom panel with a 1/3 ratio
       heights = c(1,3)
)
#layout.show(n=2)

### Plotting the top panel
par(mar= c(0,2,1,1), bty = "7", cex = 0.75, lwd = 0.75) #setting the graphic margins

plot(colnames(EU.inc),
     EU.inc[1,],
     type = "n",
     ylab = "",
     ylim = c(2.7,20),
     xlab = "",
     xaxt='n'
)

for (i in seq(nrow(EU.inc))) {
  lines(colnames(EU.inc), EU.inc[i,], type = "p", col = color.age[i], pch = 16)
  lines(lowess(colnames(EU.inc), EU.inc[i,], f = 0.75), col = color.age[i], type = "b")
}

legend("topright", legend = rownames(EU.inc), col=color.age, pch=16, bty = "n", x.intersp = 0.5, y.intersp = 0.80, cex = 0.90)

### Plotting the bottom panel
par(mar= c(3,2,1,1), bty = "u", cex = 0.75, lwd = 0.75) #setting the graphic margins

plot(colnames(EU.inc),
     EU.inc[1,],
     type = "n",
     ylab = "",
     ylim = c(0,1.5),
     xlab = ""
)

for (i in seq(nrow(EU.inc))) {
  lines(colnames(EU.inc), EU.inc[i,], type = "p", col = color.age[i], pch = 16)
  lines(lowess(colnames(EU.inc), EU.inc[i,], f = 0.75), col = color.age[i], type = "b")
}

### Add legend and graphical items
mtext(side=2, line=1, cex=1, adj = 0.5, outer=T, "Incidence (cases/100,000/year)")
mtext(paste0("All cases"), side = 3, line = 0.5, cex = 0.85, outer=T)

### Formating the mulpitle frame
pdf("Figure 3_all cases_incidence.pdf", width=8.27, height=7, paper = "a4", title = "Rplot")

dev.off()



# Supplemental Figure 2 ---------------------------------------------------

    ## SUPPLEMENTAL FIGURE 2 - Evolution of serogroup distribution in EU/EEA (2008-2017)
        ### Setting the graphical areas    
par(mfrow=c(1,1),
    oma = c(0, 0, 0, 0),
    mar=c(3, 4, 2, 1)+0.1,
    cex = 0.85,
    lwd = 1
)
        ### Data selection
df.cases <- filter(df.ecdc,
                  Age == "all",
                  Variable == "cases",
                  By  %in% c("serogroup", "overall"),
                  Year >= 2008
                  )

split.df.cases <- split(df.cases, df.cases$Year) #split of data by year
spread.df.cases <- lapply(split.df.cases, function(i) {spread(i[c(1, 5,10)], key = Serogroup, value = Value)}) #spread of data by serogroup for each year
spread.df.cases.NA <- lapply(spread.df.cases, function(i) {drop_na(i)}) #removing countries having NA value
sum.cases <- sapply(spread.df.cases.NA, function(i) {apply(i[-1], 2, function(j) {sum(j, na.rm = T)})}) #sum of cases
prop.sg <- apply(sum.cases, 2, function(i) {(i[-1]/i[1])*100}) #calculation of proportion of serogroups

        ### Plotting
plot(colnames(prop.sg), #plot the empty graph
     prop.sg[1,],
     type = "n",
     xlab = NA,
     ylab = "serogroup distribution (%)",
     ylim = c(0,100)
     )

abline(h=c(10,20,40,60,80), col= "grey", lty=4) #add reference line at 10%, 20%, 40%, 60% and 80%

####add the point and fitting line for each serogorup
for (i in seq(nrow(prop.sg))) { 
    lines(colnames(prop.sg), prop.sg[i,], type = "p", col = color.sg[-1][i], pch = 16)
    lines(lowess(colnames(prop.sg), prop.sg[i,], f = 0.75), col = color.sg[-1][i], type = "b")
}

legend("topright", row.names(prop.sg), col=color.sg[-1], pch=16, bty = "n", x.intersp = 0.5, y.intersp = 1, cex = 0.85) #add the legend

        ### Formating the mulpitle frame
pdf(paste0("Supplemental Figure 2.pdf"), width=8.27, height=4, paper = "a4", title = "Rplot")

dev.off()


# Supplemental Figure 3 ---------------------------------------------------

    ## SUPPLEMENTAL FIGURE 3 - Evolution of serogroup distribution in EU/EEA (2008-2017)
        ### Setting the graphical areas    
par(mfrow=c(1,1),
    oma = c(0, 0, 0, 0),
    mar=c(3, 4, 2, 1)+0.1,
    cex = 0.85,
    lwd = 1
)
        ### Data selection
#### Calulate the population in EU/EEA
df.pop <- filter(
    df.population,
    Age == "all",
    Variable == "population",
    By  == "overall",
    Year >= 2008
)

split.pop <- split(df.pop, df.pop$Year) #create a list of values split by year
select.pop <- lapply(split.pop, function(i) {select(i, 1,10)}) #select the needed columns

#### Calulate the number of cases in EU/EEA
df.cases <- filter(
    df.ecdc,
    Age == "all",
    Variable == "cases",
    By  %in% c("serogroup", "overall"),
    Year >= 2008
)

split.cases <- split(df.cases, df.cases$Year) #create a list of values split by year
spread.cases <- lapply(split.cases, function(i) {spread(i[,c(1,5,10)], key = Serogroup, value = Value)}) #spread the values by serogroup into a dataframe for each year of the list
spread.cases.other <- lapply(spread.cases, function(i) {i$Other <- apply(i[-1], 1, function(j) {j[1]-sum(j[2:5])}); return(i)}) #calculate the incidence for 'other" and add the values to the dataframes

#### Calulate the incidence by taking in consideration the NA value in the cases dataframe
select.pop.ordered <- mapply(function(i,j) {left_join(i[1], j, by = "Country")}, spread.cases.other, select.pop, SIMPLIFY = F) #align the population data with the order use in the case data frame
EU.inc <- mapply(function(i,j) {apply(i[-1],
                                      2,
                                      function(k) {(sum(k, na.rm = T)/sum(j[,2][!is.na(k)], na.rm = T))*100000}
                                      )},
                 spread.cases.other,
                 select.pop.ordered,
                 SIMPLIFY = T
                 )

#### Calculate the evolution
EU.evo <- EU.inc/EU.inc[,1]
EU.evo[EU.evo<1]<- 1-1/EU.evo[EU.evo<1]
EU.evo[EU.evo>=1]<- EU.evo[EU.evo>=1]-1

EU.evo <- EU.evo[-6,]

plot(colnames(EU.evo),
     EU.evo[1,],
     type = "n",
     ylab = "fold change compared to 2008",
     ylim = c(ifelse(min(EU.evo)<=-1,floor(min(EU.evo)),-1), ifelse(max(EU.evo)>=1,ceiling(max(EU.evo)),1)),
     yaxt = 'n',
     xlab = ""
     )

abline(h=0, col= "grey", lty=4)

seq <- seq(floor(min(EU.evo)),ceiling(max(EU.evo)))

axis(side = 2,
     at=seq,
     labels=c(paste0("1/", -(seq[seq<0]-1)), 1, (seq[seq>0])+1)
     )

for (i in seq(nrow(EU.evo))) {
    lines(colnames(EU.evo), EU.evo[i,], type = "p", col = color.sg[i], pch = 16)
    lines(lowess(colnames(EU.evo), EU.evo[i,], f = 0.75), col = color.sg[i], type = "b")
}

legend("topleft", row.names(EU.evo), col=color.sg, pch=16, bty = "n", x.intersp = 0.5, y.intersp = 1, cex = 0.85) #add the legend

        ### Formating the mulpitle frame
pdf(paste0("Supplemental Figure 3.pdf"), width=8.27, height=4, paper = "a4", title = "Rplot")

dev.off()


# Supplemental Figure 4 ---------------------------------------------------

        ## SUPPLEMENTAL FIGURE 4 - Ridgeline plot of serogroup W
    ### Setting the graphical areas
par(mfrow=c(1,1),
    mar=c(2.1, 4.1, 2.1, 1),
    oma = c(4, 2, 2, 2),
    cex = 2,
    lwd = 2
)

    ### Data sorting
df.inc <- filter(df.ecdc,
                Age == "all",
                Variable == "incidence",
                By == "serogroup",
                Serogroup == "W",
                Year >= 2008
)

    ### plotting
ggplot(df.inc, aes(Year, Country, height = Value)) +
    geom_ridgeline(scale = 6,alpha = 0.5, size = 0.35) +
    theme(panel.background = element_rect(fill = "white"),
          axis.title=element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(colour = "grey"),
          element_line(size = 0.25),
          axis.text = element_text(colour = "black", vjust =0),
          axis.ticks.y = element_line(colour = "grey")
          ) +
    scale_x_continuous(breaks= c(2008,2010,2012,2014,2016)) +
    scale_y_discrete(limits = c("Bulgaria",
                                "Croatia",
                                "Romania",
                                "Estonia",
                                "Latvia",
                                "Lithuania",
                                "Slovenia",
                                "Slovakia",
                                "Portugal",
                                "Hungary",
                                "Austria",
                                "Finland",
                                "Greece",
                                "Italy",
                                "Germany",
                                "Czech Republic",
                                "Norway",
                                "Poland",
                                "Spain",
                                "Belgium",
                                "France",
                                "Sweden",
                                "Denmark",
                                "Ireland",
                                "Netherlands",
                                "United Kingdom"
                                )
                     )

pdf("Figure 4 (serogroup W).pdf", width=7, height=9, paper = "a4", title = "Rplot")

dev.off()


# Supplemental Figure 6 ---------------------------------------------------

    ## SUPPLEMENTAL FIGURE 6 : Evolution of number of cases per country

        ### Selection of variables
sg <- "B"
ylabel <- c("AUT", "HRV", "EST", "DEU", "IRL", "LTU", "POL", "SVK", "SWE")

### Generation of dataframe
df.cases <- filter(df.ecdc,
                  Variable == "cases",
                  By == "serogroup",
                  Age == "all",
                  Year >=2008,
                  Serogroup == sg
)

split.df.cases <- split(df.cases, df.cases$Country)

### Generation of PDF

par(mfrow=c(9,3),
    oma = c(1, 2, 1, 1),
    mar=c(2.5, 3.8, 2.5, 1),
    cex = 0.4
    )

mtext(paste("serogroup", sg), side = 3, line = 0.5, cex = 0.85)

lapply(split.df.cases, function(i) 
    if (length(i$Value) == 0 | max(i$Value, na.rm=TRUE) == "-Inf") {
        plot(x=100, y=100, type="n", yaxt='n', main = unique(i$Country), xlab = NA, ylab = if (unique(i$ISO3) %in% ylabel) "number of cases" else "", xlim=c(2008,2017), ylim=c(0,110))
        text(2012.5, 50, "data not available")
    } else {
        plot(i$Year, i$Value, type="b", main = unique(i$Country), xlab = NA, ylab = if (unique(i$ISO3) %in% ylabel) "number of cases" else "", xlim=c(2008,2017), ylim=c(0,max(i$Value, na.rm=TRUE)), pch = 16)
    }
)

pdf(paste0("Supplemental Figure 6 (serogroup ", sg, ").pdf"), width=7.44, height=10.52, paper = "a4")

dev.off()




# Additional Figures ------------------------------------------------------


    ## ADDITIONAL FIGURE 1 - Evolution of proportion of cases by age group for each serogroup
        ### Selection of variables
sg <- "W"

        ### Generation of dataframe for individual countries
df.temp.evo <- filter(df.ecdc,
                Age != "50+",
                Variable == "distribution",
                By == "age_in_serogroup",
                Serogroup == sg,
                Year >= 2008
)

df.temp.evo.l <- split(df.temp.evo, df.temp.evo$Age)
df.temp.evo.l <- lapply(df.temp.evo.l, function(i) spread(i[c(1,3,10)], key = Year, value = Value))
df.temp.stats <- t(sapply(df.temp.evo.l, function(i) apply(i[-1], 2, fun.stats)[2,])) #use mean and not EU/EEA data


plot(colnames(df.temp.stats),
     df.temp.stats[1,],
     type = "n",
     xlab = NA,
     ylab = "%",
     main = "Evolution of distribution of age group",
     ylim = c(0,50)
     )

for (i in seq(nrow(df.temp.stats))) {
    lines(colnames(df.temp.stats), df.temp.stats[i,], type = "p", col = color.age[i], pch = 16)
    lines(lowess(colnames(df.temp.stats), df.temp.stats[i,], f = 0.75), col = color.age[i], lwd=1, type = "b")
}

legend("topright", row.names(df.temp.stats), col=color.age[1:7], pch=16, bty = "n", cex=0.9, x.intersp = 0.5, y.intersp = 1, pt.cex = 1)




   
    ## EVOLUTION OF POPULATION IN EU/EEA
        ###Selection of data

par(mfrow=c(1,1),
    oma = c(0, 0, 0, 0),
    mar=c(3, 4, 2, 1)+0.1,
    cex = 0.85,
    lwd = 0.85
)

df.pop <- filter(df.population,
                 Age %in% c("<01", "01-04", "05-14", "15-24", "25-49", "50-64", "65+"),
                 Variable == "population",
                 By  == "age",
                 Year >= 1999
                 )

split.df.pop <- split(df.pop, df.pop$Age)
spread.df.pop <- lapply(split.df.pop,function(i) {spread(i[c(2, 3,10)], key = Year, value = Value)})
sum.df.pop <- sapply(spread.df.pop, function (i) {colSums(i[-1], na.rm = T)})

matplot(row.names(sum.df.pop),
        sum.df.pop/1000000,
        pch = 16,
        col = color.age,
        type = "b",
        lty = 1,
        ylab = "Population (in million)",
        xlab = ""
        )

legend(x=1999, y = 170, legend = colnames(sum.df.pop), col = color.age, pch = 16, bty = "n")

### Formating the mulpitle frame
pdf("Population_EU.pdf", width=8.27, height=7, paper = "a4", title = "Rplot")

dev.off()


#by countries
df.pop.age <- spread.df.pop$`65+`
row.names(df.pop.age) <- df.pop.age[,1]
df.pop.age.t <- t(df.pop.age[,-1])

color.country <- c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#000000')

matplot(row.names(df.pop.age.t),
        df.pop.age.t/1000000,
        pch = 16,
        type = "b",
        lty = 1,
        cex = 0.5,
        ylab = "Population (in million)",
        xlab = "",
        ylim = c(0,4),
        col = color.country
        )

mtext(paste0("65+"), side = 3, line = 1, cex = 0.85, outer=F)

text(2017, df.pop.age.t[dim(df.pop.age.t)[1],]/1000000, colnames(df.pop.age.t), cex = 0.7, pos = 4)

### Formating the mulpitle frame
pdf("Population_age_EU.pdf", paper = "a4", title = "Rplot")

dev.off()


    ## FIGURE 3_ALL CASES_COUNTRY - Evolution of country incidence per age groups for all serogroups
        ### Generation of dataframe for EU/EEA

par(mfrow=c(1,1),
    oma = c(1, 1, 1, 1),
    mar=c(3, 4, 2, 2)+0.1,
    cex = 0.75,
    lwd = 0.75
)

country <- "Spain"
pdf(paste0("Figure 3_all cases_", country,"_2008.pdf"), width=8.27, height=7, paper = "a4", title = "Rplot")

#### Calulate the number of cases by age group in EU/EEA using the total number of case and the proportion by age group for each individual country
df.inc <- filter(df.ecdc,
                       Country == country,
                       Age %in% c("<01", "01-04", "05-14", "15-24", "25-49", "50-64", "65+"),
                       Variable == "incidence",
                       By == "age",
                       Serogroup == "all",
                       Year >= 2008
)


spread.df.inc <- spread(df.inc[,c(3,4,10)], key = "Year", value = Value)

#### Calculate the evolution
EU.evo <- spread.df.inc[-1]/spread.df.inc[,2]
EU.evo[EU.evo<1]<- 1-1/EU.evo[EU.evo<1]
EU.evo[EU.evo>=1]<- EU.evo[EU.evo>=1]-1
row.names(EU.evo) <- spread.df.inc$Age

plot(colnames(EU.evo),
     EU.evo[1,],
     type = "n",
     ylab = "Fold change compared to 2008",
     ylim = c(ifelse(min(EU.evo)<=-1,floor(min(EU.evo)),-1), ifelse(max(EU.evo)>=1,ceiling(max(EU.evo)),1)),
     yaxt='n',
     xlab = ""
)

seq <- seq(floor(min(EU.evo)),ceiling(max(EU.evo)))

axis(side = 2,
     at=seq,
     labels=c(paste0("1/", -(seq[seq<0]-1)), 1, (seq[seq>0])+1)
)

abline(h=0, col= "grey", lty=4)

for (i in seq(nrow(EU.evo))) {
    lines(colnames(EU.evo), EU.evo[i,], type = "p", col = color.age[i], pch = 16)
    lines(lowess(colnames(EU.evo), EU.evo[i,], f = 0.75), col = color.age[i], type = "b")
}

legend("bottomleft", legend = rownames(EU.evo), col=color.age, pch=16, bty = "n", x.intersp = 0.5, y.intersp = 0.80, cex = 0.90)
mtext(paste0(country), side = 3, line = 0.5, cex = 0.85)

dev.off()


    ## ADDITIONAL FIGURE - Evolution of incidence of 65+ by country

age <- "65+"
sg <- "all"

df.inc <- filter(df.ecdc,
                 Age == age,
                 Variable == "incidence",
                 By %in% c("age", "age*serogroup"),
                 Serogroup == sg,
                 Year >= 2008
)


spread.df.inc <- spread(df.inc[,c(2,3,10)], key = "Year", value = Value)

row.names(spread.df.inc) <- spread.df.inc[,1]
df.plot <- spread.df.inc[-1]

plot(colnames(df.plot),
     df.plot[1,],
     type = "n",
     ylab = "",
     ylim = c(0,max(df.plot, na.rm = T)),
     xlab = ""
)

for (i in seq(nrow(df.plot))) {
    lines(colnames(df.plot), df.plot[i,], type = "p", col = color.country[i], pch = 16)
    lines(lowess(colnames(df.plot), df.plot[i,], f = 0.6),col = color.country[i], type = "b")
}

text(2017, df.plot$`2017`, row.names(df.plot), cex = 0.7, pos = 4)

df.plot.evo <- df.plot/df.plot[,1]
df.plot.evo.na <-drop_na(df.plot.evo)
df.plot.evo[df.plot.evo<1]<- 1-1/df.plot.evo[df.plot.evo<1]
df.plot.evo[df.plot.evo>=1]<- df.plot.evo[df.plot.evo>=1]-1

plot(colnames(df.plot.evo),
     df.plot.evo[1,],
     type = "n",
     ylab = "",
     ylim = c(-3,5),
     xlab = ""
)

for (i in seq(nrow(df.plot.evo))) {
    #lines(colnames(df.plot), df.plot[i,], type = "p", col = color.age[i], pch = 16)
    lines(lowess(colnames(df.plot.evo), df.plot.evo[i,], f = 0.75),col = color.country[i], type = "b")
}
abline(h=0, col= "red", lty=4)


### CORRELATION Pop vs incidence

df.pop <- filter(df.population,
                 Age == "65+",
                 Variable == "population",
                 By  == "age",
                 Year >= 2007
                 )

spread.df.pop <- spread(df.pop[c(1, 3,10)], key = Year, value = Value)
evo.df.pop <- spread.df.pop[,-1]/spread.df.pop[,2]
row.names(evo.df.pop) <- spread.df.pop[,1]

df.inc <- filter(df.ecdc,
                 Age == "65+",
                 Variable == "incidence",
                 By == "age",
                 Serogroup == "all",
                 Year >= 2007
)


spread.df.inc <- spread(df.inc[,c(1,3,10)], key = "Year", value = Value)
evo.df.inc <- spread.df.inc[,-1]/spread.df.inc[,2]
row.names(evo.df.inc) <- spread.df.inc[,1]

plot(evo.df.pop[,8],
     evo.df.inc[,8]
     )



## SUPPLEMENTAL FIGURE 4 - Incidence by serogroup for the different age groups
### Setting the graphical areas   
par(mfrow=c(4,2),
    oma = c(4, 2, 2, 2),
    mar=c(2.1, 4.1, 2.1, 1),
    cex = 0.60,
    lwd = 0.60
)
### Selection of variables
year <- 2017
age <- "<01"

### For loop selecting the age groups
for (i in c("<01", "01-04", "05-14", "15-24", "25-49", "50-64", "65+")) {
  age <- i
  
  ### Data sorting
  df.inc <- filter(df.ecdc,
                   Age == age,
                   Variable == "incidence",
                   By %in% c("age*serogroup"),
                   Serogroup %in% c("B", "C", "W", "Y"),
                   Year == year
  )
  
  df.inc <- arrange(df.inc, Serogroup)
  
  ### Generating the box plot
  boxplot(Value~as.factor(Serogroup),
          data=df.inc,
          outline = FALSE,
          #ylim=c(0,ifelse(age %in% c("05-14", "15-24", "25-49", "50-64", "65+"),1, max(df.inc$Value, na.rm=TRUE)*1.05)),
          ylim=c(0,max(df.inc$Value, na.rm=TRUE)*1.05),
          ylab="Incidence (cases/100,000/year)",
          border = "gray70"
  )
  
  ### Generating the points layover
  beeswarm(Value~as.factor(Serogroup),
           data=df.inc,
           add = TRUE,
           pch = 20,
           col = "black",
           method = "swarm",
           spacing = 0.80
  )
  
  #### Calulate the number of cases by age group in EU/EEA using the total number of case and the proportion by age group for each individual country
  df.cases.all <- filter(df.ecdc,
                         Age == "all",
                         Variable == "cases",
                         By == "serogroup",
                         Serogroup %in% c("B", "C", "W", "Y"),
                         Year == year
  )
  
  df.age.prop <- filter(df.ecdc,
                        Age == age,
                        Variable == "distribution",
                        By == "age_in_serogroup",
                        Serogroup %in% c("B", "C", "W", "Y"),
                        Year == year
  )
  
  spread.df.cases.all <- spread(df.cases.all[,c(1,5,10)], key = "Serogroup", value = Value)
  spread.df.age.prop <- spread(df.age.prop[,c(1,5,10)], key = "Serogroup", value = Value)
  spread.df.age.prop[spread.df.cases.all == 0] <- 0
  df.age.cases <- data.frame(Country = spread.df.cases.all[1], (spread.df.cases.all[-1]*spread.df.age.prop[-1])/100)
  df.sum.cases <- apply(df.age.cases[2:5], 2, function(i) sum(i, na.rm = TRUE)) #calculate the sum of the number of cases
  
  #### Calulate the population by age group in EU/EEA taking in consideration the NA value in the cases dataframe
  df.pop <- filter(
    df.population,
    Age == age,
    By  == "age",
    Year == year
  )
  
  spread.pop.ordered <- left_join(df.age.cases[1], df.pop[c(1,10)], by = "Country") #align the population data with the order use in the case data frame
  spread.pop.ordered[2][is.na(df.age.cases[2])] <- NA #remove the population values for countries without IMD cases data
  df.sum.pop <- sum(spread.pop.ordered[-1], na.rm = TRUE) #calculate the sum of the population
  
  #### Calculate the incidence by age group
  EU.inc <- (df.sum.cases/df.sum.pop)*100000
  
  #### plot the points and text
  points(1:4, EU.inc, col = "red", pch = 18, cex = 1)
  text(1:4, 0, pos = 1, offset = 0.25, labels = paste0("n=", round(df.sum.cases)), col = "red", cex = 0.70)
  
  ### Plotting the name of outliers
  split.outliers <- split(df.inc, df.inc[,c("Serogroup")])
  q.sg <- quantile(df.inc$Value, na.rm = T)
  offset <- (q.sg[5]-q.sg[1])*0.003 #use a threshold for printing the name of outliers
  df.outliers <- lapply(split.outliers,
                        function(i) {
                          q <- quantile(i$Value, na.rm = TRUE)
                          filter(i, Value > (q[4]+1.5*(q[4]-q[2]))+offset)
                        })
  
  if(nrow(df.outliers$`B`) != 0) {text(1, df.outliers$`B`$Value, df.outliers$`B`$ISO3, pos = 4, cex = 0.75)}
  if(nrow(df.outliers$`C`) != 0) {text(2, df.outliers$`C`$Value, df.outliers$`C`$ISO3, pos = 4, cex = 0.75)}
  if(nrow(df.outliers$`W`) != 0) {text(3, df.outliers$`W`$Value, df.outliers$`W`$ISO3, pos = 4, cex = 0.75)}
  if(nrow(df.outliers$`Y`) != 0) {text(4, df.outliers$`Y`$Value, df.outliers$`Y`$ISO3, pos = 4, cex = 0.75)}
  
  legend("topright", col = c("black", "red"), pch = c(20,18), legend = c("individual county", "EU/EEA"), bty = "n", cex = 0.90)
  mtext(side=3, line=0.5, cex=0.85, adj = 0.5, outer=F, age)
}

### Formating the mulpitle frame
pdf(paste0("Supplemental Figure 4.pdf"), width=8.27, height=11.69, paper = "special", title = "Rplot")

dev.off()



## Supplemental Figure 5 - Evolution of proportion of cases for specific age group  
### Setting the graphical areas   

#k <- "all"
list.prop.df.age <- list()

for (k in c("all", "B", "C", "W", "Y")) {
  
  ### Generation of dataframe for individual countries
  df.cases.all <- filter(df.ecdc,
                         Age == "all",
                         Variable == "cases",
                         By %in% c("overall", "serogroup"),
                         Serogroup == k,
                         Year >= 2008
  )
  
  df.prop.age <- filter(df.ecdc,
                        Age %in% c("<01", "01-04", "05-14", "15-24", "25-49", "50-64", "65+"),
                        Variable == "distribution",
                        By %in% c("age_in_overall", "age_in_serogroup"),
                        Serogroup == k,
                        Year >= 2008
  )
  spread.df.cases.all <- spread(df.cases.all[c(1,3,10)], key = Year, value = Value)
  
  
  split.df.prop.age <- split(df.prop.age, df.prop.age$Age)
  spread.df.prop.age <- lapply(split.df.prop.age, function(i) spread(i[c(1,3,10)], key = Year, value = Value))
  df.age <- lapply(spread.df.prop.age, function(i) {(i[-1]*spread.df.cases.all[-1])/100})
  sum.df.age <- as.data.frame(sapply(df.age, function(i) {apply(i, 2, function(j) {sum(j, na.rm = T)})}))
  sum.df.age$Sum <- rowSums(sum.df.age)
  prop.df.age <- (sum.df.age/sum.df.age$Sum)*100
  
  list.prop.df.age[[k]] <- t(prop.df.age)
}


par(mfrow=c(4,1),
    mar=c(2.1, 4.1, 3.1, 1),
    oma = c(1, 1, 1, 1),
    cex = 0.75,
    lwd = 0.75
)

for (k in c("<01", "01-04", "15-24", "65+")) {
  
  df.plot <- t(sapply(list.prop.df.age, function(i) {i[k,]}))
  
  plot(colnames(list.prop.df.age$all),
       list.prop.df.age$all[1,],
       type = "n",
       xlab = NA,
       ylab = "proportion",
       ylim = c(0,50)
  )
  
  for (i in seq(nrow(df.plot))) {
    lines(colnames(df.plot), df.plot[i,], type = "p", col = color.sg[i], pch = 16)
    lines(lowess(colnames(df.plot), df.plot[i,], f = 0.75), col = color.sg[i], lwd=1, type = "b")
  }
  
  legend("topright", row.names(df.plot), col=color.sg[1:5], pch=16, bty = "n", cex=0.9, x.intersp = 0.5, y.intersp = 1, pt.cex = 1)
  mtext(side=3, line=0.5, cex=0.85, adj = 0.5, outer=F, k)
}

### Formating the mulpitle frame
pdf(paste0("Supplemental Figure 5.pdf"), width=8.27, height=11.69, paper = "special", title = "Rplot")

dev.off()



### LINEAR MODEL

summary(lm(EU.cases[4:7] ~ as.numeric(names(EU.cases[4:7]))))
