####-Importing Libraries and Set Up ----

# Set the working 

# Loading Packages 
library(tidyverse)
library(ggplot2)
library(rentrez)
library(seqinr)
library(dplyr)
library(BiocManager)
library(Biostrings)
library(randomForest)
library(e1071)
library(gbm)

####Section1: Data Acquisition and Cleaning ----

#####_1.Data Acquisition ----

# Data was taken from the NCBI nucleotide database on March 15, 2025
# The data was extracted for each Genus separately so that they were sampled equally.
# The search was specified to only take COI gene sequence's. Amblyomma had the least number of samples (396) so every genus was set to have this maximum number of samples to ensure a balanced data set.

# Ixodes Genus 
Ixodes_data <- entrez_search(db = "nuccore", term = "(Ixodes[ORGN] AND COI[GENE])", retmax = 396, use_history=TRUE)

# Haemaphysalis Genus
Haem_data <- entrez_search(db = "nuccore", term = "(Haemaphysalis[ORGN] AND COI[GENE])", retmax = 396, use_history=TRUE)

# Rhipicephalus Genus
Rhip_data <- entrez_search(db = "nuccore", term = "(Rhipicephalus[ORGN] AND COI[GENE])", retmax = 396, use_history=TRUE)

# Amblyomma Genus
Amb_data <- entrez_search(db = "nuccore", term = "(Amblyomma[ORGN] AND COI[GENE])", retmax = 396, use_history=TRUE)

#####_2.Organizing Data----

# Ixodes
# The fasta sequences were extracted for each unique ID
Ixodes_fetch <- entrez_fetch(db = "nuccore", web_history = Ixodes_data$web_history, rettype = "fasta", retmax = 396)
# This file is saved to the working directory set 
write(Ixodes_fetch, "Ixodes_fetch.fasta", sep = "\n") 
# The fasta file was imported from the data file:
stringSet1 <- readDNAStringSet("Ixodes_fetch.fasta")
# The fasta sequence were made into a DNA strings set
df_Ixodes <- data.frame(Identifier = names(stringSet1), Nucelotide_Sequence = paste(stringSet1))
# The species name was taken from each ID and put into a new column 
df_Ixodes$Species_Name <- word(df_Ixodes$Identifier, 2L, 3L)
df_Ixodes <- df_Ixodes[ , c("Identifier", "Species_Name", "Nucelotide_Sequence")]
# Before the identifier names were cleaned, each was checked to make sure that the sequences were only from COI gene and the mitochondria 
# Then the identifier column was edited to only keep the unique identifier that will be needed for the ML models later in order to separate samples so they are not repeated between the training and validation set.
# The names were edited to only keep the Genus names
df_Ixodes$Identifier <- sub("^(\\S+).*", "\\1", df_Ixodes$Identifier)
df_Ixodes$Species_Name <- sub("^Ixodes.*", "Ixodes", df_Ixodes$Species_Name)
# Then the column name was changed to Genus name instead of species name
names(df_Ixodes)[names(df_Ixodes) == "Species_Name"] <- "Genus_Name"

# The same process was repeated for the otehr Genus:

# Haemaphysalis
Haem_fetch <- entrez_fetch(db = "nuccore", web_history = Haem_data$web_history, rettype = "fasta", retmax = 396)
write(Haem_fetch, "Haem_fetch.fasta", sep = "\n") 
stringSet2 <- readDNAStringSet("Haem_fetch.fasta")
df_Haem <- data.frame(Identifier = names(stringSet2), Nucelotide_Sequence = paste(stringSet2))
df_Haem$Species_Name <- word(df_Haem$Identifier, 2L, 3L)
df_Haem <- df_Haem[ , c("Identifier", "Species_Name", "Nucelotide_Sequence")]
df_Haem$Identifier <- sub("^(\\S+).*", "\\1", df_Haem$Identifier)
df_Haem$Species_Name <- sub("^Haem.*", "Haem", df_Haem$Species_Name)
names(df_Haem)[names(df_Haem) == "Species_Name"] <- "Genus_Name"

# Rhipicephalus
Rhip_fetch <- entrez_fetch(db = "nuccore", web_history = Rhip_data$web_history, rettype = "fasta", retmax = 396)
write(Rhip_fetch, "Rhip_fetch.fasta", sep = "\n") 
stringSet3 <- readDNAStringSet("Rhip_fetch.fasta")
df_Rhip <- data.frame(Identifier = names(stringSet3), Nucelotide_Sequence = paste(stringSet3))
df_Rhip$Species_Name <- word(df_Rhip$Identifier, 2L, 3L)
df_Rhip <- df_Rhip[ , c("Identifier", "Species_Name", "Nucelotide_Sequence")]
df_Rhip$Identifier <- sub("^(\\S+).*", "\\1", df_Rhip$Identifier)
df_Rhip$Species_Name <- sub("^Rhip.*", "Rhip", df_Rhip$Species_Name)
names(df_Rhip)[names(df_Rhip) == "Species_Name"] <- "Genus_Name"

# Amblyomma
Amb_fetch <- entrez_fetch(db = "nuccore", web_history = Amb_data$web_history, rettype = "fasta", retmax = 396)
write(Amb_fetch, "Amb_fetch.fasta", sep = "\n") 
stringSet4 <- readDNAStringSet("Amb_fetch.fasta")
df_Amb <- data.frame(Identifier = names(stringSet4), Nucelotide_Sequence = paste(stringSet4))
df_Amb$Species_Name <- word(df_Amb$Identifier, 2L, 3L)
df_Amb <- df_Amb[ , c("Identifier", "Species_Name", "Nucelotide_Sequence")]
df_Amb$Identifier <- sub("^(\\S+).*", "\\1", df_Amb$Identifier)
df_Amb$Species_Name <- sub("^Amb.*", "Amb", df_Amb$Species_Name)
names(df_Amb)[names(df_Amb) == "Species_Name"] <- "Genus_Name"

# Then the 4 genus data frames were merged vertically so all the rows for each variable were put together 
full_df <- bind_rows(df_Rhip, df_Ixodes, df_Haem, df_Amb)
unique(full_df$Genus_Name)

#####_3.Cleaning Sequences----

# Calculate the Sequence statistics 

# Count number of un-identified bases marked as N
full_df <- full_df %>% mutate(N_count = str_count(full_df$Nucelotide_Sequence, 'N'))

# The sequence length 
full_df <- full_df %>% mutate(Length = nchar(full_df$Nucelotide_Sequence))

# The percentage of un-identified bases
full_df <- full_df %>% mutate(N_count_percentage = (N_count/Length) *100)

# The number of missing nucletoides
full_df <- full_df %>% mutate(Missing_data = str_count(full_df$Nucelotide_Sequence, '-'))

# Remove samples that had over 1% of their sequence with 
full_df <- full_df %>% filter(N_count_perentage < 0.1)

# Looking at the sequence length
summary(full_df$Length)
ggplot(full_df, aes(x=Length)) + geom_histogram()

# Look at length distributions of each genus 
boxplot(Length ~ Genus_Name, data= full_df, main = "Distribution Sequnce Length", ylab = "Length (bp)", xlab = 'Genus Name', col = c( "lightblue", "lightgreen", "salmon", "yellow"))

# Ticks typically have a COI gene of 580 base pairs (bp). The sequence lengths of each genus were examined in Figure.1 The sequences were filtered to only include those between 500 and 700 bp to account for the variation in genus length and differences in PCR amplification of the gene. 
full_df <- full_df %>% filter(Length >= 550 & Length <= 700)

# Check filtering worked 
table(full_df$Genus_Name)

####Section2: Feature Extraction ----
# To separate the genus based on sequences, the properties of the sequences will need to first be determined. The 3 K-mer frequency will be calculated. This calculates the frequency of 3 nucleotide’s appearing in a row within a sequences.

#First the nucleotide's need to be converted into biostrings so that not only character functions can be used. 
full_df <- as.data.frame(full_df)
full_df$Nucelotide_Sequence <- DNAStringSet(full_df$Nucelotide_Sequence)
class(full_df$Nucelotide_Sequence)

# Calculating the 3 K-Mer Frequency, as.prob = TRUE was used so it accounts for the variability in sequence length. 
full_df <- cbind(full_df, as.data.frame(trinucleotideFrequency(full_df$Nucelotide_Sequence, as.prob = TRUE)))
full_df$Nucelotide_Sequence <- as.character(full_df$Nucelotide_Sequence)

####Section3: Generating Data sets for ML Models ----

# Defining the data sets for the ML models 

# First a validation data set must be made, this data will be remain hidden during the training process so that the algorithm cannot see the sequences and they will be used to validate the correct classification.  It is important that a random sample of the same number of sequences are taken from each genus. Around 25% of sequences for each genus will be set aside for the validation set. The data will grouped by genus and 75 random samples from each will be taken.  The seed will be set so that this can be reproduced and the same samples will always be taken. 
# Validation 
set.seed(217)
df_Validation <- full_df %>%
  group_by(Genus_Name) %>%
  sample_n(75)

# Check all 4 Genus have the same number of samples 
table(df_Validation$Genus_Name)

# Now the training data  will be made, for this it is crucial that the data that was already taken for the validation set is not included in the training. The data is filtered with "!" so that it knows not to choose ones already in the validation data. This is why the unique NCBI ID was kept for each sequence so they could be separated. The training data consisted of 235 samples from each Genus 
set.seed(13)
df_Training <- full_df %>%
  filter(!Identifier %in% df_Validation$Identifier) %>%
  group_by(Genus_Name) %>%
  sample_n(235)

## Check all 4 Genus have the same number of samples 
table(df_Training$Genus_Name)

####Section4: ML Models ----

#####_1.Random Forest ----

#The random forest (RF) ML algorithm creates multiple decision trees and merges them together to vote on how to classify a data point, in this case by a genus. The algorithm takes a random subset of the data to build each tree to avoid specific or bias learning and keeps the model generalizable. As the tress grow and learn the patterns it will determine which features (3-Mer’s) are the most influential.

# There are 63 predictors or K-mer frequencies being used in this model. This is a relatively high number and in turn the model was set to make 130 decision trees [14]. A greater number of descension trees promotes accurate predictions and ensures that each sample is left out of the tress enough times so that it has the chance to be predicted. 

# Training the model using the K-Mer frequencies in the training data set 
genus_classifier_RF1 <- randomForest::randomForest(x = df_Training[, 8:71], y = as.factor(df_Training$Genus_Name), ntree = 130, importance = TRUE) 

# Looking at how the model trained first 
genus_classifier_RF1$confusion
# The rows show the true classes and the columns shows prediction, in this training each sample was assigned to the correct genus.
# There was a 0.96% estimated error rate

summary(genus_classifier_RF1$oob.times)
# This shows the amount of times each sample was left out when building a decision tree. This increases the generality of the model.  The sample can only be predicted when it is not included. The average time each element could be predicted is 66 times. 

genus_classifier_RF1$votes
# Not every row has agreement for only one genus meaning that not all trees had 100% agreement when trying to classify each observation. Despite this the classifier was able to separate the training data samples correctly in each genus close to perfect. 

# Now the trained model was tested on the validation set 
predictValidation_RF1 <- predict(genus_classifier_RF1, df_Validation[, c(2, 8:71)])
predictValidation_RF1

# Looking at the results, there was one miss-classification 
table(observed = df_Validation$Genus_Name, predicted = predictValidation_RF1)

# Now a new model was trained but with a higher number of tress 
genus_classifier_RF2 <- randomForest::randomForest(x = df_Training[, 8:71], y = as.factor(df_Training$Genus_Name), ntree = 180, importance = TRUE) 

genus_classifier_RF2$confusion
# Saw the error rate during training was reduced to 0.85%. 
summary(genus_classifier_RF2$oob.times)

genus_classifier_RF2$votes

# Now the new trained model was tested on the validation set 
predictValidation_RF2 <- predict(genus_classifier_RF2, df_Validation[, c(2, 8:71)])

# Again there was one incorrect classification using more trees 
table(observed = df_Validation$Genus_Name, predicted = predictValidation_RF2)

confusion_matrix <- table(observed = df_Validation$Genus_Name, predicted = predictValidation_RF2)
confusion_df <- as.data.frame(confusion_matrix)

# A heat plot is made for the findings 
ggplot(data = confusion_df, aes(x = predicted, y = observed)) +
  geom_tile(aes(fill = Freq), color = "white") +
  scale_fill_gradient(low = "lightblue", high = "dodgerblue") +
  geom_text(aes(label = Freq), color = "white") +
  labs(title = "Confusion Matrix of Models\n on Validation Data", x = "Predicted", y = "Observed") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))

# RF preformed the same with 99.6% accuracy with both number of tress. Increasing the number of trees only minimized the confusion during the training process. Overall training and validating the model only takes a few seconds. 

#####_2.Support Vector Machine ----

# In support vector machine (SVM), the algorithm separates the data across a boundary in 3D space to find the greatest distance to separate these groups. The algorithm will find the optimal boundary to separate the points that lie close together and form different clusters

# The radial shape is standard for data that is not linear so the decision boundary can find its way through the data without having to take on the shape of a straight line

# The x values in the model were defined as the 3-mer values in columns 8-71.

# y values were the factored-out genus names.

# The gamma value defines how much influence each data point has on the model during the training process. This value was set to 0.01 to ensure the decision boundary was not too specific to these data points and keep it generalizable. A lower value also lowers the complexity of the decision boundary and in turn lowers the number of points that lie on the decision boundary that are difficult to separate 

# The cost value represents the deduction implemented when a data point is classified incorrectly across the decision boundary during training

genus_classifier_SVM1 <- svm(
  x = df_Training[, 8:71],   
  y = as.factor(df_Training$Genus_Name), 
  kernel = "radial",   # A non linear shape           
  cost = 10,     # Prioritize high accuracy                         
  gamma = 0.01) # Keep it generalizable                             

genus_classifier_SVM1
# There was 156 data points lying on the decision plane 

# Now use the model on the validation set 
predictValidation_SVM1 <- predict(genus_classifier_SVM1, df_Validation[, 8:71])

# View the predictions
predictValidation_SVM1
confusion_matrix_SVM1 <- table(Observed = df_Validation$Genus_Name, Predicted = predictValidation_SVM)
print(confusion_matrix_SVM1)

# There was one miss-classification using these parameters - 99% accuarte 

# Now a new model was made changing the gamma parameter to 0.1 which makes the decision boundary more complex and there will be more emphases on each point, 
genus_classifier_SVM2 <- svm(
  x = df_Training[, 8:71],             
  y = as.factor(df_Training$Genus_Name),
  kernel = "radial", 
  cost = 10, 
  gamma = 0.1)                       

genus_classifier_SVM2
# There is now 249 data points laying on the decision boundary 

# Predict on validation set
predictValidation_SVM2 <- predict(genus_classifier_SVM2, df_Validation[, 8:71])
predictValidation_SVM2

confusion_matrix_SVM2 <- table(Observed = df_Validation$Genus_Name, Predicted = predictValidation_SVM2)
print(confusion_matrix_SVM2)

# There was 2 miss-classifications bringing the accuracy down to 97%.

# Now the gamma value was returned to 0.01 that showed better performance and the cost value was lowered to 1. 
genus_classifier_SVM3 <- svm(
  x = df_Training[, 8:71],             
  y = as.factor(df_Training$Genus_Name), 
  kernel = "radial",               
  cost = 1,   # This puts less importance on each data point                         
  gamma = 0.01)         

genus_classifier_SVM3

predictValidation_SVM3 <- predict(genus_classifier_SVM3, df_Validation[, 8:71])
confusion_matrix_SVM3 <- table(Observed = df_Validation$Genus_Name, Predicted = predictValidation_SVM3)
print(confusion_matrix_SVM3)

# There was 1 miss-classification, bring the accuray to 99%. Changing the cost vaue did not impact model performance like the gamma value. 

# Now the model shape was changed to polynomial 
genus_classifier_SVM4 <- svm(
  x = df_Training[, 8:71],             
  y = as.factor(df_Training$Genus_Name), 
  kernel = "polynomial",                
  cost = 10,                             
  gamma = 0.01)         

predictValidation_SVM4 <- predict(genus_classifier_SVM4, df_Validation[, 8:71])
confusion_matrix_SVM4 <- table(Observed = df_Validation$Genus_Name, Predicted = predictValidation_SVM4)
print(confusion_matrix_SVM4)

# This resulted in a lot of miss-classification of sequences from the Rhipicephalus genus. 

# A heat map is made for the optimal results using the first model. 
confusion_df2 <- as.data.frame(confusion_matrix_SVM1)
ggplot(data = confusion_df, aes(x = Predicted, y = Observed)) +
  geom_tile(aes(fill = Freq), color = "white") +
  scale_fill_gradient(low = "lightcoral", high = "firebrick") +
  geom_text(aes(label = Freq), color = "white") +
  labs(title = "SVM Model", x = "Predicted", y = "Observed") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))

# The gamma value at 0.01 gave optimal results for this model. The cost value did not have an impact on the model when both lowered or raised from the value 1. The effects of the parameters are summarized in Table 2 of the report. This model preformed best when the decision boundary was kept simple and could be trained and validated quickly. 

#####_3.Gradient Boosting Machine ----

# The gradient boosting machine (GBM) algorithm like FR builds a series of decision trees. However this model builds a tree then attempts to make a prediction, the following tree built will learn and correct for the errors in the prior prediction. This self correcting tool allows the model to learn as it goes providing high accuracy. GBM also offers cross validation during the training stage where it splits the data into sets to train and test its predictability 

# Training the GBM Model 
df_Training$Genus_Name <- as.factor(df_Training$Genus_Name)

# The Models paramters 

# Multinomial wadue to the 4 different genus it will be classifying

# The number of tress was set to 593 as suggested by the algorithm as this was the optimal value, any trees created past this value were not significant and could create overfitting to the data. 

# The interaction depth was set to 5, specifying that a maximum of 5 features can interact in each tree. A moderate value was chosen as higher values increase the complexity of interactions and can again lead to over fitting.

# The parameter n.minobsinnode represents the number of observations that must be found at the last node of each tree. 10 observations were used as opposed to the lower default to generate more simple trees that would not be specific to a few samples. 

# Shrinkage refers to the time each tree will take to learn, a smaller value often yields more accurate results. This value was set to a small value 0.01, although typically this would require more trees in the model. 

# Finally the cross validation cv.fold was set to 5. This mean the data was trained on 4 subsets of the data and tested on the remaining subset, this was then repeated 4 more times. 

# Train the GBM model using only selected columns
genus_classifier_GBM1 <- gbm(
  formula = Genus_Name ~ .,               # Target variable: Genus_Name
  data = df_Training[, c(8:71, 2)],       # Use columns 8:71 for features and column 2 for labels
  distribution = "multinomial",           # Multiclass classification
  n.trees = 593,                         # Number of boosting iterations (trees)
  interaction.depth = 5,                  # Depth of each tree
  shrinkage = 0.01,                      # Learning rate 
  cv.folds = 5,                          # Cross-validation folds
  n.minobsinnode = 10                    # Minimum samples per leaf node
)

# View the model summary
summary(genus_classifier_GBM1)
genus_classifier_GBM1

# Use the trained model to make predictions 
predictValidation_GBM1 <- predict(genus_classifier_GBM, df_Validation[, 8:71], n.trees = 593, type = "response")
predictValidation_GBM1

# Convert probabilities to class labels
predicted_class_GBM1 <- colnames(predictValidation_GBM1)[apply(predictValidation_GBM1, 1, which.max)]

conf_matrix_GBM <- table(Observed = df_Validation$Genus_Name, Predicted = predicted_class_GBM1)

# Print confusion matrix
print(conf_matrix_GBM1)

# 1 miss classification, overall 99% accuracy. 

# To see if perfect classification could be obtained the cross validation value was changed to 10.
genus_classifier_GBM2 <- gbm(
  formula = Genus_Name ~ .,               
  data = df_Training[, c(8:71, 2)],    
  distribution = "multinomial",           
  n.trees = 593,                         
  interaction.depth = 5,                 
  shrinkage = 0.01,                     
  cv.folds = 10,                          
  n.minobsinnode = 10                    
)

# View the model summary
summary(genus_classifier_GBM2)
genus_classifier_GBM2

# Use the trained model to make predictions 
predictValidation_GBM2 <- predict(genus_classifier_GBM2, df_Validation[, 8:71], n.trees = 593, type = "response")
predictValidation_GBM2

# Convert probabilities to class labels
predicted_class_GBM2 <- colnames(predictValidation_GBM2)[apply(predictValidation_GBM2, 1, which.max)]
conf_matrix_GBM2 <- table(Observed = df_Validation$Genus_Name, Predicted = predicted_class_GBM2)

# Print confusion matrix
print(conf_matrix_GBM2)

# Increasing the number of cross validations did not change the accuracy of the model but increased the time it took to train the model.

# The cross validation value was changed back to 5 and the shrinkage was reduced to 0.001 to allow for a longer leaning time.
genus_classifier_GBM3 <- gbm(
  formula = Genus_Name ~ .,               
  data = df_Training[, c(8:71, 2)],    
  distribution = "multinomial",           
  n.trees = 593,                         
  interaction.depth = 5,                 
  shrinkage = 0.001,                     
  cv.folds = 5,                          
  n.minobsinnode = 10                    
)

# View the model summary
summary(genus_classifier_GBM3)
genus_classifier_GBM3

# Use the trained model to make predictions 
predictValidation_GBM3 <- predict(genus_classifier_GBM2, df_Validation[, 8:71], n.trees = 593, type = "response")
predictValidation_GBM3

# Convert probabilities to class labels
predicted_class_GBM3 <- colnames(predictValidation_GBM3)[apply(predictValidation_GBM3, 1, which.max)]
conf_matrix_GBM3 <- table(Observed = df_Validation$Genus_Name, Predicted = predicted_class_GBM3)

# Print confusion matrix
print(conf_matrix_GBM3)

# Lowering the shrinkage value from 0.01 to 0.001 to increase the learning time lowered the accuracy only slightly to 99.3% contrary to what usually occurs with this parameter.

# This final matrix was made using the results in model 3. 
confusion_df3 <- as.data.frame(conf_matrix_GBM)
ggplot(data = confusion_df2, aes(x = Predicted, y = Observed)) +
  geom_tile(aes(fill = Freq), color = "white") +
  scale_fill_gradient(low = "plum", high = "mediumpurple") +
  geom_text(aes(label = Freq), color = "white") +
  labs(title = "GBM Model", x = "Predicted", y = "Observed") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))

# The effects of the parameters for this model are summarized in Table 3. This model took much longer to train and validate compared to random forest. 

####Section5: Conclusions ----

# An in depth comparison of the models can be found in the report. Ultimately each model preformed the same with 99% accuracy. Each ML model was trained on the data with parameters to increase generalizability and to reduce over-fitting to this specific data set. This was to ensure that these models could be used as a classification tool outside of this data and in the field. Although no model can be chosen to produce more accurate results, each had their unique strengths and weaknesses
