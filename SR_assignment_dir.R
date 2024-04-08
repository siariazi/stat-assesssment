# Loading the libraries 
# this is for personal use and has directories 
library(dplyr)
library(tidyr)
library(lubridate)
library(broom)
library(epitools)
library(reshape2)
library(ggplot2)
library(lme4)
library(car)
library(mgcv)
library(survival)
library(survminer)
library(scales)
library(gridExtra)
library(survminer)

# Specifying the working directory
directory <- "/Users/siavashriazi/Personal/Job/Bohdan\ Nosyk/Bohdan_data"
setwd(directory)

# A directory for saving plots
save_dir <- "/Users/siavashriazi/Dropbox/Apps/Overleaf/Assignment 2/Figures"

# Load the data
interview_data <- read.csv("assess-interviews.csv")
pvl_data <- read.csv("assess-pvl.csv")

# Correct the date conversion
interview_data$dt.int <- as.Date(interview_data$dt.int, format = "%m/%d/%Y")
interview_data$dt.dob <- as.Date(interview_data$dt.dob, format = "%m/%d/%Y")
pvl_data$dt.tst <- as.Date(pvl_data$dt.tst, format = "%Y-%m-%d")

# Calculate age at the time of interview
interview_data$age <- interval(interview_data$dt.dob, interview_data$dt.int) / years(1)

# Get the baseline (first) interview for each participant
baseline_interviews <- interview_data %>%
  arrange(id, dt.int) %>%
  group_by(id) %>%
  slice(1) %>%
  ungroup()

# Merge the baseline interviews with the pvl data to include plasma HIV-1 RNA viral load
baseline_pvl <- pvl_data %>%
  arrange(id, desc(dt.tst)) %>%
  group_by(id) %>%
  slice(1) %>%
  ungroup()

# Merge datasets on participant ID
baseline_data <- merge(baseline_interviews, baseline_pvl[, c("id", "p")], by = "id")

# Rename 'p' to 'plasma_viral_load' for clarity
names(baseline_data)[names(baseline_data) == "p"] <- "plasma_viral_load"

# First, ensure that 'homeless' is a factor for clearer plots
baseline_data$homeless <- factor(baseline_data$homeless, levels = c(0, 1), labels = c("Housed", "Homeless"))

# Reshape the data to have a long format for easier plotting with ggplot2
baseline_data_long <- melt(baseline_data, id.vars = c("id", "homeless"), measure.vars = c("age", "c", "plasma_viral_load"))

# Rename the columns for clarity
colnames(baseline_data_long) <- c("ID", "Residence", "Variable", "Value")

# Assuming baseline_data_long is already created without log transformation for RNA
baseline_data_long$Variable <- factor(baseline_data_long$Variable,
                                      levels = c("age", "c", "plasma_viral_load"),
                                      labels = c("Age", "CD4+", "Plasma HIV-1 RNA"))

# Plotting housing Age, CD4 and HIV RNA across housing status with y-axis on a logarithmic scale for RNA
plot_housing <- ggplot(baseline_data_long, aes(x = Residence, y = Value, fill = factor(Residence))) + 
  geom_boxplot(alpha = 0.5) + 
  facet_wrap(~ Variable, scales = "free_y") + 
  labs(y = "", 
       x = "") + 
  theme_minimal() + 
  scale_fill_manual(values = c('red', 'blue')) + 
  theme(legend.position = "none") +
  scale_y_continuous(trans = 'log10', labels = label_comma()) + # Log scale for y-axis
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Improve x-axis label readability

show(plot_housing)
ggsave(plot_housing, filename = "plot_housing_status.png", path = save_dir, width = 5, height = 4, dpi = 300)

# Extract year from the interview date
baseline_data$interview_year <- format(baseline_data$dt.int, "%Y") %>% as.numeric()

# Calculate median and IQR for continuous variables
continuous_summary <- baseline_data %>%
  group_by(homeless) %>%
  summarise(
    Age_Median = round(median(age, na.rm = TRUE), 2),
    Age_IQR = round(IQR(age, na.rm = TRUE), 2),
    CD4_Median = round(median(c, na.rm = TRUE), 2),
    CD4_IQR = round(IQR(c, na.rm = TRUE), 2),
    RNA_Median = round(median(plasma_viral_load, na.rm = TRUE), 2),
    RNA_IQR = round(IQR(plasma_viral_load, na.rm = TRUE), 2),
    Year_Interview_Median = median(interview_year, na.rm = TRUE),
    Year_Interview_IQR = IQR(interview_year, na.rm = TRUE)
  )

print(continuous_summary)

# Save the continuous_summary DataFrame to a CSV file
write.csv(continuous_summary, "continuous_summary.csv", row.names = FALSE)

# Define the binary variables to analyze
binary_vars <- c("female", "coke", "crac", "hero", "mmt")

# Initialize an empty data frame for binary variable summaries
binary_summary_df <- data.frame(
  Variable = character(),
  Housed_Count = integer(),
  Homeless_Count = integer(),
  Odds_Ratio = numeric(),
  CI_Lower = numeric(),
  CI_Upper = numeric(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each binary variable to calculate statistics and fill the dataframe
for (var in binary_vars) {
  # Create contingency table
  tbl <- table(baseline_data$homeless, baseline_data[[var]])
  
  # Calculate counts for non-homeless (0) and homeless (1)
  housed_count <- tbl[1, 2]  # Assuming 1st row is '0' for non-homeless, 2nd column is '1' for presence of variable
  homeless_count <- tbl[2, 2]      # Assuming 2nd row is '1' for homeless, 2nd column is '1' for presence of variable
  
  # Perform Fisher's Exact Test for Odds Ratio, CI, and P-Value
  fisher_result <- fisher.test(tbl)
  
  # Append to the dataframe
  binary_summary_df <- rbind(binary_summary_df, data.frame(
    Variable = var,
    Housed_Count = housed_count,
    Homeless_Count = homeless_count,
    Odds_Ratio = round(fisher_result$estimate, 3),
    CI_Lower = round(fisher_result$conf.int[1], 3),
    CI_Upper = round(fisher_result$conf.int[2], 3),
    P_Value = fisher_result$p.value
  ))
}

# Print the binary variable summary DataFrame
print(binary_summary_df)

# Create contingency table for 'male'
tbl_male <- table(baseline_data$homeless, baseline_data$female == 0)

# Calculate counts for non-homeless and homeless males
housed_male_count <- tbl_male[1, 2]  # Presence of males in non-homeless group
homeless_male_count <- tbl_male[2, 2]      # Presence of males in homeless group

# Perform Fisher's Exact Test for males
fisher_result_male <- fisher.test(tbl_male)

# Append the results for males to the dataframe
binary_summary_df <- rbind(binary_summary_df, data.frame(
  Variable = "male",
  Housed_Count = housed_male_count,
  Homeless_Count = homeless_male_count,
  Odds_Ratio = round(fisher_result_male$estimate, 3),
  CI_Lower = round(fisher_result_male$conf.int[1], 3),
  CI_Upper = round(fisher_result_male$conf.int[2], 3),
  P_Value = fisher_result_male$p.value
))

# Print the updated binary variable summary DataFrame
print(binary_summary_df)

# Re-ordering the rows in dataframe to have male after female
# Extract the last row (males)
male_row <- binary_summary_df[nrow(binary_summary_df), ]

# Remove the last row from the DataFrame
binary_summary_df <- binary_summary_df[-nrow(binary_summary_df), ]

# Re-insert the male row as the second row
binary_summary_df <- rbind(binary_summary_df[1, ], male_row, binary_summary_df[-1, ])

# Reset row names for clarity
rownames(binary_summary_df) <- NULL

# Save the binary_summary_df DataFrame to a CSV file
write.csv(binary_summary_df, "binary_summary_df.csv", row.names = FALSE)

###############################
# Question 2

# Merge the datasets on the 'id' column
merged_data <- merge(interview_data, pvl_data, by = "id", all = TRUE)

# Creating a year column in the merged data 
merged_data$year <- format(as.Date(merged_data$dt.tst), "%Y")

# Creating a year column in the merged data 
pvl_data$year <- format(as.Date(pvl_data$dt.tst), "%Y")

# Calculate mean, sd, se, and CI
summary_stats <- merged_data %>%
  group_by(year) %>%
  summarise(mean_pl = mean(p, na.rm = TRUE),
            sd_pl = sd(p, na.rm = TRUE),
            n = n(),
            se_pl = sd_pl / sqrt(n),
            ci_upper = mean_pl + 1.96 * se_pl,
            ci_lower = mean_pl - 1.96 * se_pl)

# Plotting
plot_viral <- ggplot(summary_stats, aes(x = year, y = mean_pl)) +
  geom_line(color = "black") +  # Line connecting the mean values
  geom_point(color = "red") +  # Points for the mean values
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, color = "black") +  # Error bars for 95% CI
  labs(x = "Year",
       y = "Plasma Viral Load") +
  theme_minimal()

show(plot_viral)
ggsave(plot_viral, filename = "plot_viral_load.png", path = save_dir, width = 5, height = 4, dpi = 300)


# Calculating mean, standard error, and confidence interval
data_summary <- merged_data %>%
  group_by(year, female, homeless) %>%
  summarise(mean_p = mean(p, na.rm = TRUE),
            se = sd(p, na.rm = TRUE) / sqrt(n()),
            ci_lower = mean_p - 1.96 * se,
            ci_upper = mean_p + 1.96 * se) %>%
  ungroup()

colnames(data_summary)[c(2,3)] <- c("gender","residence")
data_summary <- data_summary %>%
  filter(!is.na(ci_lower) & ci_lower >= 0)

data_summary$year <- as.numeric(data_summary$year)

# Convert 'gender' and 'residence' to factor with more descriptive labels
data_summary$gender <- factor(data_summary$gender, levels = c(0, 1), labels = c("Male", "Female"))
data_summary$residence <- factor(data_summary$residence, levels = c(0, 1), labels = c("Housed", "Homeless"))

# Convert categorical variables to factors
merged_data$female <- as.factor(merged_data$female)
merged_data$homeless <- as.factor(merged_data$homeless)
merged_data$coke <- as.factor(merged_data$coke)
merged_data$crac <- as.factor(merged_data$crac)
merged_data$hero <- as.factor(merged_data$hero)
merged_data$mmt <- as.factor(merged_data$mmt)

# Mixed-effects model, I didn't use this model in the report 
model_mixed <- lmer(p ~ c + female + homeless + coke + crac + hero + mmt + (1|id), 
                    data = merged_data)

summary(model_mixed)

# Histogram of plasam viral load
hist(merged_data$p, breaks = 200, main = "Distribution of Plasma Viral Load", xlab = "Plasma Viral Load")

# Log transforming the viral load
merged_data$plasma_viral_load_log <- log(merged_data$p + 1)

# Histogram of log transformed viral load
hist(merged_data$plasma_viral_load_log, breaks = 30, main = "Distribution of Log-Transformed Plasma Viral Load", xlab = "Log-Transformed Plasma Viral Load")

# Using log transformed data to model
model_log_transformed <- lm(plasma_viral_load_log ~ c + female + homeless + coke + crac + hero + mmt, data = merged_data)
#model_log_transformed2 <- lm(plasma_viral_load_log ~ c + female + homeless, data = merged_data)

summary(model_log_transformed)
vif(model_log_transformed)

AIC(model_log_transformed, model_mixed)

merged_data$year <- as.numeric(as.character(merged_data$year))

# Modelling with GAM (Generalized Additive Model)
gam_model <- gam(p ~ s(year), data = merged_data)
summary(gam_model)

plot(gam_model, pages = 1, scheme = 1, all.terms = TRUE, ylab = "Estimated Plasma Viral Load", xlab = "Year")
####################
# Calculating slopes
# Predict at specific years
# Define the years of interest
years_of_interest <- data.frame(year = c(2002, 2006, 2008, max(merged_data$year,na.rm = T)))

# Predict using the model for those years
predictions <- predict(gam_model, newdata = years_of_interest, se.fit = TRUE)


# Calculate slopes between periods
slope_2002_2006 <- (predictions$fit[2] - predictions$fit[1]) / (2006 - 2002)
slope_2006_2008 <- (predictions$fit[3] - predictions$fit[2]) / (2008 - 2006)
slope_2008_max <- (predictions$fit[4] - predictions$fit[3]) / (2011 - 2008)

#####################

# Make a data frame for predictions
finite_years <- merged_data$year[is.finite(merged_data$year)]

# Now use the finite years to create the sequence
year_range <- data.frame(year = seq(min(finite_years), max(finite_years), length.out = 100))
predictions <- predict(gam_model, newdata = year_range, se = T)

# Add predictions and confidence intervals to the data frame
year_range$fit <- predictions$fit
year_range$upper <- predictions$fit + 1.96 * predictions$se.fit
year_range$lower <- predictions$fit - 1.96 * predictions$se.fit

# Plot with ggplot2
ggplot(year_range, aes(x = year)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "blue", alpha = 0.2) +
  geom_line(aes(y = fit), color = "blue") +
  labs(y = "Plasma Viral Load", x = "Year") +
  theme_minimal()

# plotting data vas fitted values

# Create a new data frame for predictions
gam_predictions <- data.frame(year = year_range$year)

# Get predictions and standard errors
pred <- predict(gam_model, newdata = gam_predictions, se.fit = TRUE)
gam_predictions$fit <- pred$fit
gam_predictions$se <- pred$se.fit

# Calculate confidence intervals
gam_predictions$ci_lower <- gam_predictions$fit - (1.96 * gam_predictions$se)
gam_predictions$ci_upper <- gam_predictions$fit + (1.96 * gam_predictions$se)

# Convert 'year' to numeric if it is not already
gam_predictions$year <- as.numeric(as.character(gam_predictions$year))
summary_stats$year <- as.numeric(as.character(summary_stats$year))

combined_plot <- ggplot() +
  geom_ribbon(data = gam_predictions, aes(x = year, ymin = ci_lower, ymax = ci_upper, fill = "Model CI"), alpha = 0.2) +
  geom_line(data = gam_predictions, aes(x = year, y = fit, color = "Model Fit"), size = 1) +
  geom_point(data = summary_stats, aes(x = year, y = mean_pl, color = "Data Point"), size = 3) +
  geom_errorbar(data = summary_stats, aes(x = year, ymin = ci_lower, ymax = ci_upper, color = "Data CI"), width = 0.2, alpha = 0.7) +
  labs(x = "Year", y = "Plasma Viral Load", color = "Legend", fill = "Legend") +
  theme_minimal() +
  scale_color_manual(name = "", 
                     values = c("Data Point" = "black", "Data CI" = "black", "Model Fit" = "blue"),
                     labels = c("Data CI", "Data Point", "Model Fit")) +
  scale_fill_manual(name = "", 
                    values = c("Model CI" = "blue"),
                    labels = c("Model CI")) +
  guides(fill = guide_legend(override.aes = list(color = NULL)), # Hide color guide for fill
         color = guide_legend(override.aes = list(fill = NULL))) # Hide fill guide for color

# Print the plot
print(combined_plot)

ggsave(combined_plot, filename = "plot_combined.png", path = save_dir, width = 7, height = 4, dpi = 300)

#######################
# Question 3
# Plot with facets for both residence status and gender
ggplot(data_summary, aes(x = year, y = mean_p)) +
  geom_line() +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, fill = "grey50") +
  facet_grid(residence ~ gender) +  # Faceting by both homelessness status and gender
  labs(x = "Year", y = "Mean Plasma Viral Load",
       fill = "Gender") +
  theme_minimal()
gam_model_home <- gam(p ~ s(year) + homeless + s(id, bs="re"), data = merged_data)
summary(gam_model_home)
########################

# Question 4

# Merge the datasets on the 'id' column
merged_data <- merge(interview_data, pvl_data, by = "id", all = TRUE)

# Creating a year column in the merged data 
merged_data$year <- format(as.Date(merged_data$dt.tst), "%Y")
# Ensure dt.tst is in Date format
merged_data$dt.tst <- as.Date(merged_data$dt.tst)

# Calculate start_date for each participant
merged_data2 <- merged_data %>%
  group_by(id) %>%
  mutate(start_date = min(dt.tst)) %>%
  ungroup()

# Assuming `dt.tst` is a date and `c` is some condition
merged_data2 <- merged_data2 %>%
  mutate(event_date = if_else(c < 200, dt.tst, as.Date(NA)))

# Calculating duration days from start of the observation to the event
merged_data2 <- merged_data2 %>%
  mutate(duration_days = as.numeric(difftime(event_date, start_date, units = "days")))

# Defining event_occured: 1 if c < 200 
merged_data2 <- merged_data2 %>%
  mutate(event_occurred = if_else(!is.na(event_date), 1, 0))

summary(merged_data2$duration_days)

dev.off()
# Create a histogram of duration_days
ggplot(merged_data2, aes(x = duration_days)) +
  geom_histogram(binwidth = 1, fill = "blue", color = "black") +
  labs(title = "Histogram of Duration Days",
       x = "Duration in Days",
       y = "Frequency") +
  theme_minimal()

# creating an instance of cleaned_data for plotting and fitting to model
cleaned_data2 <- merged_data2

colnames(cleaned_data2)[c(5,7,10)] <- c("Gender","Residence","Heroin")
cleaned_data2$Residence[cleaned_data2$Residence==1] <- "Homeless"
cleaned_data2$Residence[cleaned_data2$Residence==0] <- "Housed"
cleaned_data2$Gender[cleaned_data2$Gender==1] <- "Female"
cleaned_data2$Gender[cleaned_data2$Gender==0] <- "Male"
cleaned_data2$Heroin[cleaned_data2$Heroin==0] <- "Not-used"
cleaned_data2$Heroin[cleaned_data2$Heroin==1] <- "Used"

surv_object <- Surv(time = merged_data2$duration_days, event = merged_data2$event_occurred)
surv_fit_res <- survfit(surv_object ~ homeless + female + hero, data = merged_data2, na.action = na.exclude)

# Defining surv object from events and duration days
surv_object <- Surv(time = merged_data2$duration_days, event = merged_data2$event_occurred)

# Fitting the surv object, looking at residence status
surv_fit_res <- survfit(surv_object ~ Residence, data = cleaned_data2, na.action = na.exclude)

# Plot with stratification by residence 
res_plot <- ggsurvplot(surv_fit_res, 
                       data = cleaned_data2, 
                       risk.table = TRUE, # Adds a table of numbers at risk at different times
                       pval = TRUE, # Adds p-value of the log-rank test
                       conf.int = TRUE, # Adds confidence intervals
                       palette = c("red", "blue"), # Colors for the groups
                       title = "Residence Survivial Curve",
                       xlab = "Time (Days)", 
                       ylab = "Survival Probability")

print(res_plot)
#ggsave(res_plot$plot, filename = "res_plot.png", path = save_dir, width = 5, height = 4, dpi = 300)

# Fitting the surv object, looking at gender 
surv_fit_gen <- survfit(surv_object ~ Gender, data = cleaned_data2, na.action = na.exclude)

# Plot with stratification by 'treatment_group'
gen_plot <- ggsurvplot(surv_fit_gen, 
                       data = cleaned_data2, 
                       risk.table = TRUE, # Adds a table of numbers at risk at different times
                       pval = TRUE, # Adds p-value of the log-rank test
                       conf.int = TRUE, # Adds confidence intervals
                       palette = c("red", "blue"), # Colors for the groups
                       title = "Gender Survival Curve",
                       xlab = "Time (Days)", 
                       ylab = "Survival Probability")
print(gen_plot)

# Fitting survival curve of Heroin use
surv_fit_drug <- survfit(surv_object ~ Heroin, data = cleaned_data2, na.action = na.exclude)
summary(surv_fit_drug)

# Arrange the plots side by side with one big title
grid.arrange(
  res_plot$plot, 
  gen_plot$plot, 
  ncol = 2, # Set the number of columns to 2 for side-by-side arrangement
  top = "Kaplan-Meier Survival Curves" # Big title at the top
)

combined_plots <- arrangeGrob(res_plot$plot, gen_plot$plot, ncol = 2)

# Use ggsave to save the grob object
ggsave(combined_plots, filename = "survival_plots.png", path = save_dir, width = 9, height = 4, dpi = 300)


# Update the palette to swap colors
drug_plot <- ggsurvplot(surv_fit_drug, 
                        data = cleaned_data2, 
                        risk.table = TRUE, 
                        pval = TRUE, 
                        conf.int = TRUE, 
                        palette = c("red", "blue"), # Swapped colors
                        title = "Kaplan-Meier Survival Curve",
                        xlab = "Time (Days)", 
                        ylab = "Survival Probability")

# Print the plot
print(drug_plot)

# Un-adjusted model for the effect of drug use
unadjusted_model <- coxph(surv_object ~ Heroin, data = cleaned_data2)
summary(unadjusted_model)
cleaned_data2$dt.dob <- as.Date(cleaned_data2$dt.dob, format = "%m/%d/%Y")

# Calculate age at the start of observation
cleaned_data2 <- cleaned_data2 %>%
  mutate(Age_at_start = as.numeric(difftime(start_date, dt.dob, units = "days")) / 365.25)

# Cox proportional hazards model
cox_model_with_age <- coxph(surv_object ~ Heroin + Gender + Residence + Age_at_start, data = cleaned_data2)
summary(cox_model_with_age)

##################################
# Additional Plots

# fitting survival of residence and gender together
surv_fit_res_gen <- survfit(surv_object ~ homeless + female, data = merged_data2, na.action = na.exclude)

res_gen_plot <- ggsurvplot(surv_fit_res_gen, 
                           data = merged_data2, 
                           risk.table = TRUE, # Adds a table of numbers at risk at different times
                           pval = TRUE, # Adds p-value of the log-rank test
                           conf.int = TRUE, # Adds confidence intervals
                           #palette = c("red", "blue"), # Colors for the groups
                           title = "Residence Survivial Curve",
                           xlab = "Time (Days)", 
                           ylab = "Survival Probability")

print(res_gen_plot)
#################################
# fitting survival of residence, gender and Heroin

surv_fit_all <- survfit(surv_object ~ homeless + female + hero, data = merged_data2, na.action = na.exclude)

all_plot <- ggsurvplot(surv_fit_all, 
                       data = merged_data2, 
                       risk.table = TRUE, # Adds a table of numbers at risk at different times
                       pval = TRUE, # Adds p-value of the log-rank test
                       conf.int = TRUE, # Adds confidence intervals
                       #palette = c("red", "blue"), # Colors for the groups
                       title = "Residence Survivial Curve",
                       xlab = "Time (Days)", 
                       ylab = "Survival Probability")

print(all_plot)
