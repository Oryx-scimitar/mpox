library(dplyr)
library(tidyr)
library(EpiCurve)
library(lubridate)
library(ggplot2)
library(incidence)
library(janitor)
library(infer)
library(gtsummary)
install.packages("gtsummary")
# Read in data and convert date columns to date format
mpx <- read.csv("data/MPX_Dataframe_for_R_Nov_2023.csv") %>%
  mutate_at(vars(any_sign_start:meds_other_end), as_date, format="%d/%m/%Y")

# Make dataframe into a "tidy dataframe":
mpx <- mpx %>%
    pivot_longer(!c(id, sex, status, group)) %>%
    mutate(date = case_when(grepl("start", name) ~ "start",
                          TRUE ~ "end")) %>%
    mutate(name = gsub("_start", "", name)) %>%
  mutate(name = gsub("_end", "", name)) %>%
    pivot_wider(names_from = date,
              values_from = value) %>% 
  rename(sign=name)

#create df with variable for fatalities
mortality <- mpx %>% 
  group_by(id) %>% 
  filter(sign == "death") %>% 
  mutate(fatal = if_else(is.na(start)==FALSE, "y", "n"))
#use ids of animals that died to create 'fatal' variable within main df
mpx <- mpx %>% 
  group_by(id) %>% 
  mutate(fatal = if_else(id == "Lola" | id == "Janet", "y", "n"))

##To find minimum and maximum likely date of point source exposure:
##extrapolate data from published human incubation periods (Nolen et al. 2016)
##then compare to onset date in the index case
##to manually calculate dates for use in epicurve 
#Determine most likely period of exposure
#Latest possible date of exposure 
#(date of onset of index case minus lower bound IQ range = 5d)
as.Date("2016-08-14") - 5
# = "2016-08-09"
#Earliest possible date of exposure
#(date of onset of index case minus upper bound IQ range = 13d)
as.Date("2016-08-14") - 13
# = "2016-08-01"
#Time between latest exposure date and onset of last case
as.Date("2016-09-05") - as.Date("2016-08-09")
# n = 27 days
#Latest possible exposure date plus maximum incubation period
#(postulate that cases occurring later than this due to chimp-to-chimp transmission)
as.Date("2016-08-14") - 5 + 13
# = "2016-08-22"
#Number of cases that occurred after latest possible exposure plus maximum incubation period
mpx %>% 
  filter(sign == "any_sign", start > as.Date("2016-08-22"))
#n=13
#Time between separation of groups and onset of first case in subgroup B
as.Date("2016-08-30") - as.Date("2016-08-17")
# n = 13 days
#Determine time between onset of first and last case
as.Date("2016-09-05") - as.Date("2016-08-14")
# n = 22 days

#EpiCurve with ggplot
epicurve <- mpx %>% 
  filter(sign == "any_sign", is.na(start)==FALSE) %>%
  ggplot +
  geom_histogram(mapping = aes(x = start, group = group, fill = group), 
                 binwidth = 1, 
                 color = "darkblue") +
  scale_x_date(date_breaks = "week", date_minor_breaks = "day",
               date_labels = "%d %b", expand = c(0,0),
               limits = as.Date(c("2016-07-30", NA))) +
  scale_fill_discrete(labels = c("Subgroup A", "Subgroup B"),  breaks = c("Tommy", "Bolly")) +
  geom_vline(aes(xintercept = as.Date(c("2016-08-09, NA"))),
             linetype = "dashed", color = "black") +
  geom_text (aes(x=as.Date(c("2016-08-09, NA"))), label = "B", 
             y = 2.6, angle = 0, vjust = -0, size = 3.5, nudge_x = 1) +
  geom_vline(aes(xintercept = as.Date(c("2016-08-01, NA"))),
             linetype = "dashed", color = "black") +
  geom_text (aes(x=as.Date(c("2016-08-01, NA"))), label = "A", 
             y = 2.6, angle = 0, vjust = -0, size = 3.5, nudge_x = 1) +
  geom_vline(aes(xintercept = as.Date(c("2016-08-22, NA"))),
             linetype = "dashed", color = "black") +
  geom_text (aes(x=as.Date(c("2016-08-22, NA"))), 
             label = "C", 
             y = 2.6, angle = 0, vjust = -0, size = 3.5, nudge_x = 1) +
  geom_vline(aes(xintercept = as.Date(c("2016-08-17, NA"))),
             linetype = "dashed", color = "red") +
  labs (title = "", x = "Date of onset", y = "Incident cases",
        fill = "",
        caption = "") +
  geom_text (aes(x=as.Date(c("2016-08-14, NA"))), 
            label = "*", 
            size = 7,
            y = 1, angle = 0, vjust = 0.4, size = 3.5) +
  geom_text (aes(x=as.Date(c("2016-08-16, NA"))), 
             label = "**", 
             size = 7,
             y = 2, angle = 0, vjust = 0.4, size = 3.5) +
  geom_text (aes(x=as.Date(c("2016-08-22, NA"))), 
             label = "+", 
             size = 5,
             y = 2, angle = 0, vjust = -0, hjust = 0.8, size = 3.5) +
  geom_text (aes(x=as.Date(c("2016-09-04, NA"))), 
             label = "++", 
             size = 5,
             y = 3, angle = 0, vjust = -0, size = 3.5) 
print (epicurve)

#Save and print epicurve
ggsave(filename = "results/epicurve_5-13d_Last_date_subgroup_direct_contact_rpt.png",
       plot = epicurve)

##DESCRIPTIVE STATS
#calculate no. of cases and attack rate
attack_rate <- mpx %>% 
  mutate(case = if_else(status == "c" | status == "p",
                        "case", "no_signs")) %>% 
  group_by(case) %>% 
  filter(sign == "any_sign") %>%  
  summarise(n = n()) %>% 
  mutate(percent = n / sum (n)*100)

#calculate no. and % confirmed, probable and negative
status <- mpx %>% 
  group_by(status) %>% 
  filter(sign == "any_sign") %>%  
  summarise(n = n()) %>% 
  mutate(percent = n / sum (n)*100)

#calculate no. and percentages affected by each sign
#N.B. denominator is total cases (i.e. excludes negative animals)
count <- mpx %>% 
  group_by(sign) %>% 
  filter(is.na(start)==FALSE) %>% 
  count %>% 
  mutate(percent = n / 20*100)

#calculate overall range and median
#(after converting duration <1 day to 1 day)
summary <- mpx %>% 
  group_by(sign) %>% 
  mutate(duration = ifelse(end - start <1, 1, end - start)) %>% 
  summarise(range = range(duration,na.rm = T),
            median = median(duration, na.rm = T))

#DESCRIPTIVE STATS: SURVIVORS VS FATALITIES
#calculate no. and percentages affected by each sign
#N.B. denominator is total cases (i.e. excludes negative animals)
count_by_death <- mpx %>% 
  group_by(sign, fatal) %>% 
  filter(is.na(start)==FALSE) %>% 
  count %>% 
  mutate(total_cases_by_death = if_else(fatal == "n", 18, 2)) %>% 
  mutate(percent = n/total_cases_by_death)

summary_by_death <- mpx %>% 
  group_by(sign, fatal) %>% 
  mutate(duration = ifelse(end - start <1, 1, end - start)) %>% 
  summarise(range = range(duration,na.rm = T),
            median = median(duration, na.rm = T))

  
#DESCRIPTIVE STATS: BY SUBGROUP
#N.B. denominator is total cases within relevant subgroup 
count_by_group <- mpx %>% 
  group_by(group, sign) %>% 
  filter(is.na(start)==FALSE) %>% 
  count %>% 
  mutate(total_cases_by_group = if_else(group == "Bolly", 7, 13)) %>% 
  mutate(percent = n/total_cases_by_group)

summary_by_group <- mpx %>% 
  group_by(sign, group) %>% 
  mutate(duration = ifelse(end - start <1, 1, end - start)) %>% 
  summarise(range = range(duration,na.rm = T),
            median = median(duration, na.rm = T))

#create dataset with signs listed as binary variables i.e. present vs absent
#to enable use for testing difference of proportions for presence/absence of each sign
compare <- mpx %>% 
  #if start = NA, recode "absent", if start = a date, recode "present"
  mutate(start = if_else((is.na(start) == T), "no", "yes")) %>% 
  rename(affected = start) %>% 
  #dplyr:: glimpse(diff_by_group) 
  #most variables are characters
  #convert character variables to factors
  mutate_if(is.character, as.factor) %>% 
  group_by(id, group, sign, affected, fatal, status) %>% 
  summarise
dplyr:: glimpse(compare) 

#create wider df from compare
#to enable comparison of individual clinical signs by e.g. group, survival 
compare_wide <- compare %>%
  pivot_wider(names_from = sign,
              values_from = affected) 

#for cases, calculate difference in signs by survival
compare_wide %>% 
  #filter to include only cases (remove negatives)
  filter(status != "n") %>% 
  tbl_summary(
    by = fatal,
    statistic = list(
      all_categorical() ~ "{n} / {N} ({p}%)"),
    include = analgesia_1 : ulceration
  ) %>% 
  add_p() %>% 
  add_ci

#calculate attack rate by group using tbl_summary
compare_wide %>% 
  tbl_summary(
    by = group,
    statistic = list(
      all_categorical() ~ "{n} / {N} ({p}%)"),
    include = any_sign
  ) %>% 
  add_p() %>% 
  add_ci
  
#for cases, calculate difference in signs by group
compare_wide %>% 
  #filter to include only cases (remove negatives)
  filter(status != "n") %>% 
  tbl_summary(
    by = group,
    statistic = list(
      all_categorical() ~ "{n} / {N} ({p}%)"),
    include = analgesia_1 : ulceration
  ) %>% 
  add_p() %>% 
  add_ci


#create variable for prodrome
prodrome <- mpx %>% 
  filter(status != "n") %>% 
  select(id, sign, start, status, fatal, group) %>% 
  filter(sign == "any_sign" | sign == "exanthema") %>% 
  pivot_wider(names_from = sign, values_from = start) %>% 
  mutate(duration = exanthema - any_sign) %>% 
  #as.numeric(duration)
  mutate(affected = if_else(duration == 0, "no", "yes")) %>% 
  #convert character variables to factors
  mutate_if(is.character, as.factor) %>% 
  group_by(id, status, fatal, group, duration, affected) %>% 
  summarise
#descriptive stats for prodrome
summary_prodrome <- prodrome %>% 
  ungroup %>% 
  filter(is.finite(duration == T), duration != 0) %>% 
  reframe(
    n=n(),
    range = range(duration, na.rm = T),
    median = median(duration, na.rm = T))
#by mortality
summary_prodrome_by_death <- prodrome %>% 
  ungroup %>% 
  filter(is.finite(duration == T), duration != 0) %>% 
  group_by(fatal) %>% 
  reframe(
    n=n(),
    range = range(duration, na.rm = T),
    median = median(duration, na.rm = T))
#by subgroup
summary_prodrome_by_group <- prodrome %>% 
  filter(is.finite(duration == T), duration != 0) %>% 
  group_by(group) %>% 
  reframe(
    n=n(),
    range = range(duration, na.rm = T),
    median = median(duration, na.rm = T))
#compare presence of prodrome between subgroups
prodrome %>% 
  #consider animals where exanthema was not observed as having no prodrome
  mutate(affected = if_else(affected == "no", "no", "yes", missing = "no")) %>% 
  tbl_summary(
  by = group,
  statistic = list(
    all_categorical() ~ "{n} / {N} ({p}%)"),
  include = affected
) %>% 
  add_p() %>% 
  add_ci
#compare presence of prodrome between survivors and fatalities
prodrome %>% 
  tbl_summary(
    by = fatal,
    statistic = list(
      all_categorical() ~ "{n} / {N} ({p}%)"),
    include = affected
  ) %>% 
  add_p() %>% 
  add_ci

#Create grouped clinical sign variables for signs affecting a particular body system or with
#a probable shared aetiology
edema <- mpx %>% 
  filter(status != "n") %>% 
  filter(sign == "edema_facial" | sign == "edema_laryngeal") %>% 
  group_by(id, fatal, group) %>% 
  reframe(min = min(start, na.rm = T),
          max = max(end, na.rm = T)) %>%
  mutate(edema_any = if_else((is.finite(min) == T),
                             "present", "absent"),
         duration = ifelse((is.finite(min) == T), 
                           ifelse(max - min <1, 1, max - min),
                           NA))

resp <- mpx %>% 
  filter(status != "n") %>% 
  filter(sign == "cough" | sign == "coryza" | sign == "dyspnea") %>% 
  group_by(id, fatal, group) %>% 
  reframe(min = min(start, na.rm = T),
          max = max(end, na.rm = T)) %>%
  mutate(resp_any = if_else((is.finite(min) == T),
                            "present", "absent"),
         duration = ifelse((is.finite(min) == T), 
                           ifelse(max - min <1, 1, max - min),
                           NA))

skin <- mpx %>% 
  filter(status != "n") %>% 
  filter(sign == "exanthema" | sign == "eschar" | sign == "ulceration") %>% 
  group_by(id, fatal, group) %>% 
  reframe(min = min(start, na.rm = T),
          max = max(end, na.rm = T)) %>%
  mutate(skin_any = if_else((is.finite(min) == T),
                            "present", "absent"),
         duration = ifelse((is.finite(min) == T), 
                           ifelse(max - min <1, 1, max - min),
                           NA))

#descriptive stats for grouped clinical sign variables
edema %>% 
  group_by(edema_any) %>% 
  reframe(n = n()) %>% 
  mutate(percent = n / sum (n)*100) 

edema %>% 
  group_by(fatal) %>% 
  reframe(
    range = range(duration,na.rm = T),
    median = median(duration, na.rm = T))

resp %>% 
  group_by(resp_any) %>% 
  reframe(n = n()) %>% 
  mutate(percent = n / sum (n)*100) 

resp %>% 
  group_by(fatal) %>% 
  reframe(
    range = range(duration,na.rm = T),
    median = median(duration, na.rm = T))

skin %>% 
  group_by(skin_any) %>% 
  reframe(n = n()) %>% 
  mutate(percent = n / sum (n)*100) 

skin %>% 
  group_by(fatal) %>% 
  reframe(
    range = range(duration,na.rm = T),
    median = median(duration, na.rm = T))

#DESCRIPTIVE STATS FOR GROUPED CLINICAL SIGN VARIABLES: SURVIVORS VS FATALITIES
#calculate no. and percentages affected by each sign
#N.B. denominator is total cases (i.e. excludes negative animals)
edema_count_by_death <- edema %>% 
  group_by(fatal) %>% 
  filter(edema_any == "present")  %>% 
  count %>% 
  mutate(total_edema_by_death = if_else(fatal == "n", 18, 2)) %>% 
  mutate(percent = n/total_edema_by_death)

edema_summary_by_death <- edema %>% 
  group_by(fatal) %>% 
  summarise(range = range(duration,na.rm = T),
            median = median(duration, na.rm = T))


skin_count_by_death <- skin %>% 
  group_by(fatal) %>% 
  filter(skin_any == "present") %>% 
  count %>% 
  mutate(total_skin_by_death = if_else(fatal == "n", 18, 2)) %>% 
  mutate(percent = n/total_skin_by_death)

skin_summary_by_death <- skin %>% 
  group_by(fatal) %>% 
  summarise(range = range(duration,na.rm = T),
            median = median(duration, na.rm = T))

resp_count_by_death <- resp %>% 
  group_by(fatal) %>% 
  filter(resp_any == "present")  %>% 
  count %>% 
  mutate(total_resp_by_death = if_else(fatal == "n", 18, 2)) %>% 
  mutate(percent = n/total_resp_by_death)

resp_summary_by_death <- resp %>% 
  group_by(fatal) %>% 
  summarise(range = range(duration,na.rm = T),
            median = median(duration, na.rm = T))

#DESCRIPTIVE STATS FOR GROUPED CLINICAL SIGN VARIABLES: BY SUBGROUP
#calculate no. and percentages affected by each sign
#N.B. denominator is total cases (i.e. excludes negative animals)
edema_count_by_group <- edema %>% 
  group_by(group) %>% 
  filter(edema_any == "present")  %>% 
  count %>% 
  mutate(total_edema_by_group = if_else(group == "Bolly", 7, 13)) %>% 
  mutate(percent = n/total_edema_by_group)

edema_summary_by_group <- edema %>% 
  group_by(group) %>% 
  summarise(range = range(duration,na.rm = T),
            median = median(duration, na.rm = T))
edema %>% 
  tbl_summary(
    by = group,
    statistic = list(
      all_categorical() ~ "{n} / {N} ({p}%)"),
    include = edema_any
  ) %>% 
  add_p() %>% 
  add_ci()

skin_count_by_group <- skin %>% 
  group_by(group) %>% 
  filter(skin_any == "present") %>% 
  count %>% 
  mutate(total_skin_by_group = if_else(group == "Bolly", 7, 13)) %>% 
  mutate(percent = n/total_skin_by_group)

skin_summary_by_group <- skin %>% 
  group_by(group) %>% 
  summarise(range = range(duration,na.rm = T),
            median = median(duration, na.rm = T))

skin %>% 
  tbl_summary(
    by = group,
    statistic = list(
      all_categorical() ~ "{n} / {N} ({p}%)"),
    include = skin_any
  ) %>% 
  add_p() %>% 
  add_ci()

resp_count_by_group <- resp %>% 
  group_by(group) %>% 
  filter(resp_any == "present")  %>% 
  count %>% 
  mutate(total_resp_by_group= if_else(group == "Bolly", 7, 13)) %>% 
  mutate(percent = n/total_resp_by_group)

resp_summary_by_group <- resp %>% 
  group_by(group) %>% 
  summarise(range = range(duration,na.rm = T),
            median = median(duration, na.rm = T))

resp %>% 
  tbl_summary(
    by = group,
    statistic = list(
      all_categorical() ~ "{n} / {N} ({p}%)"),
    include = resp_any
  ) %>% 
  add_p() %>% 
  add_ci()


#DESCRIPTIVE STATS: ASSOCIATIONS BETWEEN CLINICAL SIGNS
#look for association between grouped resp signs and death
resp %>% 
  tbl_summary(
    by = fatal,
    statistic = list(
      all_categorical() ~ "{n} / {N} ({p}%)"),
    include = resp_any
  ) %>% 
  add_p() %>% 
  add_ci()

#look for association between dyspnea and other clinical signs
compare_wide %>% 
  #filter to include only cases (remove negatives)
  filter(status != "n") %>% 
  tbl_summary(
    by = dyspnea,
    statistic = list(
      all_categorical() ~ "{n} / {N} ({p}%)"),
    include = analgesia_1 : ulceration
  ) %>% 
  add_p() %>% 
  add_ci

#look for association between peri-laryngeal edema and other clinical signs
compare_wide %>% 
  #filter to include only cases (remove negatives)
  filter(status != "n") %>% 
  tbl_summary(
    by = edema_laryngeal,
    statistic = list(
      all_categorical() ~ "{n} / {N} ({p}%)"),
    include = analgesia_1 : ulceration
  ) %>% 
  add_p() %>% 
  add_ci

#look for association between steroid use and survival in animals with 
#peri-laryngeal edema 
compare_wide %>% 
  #filter to include only cases (remove negatives)
  filter(status != "n") %>% 
  filter(edema_laryngeal == "yes") %>%
  tbl_summary(
    by = fatal,
    statistic = list(
      all_categorical() ~ "{n} / {N} ({p}%)"),
    include = analgesia_1 : ulceration
  ) %>% 
  add_p() %>% 
  add_ci
