library(tidyverse)
library(gtsummary)
library(gt)
library(ggalluvial)
library(usmap)

# mortality data
mset <- read_csv('data/pca_mortality_2011_2017.csv', skip = 2, na = c('', 'NA', '^'))
colnames(mset) <- c('fips', 
                    'all_race_rate_new', 
                    'all_race_count_new', 
                    'all_race_pop_new', 
                    'white_rate_new', 
                    'white_count_new', 
                    'white_pop_new', 
                    'black_rate_new', 
                    'black_count_new', 
                    'black_pop_new', 
                    'other_ai_rate_new', 
                    'other_ai_count_new', 
                    'other_ai_pop_new', 
                    'other_uns_rate_new', 
                    'other_uns_count_new', 
                    'other_uns_pop_new')
mset <- select(mset, 
               'fips', 
               'white_rate_new', 
               'white_count_new', 
               'white_pop_new', 
               'black_rate_new', 
               'black_count_new', 
               'black_pop_new')
mset$fips <- str_extract(mset$fips, '\\d{5}')
mset <- na.omit(mset)
mset$mort_ratio_new <- mset$black_rate_new / mset$white_rate_new

#######################
# mortality map
#######################
mmap <- plot_usmap(
  linewidth = 0.04,        
  color = "gray60",          
  regions = "counties",     
  data = mset,
  values = "mort_ratio_new",
  exclude = c("AK", "HI")   
) +
  scale_fill_gradient(
    name = "Black-to-White mortality ratio", 
    high = "darkblue",
    low = "white",
    limits = c(0, 6),  
    guide = guide_colorbar(
      title.position = "top",           
      title.hjust = 0.5,                
      barwidth = 15,                    
      barheight = 0.5                  
    )
  ) +
  labs(title = "County-Level Prostate Cancer Mortality Disparities (2011-2017)") +
  theme_void() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"), 
    text = element_text(size = 14)
  )

ggsave('results/mortality_map.png',
       plot = mmap,
       width = 10, height = 10, dpi = 500)

# screening data
bset <- read_csv('data/Behavioral_Risk_Factors__Selected_Metropolitan_Area_Risk_Trends__SMART__MMSA_Prevalence_Data__2011_to_Present__20240711.csv')
bset <- filter(bset, Year != 2020)
bset <- filter(bset, Year != 2018)

bset <- filter(bset, Response != "No")
bset <- aggregate(cbind(Data_value, Sample_Size) ~ Locationabbr, 
                  data = bset, 
                  FUN = function(x) mean(x, na.rm = TRUE))
bset$Data_value <- round(bset$Data_value, 3)
bset$Sample_Size <- round(bset$Sample_Size, 0)
colnames(bset) <- c('cbsa', 'screening_frequency', 'screening_sample_size')

# crosswalk
cset <- read.csv('data/cbsatocountycrosswalk.csv')
cset$fipscounty <- str_pad(cset$fipscounty, width = 5, pad = "0")
cset <- select(cset, cbsa, fipscounty)
cset <- rename(cset, fips = fipscounty)
cset <- na.omit(cset)

#######################
# screening map
#######################
tset <- inner_join(bset, cset, by = 'cbsa')
tset$freq_10 <- tset$screening_frequency / 10

smap <- plot_usmap(
  linewidth = 0.04,        
  color = "gray60",            
  regions = "counties",     
  data = tset,
  values = "freq_10",
  exclude = c("AK", "HI")   # exclude Alaska and Hawaii
) +
  scale_fill_gradient(
    name = "Screening rate (by 10%)", 
    high = "red",
    low = "yellow",
    limits = c(2, 7),  
    guide = guide_colorbar(
      title.position = "top",          
      title.hjust = 0.5,                
      barwidth = 15,                    
      barheight = 0.5                   
    )
  ) +
  labs(title = "County-Level PSA Screening Frequencies (2012, 2014, 2016)") +
  theme_void() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 17, face = "bold"), 
    text = element_text(size = 14)
  )

ggsave('results/screening_map.png',
       plot = smap,
       width = 10, height = 10, dpi = 500)

# svi data
vset1 <- read_csv('data/SVI_2014_US_county.csv')
vset1$fips <- str_pad(as.character(vset1$FIPS), width = 5, pad = '0')
vset1 <- select(vset1, 'fips', 'EP_POV', 'EP_UNEMP', 'EP_NOHSDP', 'EP_UNINSUR')

vset2 <- read_csv('data/SVI_2016_US_county.csv')
vset2$fips <- str_pad(as.character(vset2$FIPS), width = 5, pad = '0')
vset2 <- select(vset2, 'fips', 'EP_POV', 'EP_UNEMP', 'EP_NOHSDP', 'EP_UNINSUR')

vset3 <- read_csv('data/SVI_2018_US_county.csv')
vset3$fips <- str_pad(as.character(vset3$FIPS), width = 5, pad = '0')
vset3 <- select(vset3, 'fips', 'EP_POV', 'EP_UNEMP', 'EP_NOHSDP', 'EP_UNINSUR')
vset3 <- vset3 |> filter(fips != '35039')                                           # removed this county due to abnormal data in EP_POV and EP_UNEMP (both had value -999.0)

vset <- bind_rows(vset1, vset2, vset3)

vset <- group_by(vset, fips)
vset <- summarize(
  vset,
  avg_poverty = mean(EP_POV, na.rm = TRUE),
  avg_unemployment = mean(EP_UNEMP, na.rm = TRUE),
  avg_nohsdp = mean(EP_NOHSDP, na.rm = TRUE),
  avg_uninsured = mean(EP_UNINSUR, na.rm = TRUE),
  .groups = "drop"
)

# merge 
dset <- inner_join(bset, cset, by = 'cbsa')
dset <- inner_join(dset, mset, by = 'fips')
dset <- inner_join(dset, vset, by = 'fips')

# calculate variance ratios (delta method)
dset$white_variance_new <- dset$white_count_new / dset$white_pop_new ^ 2 * 1e5 ^ 2
dset$black_variance_new <- dset$black_count_new / dset$black_pop_new ^ 2 * 1e5 ^ 2
dset$variance_ratio_new <- 1/dset$black_rate_new ^ 2 * dset$white_variance_new + 1/dset$white_rate_new ^ 2 * dset$black_variance_new
dset$screening_frequency <- dset$screening_frequency / 10

# regression
model <- lm(mort_ratio_new ~ screening_frequency + avg_poverty + avg_uninsured + avg_nohsdp, 
            data = dset, 
            weights = 1 / variance_ratio_new)

t1 <- tbl_regression(model, label = list(screening_frequency ~ 'Screening Frequency', 
                                         avg_poverty ~ 'Poverty Rate',
                                         avg_uninsured ~ 'No Health Insurance Rate',
                                         avg_nohsdp ~ 'No High School Diploma Rate'))

t1
as_gt(t1) |> gtsave('results/regression_tbl_psa_era.html')

# scatter plot new
scatter_plot_1 <- ggplot(dset, aes(x = screening_frequency, y = mort_ratio_new, size = 1/variance_ratio_new)) +
  geom_point(color = "black", alpha = 0.57) + 
  geom_smooth(method = "lm", aes(weight = 1 / variance_ratio_new), se = FALSE, linewidth = 1.25) +
  labs(
    x = "Screening Frequency (rate / 10)",
    y = "Mortality Rate Ratio (B/W)",
    title = "County-Level Screening Frequencies (2012, 2014, 2016) vs. Mortality Rate Ratios (2011-2017)"
  ) +
  theme_classic() +  
  theme(
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5, margin = margin(b = 20)),   
    axis.title = element_text(size = 12, face = "bold"),                
    axis.text = element_text(size = 12),                                 
    #plot.margin = margin(t = 22, r = 10, b = 30, l = 14),
    legend.position = 'none'
  ) 

scatter_plot_1

ggsave('results/scatter_plot_psa_era.png',
       plot = scatter_plot_1,
       width = 12, height = 10, dpi = 300)

#######################
# pre psa era analysis
#######################
pset <- read_csv('data/pca_mortality_1984_1989.csv', skip = 2, na = c('', 'NA', '^'))
colnames(pset) <- c('fips', 
                    'all_race_rate_old', 
                    'all_race_count_old', 
                    'all_race_pop_old', 
                    'white_rate_old', 
                    'white_count_old', 
                    'white_pop_old', 
                    'black_rate_old', 
                    'black_count_old', 
                    'black_pop_old', 
                    'other_ai_rate_old', 
                    'other_ai_count_old', 
                    'other_ai_pop_old', 
                    'other_uns_rate_old', 
                    'other_uns_count_old', 
                    'other_uns_pop_old')
pset <- select(pset, 
               'fips', 
               'white_rate_old', 
               'white_count_old', 
               'white_pop_old', 
               'black_rate_old', 
               'black_count_old', 
               'black_pop_old')
pset$fips <- str_extract(pset$fips, '\\d{5}')
pset <- na.omit(pset)
pset$mort_ratio_old <- pset$black_rate_old / pset$white_rate_old

# merge with modern screening data
qset <- inner_join(bset, cset, by = 'cbsa')
qset <- inner_join(qset, pset, by = 'fips')

# calculate variance ratios (delta method)
qset$white_variance_old <- qset$white_count_old / qset$white_pop_old ^ 2 * 1e5 ^ 2
qset$black_variance_old <- qset$black_count_old / qset$black_pop_old ^ 2 * 1e5 ^ 2
qset$variance_ratio_old <- 1/qset$black_rate_old ^ 2 * qset$white_variance_old + 1/qset$white_rate_old ^ 2 * qset$black_variance_old
qset$screening_frequency <- qset$screening_frequency / 10

# regression
model2 <- lm(mort_ratio_old ~ screening_frequency, 
             data = qset, 
             weights = 1 / variance_ratio_old)

t2 <- tbl_regression(model2, label = list(screening_frequency ~ 'Screening Frequency'))

t2
as_gt(t2) |> gtsave('results/regression_tbl_pre_psa_era.html')

# scatter plot old
scatter_plot_2 <- ggplot(qset, aes(x = screening_frequency, y = mort_ratio_old, size = 1/variance_ratio_old)) +
  geom_point(color = "darkblue", alpha = 0.57) + 
  geom_smooth(method = "lm", aes(weight = 1 / variance_ratio_old), se = FALSE, linewidth = 1.25, color = "red") +
  labs(
    x = "Screening Frequency (rate / 10)",
    y = "Mortality Rate Ratio (B/W)",
    title = "County-Level Screening Frequencies (2012, 2014, 2016) vs. Mortality Rate Ratios (1984-1989)"
  ) +
  theme_classic() +  
  theme(
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5, margin = margin(b = 20)),   
    axis.title = element_text(size = 12, face = "bold"),                
    axis.text = element_text(size = 12),                                 
    #plot.margin = margin(t = 22, r = 10, b = 30, l = 14),
    legend.position = 'none'
  ) 

scatter_plot_2

ggsave('results/scatter_plot_pre_psa_era.png',
       plot = scatter_plot_2,
       width = 12, height = 10, dpi = 300)

# alluvial plot
scr <- read_csv('data/Behavioral_Risk_Factors__Selected_Metropolitan_Area_Risk_Trends__SMART__MMSA_Prevalence_Data__2011_to_Present__20240711.csv')

scr <- filter(scr, Year != 2020, Response != "No")
scr <- select(scr, Year, Locationabbr, Data_value, Sample_Size)
scr <- scr |> group_by(Year) 
scr <- scr |> mutate(quartile = ntile(Data_value, 2)) 
scr <- scr |> ungroup()
scr <- filter(scr, Year %in% c(2012, 2016))

baseline_2012 <- scr |> filter(Year == 2012) 
baseline_2012 <- baseline_2012 |> select(Locationabbr, quartile) 
baseline_2012 <- baseline_2012 |> rename(baseline = quartile)

scr <- left_join(scr, baseline_2012, by = "Locationabbr")
scr <- scr |> group_by(Locationabbr) 
scr <- scr |> filter(n_distinct(Year) == 2) 
scr <- scr |> ungroup()
scr <- select(scr, Year, Locationabbr, quartile, baseline)
scr <- na.omit(scr)

gg_theme <- function(...){
  theme_set(theme_classic())
  theme_update(axis.text=element_text(color='black'),
               axis.ticks.length=unit(0.2, 'cm'),
               strip.text=element_text(color='black', size=9),
               strip.background=element_rect(color=NA, fill=NA))
  theme_update(...)
}

scr_alluvial_plot <- function(scr){
  
  gg_theme(legend.position='none',
           axis.line.y=element_blank(),
           axis.ticks.y=element_blank(),
           axis.text=element_text(color='black', size=14),  # Increase label size
           axis.text.y=element_blank())  # Remove y-axis text
  
  scr$Year <- as.factor(scr$Year)
  scr$quartile <- as.factor(scr$quartile)
  gg <- ggplot(scr, aes(x=Year,
                        stratum=quartile,
                        alluvium=Locationabbr,
                        fill=quartile,  
                        label=quartile))
  gg <- gg + geom_alluvium(aes(fill = quartile))
  gg <- gg + geom_stratum(aes(fill = quartile), alpha = 0, size = 1, width = 0.25)
  gg <- gg + geom_text(stat='stratum', aes(label=ifelse(quartile == 1, 'Top 50%', 
                                                        ifelse(quartile == 2, 'Bottom 50%', quartile))), 
                       size=6)  
  gg <- gg + scale_x_discrete(name='', expand=c(0, 0))
  gg <- gg + scale_y_continuous(expand=c(0, 0), breaks=NULL, labels=NULL)
  gg <- gg + scale_fill_manual(values = c('#E69F00', '#480082'))  
  
  print(gg)
}

alluvial <- scr_alluvial_plot(scr)

ggsave('results/alluvial.png',
       plot = alluvial,
       width = 8, height = 5, dpi = 300)


