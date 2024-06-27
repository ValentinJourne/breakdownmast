library(gratia)
library(tidyverse)
library(here)
library(betareg)
library(broom.mixed)
library(broom)
library(mgcv)
#se correlation
correlation.spearman.se <- function (cor, n) 
{
  se.cor <- sqrt((1 - cor^2)^2 * (1 + cor^2/2)/(n - 3))
  return(se.cor)
}

#to use for beta regression, values ned to be 0<y<1
y.transf.betareg <- function(y){
  n.obs <- sum(!is.na(y))
  (y * (n.obs - 1) + 0.5) / n.obs
}

#detach package, I used to check issue coming from function redundant names 
detachAllPackages <- function() {
  
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
  
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  
  package.list <- setdiff(package.list,basic.packages)
  
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
  
}
#make calendar to associate DOY 
goingbackpastdayscalendar <-function(refday = 274, 
                                     lastdays = 1095, 
                                     yearback = 3){
  monthstart = c('-01-01')
  DATE = seq(as.Date(paste0("1946", monthstart)), as.Date(paste0("1953", monthstart)), by="days")
  MONTH =  format(as.Date(DATE, format="%Y-%m-%d"),"%m") %>% as.numeric()
  MONTHab = month.abb[MONTH]
  YEAR =  format(as.Date(DATE, format="%Y-%m-%d"),"%Y") %>% as.numeric()
  
  DOY = yday(DATE)
  dfata = data.frame(DATE,YEAR,MONTHab, DOY)
  yearperiod = 1946:1953
  sizevec = length(unique(YEAR))-yearback
  refday = refday
  vectotemp = NULL
  
  for(k in 1:sizevec){
    #year +3, becusse of the 36 month analysis 
    yearsref = yearperiod[k]
    yearrefminusOne <- yearsref-yearback #PREVIOUS = 3
    
    tt <- dfata %>% 
      filter(YEAR <= yearsref & YEAR >= yearrefminusOne) %>% 
      mutate(referenceFin = ifelse(YEAR == yearsref & DOY == refday, 1,
                                   ifelse(YEAR == yearsref & DOY > refday, NA, 0))) %>% 
      filter(!is.na(referenceFin)) %>% 
      as.data.frame()
    #create sequence going back 365 month before 
    seqDays <- seq(1,nrow(tt),1)
    newsequance <- rep(seqDays)
    
    ttup <- tt %>% 
      mutate(days.reversed = rev(newsequance))%>% 
      filter(days.reversed< lastdays )
    ttupfin = ttup %>%
      arrange(days.reversed)  %>% 
      mutate(YEAR = max(YEAR))  
    
    vectotemp <- rbind(vectotemp, ttupfin) 
  }
  
  return(vectotemp)
}

addphylogeny.family.genus = function(data = data){
  library(taxizedb)
  ids <- name2taxid(unique(data$Species), out_type="summary", db = "ncbi")
  genus <- taxa_at(ids$id,
                   rank = c('genus'), missing = "lower", db="ncbi")
  genus = do.call(rbind.data.frame, genus)
  family <- taxa_at(ids$id,
                    rank = c('family'), missing = "lower", db="ncbi")
  family = do.call(rbind.data.frame, family)
  
  taxon.list.selec = cbind(ids, genus$name, family$name) %>% 
    rename(Species = 1, genus = 3, family = 4) 
  
  phylolist = data %>% 
    dplyr::select(Species) %>% 
    distinct()%>% 
    left_join(taxon.list.selec) %>% 
    mutate(genus = ifelse(str_detect(Species, 'spp.'), 
                          str_remove(Species, " spp."), 
                          genus)) %>% 
    mutate(family = case_when(
      Species == "Nothofagus spp." ~ "Nothofagaceae",
      Species == "Quercus spp." ~ "Fagaceae",
      Species == "Abies spp." ~ "Pinaceae",
      Species == "Pinus spp." ~ "Pinaceae",
      TRUE ~ family
    ))
  return(phylolist)
}

# Function to process each year to get past moving climate 
reformat.climate.backtothepast <- function(yearsref = 2000, 
                                           climate = climate, 
                                           yearneed = 2, 
                                           refday = 274, 
                                           lastdays = 1095, 
                                           rollwin = 7, 
                                           variablemoving = 'temperature.degree') {
  # Print parameter values - to see if no shit 
  print(paste("rollwin - size window:", rollwin))
  print(paste("refday - reference day:", refday))
  print(paste("lastdays - last day:", lastdays))
  print(paste("variablemoving - variable to move:", variablemoving))
  print(paste("yearneed - adjust to maturation fruit time :", yearneed, 'years'))
  
  if (!variablemoving %in% names(climate)) {
    warning(paste("Warning: Column", variablemoving, "not found in the dataset."))
    return(NULL)
  }
  
  yearrefminusOne <- yearsref - yearneed
  tt <- climate %>%
    filter(year <= yearsref & year >= yearrefminusOne) %>%
    mutate(referenceFin = ifelse(year == yearsref & DOY == refday, 1, ifelse(year == yearsref & DOY > refday, NA, 0))) %>%
    filter(!is.na(referenceFin)) %>%
    as.data.frame()
  
  # Create sequence going back lastdays days before the reference day
  seqDays <- seq(1, nrow(tt), 1)
  newsequance <- rep(seqDays)
  
  ttup <- tt %>%
    mutate(days.reversed = rev(newsequance)) %>%
    filter(days.reversed < lastdays)
  
  #use !!sym; convert a string, here my variable name, to a symbol
  ttupfin <- ttup %>%
    arrange(days.reversed) %>%
    mutate(rolling_avg_tmean = zoo::rollmeanr(!!sym(variablemoving), k = rollwin, fill = NA, align = 'right')) %>%
    mutate(year = max(year)) %>%
    select(site_id, year, date, DOY, days.reversed, rolling_avg_tmean)
  
  return(ttupfin)
}

#function to just run correlation and beta regression model 
#define my covariation I want here for either correaltion or beta regression 
runing.movingwin.analysis = function(site = site,
                                     rolling.temperature.data = rolling.temperature.data,
                                     method = 'spearman',
                                     covariates.of.interest = 'rolling_avg_tmean',
                                     myform = formula('ScaledSeedProductionBR~rolling_avg_tmean')){
  
  #merge data seed to moving climate
  tible.sitelevel = site %>% 
    left_join(rolling.temperature.data, join_by(year)) %>% 
    drop_na(!!sym(covariates.of.interest))
  
  #define correlation - calculate correlation and se, extract also p value 
  n = tible.sitelevel %>% dplyr::select(year, site_id) %>% distinct() %>% nrow()
  
  correlation.all <- tible.sitelevel %>% 
    nest(data = -days.reversed) %>%
    mutate(correlation = map(data, ~cor.test(y=.$Value, x=.[[covariates.of.interest]], method = method)$estimate)) %>% 
    mutate(pvalue.cor = map(data, ~cor.test(y=.$Value, x=.[[covariates.of.interest]], method = method)$p.value))
  
  
  cortemp = correlation.all %>% 
    unnest(c(correlation, pvalue.cor)) %>% 
    dplyr::select(days.reversed, correlation, pvalue.cor) %>% 
    dplyr::mutate(correlation.se = correlation.spearman.se(.$correlation, n)) %>% 
    dplyr::mutate(site = unique(tible.sitelevel$site_id))
  
  #use purr for iteration 
  #here broom package used, because it is simple betareg (and not glmmTMB, need to adjust then )
  fitted_models = tible.sitelevel %>%
    nest(data = -days.reversed) %>%
    mutate(model = purrr::map(data, ~betareg(myform, data = ., na.action = na.omit)),
           tidied = purrr::map(model, tidy),
           glanced = purrr::map(model, glance),
           augmented = purrr::map(model, augment))
  
  slope = fitted_models %>%
    unnest(tidied) %>% 
    filter(str_detect(term, as.character(myform)[3])) %>% 
    dplyr::select(days.reversed,term, estimate, std.error, p.value) %>% 
    mutate(site  = unique(tible.sitelevel$site_id)) %>% 
    left_join(fitted_models %>%
                unnest(glanced) %>% 
                dplyr::select(days.reversed,pseudo.r.squared, logLik) %>% 
                mutate(site  = unique(tible.sitelevel$site_id))) %>% 
    dplyr::select(site, days.reversed, everything()) %>% 
    left_join(cortemp)
  
  return(slope)
}

#process the function for each site 
site.moving.climate.analysis = function(savefilmastree = savefilmastree,
                                        siteind = siteind,
                                        path.climate = here('climatedailyERA5/'),
                                        maturation.length = 1,
                                        refday = 274,
                                        rollwin = 7,
                                        phylolist = phylolist){
  #subset the data, 
  #scale seed production and do logit transformation (will use this logit for climwin later)
  site = savefilmastree %>% 
    dplyr::filter(sitenewname == siteind)%>% 
    dplyr::select(Alpha_Number, Species, Longitude, Latitude, sitenewname, Year, Value) %>% 
    dplyr::rename(site = sitenewname, 
           year = Year) %>% 
    mutate(ScaledSeedProduction = (Value - min(Value, na.rm = T))/(max(Value)-min(Value, na.rm = T))) %>% 
    mutate(ScaledSeedProductionBR = y.transf.betareg(ScaledSeedProduction)) %>% 
    mutate(date = as.Date(paste0(year, '-10-01')),
           logit.seed = car::logit(ScaledSeedProduction)) %>% #need to use the 0-1 data 
    filter(year < 2020) %>% 
    left_join(phylolist)
  
  #create vector site name
  #siteid = unique(site$site)
  
  #get the climate id
  climate = read_csv(paste0(path.climate, 'TemperatureData_', siteind, '.csv')) %>%
    mutate(
      temperature.degree = temperature - 273.15, #convert K to celsius 
      day = day(date),
      month = month(date),
      year = year(date),
      DOY = yday(date),
      day.month =format(date,"%m-%d")
    )
  
  # Set lastdays and yearneed based on maturation length
  if (maturation.length == 1) {
    lastdays <- 730  # 365 + 365/2
    yearneed <- 2
  } else if (maturation.length == 2) {
    lastdays <- 1095  # 365 + 365 + 365/2
    yearneed <- 3
  }

  
  # Define the year period
  yearperiod <- (min(climate$year) + yearneed):max(climate$year)
  
  # Apply the function across all years in yearperiod and combine results
  #here the map dfr will basically do same as aplly by runing the function over all the time period we want 
  rolling.temperature.data <- map_dfr(yearperiod, reformat.climate.backtothepast, 
                                      climate = climate, 
                                      yearneed = yearneed, 
                                      refday = refday, 
                                      lastdays = lastdays, 
                                      rollwin = rollwin,
                                      variablemoving = 'temperature.degree')
  
  #now run climate analysis 
  results.moving.site = runing.movingwin.analysis(site = site,
                                                  rolling.temperature.data = rolling.temperature.data,
                                                  method = 'spearman',
                                                  covariates.of.interest = 'rolling_avg_tmean',
                                                  myform = formula('ScaledSeedProductionBR~rolling_avg_tmean'))
  
  return(results.moving.site)
  
}

# Main function to process all sites
FULL.moving.climate.analysis <- function(savefilmastree = savefilmastree,
                                         path.climate = here('climatedailyERA5/'),
                                         maturation.length = 1,
                                         refday = 274,
                                         rollwin = 7,
                                         phylolist = phylolist) {
  al.sites <- unique(savefilmastree$sitenewname)
  
  results.moving <- map_dfr(al.sites, 
                            site.moving.climate.analysis, #the function I want to loop with 
                            savefilmastree = savefilmastree, 
                            path.climate = path.climate, 
                            maturation.length = maturation.length, 
                            refday = refday, 
                            rollwin = rollwin, 
                            phylolist = phylolist)
  
  return(results.moving)
}

#function from Nature Tackery paper 
#function to find concurrent window eg sequences of days
find_concurrent_period=function(temp_window,pred_C){
  
  #create dummy index
  idx=1
  #find differences between submitted dates and note whic are greater than 1
  diff_id = c(1,which(diff(temp_window)>1))
  
  if(length(diff_id)>1){if((which(diff(temp_window)>1))[1]==1){diff_id=diff_id[-1]}}
  if(length(temp_window)>1){
    #find all the  series of sequentioal dates and give unique id
    for(rv in 1:length(diff_id)){
      
      if(rv==length(diff_id)){
        idx=c(idx,rep(rv,length((diff_id[rv]+1):length(temp_window))))
      }else{
        idx=c(idx,rep(rv,length((diff_id[rv]+1):(diff_id[rv+1]))))
      }
      
    }
  }
  
  #from the estimated coefficients and the concurrent window indices (idx) find which period has the most extreme average. call that the period to use
  mx_coef = which.max(tapply(abs(pred_C$fit[temp_window]),idx,mean))
  
  return(temp_window[idx==mx_coef])
  
}

# Function to get prediction and windows based on Tackery and Simmonds et al code 
# To extract the key climate window of sensitivity from the GAMS
# we used code from Thackeray et al. 2016, which can be found at:
#https://github.com/NERC-CEH/Phenology_Climate/blob/master/Source_Code/Functions.R
# Lines 419 to 428
get_predictions_windows <- function(slope_gam, rs_gam, temporary) {
  pred_S <- predict(slope_gam, se.fit = TRUE, type = "response")
  pred_R <- predict(rs_gam, se.fit = TRUE, type = "response")
  
  coef_range_l <- which(((pred_S$fit - 1.96 * pred_S$se.fit) - 100) <= quantile(pred_S$fit - 100, 0.025))
  coef_range_u <- which(((pred_S$fit + 1.96 * pred_S$se.fit) - 100) >= quantile(pred_S$fit - 100, 0.975))
  r_sq_range <- which(((pred_R$fit + 1.96 * pred_R$se.fit) - 100) >= quantile(pred_R$fit - 100, 0.975))
  
  temp_window_l <- coef_range_l[is.element(coef_range_l, r_sq_range)]
  temp_window_u <- coef_range_u[is.element(coef_range_u, r_sq_range)]
  
  if (length(temp_window_l) == 0) temp_window_l <- coef_range_l
  if (length(temp_window_u) == 0) temp_window_u <- coef_range_u
  
  days <- find_concurrent_period(temp_window = c(temp_window_l, temp_window_u), pred_C = pred_S)
  list(pred_S = pred_S, pred_R = pred_R, days = days)
}

#function to optimize k for the gam it will provide knots value and two model one for slope and one for R2 
#run gam regression
#specify either cc because doing better with seasonal data as here or cr 
#"cc" indicates a cyclic cubic spline, which we want for the seasonal term as there should be no discontinuity between January and December.
#https://fromthebottomoftheheap.net/2014/05/09/modelling-seasonal-data-with-gam/

optimize_and_fit_gam <- function(temporary, optim.k = TRUE, plots = F) {
  if (optim.k) {
    # Function to check significance in k.check
    is_significant <- function(check) {
      p_values <- check[,'p-value']
      all(p_values < 0.05) # Returns TRUE if all p-values are significant
    }
    
    # Range of k values to try
    k_values <- seq(10, 365)  # Adjust to one per day 
    optimal_k <- NULL
    
    for (k in 1:length(k_values)) {
      # Fit the model
      slope_gam <- gam(slope ~ s(day, k = k_values[k], bs = "cr"), data = temporary)
      
      # Perform k.check
      check <- k.check(slope_gam)
      pvalue <- check[,'p-value']
      kfin <- check[,"k'"]
      
      # Check if all p-values are significant
      if (!is_significant(check)) {
        optimal_k <- k_values[k]
        break
      }
    }
    if (!is.null(optimal_k)) {
      cat("Optimal k found:", optimal_k, "\n")
      slope_gam <- gam(slope ~ s(day, k = optimal_k, bs = "cr"), data = temporary)
      rs_gam <- gam(r_s ~ s(day, k = optimal_k, bs = "cr"), data = temporary)
      k <- optimal_k
    } else {
      cat("No optimal k found within the specified range. Using default k = -1.\n")
      slope_gam <- gam(slope ~ s(day, k = -1, bs = "cr"), data = temporary)
      rs_gam <- gam(r_s ~ s(day, k = -1, bs = "cr"), data = temporary)
      k <- -1
    }
  } else {
    slope_gam <- gam(slope ~ s(day, k = -1, bs = "cr"), data = temporary)
    rs_gam <- gam(r_s ~ s(day, k = -1, bs = "cr"), data = temporary)
    k <- -1
  }
  
  if (plots) {
    results <- cowplot::plot_grid(
      draw(slope_gam, residuals = TRUE) + ylab('Slope (partial effect)'), 
      draw(rs_gam, residuals = TRUE) + ylab('R2 (partial effect)')
    )
    print(results)
  }
  
  list(slope_gam = slope_gam, rs_gam = rs_gam, k = k)
  
}

#create sequence, because some shit species does not follow one sequence (days for example going at 1 and then some much higher)
extract_consecutive_sequences <- function(values, keep_all = FALSE) {
  # Initialize variables
  sequences <- list()
  current_sequence <- integer(0)
  
  # Iterate through the values
  #if does not match afte +1 then it will create new sequence length 
  for (i in seq_along(values)) {
    if (i == 1 || values[i] == values[i - 1] + 1) {
      current_sequence <- c(current_sequence, values[i])
    } else {
      sequences <- c(sequences, list(current_sequence))
      current_sequence <- values[i]  # Start new sequence
    }
  }
  
  # Add the last sequence
  sequences <- c(sequences, list(current_sequence))
  
  if (!keep_all) {
    # Find the longest sequence
    longest_sequence <- sequences[[which.max(lengths(sequences))]]
    return(longest_sequence)
  } else {
    return(sequences)
  }
}
####################################
#MAYBE NEED TO THINK MATURATION FRUIT LENGTH to adapt the window
#get taxo for phylo maturation cycle 
phylolist = addphylogeny.family.genus(savefilmastree)

#get calendar data 
refday = 274
lastdays = 730
calendarFormat = goingbackpastdayscalendar(refday = refday, lastdays) %>% 
  filter(YEAR == 1950) %>% 
  dplyr::select(MONTHab, DOY, days.reversed) %>% 
  mutate(datefake  = as.Date(DOY, origin = "1948-01-01"),
         day.month = format(datefake,"%m-%d"), 
         year = case_when(
           days.reversed < 366 ~ 1,
           days.reversed >= 366 & days.reversed < 730 ~ 2,
           days.reversed >= 730 & days.reversed < 1095 ~ 3,
           TRUE ~ 4))

#let's think for a loop 
#maturation.length <- 1 #is it like a fagus species or a pinus species, need to think of that
# Define rollwing as it was used but not specified
#take less than 1hours to run 
savefilmastree.sub = savefilmastree %>% filter(sitenewname %in%c('6184_009_QUECER_1_7', '6184_010_QUECER_1_7'))
data.moving.roll = FULL.moving.climate.analysis(savefilmastree = savefilmastree,
                                         path.climate = here('climatedailyERA5/'),
                                         maturation.length = 1,
                                         refday = 274,
                                         rollwin = 7,
                                         phylolist = phylolist)

#qs::qsave(data.moving.roll,
#          here("data.moving.roll.7days.qs"))



data.moving.roll %>% 
  mutate(sign = ifelse(estimate >0, 'p', 'n'))%>% 
  ggplot(aes(x=days.reversed, y = estimate, group = days.reversed, col = sign, fill = sign))+
  geom_bar( stat="identity" , alpha = .9)

data.moving.roll %>% 
  ggplot(aes(x = estimate,
             y = correlation,
             size = pseudo.r.squared))+
  geom_point()

i  = 2

Results_CSP[[i]] %>% 
  mutate(sign = ifelse(estimate >0, 'p', 'n')) %>% 
  left_join(calendarFormat)%>% 
  ggplot(aes(x=days.reversed, y = estimate))+
  geom_bar(aes(col = sign, fill = sign), stat="identity" , alpha = .1)+
  theme(legend.position = 'none')+
  ylab('Spearman correlation rolling temperature')+
  xlab("Previous days")+
  scale_color_manual(values = rev(c("#00AFBB", "#FC4E07")))+
  facet_grid(.~year, scales = 'free')

Results_CSP[[i]] %>% 
  mutate(sign = ifelse(estimate >0, 'p', 'n')) %>% 
  left_join(calendarFormat)%>% 
  ggplot(aes(x=days.reversed, y = pseudo.r.squared))+
  geom_bar(aes(col = sign, fill = sign), stat="identity" , alpha = .1)+
  theme(legend.position = 'none')+
  ylab('Spearman correlation rolling temperature')+
  xlab("Previous days")+
  scale_color_manual(values = rev(c("#00AFBB", "#FC4E07")))+
  facet_grid(.~year, scales = 'free')


data.moving.roll %>% 
  ggplot(aes(x=days.reversed, y = pseudo.r.squared, group = days.reversed))+
  geom_bar( stat="identity" , alpha = .9)

#create a list for each of them 
Results_CSP = data.moving.roll %>%
  group_by(site) %>%
  group_split()
  
nameCSP =   data.moving.roll %>% 
  group_by(site) %>%
  group_keys()

# Set names of the list based on group keys
names(Results_CSP) <- apply(nameCSP, 1, paste)

output_fit_summary = NULL
for(i in 1:length(Results_CSP)){
  #now make subset and run regression
  list_slope <- as.list(Results_CSP[[i]]$estimate)
  list_rs <- as.list(Results_CSP[[i]]$pseudo.r.squared)
  
  day = seq(rollwin,(lastdays-1),1)
  slope = list_slope
  r_s = list_rs
  #k = nrow(site)
  temporary <- data.frame(slope = unlist(slope), 
                          day = day, 
                          r_s = unlist(r_s))
  #use the function to optimize k 
  results <- optimize_and_fit_gam(temporary, optim.k = T, plots = F)
  #get days windows 
  #need to add rollwin because it start at 1... 
  days = get_predictions_windows(slope_gam = results[[1]], 
                                 rs_gam = results[[2]], 
                                 temporary)$days+rollwin
  #subset the data, 
  #scale seed production and do logit transformation (will use this logit for climwin later)
  site = savefilmastree %>% 
    dplyr::filter(sitenewname == unique(Results_CSP[[i]]$site))%>% 
    dplyr::select(Alpha_Number, Species, Longitude, Latitude, sitenewname, Year, Value) %>% 
    dplyr::rename(site = sitenewname, 
                  year = Year) %>% 
    mutate(ScaledSeedProduction = (Value - min(Value, na.rm = T))/(max(Value)-min(Value, na.rm = T))) %>% 
    mutate(ScaledSeedProductionBR = y.transf.betareg(ScaledSeedProduction)) %>% 
    mutate(date = as.Date(paste0(year, '-10-01')),
           logit.seed = car::logit(ScaledSeedProduction)) %>% #need to use the 0-1 data 
    filter(year < 2020) %>% 
    left_join(phylolist)
  
  #get the climate id
  climate = read_csv(paste0(path.climate, 'TemperatureData_', unique(Results_CSP[[i]]$site), '.csv')) %>%
    mutate(
      temperature.degree = temperature - 273.15, #convert K to celsius 
      day = day(date),
      month = month(date),
      year = year(date),
      DOY = yday(date),
      day.month =format(date,"%m-%d")
    )
  
  # Set lastdays and yearneed based on maturation length
  if (maturation.length == 1) {
    lastdays <- 730  # 365 + 365/2
    yearneed <- 2
  } else if (maturation.length == 2) {
    lastdays <- 1095  # 365 + 365 + 365/2
    yearneed <- 3
  }
  
  
  # Define the year period
  yearperiod <- (min(climate$year) + yearneed):max(climate$year)
  
  # Apply the function across all years in yearperiod and combine results
  #here the map dfr will basically do same as aplly by runing the function over all the time period we want 
  rolling.temperature.data <- map_dfr(yearperiod, 
                                      reformat.climate.backtothepast, 
                                      climate = climate, 
                                      yearneed = yearneed, 
                                      refday = refday, 
                                      lastdays = lastdays, 
                                      rollwin = rollwin,
                                      variablemoving = 'temperature.degree')
  #now do the filtering based on windows found 
  windowsopen = days[1]
  windowsclose = tail(days, n=1) #add size of the moving window 
  moving.win = 'optimum'
  
  #keep only he window of interest 
  climate.windows.best = rolling.temperature.data %>% 
    dplyr::filter(days.reversed>=windowsopen & days.reversed<=windowsclose) %>% 
    group_by(site_id, year) %>% 
    summarise(mean.temperature = mean(rolling_avg_tmean, na.rm = T)) %>% 
    ungroup()
  
  if(nrow(climate.windows.best)==0){
    #it is possible that we might have shit because of gam identification (will identify not only one but additional sequences)
    sequences <- extract_consecutive_sequences(days, keep_all = F)
    windowsopen = sequences[1]
    windowsclose = tail(sequences, n=1) 
    
    climate.windows.best = rolling.temperature.data %>% 
      dplyr::filter(days.reversed>=windowsopen & days.reversed<=windowsclose) %>% 
      group_by(site_id, year) %>% 
      summarise(mean.temperature = mean(rolling_avg_tmean, na.rm = T)) %>% 
      ungroup()
    moving.win = 'longest.sequence'
  }
  
  #fit new beta reg 
  fit.best.model = site %>% 
    left_join(climate.windows.best) %>%
    betareg(ScaledSeedProductionBR~mean.temperature, data = .)
  
  intercept = tidy(fit.best.model)$estimate[1]
  intercept.se = tidy(fit.best.model)$std.error[1]
  estimate.model = tidy(fit.best.model)$estimate[2]
  estimate.se.model = tidy(fit.best.model)$std.error[2]
  pvalue.model = tidy(fit.best.model)$p.value[2]
  pseudo.r2 = glance(fit.best.model)$pseudo.r.squared
  AIC = glance(fit.best.model)$AIC
  nobs = glance(fit.best.model)$nobs
  
  #make summary talbe 
  output_fit_summary.temp <- data.frame(site = unique(Results_CSP[[i]]$site),
                                   reference.day = refday,
                                   windows.size = rollwin,
                                   knots.number = results$k,
                                   type.moving = moving.win,
                                   window.open = windowsopen,
                                   window.close = windowsclose,
                                   intercept = intercept,
                                   intercept.se = intercept.se,
                                   estimate = estimate.model,
                                   estimate.se = estimate.se.model,
                                   pvalue = pvalue.model,
                                   pseudo.R2 = pseudo.r2,
                                   AIC = AIC,
                                   nobs = nobs)
  output_fit_summary = rbind(output_fit_summary, output_fit_summary.temp)
  
}

output_fit_summary %>% 
  ggplot(aes(x = pseudo.R2,
             y = estimate,
             col = pvalue))+geom_point()

table(output_fit_summary$type.moving)

#ok now try with clinmwin 
#detachAllPackages()

library(climwin)
range = c(lastdays, 0)
#range = c(104, 0)
#52 weeks is one year 

optionwindows = 'absolute'
referenceclimwin = c(01, 10)
climwin_output_rw <- climwin::slidingwin(xvar = list(temperature.degree = climate$temperature.degree),
                                cdate = climate$date,
                                bdate = site$date,
                                baseline = lm(logit.seed ~ 1, data = site),
                                cinterval = "day",
                                range = range, 
                                type = optionwindows, 
                                refday = referenceclimwin, #absolute based on 01-10 1st october 
                                stat = "mean",
                                cmissing = 'method2',
                                func = "lin")

plotall(dataset = climwin_output_rw[[1]]$Dataset)
plotall(dataset = climwin_output_rw[[1]]$Dataset,
          bestmodel = climwin_output_rw[[1]]$BestModel,
          bestmodeldata = climwin_output_rw[[1]]$BestModelData)
climwin_output_rw$combos
top10windows = climwin_output_rw[[1]]$Dataset %>% 
  as_tibble() %>% 
  arrange(deltaAICc) %>% 
  slice_min(deltaAICc, n = 10)

hist(top10windows$WindowClose)
hist(top10windows$WindowOpen)

round(median(top10windows$WindowClose))
round(median(top10windows$WindowOpen))

#to be continued 
#check with climwin simpler stuff 
https://www.nature.com/articles/nature18608
and code her e
https://github.com/NERC-CEH/Phenology_Climate/blob/master/Source_Code/Functions.R
and other analysis 
https://besjournals.onlinelibrary.wiley.com/doi/10.1111/1365-2656.13038


#make small climate comparison eobs for this site 
#siteid = '6012_001_FAGSYL_1_9'
#listClimateTmean <- list.files(path = "/Users/vjourne/Library/CloudStorage/Dropbox/MastreeForemast/climateSiteDaily", full.names=T)#[1]
#climateEOBS = read_csv(listClimateTmean[106])
#test = climateEOBS %>% dplyr::select(year, month, day, Tmean) %>% 
#  mutate(month = as.numeric(month),
#         day = as.numeric(day)) %>% 
#  left_join(climate %>% dplyr::select(day, month, year, temperature.degree))
#plot(test$temperature.degree, test$Tmean)

#ok so going too far diifuclt to identify cues
#let's try to limit to 1years and half for most 1 year fruit maturation 
#and two years and half for 2 years maturation species 

#add ref
#https://www.nature.com/articles/s41467-023-43744-8
#https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0041010

