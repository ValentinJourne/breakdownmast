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
    dplyr::mutate(site = siteid)
  
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
    mutate(site  = siteid) %>% 
    left_join(fitted_models %>%
                unnest(glanced) %>% 
                dplyr::select(days.reversed,pseudo.r.squared, logLik) %>% 
                mutate(site  = siteid)) %>% 
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
    rename(site = sitenewname, year = Year) %>% 
    mutate(ScaledSeedProduction = (Value - min(Value, na.rm = T))/(max(Value)-min(Value, na.rm = T))) %>% 
    mutate(ScaledSeedProductionBR = y.transf.betareg(ScaledSeedProduction)) %>% 
    mutate(date = as.Date(paste0(year, '-10-01')),
           logit.seed = car::logit(ScaledSeedProduction)) %>% #need to use the 0-1 data 
    filter(year < 2020) %>% 
    left_join(phylolist)
  
  #create vector site name
  siteid = unique(site$site)
  
  #get the climate id
  climate = read_csv(paste0(path.climate, 'TemperatureData_', siteid, '.csv')) %>%
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
  
  #get calendar data 
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
  results.moving.site = runing.movingwin.analysis(site = site,rolling.temperature.data = rolling.temperature.data,
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


####################################
#MAYBE NEED TO THINK MATURATION FRUIT LENGTH to adapt the window
#get taxo for phylo maturation cycle 
phylolist = addphylogeny.family.genus(savefilmastree)

#let's think for a loop 
#maturation.length <- 1 #is it like a fagus species or a pinus species, need to think of that
# Define rollwing as it was used but not specified
#take less than 1hours to run 
data.moving.roll = FULL.moving.climate.analysis(savefilmastree = savefilmastree,
                                         path.climate = here('climatedailyERA5/'),
                                         maturation.length = 1,
                                         refday = 274,
                                         rollwin = 7,
                                         phylolist = phylolist)

qs::qsave(data.moving.roll,
          here("data.moving.roll.7days.qs"))



data.moving.roll %>% 
  mutate(sign = ifelse(estimate >0, 'p', 'n'))%>% 
  ggplot(aes(x=days.reversed, y = estimate, group = days.reversed, col = sign, fill = sign))+
  geom_bar( stat="identity" , alpha = .9)

data.moving.roll %>% 
  ggplot(aes(x = estimate,
             y = correlation,
             size = pseudo.r.squared))+
  geom_point()

tt = cortemp %>% 
  group_by(days.reversed) %>% 
  summarize_at(vars(correlation), list(mean.cor = mean), na.rm = TRUE) %>% 
  mutate(sign = ifelse(mean.cor >0, 'p', 'n')) %>% 
  left_join(calendarFormat)%>% 
  mutate(yearfac = as_factor(year)) 

tt %>% 
  ggplot(aes(x=days.reversed, y = mean.cor))+
  geom_bar(aes(col = sign, fill = sign), stat="identity" , alpha = .1)+
  theme(legend.position = 'none')+
  ylab('Spearman correlation rolling temperature')+
  xlab("Previous days")+
  scale_color_manual(values = rev(c("#00AFBB", "#FC4E07")))+
  facet_grid(.~year, scales = 'free')


r2temp %>% 
  ggplot(aes(x=days.reversed, y = pseudo.r.squared, group = days.reversed))+
  geom_bar( stat="identity" , alpha = .9)




Results_CSP = slope %>% left_join(cortemp) %>% left_join(r2temp)
list_slope <- as.list(Results_CSP$estimate)
list_rs <- as.list(Results_CSP$pseudo.r.squared)

day = seq(1,(lastdays-rollwin),1)
slope = list_slope
r_s = list_rs

k = nrow(site)
temporary <- data.frame(slope = unlist(slope), day = day, r_s = unlist(r_s))

slope_gam <- gam(slope ~ s(day, k = 1), data=temporary) #, k = k
rs_gam <- gam(r_s ~ s(day), data=temporary)


draw(slope_gam)
draw(rs_gam)
#plot(rs_gam, ylab = "Coefficient", xlab = "Day prior to 20th May", rug=F, las=1)
#plot(slope_gam, ylab = "Coefficient", xlab = "Day prior to 20th May", rug=F, las=1)

# create predictions from those gams
pred_S<-predict(slope_gam,se.fit=TRUE,type="response")
pred_R<-predict(rs_gam,se.fit=TRUE,type="response")

# To extract the key climate window of sensitivity from the GAMS
# we used code from Thackeray et al. 2016, which can be found at:
#https://github.com/NERC-CEH/Phenology_Climate/blob/master/Source_Code/Functions.R
# Lines 419 to 428
##find dates which exceed the extreme values defined by the upper and lower quantiles of the estimated coefficients
coef_range_l <- which(((pred_S$fit-1.96*pred_S$se.fit)-100)<=quantile(pred_S$fit-100,0.025))
coef_range_u <- which(((pred_S$fit+1.96*pred_S$se.fit)-100)>=quantile(pred_S$fit-100,0.975))

#simmilarly find which dates correspond to the high r squared values
r_sq_range <- which(((pred_R$fit+1.96*pred_R$se.fit)-100)>=quantile(pred_R$fit-100,0.975))

#windows are defined as where the extreme coefficient dates overlap the extreme r squared dates 
temp_window_l <- coef_range_l[is.element(coef_range_l,r_sq_range)]
temp_window_u <- coef_range_u[is.element(coef_range_u,r_sq_range)]

##if no dates overlap, then the extreme coefficient dates are chosen
if(length(temp_window_l)==0){temp_window_l = coef_range_l}
if(length(temp_window_u)==0){temp_window_u = coef_range_u}

##ensure that a concurrent period is selected. 
days = find_concurrent_period(temp_window = c(temp_window_l, temp_window_u),
                              pred_C = pred_S)
days
# This produces a vector of the days that are included in the critical window
# days are in fact days before 140 so window open = days[1]


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

