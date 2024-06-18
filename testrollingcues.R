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

#make calendar to associate DOY 
goingbackpastdayscalendar <-function(refday = 244, lastdays = 1095){
  
  DATE = seq(as.Date("1948-01-01"), as.Date("1953-01-01"), by="days")
  MONTH =  format(as.Date(DATE, format="%Y-%m-%d"),"%m") %>% as.numeric()
  MONTHab = month.abb[MONTH]
  YEAR =  format(as.Date(DATE, format="%Y-%m-%d"),"%Y") %>% as.numeric()
  
  DOY = yday(DATE)
  dfata = data.frame(DATE,YEAR,MONTHab, DOY)
  yearperiod = 1948:1953
  sizevec = length(unique(YEAR))-2
  refday = refday
  vectotemp = NULL
  
  for(k in 1:sizevec){
    #year +3, becusse of the 36 month analysis 
    yearsref = yearperiod[k]
    yearrefminusOne <- yearsref-2 #PREVIOUS = 3
    
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
      mutate(days.reversed = rev(newsequance))  %>% #because climwin start at 0 and not 1
      filter(days.reversed< lastdays )
    ttupfin = ttup %>%
      arrange(days.reversed)  %>% 
      mutate(YEAR = max(YEAR))      
    vectotemp <- rbind(vectotemp, ttupfin) 
  }
  
  return(vectotemp)
}


site = savefilmastree %>% 
  filter(sitenewname == '2543_002_FAGSYL_3_7') %>% 
  dplyr::select(Alpha_Number, Species, Longitude, Latitude, sitenewname, Year, Value) %>% 
  rename(site = sitenewname, year = Year) %>% 
  mutate(ScaledSeedProduction = (Value - min(Value, na.rm = T))/(max(Value)-min(Value, na.rm = T))) %>% 
  mutate(ScaledSeedProductionBR = y.transf.betareg(ScaledSeedProduction))

plot(site$Value, site$ScaledSeedProductionBR)

siteid = unique(site$site)

climate = read_csv(here('climatedailyERA5/TemperatureData_2543_002_FAGSYL_3_7.csv')) %>%
  mutate(
    temperature.degree = temperature - 273.15,
    day = day(date),
    month = month(date),
    year = year(date),
    DOY = yday(date),
    day.month =format(date,"%m-%d")
  )

yearneed = 2
lastdays = 730
refday = 274
rollwin = 1
method = 'spearman'

yearperiod = (min(climate$year)+yearneed):max(climate$year)
#keep vector for formatted data
datarollingTemp = NULL
for(k in 1:length(yearperiod)){
  yearsref = yearperiod[k]
  yearrefminusOne <- yearsref-yearneed
  tt <- climate %>% 
    dplyr::filter(year <= yearsref & year >= yearrefminusOne) %>% 
    dplyr::mutate(referenceFin = ifelse(year == yearsref & DOY == refday, 1,
                                 ifelse(year == yearsref & DOY > refday, NA, 0))) %>% 
    dplyr::filter(!is.na(referenceFin)) %>% 
    as.data.frame()
  #create sequence going back 365 month before 
  seqDays <- seq(1,nrow(tt),1)
  newsequance <- rep(seqDays)
  ttup <- tt %>% 
    mutate(days.reversed = rev(newsequance))  %>% 
    filter(days.reversed< lastdays )
  ttupfin = ttup %>%
    arrange(days.reversed) %>% 
    mutate(rolling_avg_tmean = zoo::rollmeanr(temperature.degree, k=rollwin, fill=NA, align='right')) %>% 
    mutate(year = max(year)) %>% 
    dplyr::select(site_id, year, date, DOY, days.reversed, 
                  rolling_avg_tmean)      
  datarollingTemp <- rbind(datarollingTemp, ttupfin) 
}


tible.sitelevel = site %>% 
  left_join(datarollingTemp, join_by(year)) %>% 
  drop_na(rolling_avg_tmean)

#define correlation
n = tible.sitelevel %>% dplyr::select(year, site_id) %>% distinct() %>% nrow()
correlation.all <- tible.sitelevel %>% 
  nest(data = -days.reversed) %>%
  mutate(correlation = map(data, ~cor.test(y=.$Value, x=.$rolling_avg_tmean, method = method)$estimate)) %>% 
  mutate(pvalue.cor = map(data, ~cor.test(y=.$Value, x=.$rolling_avg_tmean, method = method)$p.value))

cortemp = correlation.all %>% 
  unnest(c(correlation, pvalue.cor)) %>% 
  dplyr::select(days.reversed, correlation, pvalue.cor) %>% 
  dplyr::mutate(correlation.se = correlation.spearman.se(.$correlation, n)) %>% 
  dplyr::mutate(site = siteid)

calendarFormat = goingbackpastdayscalendar(refday = 244, lastdays) %>% 
  filter(YEAR == 1950) %>% 
  dplyr::select(MONTHab, DOY, days.reversed) %>% 
  mutate(datefake  = as.Date(DOY, origin = "1948-01-01"),
         day.month = format(datefake,"%m-%d")) %>% 
  mutate(year = ifelse(days.reversed < 366, 1, ifelse(days.reversed >= 366 & days.reversed<730 , 2, 3)))

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

#Celsius = Kelvin – 273.15
#ggplot(climate %>% filter(year < 1990), 
#       aes(y = temperature.degree , x = DOY, group = year, col = year))+
#  geom_point()+geom_line()

library(betareg)
library(broom.mixed)
library(broom)

myform = formula('ScaledSeedProductionBR~rolling_avg_tmean')

#use purr for iteration 
fitted_models <- tible.sitelevel %>%
  nest(data = -days.reversed) %>%
  mutate(model = purrr::map(data, ~betareg(myform, data = ., na.action = na.omit)),
         tidied = purrr::map(model, tidy),
         glanced = purrr::map(model, glance),
         augmented = purrr::map(model, augment))

r2temp = fitted_models %>%
  unnest(glanced) %>% 
  dplyr::select(days.reversed,pseudo.r.squared, logLik) %>% 
  mutate(site  = siteid)

r2temp %>% 
  ggplot(aes(x=days.reversed, y = pseudo.r.squared, group = days.reversed))+
  geom_bar( stat="identity" , alpha = .9)
  
slope = fitted_models %>%
  unnest(tidied) %>% 
  filter(str_detect(term, as.character(myform)[3])) %>% 
  dplyr::select(days.reversed,term, estimate, std.error, p.value) %>% 
  mutate(site  = siteid)

slope %>% 
  ggplot(aes(x=days.reversed, y = estimate, group = days.reversed))+
  geom_bar( stat="identity" , alpha = .9)

slope %>% left_join(cortemp) %>% left_join(r2temp) %>% 
  ggplot(aes(x = estimate,
             y = correlation,
             size = pseudo.r.squared))+
  geom_point()

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


Results_CSP = slope %>% left_join(cortemp) %>% left_join(r2temp)
list_slope <- as.list(Results_CSP$estimate)
list_rs <- as.list(Results_CSP$pseudo.r.squared)

day = seq(1,(lastdays-rollwin),1)
slope = list_slope
r_s = list_rs

temporary <- data.frame(slope = unlist(slope), day = day, r_s = unlist(r_s))
slope_gam <- gam(slope ~ s(day, k = 18), data=temporary)
rs_gam <- gam(r_s ~ s(day, k = 18), data=temporary)

library(gratia)
draw(slope_gam)
draw(rs_gam)
plot(rs_gam, ylab = "Coefficient", xlab = "Day prior to 20th May", rug=F, las=1)
plot(slope_gam, ylab = "Coefficient", xlab = "Day prior to 20th May", rug=F, las=1)

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
actual_days <- refday-days # dates in day of year
