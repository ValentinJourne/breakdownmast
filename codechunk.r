##### Rcode by Valentin and Maciej - we investigate masting breakdown based on data from MASTREE+ 
# the project started in june 2024 and is shared in Val's GitHub

########## the first step - loading libraries
## here for finding the directory
## tidyverse for manipulating data and cleaning code 
## zoo for computing rolling windows
## ggplot for plotting data
## mgcv for running GAM models

library(here)
library(tidyverse)
library(zoo)
library(ggplot)
library(mgcv)


########## the second step - getting working directory and loading data from MASTREEplus (stored in github)

here()
#here('/Users/vjourne/Documents/GITprojects/')

data = read_csv(here('MASTREEplus_2024-02-20_V2.csv'))

## show the data and explore the variables
str(data)
table(data$Spatial_unit)
table(data$Variable)

########## the third step - writing a function to extract the data for coding = cleaning and filtering dataset
## we filter the dataset using the records starting in second half of XX century (year = 1961) and measuring seed production for at least 14 years in one site
## we create a new column which holds the unique name combined with:
# Alpha_Number - an unique number of the study; Site_Number - an unique number of location; Species_code - species name; 
# VariableFactor - type of seed data as factor; UnitFactor - 

determineSiteNew <- function(data = data , yearselection = 1961, yearending = 2021, filterobs = 14){
  
  dfupdate <- data %>% 
    mutate(sitenewname= paste0(Alpha_Number, "_",Site_number, "_", Species_code, "_", VariableFactor, "_", UnitFactor)) %>% 
    filter(Year > yearselection & Year < yearending) %>% 
    group_by(sitenewname) %>% 
    mutate(newlengthdate = n()) %>% 
    filter(newlengthdate>filterobs) %>% 
    ungroup()
  
  return(dfupdate)
}  


## alternatively load the function saved in github 
#source("determineSiteNew.R")


savefilmastree = determineSiteNew(data = data %>% mutate(VariableFactor = as.character(as.numeric(as_factor(Variable))),
                                                         UnitFactor = as.character(as.numeric(as_factor(Unit))), 
                                                         yearselection = 1979,
                                                         yearending = 2021,
                                                         filterobs = 19))%>%
  filter(!Variable == "pollen" & !Spatial_unit == 'super-region' & VarType == "C") %>% 
  mutate(hemisphere = ifelse(Latitude > 0, "North", "South")) %>% 
  mutate(TropNoTrops = ifelse(Latitude < -23.44 | Latitude > 23.44, "Non-Tropics", 'Tropics'))

##########################
#for Val extraction climate
#convert to shp with WGS84
site.mastree.forclim = savefilmastree %>% ungroup() %>% 
  dplyr::filter(TropNoTrops == "Non-Tropics") %>% 
  dplyr::select(sitenewname, Longitude, Latitude) %>% distinct() %>% 
  sf::st_as_sf( coords = c("Longitude", "Latitude"), crs = 4326)
sf::st_write(site.mastree.forclim, 
             dsn = "shp.forclimate/output_mastreesite_climatelist.shp",
             append=FALSE,
             driver = "ESRI Shapefile")
#write_csv2(site.mastree.forclim, 'listname.csv')
#plot(site.mastree.forclim)
##########################

ggplot(savefilmastree%>%
         filter(Species_code == 'FAGSYL' & Year > 1950), 
       aes( y = Value, x = Year, group = sitenewname))+geom_point()+geom_line()+
  facet_wrap(.~VarType, scales = 'free')

#creat subset only for temperate sites pop 
temperate.mastree = savefilmastree  %>% 
  filter(TropNoTrops == "Non-Tropics") %>% 
  dplyr::select(sitenewname, Species, Year, Value, Collection_method, Variable, Latitude, Longitude)

zoo::rollmeanr
zoo::rollapply

lengthtime = 10 

test  = temperate.mastree %>% 
  filter(sitenewname == '6210_005_FAGSYL_1_7')%>%
  group_by(sitenewname, Species, Collection_method, Variable, Latitude, Longitude)%>%
  mutate(meanseed = mean(Value),
         sdseed = sd(Value))%>%
  mutate(seednew = (Value - meanseed)/sdseed)

ggplot(test , aes(x= Year, y = seednew))+geom_line()
ggplot(test , aes(x= Year, y = Value))+geom_line()

lengthtime = 20
data.forroll = temperate.mastree %>% 
  group_by(sitenewname, Species, Collection_method, Variable, Latitude, Longitude)%>%
  mutate(meanseed = mean(Value), sdseed = sd(Value))%>%
  mutate(seednew = (Value - meanseed)/sdseed,
         scale.seed = scale(Value))%>%
  complete(Year = full_seq(Year, period = 1), fill = list(seednew = NaN))


vector10year = data.forroll %>%
  group_by(sitenewname, Species, Collection_method, Variable, Latitude, Longitude)%>% 
  mutate(row_id=row_number()) %>% 
  filter(row_id<=10) %>% 
  summarise(firstmean10years = mean(scale.seed),
            firstcv10years.scaled = raster::cv(scale.seed)/100,
            firstcv10years = raster::cv(Value)/100)

plot(vector10year$firstcv10years.scaled, vector10year$firstcv10years)

rollingtemperate = data.forroll %>%
  group_by(sitenewname, Species, Collection_method, Variable, Latitude, Longitude)%>% 
  mutate(row_id=row_number()) %>% 
  left_join(vector10year) %>% 
  #mutate(firstmean10years = ifelse(row_id <= 10, mean(dplyr::slice(., 1:10)$scale.seed, na.rm = TRUE), NA)) %>% 
  #tidyr::fill(firstmean10years, .direction = 'down') %>% 
  group_by(sitenewname, Species, Collection_method, Variable, Latitude, Longitude)%>% 
  mutate(meanrolling = zoo::rollmeanr(seednew, k=lengthtime, fill=NA, align='right'),
         sdrolling = RcppRoll::roll_sd(Value, fill=NA, n=lengthtime, align='right'),
         cvroll = zoo::rollapply(Value, FUN = raster::cv, width=lengthtime, by = 1, fill=NA, align='right')/100,
         cv.manual = sdrolling/zoo::rollmeanr(Value, k=lengthtime, fill=NA, align='right'),
         anomaliesMeanSeed = meanrolling-firstmean10years,
         anomaliesCVpSeed = cvroll - firstcv10years)%>% #filter(Species == 'Fagus sylvatica') %>% 
          filter(Year > 1950)


plot(rollingtemperate$anomaliesMeanSeed, rollingtemperate$anomaliesCVpSeed)



rollingtemperate %>% #filter(Species == 'Fagus sylvatica') %>% 
  filter(Year > 1950) %>% 
  ggplot(aes(x = Year, y = anomaliesMeanSeed, group = sitenewname, col = Species))+
  geom_point()+geom_line()+theme(legend.position = 'none')


small.model = mgcv::gam(anomaliesMeanSeed~s(Year)+s(species,bs="re"),
                   data = rollingtemperate %>% mutate(species = as_factor(Species)))

summary(small.model)
library(gratia)
draw(small.model)

small.model.lm =glmmTMB::glmmTMB(cvroll~Year,
                        data = rollingtemperate %>% mutate(species = as_factor(Species))%>% #filter(Species == 'Fagus sylvatica') %>% 
                          filter(Year > 1980 & Species != ' Fagus sylvatica'))
summary(small.model.lm)

pred = ggeffects::ggpredict(small.model.lm)
plot(pred)

small.model.lmmean =glmmTMB::glmmTMB(meanrolling~Year,
                                 data = rollingtemperate %>% mutate(species = as_factor(Species))%>% #filter(Species == 'Fagus sylvatica') %>% 
                                   filter(Year > 1980 & Species != ' Fagus sylvatica'))

summary(small.model.lmmean)
pred2 = ggeffects::ggpredict(small.model.lmmean)
plot(pred2)


plot(rollingtemperate$seednew, rollingtemperate$scale.seed)
plot(rollingtemperate$cv.manual, rollingtemperate$cvroll)
rollingtemperate %>% ggplot(aes(cvroll, meanrolling))+geom_point()
unique(rollingtemperate$Species)
hist(rollingtemperate$cvroll)
rollingtemperate %>% filter(Species == 'Fagus sylvatica') %>% filter(Year > 1950) %>% 
  #ggplot(aes(x = Year, y = seednew, group = sitenewname, col = Species))+geom_point()+geom_line()+theme(legend.position = 'none')
  ggplot(aes(x = Year, y = cv.manual*100, group = sitenewname, col = Species))+geom_point()+geom_line()+theme(legend.position = 'none')

hist(rollingtemperate$cvroll)

test  = temperate.mastree %>% 
  group_by(sitenewname, Species, Collection_method, Variable, Latitude, Longitude)%>%
  mutate(meanseed = mean(Value),
         sdseed = sd(Value))%>%
  mutate(seednew = (Value - meanseed)/sdseed)
