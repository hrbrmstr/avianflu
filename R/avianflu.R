library(rvest)     # scraping
library(stringr)   # string manipulation
library(lubridate) # date conversion
library(dplyr)     # data mjnging
library(zoo)       # for locf
library(ggplot2)   # plotting
library(rgdal)     # map stuff
library(rgeos)     # map stuff

#' USDA site containing the avian flu data

url <- "http://www.aphis.usda.gov/wps/portal/aphis/ourfocus/animalhealth/sa_animal_disease_information/sa_avian_health/ct_avian_influenza_disease/!ut/p/a1/lVJbb4IwFP41e1qwFZDLI-oUnGgyswm8kAMUaAaFQNG4X7-ibnEPYtakDz3nO_kupyhAHgoYHGgGnFYMiv4daOFqa8vjKZad5c58wc7mY-Eaa13Z2qoA-AKA7xwL_53fvjpaP_-Gp_Z8jHcK2qMABTHjNc-RD3VO2zCuGCeMhwWNGmhOT7iFsOqaMK3irj2_gNESijAnUPD8tpLQlkBLQsrSqinPJi7tAwX2i4_5tSBgRUfYF_wM9mLqmCbIj2QzxZpMJMUYg6TGkSLBBCaSPEnSJIljXVH0q_kBdw_CO5sXkNnSslV9LQJTDRk7czGumy7GjnYFDOTrCw36XRJTRbt_mlo9VD1FgfuctqrVByA37szNBAPwXOpzR97gPi7tm30gb2AfQkxWVJH4ifvZLavFIsUQrA1JSUOaUV61HHnH43HUtQmMsuqA6vK9NJST9JluNlKwyPr7DT6YvRs!/?1dmy&urile=wcm%3apath%3a%2Faphis_content_library%2Fsa_our_focus%2Fsa_animal_health%2Fsa_animal_disease_information%2Fsa_avian_health%2Fsa_detections_by_states%2Fct_ai_pacific_flyway"

#' read in the data, extract the table and clean up the fields
#' also clean up the column names to since they are fairly nasty

pg <- html(url)

pg %>%
  html_nodes("table") %>%
  magrittr::extract2(1) %>%
  html_table(header=TRUE) %>%
  filter(`Flock size`!="pending") %>%
  mutate(Species=str_replace(tolower(Species), "s$", ""),
         `Avian influenza subtype*`=str_replace_all(`Avian influenza subtype*`, " ", ""),
         `Flock size`=as.numeric(str_replace_all(`Flock size`, ",", "")),
         `Confirmation date`=as.Date(mdy(`Confirmation date`))) %>%
  rename(state=State, county=County, flyway=Flyway, flock_type=`Flock type`,
         species=Species, subtype=`Avian influenza subtype*`, date=`Confirmation date`,
         flock_size=`Flock size`) -> birds

#' save out our data frame

write.csv(birds, "data/birds.csv", row.names=FALSE)

#' setup our hex map

us <- readOGR("data/us_states_hexgrid.geojson", "OGRGeoJSON")
centers <- cbind.data.frame(data.frame(gCentroid(us, byid=TRUE), id=us@data$iso3166_2), stringsAsFactors=FALSE)
us_map <- fortify(us, region="iso3166_2")

#' get our data in the right format. we want weekly cumulative totals by state
#' but we need them for each week, so we have to "fill in" values (a.k.a.
#' last observation carried forward). i also kinda cherry picked ranges for
#' the color breaks. doing it more scientifically is an exercise left to the
#' reader.

birds %>%
  mutate(week=as.numeric(format(birds$date, "%Y%U"))) %>%
  arrange(week) %>%
  group_by(week, state) %>%
  tally(flock_size) %>%
  group_by(state) %>%
  mutate(cum=cumsum(n)) %>%
  ungroup %>%
  select(week, state, cum) %>%
  mutate(week=as.Date(paste(week, 1), "%Y%U %u")) %>%
  left_join(tidyr::expand(., week, state), .) %>%
  group_by(state) %>%
  do(na.locf(.)) %>%
  mutate(state_abb=state.abb[match(state, state.name)],
         cum=as.numeric(ifelse(is.na(cum), 0, cum)),
         brks=cut(cum,
                  breaks=c(0, 200, 50000, 1000000, 10000000, 50000000),
                  labels=c("1-200", "201-50K", "50k-1m",
                           "1m-10m", "10m-50m"))) -> by_state_and_week

i <- 0

for (wk in unique(by_state_and_week$week)) {

  # filter by week

  by_state_and_week %>% filter(week==wk) -> this_wk

  # hack to let us color the state labels in white or black depending on
  # the value of the fill

  this_wk %>%
    filter(brks %in% c("1m-10m", "10m-50m")) %>%
    .$state_abb %>%
    unique -> white_states

  centers %>%
    mutate(txt_col="black") %>%
    mutate(txt_col=ifelse(id %in% white_states, "white", "black")) -> centers

  # setup the plot

  gg <- ggplot()
  gg <- gg + geom_map(data=us_map, map=us_map,
                      aes(x=long, y=lat, map_id=id),
                      color="white", fill="#dddddd", size=2)
  gg <- gg + geom_map(data=this_wk, map=us_map,
                      aes(fill=brks, map_id=state_abb),
                      color="white", size=2)
  gg <- gg + geom_text(data=centers,
                       aes(label=id, x=x, y=y, color=txt_col), size=4)
  gg <- gg + scale_color_identity()
  gg <- gg + scale_fill_brewer(name="Combined flock size\n(all types)",
                               palette="RdPu", na.value="#dddddd", drop=FALSE)
  gg <- gg + guides(fill=guide_legend(override.aes=list(colour=NA)))
  gg <- gg + coord_map()
  gg <- gg + labs(x=NULL, y=NULL,
                  title=sprintf("U.S. Avian Flu Total Impact as of %s\n", wk))
  gg <- gg + theme_bw()
  gg <- gg + theme(plot.title=element_text(face="bold", hjust=0, size=24))
  gg <- gg + theme(panel.border=element_blank())
  gg <- gg + theme(panel.grid=element_blank())
  gg <- gg + theme(axis.ticks=element_blank())
  gg <- gg + theme(axis.text=element_blank())
  gg <- gg + theme(legend.position="bottom")
  gg <- gg + theme(legend.direction="horizontal")
  gg <- gg + theme(legend.title.align=1)

  # save the image

  png(sprintf("output/%03d.png", i), width=800, height=500, type="quartz")
  print(gg)
  dev.off()

  i <- i + 1

}

system("convert -delay 60 -loop 1 output/*png output/avian.gif")
