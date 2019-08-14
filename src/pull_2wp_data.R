library(tidyverse)
library(sf)

#### Match reaches to sites ####

# Read in the Delaware PRMS-stream_temp stream network and the inventory of
# temperature observation sites
del_prms <- sf::read_sf('data/delaware_stream_temp_by_segment/delaware_segments/delaware_segments.shp')
wqp_partitions <- feather::read_feather('../2wp-temp-observations/1_wqp_pull/inout/wqp_pull_partitions.feather')
wqp_sites <- feather::read_feather('../2wp-temp-observations/1_wqp_pull/inout/wqp_inventory.feather') %>%
  left_join(wqp_partitions, by='MonitoringLocationIdentifier') %>%
  filter(ResolvedMonitoringLocationTypeName %in% c('Spring', 'Stream'))

# Match reaches to sites
source('src/subset_closest.R') #need sp, SearchTrees packages
system.time({
  crosswalk_site_reach <- subset_closest(reaches=del_prms, sites=wqp_sites) %>%
    sf::st_transform(crs=sf::st_crs(del_prms))
})
# ought to filter here to sites that are within some max distance of a segment;
# right now we're getting everything in the bbox
site_points <- crosswalk_site_reach %>%
  sf::st_drop_geometry() %>%
  sf::st_as_sf(coords=c('site_lon', 'site_lat')) %>%
  sf::st_set_crs(4326) %>%
  sf::st_transform(crs=sf::st_crs(del_prms)) %>%
  left_join(dplyr::select(wqp_sites, MonitoringLocationIdentifier, resultCount), by='MonitoringLocationIdentifier') %>%
  mutate(
    log10ResultCount = log10(resultCount),
    resultCountBin = base::cut(resultCount, breaks=c(min(resultCount), 10, 100, 1000, max(resultCount)), right=FALSE))
levels(site_points$resultCountBin) <- c('1-10', '10-100', '100-1000', '1000+')
reach_points <- crosswalk_site_reach %>%
  sf::st_drop_geometry() %>%
  sf::st_as_sf(coords=c('reach_median_lon', 'reach_median_lat')) %>%
  sf::st_set_crs(4326) %>%
  sf::st_transform(crs=sf::st_crs(del_prms))

#### Plot inventory and geospatial overviews ####

graphics.off()

# Guidelines for sf plotting:
# Adding to a plot of an sf object only works when using reset=FALSE in the first plot.
# If the first plot is geometry, reset can (and should) be TRUE.

# plot reaches and the reach medians
png('out/reach_centroids.png', width=600, height=960)
plot(crosswalk_site_reach['seg_id_nat'], reset=FALSE) 
plot(reach_points['seg_id_nat'], pch=20, add=T)
dev.off()

# plot 
png('out/dist_to_reach_hist.png', width=480, height=400)
ggplot(site_points, aes(x=dist_to_reach_km)) + geom_histogram(binwidth=1)
dev.off()

# plot sites colored by their distance to a matched reach
png('out/dist_to_reach.png', width=600, height=960)
plot(filter(site_points['dist_to_reach_km'], dist_to_reach_km <= 5), pch=20, reset=FALSE)
plot(sf::st_geometry(crosswalk_site_reach), add=TRUE)
dev.off()

# plot sites colored by their number of reported observations
site_cex <- function(site_points) {
  as.numeric(site_points$resultCountBin)/2
}
png('out/log10_result_count.png', width=600, height=960)
par(omi=c(0, 0, 0, 1))
filter(site_points, dist_to_reach_km <= 5) %>%
  { plot(.['resultCountBin'], pch=20, cex=site_cex(.), reset=FALSE) }
plot(sf::st_geometry(crosswalk_site_reach), col='lightgray', add=TRUE) 
filter(site_points, dist_to_reach_km <= 5, resultCount > 100) %>%
  { plot(.['resultCountBin'], pch=20, cex=site_cex(.), add=TRUE) }
dev.off()

#### Collect actual observations ####

# Identify sites
del_wqp_sites <- left_join(
  filter(site_points, dist_to_reach_km <= 5),
  dplyr::select(wqp_sites, -resultCount, -latitude, -longitude), by='MonitoringLocationIdentifier')

# Read in the observation files.
#   These files first need to be manually downloaded from
# https://drive.google.com/drive/folders/1RqKhyYi4tHvD1fli6LKocw8xHwm1PYmO to
# data/WQPTempObs190808
#   This wasn't working for me:
# library(googledrive)
# library(tidyverse)
# drive_ls('1RqKhyYi4tHvD1fli6LKocw8xHwm1PYmO')
# gd_1_wqp_pull_out <- scipiper:::gd_locate_file('1_wqp_pull/out') %>% pull(id) %>% tail(1)
# out_files <- drive_ls(gd_1_wqp_pull_out)
obs_files_avail <- dir('data/WQPTempObs190412203340/', full.names = TRUE)
obs_files_wanted <- sort(unique(file.path('data/WQPTempObs190412203340', sprintf('%s.feather', del_wqp_sites$PullTask))))
obs_files_missing <- setdiff(obs_files_wanted, obs_files_avail)
if(length(obs_files_missing) > 0) {
  warning('These wanted data files are missing: ', paste0('\n  ', basename(obs_files_missing), collapse=''))
}
obs_files_use <- intersect(obs_files_wanted, obs_files_avail)
obs_data <- bind_rows(lapply(obs_files_use, function(obs_file) {
  message(basename(obs_file))
  tryCatch({
    feather::read_feather(obs_file) %>%
      dplyr::filter(MonitoringLocationIdentifier %in% del_wqp_sites$MonitoringLocationIdentifier) %>%
      dplyr::select(MonitoringLocationIdentifier, ActivityStartDate, ResultMeasureValue)
  }, error=function(e) {
    message(e$message, appendLF = TRUE)
    NULL
  })
}))
saveRDS(obs_data, 'out/obs_data.rds')
obs_data <- readRDS('out/obs_data.rds')
# Compare the number of results we expected to those we were actually able to obtain
plot_dates <- as.Date(c('2007-01-01', '2010-01-01'))
obs_summary <- obs_data %>%
  group_by(MonitoringLocationIdentifier) %>%
  summarize(
    resultCountEmpirical = n(), # ResultMeasureValue seems never to be NA
    windowResultCount = length(which(dplyr::between(unique(ActivityStartDate), plot_dates[1], plot_dates[2])))) %>%
  left_join(dplyr::rename(sf::st_drop_geometry(del_wqp_sites), resultCountExpected = resultCount), by='MonitoringLocationIdentifier')
# We're missing some data for some sites:
nrow(obs_summary) # 3376 sites with some observations available
obs_summary %>% filter(resultCountEmpirical != resultCountExpected) %>% nrow() # 705 sites with some observations missing
sum(obs_summary$resultCountEmpirical - obs_summary$resultCountExpected) # 25872 observations missing

#### Plot PRMS predictions vs WQP obs ####

# Pick some sites to visualize
focus_sites <- filter(obs_summary, windowResultCount > 100, dist_to_reach_km < 1)
focus_obs <- filter(obs_data, MonitoringLocationIdentifier %in% focus_sites$MonitoringLocationIdentifier) %>%
  dplyr::select(MonitoringLocationIdentifier, ActivityStartDate, temp_obs=ResultMeasureValue) %>%
  left_join(dplyr::select(focus_sites, MonitoringLocationIdentifier, seg_id_nat), by='MonitoringLocationIdentifier') %>%
  group_by(MonitoringLocationIdentifier, seg_id_nat, ActivityStartDate) %>%
  summarize(n_obs = n(), temp_obs_sd = sd(temp_obs), temp_obs = mean(temp_obs)) %>%
  ungroup() %>%
  dplyr::select(date = ActivityStartDate, seg_id_nat, temp_obs, temp_obs_sd, n_obs)
dplyr::filter(focus_obs, dplyr::between(date, plot_dates[1], plot_dates[2])) %>%
  group_by(seg_id_nat) %>%
  filter(!is.na(temp_obs)) %>%
  summarize(n())

# Read the PRMS temperature predictions
prms_preds <- read_csv('data/delaware_stream_temp_by_segment/seg_tave_water.csv') %>%
  gather("seg_id_nat", "temp_pred", -date) %>% 
  mutate(
    seg_id_nat = as.integer(seg_id_nat),
    temp_pred = ifelse(temp_pred %in% c(-99.9, -98.9), NA, temp_pred))
focus_preds <- filter(prms_preds, seg_id_nat %in% focus_sites$seg_id_nat)

# Match up predictions and observations
focus_pred_obs <- full_join(focus_preds, focus_obs, by=c('seg_id_nat', 'date'))

# Plot preds and obs vs time
dplyr::filter(focus_pred_obs, dplyr::between(date, plot_dates[1], plot_dates[2])) %>%
  ggplot(aes(x=date)) +
  geom_line(aes(y=temp_pred), color='#7b3294') +
  geom_point(aes(y=temp_obs), color='#008837', size=1, alpha=0.5, na.rm=TRUE) +
  theme_bw() +
  facet_wrap(~ seg_id_nat) +
  xlab('Date') + ylab('Temperature (deg C)')
ggsave('out/predobs_v_time.png', width=7, height=7)

# Plot preds vs obs
focus_pred_obs %>%
  mutate(Date = ifelse(dplyr::between(date, plot_dates[1], plot_dates[2]), paste(format(plot_dates, '%Y'), collapse='-'), 'Other')) %>%
  ggplot(aes(x=temp_obs, y=temp_pred)) +
  geom_abline() +
  geom_point(aes(color=Date), size=0.8, alpha=0.5, na.rm=TRUE) +
  theme_bw() +
  facet_wrap(~ seg_id_nat) +
  xlab('Observed (deg C)') + ylab('Predicted (deg C)')
ggsave('out/pred_v_time.png', width=8, height=7)
