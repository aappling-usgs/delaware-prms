library(tidyverse)

# Load DRB spatial data
crosswalk_site_reach <- readRDS('out/crosswalk_site_reach.rds')
drb_net <- readRDS('out/network_full.rds')
umn_net <- readRDS('out/network_subset.rds')
data_snapshot_path <- 'data/190926-2wp-temp-observations'

# Subset the crosswalk to sites with decent matches to reaches (see
# map_sites_to_reaches for exploration). Also simplify because we don't need all
# this information below
drb_sites <- crosswalk_site_reach %>%
  mutate(
    bird_filtered = bird_dist_to_subseg_m < 250,
    fish_filtered = abs(fish_dist_to_outlet_m) < 5000) %>%
  filter(bird_filtered, fish_filtered) %>%
  st_drop_geometry() %>%
  select(site_id, subseg_id, seg_id_nat) %>%
  distinct() # slightly alarming that we need this. st_drop_geometry(crosswalk_site_reach) without select() has no duplicates
umn_sites <- filter(drb_sites, subseg_id %in% umn_net$edges$subseg_id)

# Download the inventories of temperature observation sites
# library(googledrive)
# drive_ls(as_id('1RqKhyYi4tHvD1fli6LKocw8xHwm1PYmO'), pattern='wqp_data.rds') # 1_wqp_pull/out contents
# drive_download(as_id('1zmqQe4HJfQylB7BPbZCffxTCE7MTJB2s'), file.path(data_snapshot_path, '1_wqp_pull/out/wqp_data.rds'), overwrite=TRUE)
# drive_ls(as_id('15bS-vSTwFAuqkKZplZ7ljxKW2BgNeZ_K')) # 2_nwis_pull/out contents
# drive_download(as_id('1irKTun7NdVjZQeP0aHamWEEAhmltjBii'), file.path(data_snapshot_path, '2_nwis_pull/out/nwis_dv_data.rds'), overwrite=TRUE)
# drive_download(as_id('1im_Xjnw0b3ajfTLGjskMTMTPyYgDuD7C'), file.path(data_snapshot_path, '2_nwis_pull/out/nwis_uv_data.rds'), overwrite=TRUE)

#### Collect actual observations ####

dat_wqp <- readRDS(file.path(data_snapshot_path, '1_wqp_pull/out/wqp_data.rds'))
dat_nwis_dv <- readRDS(file.path(data_snapshot_path, '2_nwis_pull/out/nwis_dv_data.rds'))
dat_nwis_uv <- readRDS(file.path(data_snapshot_path, '2_nwis_pull/out/nwis_uv_data.rds'))

# Do initial filter (to full list of acceptable DRB sites) and combine
dat_wqp_drb <- dat_wqp %>%
  select(site_id = MonitoringLocationIdentifier, date = ActivityStartDate, temp_C = ResultMeasureValue, `ResultMeasure/MeasureUnitCode`) %>%
  filter(site_id %in% drb_sites$site_id) %>%
  filter(trimws(`ResultMeasure/MeasureUnitCode`) == 'deg C', temp_C >= 0, temp_C < 35)
dat_nwis_dv_drb <- dat_nwis_dv %>% 
  mutate(site_id = sprintf('USGS-%s', site_no), date=as.Date(dateTime)) %>%
  filter(site_id %in% drb_sites$site_id) %>%
  select(site_id, date, temp_C=temp_value) %>%
  filter(temp_C >= 0, temp_C < 35)
dat_nwis_uv_drb <- dat_nwis_uv %>% 
  mutate(site_id = sprintf('USGS-%s', site_no), date=as.Date(dateTime)) %>%
  select(site_id, date, temp_C=temp_value) %>%
  filter(site_id %in% drb_sites$site_id) %>%
  filter(temp_C >= 0, temp_C < 35) %>%
  group_by(site_id, date) %>%
  summarize(temp_C = mean(temp_C, na.rm=TRUE)) %>%
  ungroup()
dat_drb_dups <- bind_rows(dat_wqp_drb, dat_nwis_dv_drb, dat_nwis_uv_drb)

# check for and handle obvious site-date duplicates
dat_drb_dups %>% select(site_id, date) %>% nrow # 383038
dat_drb_dups %>% select(site_id, date) %>% distinct %>% nrow # 347005
dat_drb_sites <- dat_drb_dups %>%
  group_by(site_id, date) %>%
  summarize(temp_C = mean(temp_C)) %>%
  ungroup()

# map sites to reaches and collapse again
dat_drb <- dat_drb_sites %>%
  left_join(drb_sites, by='site_id') %>%
  group_by(subseg_id, seg_id_nat, date) %>%
  summarize(temp_C = mean(temp_C)) %>%
  ungroup()

dat_umn <- dat_drb %>%
  filter(subseg_id %in% umn_net$edges$subseg_id)
# confirm that there are no subsegs that lack a seg_id_nat (those I added to the
# network), which means we can use seg_id_nat as a unique reach identifier in
# this subnetwork
dat_umn %>% select(subseg_id, seg_id_nat) %>% distinct %>% group_by(seg_id_nat) %>% summarize(n_subseg=length(subseg_id)) %>% filter(n_subseg != 1)
dat_umn <- dat_umn %>%
  select(-subseg_id)

# Save the temperature data for each of the two networks. I will save as csv
# rather than a python format because it's a mixed-medium data.frame, so pandas,
# and I don't know the syntax at the moment
readr::write_csv(dat_drb, path='out/obs_temp_full.csv')
readr::write_csv(dat_umn, path='out/obs_temp_subset.csv')

# Read in the discharge data
ddat_wide <- readr::read_csv('data/drb_discharge_daily_dv.csv', col_types=cols(.default=col_double(), datetime=col_date(format='')))
ddat_drb <- ddat_wide %>%
  gather(site_no, discharge_cfs, -datetime) %>%
  mutate(
    site_id = sprintf('USGS-%s', site_no),
    discharge_cms = discharge_cfs / 35.314666) %>%
  select(-site_no, -discharge_cfs) %>%
  rename(date=datetime) %>%
  filter(!is.na(discharge_cms)) %>%
  left_join(drb_sites, by='site_id') %>%
  group_by(seg_id_nat, subseg_id, date) %>% 
  summarize(discharge_cms = mean(discharge_cms)) %>%
  ungroup()

ddat_umn <- ddat_drb %>%
  filter(subseg_id %in% umn_net$edges$subseg_id) %>%
  select(-subseg_id)

readr::write_csv(ddat_drb, path='out/obs_flow_full.csv')
readr::write_csv(ddat_umn, path='out/obs_flow_subset.csv')
