# Label variables (youth), (parent)
# peta: 'growth in height' (pds_ht2_y), (pds_1_p)
# petb: 'growth of body hair' (pds_bdyhair_y), (pds_2_p)
# petc: 'noticed skin changes' (pds_skin2_y), (pds_3_p)
# petd: 'breasts begun to grow / deepening of voice' (pds_f4_2_y / pds_m4_y), (pds_f4_p / pds_m4_p)
# mpete: 'male grow hair on face' (pds_m5_y), (pds_m4_p)
# fpete: 'female begun to menstruate' (pds_f5_y), (pds_f5b_p)

select_participant_by_obs <- function(data) {
  # step 1: count the number of observations for each participant ('id')
  data_with_counts <- data %>%
    add_count(id, name = "n_obs")
  
  # step 2: identify the single participant ('id') to keep from each family
  # this is the one with the maximum number of observations.
  # slice_max handles ties by picking the first one it encounters.
  ids_to_keep <- data_with_counts %>%
    # get unique combinations to avoid redundant rows
    distinct(rel_family_id, id, n_obs) %>%
    # for each family...
    group_by(rel_family_id) %>%
    # ...find the row with the highest 'n_obs'
    slice_max(order_by = n_obs, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    # pull out just the participant ids
    pull(id)
  
  # step 3: filter the original dataframe to keep all rows for the selected ids
  data %>%
    filter(id %in% ids_to_keep)
}


# recode variables for females
raw_PDS_f <- raw_PDS_f %>%
  mutate(petbf = ifelse(petb == 1, 1,
                        ifelse(petb == 2, 2,
                               ifelse(petb == 3, 4,
                                      ifelse(petb == 4, 5, NA))))) %>%
  mutate(petcf = ifelse(petc == 1, 1,
                        ifelse(petc == 2, 2,
                               ifelse(petc == 3, 4,
                                      ifelse(petc == 4, 5, NA)))))

raw_PDS_f$adrenf <- rowMeans(cbind(raw_PDS_f$petbf, raw_PDS_f$petcf), na.rm = TRUE)
raw_PDS_f$adrenf2 <- raw_PDS_f$adrenf
raw_PDS_f$adrenf2[raw_PDS_f$petb == 1 & raw_PDS_f$adrenf == 1.5] <- 1
raw_PDS_f$adrenf2[raw_PDS_f$petb == 2 & raw_PDS_f$adrenf == 1.5] <- 2
raw_PDS_f$adrenf2[raw_PDS_f$adrenf == 2.5] <- 3
raw_PDS_f$adrenf2[raw_PDS_f$adrenf == 3.5] <- 4
raw_PDS_f$adrenf2[raw_PDS_f$adrenf == 4.5] <- 5
raw_PDS_f$adrenf2[raw_PDS_f$adrenf == 5.5] <- 5

raw_PDS_f <- raw_PDS_f %>%
  mutate(petaf = ifelse(peta == 1, 1,
                        ifelse(peta == 2, 2,
                               ifelse(peta == 3, 3,
                                      ifelse(peta == 4, 5, NA))))) %>%
  mutate(petdf = ifelse(petd == 1, 1,
                        ifelse(petd == 2, 3,
                               ifelse(petd == 3, 4,
                                      ifelse(petd == 4, 5, NA))))) %>%
  mutate(petef = ifelse(fpete == 0, 1,
                        ifelse(fpete == 1, 5, NA)))

raw_PDS_f$gonadf <- rowMeans(cbind(raw_PDS_f$petaf, raw_PDS_f$petdf), na.rm = TRUE)
raw_PDS_f$gonadf2 <- raw_PDS_f$gonadf
raw_PDS_f$gonadf2[raw_PDS_f$gonadf == 1 & raw_PDS_f$petef == 1] <- 1
raw_PDS_f$gonadf2[raw_PDS_f$gonadf == 1.5 & raw_PDS_f$petef == 1] <- 1
raw_PDS_f$gonadf2[raw_PDS_f$gonadf == 2 & raw_PDS_f$petef == 1] <- 2
raw_PDS_f$gonadf2[raw_PDS_f$gonadf == 2.5 & raw_PDS_f$petef == 1] <- 2
raw_PDS_f$gonadf2[raw_PDS_f$gonadf == 3 & raw_PDS_f$petef == 1] <- 3
raw_PDS_f$gonadf2[raw_PDS_f$gonadf == 3.5 & raw_PDS_f$petef == 1] <- 3
raw_PDS_f$gonadf2[raw_PDS_f$gonadf == 4 & raw_PDS_f$petef == 1] <- 3
raw_PDS_f$gonadf2[raw_PDS_f$gonadf == 4.5 & raw_PDS_f$petef == 1] <- 4
raw_PDS_f$gonadf2[raw_PDS_f$gonadf == 5 & raw_PDS_f$petef == 1] <- 4
raw_PDS_f$gonadf2[raw_PDS_f$gonadf == 1 & raw_PDS_f$petef == 5] <- 2
raw_PDS_f$gonadf2[raw_PDS_f$gonadf == 1.5 & raw_PDS_f$petef == 5] <- 3
raw_PDS_f$gonadf2[raw_PDS_f$gonadf == 2 & raw_PDS_f$petef == 5] <- 4
raw_PDS_f$gonadf2[raw_PDS_f$gonadf == 2.5 & raw_PDS_f$petef == 5] <- 4
raw_PDS_f$gonadf2[raw_PDS_f$gonadf == 3 & raw_PDS_f$petef == 5] <- 4
raw_PDS_f$gonadf2[raw_PDS_f$gonadf == 3.5 & raw_PDS_f$petef == 5] <- 5
raw_PDS_f$gonadf2[raw_PDS_f$gonadf == 4 & raw_PDS_f$petef == 5] <- 5
raw_PDS_f$gonadf2[raw_PDS_f$gonadf == 4.5 & raw_PDS_f$petef == 5] <- 5
raw_PDS_f$gonadf2[raw_PDS_f$gonadf == 5 & raw_PDS_f$petef == 5] <- 5

raw_PDS_f$PDSS <- rowMeans(cbind(raw_PDS_f$gonadf2, raw_PDS_f$adrenf2), na.rm = TRUE)

# recode variables for males
raw_PDS_m <- raw_PDS_m %>%
  mutate(petbm = ifelse(petb == 1, 1,
                        ifelse(petb == 2, 2,
                               ifelse(petb == 3, 4,
                                      ifelse(petb == 4, 5, NA))))) %>%
  mutate(petcm = ifelse(petc == 1, 1,
                        ifelse(petc == 2, 2,
                               ifelse(petc == 3, 3,
                                      ifelse(petc == 4, 4, NA)))))

raw_PDS_m$adrenm <- rowMeans(cbind(raw_PDS_m$petbm, raw_PDS_m$petcm), na.rm = TRUE)
raw_PDS_m$adrenm2 <- raw_PDS_m$adrenm
raw_PDS_m$adrenm2[raw_PDS_m$adrenm == 1] <- 1
raw_PDS_m$adrenm2[raw_PDS_m$adrenm == 1.5 & raw_PDS_m$petcm == 1] <- 1
raw_PDS_m$adrenm2[raw_PDS_m$adrenm == 1.5 & raw_PDS_m$petcm == 2] <- 2
raw_PDS_m$adrenm2[raw_PDS_m$adrenm == 2.5 & raw_PDS_m$petbm != 4] <- 2
raw_PDS_m$adrenm2[raw_PDS_m$adrenm == 2.5 & raw_PDS_m$petbm == 4] <- 3
raw_PDS_m$adrenm2[raw_PDS_m$adrenm == 3.5] <- 4
raw_PDS_m$adrenm2[raw_PDS_m$adrenm == 4.5] <- 5
raw_PDS_m$adrenm2[raw_PDS_m$adrenm == 5.5] <- 5

raw_PDS_m <- raw_PDS_m %>%
  mutate(petam = ifelse(peta == 1, 1,
                        ifelse(peta == 2, 3,
                               ifelse(peta == 3, 4,
                                      ifelse(peta == 4, 5, NA))))) %>%
  mutate(petdm = ifelse(petd == 1, 1,
                        ifelse(petd == 2, 2,
                               ifelse(petd == 3, 3,
                                      ifelse(petd == 4, 5, NA)))))

raw_PDS_m$gonadm <- rowMeans(cbind(raw_PDS_m$petam, raw_PDS_m$petdm), na.rm = TRUE)
raw_PDS_m$gonadm2 <- raw_PDS_m$gonadm
raw_PDS_m$gonadm2[raw_PDS_m$gonadm == 1 & raw_PDS_m$mpete == 1] <- 1
raw_PDS_m$gonadm2[raw_PDS_m$gonadm == 1 & raw_PDS_m$mpete > 1] <- 2
raw_PDS_m$gonadm2[raw_PDS_m$gonadm == 1.5 & raw_PDS_m$mpete == 1] <- 1
raw_PDS_m$gonadm2[raw_PDS_m$gonadm == 1.5 & raw_PDS_m$mpete > 1] <- 2
raw_PDS_m$gonadm2[raw_PDS_m$gonadm == 2 & raw_PDS_m$mpete == 1 & raw_PDS_m$petd == 1] <- 1
raw_PDS_m$gonadm2[raw_PDS_m$gonadm == 2 & raw_PDS_m$mpete == 1 & raw_PDS_m$petd > 1] <- 2
raw_PDS_m$gonadm2[raw_PDS_m$gonadm == 2 & raw_PDS_m$mpete > 1] <- 3
raw_PDS_m$gonadm2[raw_PDS_m$gonadm == 2.5 & raw_PDS_m$mpete == 1] <- 2
raw_PDS_m$gonadm2[raw_PDS_m$gonadm == 2.5 & raw_PDS_m$mpete > 1] <- 3
raw_PDS_m$gonadm2[raw_PDS_m$gonadm == 3] <- 3
raw_PDS_m$gonadm2[raw_PDS_m$gonadm == 3.5 & raw_PDS_m$mpete == 2] <- 4
raw_PDS_m$gonadm2[raw_PDS_m$gonadm == 3.5 & raw_PDS_m$mpete == 1] <- 4
raw_PDS_m$gonadm2[raw_PDS_m$gonadm == 3.5 & raw_PDS_m$mpete > 2] <- 5
raw_PDS_m$gonadm2[raw_PDS_m$gonadm == 4 & raw_PDS_m$mpete == 1] <- 4
raw_PDS_m$gonadm2[raw_PDS_m$gonadm == 4 & raw_PDS_m$mpete == 2] <- 4
raw_PDS_m$gonadm2[raw_PDS_m$gonadm == 4 & raw_PDS_m$mpete > 2] <- 5
raw_PDS_m$gonadm2[raw_PDS_m$gonadm > 4] <- 5

## creating scored dataframe

PDS_f <- raw_PDS_f
PDS_m <- raw_PDS_m

PDS_f$PDSS <- rowMeans(cbind(PDS_f$gonadf2, PDS_f$adrenf2), na.rm = TRUE)
PDS_m$PDSS <- rowMeans(cbind(PDS_m$gonadm2, PDS_m$adrenm2), na.rm = TRUE)


track <- read.csv(paste0(data_root, "abcd-general/ab_g_dyn.csv")) %>% 
  rename("id" = "participant_id",
         "wave" = "session_id",
         "rel_family_id" = "ab_g_stc__design_id__fam") %>% 
  select(id, wave, site, rel_family_id)

PDS_f <- left_join(PDS_f, track, by = c("id", "wave"))
PDS_m <- left_join(PDS_m, track, by = c("id", "wave"))

PDS_f <- PDS_f %>%
  group_by(id) %>%
  mutate(rel_family_id = ifelse(is.na(rel_family_id), first(rel_family_id), rel_family_id)) %>%   
  ungroup()

PDS_m <- PDS_m %>%
  group_by(id) %>%
  mutate(rel_family_id = ifelse(is.na(rel_family_id), first(rel_family_id), rel_family_id)) %>%   
  ungroup()

write.csv(PDS_f, file = paste0(data_root, "physical-health/puberty/parentyouth_tannerstages_f.csv"))
write.csv(PDS_m, file = paste0(data_root, "physical-health/puberty/parentyouth_tannerstages_m.csv"))

# now, apply the function to each dataframe
PDS_f_filtered <- select_participant_by_obs(PDS_f)
PDS_m_filtered <- select_participant_by_obs(PDS_m)

write.csv(PDS_f_filtered, file = paste0(data_root, "physical-health/puberty/filtered_combined_tannerstages_f.csv"))
write.csv(PDS_m_filtered, file = paste0(data_root, "physical-health/puberty/filtered_combined_tannerstages_m.csv"))

# recode variables for females
raw_PDS_f_y <- raw_PDS_f_y %>%
  mutate(petbf = ifelse(petb == 1, 1,
                        ifelse(petb == 2, 2,
                               ifelse(petb == 3, 4,
                                      ifelse(petb == 4, 5, NA))))) %>%
  mutate(petcf = ifelse(petc == 1, 1,
                        ifelse(petc == 2, 2,
                               ifelse(petc == 3, 4,
                                      ifelse(petc == 4, 5, NA)))))

raw_PDS_f_y$adrenf <- rowMeans(cbind(raw_PDS_f_y$petbf, raw_PDS_f_y$petcf), na.rm = TRUE)
raw_PDS_f_y$adrenf2 <- raw_PDS_f_y$adrenf
raw_PDS_f_y$adrenf2[raw_PDS_f_y$petb == 1 & raw_PDS_f_y$adrenf == 1.5] <- 1
raw_PDS_f_y$adrenf2[raw_PDS_f_y$petb == 2 & raw_PDS_f_y$adrenf == 1.5] <- 2
raw_PDS_f_y$adrenf2[raw_PDS_f_y$adrenf == 2.5] <- 3
raw_PDS_f_y$adrenf2[raw_PDS_f_y$adrenf == 3.5] <- 4
raw_PDS_f_y$adrenf2[raw_PDS_f_y$adrenf == 4.5] <- 5
raw_PDS_f_y$adrenf2[raw_PDS_f_y$adrenf == 5.5] <- 5

raw_PDS_f_y <- raw_PDS_f_y %>%
  mutate(petaf = ifelse(peta == 1, 1,
                        ifelse(peta == 2, 2,
                               ifelse(peta == 3, 3,
                                      ifelse(peta == 4, 5, NA))))) %>%
  mutate(petdf = ifelse(petd == 1, 1,
                        ifelse(petd == 2, 3,
                               ifelse(petd == 3, 4,
                                      ifelse(petd == 4, 5, NA))))) %>%
  mutate(petef = ifelse(fpete == 0, 1,
                        ifelse(fpete == 1, 5, NA)))

raw_PDS_f_y$gonadf <- rowMeans(cbind(raw_PDS_f_y$petaf, raw_PDS_f_y$petdf), na.rm = TRUE)
raw_PDS_f_y$gonadf2 <- raw_PDS_f_y$gonadf
raw_PDS_f_y$gonadf2[raw_PDS_f_y$gonadf == 1 & raw_PDS_f_y$petef == 1] <- 1
raw_PDS_f_y$gonadf2[raw_PDS_f_y$gonadf == 1.5 & raw_PDS_f_y$petef == 1] <- 1
raw_PDS_f_y$gonadf2[raw_PDS_f_y$gonadf == 2 & raw_PDS_f_y$petef == 1] <- 2
raw_PDS_f_y$gonadf2[raw_PDS_f_y$gonadf == 2.5 & raw_PDS_f_y$petef == 1] <- 2
raw_PDS_f_y$gonadf2[raw_PDS_f_y$gonadf == 3 & raw_PDS_f_y$petef == 1] <- 3
raw_PDS_f_y$gonadf2[raw_PDS_f_y$gonadf == 3.5 & raw_PDS_f_y$petef == 1] <- 3
raw_PDS_f_y$gonadf2[raw_PDS_f_y$gonadf == 4 & raw_PDS_f_y$petef == 1] <- 3
raw_PDS_f_y$gonadf2[raw_PDS_f_y$gonadf == 4.5 & raw_PDS_f_y$petef == 1] <- 4
raw_PDS_f_y$gonadf2[raw_PDS_f_y$gonadf == 5 & raw_PDS_f_y$petef == 1] <- 4
raw_PDS_f_y$gonadf2[raw_PDS_f_y$gonadf == 1 & raw_PDS_f_y$petef == 5] <- 2
raw_PDS_f_y$gonadf2[raw_PDS_f_y$gonadf == 1.5 & raw_PDS_f_y$petef == 5] <- 3
raw_PDS_f_y$gonadf2[raw_PDS_f_y$gonadf == 2 & raw_PDS_f_y$petef == 5] <- 4
raw_PDS_f_y$gonadf2[raw_PDS_f_y$gonadf == 2.5 & raw_PDS_f_y$petef == 5] <- 4
raw_PDS_f_y$gonadf2[raw_PDS_f_y$gonadf == 3 & raw_PDS_f_y$petef == 5] <- 4
raw_PDS_f_y$gonadf2[raw_PDS_f_y$gonadf == 3.5 & raw_PDS_f_y$petef == 5] <- 5
raw_PDS_f_y$gonadf2[raw_PDS_f_y$gonadf == 4 & raw_PDS_f_y$petef == 5] <- 5
raw_PDS_f_y$gonadf2[raw_PDS_f_y$gonadf == 4.5 & raw_PDS_f_y$petef == 5] <- 5
raw_PDS_f_y$gonadf2[raw_PDS_f_y$gonadf == 5 & raw_PDS_f_y$petef == 5] <- 5

raw_PDS_f_y$PDSS <- rowMeans(cbind(raw_PDS_f_y$gonadf2, raw_PDS_f_y$adrenf2), na.rm = TRUE)

# recode variables for males
raw_PDS_m_y <- raw_PDS_m_y %>%
  mutate(petbm = ifelse(petb == 1, 1,
                        ifelse(petb == 2, 2,
                               ifelse(petb == 3, 4,
                                      ifelse(petb == 4, 5, NA))))) %>%
  mutate(petcm = ifelse(petc == 1, 1,
                        ifelse(petc == 2, 2,
                               ifelse(petc == 3, 3,
                                      ifelse(petc == 4, 4, NA)))))

raw_PDS_m_y$adrenm <- rowMeans(cbind(raw_PDS_m_y$petbm, raw_PDS_m_y$petcm), na.rm = TRUE)
raw_PDS_m_y$adrenm2 <- raw_PDS_m_y$adrenm
raw_PDS_m_y$adrenm2[raw_PDS_m_y$adrenm == 1] <- 1
raw_PDS_m_y$adrenm2[raw_PDS_m_y$adrenm == 1.5 & raw_PDS_m_y$petcm == 1] <- 1
raw_PDS_m_y$adrenm2[raw_PDS_m_y$adrenm == 1.5 & raw_PDS_m_y$petcm == 2] <- 2
raw_PDS_m_y$adrenm2[raw_PDS_m_y$adrenm == 2.5 & raw_PDS_m_y$petbm != 4] <- 2
raw_PDS_m_y$adrenm2[raw_PDS_m_y$adrenm == 2.5 & raw_PDS_m_y$petbm == 4] <- 3
raw_PDS_m_y$adrenm2[raw_PDS_m_y$adrenm == 3.5] <- 4
raw_PDS_m_y$adrenm2[raw_PDS_m_y$adrenm == 4.5] <- 5
raw_PDS_m_y$adrenm2[raw_PDS_m_y$adrenm == 5.5] <- 5

raw_PDS_m_y <- raw_PDS_m_y %>%
  mutate(petam = ifelse(peta == 1, 1,
                        ifelse(peta == 2, 3,
                               ifelse(peta == 3, 4,
                                      ifelse(peta == 4, 5, NA))))) %>%
  mutate(petdm = ifelse(petd == 1, 1,
                        ifelse(petd == 2, 2,
                               ifelse(petd == 3, 3,
                                      ifelse(petd == 4, 5, NA)))))

raw_PDS_m_y$gonadm <- rowMeans(cbind(raw_PDS_m_y$petam, raw_PDS_m_y$petdm), na.rm = TRUE)
raw_PDS_m_y$gonadm2 <- raw_PDS_m_y$gonadm
raw_PDS_m_y$gonadm2[raw_PDS_m_y$gonadm == 1 & raw_PDS_m_y$mpete == 1] <- 1
raw_PDS_m_y$gonadm2[raw_PDS_m_y$gonadm == 1 & raw_PDS_m_y$mpete > 1] <- 2
raw_PDS_m_y$gonadm2[raw_PDS_m_y$gonadm == 1.5 & raw_PDS_m_y$mpete == 1] <- 1
raw_PDS_m_y$gonadm2[raw_PDS_m_y$gonadm == 1.5 & raw_PDS_m_y$mpete > 1] <- 2
raw_PDS_m_y$gonadm2[raw_PDS_m_y$gonadm == 2 & raw_PDS_m_y$mpete == 1 & raw_PDS_m_y$petd == 1] <- 1
raw_PDS_m_y$gonadm2[raw_PDS_m_y$gonadm == 2 & raw_PDS_m_y$mpete == 1 & raw_PDS_m_y$petd > 1] <- 2
raw_PDS_m_y$gonadm2[raw_PDS_m_y$gonadm == 2 & raw_PDS_m_y$mpete > 1] <- 3
raw_PDS_m_y$gonadm2[raw_PDS_m_y$gonadm == 2.5 & raw_PDS_m_y$mpete == 1] <- 2
raw_PDS_m_y$gonadm2[raw_PDS_m_y$gonadm == 2.5 & raw_PDS_m_y$mpete > 1] <- 3
raw_PDS_m_y$gonadm2[raw_PDS_m_y$gonadm == 3] <- 3
raw_PDS_m_y$gonadm2[raw_PDS_m_y$gonadm == 3.5 & raw_PDS_m_y$mpete == 2] <- 4
raw_PDS_m_y$gonadm2[raw_PDS_m_y$gonadm == 3.5 & raw_PDS_m_y$mpete == 1] <- 4
raw_PDS_m_y$gonadm2[raw_PDS_m_y$gonadm == 3.5 & raw_PDS_m_y$mpete > 2] <- 5
raw_PDS_m_y$gonadm2[raw_PDS_m_y$gonadm == 4 & raw_PDS_m_y$mpete == 1] <- 4
raw_PDS_m_y$gonadm2[raw_PDS_m_y$gonadm == 4 & raw_PDS_m_y$mpete == 2] <- 4
raw_PDS_m_y$gonadm2[raw_PDS_m_y$gonadm == 4 & raw_PDS_m_y$mpete > 2] <- 5
raw_PDS_m_y$gonadm2[raw_PDS_m_y$gonadm > 4] <- 5

## creating scored dataframe

PDS_f_y <- raw_PDS_f_y
PDS_m_y <- raw_PDS_m_y

PDS_f_y$PDSS <- rowMeans(cbind(PDS_f_y$gonadf2, PDS_f_y$adrenf2), na.rm = TRUE)
PDS_m_y$PDSS <- rowMeans(cbind(PDS_m_y$gonadm2, PDS_m_y$adrenm2), na.rm = TRUE)


track <- read.csv(paste0(data_root, "abcd-general/ab_g_dyn.csv")) %>% 
  rename("id" = "participant_id",
         "wave" = "session_id",
         "rel_family_id" = "ab_g_stc__design_id__fam") %>% 
  select(id, wave, site, rel_family_id)

PDS_f_y <- left_join(PDS_f_y, track, by = c("id", "wave"))
PDS_m_y <- left_join(PDS_m_y, track, by = c("id", "wave"))

PDS_f_y <- PDS_f_y %>%
  group_by(id) %>%
  mutate(rel_family_id = ifelse(is.na(rel_family_id), first(rel_family_id), rel_family_id)) %>%   
  ungroup()

PDS_m_y <- PDS_m_y %>%
  group_by(id) %>%
  mutate(rel_family_id = ifelse(is.na(rel_family_id), first(rel_family_id), rel_family_id)) %>%   
  ungroup()

write.csv(PDS_f_y, file = paste0(data_root, "physical-health/puberty/youth_tannerstages_f.csv"))
write.csv(PDS_m_y, file = paste0(data_root, "physical-health/puberty/youth_tannerstages_m.csv"))

# now, apply the function to each dataframe
PDS_f_y_filtered <- select_participant_by_obs(PDS_f_y)
PDS_m_y_filtered <- select_participant_by_obs(PDS_m_y)

write.csv(PDS_f_y_filtered, file = paste0(data_root, "physical-health/puberty/filtered_youth_tannerstages_f.csv"))
write.csv(PDS_m_y_filtered, file = paste0(data_root, "physical-health/puberty/filtered_youth_tannerstages_m.csv"))

# recode variables for females
raw_PDS_f_p <- raw_PDS_f_p %>%
  mutate(petbf = ifelse(petb == 1, 1,
                        ifelse(petb == 2, 2,
                               ifelse(petb == 3, 4,
                                      ifelse(petb == 4, 5, NA))))) %>%
  mutate(petcf = ifelse(petc == 1, 1,
                        ifelse(petc == 2, 2,
                               ifelse(petc == 3, 4,
                                      ifelse(petc == 4, 5, NA)))))

raw_PDS_f_p$adrenf <- rowMeans(cbind(raw_PDS_f_p$petbf, raw_PDS_f_p$petcf), na.rm = TRUE)
raw_PDS_f_p$adrenf2 <- raw_PDS_f_p$adrenf
raw_PDS_f_p$adrenf2[raw_PDS_f_p$petb == 1 & raw_PDS_f_p$adrenf == 1.5] <- 1
raw_PDS_f_p$adrenf2[raw_PDS_f_p$petb == 2 & raw_PDS_f_p$adrenf == 1.5] <- 2
raw_PDS_f_p$adrenf2[raw_PDS_f_p$adrenf == 2.5] <- 3
raw_PDS_f_p$adrenf2[raw_PDS_f_p$adrenf == 3.5] <- 4
raw_PDS_f_p$adrenf2[raw_PDS_f_p$adrenf == 4.5] <- 5
raw_PDS_f_p$adrenf2[raw_PDS_f_p$adrenf == 5.5] <- 5

raw_PDS_f_p <- raw_PDS_f_p %>%
  mutate(petaf = ifelse(peta == 1, 1,
                        ifelse(peta == 2, 2,
                               ifelse(peta == 3, 3,
                                      ifelse(peta == 4, 5, NA))))) %>%
  mutate(petdf = ifelse(petd == 1, 1,
                        ifelse(petd == 2, 3,
                               ifelse(petd == 3, 4,
                                      ifelse(petd == 4, 5, NA))))) %>%
  mutate(petef = ifelse(fpete == 0, 1,
                        ifelse(fpete == 1, 5, NA)))

raw_PDS_f_p$gonadf <- rowMeans(cbind(raw_PDS_f_p$petaf, raw_PDS_f_p$petdf), na.rm = TRUE)
raw_PDS_f_p$gonadf2 <- raw_PDS_f_p$gonadf
raw_PDS_f_p$gonadf2[raw_PDS_f_p$gonadf == 1 & raw_PDS_f_p$petef == 1] <- 1
raw_PDS_f_p$gonadf2[raw_PDS_f_p$gonadf == 1.5 & raw_PDS_f_p$petef == 1] <- 1
raw_PDS_f_p$gonadf2[raw_PDS_f_p$gonadf == 2 & raw_PDS_f_p$petef == 1] <- 2
raw_PDS_f_p$gonadf2[raw_PDS_f_p$gonadf == 2.5 & raw_PDS_f_p$petef == 1] <- 2
raw_PDS_f_p$gonadf2[raw_PDS_f_p$gonadf == 3 & raw_PDS_f_p$petef == 1] <- 3
raw_PDS_f_p$gonadf2[raw_PDS_f_p$gonadf == 3.5 & raw_PDS_f_p$petef == 1] <- 3
raw_PDS_f_p$gonadf2[raw_PDS_f_p$gonadf == 4 & raw_PDS_f_p$petef == 1] <- 3
raw_PDS_f_p$gonadf2[raw_PDS_f_p$gonadf == 4.5 & raw_PDS_f_p$petef == 1] <- 4
raw_PDS_f_p$gonadf2[raw_PDS_f_p$gonadf == 5 & raw_PDS_f_p$petef == 1] <- 4
raw_PDS_f_p$gonadf2[raw_PDS_f_p$gonadf == 1 & raw_PDS_f_p$petef == 5] <- 2
raw_PDS_f_p$gonadf2[raw_PDS_f_p$gonadf == 1.5 & raw_PDS_f_p$petef == 5] <- 3
raw_PDS_f_p$gonadf2[raw_PDS_f_p$gonadf == 2 & raw_PDS_f_p$petef == 5] <- 4
raw_PDS_f_p$gonadf2[raw_PDS_f_p$gonadf == 2.5 & raw_PDS_f_p$petef == 5] <- 4
raw_PDS_f_p$gonadf2[raw_PDS_f_p$gonadf == 3 & raw_PDS_f_p$petef == 5] <- 4
raw_PDS_f_p$gonadf2[raw_PDS_f_p$gonadf == 3.5 & raw_PDS_f_p$petef == 5] <- 5
raw_PDS_f_p$gonadf2[raw_PDS_f_p$gonadf == 4 & raw_PDS_f_p$petef == 5] <- 5
raw_PDS_f_p$gonadf2[raw_PDS_f_p$gonadf == 4.5 & raw_PDS_f_p$petef == 5] <- 5
raw_PDS_f_p$gonadf2[raw_PDS_f_p$gonadf == 5 & raw_PDS_f_p$petef == 5] <- 5

raw_PDS_f_p$PDSS <- rowMeans(cbind(raw_PDS_f_p$gonadf2, raw_PDS_f_p$adrenf2), na.rm = TRUE)

# recode variables for males
raw_PDS_m_p <- raw_PDS_m_p %>%
  mutate(petbm = ifelse(petb == 1, 1,
                        ifelse(petb == 2, 2,
                               ifelse(petb == 3, 4,
                                      ifelse(petb == 4, 5, NA))))) %>%
  mutate(petcm = ifelse(petc == 1, 1,
                        ifelse(petc == 2, 2,
                               ifelse(petc == 3, 3,
                                      ifelse(petc == 4, 4, NA)))))

raw_PDS_m_p$adrenm <- rowMeans(cbind(raw_PDS_m_p$petbm, raw_PDS_m_p$petcm), na.rm = TRUE)
raw_PDS_m_p$adrenm2 <- raw_PDS_m_p$adrenm
raw_PDS_m_p$adrenm2[raw_PDS_m_p$adrenm == 1] <- 1
raw_PDS_m_p$adrenm2[raw_PDS_m_p$adrenm == 1.5 & raw_PDS_m_p$petcm == 1] <- 1
raw_PDS_m_p$adrenm2[raw_PDS_m_p$adrenm == 1.5 & raw_PDS_m_p$petcm == 2] <- 2
raw_PDS_m_p$adrenm2[raw_PDS_m_p$adrenm == 2.5 & raw_PDS_m_p$petbm != 4] <- 2
raw_PDS_m_p$adrenm2[raw_PDS_m_p$adrenm == 2.5 & raw_PDS_m_p$petbm == 4] <- 3
raw_PDS_m_p$adrenm2[raw_PDS_m_p$adrenm == 3.5] <- 4
raw_PDS_m_p$adrenm2[raw_PDS_m_p$adrenm == 4.5] <- 5
raw_PDS_m_p$adrenm2[raw_PDS_m_p$adrenm == 5.5] <- 5

raw_PDS_m_p <- raw_PDS_m_p %>%
  mutate(petam = ifelse(peta == 1, 1,
                        ifelse(peta == 2, 3,
                               ifelse(peta == 3, 4,
                                      ifelse(peta == 4, 5, NA))))) %>%
  mutate(petdm = ifelse(petd == 1, 1,
                        ifelse(petd == 2, 2,
                               ifelse(petd == 3, 3,
                                      ifelse(petd == 4, 5, NA)))))

raw_PDS_m_p$gonadm <- rowMeans(cbind(raw_PDS_m_p$petam, raw_PDS_m_p$petdm), na.rm = TRUE)
raw_PDS_m_p$gonadm2 <- raw_PDS_m_p$gonadm
raw_PDS_m_p$gonadm2[raw_PDS_m_p$gonadm == 1 & raw_PDS_m_p$mpete == 1] <- 1
raw_PDS_m_p$gonadm2[raw_PDS_m_p$gonadm == 1 & raw_PDS_m_p$mpete > 1] <- 2
raw_PDS_m_p$gonadm2[raw_PDS_m_p$gonadm == 1.5 & raw_PDS_m_p$mpete == 1] <- 1
raw_PDS_m_p$gonadm2[raw_PDS_m_p$gonadm == 1.5 & raw_PDS_m_p$mpete > 1] <- 2
raw_PDS_m_p$gonadm2[raw_PDS_m_p$gonadm == 2 & raw_PDS_m_p$mpete == 1 & raw_PDS_m_p$petd == 1] <- 1
raw_PDS_m_p$gonadm2[raw_PDS_m_p$gonadm == 2 & raw_PDS_m_p$mpete == 1 & raw_PDS_m_p$petd > 1] <- 2
raw_PDS_m_p$gonadm2[raw_PDS_m_p$gonadm == 2 & raw_PDS_m_p$mpete > 1] <- 3
raw_PDS_m_p$gonadm2[raw_PDS_m_p$gonadm == 2.5 & raw_PDS_m_p$mpete == 1] <- 2
raw_PDS_m_p$gonadm2[raw_PDS_m_p$gonadm == 2.5 & raw_PDS_m_p$mpete > 1] <- 3
raw_PDS_m_p$gonadm2[raw_PDS_m_p$gonadm == 3] <- 3
raw_PDS_m_p$gonadm2[raw_PDS_m_p$gonadm == 3.5 & raw_PDS_m_p$mpete == 2] <- 4
raw_PDS_m_p$gonadm2[raw_PDS_m_p$gonadm == 3.5 & raw_PDS_m_p$mpete == 1] <- 4
raw_PDS_m_p$gonadm2[raw_PDS_m_p$gonadm == 3.5 & raw_PDS_m_p$mpete > 2] <- 5
raw_PDS_m_p$gonadm2[raw_PDS_m_p$gonadm == 4 & raw_PDS_m_p$mpete == 1] <- 4
raw_PDS_m_p$gonadm2[raw_PDS_m_p$gonadm == 4 & raw_PDS_m_p$mpete == 2] <- 4
raw_PDS_m_p$gonadm2[raw_PDS_m_p$gonadm == 4 & raw_PDS_m_p$mpete > 2] <- 5
raw_PDS_m_p$gonadm2[raw_PDS_m_p$gonadm > 4] <- 5

## creating scored dataframe

PDS_f_p <- raw_PDS_f_p
PDS_m_p <- raw_PDS_m_p

PDS_f_p$PDSS <- rowMeans(cbind(PDS_f_p$gonadf2, PDS_f_p$adrenf2), na.rm = TRUE)
PDS_m_p$PDSS <- rowMeans(cbind(PDS_m_p$gonadm2, PDS_m_p$adrenm2), na.rm = TRUE)


track <- read.csv(paste0(data_root, "abcd-general/ab_g_dyn.csv")) %>% 
  rename("id" = "participant_id",
         "wave" = "session_id",
         "rel_family_id" = "ab_g_stc__design_id__fam") %>% 
  select(id, wave, site, rel_family_id)

PDS_f_p <- left_join(PDS_f_p, track, by = c("id", "wave"))
PDS_m_p <- left_join(PDS_m_p, track, by = c("id", "wave"))

PDS_f_p <- PDS_f_p %>%
  group_by(id) %>%
  mutate(rel_family_id = ifelse(is.na(rel_family_id), first(rel_family_id), rel_family_id)) %>%   
  ungroup()

PDS_m_p <- PDS_m_p %>%
  group_by(id) %>%
  mutate(rel_family_id = ifelse(is.na(rel_family_id), first(rel_family_id), rel_family_id)) %>%   
  ungroup()

write.csv(PDS_f_p, file = paste0(data_root, "physical-health/puberty/parent_tannerstages_f.csv"))
write.csv(PDS_m_p, file = paste0(data_root, "physical-health/puberty/parent_tannerstages_m.csv"))

# now, apply the function to each dataframe
PDS_f_p_filtered <- select_participant_by_obs(PDS_f_p)
PDS_m_p_filtered <- select_participant_by_obs(PDS_m_p)

write.csv(PDS_f_p_filtered, file = paste0(data_root, "physical-health/puberty/filtered_parent_tannerstages_f.csv"))
write.csv(PDS_m_p_filtered, file = paste0(data_root, "physical-health/puberty/filtered_parent_tannerstages_m.csv"))
