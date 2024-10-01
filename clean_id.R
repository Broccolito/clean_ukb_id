library(dplyr)
library(data.table)
library(purrr)

sequencing = fread("ukb_wgs_ids.fam") %>%
  filter(V6 == "-9")
names(sequencing)[c(1,2)] = c("FID", "IID")
sequencing = sequencing %>%
  select(FID, IID) %>%
  mutate(fid_iid = paste(FID, IID, sep = "_"))

afr = fread("UKBiobank_genoQC_AFR_unrelated_list_kinship_0.0442.txt") %>%
  mutate(ancestry = "afr")
eur = fread("UKBiobank_genoQC_EUR_unrelated_list_kinship_0.0442.txt") %>%
  mutate(ancestry = "eur")
eas = fread("UKBiobank_genoQC_EAS_unrelated_list_kinship_0.0442.txt") %>%
  mutate(ancestry = "eas")
sas = fread("UKBiobank_genoQC_SAS_unrelated_list_kinship_0.0442.txt") %>%
  mutate(ancestry = "sas")

sas_sri = fread("UKBiobank_genoQC_SAS_Sri_Lankan_unrelated_list_kinship_0.0442.txt") %>%
  mutate(ancestry = "sas") %>%
  mutate(subancestry = "sas_sri_lankan") %>%
  mutate(fid_iid = paste(FID, IID, sep = "_"))
sas_ind = fread("UKBiobank_genoQC_SAS_Indian_unrelated_list_kinship_0.0442.txt") %>%
  mutate(ancestry = "sas") %>%
  mutate(subancestry = "sas_indian") %>%
  mutate(fid_iid = paste(FID, IID, sep = "_"))
sas_ban = fread("UKBiobank_genoQC_SAS_Bangladeshi_unrelated_list_kinship_0.0442.txt") %>%
  mutate(ancestry = "sas") %>%
  mutate(subancestry = "sas_bangladeshi") %>%
  mutate(fid_iid = paste(FID, IID, sep = "_"))
sas_pak = fread("UKBiobank_genoQC_SAS_Pakistani_unrelated_list_kinship_0.0442.txt") %>%
  mutate(ancestry = "sas") %>%
  mutate(subancestry = "sas_Pakistani") %>%
  mutate(fid_iid = paste(FID, IID, sep = "_"))
afr_car = fread("UKBiobank_genoQC_AFR_Caribbean_unrelated_list_kinship_0.0442.txt") %>%
  mutate(ancestry = "afr") %>%
  mutate(subancestry = "afr_caribbean") %>%
  mutate(fid_iid = paste(FID, IID, sep = "_"))
afr_noc = fread("UKBiobank_genoQC_AFR_Non_Caribbean_unrelated_list_kinship_0.0442.txt") %>%
  mutate(ancestry = "afr") %>%
  mutate(subancestry = "afr_non_caribbean") %>%
  mutate(fid_iid = paste(FID, IID, sep = "_"))

genotyping = rbind.data.frame(
  afr, eur, eas, sas
) %>%
  mutate(fid_iid = paste(FID, IID, sep = "_"))

unrelated_common_set = inner_join(sequencing, genotyping, by = "fid_iid", 
                                  suffix = c("", "_geno")) %>%
  select(FID, IID, ancestry, fid_iid) %>%
  mutate(subancestry = ancestry) %>%
  mutate(subancestry = ifelse(
    fid_iid %in% sas_sri$fid_iid, "sas_sri_lankan", ifelse(
      fid_iid %in% sas_ind$fid_iid, "sas_indian", ifelse(
        fid_iid %in% sas_ban$fid_iid, "sas_bangladeshi", ifelse(
          fid_iid %in% sas_pak$fid_iid, "sas_Pakistani", ifelse(
            fid_iid %in% afr_car$fid_iid, "afr_caribbean", ifelse(
              fid_iid %in% afr_noc$fid_iid, "afr_non_caribbean", subancestry
            )
          )
        )
      )
    )
  )) %>%
  filter(subancestry != "afr") %>%
  select(FID, IID, ancestry, subancestry)


