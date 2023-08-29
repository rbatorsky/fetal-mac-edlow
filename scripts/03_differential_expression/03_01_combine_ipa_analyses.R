
library(tidyverse)
library(openxlsx)
files <- Sys.glob("~/analysis/ipa/mf_together/Canonical Pathways/export_each_cluster_cp_abslfc_0.25/*xlsx")
df =NULL
files
for (f in files) {
  print(f)
  #new=read_excel(f, skip = 2)

  new$cluster = gsub("_MF.xlsx","",basename(f))
  if(nrow(new)>0)
    if(is.null(df))
      df <- new
    else
      df = rbind(df, new)
}

write.xlsx(df, "~/analysis/ipa/mf_together/Canonical Pathways/export_cp_mftogether_abslfc_0.25.xlsx")


files <- Sys.glob("~/analysis/ipa/mf_together/Canonical Pathways/export_each_cluser_cp_abslfc_0.5/*txt")
df =NULL
files
for (f in files) {
  print(f)
  new=read.table(f, skip = 2, header=T, sep="\t")

  new$cluster = gsub("_MF_0.5.txt","",basename(f))
  if(nrow(new)>0)
    if(is.null(df))
      df <- new
  else
    df = rbind(df, new)
}

write.xlsx(df, "~/analysis/ipa/mf_together/Canonical Pathways/export_cp_mftogether_abslfc_0.5.xlsx")


files <- Sys.glob("~/analysis/ipa/mf_together/Canonical Pathways/export_cp_mftogether_abslfc_0.3/*.xlsx")
df =NULL
files
for (f in files) {
  print(f)
  new=read_excel(f, skip = 2)

  new$cluster = gsub("_MF_0.3.xlsx","",basename(f))

  if(nrow(new)>0)
    print(new %>% dplyr::filter(grepl("CTL",`Ingenuity Canonical Pathways`)))
    if(is.null(df))
      df <- new
  else
    df = rbind(df, new)
}
head(df)

df %>% dplyr::filter(grepl("CTL",`Ingenuity Canonical Pathways`))
write.xlsx(df, "~/analysis/ipa/mf_together/Canonical Pathways/export_cp_mftogether_abslfc_0.3.xlsx")

files <- Sys.glob("analysis/ipa/mf_together/df/*.xlsx")
df =NULL
files
for (f in files) {
  print(f)
  
  new=read_excel(f)
  name=gsub(".xlsx","",basename(f))
  print(name)
  new$cluster = name
  print(head(new))
  
  if(is.null(df))
    df <- new
  else
    df = rbind(df, new)
}

write.xlsx(df, "analysis/ipa/mf_together/diseases_functions_lfc0.25.xlsx")


# need to add this one MgYSI_F_cp.txt
files <- Sys.glob("analysis/ipa/mf_separate_combined_cluster/cp/*xlsx")
df =NULL
files
for (f in files) {
  print(f)
  new=read_excel(f, skip=1)
  name=gsub("_cp.xlsx","",basename(f))
  print(name)
  new$cluster = name
  print(head(new))
  
  #if(nrow(new)>0)
  #  print(new %>% dplyr::filter(grepl("CTL",`Ingenuity Canonical Pathways`)))
  if(is.null(df))
    df <- new
  else
    df = rbind(df, new)
}
head(df)
write.xlsx(df, "analysis/ipa/mf_separate_combined_cluster/canonical_pathways_separate_lfc0.25.xlsx")


files <- Sys.glob("analysis/ipa/mf_separate_combined_cluster/df/*.xlsx")
df =NULL
files
for (f in files) {
  #print(f)
  new=read_excel(f, skip=1)
  name=gsub("_df.xlsx","",basename(f))
  print(name)
  new$cluster = name
  print(head(new))
  
  #if(nrow(new)>0)
  #  print(new %>% dplyr::filter(grepl("CTL",`Ingenuity Canonical Pathways`)))
  if(is.null(df))
    df <- new
  else
    df = rbind(df, new)
}
head(df)
write.xlsx(df, "analysis/ipa/mf_separate_combined_cluster/diseases_functions_separate_lfc0.25.xlsx")
