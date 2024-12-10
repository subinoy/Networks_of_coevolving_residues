## 1. FEN1
fen1_intra <- read_csv("data/intra_FEN1.csv")
fen1_intra
## Applying function to process the data frame
fen1_intra <- df_intra_process(fen1_intra)


## Generating frequency table
fen1_intra_cov_freq <- prot_intra_cov_freq(fen1_intra,"FEN1")
## Plotting bar plot
#draw_bar_plot(fen1_intra_cov_freq )

intra_draw_bar_plot(fen1_intra_cov_freq, "FEN1")

# ggplot(fen1_intra_cov_freq, aes(AA1, y=freq, fill=protein)) +
#     geom_bar(stat="identity")

# **************************************************************
## 2. PCNA
pcna_intra <- read_csv("data/intra_PCNA.csv")
pcna_intra
## Applying function to process the data frame
pcna_intra <- df_intra_process(pcna_intra)
str(pcna_intra)
head(pcna_intra)
## Generating frequency table
pcna_intra_cov_freq <- prot_intra_cov_freq(pcna_intra,"PCNA")
## Plotting bar plot
#draw_bar_plot(fen1_intra_cov_freq )

intra_draw_bar_plot(pcna_intra_cov_freq, "PCNA")

# **************************************************************
## 3. POLD2
pold2_intra <- read_csv("data/intra_POLD2.csv")
pold2_intra
## Applying function to process the data frame
pold2_intra <- df_intra_process(pold2_intra)


## Generating frequency table
pold2_intra_cov_freq <- prot_intra_cov_freq(pold2_intra,"POLD2")
## Plotting bar plot
#draw_bar_plot(fen1_intra_cov_freq )

intra_draw_bar_plot(pold2_intra_cov_freq, "POLD2")


# **************************************************************
## 4.  RFC3
rfc3_intra <- read_csv("data/intra_RFC3.csv")
rfc3_intra
## Applying function to process the data frame
rfc3_intra <- df_intra_process(rfc3_intra)


## Generating frequency table
rfc3_intra_cov_freq <- prot_intra_cov_freq(rfc3_intra,"RFC3")
## Plotting bar plot
#draw_bar_plot(fen1_intra_cov_freq )

intra_draw_bar_plot(rfc3_intra_cov_freq, "RFC3")

# *************************************************************
## 5.  RFC1
rfc1_intra <- read_csv("data/intra_RFC1.csv")
rfc1_intra
## Applying function to process the data frame
rfc1_intra <- df_intra_process(rfc1_intra)


## Generating frequency table
rfc1_intra_cov_freq <- prot_intra_cov_freq(rfc1_intra,"RFC1")
## Plotting bar plot
#draw_bar_plot(fen1_intra_cov_freq )

intra_draw_bar_plot(rfc1_intra_cov_freq, "RFC1")

# _____________________________________________________________
# *************************************************************
## 6.  RFC4
rfc4_intra <- read_csv("data/intra_RFC4.csv")
rfc4_intra
## Applying function to process the data frame
rfc4_intra <- df_intra_process(rfc4_intra)


## Generating frequency table
rfc4_intra_cov_freq <- prot_intra_cov_freq(rfc4_intra,"RFC4")
## Plotting bar plot
#draw_bar_plot(fen1_intra_cov_freq )

intra_draw_bar_plot(rfc4_intra_cov_freq, "RFC4")

# _____________________________________________________________

# **************************************************************
## 7. POLD1
pold1_intra <- read_csv("data/intra_POLD1.csv")
pold1_intra
## Applying function to process the data frame
pold1_intra <- df_intra_process(pold1_intra)


## Generating frequency table
pold1_intra_cov_freq <- prot_intra_cov_freq(pold1_intra,"POLD1")
## Plotting bar plot
#draw_bar_plot(fen1_intra_cov_freq )

intra_draw_bar_plot(pold1_intra_cov_freq, "POLD1")


# **************************************************************

uniq_fen1 <- read.csv("data/total_uniq_FEN1", header=T, row.names = 1)
head(uniq_fen1)
dim(uniq_fen1)
write_csv(uniq_fen1, path = "data/uniq_FEN1_intra_table.csv")

save.image(file="data/Sept_25th_upto_RFC4.RData")
