library(ggplot2)
library(ggrepel)
library(stringr)
library(dplyr)
library(tidyr)
library(gridExtra)
library(grid)
library(corrplot)
library(MultinomialCI)

dt_raw <- read.csv("scopus_CodedFINAL_Rinput.csv", header = T)

#clean-up raw data
#step1: extract studies that are included for the assessment, 
dt <- dt_raw[which(startsWith(dt_raw$included, "Yes")),]


#recoding function to be used for all indicators (replace original name by a new shortname [newcode_df])
recoding <- function(raw_df, raw_df_col_N, newcode_df)  {
  for (i in 1:nrow(raw_df)){
    if (is.na(raw_df[,raw_df_col_N][i])){
      raw_df[,raw_df_col_N][i] <- "N/A"
    } else {
      for (j in 1:nrow(newcode_df)){
        if(raw_df[,raw_df_col_N][i] == newcode_df[j,1]){
          raw_df[,raw_df_col_N][i] <- newcode_df[j,2]
        } 
      }
    }
  }
  return(raw_df)
}




#step2: assign new shortname to each indicator
#col = 7 for "whether or not primary LCI data are available?", with four initial levels:
pri_avl_raw <- levels(factor(dt$primary_available))
pri_avl_new <- c("wo_statement", "others", "wh_statement_available",  "wh_statement_not_available")
pri_avl_final <- cbind(pri_avl_raw,pri_avl_new)

#col = 8 for "How does the study indicate that primary LCI data are available", with three initial levels, excluding NA
pri_how_raw <- levels(factor(dt$primary_how ))
pri_how_new <- c("3rd_party", "in_paper", "journal_SI")
pri_how_final <- cbind(pri_how_raw,pri_how_new)

#col = 9 for "whether or not secondary LCI data are available?", with six initial levels:
bk_avl_raw <- levels(factor(dt$secondary_available))
bk_avl_new <- c("wo_statement", "others", "source_only_no_data",  "available&reproducible", "not_available", "available_not_reproducible")
bk_avl_final <- cbind(bk_avl_raw,bk_avl_new)

#col = 10 for "Secondary LCI source", various sources, split all source:
BKDT_split <- unlist(strsplit(dt$secondary_source, ";")) 
for (i in 1:length(BKDT_split)) {
  BKDT_split[i] <- gsub(paste("and", collapse = "|"), "", BKDT_split[i])   #delete all "and"
  if(grepl(" ",BKDT_split[i],  fixed = TRUE)) {
    BKDT_split[i] <- trimws(BKDT_split[i], which = c("both"))              #delete all space
  }
  if(grepl("ecoinvent",BKDT_split[i],  fixed = TRUE) | grepl("Ecoinvent",BKDT_split[i],  fixed = TRUE)) {
    BKDT_split[i] = "Ecoinvent"                        #to Ecoinvent (regardless of version)
  }
  if(grepl("Agri-footprint",BKDT_split[i],  fixed = TRUE )){
    BKDT_split[i] = "Agri-footprint"    
  }
  if(grepl("Gabi",BKDT_split[i], fixed = TRUE)  | grepl("GaBi",BKDT_split[i],  fixed = TRUE) | grepl("Gabi",BKDT_split[i],  fixed = TRUE)  ){
    BKDT_split[i] = "Gabi"
  }
}
#delete all "various sources" as only focus on LCI database: 
BKDT_split <- BKDT_split[-which(BKDT_split == 'various sources')]
as.data.frame(table(BKDT_split))

length(BKDT_split)


#col = 11 for "How does the study indicate that secondary LCI data are available?", with 5 initial levels (excluding NA):
bk_how_raw <- levels(factor(dt$secondary_how))
bk_how_new <- c("author_website", "3rd_party", "in_paper", "within_analysis_script", "journal_SI")
bk_how_final <- cbind(bk_how_raw,bk_how_new)


#col = 14 for "Analysis tool used", e.g., SimaPro, Gabi, openLCA, keep the raw text
LCA_tool_raw <- levels(factor(dt$analysis_tool))
LCA_tool_final <- LCA_tool_raw
#too many tool to plot, choose the top 10 tools

#col = 15 for "whether or not analysis scripts are available?" 
LCA_ANAL_raw <- levels(factor(dt$analysis_available))
LCA_ANAL_new <- c("wo_statement", "others", "not_available", "available")
LCA_ANAL_final <- cbind(LCA_ANAL_raw,LCA_ANAL_new)

#col = 16 for "How does the statement indicate the analysis scripts are available"
ANAL_how_raw <- levels(factor(dt$analysis_how))
ANAL_how_new <- c("3rd_party","journal_SI")
ANAL_how_final <- cbind(ANAL_how_raw,ANAL_how_new)



#step3: plotting 
#bar_plot function to be applied for each var. 
bplot <- function (df_plot,legent_title) {
  if (nrow(df_plot) > 9) {               #only extract the first 10 rows for plotting if too many attributes (e.g., databases; LCA tool)
    df_plot <- df_plot[order(df_plot$Freq, decreasing = T), ][1:10,]    
  } 
  p <- ggplot(data=df_plot, aes(na.rm = TRUE, x="", y=Freq, fill=df_plot[,1])) +  
    geom_bar(stat="identity") +
    coord_flip() + 
    labs(subtitle = legent_title, caption = paste("N =", sum(df_plot$Freq)), fill="") +  
    geom_text(aes(label = Freq), size = 2, angle = 0, hjust = 0.5, position = position_stack(vjust = 0.5) ) + 
    scale_fill_brewer(palette="Paired") +
    theme(legend.position="top", legend.text = element_text(size=7) ) +   
    xlab("") + ylab("") 
  
  return(p)
}





lab_text <- c("Primary data availability", "Primary data presentation", "Secondary data availability", "Secondary data presentation", 
              "LCA analysis script availability", "LCA analysis script presentation", "Secondary data source", "LCA analysis tool")  

#too many sources for LCI databases and LCA tools, select most 10 frequent, then combine others  
LCIDB_plot <- as.data.frame(table(BKDT_split))[!as.data.frame(table(BKDT_split))$BKDT_split == "N/A", ]
LCIDB_top9 = LCIDB_plot[order(LCIDB_plot$Freq, decreasing = T),][1:9,]
LCIDB_other <- cbind("other_databases", 
                     sum(as.data.frame(table(BKDT_split))[!as.data.frame(table(BKDT_split))$BKDT_split == "N/A", ]$Freq) - sum(LCIDB_top9$Freq))
colnames(LCIDB_other)<-colnames(LCIDB_plot)
LCIDB_plot <- rbind(LCIDB_top9,LCIDB_other)
LCIDB_plot$Freq <- as.numeric(LCIDB_plot$Freq)

tool <- as.data.frame(table(dt$analysis_tool))[!as.data.frame(table(dt$analysis_tool))$Var1 == "N/A", ]
tool_top9 = tool[order(tool$Freq, decreasing = T),][1:9,]
tool_other <- cbind("other_tools", 
                    sum(as.data.frame(table(dt$analysis_tool))[!as.data.frame(table(dt$analysis_tool))$Var1 == "N/A", ]$Freq) - sum(tool_top9$Freq))
colnames(tool_other)<-colnames(tool)
tool_plot <- rbind(tool_top9,tool_other)
tool_plot$Freq <- as.numeric(tool_plot$Freq)
  
### delete all "N/A" not applicable in ploting
table_list <- list( 
  as.data.frame(table(recoding(dt,7,pri_avl_final)[,7])),
  as.data.frame(table(recoding(dt,8,pri_how_final)[,8]))[!as.data.frame(table(recoding(dt,8,pri_how_final)[,8]))$Var1 == "N/A", ],
  as.data.frame(table(recoding(dt,9,bk_avl_final)[,9])),
  as.data.frame(table(recoding(dt,11,bk_how_final)[,11]))[!as.data.frame(table(recoding(dt,11,bk_how_final)[,11]))$Var1 == "N/A", ],
  as.data.frame(table(recoding(dt,15,LCA_ANAL_final)[,15])),
  as.data.frame(table(recoding(dt,16,ANAL_how_final)[,16]))[!as.data.frame(table(recoding(dt,16,ANAL_how_final)[,16]))$Var1 == "N/A", ],
  LCIDB_plot,
  #as.data.frame(table(BKDT_split))[!as.data.frame(table(BKDT_split))$BKDT_split == "N/A", ],
  tool_plot
  #as.data.frame(table(dt$analysis_tool))[!as.data.frame(table(dt$analysis_tool))$Var1 == "N/A", ]
)


length(lab_text) == length(table_list)

for (n in 1:8){
  assign(paste("x",n, sep=""), bplot(table_list[[n]],lab_text[n]) )
}
grid.arrange(x1, x2, x3, x4, x5, x6, x7, x8, ncol=2)
