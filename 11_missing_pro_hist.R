#load 0505
colnames(protome_all)
missing = data.frame(missing = rowSums(is.na(protome_all[16:69])))
missing$completed = 54-missing$missing

g = 
  ggplot(missing,aes(x=completed))+  geom_histogram(color="black",fill="black",binwidth=1)+ theme_classic() +
  scale_color_npg()+
  #theme(text = element_text( family = 'Arial'))+
  xlab("Number of Measurements in Set1 & Set2")+
  ylab("Number of Proteins")+
  ggtitle("Measurement Completeness in \nProteomics Data")
g
ggsave("11_9-missing_pro.pdf",height=4,width=4,plot=g)

missing$percentage = missing$completed/54 * 100

# Load the necessary library
library(ggplot2)

# Example dataframe
df = missing
# Bin the percentages into groups
df$group <- cut(df$percentage, 
                breaks = c(-Inf, 20, 40, 60, 80, 100,Inf), 
                labels = c("0-20%", "20-40%", "40-60%", "60-80%", "80-100%","100%"),
                right = FALSE)

# Count the number of proteins in each group
group_counts <- as.data.frame(table(df$group))

# Calculate the percentage for each group
group_counts$percentage <- group_counts$Freq / sum(group_counts$Freq) * 100

# Create a label column for the pie chart
group_counts$label <- paste(group_counts$Var1, "\n", group_counts$Freq, " proteins")

# Create the pie chart
pie_chart <- ggplot(group_counts, aes(x = "", y = percentage, fill = Var1)) +
  scale_fill_npg()+
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0) +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3) +
  theme_void() +
  labs(fill = "Percentage Range", 
       title = "Distribution of Proteins by Percentage Presence in Samples") +
  
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 8))

pie_chart
ggsave("missing_pro_new.pdf",height=6,width=6,plot=pie_chart)

