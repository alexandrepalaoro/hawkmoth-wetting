require(ggplot2)
require(ggpubr)

labv = read.csv(file = "data/hawkmoth-labviewimagej.csv", h = T, sep = ';')

labv$species = gsub("_", " ", labv$species)

png(filename = "figures/labvimagej.png", res = 600, units = 'mm',
    h = 160, w = 220)
ggplot(labv, aes(x=as.factor(species), y=avg, fill=type)) +
  geom_bar(position=position_dodge(), stat="identity", colour='black') +
  geom_errorbar(aes(ymin=avg-sd, ymax=avg+sd), width=.2,position=position_dodge(.9)) +
  theme_bw() + xlab("Species") + ylab("Average contact angle") + 
  scale_fill_discrete(name = "Type of calculation", labels = c("ImageJ","LabVIEW")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


sides = read.csv(file = "data/hawkmoth-rightandleft.csv", h = T, sep = ';')

sides$species = gsub("_", " ", sides$species)

png(filename = "figures/rightandleft.png", res = 600, units = 'mm',
    h = 160, w = 260)
ggplot(sides, aes(x = as.factor(species), y = avg, fill = side)) +
  geom_bar(position=position_dodge(), stat="identity", colour='black') +
  geom_errorbar(aes(ymin=avg-sd, ymax=avg+sd), width=.2,position=position_dodge(.9)) +
  theme_bw() + xlab("Species") + ylab("Average contact angle") + 
  scale_fill_discrete(name = "Galea measured", labels = c("Left","Right")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
