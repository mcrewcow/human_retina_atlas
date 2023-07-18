midpoint <-sum(atlas$Percentage) -  cumsum(atlas$Percentage) + atlas$Percentage/2

ggplot(atlas, aes(x = "", y=Percentage, fill = Type)) +
  geom_bar(width = 1, stat = "identity") +
  scale_y_continuous(breaks=midpoint, labels=atlas$Type) +
  scale_fill_brewer(palette = "Blues", direction = -1) +
  labs(fill="Bachelor/ Vordiplom",x=NULL,y=NULL,title="",caption="")  +
  geom_text_repel(aes(label = round(Percentage,2)),
                  position = position_stack(vjust = 0.5),
                  show.legend = FALSE) +
  coord_polar(theta = "y", start=0) +
  theme(axis.ticks = element_blank(), panel.grid  = element_blank(), axis.line = element_blank(), strip.background = element_blank(), panel.background = element_blank())
