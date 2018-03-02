#library(ggplot2)
#library(scales)
#Gallaecimonas
Gallaecimonas_abs_water_plot<-subset(abs_abu_molten_tax, habitat == "water" & grepl('Gallaecimonas',genus) & days > 40 & treatment =="glyph")
species_title_Gallaecimonas<-expression(paste(,italic("Gallaecimonas")," sp."))
Gallaecimonas_abs_water_plot$treatment2<-factor(Gallaecimonas_abs_water_plot$treatment,labels="Treatment")



ggplot(data=Gallaecimonas_abs_water_plot, aes(x=days-69,y=value,colour=nucleic_acid,group=nucleic_acid,lty=nucleic_acid))+ 
	geom_vline(data=subset(Gallaecimonas_abs_water_plot,treatment=="glyph"),aes(xintercept=0),linetype="dashed", size=1.2)+
	geom_point(aes(),colour="black")+
	stat_summary(aes(colour=nucleic_acid),colour="black",fun.y="mean", geom="line", size=1.5)+
	scale_linetype_manual(values=c("dna"=1,"cdna"=4),
						name="Nucleic acid  ",
						breaks=c("cdna","dna"),
						labels=c("16S rRNA","16S rDNA"))+
	scale_x_continuous(breaks = scales::pretty_breaks(n = 20))+
	theme_bw()+
	scale_y_continuous(label= function(x) {ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))})+
	theme(axis.title.y = element_text(angle=90,vjust=0.5))+
	theme(axis.text=element_text(size=12))+
	theme(legend.text=element_text(size=11))+
	theme(panel.grid.major=element_line(colour = NA, size = 0.2))+
	theme(panel.grid.minor=element_line(colour = NA, size = 0.5))+
	theme(legend.position="bottom")+
	theme(legend.background = element_rect(fill="grey90",linetype="solid"))+
	labs(x="Days",
			 y=expression(atop("Absolute abundance in cell equivalents",bgroup("[",relative~abundance~x~cells~mL^{-1},"]"))))

			 
sample_data<-Gallaecimonas_abs_water_plot[,c(1,2,4,5,6,14,20)]
require(ggplot2)
require(scales)
ggplot(data=sample_data, aes(x=days-69,y=value,colour=nucleic_acid,group=nucleic_acid,lty=nucleic_acid))+ 
	geom_vline(aes(xintercept=0),linetype="dashed", size=1.2)+
	geom_point(aes(),colour="black")+
	stat_summary(aes(colour=nucleic_acid),colour="black",fun.y="mean", geom="line", size=1.5)+
	scale_linetype_manual(values=c("dna"=1,"cdna"=4),
						name="Nucleic acid  ",
						breaks=c("cdna","dna"),
						labels=c("16S rRNA","16S rDNA"))+
	scale_x_continuous(breaks = scales::pretty_breaks(n = 20))+
	theme_bw()+
	scale_y_continuous(label= function(x) {ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))})+
	theme(axis.title.y = element_text(angle=90,vjust=0.5))+
	theme(axis.text=element_text(size=12))+
	theme(legend.text=element_text(size=11))+
	theme(panel.grid.major=element_line(colour = NA, size = 0.2))+
	theme(panel.grid.minor=element_line(colour = NA, size = 0.5))+
	theme(legend.position="bottom")+
	theme(legend.background = element_rect(fill="grey90",linetype="solid"))+
	labs(x="Days",
			 y=expression(atop("Absolute abundance in cell equivalents",bgroup("[",relative~abundance~x~cells~mL^{-1},"]"))))
			 
sample_data2<-structure(list(time = c(10L, 10L, 10L, 10L, 10L, 10L, 11L, 11L, 
11L, 11L, 11L, 11L, 12L, 12L, 12L, 12L, 12L, 12L, 13L, 13L, 13L, 
13L, 13L, 13L, 14L, 14L, 14L, 14L, 14L, 14L, 15L, 15L, 15L, 15L, 
15L, 15L, 16L, 16L, 16L, 16L, 16L, 16L, 17L, 17L, 17L, 17L, 18L, 
18L, 18L, 18L, 18L, 18L, 19L, 19L, 19L, 19L, 19L, 19L, 4L, 4L, 
4L, 4L, 4L, 4L, 5L, 5L, 5L, 5L, 5L, 5L, 6L, 6L, 6L, 6L, 6L, 6L, 
7L, 7L, 7L, 7L, 7L, 7L, 8L, 8L, 8L, 8L, 8L, 8L, 9L, 9L, 9L, 9L, 
9L, 9L), days = c(83L, 83L, 83L, 83L, 83L, 83L, 86L, 86L, 86L, 
86L, 86L, 86L, 91L, 91L, 91L, 91L, 91L, 91L, 98L, 98L, 98L, 98L, 
98L, 98L, 105L, 105L, 105L, 105L, 105L, 105L, 112L, 112L, 112L, 
112L, 112L, 112L, 119L, 119L, 119L, 119L, 119L, 119L, 126L, 126L, 
126L, 126L, 133L, 133L, 133L, 133L, 133L, 133L, 140L, 140L, 140L, 
140L, 140L, 140L, 44L, 44L, 44L, 44L, 44L, 44L, 62L, 62L, 62L, 
62L, 62L, 62L, 69L, 69L, 69L, 69L, 69L, 69L, 72L, 72L, 72L, 72L, 
72L, 72L, 76L, 76L, 76L, 76L, 76L, 76L, 79L, 79L, 79L, 79L, 79L, 
79L), parallel = c(3L, 1L, 2L, 2L, 3L, 1L, 2L, 3L, 3L, 2L, 1L, 
1L, 2L, 1L, 3L, 3L, 1L, 2L, 2L, 3L, 3L, 1L, 1L, 2L, 2L, 3L, 1L, 
1L, 3L, 2L, 1L, 1L, 2L, 3L, 3L, 2L, 2L, 3L, 3L, 1L, 1L, 2L, 3L, 
1L, 1L, 3L, 2L, 3L, 1L, 1L, 2L, 3L, 1L, 2L, 3L, 3L, 1L, 2L, 2L, 
3L, 3L, 1L, 1L, 2L, 2L, 3L, 1L, 1L, 3L, 2L, 1L, 2L, 3L, 3L, 1L, 
2L, 2L, 3L, 3L, 1L, 1L, 2L, 2L, 1L, 1L, 2L, 3L, 3L, 1L, 2L, 3L, 
3L, 1L, 2L), nucleic_acid = structure(c(1L, 1L, 1L, 2L, 2L, 2L, 
2L, 2L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 1L, 1L, 1L, 2L, 2L, 1L, 2L, 
1L, 1L, 2L, 2L, 2L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 2L, 2L, 1L, 1L, 
2L, 1L, 2L, 2L, 1L, 1L, 2L, 2L, 1L, 1L, 1L, 2L, 2L, 2L, 1L, 1L, 
1L, 2L, 2L, 2L, 2L, 2L, 1L, 2L, 1L, 1L, 2L, 2L, 2L, 1L, 1L, 1L, 
2L, 2L, 2L, 1L, 1L, 1L, 2L, 2L, 1L, 2L, 1L, 1L, 2L, 2L, 1L, 1L, 
1L, 2L, 1L, 1L, 1L, 2L, 2L, 2L), .Label = c("cdna", "dna"), class = "factor"), 
    habitat = structure(c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L), .Label = "water", class = "factor"), 
    value = c(5316639.62, 6402573.912, 6294710.95, 2369809.996, 
    2679661.691, 2105693.166, 2108794.224, 2487177.041, 6021765.438, 
    5524939.499, 6016021.786, 2628427.206, 3164229.113, 896068.7656, 
    2966515.364, 4436008.425, 1860580.149, 3911309.508, 888489.0268, 
    1004334.365, 1141636.992, 961140.0729, 1072009.18, 1134997.852, 
    668013.4333, 459645.1058, 645944.1129, 702293.6865, 590620.3693, 
    642136.7523, 932531.1588, 1224299.065, 1502344.5, 1545034.46, 
    1122002.798, 1411050.57, 1465061.711, 1378876.488, 810348.2823, 
    1361496.248, 1056558.288, 897876.4169, 931519.9524, 1165768.09, 
    957873.9045, 746011.7558, 624116.5603, 522209.2283, 551120.1371, 
    440096.4446, 565108.4447, 373304.8604, 266595.7171, 333767.4042, 
    185612.6681, 144899.8736, 173739.3969, 211490.827, 223815.0867, 
    296455.4243, 1278759.217, 247292.4355, 1171554.199, 1146278.577, 
    227443.8462, 233542.6719, 253224.2629, 875040.4892, 1151921.616, 
    1285744.479, 355381.9156, 110724.7928, 252238.9632, 912865.3372, 
    608269.6498, 500307.5301, 774955.9598, 1374106.94, 3121909.308, 
    1071086.757, 3033665.589, 2984567.998, 1396313.444, 1356465.773, 
    4480581.956, 4273141.231, 4957691.655, 1910056.657, 5520085.32, 
    5094686.657, 5990052.759, 2272441.566, 1513268.608, 1821716.75
    ), treatment2 = structure(c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L), .Label = "Treatment", class = "factor")), .Names = c("time", 
"days", "parallel", "nucleic_acid", "habitat", "value", "treatment2"
), class = "data.frame", row.names = c(51243L, 51244L, 51245L, 
51246L, 51247L, 51248L, 51255L, 51256L, 51257L, 51258L, 51259L, 
51260L, 51267L, 51268L, 51269L, 51270L, 51271L, 51272L, 51279L, 
51280L, 51281L, 51282L, 51283L, 51284L, 51291L, 51292L, 51293L, 
51294L, 51295L, 51296L, 51303L, 51304L, 51305L, 51306L, 51307L, 
51308L, 51315L, 51316L, 51317L, 51318L, 51319L, 51320L, 51326L, 
51327L, 51328L, 51329L, 51336L, 51337L, 51338L, 51339L, 51340L, 
51341L, 51348L, 51349L, 51350L, 51351L, 51352L, 51353L, 51360L, 
51361L, 51362L, 51363L, 51364L, 51365L, 51372L, 51373L, 51374L, 
51375L, 51376L, 51377L, 51384L, 51385L, 51386L, 51387L, 51388L, 
51389L, 51396L, 51397L, 51398L, 51399L, 51400L, 51401L, 51408L, 
51409L, 51410L, 51411L, 51412L, 51413L, 51420L, 51421L, 51422L, 
51423L, 51424L, 51425L))