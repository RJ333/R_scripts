#fancy scientific hack (funktionioert nicht ganz)

> fancy_scientific <- function(l) {
      # turn in to character string in scientific notation
     l <- format(l, scientific = TRUE)
	 l <- gsub("0e\\+00","0",l)
     # quote the part before the exponent to keep all the digits
     l <- gsub("^(.*)e", "'\\1'e", l)
     # turn the 'e+' into plotmath format
     l <- gsub("e", "%*%10^", l)
	 # remove + after exponent, if exists. E.g.: (3x10^+2 -> 3x10^2) 
	 l <- gsub("e\\+","e",l)
	 # convert 1x10^ or 1.000x10^ -> 10^ 
	 l <- gsub("\\'1[\\.0]*\\'\\%\\*\\%", "", l)
	 # return this as an expression
     parse(text=l)
}

#https://rdrr.io/r/grDevices/plotmath.html

#http://www.cookbook-r.com/Graphs/Axes_(ggplot2)/
## Self-defined formatting function for times.
timeHMS_formatter <- function(x) {
    h <- floor(x/60)
    m <- floor(x %% 60)
    s <- round(60*(x %% 1))                   # Round to nearest second
    lab <- sprintf('%02d:%02d:%02d', h, m, s) # Format the strings as HH:MM:SS
    lab <- gsub('^00:', '', lab)              # Remove leading 00: if present
    lab <- gsub('^0', '', lab)                # Remove leading 0 if present
}

bp + scale_y_continuous(label=timeHMS_formatter)


#superscript 
labs(x=expression(Production~rate~" "~mu~moles~NO[3]^{-1}-N~Kg^{-1}),
     y=expression(glyphosate concentration~mg~L^{-1}))

	 
	 
#second axis + scientific (funktioniert!)
scale_y_continuous(label= function(x) {ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))}, sec.axis=sec_axis(~./2200000,name=expression(paste("Glyphosate concentration    ",bgroup("[",mg~L^{-1},"]"))))) +