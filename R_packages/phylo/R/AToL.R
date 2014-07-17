## DEFINITIONS ##
## from David Tank: based on Soltis et al. 2011 (AJB) and Cantino et al. 2007 (Taxon)
spermatophyta_APGIII_clades=list(
 Spermatophyta=c("Angiospermae", "Acrogymnospermae"),

  Acrogymnospermae = c("Cycadales", "Gnetales"),
	Coniferae=c("Pinales", "Cupressales"),
		Cupressophyta=c("Cupressales", "Araucariales"),
		
	Cycadophyta=c("Cycas", "Zamia"),
	Gnetophyta=c("Gnetales", "Welwitschiales"),
	

  Angiospermae=c("Amborellales", "Asterales"),

	Mesangiospermae=c("Asterales", "Chloranthales"),

		Monocotyledoneae=c("Acorales", "Commelinales"),
			Nartheciidae=c("Alismatales", "Liliales"),
				Petrosaviidae=c("Liliales", "Commelinales"),
					Commelinidae=c("Asparagales", "Poales"),
					
		Magnoliidae = c("Magnoliales", "Piperales"),

#		Austrobaileyales = c("Trimeniaceae", "Schisandraceae", "Austrobaileyaceae"),

		Eudicotyledoneae = c("Ranunculales", "Fagales"),
		   Gunneridae=c("Gunnerales", "Fagales"),
			Superasteridae=c("Dilleniales", "Lamiales"),
				Pentapetalae=c("Asteridae", "Santalales"),
					Asteridae=c("Cornales", "Lamiales"),
						Lamiidae = c("Lamiales", "Garryales"),
						Campanulidae = c("Asterales", "Aquifoliales"),
			Superrosidae=c("Saxifragales","Fagales"),
				Rosidae=c("Vitales", "Fagales"),
					Fabidae = c("Celastrales", "Fagales"), 	
					Malvidae = c("Myrtales", "Malvales")
)

.fix_classification.nescent=function(csv=""){
#	e.g., csv="genus_order_lookup.11162011.csv"
	dat=read.csv(csv, as.is=TRUE)
	if(all(names(dat)==c("GeneraFromWill","genus_final","family_final","order_final","RESOURCE","FAMILY_REFERENCE","ORDER_REFERENCE","NOTES"))){
		NESCent_APGIII_classification=dat[,c("genus_final","family_final","order_final")]
		names(NESCent_APGIII_classification)=c("genus","family","order")
		attr(NESCent_APGIII_classification, "date")=Sys.Date()
		attr(NESCent_APGIII_classification, "file")=csv
		save(NESCent_APGIII_classification, file="NESCent_APGIII_classification.rda")
	} else {
		stop()
	}
}

.fix_classification.orders=function(csv=""){
#	e.g., csv="APGIII_orders.2011-11-17.csv"
	dat=read.csv(csv, as.is=TRUE)
	notes=dat$notes
	names(notes)=dat$family
	dat=dat[,c("family","order")]
	dat$order[dat$order==""]="UNPLACED"
	spermatophyta_APGIII_orders=split(dat$family, dat$order)
	attr(spermatophyta_APGIII_orders, "date")=Sys.Date()
	attr(spermatophyta_APGIII_orders, "file")=csv
	attr(spermatophyta_APGIII_orders, "notes")=notes
	save(spermatophyta_APGIII_orders, file="spermatophyta_APGIII_orders.rda")
	return(spermatophyta_APGIII_orders)
}

