## THIS IS JUST A SIMPLE CODE TO GET THE SlOPES OF A LINEAR RELATIONSHIP BETWEEN
## CONTACT ANGLE AND THE POSITION ALONG THE PROBOSCIS. THIS MEASURE WILL BE USED 
## AS A RATE OF CHANGE OF CONTACT ANGLE. VALUES WERE CALCULATE AND THEN COPIED
## AND PASTED ON THE FILE WITH SPECIES MEANS

ca_slope = read.csv(file = "data/hawkmoth-ca-slopes.csv", h = T, sep = ";")

str(ca_slope)

ca_slope$species = as.factor(ca_slope$species)
ca_slope$normalized.distance = as.numeric(ca_slope$normalized.distance)

agrius = ca_slope[ca_slope$species=="Agrius_cingulata",]

agrius.m = lm(CA ~ normalized.distance , data = agrius)
agrius.m

max(agrius$CA)/min(agrius$CA)

dolba = ca_slope[ca_slope$species=="Dolba_hyloeus",]

dolba.m = lm(CA ~ normalized.distance, data = dolba)
dolba.m

max(dolba$CA)/min(dolba$CA)

darapsa = ca_slope[ca_slope$species=="Darapsa_myron",]

darapsa.m = lm(CA ~ normalized.distance, data = darapsa)
darapsa.m

max(darapsa$CA)/min(darapsa$CA)

fasciatus = ca_slope[ca_slope$species=="Eumorpha_fasciatus",]

fasciatus.m = lm(CA ~ normalized.distance, data = fasciatus)
fasciatus.m

max(fasciatus$CA)/min(fasciatus$CA)

pandorus = ca_slope[ca_slope$species=="Eumorpha_pandorus",]

pandorus.m = lm(CA ~ normalized.distance, data = pandorus)
pandorus.m

max(pandorus$CA)/min(pandorus$CA)

enyo = ca_slope[ca_slope$species=="Enyo_lugubris",]

enyo.m = lm(CA ~ normalized.distance, data = enyo)
enyo.m

max(enyo$CA)/min(enyo$CA)

hemaris = ca_slope[ca_slope$species=="Hemaris_thysbe",]

hemaris.m = lm(CA ~ normalized.distance, data = hemaris)
hemaris.m

max(hemaris$CA)/min(hemaris$CA)

quinque = ca_slope[ca_slope$species=="Manduca_quinquemaculata",]

quinque.m = lm(CA ~ normalized.distance, data = quinque)
quinque.m

max(quinque$CA)/min(quinque$CA)

rustic = ca_slope[ca_slope$species=="Manduca_rustica",]

rustic.m = lm(CA ~ normalized.distance, data = rustic)
rustic.m

max(rustic$CA)/min(rustic$CA)

sexta = ca_slope[ca_slope$species=="Manduca_sexta",]

sexta.m = lm(CA ~ normalized.distance, data = sexta)
sexta.m

max(sexta$CA)/min(sexta$CA)

plebeian = ca_slope[ca_slope$species=="Paratrea_plebeja",]

plebeian.m = lm(CA ~ normalized.distance, data = plebeian)
plebeian.m

max(plebeian$CA)/min(plebeian$CA)

tersa = ca_slope[ca_slope$species=="Xylophanes_tersa",]

tersa.m = lm(CA ~ normalized.distance, data = tersa)
tersa.m

max(tersa$CA)/min(tersa$CA)

hyles = ca_slope[ca_slope$species=="Hyles_lineata",]

hyles.m = lm(CA ~ normalized.distance, data = hyles)
hyles.m

max(hyles$CA)/min(hyles$CA)

