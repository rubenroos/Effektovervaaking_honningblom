# Estimering av størrelsen av fire honningblom populasjonser basert på
# linje og rose design, der det på grunn av datagrunnlaget benyttes ulike
# tilnærminger for de fire populasjonene.

# Fix problems with norwegian names
Sys.setlocale(locale='no_NB.utf8')

# Last inn nødvendige R-pakker, eget funksjonsbibliotek og data.

library(boot)
library(rio)

################################################################################
# Funksjons-bibliotek
################################################################################

estim.fa.linje <- function(lengde,bredde,andeler) {
  return(lengde*bredde*andeler)
}

estim.fa.rose <- function(ua,andeler) {
  return(ua*andeler)
}

estim.fa.tetthet <- function(x) {
  x[x==0] <- NA
  return(colMeans(x,na.rm = TRUE))
}

ratioest <- function(z){
  ll <- length(z)
  y <- z[1:(ll/2)]
  x <- z[(1+(ll/2)):ll]
  if (sum(x,na.rm = T) == 0) {return(0)}
  if (sum(y,na.rm = T) == 0) {return(0)}
  
  naivest <- mean(y,na.rm = T)/mean(x,na.rm = T)
  if (is.na(naivest) | is.nan(naivest)) {return(0)}
  
  corrterm <- try(cov((y/x),x,use="pairwise.complete.obs")/mean(x,na.rm = T),silent = T)
  if (!is.numeric(corrterm) | is.na(corrterm) | is.nan(corrterm)) {return(naivest)}
  if ((corrterm > naivest) | (naivest - corrterm > 1)) {return(naivest)}
  testsig <- try(cor.test(y/x,x)[[3]],silent = T)
  if (!is.numeric(testsig) | is.na(testsig) | is.nan(testsig)) {return(naivest)}
  if (testsig > 0.05) {return(naivest)}
  return(naivest - corrterm)
}

ratioest2 <- function(z){
  ll <- length(z)
  y <- z[1:(ll/2)]
  x <- z[(1+(ll/2)):ll]
  if (sum(x,na.rm = T) == 0) {return(0)}
  if (sum(y,na.rm = T) == 0) {return(0)}
  
  naivest <- mean(y,na.rm = T)/mean(x,na.rm = T)
  if (is.na(naivest) | is.nan(naivest)) {return(0)}
  
  corrterm <- try(cov((y/x),x,use="pairwise.complete.obs")/mean(x,na.rm = T),silent = T)
  if (!is.numeric(corrterm) | is.na(corrterm) | is.nan(corrterm)) {return(naivest)}
  if ((corrterm > naivest) | (naivest - corrterm > 1)) {return(naivest)}
  testsig <- try(cor.test(y/x,x)[[3]],silent = T)
  if (!is.numeric(testsig) | is.na(testsig) | is.nan(testsig)) {return(naivest)}
  if (testsig > 0.05) {return(naivest)}
  return(naivest)
}

nbinom.likelihood <- function(x,n.transektruter,n.forekomstruter,antall,
                              n.analyseruter.per.transektruter) {
  
  # Veid maksimum likelihood funksjon negativ binomisk fordeling
  
  cond.dnbinom <- function(y,x) {
    return(dnbinom(y,size=x[1],prob=x[2]))#/(1-dnbinom(0,size=x[1],prob=x[2])))
  }
  
  x[1] <- exp(x[1])
  x[2] <- inv.logit(x[2])
  transektledd <- NA
  if (n.analyseruter.per.transektruter == 1) {
    transektledd <- (log(dnbinom(0,size=x[1],prob=x[2]))*(n.transektruter - n.forekomstruter)) + 
      ((log(1-dnbinom(0,size=x[1],prob=x[2])))*n.forekomstruter)
  }
  if (n.analyseruter.per.transektruter == 2) {
    transektledd <- (log(dnbinom(0,size=x[1],prob=x[2])^2)*(n.transektruter - n.forekomstruter)) + 
      ((log(1-(dnbinom(0,size=x[1],prob=x[2])^2)))*n.forekomstruter)
  }
  
  analyseledd <- sum(log(sapply(antall,cond.dnbinom,x)))
  
  return(-1*((n.transektruter*transektledd) + (length(antall)*analyseledd))/(n.transektruter+length(antall)))
}

param.est.andel <- function(serie = NULL,
                            ant.trans.ruter = NULL,
                            napt = 1,
                            n.forekomstruter = 0) {
  
  # Forsøker å reskalere rutestørrelsene for hvert år slik at
  # antall planter i hver analyserute er i samsvar med en
  # binomisk fordeling. Dette er nødvendig for enkelte år i Teneskjær
  # populasjonen der det er ekstremt mange planter i enkelte av analyse-
  # rutene, som medfører at ingen aktuelle fordelingsmodeller (poisson,
  # ZIP, hurdeled poisson, binomisk, ZI neg.binomisk osv.) fungerer som en
  # modell for dataene.
  
  # Reskalering av rutene innebærer å fordele de observerte
  # plantene i ei analyserute over sc.f reskalerte ruter og øke
  # antall transektruter (uten planter) med den samme faktoren sc.f.
  
  # Denne manipulasjonen antar implisitt at plantene er jevnt fordelt innenfor
  # analyserutene. Antagelsen har antagelig liten btydning sammenliknet med
  # antagelsen om uavhengighet mellom transektrutene, som er den kritiske 
  # antagelsen i denne sammenheng.
  
  sc.f <- 10^((max(serie,na.rm = TRUE) %/% c(300) > 10) +
                (max(serie,na.rm = TRUE) %/% c(300) > 0))
  
  antall <- rep(ceiling(serie/sc.f),sc.f)
  antall <- antall[!is.na(antall)]
  mini.nbinom <- nlminb(c(1,1), 
                        objective=nbinom.likelihood,
                        n.transektruter=ant.trans.ruter*sc.f,
                        n.forekomstruter=n.forekomstruter,
                        antall=antall,
                        n.analyseruter.per.transektruter=napt)
  mini.nbinom$par <- c(exp(mini.nbinom$par[1]),inv.logit(mini.nbinom$par[2]))
  
  andel <- pnbinom(0,mini.nbinom$par[1],mini.nbinom$par[2],lower.tail = FALSE)
  return(andel)
}

estim.ua.rose <- function(x)
  # Tilsvarende funksjon hos Olav, "beregn.arealvekt", beregner arealet til 
  # de ulike delarealene i rose-designet. Denne funksjonen er her endret slik 
  # den beregner delen av det undersøkte området ved lokaliteten som det ett
  # transekt i rose-designet representerer.
  
{
  r1 <- 0             # radius for indre sirkel (transektets startpunkt)
  r2 <- x             # radius for ytre sirkel (transektets endepunkt)
  return( (pi*r2^2 - pi*r1^2)/8 )
}

estim.pop.size <- function(transektruteserier = NULL,
                           analyseruteserier = NULL,
                           transects = NULL,
                           transektlengder = NULL,
                           design = NULL,
                           n.analyseruter.per.transektruter = NULL,
                           transekt.rute.lengde = NULL,
                           loc = NULL,
                           transektbredde = NULL,
                           years = NULL,
                           grense = 0) {
  
  # Estimer tetthet av individer, fertile- og vegetative planter innenfor
  # forekomstarealet som gjennomsnittlig tetthet i analyseruter med planter.
  
  # Estimer populasjonsstørrelsene for hvert år som 
  # forekomstareal * tettheten  innenfor forekomstarealet.
  
  if (design == "Linje") {
    
    
    andeler <- colMeans(transektruteserier)
    zero.fa <- andeler <= grense
    naf <- colSums(transektruteserier)
    
    for (year in years[zero.fa]) {
      # Parameter-basert estimering av forekomstarealet
      andeler[year] <- param.est.andel(serie = analyseruteserier[,year],
                                       ant.trans.ruter = dim(transektruteserier)[1],
                                       napt = n.analyseruter.per.transektruter[loc],
                                       n.forekomstruter = naf[year])
    }
    
    fa <- estim.fa.linje(lengde = dim(transektruteserier)[1] * transekt.rute.lengde[loc],
                         bredde = transektbredde,
                         andeler = andeler)
    
    popsize <- round(estim.fa.tetthet(analyseruteserier) * fa)
    popsize[is.nan(popsize)] <- 0
    
  } else { # design == "Rose"
    ua <- 0
    for (tname in unique(transects$names)) {
      transl <- transektlengder[transects$names == tname][1]
      ua <- ua + estim.ua.rose(transl)
    }
    
    andeler <- colMeans(transektruteserier)
    zero.fa <- andeler <= grense
    naf <- colSums(transektruteserier)
    for (year in years[zero.fa]) {
      # Parameter-basert estimering av forekomstarealet
      andeler[year] <- param.est.andel(serie = analyseruteserier[,year],
                                       ant.trans.ruter = dim(transektruteserier)[1],
                                       napt = n.analyseruter.per.transektruter[loc],
                                       n.forekomstruter = naf[year])
      
    }
    
    fa <- estim.fa.rose(ua = ua,
                        andeler = andeler)
    popsize <- round(estim.fa.tetthet(analyseruteserier) * fa)
    popsize[is.nan(popsize)] <- 0
  }
  return(list(popsize = popsize, fa = fa))
}

################################################################################
# Data
################################################################################

#lokalitetsdata <- import("Data/Honningblom_2022.xlsx",sheet="Lokalitetsdata",na=c("NA"))
#transektdata <- import("Data/Honningblom_2022.xlsx",sheet="transektdata",na=c("NA"," "))
#rutedata <- import("Data/Honningblom_2022.xlsx",sheet="Honningblom_data",na=c("NA"))

###### Ny kode #################
getwd()
lokalitetsdata <- import("Data/Honningblom_2024_copy.xlsx",sheet="Lokalitetsdata",na=c("NA"))
transektdata <- import("Data/Honningblom_2024_copy.xlsx",sheet="transektdata",na=c("NA"," "))
rutedata <- import("Data/Honningblom_2024_copy.xlsx",sheet="Honningblom_data",na=c("NA"))

#lokalitetsdata$År <- lokalitetsdata$Year # Går tilbake til opprinnelige variabelnavnene fra i fjor
#rutedata$Ar <- rutedata$Year
#transektdata$År <- transektdata$Year

sort(unique(rutedata$Rute.ID))

# Allokerer nye ID til ruter som mangler ID

rutedata$Rute.ID_orig <- rutedata$Rute.ID # Tar vare på opprinnelig input

rutedata$Rute.ID[rutedata$Rute == "SKI.10.A"] <- "SKI10"
rutedata$Rute.ID[rutedata$Rute == "SKI.5.A"] <- "SKI5"
rutedata$Rute.ID[rutedata$Rute == "SKX.2.A"] <- "SKJ51"
rutedata$Rute.ID[rutedata$Rute == "SKX.4.A"] <- "SKJ52"
rutedata$Rute.ID[rutedata$Rute == "SKX.5.A"] <- "SKJ53"
rutedata$Rute.ID[rutedata$Rute == "SKX.9.A"] <- "SKJ54"
rutedata$Rute.ID[rutedata$Rute == "SKY.1.A"] <- "SKJ55"
rutedata$Rute.ID[rutedata$Rute == "SKY.5.A"] <- "SKJ56"
rutedata$Rute.ID[rutedata$Rute == "SKY.6.A"] <- "SKJ57"
rutedata$Rute.ID[rutedata$Rute == "SKY.8.A"] <- "SKJ58"

rutedata$Rute_orig <- rutedata$Rute # Tar vare på opprinnelig input
rutedata$Rute <- rutedata$Rute.ID # Kopierer variabler for å unngå rekoding senere i skriptet

table(sort(rutedata$Rute)) # Sjekk

###### Slutt ny kode #################


################################################################################
# Script
################################################################################

# Les, initier og beregn designvariabler - felles kode for alle lokaliteter.

roser <- unique(lokalitetsdata$Lokalitet[lokalitetsdata$Design == "Rose"])
linjer <- unique(lokalitetsdata$Lokalitet[lokalitetsdata$Design == "Linje"])
forekomst.transekt <- rep("",length(c(roser,linjer)))
transekt.rute.lengde <- transekt.rute.bredde <- rep(0,length(c(roser,linjer)))
names(forekomst.transekt) <- names(transekt.rute.lengde) <- names(transekt.rute.bredde) <- c(roser,linjer)
forekomst.transekt[roser] <- "Ant.0,5 m.med"
forekomst.transekt[linjer] <- "Ant.m.med"
transekt.rute.lengde[roser] <- transekt.rute.bredde[roser] <- 0.5
transekt.rute.lengde[linjer] <- 1
transekt.rute.bredde[linjer] <- 0.5
analyse.rute.lengde <- transekt.rute.bredde
forekomst.avstand <- "Forekomst honningblom (m)"
n.analyseruter.per.transektruter <- transekt.rute.lengde / analyse.rute.lengde

localities <- unique(lokalitetsdata$Lokalitet)
nboot <- 2000  # Antall simuleringer
grense <- 0.05 # Når andelen transektruter med Honningblom er mindre enn grense,
# gjennomføres parametrisk estimering av forekomstarealet. Kan 
# settes lavere når en har mange transektruter, og individuelt 
# for hver populasjon. Men da må scriptet endres tilsvarende
# på linje 365 og 426.
pop <- list() # Resultat-objekt

for (loc in localities) {
  
  # Organiser transektdatasettet i ei transektrute x årMedTransektdata matrise
  # 1) Finn alle transekter og gi dem navn
  # 2) Finn transektenes lengde, utled antall ruter langs hvert transekt og
  #    navngi rutene.
  # 3) Lag matrisa og fyll den med presence/absence data
  
  transectyears <- sort(unique(lokalitetsdata$Year[lokalitetsdata$Lokalitet == loc]))
  design <- unique(lokalitetsdata[lokalitetsdata$Lokalitet==loc,"Design"])
  transects <- transektdata[transektdata$Lokalitet == loc,]
  transects$names <- paste(transects$Lokalitet,transects$Transekt,sep = "_")
  transektbredde <- lokalitetsdata[lokalitetsdata$Lokalitet==loc,"Transektbredde"] [1]
  
  transektlengder <-  as.numeric(transects$`Lengde (m)`)
  transektrutenavn <- NULL
  for (tname in unique(transects$names)) {
    rutesequence <- 1:as.integer(transektlengder[transects$names == tname][1]/
                                   transekt.rute.lengde[loc])
    transektrutenavn <- c(transektrutenavn,
                          paste(tname,
                                (rutesequence - 1)*transekt.rute.lengde[loc],
                                sep = "_"))
  }
  
  transektruteserier <- matrix(0,nrow = length(transektrutenavn),ncol = length(transectyears))
  dimnames(transektruteserier) <- list(transektrutenavn,paste(transectyears))
  
  for (tname in unique(transects$names)) {
    for (year in transectyears) {
      forekomstruter <- 
        as.numeric(unlist(sapply(transects[transects$names == tname &
                                             transects$Year == year,
                                           forekomst.avstand],strsplit,split=',')))
      if (sum(is.na(forekomstruter))==0) {
        forekomstruter <- paste(tname,forekomstruter,sep = "_")
        transektruteserier[forekomstruter,as.character(year)] <- 1
      } 
    }
  }
  
  # Organiser analysedatasettet i ei analyserute x årMedAnalyser matrise
  # Tar i utgangspunktet med alle analyserutene fra lokaliteten, fjerner
  # deretter rutene som aldri i løpet av tidsserien inngår i forekomstarealet
  # og år uten data. Denne prosedyren kan virke unødvendig, men er trolig
  # nødvendig for å inkludere analyserutene som i starten var utenfor 
  # forekomstarealet, men som senere blir en del av det.
  
  analyse.rute.data <- rutedata[rutedata$Lokalitet == loc,]
  
  ####################################################################################
  #### Her er roten til alt ondt i 2022! Koden nedenfor er en ad hoc tilpasning til fjorårets
  #### input data der ingen av variablene "Rute" og "Rute.ID" var komplette.
  #### Koden under fungerte bra i fjor, men ikke i år fordi variabelen Rute nå ikke er oppgitt
  #### for flere analyseruter i 2022-dataene. Koden nedenfor introduserer dermed i effekt 
  #### flere nye ruter som bare er observert i 2022. Dette skaper igjen problemer i 
  #### bootstrap rutinen nedenfor fordi sjansene for at bootstrapping gir datasett der
  #### alle ruter har 0 individer øker betraktelig i forhold til i fjor. Dermed oppstår det
  #### problemer når en under likelihood-estimeringen leter etter optima vha. funksjonen
  #### nlminb (linje 118).
  #### Koden på linjene 231 - 245 medfører at ingen ruter mangler navn slik at koden på
  #### linjer 343 - 344 blir overflødig. Men jeg lar den likevel stå.
  
  mangler.navn <- is.na(analyse.rute.data$Rute)
  analyse.rute.data$Rute[mangler.navn] <- analyse.rute.data$Rute.ID[mangler.navn]
  
  #### Slutt på alt ondt!
  ####################################################################################
  
  analyseruter <- sort(unique(analyse.rute.data$Rute))
  years <- sort(unique(analyse.rute.data$Year))
  
  analyseruteserier <- matrix(NA,nrow = length(analyseruter),
                              ncol=length(years))
  dimnames(analyseruteserier)[[1]] <- analyseruter
  dimnames(analyseruteserier)[[2]] <- years
  fertilitetsserier <- sterilitetsserier <- analyseruteserier
  tetthet <- round(analyse.rute.data$Ant.individer/analyse.rute.data$Rutestrl)
  fert.tetthet <- round(analyse.rute.data$Ant.fertile/analyse.rute.data$Rutestrl)
  veg.tetthet <- tetthet - fert.tetthet
  
  for (i in 1:length(tetthet)) {
    analyseruteserier[analyse.rute.data$Rute[i],
                      as.character(analyse.rute.data$Year[i])] <- tetthet[i]
    fertilitetsserier[analyse.rute.data$Rute[i],
                      as.character(analyse.rute.data$Year[i])] <- fert.tetthet[i]
    sterilitetsserier[analyse.rute.data$Rute[i],
                      as.character(analyse.rute.data$Year[i])] <- veg.tetthet[i]
  }
  
  ta.med <- rowSums(analyseruteserier,na.rm = T) > 0
  analyseruteserier <- analyseruteserier[ta.med,]
  fertilitetsserier <- fertilitetsserier[ta.med,]
  sterilitetsserier <- sterilitetsserier[ta.med,]
  notOnlyNAs <- function(x) {return(length(x[!is.na(x)]) > 0)}
  ta.med <- apply(analyseruteserier,2,notOnlyNAs)
  analyseruteserier <- analyseruteserier[,ta.med]
  fertilitetsserier <- fertilitetsserier[,ta.med]
  sterilitetsserier <- sterilitetsserier[,ta.med]
  
  years <- dimnames(analyseruteserier)[[2]]
  
  
  minsize <- as.integer(years)
  names(minsize) <- years
  minfert <- minsize
  for(year in years) {
    minsize[year] <- sum(analyse.rute.data$Ant.individer[analyse.rute.data$Year == as.integer(year)],na.rm = T)
    minfert[year] <- sum(analyse.rute.data$Ant.fertile[analyse.rute.data$Year == as.integer(year)],na.rm = T)
  }
  
  # For Skjellvik populasjonen gjøres antagelsen at transekt-datasettet fra 2021 
  # også reflekterer forekomstarealet gjennom hele perioden 2014 - 2020. 
  # Denne antagelsen er ikke rimelig for Teneskjær. Her er det 
  # dessuten ikke mulig å estimere et forekomstareal ut fra transektdatasettet.
  
  # For Skjellvik populasjonen estimeres forekomstarealet for årene 2014-2017 
  # og 2019 - 2020 fra transektdatasettet for 2021. I praksis, kopieres 
  # transektdataene for 2021 for disse 6 årene. For Filletassen, Skjellvik og 
  # Skipstadsand estimeres forekomstarealet som andel transektruter med 
  # honningblom * størrelse av undersøkt areal. Størrelse av undersøkt areal er
  # i linje-designet ved Skjellvik, Skipstadsand og Teneskjær lik 
  # transektbredden * summerte transekt-lengder. I rose-designet settes 
  # størrelsen av det undersøkte arealet lik summen av sektorenes areal, der
  # kun de reelle sektorene telles med, dvs. de som har ei transektrute.
  
  # For Teneskjær ser det ut til at forekomstarealet er blitt redusert i løpet
  # perioden mellom 2014 og 2021. Det antas at det ikke har skjedd endringer i
  # transektrutene gjennom perioden. Ei negativ binomisk fordeling tilpasses
  # datasett fra hvert år, dvs. analysedata fra det aktuelle året samt
  # transektdataene, basert på en veid likelihood funksjon. Vektene i funksjonen
  # er antall observasjoner i årets analysedatasett og i transektdataene.
  # Forekomstarealet ved Teneskjær for et aktuelt år estimeres som 
  # fat * størrelsen av det undersøkte arealet: 
  
  # fat <- pnbinom(0,par[1],par[2],lower.tail = FALSE)
  
  # Dette vil medføre at det estimerte forekomstarealet vil gå ned i løpet av 
  # perioden fordi antall individer i flere av analyserutene er redusert til 0.
  
  # For Skipstadsand-populasjonen benyttes samme tilnærming som for populasjonen
  # på Teneskjær. Ved skipstadsand er det kun to forekomstruter i transektene,
  # noe som gir svært usikre estimat av forekomstarealet. En parametrisk
  # tilnærming reduserer usikkerheten noe sammenliknet med en ikke-parametrisk
  # bootstrap tilnærming, naturlig nok, men usikkerheten er fortsatt svært stor.
  
  # Den negative binomiske fordelingen (NB) ser ut til å fungere bedre som modell for
  # den romlige fordelingen av Honningblom enn andre aktuelle fordelinger som er
  # utprøvd: poisson, hurdeled poisson, zeroinflated poisson, hurdeled NB og 
  # zeroinflated NB.
  
  if (loc %in% c("Skipstadsand","Teneskjaer","Skjellvik")) {
    years.without.transect <- years[which(!(years %in% transectyears))]
    first.year.with.transect <- min(transectyears)
    ex.transektruteserier <- matrix(transektruteserier[,as.character(first.year.with.transect)],
                                    nrow = dim(transektruteserier)[1], 
                                    ncol = length(years.without.transect))
    transektruteserier <- cbind(ex.transektruteserier,transektruteserier)
    dimnames(transektruteserier)[[2]] <- paste(c(years.without.transect,transectyears))
  }
  
  popsize <- matrix(0,nrow = nboot+1, ncol = length(years))
  dimnames(popsize)[[2]] <- years
  forekomstareal <- vegsize <- fertsize <- popsize
  
  returned.stuff <- estim.pop.size(transektruteserier = transektruteserier,
                                   analyseruteserier = analyseruteserier,
                                   transects = transects,
                                   transektlengder = transektlengder,
                                   design = design,
                                   n.analyseruter.per.transektruter = 
                                     n.analyseruter.per.transektruter,
                                   transekt.rute.lengde = transekt.rute.lengde,
                                   loc = loc,
                                   transektbredde = transektbredde,
                                   years = years,
                                   grense = grense)
  
  popsize[1,] <- round(returned.stuff$popsize)
  popsize[1,popsize[1,] < minsize] <- minsize[popsize[1,] < minsize]
  forekomstareal[1,] <- returned.stuff$fa
  
  fertsize[1,] <- round(popsize[1,] * apply(rbind(fertilitetsserier,analyseruteserier),2,ratioest))
  #fertsize[1,] <- estim.pop.size(transektruteserier = transektruteserier,
  #                                       analyseruteserier = fertilitetsserier,
  #                                       transects = transects,
  #                                       transektlengder = transektlengder,
  #                                       design = design,
  #                                       n.analyseruter.per.transektruter = 
  #                                         n.analyseruter.per.transektruter,
  #                                       transekt.rute.lengde = transekt.rute.lengde,
  #                                       loc = loc,
  #                                       transektbredde = transektbredde,
  #                                       years = years,
  #                                       grense = grense)$popsize
  fertsize[1,is.na(fertsize[1,])] <- 0
  fertsize[1,fertsize[1,] < minfert] <- minfert[fertsize[1,] < minfert]
  
  vegsize[1,] <- popsize[1,] - fertsize[1,] 
  
  # Estimer usikkerheten gjennom permutasjoner/ikke-parametrisk bootstrap.
  # For Skjellvik bootstrap transektrutetidsseriene og estimer
  # forekomstarealet for hvert bootstrapsample (og år). Tilsvarende, bootstrap
  # analysetidsseriene og estimer tetthet for hvert bootstrapsample (og år).
  # For Filletassen bootstrap transektrutetidsseriene og beregn forekomstarealet
  # for hver permutasjon og år uten å endre undersøkelsesområdet i hvert sample. 
  # Bootstrap analysetidsseriene og estimer tetthet for hvert bootstrapsample 
  # (og år).
  # For Teneskjær bootstrap analysetidsseriene og estimer forekomstarealet og 
  # tettheten for hver permutasjon og år.
  # For Skipstadsand bootstrap transektrutetidsseriene, bootstrap 
  # analysetidsseriene og estimer forekomstarealet for hvert år fra hvert par av 
  # bootstrapsampler. Benytt bootstrap sample fra analysetidsseriene til også
  # å estimere tetthet innenfor forekomstarealet for hvert år.
  # I Praksis, slik estimeringen nå er programmert, kan de samme resampling-
  # rutinene benyttes for alle populasjonene.
  #
  # For hver simulering beregn populasjonsstørrelsen for hvert år og estimer 
  # konfidens-intervallet fra fordelingen av populasjonsstørrelser for hvert år.
  # Dog med den begrensning at populasjonsstørrelsen er større eller lik 
  # antallet individer observert i analyserutene
  # 
  # Tilsvarende for estimering av forekomstareal og antall fertile planter.
  
  for (i in 1:nboot) {
    boot.names <- sample(dimnames(transektruteserier)[[1]],replace = T) 
    transektruteserier.boot <- transektruteserier[boot.names,]
    boot.names <- sample(dimnames(analyseruteserier)[[1]],replace = T)
    analyseruteserier.boot <- analyseruteserier[boot.names,]
    fertilitetsserier.boot <- fertilitetsserier[boot.names,]
    returned.stuff <- estim.pop.size(transektruteserier = transektruteserier.boot,
                                     analyseruteserier = analyseruteserier.boot,
                                     transects = transects,
                                     transektlengder = transektlengder,
                                     design = design,
                                     n.analyseruter.per.transektruter = 
                                       n.analyseruter.per.transektruter,
                                     transekt.rute.lengde = transekt.rute.lengde,
                                     loc = loc,
                                     transektbredde = transektbredde,
                                     years = years,
                                     grense = grense)
    popsize[i+1,] <- round(returned.stuff$popsize)
    popsize[i+1,popsize[i+1,] < minsize] <- minsize[popsize[i+1,] < minsize]
    forekomstareal[i+1,] <- returned.stuff$fa
    fertsize[i+1,] <- round(popsize[i+1,] * apply(rbind(fertilitetsserier.boot,analyseruteserier.boot),2,ratioest))
    #  fertsize[i+1,] <- estim.pop.size(transektruteserier = transektruteserier.boot,
    #                                   analyseruteserier = fertilitetsserier.boot,
    #                                   transects = transects,
    #                                   transektlengder = transektlengder,
    #                                   design = design,
    #                                   n.analyseruter.per.transektruter = 
    #                                     n.analyseruter.per.transektruter,
    #                                   transekt.rute.lengde = transekt.rute.lengde,
    #                                   loc = loc,
    #                                   transektbredde = transektbredde,
    #                                   years = years,
    #                                   grense = grense)$popsize
    fertsize[i+1,is.na(fertsize[i+1,])] <- 0
    fertsize[i+1,fertsize[i+1,] < minfert] <- minfert[fertsize[i+1,] < minfert]
    vegsize[i+1,] <- popsize[i+1,] - fertsize[i+1,]
    message(paste("Lokalitet:",loc,", Permutasjon:",i))
  }
  sizecofint <- round(apply(popsize[2:nboot+1,],2,quantile,probs =c(0.025,0.975)))
  fertcofint <- round(apply(fertsize[2:nboot+1,],2,quantile,probs =c(0.025,0.975)))
  vegcofint <- round(apply(vegsize[2:nboot+1,],2,quantile,probs =c(0.025,0.975)))
  forekomstcofint <- apply(forekomstareal[2:nboot+1,],2,quantile,probs =c(0.025,0.975))
  pop[[loc]][["size"]] <- rbind(popsize[1,],sizecofint)
  pop[[loc]][["fert"]] <- rbind(fertsize[1,],fertcofint)
  pop[[loc]][["veg"]] <- rbind(vegsize[1,],vegcofint)
  pop[[loc]][["forekomst"]] <- rbind(forekomstareal[1,],forekomstcofint)
  dimnames(pop[[loc]][["size"]])[[1]][1] <- "est"
  dimnames(pop[[loc]][["fert"]])[[1]][1] <- "est"
  dimnames(pop[[loc]][["veg"]])[[1]][1] <- "est"
  dimnames(pop[[loc]][["forekomst"]])[[1]][1] <- "est"
}

print(pop)
sink('pop2024.txt')
print(pop)
sink()

# Figurer

years <- sizes <- NULL
for (loc in localities) {
  years <- c(years, as.integer(dimnames(pop[[loc]]$size)[[2]]))
  sizes <- c(sizes, log(pop[[loc]]$size["97.5%",]))
}
xlimits <- c(min(years) - 0.5, max(years) + 0.5)
ylimits <- c(0, max(sizes) + 0.5)

#if(!dir.exists("Figurer")) dir.create("Figurer")
#pdf(paste("Figurer/","Alle lokaliteter log",".pdf",sep=""))

#pdf("Figures/Alle lokaliteter 2022.pdf") 
png("Figures/Alle lokaliteter 2024.png", width = 840, height = 840)
par(mfrow=c(2,2))
for (loc in localities) {
  plot(as.integer(dimnames(pop[[loc]]$size)[[2]]), 
       log(pop[[loc]]$size["est",]),
       col=1,pch=19,cex=1,type="o",
       ylim=ylimits,xlim=xlimits,
       xlab = expression(Aar), 
       ylab = "Antall individer", axes = F,cex.lab=1.25,
       main=loc)
  axis(1,
       at=c(xlimits[1],(xlimits[1]+0.5):(xlimits[2]-0.5)),
       labels = c("",as.character((xlimits[1]+0.5):(xlimits[2]-0.5))),
       pos=0,
       cex.lab=1.25, 
       cex.axis=1.25)
  axis(2,
       at=c(0,log(c(10,100,1000,10000,exp(ylimits[2])))),
       labels = c("",as.character(c(10,100,1000,10000)),""),
       pos=xlimits[1],
       cex.lab=1.25, 
       cex.axis=1.25)
  polygon(c(as.integer(dimnames(pop[[loc]]$size)[[2]]),rev(as.integer(dimnames(pop[[loc]]$size)[[2]]))),
          c(log(pop[[loc]]$size["97.5%",]),rev(log(pop[[loc]]$size["2.5%",]))),
          col=rgb(0.5,0.5,0.5,alpha=0.2),border=NA)
  points(as.integer(dimnames(pop[[loc]]$fert)[[2]]), 
         log(pop[[loc]]$fert["est",]),col=2,pch=19,cex=1,type="o")
  polygon(c(as.integer(dimnames(pop[[loc]]$fert)[[2]]),rev(as.integer(dimnames(pop[[loc]]$fert)[[2]]))),
          c(log(pop[[loc]]$fert["97.5%",]),rev(log(pop[[loc]]$fert["2.5%",]))),
          col=rgb(1,0,0,alpha=0.2),border=NA)
}

dev.off()


#### Check if it is possible to make some ggplots ####
dev.off()
library(reshape2)
library(tidyverse)
library(gganimate)
library(av)
library(gifski)

#Get the data (list format) into a dataframe
df.pop <- data.frame(melt(pop))
df.pop <- df.pop %>% 
  rename("Estimate" = Var1, 
         "Year" = Var2,
         "Value" = value,
         "Pop.variable" = L2, 
         "Locality" = L1) %>% 
  mutate(Estimate = recode(Estimate, est = "est",
                                                '2.5%' = "lower",
                                                '97.5%' = "upper")) %>% 
  mutate_if(is.numeric, round, 1)

#restructure the data so that each estimate is in its own column
df.pop <- reshape(data = df.pop,
        idvar = c("Locality", "Pop.variable", "Year"),
        v.names = c("Value"),
        timevar = "Estimate",
        direction = "wide") %>% 
  mutate(Year = as.integer(Year))

str(df.pop)

p.TEN <- df.pop %>% 
  dplyr::filter(Pop.variable == "size") %>% 
  dplyr::filter(Locality == "Teneskjaer") %>% 
  ggplot(., aes(x = Year, y = Value.est, ymin = Value.lower, ymax = Value.upper)) + 
  geom_ribbon(alpha = 0.2, fill = "lightblue") + 
  geom_point(color = "darkblue")+
  geom_line(color = "darkblue")+
  scale_x_continuous(breaks = 2014:2024, labels = 2014:2024, expand = c(0,0.1)) +
  scale_y_continuous(expand = c(0, 0))+
  theme_classic() + 
  theme(axis.text.x=element_text(angle = 45, hjust = 1, size = 20), 
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20), 
        axis.text.y = element_text(size = 20), 
        axis.title.y = element_text(size = 20)) +
  labs(y = "Estimert populasjonsstørrelse") + 
  ggtitle("Teneskjær")

ggsave("Teneskjaer_2024_static.png", width = 5.4, height = 9.6, units = c("in"), dpi = 300)

p.TEN_anim <- p.TEN + transition_reveal(Year)

#create avi file
animate(p.TEN_anim, duration = 5, fps = 25, width = 540, height = 960, renderer = av_renderer())
anim_save("Teneskjaer_2024_animated.avi")

#create gif file
animate(p.TEN_anim, duration = 5, fps = 25, width = 540, height = 960, renderer = gifski_renderer())
anim_save("Teneskjaer_2024_animated.gif")

#If wanted, we can logtransfer the the y-axis
p.TEN + scale_y_continuous(trans='log10')

## Filletassen 
p.FIL <- df.pop %>% 
  dplyr::filter(Pop.variable == "size") %>% 
  dplyr::filter(Locality == "Filletassen") %>% 
  ggplot(., aes(x = Year, y = Value.est, ymin = Value.lower, ymax = Value.upper)) + 
  geom_ribbon(alpha = 0.2, fill = "lightblue") + 
  geom_point(color = "darkblue")+
  geom_line(color = "darkblue")+
  scale_x_continuous(breaks = 2014:2024, labels = 2014:2024, expand = c(0,0.1)) +
  scale_y_continuous(expand = c(0, 0))+
  theme_classic() + 
  theme(axis.text.x=element_text(angle = 45, hjust = 1, size = 20), 
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20), 
        axis.text.y = element_text(size = 20), 
        axis.title.y = element_text(size = 20)) +
  labs(y = "Estimert populasjonsstørrelse") + 
  ggtitle("Filletassen") 

ggsave("Filletassen_2024_static.png", width = 5.4, height = 9.6, units = c("in"), dpi = 300)

p.FIL_anim <- p.FIL +   transition_reveal(Year)

#create avi file
animate(p.FIL_anim, duration = 5, fps = 25, width = 540, height = 960, renderer = av_renderer())
anim_save("Filletassen_2024_animated.avi")

#create gif file
animate(p.FIL_anim, duration = 5, fps = 25, width = 540, height = 960, renderer = gifski_renderer())
anim_save("Filletassen_2024_animated.gif")

p.FIL + scale_y_continuous(trans='log10')

## Skipstadsand
p.SKI <- df.pop %>% 
  dplyr::filter(Pop.variable == "size") %>% 
  dplyr::filter(Locality == "Skipstadsand") %>% 
  ggplot(., aes(x = Year, y = Value.est, ymin = Value.lower, ymax = Value.upper)) + 
  geom_ribbon(alpha = 0.2, fill = "lightblue") + 
  geom_point(color = "darkblue")+
  geom_line(color = "darkblue")+
  scale_x_continuous(breaks = 2014:2024, labels = 2014:2024, expand = c(0,0.1)) +
  scale_y_continuous(expand = c(0, 0))+
  theme_classic() + 
  theme(axis.text.x=element_text(angle = 45, hjust = 1, size = 20), 
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20), 
        axis.text.y = element_text(size = 20), 
        axis.title.y = element_text(size = 20)) +
  labs(y = "Estimert populasjonsstørrelse") + 
  ggtitle("Skipstadsand") 

ggsave("Skipstadsand_2024_static.png", width = 5.4, height = 9.6, units = c("in"), dpi = 300)

p.SKI_anim <- p.SKI + transition_reveal(Year)

#create avi file
animate(p.SKI_anim, duration = 5, fps = 25, width = 540, height = 960, renderer = av_renderer())
anim_save("Skipstandsand_2024_animated.avi")

#create gif file
animate(p.SKI_anim, duration = 5, fps = 25, width = 540, height = 960, renderer = gifski_renderer())
anim_save("Skipstadsand_2024_animated.gif")

p.SKI + scale_y_continuous(trans='log10')

## Skjellvik
p.SKJ <- df.pop %>% 
  dplyr::filter(Pop.variable == "size") %>% 
  dplyr::filter(Locality == "Skjellvik") %>% 
  ggplot(., aes(x = Year, y = Value.est, ymin = Value.lower, ymax = Value.upper)) + 
  geom_ribbon(alpha = 0.2, fill = "lightblue") + 
  geom_point(color = "darkblue")+
  geom_line(color = "darkblue")+
  scale_x_continuous(breaks = 2014:2024, labels = 2014:2024, expand = c(0,0.1)) +
  scale_y_continuous(expand = c(0, 0))+
  theme_classic() + 
  theme(axis.text.x=element_text(angle = 45, hjust = 1, size = 20), 
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 20), 
        axis.text.y = element_text(size = 20), 
        axis.title.y = element_text(size = 20)) +
  labs(y = "Estimert populasjonsstørrelse") + 
  ggtitle("Skjellvik")

ggsave("Skjellvik_2024_static.png", width = 5.4, height = 9.6, units = c("in"), dpi = 300)

p.SKJ_anim <- p.SKJ + transition_reveal(Year)

#create avi file
animate(p.SKJ_anim, duration = 5, fps = 25, width = 540, height = 960, renderer = av_renderer())
anim_save("Skjellvik_2024_animated.avi")

#create gif file
animate(p.SKJ_anim, duration = 5, fps = 25, width = 540, height = 960, renderer = gifski_renderer())
anim_save("Skjellvik_2024_animated.gif")

p.SKJ + scale_y_continuous(trans='log10')


