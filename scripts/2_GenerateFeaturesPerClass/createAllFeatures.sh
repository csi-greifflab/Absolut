#!/bin/bash
let "k=0"
for i in *Analyses/
do
        #let "k=k+1"
		antigenID="${i%Analyses/*}"
        echo "Processing ${antigenID}"
        echo "nr  ${k}"
        #Stupid linux needs a space after if
        #if [[ $k == 20 ]]; then
            #let "k = 0"
		./AbsolutNoLib getFeatures "${antigenID}" "${antigenID}Analyses/${antigenID}_Looser.txt"	        "${antigenID}Analyses/${antigenID}_LooserFeatures.txt"	1	TRUE &
		./AbsolutNoLib getFeatures "${antigenID}" "${antigenID}Analyses/${antigenID}_LooserExclusive.txt"	"${antigenID}Analyses/${antigenID}_LooserExclusiveFeatures.txt"	1	TRUE &
		./AbsolutNoLib getFeatures "${antigenID}" "${antigenID}Analyses/${antigenID}_Mascotte.txt"	        "${antigenID}Analyses/${antigenID}_MascotteFeatures.txt"	1	TRUE &
		./AbsolutNoLib getFeatures "${antigenID}" "${antigenID}Analyses/${antigenID}_MascotteExclusive.txt"	"${antigenID}Analyses/${antigenID}_MascotteExclusiveFeatures.txt"	1	TRUE &
		./AbsolutNoLib getFeatures "${antigenID}" "${antigenID}Analyses/${antigenID}_Heroes.txt"	        "${antigenID}Analyses/${antigenID}_HeroesFeatures.txt"	1	TRUE &
		./AbsolutNoLib getFeatures "${antigenID}" "${antigenID}Analyses/${antigenID}_HeroesExclusive.txt"	"${antigenID}Analyses/${antigenID}_HeroesExclusiveFeatures.txt"	1	TRUE &
		./AbsolutNoLib getFeatures "${antigenID}" "${antigenID}Analyses/${antigenID}_superHeroes.txt"	    "${antigenID}Analyses/${antigenID}_superHeroesFeatures.txt"	1	TRUE &
		./AbsolutNoLib getFeatures "${antigenID}" "${antigenID}Analyses/${antigenID}_500kNonMascotte.txt"	"${antigenID}Analyses/${antigenID}_500kNonMascotteFeatures.txt"	1	TRUE 

		./AbsolutNoLib getFeatures "${antigenID}" "${antigenID}Analyses/${antigenID}_LooserSlices.txt"	        "${antigenID}Analyses/${antigenID}_LooserSlicesFeatures.txt"	1	TRUE &
		./AbsolutNoLib getFeatures "${antigenID}" "${antigenID}Analyses/${antigenID}_LooserExclusiveSlices.txt"	"${antigenID}Analyses/${antigenID}_LooserExclusiveSlicesFeatures.txt"	1	TRUE &
		./AbsolutNoLib getFeatures "${antigenID}" "${antigenID}Analyses/${antigenID}_MascotteSlices.txt"	        "${antigenID}Analyses/${antigenID}_MascotteSlicesFeatures.txt"	1	TRUE &
		./AbsolutNoLib getFeatures "${antigenID}" "${antigenID}Analyses/${antigenID}_MascotteExclusiveSlices.txt"	"${antigenID}Analyses/${antigenID}_MascotteExclusiveSlicesFeatures.txt"	1	TRUE &
		./AbsolutNoLib getFeatures "${antigenID}" "${antigenID}Analyses/${antigenID}_HeroesSlices.txt"	        "${antigenID}Analyses/${antigenID}_HeroesSlicesFeatures.txt"	1	TRUE &
		./AbsolutNoLib getFeatures "${antigenID}" "${antigenID}Analyses/${antigenID}_HeroesExclusiveSlices.txt"	"${antigenID}Analyses/${antigenID}_HeroesExclusiveSlicesFeatures.txt"	1	TRUE &
		./AbsolutNoLib getFeatures "${antigenID}" "${antigenID}Analyses/${antigenID}_SuperHeroesSlices.txt"	    "${antigenID}Analyses/${antigenID}_SuperHeroesSlicesFeatures.txt"	1	TRUE &

		./AbsolutNoLib getFeatures "${antigenID}" "${antigenID}Analyses/${antigenID}_500kNonMascotteSlices.txt"	"${antigenID}Analyses/${antigenID}_500kNonMascotteSlicesFeatures.txt"	1	TRUE 
 
        #else
        #fi
done
