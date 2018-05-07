Code relating to analysis looking fro change in dominance structure using teh BioTIME Database. 
Instead we found evidence fo regulation in teh dominance structure. 

Only a subset of the BioTIME database was used. studies with 10 years or more of data, count or 
denisty data, and no unsuitable experimental datasets 

NullModel_turnover_inc0s.r - teh code where i run a null model randomly reordering species abundance. 
This code needs to be run on teh server, not a normal computer.

ActualDominance_correctsubset.r - the code running the mixed model looking for evidence fo global change 
in absolute (numerical) dominance. 

PercentDominance_correctsubset.r -  the code running the mixed model looking for evidence fo global change 
in relative (percent) dominance. 

AssSize_correctsubset.r - the code running the mixed model looking for evidence fo global change 
in assembaleg size (numerical abundance of all individuals in teh assemblge. 

AssSize_DomChange.r - looking to see if assemblage size change is related to dominance change. 
This is where the evidence for gerulation is strong 

compare_percentiles.r - comparing empirical data to teh null model by focusing on the number
of studies whos sloeps fall without the null model expectations

compareZScores_csubset.r _ comparing empirical data to teh null model by focusing on the z 
scored of the dominance data/null model results

plotlines.r - ploting the main figure of teh result, which is a pannel of two lien plots showing 
rates of change in each study and one overall slope for the model
