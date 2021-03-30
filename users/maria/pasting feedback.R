journ<-read.csv("journ.csv", stringsAsFactors = FALSE)
journ2<-c()
for (i in 1:22){
  a<-paste("matriculation number", journ$matric[i], sep = ": ")
  a<- paste(a, "The assignment was assessed according to six criteria: pitch (appropriate to the target audience), order the information was presented, interest (whether it keeps the attention of the reader), clarity, integration of wider knowledge and readibility and style. Pitch:", journ$pitched.at.right.audience[i],sep = " ")
  a<-paste(a, "Order:", journ$order.of.information[i], "Interest:", journ$interesting[i], "Clarity:", journ$clear[i], "Wider sources: ", journ$integrates.wider.knowledge[i], "Style: ", journ$style[i] ,sep = " ")
  journ2[i]<-a
  rm(a)
}

write.csv(journ2, "journ2.csv")
