setwd('D:/qbx/resistence/salt/')
library(ftrCOOL)



## ===== Step 1: AAKC =====
filePrs<-system.file("extdata/salt/pos.fasta",package="ftrCOOL")
mat<-AAKpartComposition(seqs=filePrs,k=5,normalized=FALSE)
write.csv(mat, file = "feature/AAKC_pos.csv", row.names = TRUE)


filePrs<-system.file("extdata/salt/neg.fasta",package="ftrCOOL")
mat<-AAKpartComposition(seqs=filePrs,k=5,normalized=FALSE)
write.csv(mat, file = "feature/AAKC_neg.csv", row.names = TRUE)





### ===== Step 2: AAutoCor =====
filePrs<-system.file("extdata/salt/pos.fasta",package="ftrCOOL")
mat1<-AAutoCor(seqs=filePrs,maxlag=20,threshold=0.9,
               type=c("Moran","Geary","NormalizeMBorto","AC"))

mat2<-AAutoCor(seqs=filePrs,maxlag=20,threshold=0.9,selectedAAidx=
                 list(c('CIDH920105','BHAR880101','CHAM820101','CHAM820102'),c('CHOC760101','BIGC670101')
                      ,c('CHAM810101','DAYM780201')),type=c("AC","CC","ACC"))

write.csv(mat1, file = "feature/AAKC_autoCor1_pos.csv", row.names = TRUE)
write.csv(mat2, file = "feature/AAKC_autoCor2_pos.csv", row.names = TRUE)


filePrs<-system.file("extdata/salt/neg.fasta",package="ftrCOOL")
mat1<-AAutoCor(seqs=filePrs,maxlag=20,threshold=0.9,
               type=c("Moran","Geary","NormalizeMBorto","AC"))

mat2<-AAutoCor(seqs=filePrs,maxlag=20,threshold=0.9,selectedAAidx=
                 list(c('CIDH920105','BHAR880101','CHAM820101','CHAM820102'),c('CHOC760101','BIGC670101')
                      ,c('CHAM810101','DAYM780201')),type=c("AC","CC","ACC"))

write.csv(mat1, file = "feature/AAKC_autoCor1_neg.csv", row.names = TRUE)
write.csv(mat2, file = "feature/AAKC_autoCor2_neg.csv", row.names = TRUE)





### ===== Step 3: APAAC =====
filePrs<-system.file("extdata/salt/pos.fasta",package="ftrCOOL")
mat<-APAAC(seqs=filePrs,l=2,lambda=3,threshold=1)
write.csv(mat, file = "feature/APAAC_pos.csv", row.names = TRUE)


filePrs<-system.file("extdata/salt/neg.fasta",package="ftrCOOL")
mat<-APAAC(seqs=filePrs,l=2,lambda=3,threshold=1)
write.csv(mat, file = "feature/APAAC_neg.csv", row.names = TRUE)





### ===== Step 4: ASDC =====
filePrs<-system.file("extdata/salt/pos.fasta",package="ftrCOOL")
mat<-ASDC(seqs=filePrs)
write.csv(mat, file = "feature/ASDC_pos.csv", row.names = TRUE)


filePrs<-system.file("extdata/salt/neg.fasta",package="ftrCOOL")
mat<-ASDC(seqs=filePrs)
write.csv(mat, file = "feature/ASDC_neg.csv", row.names = TRUE)





## ===== Step 5: CkSGAApair =====
filePrs<-system.file("extdata/salt/pos.fasta",package="ftrCOOL")
mat1<-CkSGAApair(seqs=filePrs,rng=2,upto=TRUE,Grp="aromatic")

mat2<-CkSGAApair(seqs=filePrs,rng=c(1,3,5),upto=FALSE,Grp=
                   list(Grp1=c("G","A","V","L","M","I","F","Y","W"),Grp2=c("K","R","H","D","E")
                        ,Grp3=c("S","T","C","P","N","Q")))
write.csv(mat1, file = "feature/CkSGAApair1_pos.csv", row.names = TRUE)
write.csv(mat2, file = "feature/CkSGAApair2_pos.csv", row.names = TRUE)


filePrs<-system.file("extdata/salt/neg.fasta",package="ftrCOOL")
mat1<-CkSGAApair(seqs=filePrs,rng=2,upto=TRUE,Grp="aromatic")

mat2<-CkSGAApair(seqs=filePrs,rng=c(1,3,5),upto=FALSE,Grp=
                   list(Grp1=c("G","A","V","L","M","I","F","Y","W"),Grp2=c("K","R","H","D","E")
                        ,Grp3=c("S","T","C","P","N","Q")))
write.csv(mat1, file = "feature/CkSGAApair1_neg.csv", row.names = TRUE)
write.csv(mat2, file = "feature/CkSGAApair2_neg.csv", row.names = TRUE)





## ===== Step 6: conjointTriad =====
filePrs<-system.file("extdata/salt/pos.fasta",package="ftrCOOL")
mat1<-conjointTriad(seqs=filePrs)
write.csv(mat1, file = "feature/conjointTriad_pos.csv", row.names = TRUE)


filePrs<-system.file("extdata/salt/neg.fasta",package="ftrCOOL")
mat1<-conjointTriad(seqs=filePrs)
write.csv(mat1, file = "feature/conjointTriad_neg.csv", row.names = TRUE)






## ===== Step 7: conjointTriadKS =====
filePrs<-system.file("extdata/salt/pos.fasta",package="ftrCOOL")
mat1<-conjointTriad(seqs=filePrs)
write.csv(mat1, file = "feature/conjointTriadKS_pos.csv", row.names = TRUE)


filePrs<-system.file("extdata/salt/neg.fasta",package="ftrCOOL")
mat1<-conjointTriad(seqs=filePrs)
write.csv(mat1, file = "feature/conjointTriadKS_neg.csv", row.names = TRUE)






## ===== Step 8: CTD =====
filePrs<-system.file("extdata/salt/pos.fasta",package="ftrCOOL")
CTDtotal<-CTD(seqs=filePrs,normalized=FALSE)
write.csv(CTDtotal, file = "feature/CTD_pos.csv", row.names = TRUE)


filePrs<-system.file("extdata/proteins.fasta",package="ftrCOOL")
CTDtotal<-CTD(seqs=filePrs,normalized=FALSE)
write.csv(CTDtotal, file = "feature/CTD_neg.csv", row.names = TRUE)






# ===== Step 9: CTDC =====
filePrs<-system.file("extdata/salt/pos.fasta",package="ftrCOOL")
mat<-CTDC(seqs=filePrs,normalized=FALSE,label=c())
write.csv(mat, file = "feature/CTDC_pos.csv", row.names = TRUE)


filePrs<-system.file("extdata/salt/neg.fasta",package="ftrCOOL")
mat<-CTDC(seqs=filePrs,normalized=FALSE,label=c())
write.csv(mat, file = "feature/CTDC_neg.csv", row.names = TRUE)






# ===== Step 10: CTDD =====
filePrs<-system.file("extdata/salt/pos.fasta",package="ftrCOOL")
mat<-CTDD(seqs=filePrs)
write.csv(mat, file = "feature/CTDD_pos.csv", row.names = TRUE)


filePrs<-system.file("extdata/salt/neg.fasta",package="ftrCOOL")
mat<-CTDD(seqs=filePrs)
write.csv(mat, file = "feature/CTDD_neg.csv", row.names = TRUE)






# ===== Step 11: CTDT =====
filePrs<-system.file("extdata/salt/pos.fasta",package="ftrCOOL")
CTD_T<-CTDT(seqs=filePrs,normalized=FALSE)
write.csv(mat, file = "feature/CTDT_pos.csv", row.names = TRUE)


filePrs<-system.file("extdata/salt/neg.fasta",package="ftrCOOL")
CTD_T<-CTDT(seqs=filePrs,normalized=FALSE)
write.csv(mat, file = "feature/CTDT_neg.csv", row.names = TRUE)





# ===== Step 12: DDE =====
filePrs<-system.file("extdata/salt/pos.fasta",package="ftrCOOL")
mat<-DDE(seqs=filePrs)
write.csv(mat, file = "feature/DDE_pos.csv", row.names = TRUE)


filePrs<-system.file("extdata/salt/neg.fasta",package="ftrCOOL")
mat<-DDE(seqs=filePrs)
write.csv(mat, file = "feature/DDE_neg.csv", row.names = TRUE)






# ===== Step 13: EVAA =====
filePrs<-system.file("extdata/salt/pos.fasta",package="ftrCOOL")
mat<-ExpectedValueAA(seqs=filePrs,k=2,normalized=FALSE)
write.csv(mat, file = "feature/EVAA_pos.csv", row.names = TRUE)


filePrs<-system.file("extdata/salt/neg.fasta",package="ftrCOOL")
mat<-ExpectedValueAA(seqs=filePrs,k=2,normalized=FALSE)
write.csv(mat, file = "feature/EVAA_neg.csv", row.names = TRUE)






# ===== Step 14: EVGAA =====
filePrs<-system.file("extdata/salt/pos.fasta",package="ftrCOOL")
mat1<-ExpectedValueGAA(seqs=filePrs,k=2,Grp="locFus")

mat2<-ExpectedValueGAA(seqs=filePrs,k=1,Grp=
                         list(Grp1=c("G","A","V","L","M","I","F","Y","W"),Grp2=c("K","R","H","D","E")
                              ,Grp3=c("S","T","C","P","N","Q")))
write.csv(mat1, file = "feature/EVGAA1_pos.csv", row.names = TRUE)
write.csv(mat2, file = "feature/EVGAA2_pos.csv", row.names = TRUE)


filePrs<-system.file("extdata/salt/neg.fasta",package="ftrCOOL")
mat1<-ExpectedValueGAA(seqs=filePrs,k=2,Grp="locFus")

mat2<-ExpectedValueGAA(seqs=filePrs,k=1,Grp=
                         list(Grp1=c("G","A","V","L","M","I","F","Y","W"),Grp2=c("K","R","H","D","E")
                              ,Grp3=c("S","T","C","P","N","Q")))
write.csv(mat1, file = "feature/EVGAA1_neg.csv", row.names = TRUE)
write.csv(mat2, file = "feature/EVGAA2_neg.csv", row.names = TRUE)




# ===== Step 15: EVGKAA =====
filePrs<-system.file("extdata/salt/pos.fasta",package="ftrCOOL")
mat1<-ExpectedValueGKmerAA(seqs=filePrs,k=2,Grp="locFus")

mat2<-ExpectedValueGKmerAA(seqs=filePrs,k=1,Grp=
                             list(Grp1=c("G","A","V","L","M","I","F","Y","W"),Grp2=c("K","R","H","D","E")
                                  ,Grp3=c("S","T","C","P","N","Q")))
write.csv(mat1, file = "feature/EVGKAA1_pos.csv", row.names = TRUE)
write.csv(mat2, file = "feature/EVGKAA2_pos.csv", row.names = TRUE)


filePrs<-system.file("extdata/salt/neg.fasta",package="ftrCOOL")
mat1<-ExpectedValueGKmerAA(seqs=filePrs,k=2,Grp="locFus")

mat2<-ExpectedValueGKmerAA(seqs=filePrs,k=1,Grp=
                             list(Grp1=c("G","A","V","L","M","I","F","Y","W"),Grp2=c("K","R","H","D","E")
                                  ,Grp3=c("S","T","C","P","N","Q")))
write.csv(mat1, file = "feature/EVGKAA1_neg.csv", row.names = TRUE)
write.csv(mat2, file = "feature/EVGKAA2_neg.csv", row.names = TRUE)






# ===== Step 16: EVKAA =====
filePrs<-system.file("extdata/salt/pos.fasta",package="ftrCOOL")
mat<-ExpectedValueKmerAA(filePrs,k=2,normalized=FALSE)
write.csv(mat, file = "feature/EVKAA_pos.csv", row.names = TRUE)


filePrs<-system.file("extdata/salt/neg.fasta",package="ftrCOOL")
mat<-ExpectedValueKmerAA(filePrs,k=2,normalized=FALSE)
write.csv(mat, file = "feature/EVKAA_neg.csv", row.names = TRUE)






# ===== Step 17: GAAKC =====
filePrs<-system.file("extdata/salt/pos.fasta",package="ftrCOOL")
mat1<-GAAKpartComposition(seqs=filePrs,k=5,Grp="aromatic")

mat2<-GAAKpartComposition(seqs=filePrs,k=3,normalized=FALSE,Grp=
                            list(Grp1=c("G","A","V","L","M","I","F","Y","W"),Grp2=c("K","R","H","D","E")
                                 ,Grp3=c("S","T","C","P","N","Q")))
write.csv(mat1, file = "feature/GAAKC1_pos.csv", row.names = TRUE)
write.csv(mat2, file = "feature/GAAKC2_pos.csv", row.names = TRUE)


filePrs<-system.file("extdata/salt/neg.fasta",package="ftrCOOL")
mat1<-GAAKpartComposition(seqs=filePrs,k=5,Grp="aromatic")

mat2<-GAAKpartComposition(seqs=filePrs,k=3,normalized=FALSE,Grp=
                            list(Grp1=c("G","A","V","L","M","I","F","Y","W"),Grp2=c("K","R","H","D","E")
                                 ,Grp3=c("S","T","C","P","N","Q")))
write.csv(mat1, file = "feature/GAAKC1_neg.csv", row.names = TRUE)
write.csv(mat2, file = "feature/GAAKC2_neg.csv", row.names = TRUE)






# ===== Step 18: KGAAC =====
filePrs<-system.file("extdata/salt/pos.fasta",package="ftrCOOL")
mat1<-CkSGAApair(seqs=filePrs,rng=2,upto=TRUE,Grp="aromatic")

mat2<-CkSGAApair(seqs=filePrs,rng=c(1,3,5),Grp=
                   list(Grp1=c("G","A","V","L","M","I","F","Y","W"),Grp2=c("K","R","H","D","E")
                        ,Grp3=c("S","T","C","P","N","Q")))
write.csv(mat1, file = "feature/KGAAC1_pos.csv", row.names = TRUE)
write.csv(mat2, file = "feature/KGAAC2_pos.csv", row.names = TRUE)


filePrs<-system.file("extdata/salt/neg.fasta",package="ftrCOOL")
mat1<-CkSGAApair(seqs=filePrs,rng=2,upto=TRUE,Grp="aromatic")

mat2<-CkSGAApair(seqs=filePrs,rng=c(1,3,5),Grp=
                   list(Grp1=c("G","A","V","L","M","I","F","Y","W"),Grp2=c("K","R","H","D","E")
                        ,Grp3=c("S","T","C","P","N","Q")))
write.csv(mat1, file = "feature/KGAAC1_neg.csv", row.names = TRUE)
write.csv(mat2, file = "feature/KGAAC2_neg.csv", row.names = TRUE)






# ===== Step 19: PSEAAC =====
filePrs<-system.file("extdata/salt/pos.fasta",package="ftrCOOL")
mat<-PSEAAC(seqs=filePrs,l=2)
write.csv(mat, file = "feature/PSEAAC_pos.csv", row.names = TRUE)


filePrs<-system.file("extdata/salt/neg.fasta",package="ftrCOOL")
mat<-PSEAAC(seqs=filePrs,l=2)
write.csv(mat, file = "feature/PSEAAC_neg.csv", row.names = TRUE)







# ===== Step 20: PSEKRAAC =====
filePrs<-system.file("extdata/salt/pos.fasta",package="ftrCOOL")
mat1<-PseKRAAC_T11(seqs=filePrs,type="gap",Grp=4,GapOrLambdaValue=3,k=2)
mat2<-PseKRAAC_T11(seqs=filePrs,type="lambda",Grp=4,GapOrLambdaValue=3,k=2)
write.csv(mat1, file = "feature/PSEKRAAC1_pos.csv", row.names = TRUE)
write.csv(mat2, file = "feature/PSEKRAAC2_pos.csv", row.names = TRUE)


filePrs<-system.file("extdata/salt/neg.fasta",package="ftrCOOL")
mat1<-PseKRAAC_T11(seqs=filePrs,type="gap",Grp=4,GapOrLambdaValue=3,k=2)
mat2<-PseKRAAC_T11(seqs=filePrs,type="lambda",Grp=4,GapOrLambdaValue=3,k=2)
write.csv(mat1, file = "feature/PSEKRAAC1_neg.csv", row.names = TRUE)
write.csv(mat2, file = "feature/PSEKRAAC2_neg.csv", row.names = TRUE)






# ===== Step 21: QSOrder =====
filePrs<-system.file("extdata/salt/pos.fasta",package="ftrCOOL")
mat<-QSOrder(seqs=filePrs,nlag=25)
write.csv(mat, file = "feature/QSOrde_pos.csv", row.names = TRUE)
filePrs<-system.file("extdata/salt/neg.fasta",package="ftrCOOL")
mat<-QSOrder(seqs=filePrs,nlag=25)
write.csv(mat, file = "feature/QSOrde_neg.csv", row.names = TRUE)














