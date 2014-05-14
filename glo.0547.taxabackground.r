#!/usr/bin/env Rscript

#perl -aln -F"\t" -E  'say $F[4] unless ($. == 1..4)' <(head -n10000 ~/GLO/out/glo.0546.input) | perl -nE 'BEGIN{%taxa}; @taxa = split(/_/); foreach (@taxa) { $taxa{$_}++ unless exists $taxa{$_}}; END{say scalar keys %taxa}' 

#Need to check the previous analysis's GLO are not GLOs but strings from base to top
library(ggplot2)
library(data.table)
library(igraph)
library(dplyr)
library(reshape2)
eval.string.dplyr = function(.data, .fun.name, ...) {
  args = list(...)
  args = unlist(args)
  code = paste0(.fun.name,"(.data,", paste0(args, collapse=","), ")")
  df = eval(parse(text=code,srcfile=NULL))
  df  
}

s_group_by = function(.data, ...) {
  eval.string.dplyr(.data,"group_by", ...)
}
  
#Importing
fullkos=list.files(path="~/GLO/out/glo.0546", pattern="K\\d{5}$")


lapply(fullkos, function(kos){
####################################################################################################

data=setNames(fread(sprintf("~/GLO/out/glo.0546/%s",kos)), c("CONTIG","GI","BASETAXA","BITSCORE","TREE"))
nodes=setNames(fread(sprintf("~/GLO/out/glo.0546/%s_nodes",kos)), c("KO","taxid"))
data2=data[complete.cases(data),]#Remove nulls ie. GI not linked to any TAXID; a Database problem

#Creating the graph

#Component1: EDGELIST
edgelist=unique(
rbindlist(lapply(unique(data2$TREE), function(tree) {
	splitted = strsplit(tree,"_")[[1]]
	if(length(splitted) >1){
	df = setNames(
	as.data.frame(rbind(
	do.call(rbind,lapply(1:(length(splitted)-1), function(i){
		matrix(c(splitted[i], splitted[i+1]), ncol=2)
		})),
	matrix(c(last(splitted), 1), ncol=2)
	)), 
	c("start","end"))
	df
cbind(df, data.frame(rank=c("base","genus","family","order","class","phylum")))
}
else{
	cbind(as.data.frame(matrix(c(last(splitted), 1), ncol=2)), 
	data.frame(rank=c("base"))
	)
}
}
)))

#Component2: Vertices
vertices = rbind(data.frame(node=edgelist$start, rank=edgelist$rank), data.frame(node="1", rank="root"))
#Build graph
g=graph.data.frame(d=edgelist,directed=F,vertices=vertices)
layoutt=layout.reingold.tilford(g, root=which(V(g)$name==1))	#6 is the highest lvl ie. root

#Print 
df(sprintf("out/glo.0547/%s_graph.pdf",kos))
par(srt=90)
    plot(g,
    vertex.size=0.5,
#    vertex.label=labels,
    vertex.label.cex=0.25,
    vertex.label.degree=pi,
    layout=layoutt
    )
dev.off()

#Part2:Heatmap
root=which(V(g)$name==1)
noi= V(g)[outnei(root, mode="out")]

#Weirdos
weirdos=V(g)$name[which(V(g)$name %in% noi[noi$rank== 'base']$name)]

#g=delete.vertices(g, which(V(g)$name %in% weirdos))
rankdf=data.frame(ranklvl=0:5,rankname=c("base","genus","family","order","class","phylum"))

tt=setNames(as.data.frame(
t(apply(data2, 1, function(x) {
    members=unlist(strsplit(x[names(x) == 'TREE'],"_"))
if(length(members) > 1){
	matrix(unlist(strsplit(x[names(x) == 'TREE'],"_")), nrow=1,ncol=6)
}else{
	matrix(c(members, rep("",5)), nrow=1,ncol=6)
}
})
)), c("base","genus","family","order","class","phylum"))
data3=cbind(data2,tt)


#6 plots across the different ranks

lapply(1:6, function(idk) {
level = rankdf[idk,1]
levelname=rankdf[idk,2]

sublayout= layoutt

#Order
orderlist=data.frame(
ID=V(g)$name[which(layoutt[,2]==level)],	#the layout retains weirdos
rank=V(g)$rank[which(layoutt[,2]==level)],	#the layout retains weirdos
loc=layoutt[which(layoutt[,2]==level),1]
)
orderlist=orderlist[!orderlist$ID %in% weirdos,]
#Adds back in the weirdos
orderlist=c(as.character(orderlist[order(orderlist$loc),"ID"]), weirdos)

sumdata= data3%.% 
group_by(CONTIG) %.%
s_group_by(levelname) %.% 
summarise(num=n())
sumdata = as.data.frame(sumdata)

wahliao=as.data.frame(
do.call(rbind,lapply(unique(sumdata$CONTIG), function(contig){
contigmatrix=matrix(rep(0,length(orderlist)),nrow=1)
colnames(contigmatrix)=orderlist
rownames(contigmatrix)=contig
bigcontig = filter(sumdata, CONTIG == contig)
apply(bigcontig,1, function(yy) { 
contigmatrix[,which(colnames(contigmatrix) == as.character(as.integer(yy[names(yy) == levelname])))]<<- as.integer(yy[names(yy)=='num'])
})
as.data.frame(contigmatrix)
})
))
wahliao$CONTIG = rownames(wahliao)
wh=setNames(melt(wahliao), c("CONTIG","base","value"))
wh$base = factor(wh$base, levels=as.character(orderlist))

wh$value[which(wh$value > 4)]<-4
wh$value=as.factor(wh$value)




geomdf=setNames(data.frame(t(apply(get.edgelist(g), 1, function(x) { 
c(layoutt[which(V(g)$name == x[[1]]),],layoutt[which(V(g)$name == x[[2]]),])
}))), 
c("xstart","ystart","xend","yend"))
ggplot(geomdf,aes(x=xstart,y=ystart, xend=xend, yend=yend)) + 
geom_segment()+
new_theme_empty


heatmap_plot=
	ggplot(wh, aes(x=base,y=CONTIG)) + 
	geom_tile(aes(fill=value))+
	scale_fill_manual(name="Occurence",values=c("#FFFFFF","#fef0d9","#fdcc8a","#fc8d59","#d7301f"),labels=c("0","1","2","3",">=4"))+
	theme(axis.text.x=element_text(angle=90), axis.text.y=element_text(size=0.2))+
	labs(x="TaxonID",y="Contigs")+
	ggtitle(as.character(levelname))

pdf(sprintf("out/glo.0547/%s_%s_%s_heatmap.pdf",level,levelname,kos),w=20,h=20)
print(
heatmap_plot
     )
dev.off()
})
#Add in the bit-scores
