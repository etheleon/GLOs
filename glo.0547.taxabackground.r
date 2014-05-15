#!/usr/bin/env Rscript

#perl -aln -F"\t" -E  'say $F[4] unless ($. == 1..4)' <(head -n10000 ~/GLO/out/glo.0546.input) | perl -nE 'BEGIN{%taxa}; @taxa = split(/_/); foreach (@taxa) { $taxa{$_}++ unless exists $taxa{$_}}; END{say scalar keys %taxa}' 

#Need to check the previous analysis's GLO are not GLOs but strings from base to top
library(gridExtra)
library(ggplot2)
library(data.table)
library(igraph)
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(RJSONIO)
library(RCurl)
args = commandArgs(T)
print(args[1])
dir.create("out/glo.0547/table")
#dir.create("out/glo.0547/occurrence")
colorlabels = brewer.pal(n=6,name="Set1")
ranklist=c("base","genus","family","order","class","phylum","root")
relate=c("base-root","base-genus","genus-family","family-order","order-class","class-phylum","phylum-root")


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
  
new_theme_empty <- theme_bw()
new_theme_empty$line <- element_blank()
new_theme_empty$rect <- element_blank()
new_theme_empty$strip.text <- element_blank()
new_theme_empty$axis.text <- element_blank()
new_theme_empty$plot.title <- element_blank()
new_theme_empty$plot.margin <- structure(c(0, 0, -1, -1), unit = "lines", valid.unit = 3L, class = "unit")

#Importing
#fullkos=list.files(path="~/GLO/out/glo.0546", pattern="K\\d{5}$")
#fullkos = 'K15512'
#Iterate through kos 
#lapply(fullkos, function(kos){
kos = args[1]
pdf(sprintf("out/glo.0547/%s_heatmap.pdf",kos),w=20,h=30)
####################################################################################################
data=setNames(fread(sprintf("~/GLO/out/glo.0546/%s",kos)), c("CONTIG","GI","BASETAXA","BITSCORE","TREE"))
nodes=setNames(fread(sprintf("~/GLO/out/glo.0546/%s_nodes",kos)), c("KO","taxid"))
data2=data[complete.cases(data),]#Remove nulls ie. GI not linked to any TAXID; a Database problem

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

#Print Full graph
#pdf(sprintf("out/glo.0547/%s_graph.pdf",kos))
par(srt=90)
    plot(g,
    vertex.size=0.5,
#    vertex.label=labels,
    vertex.label.cex=0.25,
    vertex.label.degree=pi,
    layout=layoutt
    )
#dev.off()

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
data4=data3[with(data3, order(phylum, class, order, family, genus, base)),]
write.csv(data4, file=sprintf("out/glo.0547/table/%s.table.csv",kos))

#6 plots across the different ranks "base","genus","family","order","class","phylum"
lapply(1:6, function(idk) {
level = rankdf[idk,1]
levelname=rankdf[idk,2]

if(level == 0) { 
#Move the unclassified to the bottom
layoutt[which(V(g)$name %in% weirdos),1]<-seq(from=max(layoutt[layoutt[,2]==0,1])+1, to=max(layoutt[layoutt[,2]==0,1]+1 + length(weirdos))-1, by=1)
layoutt[which(V(g)$name %in% weirdos),2]<-0
sublayoutt= layoutt
colorlabels=c("black",colorlabels)
}else{
#Delete nodes
g=delete.vertices(g, which(V(g)$rank %in% as.character(rankdf[rankdf$ranklvl %in% c(0:6)[0:6 < level],2])))
sublayoutt=layout.reingold.tilford(g, root=which(V(g)$name==1))	#6 is the highest lvl ie. root
}

#Order for the heatmap
orderlist=data.frame(
ID=V(g)$name[which(sublayoutt[,2]==0)],	#the layout retains weirdos
rank=V(g)$rank[which(sublayoutt[,2]==0)],	#the layout retains weirdos
loc=sublayoutt[which(sublayoutt[,2]==0),1]
)
orderlist=as.character(orderlist[order(orderlist$loc),"ID"])

####################################################################################################
query = 'start basetaxa=node:ncbitaxid(taxid={taxid}) return basetaxa.name'
cypherurl="http://192.168.100.1:7474/db/data/cypher"
namelist = do.call(rbind,lapply(orderlist, function(taxa) { 
post=toJSON(
	    list(
		query=query,
		params=list(taxid=taxa)
		)
	    )
result=fromJSON(
	getURL(cypherurl, 
	customrequest='POST', 
	httpheader=c('Content-Type'='application/json'), 
	postfields=post
	))
if(length(unlist(result$data)) > 0 ){
data.frame(id=taxa,name=unlist(result$data))
}
}))
####################################################################################################


#melted data.frame
sumdata= data3 %.% 
group_by(CONTIG) %.%
s_group_by(levelname) %.% 
summarise(num=n())
sumdata = as.data.frame(sumdata)

#Construct the matrix
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
contigorder = rownames(wahliao)[hclust(dist(select(wahliao, -CONTIG)))$order]

#Melt again
wh=setNames(melt(wahliao), c("CONTIG","base","value"))
wh$base = factor(wh$base, levels=as.character(orderlist))

#Change the numbering need validating
#pdf(sprintf("out/glo.0547/occurrence/%s_%s_occurence.distribution.pdf",kos,levelname),w=10,h=10)
#qplot(wh$value)+xlab("Occurence")+ylab("Frequency")+ggtitle(levelname)
#dev.off()
occurrencetable = setNames(melt(table(wh$value)), c("Occurrence","Frequency"))
occurrencetable$rank = levelname 
write.table(occurrencetable , file=sprintf("out/glo.0547/table/%s_frequency.txt",kos),append=TRUE,quote=F,row.names=F,col.names=F)

wh$value[which(wh$value > 4)]<-4
wh$value=as.factor(wh$value)
wh$base = factor(wh$base, labels=as.character(namelist[match(levels(wh$base),namelist$id),]$name))
wh$CONTIG = factor(wh$CONTIG, levels = contigorder)


#The heatmap
heatmap_plot=
	ggplot(wh, aes(x=base,y=CONTIG)) + 
	geom_tile(aes(fill=value))+
	scale_fill_manual(name="Occurence",values=c("#FFFFFF","#fef0d9","#fdcc8a","#fc8d59","#d7301f"))+
	theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0), axis.text.y=element_text(size=0.2),legend.position="none")+
	labs(x="TaxonID",y="Contigs")

#graph
geomdf=setNames(
	data.frame(t(apply(get.edgelist(g), 1, function(x) { 
		    c(sublayoutt[which(V(g)$name == x[[1]]),],sublayoutt[which(V(g)$name == x[[2]]),])
		    }))), 
	c("xstart","ystart","xend","yend"))

relation = data.frame(ystart=c(0:(6-idk+1)),yend=c(1:(6-idk+1),0))

relation$type = apply(relation, 1, function(x) { 
paste(ranklist[sort(x)[[1]]+idk], ranklist[sort(x)[[2]]+idk],sep="-")
})

binding=unname(as.matrix(relation))
binding=setNames(as.data.frame(rbind(binding,binding[,c(2,1,3)])), c("ystart","yend","type"))
geomdf=merge(geomdf, binding,all.x=T, by=c("ystart","yend"))
geomdf$type = factor(geomdf$type)
if(idk==1){
geomdf$type = factor(geomdf$type, levels = relate)
}else{
geomdf$type = factor(geomdf$type, levels = relate[(idk+1):7])
}

leaves = setNames(
data.frame(cbind(V(g)$name[which(V(g)$name %in% orderlist)], sublayoutt[which(V(g)$name %in% orderlist),])),
 c("id","x","y"))

namelist=merge(namelist,
leaves,
by="id",
all.x=T)

#The dendrogram
tree=ggplot(geomdf,aes(
	    x=xstart,
	    y=ystart, 
	    xend=xend, 
	    yend=yend
	    )) + 
geom_segment(aes(colour=type),size=2)+
geom_text(data=namelist, aes(x=as.numeric(as.character(x)),y=as.numeric(as.character(y)),label=as.character(name)),angle=90,inherit.aes=FALSE,hjust=1,vjust=0,size=3)+
scale_y_continuous(labels=ranklist[idk:(6+1)],	breaks=c(0:(6-level)),limits=c(-1,max(geomdf$yend)))+
scale_color_manual(values=colorlabels[idk:length(colorlabels)])+
new_theme_empty+
theme(legend.position="none",plot.margin=unit(c(0,-1.5,0,-1.5), "cm"))

print(
grid.arrange(tree,heatmap_plot,ncol=1)
)
})
dev.off()
#})
