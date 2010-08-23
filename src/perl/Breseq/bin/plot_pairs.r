min_x<-1
max_x<-4700000
min_y<-1
max_y<-4700000
infile<-"pairs.tab"
in_res<-100
use_res<-100;


for (e in commandArgs()) {
  ta = strsplit(e,"=",fixed=TRUE)
  if(! is.na(ta[[1]][2])) {
    temp = ta[[1]][2]
    if(substr(ta[[1]][1],nchar(ta[[1]][1]),nchar(ta[[1]][1])) == "I") {
      temp = as.integer(temp)
    }
    if(substr(ta[[1]][1],nchar(ta[[1]][1]),nchar(ta[[1]][1])) == "N") {
      temp = as.numeric(temp)
    }
    assign(ta[[1]][1],temp)
    cat("assigned ",ta[[1]][1]," the value of |",temp,"|\n")
  } else {
    assign(ta[[1]][1],TRUE)
    cat("assigned ",ta[[1]][1]," the value of TRUE\n")
  }
}

min_x<-as.numeric(min_x)
max_x<-as.numeric(max_x)
min_y<-as.numeric(min_y)
max_y<-as.numeric(max_y)


print(min_x)
print(max_x)
print(min_y)
print(max_y)

in_res<-as.numeric(in_res)
use_res<-as.numeric(use_res)
print(in_res)
print(use_res)

raw_data<-read.table(infile)


if (use_res != in_res)
{
	
	raw_data<-as.matrix(raw_data)

	#correct to actual coordinates!
	raw_data[,1] = raw_data[,1] * in_res;
	raw_data[,2] = raw_data[,2] * in_res;
	
	X<- array(0, dim=c(4, trunc(max_x/use_res), trunc(max_y/use_res)))
	
	
	for(i in 1:dim(raw_data)[1])
	{
		xslot <- trunc(raw_data[i,1] / use_res)+1;
		yslot <- trunc(raw_data[i,2] / use_res)+1;
		strand_slot <- raw_data[i,3]+1;
	
		X[strand_slot,xslot,yslot] = X[strand_slot,xslot,yslot] + raw_data[i,4];
	}

	max_count<-max(X)
} else {

	raw_data<-as.matrix(raw_data)
	#correct to actual coordinates!
	raw_data[,1] = raw_data[,1] * in_res;
	raw_data[,2] = raw_data[,2] * in_res;
	max_count<-max(raw_data[,4])
}

my_colors<-c('yellow', 'gray', 'orange', 'purple', 'blue', 'red', 'black');
chunk = trunc(max_count / length(my_colors));
chunk = 100;
offset = 100;



pdf("pair.pdf")

my_titles = c("Not Rev, Not Rev", "Rev, Not Rev", "Not Rev, Rev", "Rev, Rev");

for (i in 1:4)
{	

	if (use_res != in_res)
	{
		Z<-X[i,,]
		A<-c()
		for(j in 1:dim(Z)[1])
		{
			for(k in 1:dim(Z)[2])
			{
	
				if (Z[j,k] > 0)
				{
					A<-rbind(A, c(use_res*j, use_res*k, Z[j,k])) 
				}
			}
		}
		A<-data.frame(A);
	}
	else
	{
		A<-data.frame(raw_data);
		names(A)<-c("X1", "X2", "X3", "X4")
		A<-subset(A, X3 == i-1)
		
		A<-as.matrix(A);
		A<-data.frame(cbind(A[,1], A[,2], A[,4]))
		names(A)<-c("X1", "X2", "X3")
				
		#print (A)
				
	}
	
	my_title<-my_titles[i];
	
	for (m in 1:length(my_colors))
	{
		print (offset+chunk*(m-1))
		print (offset+chunk*m)
		
		Y<-subset(A, X3>=offset+chunk*(m-1))
		
		if (m < length(my_colors))
		{
			Y<-subset(Y, X3<offset+chunk*m)
		}
		
		Y<-subset(Y, abs(X1-X2)>2000)
		
		Y<-subset(Y, X1>=min_x)
		Y<-subset(Y, X1<=max_x)
		Y<-subset(Y, X2>=min_y)
		Y<-subset(Y, X2<=max_y)
				
		if (m==1)
		{
			plot(Y$X1, Y$X2, "p", cex=0.5, col=my_colors[m], xlim = c(min_x,max_x), ylim = c(min_y,max_y), xlab="Genome Coordinate", ylab="Genome Coordinate", main=my_title)
		}
		else
		{
			points(Y$X1, Y$X2, "p", cex=0.5, col=my_colors[m])
		}
	}
}

dev.off()