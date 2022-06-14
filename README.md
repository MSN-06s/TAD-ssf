# TAD-ssf
Scripts for calculating conserved/fused/seperated TADs among 2 samples

This is how this script works:
1. Add a column with "1" in the boudary files (column "bd")
2. Join the boudaries file with bins distribution file with the the "bd" column renamed as "bd.s1" and "bd.s2"
3. Fill the empty position in bd.s1 and bd.s2 with "0"
4. Add a new column with bd.s1-bd.s2 as the column "bd.comp" in the joined table
5. Add a new column with bd.s1*bd.s2 as the column "mul" in the joined table
6. For every chromosome, locate each line with mul=1(correlated boudary among 2 samples), and record their row number. Then extract every region between every correlated boundaries. 
7. For every region between correlated boundaries,sum the bd.comp column(recorded as "sum"), then extract the rownumber of bd.comp =1(recored as i.line.num) and =-1(recorded as i.m.line.num), calculate minimum gap between elements in i.line.num and i.line.num(recorded as min.gap)
8. if (sum<0){
      if (1 %in% test.mat[,6]){
        if (min.gap>=gap.threshold){
          fushion.shift<-fushion.shift+dis
        }
        else {
          fushion<-fushion+dis
        }
      }
      else{ 
        fushion<-fushion+dis
      }
    }
    else if (sum>0){
      if (-1 %in% test.mat[,6]){
        if (min.gap>=gap.threshold){
          seperation.shift<-seperation.shift+dis
        }
        else{
          seperation<-seperation+dis
        }
      }
      else{
        seperation<-seperation+dis
      }
    }
    else if (sum==0){
      if (-1 %in% test.mat[,6]){
        if (min.gap>=gap.threshold){
          shift<-shift+dis
        }
      }
    }
  }
}

![Shema of TADcompare](https://user-images.githubusercontent.com/53066081/173478571-0927fc3a-edf9-4018-af77-b4fc9e854909.png)
