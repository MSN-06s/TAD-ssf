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
8. Identify different type of TAD variant in following manner:

![Fushion](https://user-images.githubusercontent.com/53066081/173772734-a225e944-4cd4-4499-b4f7-0d4900e30d5e.png)
Fushion: Sum<0,no "1" presents in the bd.comp column in the corresponding region or "1" presents, but min<=gap.threshold

![Seperation](https://user-images.githubusercontent.com/53066081/173773065-e986937c-b4a3-48a4-aa07-1ec4c15db83f.png)
Seperation: Sum>0,no "-1" presents in the bd.comp column in the corresponding region or "-1" presents, but min<=gap.threshold

![Shift](https://user-images.githubusercontent.com/53066081/173773253-3158c35e-18c5-4577-a99c-76f0512b28b1.png)
Shift: Sum=0,"1" Fushion: Sum<0,no "1" presents in the bd.comp column in the corresponding region or 1 presents, but min<=gap.threshold and min>=gap.threshold

![Seperation_shift](https://user-images.githubusercontent.com/53066081/173773785-df35f13f-92a9-4fd6-9bf6-a9aba7a3c321.png)
Seperation&Shift: Sum>0, "-1" presents in the bd.comp column in the corresponding region and min>=gap.threshold

![Fushion_shift](https://user-images.githubusercontent.com/53066081/173774118-5590672e-51a8-4059-981b-1d7b2f892731.png)
Fushion&Shift: Sum<0, "1" presents in the bd.comp column in the corresponding region and min>=gap.threshold
