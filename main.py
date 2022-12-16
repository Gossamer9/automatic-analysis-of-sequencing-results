#import requests
#with open("C:/Users/Куська/Desktop/python/hey.txt",'r') as inf:

from Bio import AlignIO
from Bio.Seq import Seq

alignment = AlignIO.read(open("C:/Users/Куська/Desktop/C++/cpp/sample_002.fas"), "fasta")
size=len(alignment[1,:]);
gap=0;
start=1620;
end=2441;
startcodon=211;
num=[];

for i in range(size):    #счётчик gap в референсе num, где индексы пробелов +
    if str(alignment[1,i])== '-':
        num.append(i);
frame=startcodon%3;      #рамка считывания
for i in range(start,end,1):   #замена аминокислот
    if alignment[0,i]!=alignment[1,i] and str(alignment[0,i])!= '-':
        if i%3==frame and str(alignment[0,i+1])!= '-' and str(alignment[0,i+2])!= '-': #if первая в кодоне
            x = str();
            for el in alignment[0,i:i+3]:
                x += el;
            x=Seq(x);
            y = str();
            for el in alignment[1,i:i+3]:
                y += el;
            y=Seq(y);
            print('Замена АК', x.translate() , '-->', y.translate() )
        else:                  #if вторая в кодоне
            if i%3==(frame+1)%3 and str(alignment[0,i+1])!= '-' and str(alignment[0,i-1])!= '-':
                x = str();
                for el in alignment[0,i-1:i+2]:
                    x += el;
                x=Seq(x);
                y = str();
                for el in alignment[1,i-1:i+2]:
                    y += el;
                y=Seq(y);
                print('Замена АК', y.translate() , '-->', x.translate() )
            else:                  #if третья в кодоне
                if i%3==(frame+2)%3 and str(alignment[0,i-2])!= '-' and str(alignment[0,i-1])!= '-':
                    x = str();
                    for el in alignment[0,i-2:i+1]:
                        x += el;
                    x=Seq(x);
                    y = str();
                    for el in alignment[1,i-2:i+1]:
                        y += el;
                    y=Seq(y);
                    print('Замена АК', y.translate() , '-->', x.translate() )

 #   alig=Seq("TCC");
#print(alig.translate())
