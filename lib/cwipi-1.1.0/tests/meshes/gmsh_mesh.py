#!/usr/bin/python

import sys

lenarg = len(sys.argv)
argu = sys.argv

modelName=argu[2]
modelFile = open(modelName+"Model.geo","r")
modelString = modelFile.read()
modelFile.close()

R1=1.0
if(modelName.find("2")>-1):
  quad = 1
  h=0.2
else:
  quad = 0
  h=0.3

XCenter =0
YCenter =0 
for i in range(0,int(argu[1])):
  XCenter = XCenter + 4*R1; 
  
  if(i % 5 == 0):
    YCenter = YCenter + 4*R1
    XCenter = 0
  newFileString=modelString.replace("%h",str(h)).replace("%R1",str(R1)).replace("%quad",str(quad)).replace("%XCenter",str(XCenter)).replace("%YCenter",str(YCenter))
  newFile=open(modelName+"_part"+str(i)+".geo","w")
  newFile.write(newFileString)
  newFile.close()
