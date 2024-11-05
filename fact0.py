# -*- encoding: UTF-8 -*-
#
#from lg3alib0 import *

##
def v2i(v,i0,i1):
  a=0; f=1
  for ii in range(i0,i1): a = a+v[ii]*f; f = f*2
  return a
def i2v(a,nb=-1):
  v=[]
  if nb > -1:
    for ii in range(nb): v.append(a%2); a=a//2
  else:
    while a > 0: v.append(a%2); a=a//2
  return v

def v2s(v):
  ss = ""
  for b in v: ss = ss + str(b) 
  return ss

def s2v(s):
  v=[]
  for ii in range(len(s)):
    if s[ii]=="0": v.append(0)
    elif s[ii]=="1": v.append(1)
    elif s[ii]=="2": v.append(2)
    else: v.append(3)
  return v

def foutpara(fname,nb1,nb2,s0,sn=[]):
  if sn==[]: sn=[-1]
  f = open(fname,"wt")
  f.write("NB1 %d\nNB2 %d"%(nb1,nb2))
  f.write("\nS0 ");
  f.write(" ".join(map(str,s0)))
  f.write("\nSn ");
  f.write(" ".join(map(str,sn)))
  f.close()



###########################
#
if __name__ == "__main__":


  n51b1 = 2251799813649961;    n51b2 = 2251799813614873
  n50b1 = 1125899906812609;    n50b2 = 1125899906782517
  n49b1 = 562949953387441; n49b2 = 562949953352827
  n33b1 = 8589903197; n33b2 = 8589912691; # 6633/101, 61347->30291->17688, 594.6
  n32b1 = 4294961137; n32b2 = 4294951021; # 6274/98, 24900->24711->16447, 162
  n24b1 = 16747207; n24b2 = 16717213
  n23b1 = 8358607; n23b2 = 8328583
  n21b1 = 2067137; n21b2 = 2037151;
  n20a1 = 1048703; n20a2 = 1048837; n20b1 = 1048391; n20b2 = 1048219
  n16a1 = 65599; n16a2 = 65699; n16b1 = 65413; n16b2 = 65293


#  a1 = n50b1; a2 = n50b2; a = a1*a2; # 8.2
#  a1 = n49b1; a2 = n51b1; a = a1*a2; # 9.4
#  a1 = n33b1; a2 = n33b2; a = a1*a2; # 3.4
#  a1 = n32b1; a2 = n33b2; a = a1*a2; # 1.75

#  a1 = n24b1; a2 = n24b2; a = a1*a2; # 2.52
#  a1 = n23b1; a2 = n24b1; a = a1*a2; # 2.44
#  a1 = n21b1; a2 = n21b2; a = a1*a2; # 2.16
#  a1 = n20b1; a2 = n21b2; a = a1*a2; # 0.68
  a1 = n16b1; a2 = n20b1; a = a1*a2; # 0.83
#  a1 = n16b1; a2 = n16b2; a = a1*a2; # 0.84


  print("ANS:%d"%a)
  nb1 = len(i2v(a1)); nb2 = len(i2v(a2));
  print(nb1,nb2)

  s0 = [2]*(nb1+2)+i2v(a,nb1+nb2)

  import os
  import time
  import datetime

  fnamep = "fnamep.txt"
  os.system("gcc fact0.c -O")
  sn = [3]*nb1+[2]*(nb1+2)+[3]*nb2; sn[0]=1;sn[nb1-1]=1;sn[2*nb1+2]=1; sn[2*nb1+nb2+1]=1
  foutpara(fnamep,nb1,nb2,s0,sn)
  t0 = time.time();  print(datetime.datetime.now())
  os.system("a %s"%(fnamep))
  print("Time: %f"%(time.time()-t0))
  s=open("out.txt","rt").read()
  print(v2s(i2v(v2i(s2v(s),0,nb1))))
  print(v2s(i2v(v2i(s2v(s),nb1*2+2,nb1*2+2+nb2))))
  print("p1 = "+str(v2i(s2v(s),0,nb1)))
  print("p2 = "+str(v2i(s2v(s),nb1*2+2,nb1*2+2+nb2)))


