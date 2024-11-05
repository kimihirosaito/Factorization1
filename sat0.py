##
# Minisat is required.
# Execute on LINUX (cygwin)
#

# Logics
#
def mkand(a,b,c): # c = a and b   (3 terms)
  return "%d %d %d 0\n%d %d 0\n%d %d 0\n"%(-a,-b,c,a,-c,b,-c)

def mkor(a,b,c): # c = a or b  (3 terms)
  return "%d %d %d 0\n%d %d 0\n%d %d 0\n"%(a,b,-c,-a,c,-b,c)

def mkxor(a,b,c): # c = a xor b  (4 terms)
  return "%d %d %d 0\n%d %d %d 0\n%d %d %d 0\n%d %d %d 0\n"\
       %(a,b,-c,-a,-b,-c,-a,b,c,a,-b,c)

def mknot(a,c): # c = not a  (2 terms)
  return "%d %d 0\n%d %d 0\n"%(a,c,-a,-c)

##
# p+4,p+5 : S,C
#   p -> p+6
# (20 terms)
def mkada42(a,b,c,d,p):
  s = ""
  s = s + mkand(c,d,p)
  s = s + mkand(a,b,p+1)
  s = s + mkxor(a,b,p+2)
  s = s + mkand(p+2,p,p+3)
  s = s + mkxor(p+2,p,p+4)
  s = s + mkor(p+1,p+3,p+5)
  return s

##
# aa[n]*b+cc[n] -> mkadanoutp(p0,n):
#  p->p+1+6*n
#
def mkadan(aa,b,cc,p,dbg=False):
  n = len(aa)
  s = ""
  if not dbg: s = s+"%d 0\n"%(-p)
  s = s+mkada42(p,cc[0],aa[0],b,p+1); p=p+1
  for ii in range(1,n):
    p=p+6
    s = s+mkada42(p-1,cc[ii],aa[ii],b,p);
  return s

#
# print(mkadanoutp(5,2))
def mkadanoutp(p0,n):
  return [p0+5+6*ii for ii in range(n)]+[p0+6*n]


####
def v2s(v,idxs):
  s = ""
  for ii in idxs: s=s+str(v[ii-1])
  return s

def s2i(s):
  jj=0;
  for ii in range(len(s)): jj = jj+2**ii*int(s[ii])
  return jj

#  print(i2s(10,range(1,10)))
def i2s(ii,idxs):
  si = format(ii,"b")[::-1]; #print(si)
  n = max(len(si),len(idxs))
  s = ""
  for ii in range(len(si)):
    if ii < len(si):
      s = s+"%d 0\n"%(-idxs[ii] if int(si[ii]) <= 0 else idxs[ii])
    else: s = s+"%d 0\n"%(-idxs[ii])
  return s


# input t0, output t1
#
def runms(s,f0="t"):
  open("%s0"%f0,"wt").write(s)

  os.system("minisat %s0 %s1"%(f0,f0))
#  os.system("cat t1")

  so = open("%s1"%f0,"rt").read();  #print(so)
  so = map(int,so.split()[1:-1])
  for ii in range(len(so)): so[ii] = 0 if so[ii] < 0 else 1
  return so

#
# p0 = 2*na+nb+1
# out = [p0+(1+6*na)*(ii+1)-2 for ii in range(na)]+[p0+(1+6*na)*na-1]
#
def mkmul(na,nb):
  s = ""
  for ii in range(na): s = s+"-%d 0\n"%(na+nb+ii+1)
  p = 2*na+nb+1
  s = s+mkadan(range(1,na+1),na+1,range(na+nb+1,2*na+nb+1),p);
  for ii in range(1,nb):
    cc = mkadanoutp(p,na)[1:]; #print(cc)
    p=max(cc)+1
    s = s+mkadan(range(1,na+1),na+ii+1,cc,p);
  return s

#
def muloutp(na,nb,p0=0):
  if p0 <= 0: p0 = 2*na+nb+1
  p = [p0+5+(6*na+1)*ii for ii in range(nb)]
  p = p+[p0+5+(6*na+1)*(nb-1)+6*ii for ii in range(1,nb)]
  p = p+[p[-1]+1]
  return p


###

if __name__ == "__main__":
  import os

#  na = 49; nb = 51; nv = 1267650600112341087037729009993
#  na = 50; nb = 50; nv = 1267650600126761050094036356853
  na = 33; nb = 33; nv = 73786518486371773127
#  na = 32; nb = 33; nv = 36893341178068089667
#  na = 24; nb = 24; nv = 279966626574091
#  na = 23; nb = 24; nv = 139983321660649
#  na = 21; nb = 21; nv = 4211070206687
#  na = 20; nb = 21; nv = 2135730774041
#  na = 16; nb = 16; nv = 4271011009
#  na = 16; nb = 20; nv = 68578400483

  s = ""
  s = s + mkmul(na,nb)
  outi = muloutp(na,nb); print(outi) ;print(len(outi))
  s = s + i2s(nv,outi)
  print(len(s))

  so = runms(s,"tt")

  print(len(so))
  print(v2s(so,range(1,na+1)))
  print(v2s(so,range(na+1,na+nb+1)))
  print
  print("p1= "+str(s2i(v2s(so,range(1,na+1)))))
  print("p2= "+str(s2i(v2s(so,range(na+1,na+nb+1)))))
