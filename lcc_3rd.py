import secondQuantizationAlgebra as sqa
from tools_for_sqa import *
from tools_for_pt import *
from old_tools_for_sqa import *
sqa.options.verbose = False







# PS1 ---------------------------------------------------------------------------------------------
#  |Psi1> =         d_mu O  _mu |Psi_0>
#  <Psi1| = <Psi_0| d_mu O^T_mu
#  where:
#   - the indexes of "O_mu" and "O^T_mu" are dictated by the class itself,
#     see that they match with the rest of the sqa files!
#       (first 4 and next 4 of "indExOp_Psi1")
#   - the indexes of "d_mu" for each class are to match the unique tensors outputed by icpt!
#       (first 4 and next 4 of "indTens_Psi1")
#   - those tensors are labelled by "names" that are the classes
#     and not the correspondinng string of indexes types!
#       (last element of "indTens_Psi1")

# List of indexes for the ExOp:
#              O_mu              O^T_nu            (name)
indExOp_Psi1=[[a1 ,i1,a8 ,i8 ,   i18,a18,i15,a15,  'vvcc'],
              [a2 ,i2,a9 ,p1 ,   p9 ,a19,i9 ,a16,  'vvca'],
              [p2 ,i3,a3 ,i10,   i19,a10,i16,p10,  'ccav'],
              [a4 ,p3,a11,p11,   p23,a20,p17,a17,  'vvaa'],
              [p4 ,i4,p12,i11,   i20,p24,i17,p18,  'ccaa'],
              [a5 ,i5,p13,p5 ,   p19,p25,i12,a12,  'vaca'],
              [p14,i6,a6 ,p6 ,   p20,a13,i13,p26,  'avca'],
              [a7 ,p7,p21,p15,   p29,p31,p27,a14,  'vaaa'],
              [p16,i7,p22,p8 ,   p28,p32,i14,p30,  'caaa']]
# List of indexes for the tensors:
#              d_mu              d_nu              (name)
indTens_Psi1=[[a1 ,a8 ,i1 ,i8 ,  a15,a18,i15,i18,  'vvcc'],
              [a2 ,a9 ,i2 ,p1 ,  a16,a19,i9 ,p9 ,  'vvca'],
              [i3 ,i10,p2 ,a3 ,  i16,i19,p10,a10,  'ccav'],
              [a4 ,a11,p3 ,p11,  a17,a20,p17,p23,  'vvaa'],
              [i4 ,i11,p4 ,p12,  i17,i20,p18,p24,  'ccaa'],
              [a5 ,p13,i5 ,p5 ,  a12,p25,i12,p19,  'vaca'],
              [p14,a6 ,i6 ,p6 ,  p26,a13,i13,p20,  'avca'],
              [a7 ,p21,p7 ,p15,  a14,p31,p27,p29,  'vaaa'],
              [i7 ,p22,p16,p8 ,  i14,p32,p30,p28,  'caaa']]
for index in range(len(indExOp_Psi1)):
  if sorted(indExOp_Psi1[index][:4])!=sorted(indTens_Psi1[index][:4]):
    print "WHAT??"
    exit()
  if sorted(indExOp_Psi1[index][4:8])!=sorted(indTens_Psi1[index][4:8]):
    print "WHAT??"
    exit()
  if indExOp_Psi1[index][-1]!=indTens_Psi1[index][-1]:
    print "WHAT??"
    exit()




# H -----------------------------------------------------------------------------------------------
#   H = t E2 + w E4
#   where for "E4":
#    - the indexes of "E4" run through all possible unique 4-tuples
#        (this is done by sorted the first two indexes)
#      UPDATE: it now runs through all possible 4-tuples
#    - some of these "E4" operators belong to particular classes
#        (this is associated via the last elt of "indExOp_H4" using the "indExOpClass" dictionary)
#   and where for "w"
#    - the indexes of "w" are permutated to refer to unique tensors
#        (this is done by using the "symmTens" list)




# E4 ----------------------------------------------------------------------------------------------
# Links an ExOp of the Hamiltonian to a Class
indExOpClass={}
indExOpClass['vvcc'] =['vvcc']
indExOpClass['vvca'] =['vvac','vvca']
indExOpClass['ccav'] =['avcc','vacc']
indExOpClass['vvaa'] =['vvaa']
indExOpClass['ccaa'] =['aacc']
indExOpClass['vaca'] =['avac','vaca']
indExOpClass['avca'] =['avca','vaac']
indExOpClass['vaaa'] =['avaa','vaaa']
indExOpClass['caaa'] =['aaac','aaca']
indExOpClass['h0']   =['cc','aa','vv',\
                       'cccc','aaaa','vvvv',\
                       'acac','acca','caca','caac',\
                       'cvcv','cvvc','vcvc','vccv',\
                       'avav','avva','vava','vaav']
indExOpClass['heff2']=[]
indExOpClass['heff4']=[]
indExOpClass['heff2']=indExOpClass['heff2']+\
                      ['vc',\
                       'va',\
                       'ac']
indExOpClass['heff4']=indExOpClass['heff4']+\
                      ['cvcc','vccc',\
                       'cvca','cvac','vcac','vcca',\
                       'accc','cacc']
indExOpClass['heff4']=indExOpClass['heff4']+\
                      ['vvcv','vvvc',\
                       'vvav','vvva',\
                       'avcv','avvc','vavc','vacv']

# Construct the full Hamiltonian ExOp and assign classes.
# (The indexes run through ABCD.)
indlist=[[iH,pH,aH],[jH,qH,bH],[kH,rH,cH],[lH,sH,dH]]
factor={}
for ind0 in range(3):
  for ind1 in range(3):
    for ind2 in range(3):
      for ind3 in range(3):
        factor[code([indlist[0][ind0],indlist[1][ind1],indlist[2][ind2],indlist[3][ind3]])]=0

print '// To make up E4(AB<-CD), there is 3*3*3*3=81 possible 4-tuples'
print '// some of them are redundant are skipped.'
indExOp_H4=[]
iter=0
for ind0 in range(3):
  for ind1 in range(3):
    for ind2 in range(3):
      for ind3 in range(3):
        list=[indlist[0][ind0],indlist[1][ind1],indlist[2][ind2],indlist[3][ind3]]
        name=code(list)

        iter=iter+1
        skip=False
        if ((name[0]==name[1] and name[-2:]!=''.join(sorted(name[-2:])))
         or (''.join(sorted(name[:2]))!=name[:2])):
            skip=True
        if not skip:
          factor[name]=factor[name]+1
          line='{:2} "{:4}" {:} {:4}'.format(iter,code(list),'included as ',name)
        else:
          factor[name[1]+name[0]+name[3]+name[2]]=factor[name[1]+name[0]+name[3]+name[2]]+1
          line='{:2} "{:4}" {:} {:4}'.format(iter,code(list),'skipped (see',name[1]+name[0]+name[3]+name[2]+')')
        print '//',line
        if skip:
          continue

        newname='other'
        for cl in ClassNames+['h0','heff2','heff4']:
          if name in indExOpClass[cl]:
            newname=cl
        #for Class,Pattern in occupations.items():
        #  if Pattern == pattern(code_rearrange(list)):
        #    newname=Class
        indExOp_H4.append(list+[newname])
print '// Sum of factors:',sum(factor.values())
print '//\n//'




# W -----------------------------------------------------------------------------------------------
#print '// To find unique corresponding tensor,'
#print '// we need to translate some indexes.'
# Uses symmetry of two-el integrals
# to use as few unique tensors as possible
symmTens=[
          ['vvcc','ccvv','vccv','cvvc'],
          ['vvca','acvv','avvc','vvac','vcav','cvva','vacv','cavv'],
          ['ccav','accv','ccva','cvac','avcc','vcca','cavc','vacc'],
          ['vvaa','aavv','vaav','avva'],
          ['ccaa','aacc','acca','caac'],
          ['vaca','cava','acav','avac','cava'],
          ['avca','aacv','aavc','cvaa','vaac','acva','vcaa','caav'],
          ['vaaa','aaav','avaa','aava'],
          ['caaa','aaac','acaa','aaca'],
          #
          ['ccca','accc','cacc','ccac'],
          ['cccv','vccc','cvcc','ccvc'],
          ['vvvc','cvvv','vcvv','vvcv'],
          ['vvva','avvv','vavv','vvav'],
          ['caca','acac'],
          ['cvcv','vcvc'],
          ['avav','vava'],
          ['cacv','cvca','vcac','acvc'],
          ['avcv','vavc','vcva','cvav'],
         ]
# DOUBLE CHECK FOR ERRORS
# check that all lists indeed give the same integral
for elt in symmTens:
  comp=elt[0]
  model1=comp[0]+comp[2]
  model2=comp[1]+comp[3]
  for iteri in range(1,len(elt)):
    comp=elt[iteri]
    comp1N=comp[0]+comp[2]
    comp1R=comp[2]+comp[0]
    comp2N=comp[1]+comp[3]
    comp2R=comp[3]+comp[1]
    if((comp1N+comp2N!=model1+model2)and
       (comp1N+comp2R!=model1+model2)and
       (comp1R+comp2N!=model1+model2)and
       (comp1R+comp2R!=model1+model2)and
       (comp2N+comp1N!=model1+model2)and
       (comp2R+comp1N!=model1+model2)and
       (comp2N+comp1R!=model1+model2)and
       (comp2R+comp1R!=model1+model2)):
      print 'ERROR symmTens',elt[0], comp, comp1N, comp2N, comp1R, comp2R
      exit()

# Construct the tensors W
indTens_H4=[]
for iteri in range(len(indExOp_H4)):
  # It's E(pqrs)=adagger_p(1) adagger_q(2) a_s(2) a_r(1) -> <pq|rs>
  tensor=[indExOp_H4[iteri][i] for i in [0,1,2,3]]
  name  =code(tensor)

  for iterj in range(len(symmTens)):
    if name in symmTens[iterj] and name!=symmTens[iterj][0]:
      comp=symmTens[iterj][0]
      model1=comp[0]+comp[2]
      model2=comp[1]+comp[3]
      comp1N=name[0]+name[2]
      comp1R=name[2]+name[0]
      comp2N=name[1]+name[3]
      comp2R=name[3]+name[1]
      if(comp1N+comp2N==model1+model2):
        trans=[0,1,2,3]
      elif(comp1N+comp2R==model1+model2):
        trans=[0,3,2,1]
      elif(comp1R+comp2N==model1+model2):
        trans=[2,1,0,3]
      elif(comp1R+comp2R==model1+model2):
        trans=[2,3,0,1]
      elif(comp2N+comp1N==model1+model2):
        trans=[1,0,3,2]
      elif(comp2R+comp1N==model1+model2):
        trans=[3,0,1,2]
      elif(comp2N+comp1R==model1+model2):
        trans=[1,2,3,0]
      elif(comp2R+comp1R==model1+model2):
        trans=[3,2,1,0]
      #print '//',name,'=',\
      #      ''.join([name[A]       +str(A) for A in range(4)]),'->',\
      #      ''.join([name[trans[A]]+str(trans[A]) for A in range(4)]),'=',\
      #      symmTens[iterj][0]
      tensor=[tensor[iterk] for iterk in trans]
      name=''.join([name[iterk] for iterk in trans])
      if name!=symmTens[iterj][0]:
        print 'ERROR symmTens',name,symmTens[iterj][0]
        exit()

  indTens_H4.append(tensor+[name])




# E2 ----------------------------------------------------------------------------------------------
# Construct the full Hamiltonian ExOp
# (The indexes run through AB)
indlist=[[iH,pH,aH],[jH,qH,bH]]
indExOp_H2=[]
for ind0 in range(3):
  for ind1 in range(3):
    list=[indlist[0][ind0],indlist[1][ind1]]
    name=code(list)

    newname='other'
    for cl in ClassNames+['h0','heff2','heff4']:
      if name in indExOpClass[cl]:
        newname=cl
    #for Class,Pattern in occupations.items():
    #  if Pattern == pattern(code_rearrange(list)):
    #    newname=Class
    indExOp_H2.append(list+[newname])




# T -----------------------------------------------------------------------------------------------
# Construct the tensors t
# (again we use the symmetries of the one-el integral
#  to minimize the number of different tensors needed.)
indTens_H2=[]
for iteri in range(len(indExOp_H2)):
  tensor=[indExOp_H2[iteri][i] for i in [0,1]]
  name  =code(tensor)
  if name=='ac' or name=='av' or name=='cv':
    tensor=[tensor[1],tensor[0]]
    name=name[1]+name[0]
  indTens_H2.append(tensor+[name])




# Print out stuff ---------------------------------------------------------------------------------
print '//'
print '// With the unique tensors found, |Psi1> is:'
for cl in ClassNames+['h0','heff2','heff4']:
  for iteri in range(len(indExOp_Psi1)):
    if indExOp_Psi1[iteri][-1]==cl:
      print '// d_{:5} E_{:2} E_{:2}        Class:{:6} OccPattern:{:9}'.format(\
                            code(indTens_Psi1[iteri][:4]),\
                            code(indExOp_Psi1[iteri][:2]),\
                            code(indExOp_Psi1[iteri][2:4]),\
                            cl,\
                            output_pattern(code(indExOp_Psi1[iteri][:4])))
print '//'
print '// With the unique tensors found, <Psi1| is:'
for cl in ClassNames+['h0','heff2','heff4']:
  for iteri in range(len(indExOp_Psi1)):
    if indExOp_Psi1[iteri][-1]==cl:
      print '// d_{:5} E_{:2} E_{:2}        Class:{:6} OccPattern:{:9}'.format(\
                            code(indTens_Psi1[iteri][4:8]),\
                            code(indExOp_Psi1[iteri][4:6]),\
                            code(indExOp_Psi1[iteri][6:8]),\
                            cl,\
                            output_pattern(code(indExOp_Psi1[iteri][4:8])))
print '//'
print '//'
print '// With the unique tensors found, the Hamiltonian is:'
for cl in ClassNames+['h0','heff2','heff4']:
  for iteri in range(len(indExOp_H4)):
    if indExOp_H4[iteri][-1]==cl:
      print '// w_{:5} E_{:7} ({:1})    Class:{:6} OccPattern:{:9}'.format(\
                            code(indTens_H4[iteri][:4]),\
                            code(indExOp_H4[iteri][:4]),\
                            factor[code(indExOp_H4[iteri][:4])],\
                            cl,\
                            output_pattern(code_rearrange(indExOp_H4[iteri][:4])))
  for iteri in range(len(indExOp_H2)):
    if indExOp_H2[iteri][-1]==cl:
      print '// t_{:5} E_{:7} ({:1})    Class:{:6} OccPattern:{:9}'.format(\
                            code(indTens_H2[iteri][:2]),\
                            code(indExOp_H2[iteri][:2]),\
                            1,\
                            cl,
                            output_pattern(code_rearrange(indExOp_H2[iteri][:2])))
for iteri in range(len(indExOp_H4)):
  if indExOp_H4[iteri][-1]=='other':
      print '// w_{:5} E_{:7} ({:1})    Class:{:6} OccPattern:{:9}'.format(\
                            code(indTens_H4[iteri][:4]),\
                            code(indExOp_H4[iteri][:4]),\
                            factor[code(indExOp_H4[iteri][:4])],\
                            'other',\
                            output_pattern(code_rearrange(indExOp_H4[iteri][:4])))
for iteri in range(len(indExOp_H2)):
  if indExOp_H2[iteri][-1]=='other':
      print '// t_{:5} E_{:7} ({:1})    Class:{:6} OccPattern:{:9}'.format(\
                            code(indTens_H2[iteri][:2]),\
                            code(indExOp_H2[iteri][:2]),\
                            1,\
                            'other',
                            output_pattern(code_rearrange(indExOp_H2[iteri][:2])))



# CALCULATIONS ------------------------------------------------------------------------------------
print '//'
print '// The actual calculation yields:'
global_result=[]
log_H4={}
log_H2={}

# Loop over Classes
for this_class in ClassNames+['other']:
  if this_class=='avca':
    continue
  print '// for V of class '+this_class
  result=[]

  # Loop over all H4tau of 'd_mu d_nu w_tau <O^T_nu V_tau O_mu>'
  for tau in range(len(indExOp_H4)):
    listtau=indExOp_H4[tau][:4]

    # This is to include "heff" etc... into a class,
    # so long as the pattern matches...
    if this_class!='other':
      if not (pattern(code_rearrange(listtau))==occupations[this_class]):
        continue
    # ... except is the class is 'other',
    # in this case the pattern is diverse,
    # and we include all "other" signatures
    elif indExOp_H4[tau][-1]!='other':
      continue

    # Loop over all Omu and Onu of 'd_mu d_nu w_tau <O^T_nu V_tau O_mu>'
    for mu in range(len(indExOp_Psi1)):
      for nu in range(len(indExOp_Psi1)):

        #if mu==nu:
        #  continue
        #if ((indExOp_Psi1[mu][-1]=='vaca' and indExOp_Psi1[nu][-1]=='avca')\
        # or (indExOp_Psi1[mu][-1]=='avca' and indExOp_Psi1[nu][-1]=='vaca')):
        #  continue

        listmu=indExOp_Psi1[mu][:4]
        listnu=indExOp_Psi1[nu][4:8]
        if is_non_zero(code(listnu)+code_rearrange(listtau)+code(listmu)):
          log_H4[tau]='evaluated'

          # V = Atau[0:4] Htau[0:4]
          if indExOp_H4[tau][-1]=='other' or indExOp_H4[tau][-1]=='h0':
            name='AO'+indTens_H4[tau][-1]
          elif indExOp_H4[tau][-1]=='heff4':
            name='AE'+indTens_H4[tau][-1]
          else:
            name='AC'+indTens_H4[tau][-1]
          Atau = sqa.tensor(name, indTens_H4[tau][:4], symm(indTens_H4[tau][:4]))
          Vtau = sqa.term( factor[code(listtau)]*0.5, [''],  [Atau, sqa.sfExOp(listtau)])

          # |Psi1> = dmu[0:4] Omu[01] Omu[23]
          Dmu    = sqa.tensor('D'+indTens_Psi1[mu][-1], indTens_Psi1[mu][:4], [])#symm(indTens_Psi1[mu][:4]))
          Omu    = sqa.term( 1.0, [''],  [Dmu, sqa.sfExOp(listmu[:2]), sqa.sfExOp(listmu[2:4])])

          # <Psi1| = dnu[4:8] Onu[45] Onu[67]
          Dnu    = sqa.tensor('D'+indTens_Psi1[nu][-1], indTens_Psi1[nu][4:8], [])#symm(indTens_Psi1[nu][4:8]))
          Onu    = sqa.term( 1.0, [''],  [Dnu, sqa.sfExOp(listnu[:2]), sqa.sfExOp(listnu[2:4])])

          # Define the normal-ordered:   A_nu = D^mu D^nu <O^T_nu V_tau O_mu>
          term=sqa.multiplyTerms(sqa.multiplyTerms(Onu, Vtau),  Omu)
          result=sqa.normalOrder(term)
          result=simplify_all(result,[deltaC, deltaA, deltaV])

          # Write out a message about this contribution
          if (len(result)!=0):
            global_result+=result
            print '// <Psi1_{:5} | V_{:5} | Psi1_{:5}> ={:5}    Class:{:6}'.format(\
                        indExOp_Psi1[nu][-1],\
                        code(listtau),\
                        indExOp_Psi1[mu][-1],\
                        len(result),\
                        indExOp_H4[tau][-1])




  # Loop over all H2tau of 'd_mu d_nu w_tau <O^T_nu V_tau O_mu>'
  for tau in range(len(indExOp_H2)):
    listtau=indExOp_H2[tau][:2]

    # This is to include "heff" etc... into a class,
    # so long as the pattern matches...
    if this_class!='other':
      if not (pattern(code_rearrange(listtau))==occupations[this_class]):
        continue
    # ... except is the class is 'other',
    # in this case the pattern is diverse,
    # and we include all "other" signatures
    elif indExOp_H2[tau][-1]!='other':
      continue

    # Loop over all Omu and Onu of 'd_mu d_nu w_tau <O^T_nu V_tau O_mu>'
    for mu in range(len(indExOp_Psi1)):
      for nu in range(len(indExOp_Psi1)):

        #if mu==nu:
        #  continue
        #if ((indExOp_Psi1[mu][-1]=='vaca' and indExOp_Psi1[nu][-1]=='avca')\
        # or (indExOp_Psi1[mu][-1]=='avca' and indExOp_Psi1[nu][-1]=='vaca')):
        #  continue

        listmu=indExOp_Psi1[mu][:4]
        listnu=indExOp_Psi1[nu][4:8]
        if is_non_zero(code(listnu)+code_rearrange(listtau)+code(listmu)):
          log_H2[tau]='evaluated'

          # V = Atau[0:2] Htau[0:2]
          if indExOp_H2[tau][-1]=='other' or indExOp_H2[tau][-1]=='h0':
            name='AO'+indTens_H2[tau][-1]
          elif indExOp_H2[tau][-1]=='heff2':
            name='AE'+indTens_H2[tau][-1]
          else:
            name='AC'+indTens_H2[tau][-1]
          Atau = sqa.tensor(name, indTens_H2[tau][:2], symm(indTens_H2[tau][:2]))
          Ttau   = sqa.term( 1.0, [''],  [Atau, sqa.sfExOp(listtau)])

          # |Psi1> = dmu[0:4] Omu[01] Omu[23]
          Dmu    = sqa.tensor('D'+indTens_Psi1[mu][-1], indTens_Psi1[mu][:4], [])#symm(indTens_Psi1[mu][:4]))
          Omu    = sqa.term( 1.0, [''],  [Dmu, sqa.sfExOp(listmu[:2]), sqa.sfExOp(listmu[2:4])])

          # <Psi1| = dnu[4:8] Onu[45] Onu[67]
          Dnu    = sqa.tensor('D'+indTens_Psi1[nu][-1], indTens_Psi1[nu][4:8], [])#symm(indTens_Psi1[nu][4:8]))
          Onu    = sqa.term( 1.0, [''],  [Dnu, sqa.sfExOp(listnu[:2]), sqa.sfExOp(listnu[2:4])])

          # Define the normal-ordered:   A_nu = D^mu D^nu <O^T_nu T_tau O_mu>
          term=sqa.multiplyTerms(sqa.multiplyTerms(Onu, Ttau),  Omu)
          result=sqa.normalOrder(term)
          result=simplify_all(result,[deltaC, deltaA, deltaV])

          # Write out a message about this contribution
          if (len(result)!=0):
            global_result+=result
            print '// <Psi1_{:5} | T_{:5} | Psi1_{:5}> ={:5}    Class:{:6}'.format(\
                        indExOp_Psi1[nu][-1],\
                        code(listtau),\
                        indExOp_Psi1[mu][-1],\
                        len(result),\
                        indExOp_H2[tau][-1])




# OUTPUT ------------------------------------------------------------------------------------------
# Write out the tensors and the equations
#global_result=simplify_all(global_result,[deltaC, deltaA, deltaV])

ACEOtensor=[]
Wktensor=[]
for t in global_result:
  for i in t.tensors:
    if i.name[0]=='A':
      if i.name not in [this[0] for this in ACEOtensor]:
        ACEOtensor.append([i.name,         i.name[:2], i.name[2:].replace('v','e'), 'A'])
        if len(i.name)==6 and 'W'+i.name[2:] not in [this[0] for this in Wktensor]:
          Wktensor.append(['W'+i.name[2:], 'W',        i.name[2:].replace('v','e'), 'A'])
        elif len(i.name)==4 and 'k'+i.name[2:] not in [this[0] for this in Wktensor]:
          Wktensor.append(['k'+i.name[2:], 'k',        i.name[2:].replace('v','e'), 'A'])
sort_indACEO=sorted(xrange(len(ACEOtensor)), key=lambda ix: ACEOtensor[ix][0])
sort_indWk  =sorted(xrange(len(Wktensor)), key=lambda ix: Wktensor[ix][0])
tensors=[ACEOtensor[i] for i in sort_indACEO]+[Wktensor[i] for i in sort_indWk]+tensors

if len(global_result)!=0:
  AllTensors,nbrTensors,nbrEqs = write_tensors_and_equations('MRLCC3',global_result,tensors)

  # Write out GetMethodInfo
  print '  static void GetMethodInfo(FMethodInfo &Out) {'
  print '    Out = FMethodInfo();'
  print '    Out.pName = "MRLCC3";'
  print '    Out.pSpinClass = "restricted";'
  print '    Out.pTensorDecls = &TensorDecls[0];'
  print '    Out.nTensorDecls = %i;'%(nbrTensors)
  print '    Out.EqsRes = FEqSet(&EqsRes[0], %i, "MRLCC3/Res");'%(len(global_result))
  print '  };'
  print '};'




# LOG ---------------------------------------------------------------------------------------------
print ''
print '// Log:'
for tau in range(len(indExOp_H4)):
  if tau in log_H4.keys():
    print '// {:4} Class:{:5}    {:}'.format(code(indExOp_H4[tau][:4]),indExOp_H4[tau][-1],log_H4[tau])
  else:
    print '// {:4} Class:{:5}    {:}'.format(code(indExOp_H4[tau][:4]),indExOp_H4[tau][-1],'passed')
for tau in range(len(indExOp_H2)):
  if tau in log_H2.keys():
    print '// {:4} Class:{:5}    {:}'.format(code(indExOp_H2[tau][:2]),indExOp_H2[tau][-1],log_H2[tau])
  else:
    print '// {:4} Class:{:5}    {:}'.format(code(indExOp_H2[tau][:2]),indExOp_H2[tau][-1],'passed')




