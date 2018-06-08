import secondQuantizationAlgebra as sqa
from tools_for_sqa import *
from tools_for_pt import *
sqa.options.verbose = False







# PS1 ---------------------------------------------------------------------------------------------
#  |Psi1> =         d_mu O  _mu |Psi_0>
#  <Psi1| = <Psi_0| d_mu O^T_mu
#  where:
#   - the indexes of the basis are dictated by the class itself,
#     see that they match with the rest of the sqa files!
#   - the indexes of the tensor for each class are to match the unique tensors outputed by icpt!
#     (this is taken care of by the "unique_tensor" by construction)
#   - those tensors are labelled by "names" that are the classes
#     and not the correspondinng string of indexes types!
names=['vvcc','vvca','ccav','vvaa',\
       'ccaa','vaca','avca','vaaa','caaa']
basis=[[gimme('vc'),gimme('vc')],\
       [gimme('vc'),gimme('va')],\
       [gimme('ac'),gimme('vc')],\
       [gimme('va'),gimme('va')],\
       [gimme('ac'),gimme('ac')],\
       [gimme('vc'),gimme('aa')],\
       [gimme('ac'),gimme('va')],\
       [gimme('va'),gimme('aa')],\
       [gimme('ac'),gimme('aa')]]
Psi1Right=[]
Psi1Left =[]
for i in range(9):
  basis1=sqa.sfExOp(basis[i][0])
  basis2=sqa.sfExOp(basis[i][1])
  tensor=sqa.tensor('D'+names[i], unique_tensor(from_basis(basis1,basis2)), [])
  print '1:',code(unique_tensor (from_basis(basis1,basis2))),\
             code(unique_tensor_old(from_basis(basis1,basis2)))
  Psi1Right.append(sqa.term(1.0,[''],[tensor,basis1,basis2]))

  basis1=sqa.sfExOp([basis[i][1][t] for t in [1,0]])
  basis2=sqa.sfExOp([basis[i][0][t] for t in [1,0]])
  tensor=sqa.tensor('D'+names[i], unique_tensor(from_basis(basis1,basis2)), [])
  print '1:',code(unique_tensor (from_basis(basis1,basis2))),\
             code(unique_tensor_old(from_basis(basis1,basis2)))
  Psi1Left.append( sqa.term(1.0,[''],[tensor,basis1,basis2]))




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
# Construct the full Hamiltonian ExOp
# (The indexes run through ABCD.)
indlist=[gimme('cav'),gimme('cav'),gimme('cav'),gimme('cav')]

factor={}
for ind0 in range(3):
  for ind1 in range(3):
    for ind2 in range(3):
      for ind3 in range(3):
        factor[code([indlist[0][ind0],indlist[1][ind1],indlist[2][ind2],indlist[3][ind3]])]=0
    factor[code([indlist[0][ind0],indlist[1][ind1]])]=1

print '// To make up E4(AB<-CD), there is 3*3*3*3=81 possible 4-tuples'
print '// some of them are redundant are skipped.'
H=[]
iter=0
for ind0 in range(3):
  for ind1 in range(3):
    for ind2 in range(3):
      for ind3 in range(3):
        list=[indlist[0][ind0],indlist[1][ind1],indlist[2][ind2],indlist[3][ind3]]

        iter=iter+1
        name=code(list)
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

        exop  =sqa.sfExOp(list)
        tensor=sqa.tensor('AX'+code(unique_tensor(list)), unique_tensor(list), [])
        print '1:',code(unique_tensor (list)),\
                   code(unique_tensor_old(list))
        H.append(sqa.term(0.5,[''],[tensor,exop]))

    list=[indlist[0][ind0],indlist[1][ind1]]
    exop  =sqa.sfExOp(list)
    tensor=sqa.tensor('AX'+code(unique_tensor(list)), unique_tensor(list), [])
    print '1:',code(unique_tensor (list)),\
               code(unique_tensor_old(list))
    H.append(sqa.term(1.0,[''],[tensor,exop]))

for i in range(len(H)):
  el=H[i]
  string=code([t.indices for t in el.tensors if isinstance(t,sqa.sfExOp)])
  H[i]=sqa.term(el.numConstant*factor[string],\
                el.constants,\
                el.tensors)

print '// Sum of factors:',sum(factor.values())
print '//\n//'




# Print out stuff ---------------------------------------------------------------------------------
print '//'
print '// With the unique tensors found, |Psi1> is:'
for el in Psi1Right:
  print '// d_{:5} E_{:2} E_{:2}        Class:{:6} OccPattern:{:9}'.format(\
                        code([i.indices for i in el.tensors if not isinstance(i,sqa.sfExOp)][0]),\
                        code([i.indices for i in el.tensors if isinstance(i,sqa.sfExOp)][0]),\
                        code([i.indices for i in el.tensors if isinstance(i,sqa.sfExOp)][1]),\
                        [i.name for i in el.tensors if not isinstance(i,sqa.sfExOp)][0][1:],\
                        '%3i%3i%3i'%(tuple(pattern([el]))))
print '//'
print '// With the unique tensors found, <Psi1| is:'
for el in Psi1Left:
  print '// d_{:5} E_{:2} E_{:2}        Class:{:6} OccPattern:{:9}'.format(\
                        code([i.indices for i in el.tensors if not isinstance(i,sqa.sfExOp)][0]),\
                        code([i.indices for i in el.tensors if isinstance(i,sqa.sfExOp)][0]),\
                        code([i.indices for i in el.tensors if isinstance(i,sqa.sfExOp)][1]),\
                        [i.name for i in el.tensors if not isinstance(i,sqa.sfExOp)][0][1:],\
                        '%3i%3i%3i'%(tuple(pattern([el]))))
print '//'
print '//'
print '// With the unique tensors found, the Hamiltonian is:'
ClassNames=['vvcc','vvca','ccav','vvaa','ccaa','vaca','vaaa','caaa']
for cl in ClassNames+['h0']:#,'heff2','heff4']:
  for el in H:
    if pattern([el])==pattern_of[cl]:
      print '// w_{:5} E_{:7} ({:1})    Class:{:6} OccPattern:{:9}'.format(\
                        code([i.indices for i in el.tensors if not isinstance(i,sqa.sfExOp)][0]),\
                        code([i.indices for i in el.tensors if isinstance(i,sqa.sfExOp)][0]),\
                        factor[code([i.indices for i in el.tensors if isinstance(i,sqa.sfExOp)][0])],\
                        cl,\
                        '%3i%3i%3i'%(tuple(pattern([el]))))
for el in H:
  if pattern([el]) not in pattern_of.values():
      print '// w_{:5} E_{:7} ({:1})    Class:{:6} OccPattern:{:9}'.format(\
                        code([i.indices for i in el.tensors if not isinstance(i,sqa.sfExOp)][0]),\
                        code([i.indices for i in el.tensors if isinstance(i,sqa.sfExOp)][0]),\
                        factor[code([i.indices for i in el.tensors if isinstance(i,sqa.sfExOp)][0])],\
                        'other',\
                        '%3i%3i%3i'%(tuple(pattern([el]))))



# CALCULATIONS ------------------------------------------------------------------------------------
print '//'
print '// The actual calculation yields:'
global_result=[]
log_H={}

# Loop over Classes
for this_class in ClassNames+['other']:
  print '// for V of class '+this_class
  result=[]

  # Loop over all H4tau of 'd_mu d_nu w_tau <O^T_nu V_tau O_mu>'
  for hTerm in H:
    # This is to include "heff" etc... into a class,
    # so long as the pattern matches...
    if this_class!='other':
      if pattern([hTerm])!=pattern_of[this_class]:
        continue
    # ... except is the class is 'other',
    # in this case the pattern is diverse,
    # and we include all "other" signatures
    elif pattern([hTerm]) in pattern_of.values():
      continue

    # Loop over all Omu and Onu of 'd_mu d_nu w_tau <O^T_nu V_tau O_mu>'
    for leftTerm in Psi1Left:
      for rightTerm in Psi1Right:

        if is_non_zero([leftTerm,hTerm,rightTerm]):
          log_H[H.index(hTerm)]='evaluated'

          ## V = Atau[0:4] Htau[0:4]
          #if indExOp_H4[tau][-1]=='other' or indExOp_H4[tau][-1]=='h0':
          #  name='AO'+indTens_H4[tau][-1]
          #elif indExOp_H4[tau][-1]=='heff4':
          #  name='AE'+indTens_H4[tau][-1]
          #else:
          #  name='AC'+indTens_H4[tau][-1]
          #Atau = sqa.tensor(name, indTens_H4[tau][:4], symm(indTens_H4[tau][:4]))
          #Vtau = sqa.term( factor[code(listtau)]*0.5, [''],  [Atau, sqa.sfExOp(listtau)])

          # Define the normal-ordered:   A_nu = D^mu D^nu <O^T_nu V_tau O_mu>
          term=sqa.multiplyTerms(sqa.multiplyTerms(leftTerm, hTerm),  rightTerm)
          result=sqa.normalOrder(term)
          result=simplify_all(result,[deltaC, deltaA, deltaV])

          # Write out a message about this contribution
          if (len(result)!=0):
            global_result+=result
            print '// <Psi1_{:5} | V_{:5} | Psi1_{:5}> ={:5}    Class:{:6}'.format(\
                        code([i.indices for i in leftTerm.tensors  if isinstance(i,sqa.sfExOp)]),\
                        code([i.indices for i in hTerm.tensors     if isinstance(i,sqa.sfExOp)]),\
                        code([i.indices for i in rightTerm.tensors if isinstance(i,sqa.sfExOp)]),\
                        len(result),\
                        this_class)




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



exit()

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




