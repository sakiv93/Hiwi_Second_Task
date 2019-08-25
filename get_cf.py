from odbAccess import *
from time import *

import sys

lt = localtime()
print 'start get_cf' + strftime("   %d.%m.%Y, %H:%M:%S", lt)

odb_name = sys.argv[1]

odb = openOdb(path=odb_name+'.odb', readOnly=False)
#
#
##### Parameters for gauss integration####
# coordinate of integration point
#$ipcoor = 0.577350269189626
ipcoor= 0.7745966692414834
# gauss weight
#$gweight = 1.
gweight_1 = 0.8888888888888889 #gauss weight for coordinate 0
gweight_2 = 0.5555555555555556 #gauss weight for coordinate 0.7745966692414834
##########################################
# initialize local coordinates of integration points
#$allxi =  [ -ipcoor, ipcoor, -ipcoor, ipcoor, -ipcoor, ipcoor, -ipcoor, ipcoor]
#$alleta =  [ -ipcoor, -ipcoor, ipcoor, ipcoor, -ipcoor, -ipcoor, ipcoor, ipcoor]
#$allzeta = [ -ipcoor, -ipcoor, -ipcoor, -ipcoor, ipcoor, ipcoor, ipcoor, ipcoor]
allxi =  [ -ipcoor, 0, ipcoor, -ipcoor, 0, ipcoor, -ipcoor, 0, ipcoor, -ipcoor, 0, ipcoor, -ipcoor, 0, ipcoor, -ipcoor, 0, ipcoor, -ipcoor, 0, ipcoor, -ipcoor, 0, ipcoor, -ipcoor, 0,, ipcoor]
alleta =  [ -ipcoor, -ipcoor, -ipcoor, 0, 0, 0, ipcoor, ipcoor, ipcoor, -ipcoor, -ipcoor, -ipcoor, 0, 0, 0, ipcoor, ipcoor, ipcoor, -ipcoor, -ipcoor, -ipcoor, 0, 0, 0, ipcoor, ipcoor, ipcoor]
allzeta = [ -ipcoor, -ipcoor, -ipcoor, -ipcoor, -ipcoor, -ipcoor, -ipcoor, -ipcoor, -ipcoor, 0, 0, 0, 0, 0, 0, 0, 0, 0, ipcoor, ipcoor, ipcoor, ipcoor, ipcoor, ipcoor, ipcoor, ipcoor, ipcoor]

# initialize local coordinates of nodes
#$allnodexi = [ -1., 1., 1., -1., -1., 1., 1., -1.]
#$allnodeeta = [ -1., -1., 1., 1., -1., -1., 1., 1. ]
#$allnodezeta = [ -1., -1., -1., -1., 1., 1., 1., 1. ]
allnodexi = [ -1., 1., 1., -1., -1., 1., 1., -1., 0., 1., 0., -1., 0., 1., 0., -1., -1., 1., 1., -1. ]
allnodeeta = [ -1., -1., 1., 1., -1., -1., 1., 1., -1., 0., 1., 0., -1., 0., 1., 0., -1., -1., 1., 1. ]
allnodezeta = [ -1., -1., -1., -1., 1., 1., 1., 1., -1., -1., -1., -1., 1., 1., 1., 1., 0., 0., 0., 0. ]
# Kronecker delta 
kronkj = [[ 1., 0., 0.],[0., 1., 0.],[0., 0., 1.]]
#
for partkey in odb.parts.keys():
    part=odb.parts[partkey]
    if len(part.nodes)>0:
# make list "elatnode" containing element labels of all elements at an specific node
        elatnode=[]
        nodelabels=[]
        for i in range(len(part.nodes)):
            elatnode.append([])
            nodelabels.append(i+1)
        elpos=0
        for elem in part.elements:
            for node in elem.connectivity:
                elatnode[node-1].append(elpos)
            elpos=elpos+1
#
###############################################################################################
#
        firsttime = 0
        indexfield = []
        for elem in part.elements:
            #$indexfield.append([0,0,0,0,0,0,0,0])
            indexfield.append([0,0,0,0,0,0,0,0])
        for s in odb.steps.keys():
	    step = odb.steps[s]
	    for frame in step.frames:
# If first time get positions of elements and gauss points in field output
                if (firsttime == 0) :
                    print 'start get positions of elements and gauss points in field output'
                    for index in range(len(part.elements)*8):
                        indexelement = frame.fieldOutputs['S'].values[index].elementLabel - 1
                        indexgauss = frame.fieldOutputs['S'].values[index].integrationPoint - 1
                        indexfield[indexelement][indexgauss] = index
#                        print index
                    print 'get positions of elements and gauss points in field output completed'
                firsttime = 1
                lt = localtime()
                print 'start processing increment ' + str(frame.incrementNumber) + strftime("   %d.%m.%Y, %H:%M:%S", lt)
# Gnode - configurational force at node
                Gnode=[]
# Gelnode - contribution of element to configurational force at node
                Gelnode=[]
                for elem in part.elements:
                    Gelnode.append([[0., 0., 0.], [0., 0., 0.], [0., 0., 0.], [0., 0., 0.], [0., 0., 0.], [0., 0., 0.], [0., 0., 0.],[0., 0., 0.]])
                for elem in part.elements:
#                    print elem.label
                    #$for Gp in range(8):
                    for Gp in range(27):
# Get values at integration point
                        index = indexfield[elem.label-1][Gp]
                        stress = frame.fieldOutputs['S'].values[index].data
                        strain = frame.fieldOutputs['E'].values[index].data
                        eirr = frame.fieldOutputs['EIRR'].values[index].data
                        Pirr = frame.fieldOutputs['POL'].values[index].data
                        epg = frame.fieldOutputs['EPG'].values[index].data
                        eflx = frame.fieldOutputs['EFLX'].values[index].data
                        elin = [0., 0., 0., 0., 0., 0.]
                        DplusP = [0., 0., 0.]
# Calculate electric enthalpy density h
                        stresselin = 0.
                        for i in range( 6):
                            elin[i] = strain[i] - eirr[i]
                            stresselin = stresselin + (stress[i] * elin[i])
                        DplusPepg = 0.
                        for i in range(3):
                            DplusP[i] = eflx[i] + Pirr[i]
                            DplusPepg = DplusPepg + (DplusP[i] * epg[i])
                        h = 0.5 * stresselin - 0.5 * DplusPepg
#
###############################################################################################
#
# global derivatives of shape funktions at integration point
# 1. local coordinates of integration point
                        xi = allxi[Gp]
                        eta = alleta[Gp]
                        zeta = allzeta[Gp]
#                        print xi
#                        print etaindexfield.append([0,0,0,0,0,0,0,0])
#                        print zeta
# 2. local coordinates of nodes and local partial derivatives of shape functions
                        Nxi = [0., 0., 0., 0., 0., 0., 0.,0.,0., 0., 0., 0., 0., 0., 0.,0.,0., 0., 0., 0.]
                        Neta = [0., 0., 0., 0., 0., 0., 0.,0.,0., 0., 0., 0., 0., 0., 0.,0.,0., 0., 0., 0.]
                        Nzeta = [0., 0., 0., 0., 0., 0., 0.,0.,0., 0., 0., 0., 0., 0., 0.,0.,0., 0., 0., 0.]
                        for node in range(20):
                            nodexi = allnodexi[node]
                            nodeeta = allnodeeta[node]
                            nodezeta = allnodezeta[node]
#                            print 'node coordinates'
#                            print nodexi
#                            print nodeeta
#                            print nodezeta
#                            print (1./8.)*nodexi*(1.+nodeeta*eta)*(1.+nodezeta*zeta)
                            if(nodexi!=0 and nodeeta!=0 and nodezeta!=0):
                                #$Nxi[node] = (1./8.)*nodexi*(1.+nodeeta*eta)*(1.+nodezeta*zeta)
                                Nxi[node]= (1./8.)*nodexi*(1.+nodeeta*eta)*(1.+nodezeta*zeta)*(nodexi*xi+nodeeta*eta+nodezeta*zeta-2) + (1./8.)*(1.+nodexi*xi)*(1.+nodeeta*eta)*(1.+nodezeta*zeta)*(nodexi)
                                #$Neta[node] = (1./8.)*(1.+nodexi*xi)*nodeeta*(1.+nodezeta*zeta)
                                Nxi[node]= (1./8.)*(1.+nodexi*xi)*(nodeeta)*(1.+nodezeta*zeta)*(nodexi*xi+nodeeta*eta+nodezeta*zeta-2) + (1./8.)*(1.+nodexi*xi)*(1.+nodeeta*eta)*(1.+nodezeta*zeta)*(nodeeta)
                                #$Nzeta[node] = (1./8.)*(1.+nodexi*xi)*(1.+nodeeta*eta)*nodezeta
                                Nxi[node]= (1./8.)*(1.+nodexi*xi)*(1.+nodeeta*eta)*nodezeta*(nodexi*xi+nodeeta*eta+nodezeta*zeta-2) + (1./8.)*(1.+nodexi*xi)*(1.+nodeeta*eta)*(1.+nodezeta*zeta)*(nodezeta)
                            elif(nodexi==0):
                                Nxi[node] = (-1./2.)*nodexi*(1.+nodeeta*eta)*(1.+nodezeta*zeta)
                            elif(nodeeta==0):
                                Nxi[node]= (-1./2.)*(1.+nodexi*xi)*(nodeeta)*(1.+nodezeta*zeta)
                            elif(nodeeta==0):
                                Nxi[node]= (-1./2.)*(1.+nodexi*xi)*(1.+nodeeta*eta)*nodezeta


#                            print 'shape functions'
#                            print Nxi[node]
#                            print Neta[node]
#                            print Nzeta[node]
# 3. jacobian matrix
                        J=[[ 0., 0.,0.],[ 0., 0.,0.],[ 0., 0.,0.]]
                        for node in range(8):
                            nodelab=elem.connectivity[node]
                            x = part.nodes[nodelab-1].coordinates[0]
                            y = part.nodes[nodelab-1].coordinates[1]
                            z = part.nodes[nodelab-1].coordinates[2]
                            J[0][0] = J[0][0] + Nxi[node]*x
                            J[0][1] = J[0][1] + Nxi[node]*y
                            J[0][2] = J[0][2] + Nxi[node]*z
                            J[1][0] = J[1][0] + Neta[node]*x
                            J[1][1] = J[1][1] + Neta[node]*y
                            J[1][2] = J[1][2] + Neta[node]*z
                            J[2][0] = J[2][0] + Nzeta[node]*x
                            J[2][1] = J[2][1] + Nzeta[node]*y
                            J[2][2] = J[2][2] + Nzeta[node]*z
#                            print Nxi[node]
#                            print Neta[node]
#                            print Nzeta[node]
#                            print x
#                            print y
#                            print z
#                            print J
# 4. determinant of jacobian
                        detJ = J[0][0] * J[1][1] * J[2][2]\
                             + J[0][1] * J[1][2] * J[2][0]\
                             + J[0][2] * J[1][0] * J[2][1]\
                             - J[0][0] * J[1][2] * J[2][1]\
                             - J[0][1] * J[1][0] * J[2][2]\
                             - J[0][2] * J[1][1] * J[2][0]
#                        print detJ
# 5. inverse of jacobian
                        Jinv = [[ 0., 0.,0.],[ 0., 0.,0.],[ 0., 0.,0.]]
                        Jinv[0][0] = (1./detJ) * ( J[1][1] * J[2][2] - J[1][2] * J[2][1])
                        Jinv[0][1] = (1./detJ) * ( J[0][2] * J[2][1] - J[0][1] * J[2][2])
                        Jinv[0][2] = (1./detJ) * ( J[0][1] * J[1][2] - J[0][2] * J[1][1])
                        Jinv[1][0] = (1./detJ) * ( J[1][2] * J[2][0] - J[1][0] * J[2][2])
                        Jinv[1][1] = (1./detJ) * ( J[0][0] * J[2][2] - J[0][2] * J[2][0])
                        Jinv[1][2] = (1./detJ) * ( J[0][2] * J[1][0] - J[0][0] * J[1][2])
                        Jinv[2][0] = (1./detJ) * ( J[1][0] * J[2][1] - J[1][1] * J[2][0])
                        Jinv[2][1] = (1./detJ) * ( J[0][1] * J[2][0] - J[0][0] * J[2][1])
                        Jinv[2][2] = (1./detJ) * ( J[0][0] * J[1][1] - J[0][1] * J[1][0])
# 6. get global partial derivatives of shape functions
                        Nj = [[0., 0., 0., 0., 0., 0., 0., 0.],[0., 0., 0., 0., 0., 0., 0., 0.],[0., 0., 0., 0., 0., 0., 0., 0.]]
                        for node in range(8):
                            Nj[0][node] = Jinv[0][0] * Nxi[node] +Jinv[0][1] * Neta[node] +Jinv[0][2] * Nzeta[node]
                            Nj[1][node] = Jinv[1][0] * Nxi[node] +Jinv[1][1] * Neta[node] +Jinv[1][2] * Nzeta[node]
                            Nj[2][node] = Jinv[2][0] * Nxi[node] +Jinv[2][1] * Neta[node] +Jinv[2][2] * Nzeta[node]
# calculate contribution of element to configurational force at node Gelnode
# 1. calculate eshelby stress tensor at integration point ESTkj
                        ESTkj = [[0., 0., 0.],[0., 0., 0.],[0., 0., 0.]]
# 1.a get global derivatives of displacement ui,k
                        uik = [[0., 0., 0.],[0., 0., 0.],[0., 0., 0.]]
                        for node in range(8):
                            nodeindex = elem.connectivity[node]-1
#                            print nodeindex
                            u = frame.fieldOutputs['U'].values[nodeindex].data
                            for i in range(3):
                                for k in range(3):
                                    uik[i][k] = uik[i][k] + Nj[k][node]*u[i]
# 1.b define stress in standard notation
                        stressij = [[ stress[0], stress[3], stress[4]],[ stress[3], stress[1], stress[5]],[ stress[4], stress[5], stress[2]]]
# 1.c eshelby stress tensor
                        for k in range(3):
                            for j in range(3):
                                ESTkj[k][j] = h * kronkj[k][j]\
                                            - stressij[0][j] * uik[0][k]\
                                            - stressij[1][j] * uik[1][k]\
                                            - stressij[2][j] * uik[2][k]\
                                            + eflx[j] * epg[k]
# 2. calculate Gelnode
                        for node in range(8):
                            for k in range(3):
                                labelelem = elem.label-1
                                Gelnode[labelelem][node][k] = Gelnode[labelelem][node][k]\
                                                               + ESTkj[k][0] * Nj[0][node] * gweight * detJ\
                                                               + ESTkj[k][1] * Nj[1][node] * gweight * detJ\
                                                               + ESTkj[k][2] * Nj[2][node] * gweight * detJ
#
###############################################################################################
#
# Gnode = sum Gelnode
                for node in range(len(part.nodes)):
                    Gnode.append([0., 0., 0.])
                    for elempos in elatnode[node]:
                        elemen = part.elements[elempos]
                        elnodeposition=0
                        for elnode in elemen.connectivity:
                            if elnode==node+1:
                                for i in range(3):
                                    Gnode[node][i]=Gnode[node][i]+Gelnode[elemen.label-1][elnodeposition][i]
                            elnodeposition = elnodeposition + 1
# update odb
                lt = localtime()
                print 'start update odb' + strftime("   %d.%m.%Y, %H:%M:%S", lt)
                GField = frame.FieldOutput(name='G',description='Configurational Forces', type=VECTOR)
                GField.addData(position=NODAL, instance=odb.rootAssembly.instances[instancename], labels=nodelabels, data=Gnode)
lt = localtime()
print 'get_cf completed' + strftime("   %d.%m.%Y, %H:%M:%S", lt)

