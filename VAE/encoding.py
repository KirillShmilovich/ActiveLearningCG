import numpy as np
from martini22_ff import martini22

aminoAcidAlphabet = 'arndcqeghilkmfpstwyv'
#aminoAcidAlphabet = 'ancqeghilmfpstwyv'
aminoAcidDict = { 'g':'GLY','a':'ALA','l':'LEU',
                 'm':'MET','f':'PHE','w':'TRP',
                 'k':'LYS','q':'GLN','e':'GLU',
                 's':'SER','p':'PRO','v':'VAL',
                 'i':'ILE','y':'TYR','h':'HIS',
                 'r':'ARG','n':'ASN','d':'ASP','t':'THR','c':'CYS'}
with open('martini.itp','r') as f:
    martiniItp = f.read().splitlines()
# Dict of all the possible interactions between beads
interactionDict = { 'supraattractive':2,
                    'attractive':3,
                    'almostattractive':4,
                    'semiattractive':5,
                    'intermediate':6,
                    'almostintermediate':7,
                    'semirepulsive':8,
                    'almostrepulsive':9,
                    'repulsive':10,
                    'superrepulsive':11}
# Maximum number of beads in one martini residue
beadDict = { 'P1' :1,
             'P4' :2,
             'P5' :3,
             'SC4':4,
             'SC5':5,
             'SP1':6,
             'SNd':7,
             'N0' :8,
             'Qd' :9,
             'Qa' :10,
             'C3' :11,
             'C5' :12,
             'AC1':13,
             'AC2':14 }
NUM_BEADS = 5

class Processing():

    def __init__(self, input_wing):

        # Ensure input is lowercase, for indexing
        input_wing = input_wing.lower()
        wing_len = len(input_wing)
        aminoAcidString = [aminoAcidDict[char] for char in input_wing] 
        self.scaledMatrix = list()
        self.beadArray = np.zeros(wing_len*NUM_BEADS)
        self.mols = 0
        tempMatrix = np.zeros((wing_len*NUM_BEADS,wing_len*NUM_BEADS))
        #tempMatrix = np.zeros((len(beadDict),wing_len*NUM_BEADS,wing_len*NUM_BEADS))
        k=0
        self.P1Bool = 0
        self.P3Bool = 0
        self.P4Bool = 0
        self.P5Bool = 0
        self.SC4Bool = 0
        self.SC5Bool = 0
        self.SP1Bool = 0
        self.SNdBool = 0
        self.SQdBool = 0
        self.N0Bool = 0
        self.QdBool = 0
        self.QaBool = 0
        self.C3Bool = 0
        self.C5Bool = 0
        self.AC1Bool = 0
        self.AC2Bool = 0
        for residue in aminoAcidString:
            #tempMatrix = np.zeros((NUM_BEADS,NUM_BEADS))
            try:
                bonds = martini22().connectivity[residue][0]
                beads = [martini22().bbGetBead(residue)]+martini22().sidechains[residue][0]
                molNum = np.amax(bonds)+1
            except IndexError:
                bonds = []
                beads = [martini22().bbGetBead(residue)]
                molNum = 1
            self.mols += molNum
            if 'P1' in beads:
                self.P1Bool = 1
            if 'P3' in beads:
                self.P3Bool = 1
            if 'P4' in beads:
                self.P4Bool = 1
            if 'P5' in beads:
                self.P5Bool = 1
            if 'SC4' in beads:
                self.SC4Bool = 1
            if 'SC5' in beads:
                self.SC5Bool = 1
            if 'SP1' in beads:
                self.SP1Bool = 1
            if 'SNd' in beads:
                self.SNdBool = 1
            if 'SQd' in beads:
                self.SQdBool = 1
            if 'N0' in beads:
                self.N0Bool = 1
            if 'Qd' in beads:
                self.QdBool = 1
            if 'Qa' in beads:
                self.QaBool = 1
            if 'C3' in beads:
                self.C3Bool = 1
            if 'C5' in beads:
                self.C5Bool = 1
            if 'AC1' in beads:
                self.AC1Bool = 1
            if 'AC2' in beads:
                self.AC2Bool = 1
            for i in range(molNum):
                self.beadArray[i+k*NUM_BEADS]=beadDict[beads[i]]
                for j in range(molNum):
                    if i==j: pass 
                        #tempMatrix[i+k*NUM_BEADS][j+k*NUM_BEADS]=beadDict[beads[i]]
                        #tempMatrix[i+k*NUM_BEADS][j+k*NUM_BEADS]=1
                    else: pass
                        #tempMatrix[i+k*NUM_BEADS][j+k*NUM_BEADS]=1
            for i in bonds:
                #bead_1 = beads[i[0]]
                #bead_2 = beads[i[1]]
                tempMatrix[i[0]+k*NUM_BEADS][i[1]+k*NUM_BEADS]=1
                tempMatrix[i[1]+k*NUM_BEADS][i[0]+k*NUM_BEADS]=1
            k+=1
        if wing_len>1:
            for i in range(wing_len-1):
                tempMatrix[i*NUM_BEADS][(i+1)*NUM_BEADS]=1
                tempMatrix[(i+1)*NUM_BEADS][i*NUM_BEADS]=1
        self.scaledMatrix = np.array(tempMatrix)
        self.scaledArray = list()
        self.bianaryArray = list()
        for i in range(wing_len*NUM_BEADS):
            for j in range(i,wing_len*NUM_BEADS):
                self.scaledArray.append(self.scaledMatrix[i,j])
                if (i==j):
                    oneHotVector = [0]*(1)
                else:
                    oneHotVector = [0]*(1)
                if int(self.scaledMatrix[i,j])==0:
                    self.bianaryArray = self.bianaryArray + oneHotVector
                else:
                    oneHotVector[0]=1
                    self.bianaryArray = self.bianaryArray + oneHotVector
        for i in range(wing_len*NUM_BEADS):
            oneHotVector = [0]*(len(beadDict))
            if int(self.beadArray[i])!=0:
                oneHotVector[int(self.beadArray[i])-1]=1
            self.bianaryArray = self.bianaryArray + oneHotVector
        self.scaledArray = np.array(self.scaledArray)
        self.bianaryArray = np.array(self.bianaryArray)
    def unscale(self, inputMatrix):
        return inputMatrix
    def toMatrix(self, bianaryArray):
        k = 0
        #CHANGE THIS TO BE DYNAMIC
        outMatrix = np.zeros((3,NUM_BEADS,NUM_BEADS)) 
        for l in range(3):
            for i in range(NUM_BEADS):
                for j in range(i,NUM_BEADS):
                    if i==j:
                        tempVec = bianaryArray[k:k+(len(beadDict))]
                        k+=len(beadDict)
                        if sum(tempVec)>1:
                            print('We got an issue here')
                            print(tempVec)
                        if sum(tempVec)==1:
                           outMatrix[l,i,j]=np.where(tempVec==1)[0][0]+1 
                           outMatrix[l,j,i]=np.where(tempVec==1)[0][0]+1 
                    else:
                        tempVec = bianaryArray[k:k+len(interactionDict)+1]
                        k += len(interactionDict)+1
                        if sum(tempVec)>1:
                            print('We got an issue here')
                            print(tempVec)
                        if sum(tempVec)==1:
                           outMatrix[l,i,j]=np.where(tempVec==1)[0][0]+1 
                           outMatrix[l,j,i]=np.where(tempVec==1)[0][0]+1 
        return outMatrix
def main(): pass 
if __name__ == '__main__':
   main() 
