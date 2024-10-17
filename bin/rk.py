'''
Author: Zhuorui, Mintao
Version: 0.2
Date: 2022, Nov, 11
'''
import math
from pathlib import Path
import json

class motifScore(object):
    def __init__(self, boxSize, scorePath, JSONpath):
        #boxSize is the length for motif, scorePath is the path to the scoreMatrix; JSONpath is to the json file saving the hashcode version matrix
        self.boxSize = boxSize
        self.RBPmatrix = {}
        self.curSeq = 0
        self.HashBase = 10 ** boxSize
        self.scoreDic = {"A": 1, "C": 2, "G": 3, "T": 4, "N": 4,"a": 1, "c": 2, "g":3, "t":4, "n":4} 
        #for RBPmatrix, the key word is RBP 
        #the value is a dic, for this dic, key is sequence(in hashcode) whose size is boxSize, eg ATCGACA(size == 7) => in hash mode(1234131), value is score
        #A == 1  C == 2 G == 3 T == 4
        JSON = Path(JSONpath)
        if JSON.exists():
            #read JSON file
            with open(JSONpath) as json_file:
                self.RBPmatrix = json.load(json_file)
        else:
            self.geneJSON(scorePath, JSONpath)
    
    def geneJSON(self, scorePath, JSONpath):
        def geneHash(curL: int, matrix: dict, curScore: float, curSeq: str, RBPname: str):
            if curL == len(matrix):
                self.RBPmatrix[RBPname][curSeq] = curScore
                return
            geneHash(curL + 1, matrix, curScore + math.log2(4 * matrix[curL][0]), curSeq + "1", RBPname)
            geneHash(curL + 1, matrix, curScore + math.log2(4 * matrix[curL][1]), curSeq + "2", RBPname)
            geneHash(curL + 1, matrix, curScore + math.log2(4 * matrix[curL][2]), curSeq + "3", RBPname)
            geneHash(curL + 1, matrix, curScore + math.log2(4 * matrix[curL][3]), curSeq + "4", RBPname)
        
        matrix = {}
        human = -1
        with open(scorePath, "r") as f:
            lines = f.readlines()
            for line in lines:
                #title line
                if line.startswith('>'):
                    #in case it is not a human motif
                    if "Homo" in line:
                        human = 1
                    else:
                        human = -1
                        continue
                    index1 = line.find('(')
                    RBPname = line[13:index1]
                    matrix[RBPname] = []
                elif human == 1:
                    str_list = line.split()
                    float_list = []
                    for i in str_list:
                        float_list.append(float(i))
                    matrix[RBPname].append(float_list)
        
        #print(matrix['RBFOX1'])
        for RBPname in matrix:
            self.RBPmatrix[RBPname] = {}
            geneHash(0, matrix[RBPname], 0, "", RBPname)
            #print(self.RBPmatrix[RBPname])
        
        # Serializing json
        json_object = json.dumps(self.RBPmatrix, indent=4)
 
        # Writing to sample.json
        with open(JSONpath, "w") as outfile:
            outfile.write(json_object)
    
    def askScore(self, newBase, boxSize: int, RBPname):
        #insert a new gene Base
        #return: None, current seq is not a motif in this boxSize; Score for the probability of this RBP is motivated by this seq
        assert newBase in self.scoreDic
        assert boxSize == self.boxSize
        
        self.curSeq = self.curSeq * 10 % self.HashBase + self.scoreDic[newBase]
        #print(self.curSeq)
        if str(self.curSeq) not in self.RBPmatrix[RBPname]:
            return None
        return self.RBPmatrix[RBPname][str(self.curSeq)]

    def getRBPlist(self):
        return list(self.RBPmatrix.keys())

    def initSeq(self):
        self.curSeq = 0

    



if __name__ == "__main__":
    test = motifScore(7, "/restricted/projectnb/casa/bu_brain_rnaseq/hippo_rawdata/circrna/homer/data/knownTFs/known.rna.motifs", "./tmpJSON.json")
    testseq = "AGCATGT"
    for base in testseq:
        print(test.askScore(base, 7, 'RBFOX1'))
    print(test.getRBPlist())