import math
import matplotlib.pyplot as plt


def main():
    protein = ReadProteinFromPDB('proteins/e3e59A1_AB.pdb')
    coordinatesToPlot = []

    for chain in protein:
        for i in range(1, len(chain)-1):
            currPhiAngle = CalculateAngleByAminoAcid(chain[i-1], chain[i], 1)
            currPsiAngle = CalculateAngleByAminoAcid(chain[i+1], chain[i], 0)
            coordinateToPush = (currPhiAngle, currPsiAngle)
            coordinatesToPlot.append(coordinateToPush)
    for coordinate in coordinatesToPlot:
        plt.scatter(coordinate[0], coordinate[1], 12)
    plt.xlim(-180, 180)
    plt.ylim(-180, 180)
    plt.title('Ramachandran Plot for e3e59A1_AB.pdb')
    plt.ylabel('Psi')
    plt.xlabel('Phi')
    plt.show()


def CalculateAngleByAminoAcid(aminoAcid1, aminoAcid2, isPhi):
    if isPhi:
        p0 = aminoAcid1[2]  # C of prev
        p1 = aminoAcid2[0]  # N
        p2 = aminoAcid2[1]  # CA
        p3 = aminoAcid2[2]  # C
    else:
        p0 = aminoAcid2[0]  # N
        p1 = aminoAcid2[1]  # CA
        p2 = aminoAcid2[2]  # C
        p3 = aminoAcid1[0]  # N of next

    r1 = SubtractVectors3D(p0, p1)
    r2 = SubtractVectors3D(p2, p1)
    r3 = SubtractVectors3D(p1, p2)
    r4 = SubtractVectors3D(p3, p2)

    n1 = CalcVectorCrossProduct3D(r1, r2)
    n2 = CalcVectorCrossProduct3D(r3, r4)

    n1Normalized = NormalizeVector3D(n1)
    n2Normalized = NormalizeVector3D(n2)

    dot_n1n2 = CalcVectorDotProduct3D(n1Normalized, n2Normalized)

    if CalcVectorDotProduct3D(CalcVectorCrossProduct3D(n1, n2), r3) > 0:
        return -math.degrees(math.acos(dot_n1n2))

    return math.degrees(math.acos(dot_n1n2))


def CalcVectorDotProduct3D(vector1, vector2):
    result = vector1[0]*vector2[0]+vector1[1]*vector2[1]+vector1[2]*vector2[2]
    return result


def CalcVectorCrossProduct3D(vector1, vector2):
    newVector = (vector1[1] * vector2[2] - vector1[2] * vector2[1],
                 vector1[2] * vector2[0] - vector1[0] * vector2[2],
                 vector1[0] * vector2[1] - vector1[1] * vector2[0])
    return newVector


def SubtractVectors3D(vector1, vector2):
    newVector = (vector1[0] - vector2[0], vector1[1] - vector2[1], vector1[2] - vector2[2])
    return newVector


def NormalizeVector3D(vector):
    magnitude = math.sqrt(vector[0] ** 2 + vector[1] ** 2 + vector[2] ** 2)
    newVector = (vector[0] / magnitude, vector[1] / magnitude, vector[2] / magnitude)
    return newVector


def ReadProteinFromPDB(fileName):
    protein = []
    chain = []
    backbone = {}

    pdbFile = open(fileName)
    for line in pdbFile:
        if line[:3] == "TER":
            # Chain break
            if chain:
                protein.append(chain)
                chain = []
            continue
        if line[:4] == "ATOM":
            atomName = line[12:16].strip()
            if atomName in ["N", "CA", "C"]:

                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])

                backbone[atomName] = (x, y, z)

                if len(backbone) == 3:
                    aminoAcid = [backbone["N"], backbone["CA"], backbone["C"]]
                    chain.append(aminoAcid)
                    backbone = {}
    if len(chain):
        protein.append(chain)
    return protein


if __name__ == '__main__':
    main()


