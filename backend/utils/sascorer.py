#
# calculation of synthetic accessibility score as described in:
#
# Estimation of Synthetic Accessibility Score of Drug-like Molecules based on Molecular Complexity and Fragment Contributions
# Peter Ertl and Ansgar Schuffenhauer
# Journal of Cheminformatics 1:8 (2009)
# http://www.jcheminf.com/content/1/1/8
#

import math
import pickle
import gzip
import os
import os.path as op
from rdkit import Chem, RDConfig
from rdkit.Chem import rdFingerprintGenerator, rdMolDescriptors

_fscores = None
mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2)


def readFragmentScores(name=None):
    global _fscores
    if name is None:
        local_path = op.join(op.dirname(__file__), 'fpscores.pkl.gz')
        if op.exists(local_path):
            name = local_path
        elif hasattr(RDConfig, 'RDContribDir'):
            contrib_path = op.join(RDConfig.RDContribDir, 'SA_Score', 'fpscores.pkl.gz')
            if op.exists(contrib_path):
                name = contrib_path
            else:
                name = local_path
        else:
            name = local_path

    with gzip.open(name, 'rb') as f:
        data = pickle.load(f)
    outDict = {}
    for i in data:
        for j in range(1, len(i)):
            outDict[i[j]] = float(i[0])
    _fscores = outDict


def numBridgeheadsAndSpiro(mol, ri=None):
    nSpiro = rdMolDescriptors.CalcNumSpiroAtoms(mol)
    nBridgehead = rdMolDescriptors.CalcNumBridgeheadAtoms(mol)
    return nBridgehead, nSpiro


def calculateScore(m):
    if m is None or not m.GetNumAtoms():
        return None

    if _fscores is None:
        readFragmentScores()

    # Fragment score
    sfp = mfpgen.GetSparseCountFingerprint(m)

    score1 = 0.0
    nf = 0
    nze = sfp.GetNonzeroElements()
    for fid, count in nze.items():
        nf += count
        score1 += _fscores.get(fid, -4) * count

    if nf == 0:
        return None
    score1 /= nf

    # Features score
    nAtoms = m.GetNumAtoms()
    nChiralCenters = len(Chem.FindMolChiralCenters(m, includeUnassigned=True))
    ri = m.GetRingInfo()
    nBridgeheads, nSpiro = numBridgeheadsAndSpiro(m, ri)
    nMacrocycles = 0
    for x in ri.AtomRings():
        if len(x) > 8:
            nMacrocycles += 1

    sizePenalty = nAtoms**1.005 - nAtoms
    stereoPenalty = math.log10(nChiralCenters + 1)
    spiroPenalty = math.log10(nSpiro + 1)
    bridgePenalty = math.log10(nBridgeheads + 1)
    macrocyclePenalty = 0.0
    if nMacrocycles > 0:
        macrocyclePenalty = math.log10(2)

    score2 = 0.0 - sizePenalty - stereoPenalty - spiroPenalty - bridgePenalty - macrocyclePenalty

    score3 = 0.0
    numBits = len(nze)
    if nAtoms > numBits and numBits > 0:
        score3 = math.log(float(nAtoms) / numBits) * 0.5

    sascore = score1 + score2 + score3

    # Transform raw value into scale between 1 and 10
    min_val = -4.0
    max_val = 2.5
    sascore = 11.0 - (sascore - min_val + 1) / (max_val - min_val) * 9.0

    if sascore > 8.0:
        try:
            sascore = 8.0 + math.log(sascore + 1.0 - 9.0)
        except ValueError:
            pass
    if sascore > 10.0:
        sascore = 10.0
    elif sascore < 1.0:
        sascore = 1.0

    return sascore
