reset()

from sage.crypto.boolean_function import BooleanFunction
import sys
load("aux_function.sage")
load("MILP_function.sage")

def computeInequalitiesMC(M):
	#Precompute inequalities for MC
	M_nrows = M.nrows()
	M_ncols = M.ncols()

	V = vectorF2n(M_nrows)
	L = [M*vector(x) for x in V]
	MCbox = [None for _ in range(1 << M_nrows)]
	for i in range(1 << M_nrows):
		x = 0
		for j in range(M_nrows):
			if L[i][j] == 1:
				x += 2^j
		MCbox[i] = x

	(BPR,anf) = anfSbox(MCbox)
	TMC = divPropANFBinTable(anf)
	return sboxReducedInequalities(TMC)

def applySbox(m, xr, yr, ineq, sboxSize, nbSbox):
	#Add the constraints for the Sbox layer from xr to yr, using the inequalities in ineq
	#Constraints for Sbox

	for j in range(nbSbox):
		inputvar = [xr[sboxSize*j+i] for i in range(sboxSize)]
		outputvar = [yr[sboxSize*j+i] for i in range(sboxSize)]
		addSboxConstr(m, ineq, inputvar, outputvar)

def applyMCAsLbox(m, wr, xrp1, ineqMC, M_ncols, nbCol, sboxSize):
	colSize = M_ncols*sboxSize
	for col in range(nbCol):
		for offset in range(sboxSize):
			inputvar = [wr[col*colSize + offset + j*sboxSize] for j in range(M_ncols)]
			outputvar = [xrp1[col*colSize + offset + j*sboxSize] for j in range(M_ncols)]
			addSboxConstr(m,ineqMC,inputvar,outputvar)

def applyMCAsCopyXor(m, wr, xrp1, M, nbCol, r):
	M_nrows = M.nrows()
	M_ncols = M.ncols()
	for col in range(nbCol):
		#First, create the temporary variables for the Copy + XOR modeling
		t = [[0 for __ in range(M_ncols)] for _ in range(M_nrows)]
		for i in range(M_nrows):
			for j in range(M_ncols):
				if M[i][j] == 1:
					t[i][j] = m.addVar(vtype=GRB.BINARY, name="t_"+str(r)+"_"+str(col)+"_"+str(i)+"_"+str(j))

		#Copy constraints
		for j in range(M_ncols):
			m.addGenConstrOr(wr[col*M_ncols + j], [t[i][j] for i in range(M_nrows) if M[i][j] == 1])

		#XOR constraints
		for i in range(M_nrows):
			m.addConstr(xrp1[col*M_ncols + i] == quicksum(t[i][j] for j in range(M_ncols) if M[i][j] == 1))


rMax = 13

S    = [0xc,0x5,0x6,0xb,0x9,0x0,0xa,0xd,0x3,0xe,0xf,0x8,0x4,0x7,0x1,0x2] #normal Sbox
Sin  = [0x8,0x5,0x0,0x2,0xa,0x9,0x4,0xc,0xe,0x7,0xb,0xd,0x6,0x1,0xf,0x3] #Sbox for first round
Sout = [0x6,0x7,0xe,0x3,0x1,0x0,0x8,0xd,0x9,0x4,0xf,0xa,0xc,0x5,0xb,0x2] #sbox for last round

P = [0,4,8,12,1,5,9,13,2,6,10,14,3,7,11,15]
M = matrix(GF(2),[
    [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,],
    [0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,],
    [0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,],
    [0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,],
    [0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,],
    [0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,],
    [0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,],
    [0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,],
    [0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,],
    [0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,],
    [0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,],
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,],
    [0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,],
    [0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,],
    [0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,],
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,],
    ])
name = "PRESENT_13r"
linAsSbox = False
noLastMC = True


#Precompute inequalities for the sboxes
(BPR,anf) = anfSbox(S)
T = divPropANFBinTable(anf)
ineqSbox = sboxReducedInequalities(T)
(BPR,anfin) = anfSbox(Sin)
Tin = divPropANFBinTable(anfin)
ineqSboxin = sboxReducedInequalities(Tin)
(BPR,anfout) = anfSbox(Sout)
Tout = divPropANFBinTable(anfout)
ineqSboxout = sboxReducedInequalities(Tout)

M_nrows = M.nrows()
M_ncols = M.ncols()
if linAsSbox:
	#Precompute inequalities for MC
	ineqMC = computeInequalitiesMC(M)

#Some useful constants
sboxSize = anf[0].parent().n_variables() #size of the Sbox
blockSize = sboxSize*len(P) #block size
nbSbox = len(P) #number of sboxes
nbCol = blockSize//M_nrows #Number of column in the state if linAsSbox = false
if linAsSbox:
	nbCol = nbSbox//M_nrows
colSize = M_ncols*sboxSize

#Compute the permutation extended at the bit level
expandP = [0 for i in range(blockSize)]
for isbox in range(nbSbox):
	for j in range(sboxSize):
		expandP[isbox*sboxSize+j] = P[isbox]*sboxSize+j

#Create the models and the variables
m = Model(name)

if noLastMC:
	#Variables are x[r] --S--> y[r] --P--> z[r] --MC--> w[r] --ARK--> x[r+1]
	#No MC on last round, stops at y[rMax-1]
	x = [[m.addVar(vtype=GRB.BINARY, name="x_"+str(i)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax)]
	y = [[m.addVar(vtype=GRB.BINARY, name="y_"+str(i)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax)]
	z = [[m.addVar(vtype=GRB.BINARY, name="z_"+str(i)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax-1)]
	w = [[m.addVar(vtype=GRB.BINARY, name="w_"+str(i)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax-1)]
	k = [[m.addVar(vtype=GRB.BINARY, name="k_"+str(i)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax-1)]
else:
	#Variables are x[r] --S--> y[r] --P--> z[r] --MC--> w[r] --ARK--> x[r+1]
	x = [[m.addVar(vtype=GRB.BINARY, name="x_"+str(i)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax+1)]
	y = [[m.addVar(vtype=GRB.BINARY, name="y_"+str(i)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax)]
	z = [[m.addVar(vtype=GRB.BINARY, name="z_"+str(i)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax)]
	w = [[m.addVar(vtype=GRB.BINARY, name="w_"+str(i)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax)]
	k = [[m.addVar(vtype=GRB.BINARY, name="k_"+str(i)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax)]

m.update()

#Round Function constraints (except last round)
for r in range(rMax-1):
	#Constraints for Sbox
	if r == 0:
		applySbox(m, x[r], y[r], ineqSboxin, sboxSize, nbSbox)
	else:
		applySbox(m, x[r], y[r], ineqSbox, sboxSize, nbSbox)

	#Constraints for Permutation
	for i in range(blockSize):
		m.addConstr(z[r][expandP[i]] == y[r][i])

	#Constraints for MixColumn
	if linAsSbox:
		applyMCAsLbox(m, z[r], w[r], ineqMC, M_ncols, nbCol, sboxSize)

	else:
		applyMCAsCopyXor(m, z[r], w[r], M, nbCol, r)

	#Constraints for ARK
	for i in range(blockSize):
		m.addConstr(w[r][i] + k[r][i] == x[r+1][i])

#last round
r = rMax-1
#Constraints for Sbox
applySbox(m, x[r], y[r], ineqSboxout, sboxSize, nbSbox)

m.update()

m.write(name+".mps")
