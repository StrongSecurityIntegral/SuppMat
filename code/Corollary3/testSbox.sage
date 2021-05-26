def greater(a,b):
	#return True if a[i] >= b[i] for all i
	#False otherwise
	for i in range(len(a)):
		if a[i] < b[i]:
			return False
	return True

def vecToInt(x):
	vx = 0
	for i in range(len(x)):
		if x[i] == 1:
			vx += (1 << i)

	return vx
 
def custom_ANF(B,SBOX,nInput,nOutput):
	"""
	Compute the ANF of S in the ring B
	Quick and dirty, just because of a memory leak from Sage... (Ticket 21892)
	Should be much better with a Moebius transform, but works fine for small sboxes
	"""

	#Get the bit representation of each value
	SBOX_bits = [Integer(x).bits() + [0 for i in range(nOutput-len(Integer(x).bits()))] for x in SBOX]
	xvar = B.gens()[:nInput]

	ANF = [B(0) for _ in range(nOutput)]
	F2n = [tuple(Integer(c).bits() + [0 for i in range(nInput-Integer(c).nbits())]) for c in range(1 << nInput)]

	for i in range(nOutput):
		#For each component of the ANF
		for u in F2n:
			#compute a_u for ANF[i]
			a_u = 0
			for x in F2n:
				if greater(u,x):
					vx = vecToInt(x)
					a_u ^^= SBOX_bits[vx][i]
			if a_u == 1:
				xu = B(prod(xvar[j] for j in range(nInput) if u[j] == 1))
				ANF[i] += xu

	return ANF

def computeANFWithKey(S, BPR, nInput, nOutput):
	"""
	Compute the ANF of S(x+k), S being an Sbox in nInput variables with nOutput bits
	BPR must be a BooleanPolynomialRing in 2*nInput variables to handle a memory leak from Sage
	"""

	#ANF of S
	anfS = custom_ANF(BPR,S,nInput,nOutput)
	x = BPR.gens()[:nInput]
	k = BPR.gens()[nInput:]

	#ANF of S(x+k)
	d = {x[i] : x[i] + k[i] for i in range(nInput)}
	anf = [BPR(anfS[i].subs(d)) for i in range(nOutput)]

	return anf

def getKeyPolynomials(f,mapMonToXK,mapMonToExp,n):
	"""
	Return a list L such that L[u] contains the polynomials p_u(k) form f
	"""

	L = [f.parent()(0) for u in range(1 << n)]
	for mon in f.monomials():
		xu,kv = mapMonToXK[mon]
		L[mapMonToExp[xu]] += kv

	return L

def is_in_span(M,v):
	"""
	Return True is v is the the span of the rows of M
	"""
	try:
		M.solve_left(v)
		return True
	except ValueError:
		return False

def preserveDegree(f, mapMonToXK,mapMonToExp,n):
	"""
	Given a boolean function f (i.e. one cordinate function of S(x+k)), return a vector v such that v[i] = True if f preserves the degree of the i-th input
	"""

	#compute the p_u for f
	p = getKeyPolynomials(f, mapMonToXK,mapMonToExp,n)

	#compute the matrix M s.t. M[u][v] = lamda_u,v, i.e. M[u][v] = 1 iif k^v belongs to p_u
	M = matrix(GF(2), 1 << n, 1 << n)
	for u in range(1 << n):
		for mon in p[u].monomials():
			M[u,mapMonToExp[mon]] = 1

	v = [None for i in range(n)]
	for i in range(n):
		#Matrix without the row foe e_i and row 0
		T = M.matrix_from_rows_and_columns([j for j in range(1 << n) if j != (1 << i) and j != 0], range(1 << n))
		v[i] = not is_in_span(T,vector(M[(1 << i)]))

	return v

def checkPreserveAlgDegree(S,BPR,mapMonToXK,mapMonToExp,nInput,nOutput):
	anf = computeANFWithKey(S,BPR,nInput,nOutput)

	A = [[None for i in range(nInput)] for j in range(nOutput)]
	#A[j][i] = True means that the coordinate function S_j preserves the degree of the i-th input

	for j in range(nOutput):
		A[j] = preserveDegree(anf[j],mapMonToXK,mapMonToExp,nInput)

	#Check if S preserves the algebraic degree
	#This means that for all i, there is a coordinate preserving the degree of the i-th input
	#i.e. for all i, there exists j such that A[j][i] = True
	check = [None for i in range(nInput)] #check[i] = True if degree of input i is preserved

	for i in range(nInput):
		checki = False
		for j in range(nOutput):
			if A[j][i]: #there is a degree preserving coordinate for input i
				check[i] = True
				checki = True
				break
		if not checki: #no degree preserving coordinate for input i
			# print("No degree preserving coordinate for input " + str(i))
			check[i] = False

	return check

def checkDegreeListSbox(listNameSbox, nInput,nOutput, printCounter=False):
	"""
	Given a list #listNameSbox of Sboxes with an input/output size of #nInput and #nOutput
	"""


	BPR = BooleanPolynomialRing(2*nInput,["x"+str(i) for i in range(nInput)] + ["k"+str(i) for i in range(nInput)])
	x = BPR.gens()[:nInput]
	k = BPR.gens()[nInput:]

	#Sage multivariate polynomials are not the best to check which monomials appears
	#Prepare some precomputed maps that will be used to facilitate that

	#mapMonToXK[m] = [m_x, m_k], with m = m_x*m_k, m containing only x variables, m_k containing only k variables
	#mapMonToExp[m] = u, where m is either a monomial in x or k, u is the interget such that m = x^u (or m = k^u)
	mapMonToXK = dict()
	mapMonToExp = dict()
	for u in range(1 << nInput):
		xu = BPR(prod(x[i] for i in range(nInput) if ((u >> i)&1) == 1))
		ku = BPR(prod(k[i] for i in range(nInput) if ((u >> i)&1) == 1))
		mapMonToExp[xu] = u
		mapMonToExp[ku] = u

		for v in range(1 << nInput):
			kv = BPR(prod(k[i] for i in range(nInput) if ((v >> i)&1) == 1))
			mapMonToXK[xu*kv] = [xu,kv]

	#Counter for each case
	ctr = dict()

	for name, S in listNameSbox:
		checkAlgDegree = checkPreserveAlgDegree(S,BPR,mapMonToXK,mapMonToExp,nInput,nOutput)

		checkAllInput = True
		for c in checkAlgDegree:
			if not c:
				checkAllInput = False
				break

		print(name + " : " + str(checkAlgDegree) + ", preserves all input : " + str(checkAllInput))
		# print(name + " algDegree: " + str(checkAlgDegree), flush=True)
		if str(checkAlgDegree) in ctr:
			ctr[str(checkAlgDegree)] += 1
		else:
			ctr[str(checkAlgDegree)] = 1

	if printCounter:
		for k in ctr:
			print(str(k) + " : " + str(ctr[k]))


# # ---- Check all 302 4-bit affine equivalence classes ----
# from sage.crypto.sboxes import affine_equiv_classes
# listSbox = affine_equiv_classes[4]
# listNameSbox = []
# for S in listSbox:
# 	hexrep = ""
# 	for tmp in S:
# 		hexrep += hex(tmp).lstrip("0x")
# 	listNameSbox.append([hexrep, S])
# n = 4
# checkDegreeListSbox(listNameSbox,n,n,printCounter=True)


# #---- Check for all Sboxes in the sboxes module ----
from sage.crypto.sboxes import sboxes
for n in range(3,9):
	print("Sboxes of size " + str(n))
	listNameSbox = []
	for name in sboxes:
		if sboxes[name].input_size() == n and sboxes[name].output_size() == n:
			listNameSbox.append([name,sboxes[name]])
	checkDegreeListSbox(listNameSbox,n,n,printCounter=True)	


# #---- Check for all inverse mappings ----
# for n in range(3,9):
# 	# print("Inversion " + str(n)+"-bit")
# 	#Compute the Sbox of x^-1 in GF(2^n)
# 	F = GF(2^n, "alpha", modulus="primitive", repr="int")
# 	a = F.gen()
# 	S = [None for _ in range(1 << n)]
# 	S[0] = 0

# 	for x in range(1,1 << n):
# 		e = F(0)
# 		for i in range(n):
# 			if ((x >> i)&1) == 1:
# 				e += a^i
# 		inve = e^(-1)
# 		S[x] = inve.integer_representation()

# 	listNameSbox = [["inversion"+str(n), S]]
# 	checkDegreeListSbox(listNameSbox,n,n)	


# #---- For a single Sbox ----
# from sage.crypto.sboxes import sboxes
# name = "PRESENT"
# listNameSbox = [[name, sboxes[name]]]
# n = sboxes[name].input_size()
# checkDegreeListSbox(listNameSbox,n,n)

# #---- Check all 5-bit boolean functions (up to equivalence)----
# #For all balanced 5-bit functions, replace the file by F5b.txt
# #For all balanced 6-bit functions, replace the file by F6b.txt and nInput = 6

# with open("F5.txt","r") as inputfile:
# 	nInput = 5
# 	nOutput = 1
# 	listNameSbox = []
# 	for l in inputfile:
# 		s = l.rstrip("\n")
# 		S = [int(x) for x in s]
# 		listNameSbox.append([s,S])
# 	checkDegreeListSbox(listNameSbox,nInput,nOutput,printCounter=True)