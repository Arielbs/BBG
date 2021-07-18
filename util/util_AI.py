
import numpy as np
from scipy.spatial.transform import Rotation as R


#### Creating 4letters words
# Due to the 3.8 Ang constrains each word can be described by 3 angular parameters: theta3, theta4, phi4. 
# class object W4L includes all the functions to move between euclidian and angular description, as well as 
# obtaining a word state (position and orientation). 

def XYZ2tf():
    return 1
def tf2XYZ():
    return 1
def buildLetter4DB():
    return 1
def w4l_canonical(w4l):
    '''
    returns the canonical euclidian word description 
    input: angular description (list(len==3)), euclidian minimal (np.array shape(5,)) always canonical
           euclidian (np.array((4,3))) 
           returns np.array at the form of ([0,0,0],[0,0,R],[0,Y3,Z3],[X4,Y4,Z4])
           or the minimal form [Y3,Z3,X4,Y4,Z4]
    method: get the vector v12 (l1 to l2), align np.array with v12 to [0,0,1]
    '''
    return 1

def wordXYZ2angle():
    '''
    Takes the XYZ coordinate of a work np.array(4X3) 
    returns minimal word description (theta3,theta4,phi4)
    '''
    return 1
def wordAngle2XYZ(ang: list ,R: float=3.8, description: str='2d'):
    '''
    Takes word angular description (theta3,theta4,phi4) 
    returns the XYZ coordinate of a work np.array(4X3) in its canonical form if full description is requireed
    returns flatten version of the last 5 positions [Y3,Z3,X4,Y4,Z4]
    
    '''
    theta3 = ang[0] ; theta4 = ang[1] ; phi4 = ang[2] 
    w4l = np.zeros((4,3)) ; w4l[1,2] = R
    w4l[2,1]=w4l[1,1] + R*np.sin(theta3) ;  w4l[2,2]=w4l[1,2] + R*np.cos(theta3)
    w4l[3,0]=w4l[2,0] + R*np.sin(theta4)*np.cos(phi4) ; w4l[3,1]=w4l[2,1] + R*np.sin(theta4)*np.sin(phi4) ; w4l[3,2]=w4l[2,2] + R*np.cos(theta4)
    if description=='2d':
        return w4l
    if description=='1d':
        return w4l[2:].flatten()[1:] 


 ###### functions to build equaly distributed points in spheres
def get_letter3(d0,R=3.8,XYZ=[0,0,3.8]):
    '''
    This function resturns a NX3 size numpy arraya N number of points as determined by d0
    do is the distance between points on the curve, R radi of interaction set by default to 3.8 Ang
    XYZ is the coordinate of letter #2 (default at [0,0,R])
    returns: 
    theta3array NX3 np.arrat [theta3, Z0,Y0] note theta is the angle on the ZX plane (by arbitrary selection)
    '''
    thetaRange = [0,2/3*np.pi]  # range in angles [0,120] any value larger the 120O would violate the min distance restriction 
    d_theta = 2*np.arcsin(d0/(2*R))  # distance between each two points
    thetaList = np.expand_dims(np.arange(thetaRange[0],thetaRange[1],d_theta),axis=1)
    theta3array = np.concatenate((thetaList, R*np.cos(thetaList)+XYZ[2], R*np.sin(thetaList)+XYZ[1] ), axis=1)
    ### to test use np.linalg.norm(ZYtheta[k,1:]-ZYtheta[k+1,1:])
    return theta3array



def fibonacci_sphere(num_points: int, rad=3.8):
    ga = (3 - np.sqrt(5)) * np.pi # golden angle                                                                             

    # Create a list of golden angle increments along tha range of number of points                                           
    theta = ga * np.arange(num_points)

    # Z is a split into a range of -1 to 1 in order to create a unit circle                                                  
    z = np.linspace(0.5/num_points-1, 1-1/num_points, num_points)
    # a list of the radii at each height step of the unit circle                                                             
    radius = np.sqrt(1 - z * z)
    # Determine where xy fall on the sphere, given the azimuthal and polar angles                                            
    y = radius * np.sin(theta)
    x = radius * np.cos(theta)
    XYZ = np.array((x,y,z)).T
    return XYZ*rad


def get_fibonacci_spacing(n: int,R: float):
    XYZ = fibonacci_sphere(n,R)
    binN = int(n/2)
    distMap = sc.spatial.distance.cdist(XYZ, XYZ, metric='euclidean')
    flatDistMap = distMap.flatten()
    counts, bin_edges = np.histogram(flatDistMap, bins=binN)
    countsLocList = [N for N,x in enumerate(counts) if x>(n+2)]
    return bin_edges[countsLocList[0]]


def get_fibo_n(target_d: float,n=100,R=3.8):
    '''
    This function should iteratively find thw number of points needed to be generated to obtain a spherically 
    evenly distributed points (sorry this is brute force and could better done analytically)
    Input: target_d: target distance [float], nstart [int] = number of points to start the search from, R - sphere radius
    Output: nSelect [int]
    '''
    fibo_d = get_fibonacci_spacing(n,R)  # distance between two close points on the sphere 
    ratio = fibo_d/target_d  # if >1 than more points should be created 
    while np.abs(ratio-1)>0.05:
        n = int(n*ratio*ratio)
        print("n",n)
        fibo_d = get_fibonacci_spacing(n,R)
        print("fibo_d: ",fibo_d)
        ratio = fibo_d/target_d
        print("ratio: ",ratio)
    return n, fibo_d






def getAllPDBsID():
	"""Return a list of all PDB entries currently in the RCSB Protein Data Bank
	Returns
	-------
	out : list of str
	    A list of all of the PDB IDs currently in the RCSB PDB
	Examples
	--------
	>>> print(get_all()[:10])
	['100D', '101D', '101M', '102D', '102L', '102M', '103D', '103L', '103M', '104D']
	"""
	url = 'http://www.rcsb.org/pdb/rest/getCurrent'
	response = requests.get(url)
	if response.status_code == 200:
		pass
	else:
		warnings.warn("Retrieval failed, returning None")
		return
	result  = str(response.text)
	return re.findall(r'"(.*?)"', result)




def pdb2alpha(p):
	'''
	Generate an alternative pdb file which includes info only on the Calpha
	get path to pdb [str] and return DB [Dataframe] of the ordered Calphas of chain A and only if it's well ordered 
	'''
	pF = open(p,'r')
	lines = [line.split() for line  in pF.readlines() if (line[:4]=="ATOM" and line[13:15]=="CA")]
	pF.close()
	listOfLists = [[line[5],line[6],line[7],line[8],line[3]] for line in lines if (line[-3]=="1.00" and line[4]=="A") ]
	return pd.DataFrame(listOfLists,columns=["n","x","y","z","res"])



def checkConsecutive(l):
	'''
	naive method to make sure residue sequnces is concecutive.
	get list of integers and check if all are concecutive 
	'''
	if np.sum([0 if l2-l1==1 else 1  for l1,l2 in zip(l[:-1],l[1:])])==0:
		return 1  # list is good 
	else:
		return 0 # don't use this list and continue to next step 



def get_rotation_matrix(i_v, unit=None):
	'''# returns a rotation matrix between 3D vector i_v and "unit" vector '''
	# From http://www.j3d.org/matrix_faq/matrfaq_latest.html#Q38
	if unit is None:
		unit = [1.0, 0.0, 0.0]
	# Normalize vector length
	i_v /= np.linalg.norm(i_v)
	# Get axis
	uvw = np.cross(i_v, unit)
	# compute trig values - no need to go through arccos and back
	rcos = np.dot(i_v, unit)
	rsin = np.linalg.norm(uvw)
	#normalize and unpack axis
	if not np.isclose(rsin, 0):
		uvw /= rsin
	u, v, w = uvw
	# Compute rotation matrix - re-expressed to show structure
	return (
		rcos * np.eye(3) +
		rsin * np.array([
			[ 0, -w,  v],
			[ w,  0, -u],
			[-v,  u,  0]
		]) +
		(1.0 - rcos) * uvw[:,None] * uvw[None,:]
	)


def rotateMat(rotMat,vecS):
	'''
	recieves: rotMat - a calid rotation matrix, and an np.array(N,3)
	returns np.array(N,3) after rotation operation around the [0,0,0] point
	'''
	newVec = np.zeros((vecS.shape))
	for n,ele in enumerate(vecS):
		if n >0:
			newVec[n,:] =  np.dot(rotMat,ele)
	return newVec # this is the rotated set of arrays 


def resDist(A):
	'''
	recieves np.array(N,3) representing N (x,y,z) positions 
	returns np.array(N-1) for the distances between sequential set of positions 
	'''
	distVec = np.zeros((A.shape[0]-1))
	for n,ele in enumerate(A):
		if n >0:
			distVec[n-1] = np.linalg.norm(ele-A[n-1,:])  
	return distVec


def orientQuadruplets(A,check=[0,0]):
	'''
	# receives numpy array of size (4,3) representing 4 sequntial C alpha x,y,z positions 
	# returns a canonical form of the positions Quadruplets in the form of:
	# np.array([0,0,0],[x1,0,0],[x2,y2,0],[x3,y3,z3])
	'''
	if check[0]:
		distVec = resDist(A)
		if (max(distVec)>4 or np.min(distVec)<3.6):
			print("Distances within the original date is not compatible with C-C bond", distVec)
			return np.zeros((4,3))
	#### begin transformations 
	A0 = A-A[0]   # translating atom 1 to [0,0,0]
	val = np.linalg.norm(A0[1])
	rotMat = get_rotation_matrix(A0[1], unit=None) ; A0[1] = A0[1]*val  # geenrate rotation map to align second atom with axis X
	A1 = rotateMat(rotMat,A0)  # align the  Quadruplets such that the second atom is on the x axis        
	### calc. angle between the third atom and axis y
	cosine_angle = np.dot(A1[2,1:], [1.,0.]) / (np.linalg.norm(A1[2,1:]) * np.linalg.norm([1.,0.]))
	angle = np.arccos(cosine_angle) # get the angle
	rb = R.from_euler('xyz', [-angle,0,0], degrees=False) # rotation matrix to rotate around the x axis   
	A2 = rotateMat(rb.as_dcm(),A1)  #   align the  Quadruplets such that the third atom is on the y axis 
	if check[1]:
		distVec = resDist(A2)
		if (max(distVec)>4 or np.min(distVec)<3.6):
			print("Distances within the result data is not compatible with C-C bond", distVec)
			return np.zeros((4,3))
	return A2


def calc_angle(u1, u2, u3):
	""" Calculate angle between three points.
	The angle is in [-pi, pi] and rad.
	"""
	V1 = u1-u2 ; V2 = u3-u2 
	cosine_angle1 = np.dot(V1,V2) / (np.linalg.norm(V1) * np.linalg.norm(V2))
	return np.arccos(cosine_angle1) # get the psi ()


def calc_dihedral(u1, u2, u3, u4):
	""" Calculate dihedral angle method. From bioPython.PDB
	(adapted to np.array)
	Calculate the dihedral angle between 4 vectors
	representing 4 connected points. The angle is in
	[-pi, pi].
	"""
	a1 = u2 - u1
	a2 = u3 - u2
	a3 = u4 - u3
	
	v1 = np.cross(a1, a2)
	v1 = v1 / (v1 * v1).sum(-1)**0.5
	v2 = np.cross(a2, a3)
	v2 = v2 / (v2 * v2).sum(-1)**0.5
	porm = np.sign((v1 * a3).sum(-1))
	rad = np.arccos((v1*v2).sum(-1) / ((v1**2).sum(-1) * (v2**2).sum(-1))**0.5)
	if not porm == 0:
		rad = rad * porm
	return rad   


def quadruplets_Cartesian2Angles(A):
	'''
	Recieves A - np.array((12)) of 4 atom positions 
	Returns np.array((3)) [theta1,theta2,dihedral], three angles for the three degrees of freedom in the system
	
	psi is the dihedral angle on the Calpha-C bond. 
	phi is the dihedral angle on the C-N bond
	Here we work only with the Calphas of four residues and the representation has 
	to be based on 2 angle and 1 dihedral angle:
	theta1 - angle between first three atoms, 
	Theta2 - angle between last three atoms
	dihedral - torsion angle between the second and third atom
	'''
	A = A.reshape((4,3))
	output = np.zeros((3))
	output[0] = calc_angle(A[0],A[1],A[2])         # theta1
	output[1] = calc_angle(A[1],A[2],A[3])         # theta2
	output[2] = calc_dihedral(A[0],A[1],A[2],A[3]) # dihedral
	return output








