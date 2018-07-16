#Converting below from C++ to Python
"""NOTES:
My comments all start with ì#TO DO:î to differentiate between my comments and the existing comments from the original program.

Overall to do:
Replacing the pointers to arrays (in Python arrays are mutable)

Including the header file (The header file looks like it was used to pre-define the variable types)

ifstream/ofstream input/output equivalents in Python (there are some writeinput...etc. not sure what they are for?)

Meshout and other output formatting (Still trying to find the equivalent in Python )
"""
#PROGRAM STARTS:





"""***************** 
* CANOE CREATOR *
*****************

Version 1.02

Created for the University of Toronto Concrete Canoe Team by David Ruggiero

""" 

import math

#TO DO: need to include header Creator.h?, not sure how to do this

counter = 0
bowpower = 0			# these shouldn't be global variables but I'm tired and it's 2:15 AM
sternpower = 0

#TO DO: predef functions:
def Spline( A, B, C, D, t): #TO DO: done
	return (A*t*t*t + B*t*t + C*t + D)

def SplineArea( P0y, P0z, P1y, P1z, P2y, P2z, P3y, P3z, ti, tf): #TO DO: done
	A1 = P0y - 3*P1y + 3*P2y - P3y
	B1 = 3*P1y - 6*P2y + 3*P3y
	C1 = 3*P2y - 3*P3y;
	D1 = P3y;

	A2 = 0;
	B2 = 3*P0z - 9*P1z + 9*P2z - 3*P3z;
	C2 = 6*P1z - 12*P2z + 6*P3z;
	D2 = 3*P2z - 3*P3z;

	E = [0]
	E = E*7
	E[0] = D1*D2
	E[1] = C1*D2 + D1*C2
	E[2] = B1*D2 + C1*C2 + D1*B2
	E[3] = A1*D2 + B1*C2 + C1*B2 + D1*A2
	E[4] = A1*C2 + B1*B2 + C1*A2
	E[5] = A1*B2 + B1*A2
	E[6] = A1*A2
	
	area = 0
	for i in range (7):
		area += E[i]/(i+1) * (pow(tf,(i+1)) - pow(ti,(i+1)))
		
	return area 

def SplineMoment(P0y, P0z, P1y, P1z, P2y, P2z, P3y, P3z, ti, tf): #TO DO: done

	A1 = P0y - 3*P1y + 3*P2y - P3y
	B1 = 3*P1y - 6*P2y + 3*P3y
	C1 = 3*P2y - 3*P3y
	D1 = P3y

	A2 = 0
	B2 = 3*P0z - 9*P1z + 9*P2z - 3*P3z
	C2 = 6*P1z - 12*P2z + 6*P3z
	D2 = 3*P2z - 3*P3z

	A3 = P0z - 3*P1z + 3*P2z - P3z
	B3 = 3*P1z - 6*P2z + 3*P3z
	C3 = 3*P2z - 3*P3z
	D3 = P3z

	A4 = 0
	B4 = 3*P0y - 9*P1y + 9*P2y - 3*P3y
	C4 = 6*P1y - 12*P2y + 6*P3y
	D4 = 3*P2y - 3*P3y

	E = [0]
	E = E*7
	E[6] = A1*A2
	E[5] = A1*B2 + B1*A2
	E[4] = A1*C2 + B1*B2 + C1*A2
	E[3] = A1*D2 + B1*C2 + C1*B2 + D1*A2
	E[2] = B1*D2 + C1*C2 + D1*B2
	E[1] = C1*D2 + D1*C2
	E[0] = D1*D2

	L = [0]
	L = L*10
	L[9] = A3*E[6]
	L[8] = B3*E[6] + A3*E[5]
	L[7] = C3*E[6] + B3*E[5] + A3*E[4]
	L[6] = D3*E[6] + C3*E[5] + B3*E[4] + A3*E[3]
	L[5] = D3*E[5] + C3*E[4] + B3*E[3] + A3*E[2]
	L[4] = D3*E[4] + C3*E[3] + B3*E[2] + A3*E[1]
	L[3] = D3*E[3] + C3*E[2] + B3*E[1] + A3*E[0]
	L[2] = D3*E[2] + C3*E[1] + B3*E[0]
	L[1] = D3*E[1] + C3*E[0]
	L[0] = D3*E[0]

	Ey = [0]
	Ey = Ey*7
	Ey[6] = A3*A4
	Ey[5] = A3*B4 + B3*A4
	Ey[4] = A3*C4 + B3*B4 + C3*A4
	Ey[3] = A3*D4 + B3*C4 + C3*B4 + D3*A4
	Ey[2] = B3*D4 + C3*C4 + D3*B4
	Ey[1] = C3*D4 + D3*C4
	Ey[0] = D3*D4

	Ly = [0]
	Ly = Ly*10
	Ly[9] = A1*Ey[6]
	Ly[8] = B1*Ey[6] + A1*Ey[5]
	Ly[7] = C1*Ey[6] + B1*Ey[5] + A1*Ey[4]
	Ly[6] = D1*Ey[6] + C1*Ey[5] + B1*Ey[4] + A1*Ey[3]
	Ly[5] = D1*Ey[5] + C1*Ey[4] + B1*Ey[3] + A1*Ey[2]
	Ly[4] = D1*Ey[4] + C1*Ey[3] + B1*Ey[2] + A1*Ey[1]
	Ly[3] = D1*Ey[3] + C1*Ey[2] + B1*Ey[1] + A1*Ey[0]
	Ly[2] = D1*Ey[2] + C1*Ey[1] + B1*Ey[0]
	Ly[1] = D1*Ey[1] + C1*Ey[0]
	Ly[0] = D1*Ey[0]

	momentareaz = 0
	momentareay = 0
	for i in range (10):
		momentareaz +=  L[i]/(i+1) * (pow(tf,(i+1)) - pow(ti,(i+1)))
		momentareay += Ly[i]/(i+1) * (pow(tf,(i+1)) - pow(ti,(i+1)))

	Piy = Spline(A1,B1,C1,D1,ti) 
	Piz = Spline(A3,B3,C3,D3,ti)
	Pfy = Spline(A1,B1,C1,D1,tf)
	Pfz = Spline(A3,B3,C3,D3,tf)

	momentareay += Piy*Piy*Piz/2 - Pfy*Pfy*Pfz/2	# extra rectangular regions

	momentareaz = momentareaz		# I think the negatives should be here...
	momentareay = -momentareay			# But not here?

	return (momentareay, momentareaz)


def SolveCubic(a, b, c, d): #TO DO: done
	""" Only works if there's 1 real root: 
	
	q = (9*a*b*c - 27*a*a*d - 2*b*b*b)/(54*a*a*a)
	r = sqrt(pow(((3*a*c - b*b)/(9*a*a)),3) + q*q)
	s = pow((q+r),1./3.)
	if (q>r):
		t = pow((q-r),1./3.)
	else:
		t = -pow((r-q),1./3.)

	return (s + t - b/(3.*a))
	"""
	ap = b/a
	bp = c/a
	cp = d/a

	p = bp - ap*ap/3.
	q = cp + (2*ap*ap*ap - 9*ap*bp)/27.

	if (p == 0 and q == 0):
		return -ap/3.
	if (p == 0):
		magg = q
		argg = 0
		
	else:
		f = q*q/4. + p*p*p/27.
		if (f < 0) :
			rootfr = 0
			rootfi = sqrt(-f)
	
		else :
			rootfr = sqrt(f)
			rootfi = 0

		gr = rootfr + q/2.
		gi = rootfi

		magg = sqrt(gr*gr+gi*gi)
		argg = atan2(gi,gr) + (2*math.pi if rootfi < 0 else 0)
	

	magu1 = pow(magg,1./3.)
	argu1 = argg/3.

	u1r = magu1*cos(argu1)
	u1i = magu1*sin(argu1)

	u2r = -0.5*u1r - sqrt(3)/2.*u1i
	u2i = -0.5*u1i + sqrt(3)/2.*u1r

	u3r = -0.5*u1r + sqrt(3)/2.*u1i
	u3i = -0.5*u1i - sqrt(3)/2.*u1r

	x1r = p/3.*u1r/(u1r*u1r+u1i*u1i) - u1r - ap/3.
	x1i = -p/3.*u1i/(u1r*u1r+u1i*u1i) - u1i

	x2r = p/3.*u2r/(u2r*u2r+u2i*u2i) - u2r - ap/3.
	x2i = -p/3.*u2i/(u2r*u2r+u2i*u2i) - u2i

	x3r = p/3.*u3r/(u3r*u3r+u3i*u3i) - u3r - ap/3.
	x3i = -p/3.*u3i/(u3r*u3r+u3i*u3i) - u3i

	
	#	print(x1r, " + ", x1i, "i\n")
	#	print(x2r, " + ", x2i, "i\n")
	#	print(x3r, " + ", x3i, "i\n")
	
	
	# determine which real root falls between 0 and 1 - if none do, return -1
	
	threshold = 1e-10		# zero threshold
	ROUNDOFF = threshold*1000	# this one is killing me
	if (x1r <= 1+ROUNDOFF and x1r >= 0-ROUNDOFF and fabs(x1i)<threshold):
		return x1r
	elif (x2r <= 1+ROUNDOFF and x2r >= 0-ROUNDOFF and fabs(x2i)<threshold):
		return x2r
	elif (x3r <= 1+ROUNDOFF and x3r >= 0-ROUNDOFF and fabs(x3i)<threshold):
		return x3r

	return -1

def SectionResultant(P0y, P0z, P1y, P1z, P2y, P2z, P3y, P3z, d, theta):  # d is freeboard #TO DO: done, used tuples
	
	if (d > -P3z):
		magnitude = 0
		yloc = 0
		zloc = 0
		return (magnitude, yloc, zloc, 0)

	flag = 0

	A1 = P0y - 3*P1y + 3*P2y - P3y
	B1 = 3*P1y - 6*P2y + 3*P3y
	C1 = 3*P2y - 3*P3y
	D1 = P3y

	A3 = P0z - 3*P1z + 3*P2z - P3z
	B3 = 3*P1z - 6*P2z + 3*P3z
	C3 = 3*P2z - 3*P3z
	D3 = P3z

	ti = tfr = tfl= 0
	if (d == 0 and theta == 0):
		ti = 0
		tfr = 1
		tfl = 1
	else:
		ti = 0						# lower bound
		tfr = SolveCubic(A3-tan(theta)*A1, B3-tan(theta)*B1, C3-tan(theta)*C1, D3-tan(theta)*D1+d-P0z)	# upper bound right side
		if (tfr < 0):
			flag = -1          			# leak point
		tfl = SolveCubic(A3+tan(theta)*A1, B3+tan(theta)*B1, C3+tan(theta)*C1, D3+tan(theta)*D1+d-P0z)	# upper bound left side
		if (tfl < 0):
			flag = -1

	arearspline = SplineArea(P0y,P0z,P1y,P1z,P2y,P2z,P3y,P3z,ti,tfr)
	arealspline = SplineArea(P0y,P0z,P1y,P1z,P2y,P2z,P3y,P3z,ti,tfl)

	momyrspline, momzrspline = SplineMoment(P0y,P0z,P1y,P1z,P2y,P2z,P3y,P3z,ti,tfr) #returns momy first, then momz
	momylspline, momzlspline = SplineMoment(P0y,P0z,P1y,P1z,P2y,P2z,P3y,P3z,ti,tfl)

	momylspline = - momylspline

	ytfr = Spline(A1,B1,C1,D1,tfr);
	ztfr = Spline(A3,B3,C3,D3,tfr);
	ytfl = -Spline(A1,B1,C1,D1,tfl);
	ztfl = Spline(A3,B3,C3,D3,tfl);

	areartri = -ytfr*((P0z-d)-ztfr)/2;
	arealtri = -ytfl*((P0z-d)-ztfl)/2;

	arear = arearspline - areartri;		# subtract triangular portion - there's an explicit negative here
	areal = arealspline + arealtri;		# add triangular portion

	magnitude = arear + areal

	ymoment = momyrspline - 1./3.*ytfr*areartri + momylspline + 1./3.*ytfl*arealtri
	zmoment = momzrspline - (ztfr - 1./3.*(ztfr - (P0z-d)))*areartri+ momzlspline + (P0z-d - 2./3.*((P0z-d) - ztfl))*arealtri

	yloc = 0 if magnitude == 0 else ymoment/magnitude	
	zloc = 0 if magnitude == 0 else zmoment/magnitude

	if (flag == -1):
		return (magnitude, yloc, zloc, -1)
	return (magnitude, yloc, zloc, 0)



def Canoe_ControlPoints(station, P0y, P0z, P1y, P1z, P2y, P2z, P3y, P3z): #TO DO: arrays...syntax --> Python
	#x = station*increment                                            
	w = wlvalues[station]                                              
	d = klvalues[station]
	flareangle = 0
	fmax = asin(w/(2.*d))
	if (fmax < flare):                                                  #TO DO: wlvalues, klvalues, d, flare not defined???
		flareangle = fmax
	else:
		flareangle = flare

	u1max = d/cos(flareangle)
	u2max = w - d*tan(flareangle)
	u1 = u1max*(0.25+0.75*shapeparam)                                   #TO DO: shapearam not defined???
	u2 = u2max*(0.5+0.5*shapeparam)

	P0y = w
	P0z = 0
	P1y = P0y - u1*sin(flareangle)
	P1z = P0z - u1*cos(flareangle)
	P3y = 0
	P3z = -d
	P2y = P3y + u2
	P2z = P3z

	return (P0y, P0z, P1y, P1z, P2y, P2z, P3y, P3z)

# OVERLOADED FUNCTION
# This version works for any x value, not just at predefined stations. This one is more computationally intensive,
# so it should be used only sparingly
def Canoe_ControlPoints(x, P0y, P0z, P1y, P1z, P2y, P2z, P3y, P3z): #TO DO: why is this "expensive" function named the same as the one above...
	w = Waterline(x)                                            #TO DO: some functions might need to be moved above this one...undefined variables again
	d = Keelline(x)
	flareangle = 0
	fmax = acos(w/(2.*d))
	if (fmax < flare):
		flareangle = fmax
	else:
		flareangle = flare

	u1max = d/cos(flareangle)
	u2max = w - d*tan(flareangle)
	u1 = u1max*(0.25+0.75*shapeparam)
	u2 = u2max*(0.5+0.5*shapeparam)

	P0y = w
	P0z = 0
	P1y = P0y - u1*sin(flareangle)
	P1z = P0z - u1*cos(flareangle)
	P3y = 0
	P3z = -d
	P2y = P3y + u2
	P2z = P3z

	return (P0y, P0z, P1y, P1z, P2y, P2z, P3y, P3z)


def Ramp(x, t): #TO DO: done
	return (t/2. * log(0.5*exp(x/t) + 0.5*exp(-x/t)) + x/2. + t*log(2.)/2.)

def Canoe_Waterline(x): # TO DO: some variables are undefined like width lfirst, length, lpaddler....stuff like that but otherwise DONE
	m1 = width/(2.*lfirst)
	m2 = width/(2.*(length - lpaddler - lfirst))
	rbow = width/2. - m1*Ramp(width/(2.*m1),smooth1) - m2*Ramp(-length + width/(2.*m2),smooth2)
	rstern = width/2. - m1*Ramp(width/(2.*m1)-length,smooth1) - m2*Ramp(width/(2.*m2),smooth2)

	return (width/2. - m1*Ramp(width/(2.*m1)-x,smooth1) - m2*Ramp(x-length+width/(2.*m2),smooth2) - (rstern-rbow)*x/length - rbow)

def Canoe_Keelline(x): #TO DO: more undefined variables 
	if (x < ldeepest):
		return (depth - brocker*fabs(pow((1.-x/ldeepest),bowpower)))
	else:
		return (depth - srocker*fabs((pow((ldeepest-x)/(length-ldeepest),sternpower))))
	return 0

def Canoe_InitializeCanoe(L, Lp, Ld, Lf, W, t1, t2, d, h, b, s, f, n, density): # TO DO: k not sure what this was suppose to return also array stuff below 
	numstations = 101
	increment = L/(numstations-1)

	# Default tolerance values
	wltol = 1e-4
	tippingtol = 1e-2
	bowtol = 1e-4
	wlmaxit = 50
	tippingmaxit = 50

	wlvalues = [0]
	wlvalues = wlvalues*numstations
	klvalues = [0]
	klvalues = klvalues*numstations

	length = L
	lpaddler = Lp
	ldeepest = Ld
	lfirst = Lf
	width = W
	smooth1 = t1
	smooth2 = t2
	depth = d
	hflange = h
	brocker = b
	srocker = s
	flare = f
	shapeparam = n
	this->density = density #TO DO: dereferencing pointer I think...take this out and replace with something?

	for i in range (numstations):
		wlvalues[i] = Waterline(i*increment)
		klvalues[i] = Keelline(i*increment)
	
	# 3 men load case - turning is important       #TO DO: POINTERS (take away, just pass in array)...also these are structs...???
	loadcase[0].numpaddlers = 3
	loadcase[0].paddlerweight = 85		# about 185 pounds
	loadcase[0].paddlercm = .40		# 40 cm (rough estimate)

	# 2+2 load case - speed is important
	loadcase[1].numpaddlers = 4
	loadcase[1].paddlerweight = 72.5			# about 135 pounds for women
	loadcase[1].paddlercm = .35			# about 30 cm for women

	return 0


def Canoe_Destruct():   #TO DO: python is dynamic, might not need this function actually anymore
	delete[] wlvalues
	delete[] klvalues
	return 0

def Canoe_Canoe(L, Lp, Ld, Lf, W, t1, t2, d, h, b, s, f, n): #TO DO: were these values suppose to be returned by reference? there was no double&...
	length = L
	lpaddler = Lp
	ldeepest = Ld
	lfirst = Lf
	width = W
	smooth1 = t1
	smooth2 = t2
	depth = d
	hflange = h
	brocker = b
	srocker = s
	flare = f
	shapeparam = n
	
	#TO DO: no return????

def Canoe_Canoe(): #TO DO: what was this suppose to be?

def Canoe_Analyze(d, theta, trim, volume, cenx, ceny, cenz): # TO DO: POINTERS (take away, just pass in array)
	double* areas = new double[numstations];   
	double* centroidx = new double[numstations];
	double* centroidy = new double[numstations];
	double* centroidz = new double[numstations];

	x = 0

	P0y = 0
	P0z = 0
	P1y = 0
	P1z = 0
	P2y = 0
	P2z = 0
	P3y = 0
	P3z = 0
	flag = 0
	
	for i in range (numstations): #TO DO: numstations not defined
		
		P0y, P0z, P1y, P1z, P2y, P2z, P3y, P3z = Canoe_ControlPoints(i, P0y,P0z,P1y,P1z,P2y,P2z,P3y,P3z) 
		# freeboard at distance x is given by d - trim/length*x
		if (SectionResultant(P0y,P0z,P1y,P1z,P2y,P2z,P3y,P3z, (d if trim == 0 else d-trim/length*x),theta, areas[i],centroidy[i],centroidz[i])[3] == -1):
			#print ("LEAK\n")			# let it pass for now 
			flag = 1
			#return 1
			
		centroidx[i] = x

		x += increment               #TO DO: increment not defined

	volume = 0
	for i in range (numstations):                      #numerically integrate using parallelograms
		volume += (areas[i]+areas[i-1])/2
	
	volume *= increment

	cenx = 0
	ceny = 0
	cenz = 0

	for in range (numstations):
		cenx += centroidx[i]*areas[i]*increment			# this is not 100% correct... but it's not too far wrong
		ceny += centroidy[i]*areas[i]*increment
		cenz += centroidz[i]*areas[i]*increment

	cenx /= volume
	ceny /= volume
	cenz /= volume

	delete[] areas                 #TO DO: recall: python --> dynamic, don't need these anymore I think?
	delete[] centroidx
	delete[] centroidy
	delete[] centroidz

	if (flag == 1):
		return (volume, cenx, ceny, cenz, 1)

	return (volume, cenx, ceny, cenz, 0)

def Canoe_FindWLine(vtarget, theta, trim):
	v= 0
	cx = 0
	cy = 0
	cz = 0
	
	MAXITERATIONS = wlmaxit   #TO DO: wlmaxit, wltotal not defined
	tolerance = wltol

	testd = depth/2
	width = depth/2

	flag = 0
	# Check if it's possible to attain equilibrium (only valid for the zero trim case)
	if (trim == 0):
		width = GetBWL(0)/2.  #TO DO: undefined function...may have to reorganize order functions are defined
		dcheck = width*tan(theta)
		v, cx, cy, cz , flag= Canoe_Analyze(dcheck,theta,trim,v,cx,cy,cz) 
		if (v < vtarget):
			return -1
	

	# binary-type algorithm
	for i in range (MAXITERATIONS):
		flag = Canoe_Analyze(testd,theta,trim,v,cx,cy,cz)[4]
		if (fabs(v-vtarget) < tolerance):
			return testd
#		if (v < vtarget and flag == 1):
#			return -1						# equilibrium not possible -- leak
		if (v > vtarget or flag == 1):
			testd += width/2.
		else:
			testd -= width/2.
		width /= 2.
	
#	print ("WARNING: LOSS OF PRECISION\n")                           
	return -1


def Moment( y1, z1, y2, z2, weight, theta): # TO DO: done
	return (-(cos(theta)*(y2-y1) - sin(theta)*(z2-z1))*weight)


# NOTE - cmheight is the centre of mass height taken ABOVE THE TOP OF THE GUNWALE
# The tipping angle calculation assumes the leak happens at the widest point (which is always true for zero trim)
def Canoe_CrossCurve(weight, cmheight, tangent, tipping): # TO DO: return tangent and tipping by REF????
	theta = 0
	tinc = 0.01
	v = 0
	cx = 0
	cy = 0
	cz = 0
	
	# roughly get tangent	
	theta1 = theta
	theta2 = theta+tinc
	d1 = FindWLine(weight, theta1, 0)
	v, cx, cy, cz , flag = Canoe_Analyze(d1,theta1,0,v,cx,cy,cz)                            
	moment1 = Moment(cy,cz,0,cmheight,weight,theta)
	d2 = FindWLine(weight, theta2, 0)
	v, cx, cy, cz , flag = Canoe_Analyze(d2,theta2,0,v,cx,cy,cz)
	moment2 = Moment(cy,cz,0,cmheight,weight,theta)
	tangent = (moment2-moment1)/(theta2-theta1)

	# find critical tipping angle
	tolerance = tippingtol
	width = GetBWL(0)/2.
	dtest = d1					# start off with d for 0 heel
	dmax = depth					# initial maximum value of d is depth of canoe
	dmin = 0.						# initial minimum value of d is 0


	theta = atan2(dtest,width)    
	v, cx, cy, cz , flag = Canoe_Analyze(dtest,theta,0,v,cx,cy,cz)
	MAXITERATIONS = tippingmaxit
	count = 0
	flag = 0
	while (fabs(weight-v) > tolerance):
		if (count >= MAXITERATIONS):
			flag = 1
			break
		
		count += 1                   
		if (weight < v):
			dmin = dtest
			dtest = (dtest+dmax)/2.
		else:
			dmax = dtest
			dtest = (dtest+dmin)/2.
		
		theta = atan2(dtest,width)
		v, cx, cy, cz , flag = Canoe_Analyze(dtest,theta,0,v,cx,cy,cz)

	tipping = theta

	return flag


""" OLD algorithm - it works too but it is very slow
# NOTE - cmheight is the centre of mass height taken ABOVE THE TOP OF THE GUNWALE
def Canoe_CrossCurve(weight, cmheight, tangent, tipping): #return tangent and tipping REF
	theta = 0
	tinc = 0.01
	v = 0
        cx = 0
        cy = 0
        cz = 0
	
	# roughly get tangent	
	d = FindWLine(weight, theta, 0)
	theta1 = theta
	theta2 = theta+tinc
	v, cx, cy, cz , flag = Canoe_Analyze(d,theta1,0,v,cx,cy,cz)
	moment1 = Moment(cy,cz,0,cmheight,weight,theta)
	v, cx, cy, cz , flag = Canoe_Analyze(d,theta2,0,v,cx,cy,cz)
	moment2 = Moment(cy,cz,0,cmheight,weight,theta)
	tangent = (moment2-moment1)/(theta2-theta1)

	while (d >= 0):
		v, cx, cy, cz , flag = Canoe_Analyze(d,theta,0,v,cx,cy,cz)
		
		output << theta << '\t' << Moment(cy,cz,0,cmheight,weight,theta) << endl #TO DO: output is not the same as cout??

		theta += tinc

		d = FindWLine(weight, theta, 0)

	tipping = theta - tinc

	return 0

""" 

def Canoe_DisplacementCurve(): 
	d = depth
	dinc = 0.01
	v = 0
	cx = 0
	cy = 0
	cz = 0
	
	while (d >= 0):
		v, cx, cy, cz , flag = Canoe_Analyze(d,0,0,v,cx,cy,cz)

		print( d , '\t' , v )

		d -= dinc
	
	return 0


def SplineLength(P0y, P0z, P1y, P1z, P2y, P2z, P3y, P3z, d): #TO DO: done
	if (d > -P3z):
		return (0, 0)
	

	A1 = P0y - 3*P1y + 3*P2y - P3y
	B1 = 3*P1y - 6*P2y + 3*P3y
	C1 = 3*P2y - 3*P3y
	D1 = P3y

	A3 = P0z - 3*P1z + 3*P2z - P3z
	B3 = 3*P1z - 6*P2z + 3*P3z
	C3 = 3*P2z - 3*P3z
	D3 = P3z

	ti = 0							# lower bound
	tf = SolveCubic(A3, B3, C3, D3+d-P0z)	# upper bound
	if (tf < 0):
		tf = 0    	# if there's a leak we ignore the section - presumably there should not be any leaks in this routine
	if (d == 0):
		tf = 1

	t = 0
	tinc = 0.01	# UNDEFINED CONSTANT
	y1 = 0
	z1 = 0
	y2 = 0
	z2 = 0
	l = 0

	y1 = Spline(A1,B1,C1,D1,ti)
	z1 = Spline(A3,B3,C3,D3,ti)

	zmoment = 0

	for t in range (ti+tinc, tf + 1, tinc):
		y2 = Spline(A1,B1,C1,D1,t)
		z2 = Spline(A3,B3,C3,D3,t)
		segment = sqrt((y2-y1)*(y2-y1) + (z2-z1)*(z2-z1))
		l += segment
		zmoment += segment*(z1+z2)/2.		# moment = length * average height

		y1 = y2
		z1 = z2
	

	cmz = zmoment/l 	#z centre of mass

	return (cmz, l)

def Canoe_SurfaceArea(d, cx, cz): #TO DO: POINTERS (take away, just pass in array)
	area = 0
	areamoment = 0
	areazmoment = 0
	
	double* lengths = new double[numstations] #TO DO: POINTERS (take away, just pass in array)
	double* cenzs = new double[numstations]
	P0y = 0
	P0z = 0
	P1y = 0
	P1z = 0
	P2y = 0
	P2z = 0
	P3y = 0
	P3z = 0
#	x = 0
	cmz = 0
	
	for i in range (numstations):
		P0y, P0z, P1y, P1z, P2y, P2z, P3y, P3z = Canoe_ControlPoints(i, P0y,P0z,P1y,P1z,P2y,P2z,P3y,P3z)
		lengths[i] = SplineLength(P0y,P0z,P1y,P1z,P2y,P2z,P3y,P3z,d)[1] 
		#x += increment;
		cmz = SplineLength(P0y,P0z,P1y,P1z,P2y,P2z,P3y,P3z,d)[0]
		cenzs[i] = cmz 
	
	for i in range (numstations): 							# numerically integrate using parallelograms
		area += (lengths[i]+lengths[i-1])/2
		areamoment += (lengths[i]*i*increment)
		areazmoment += cenzs[i]*lengths[i]

	cx = areamoment/area
	cz = areazmoment/area

	area *= increment #TO DO: change the pointer stuff again
	area *= 2.

	delete[] lengths;
	delete[] cenzs;

	return ( cx, cz, area) 
    
# John Winters' KAPER formula - disp is in long tonnes, lengths in feet, speed in knots, LCB is a fraction of length
def Canoe_Kaper(BWL, EWL, WS, Cp, Cv, LCB, disp, le, v): #TO DO: done
	vtol = v/sqrt(EWL)

	c1 = 0.002*sqrt(BWL/EWL)*pow((4*vtol),4)
	c2 = 0.005*(sin(le*360/2./PI))*(4*vtol)*(4*vtol)
	c3 = (1 if vtol < 1.5 else 0.7)*0.8*cos(3.65*vtol+0.07) + (1 if vtol < 1.5 else 0.96)*1.2
	c4 = GetC4(Cv,Cp,vtol)
	c5 = pow((0.5/LCB),0.35)
	c6 = 1.0

	Rr = (4*c5*disp*pow(vtol,4) + c1 + c2)*(c3*c4*c5)
	if (v >= 1.6 and v < 3):
		Rr *= 0.75
	if (v >= 1.4 and v < 1.6):
		Rr *= 0.85

	reynold = v*1.6889*EWL/1.2791*100000
	Cf = 0.075/pow((log10(reynold)-2),2)
	Rf = 0.99525*Cf*WS*(v*1.6889)*(v*1.6889)		# speed in feet/sec

	return (Rf+Rr)


def ReadC4Table(ifstream& c4): #TO DO: what is IFSTREAM&?
	
	for i in range (3):
		for j in range (17):
			for k in range (14):
				c4 >> C4TABLE[i][j][k] #TO DO: what is c4? is that just reading input again with the word c4?
		
	return 0


# Look up c4 in the table
def GetC4(Cv, Cp, vtol): #TO DO: done
	Cvind1 = int(((Cv - .001)/.0005))
	Cvind2 = int(((Cv - .001)/.0005) + 1)
	if (Cvind1 < 0):
		Cvind1 = 0
	if (Cvind1 > 2):
		Cvind1 = 2
	if (Cvind2 < 0):
		Cvind2 = 0
	if (Cvind2 > 2):
		Cvind2 = 2

	Cpind1 = int(((Cp - 0.48)/0.01))
	Cpind2 = int(((Cp - 0.48)/0.01) + 1)
	if (Cpind1 < 0):
		Cpind1 = 0
	if (Cpind1 > 16) :
		Cpind1 = 16
	if (Cpind2 < 0):
		Cpind2 = 0
	if (Cpind2 > 16):
		Cpind2 = 16

	vtolind1 = int(((vtol - 0.4)/0.1))
	vtolind2 = int(((vtol - 0.4)/0.1) + 1)
	if (vtolind1 < 0):
		vtolind1 = 0
	if (vtolind1 > 13):
		vtolind1 = 13
	if (vtolind2 < 0):
		vtolind2 = 0
	if (vtolind2 > 13):
		vtolind2 = 13

	Cvratio = 0
	Cpratio = 0
	vtolratio = 0
	if (Cvind1 == Cvind2):
		Cvratio = 0
	else:
		Cvratio = (Cv - (Cvind1*.0005 + .001))/(Cvind2*.0005 - Cvind1*.0005)
	if (Cpind1 == Cpind2):
		Cpratio = 0
	else:
		Cpratio = (Cp - (Cpind1*.01 + .48))/(Cpind2*.01 - Cpind1*.01)
	if (vtolind1 == vtolind2):
		vtolratio = 0
	else:
		vtolratio = (vtol - (vtolind1*.1 + .4))/(vtolind2*.1 - vtolind1*.1)

	# 8 corners of the prism to be interpolated
	corners = [[[0 for x in range(2)] for y in range(2)] for z in range(2)] 	

	for i in range (Cvind1, Cvind2+1, 1):
		for j in range (Cpind1, Cpind2, 1):
			for k in range (vtolind1, vtolind2+1, 1):
				corners[i-Cvind1][j-Cpind1][k-vtolind1] = C4TABLE[i][j][k]

	actual = [[0 for x in range(2)] for y in range(2)]
	for i in range (2): 
		for j in range (2):
			actualCv[i][j] = Cvratio*corners[1][i][j] + (1-Cvratio)*corners[0][i][j]

	actualCp = [0 for x in range 2]
	for i in range (2):
		actualCp[i] = Cpratio*actualCv[1][i] + (1-Cpratio)*actualCv[0][i]

	c4 = vtolratio*actualCp[1] + (1-vtolratio)*actualCp[0]
	
	print ("Cv from " , Cvind1 , " to " , Cvind2 )                             
	print ("Cp from " , Cpind1 , " to " , Cpind2 )
	print ("VtoL from ", vtolind1 , " to " , vtolind2 )

	print ("Cv ratio: " , Cvratio )
	print ("Cp ratio: " , Cpratio )
	print ("VtoL ratio: " , vtolratio)

	print ("C4 = " , c4 )
	

	return c4


# Get "Effective waterline length" (metres)
def Canoe_GetEWL(d, theta, trim): #TO DO: pointers
	double* areas = new double[numstations] #TO DO: POINTERS (take away, just pass in array)

	x = 0
	dummy1 = dummy2 = 0
	areamax = 0
	xmax = 0

	P0y = P0z = P1y = P1z = P2y = P2z = P3y = P3z = 0
	
	for i in range (numstations):
		P0y, P0z, P1y, P1z, P2y, P2z, P3y, P3z = Canoe_ControlPoints(i, P0y,P0z,P1y,P1z,P2y,P2z,P3y,P3z)
		# freeboard at distance x is given by d - trim/length*x
		if (SectionResultant(P0y,P0z,P1y,P1z,P2y,P2z,P3y,P3z, d-trim/length*x,theta, areas[i],dummy1,dummy2)[3] == -1):
			print("LEAK\n") 
			delete[] areas #TO DO: recall python -->dynamic memory, don't need this anymore
			return 1

		if (areas[i] > areamax): 
			areamax = areas[i]
			xmax = i

		x += increment

	"""/*
	output << "\n\nCurve of Areas\n";
	for i in range (numstations): 
		output << i*increment << '\t' << areas[i] << endl;
	
	*/""" 

	for xstart in range (xmax):
		flag = false
		for i in range (xmax):
			if ((i-xstart)/(xmax-xstart)*areamax > areas[i]):			# Find (approx) tangent to the curve
				flag = true
				break
			
		if (not flag):
			break	# We now have xstart for EWL

	for  xend in range(numstations-1, xmax, -1):
		flag = false
		for i in range(xend, xmax, -1):
			if ((xend-i)/(xend-xmax)*areamax > areas[i]):		# Find (approx) tangent to the curve
				flag = true
				break
			
		if (not flag):
			break	# We now have xend for EWL

	EWL = (xend - xstart)*increment
#	print("EWL = " , EWL)
	delete[] areas #TO DO:  python is dynamic 

	return EWL

# Get waterline length (metres)
def Canoe_GetLWL(d): #TO DO: done

	flag = false
	starti = endi = 0

	for i in range (numstations):
		if (klvalues[i] > d):
			flag = true
			starti = i
			break

	if (not flag):
		return length

	flag = false
	for i in range(numstations-1, -1, -1):
		if (klvalues[i] > d):
			flag = true
			endi = i
			break

	return ((endi-starti)*increment)

# Get waterline beam (metres)
def Canoe_GetBWL(d): #TO DO: done
	BWL = 0
	P0y = P0z = P1y = P1z = P2y = P2z = P3y = P3z = 0

	width = 0

	for i in range (numstations):
		P0y, P0z, P1y, P1z, P2y, P2z, P3y, P3z = Canoe_ControlPoints(i, P0y,P0z,P1y,P1z,P2y,P2z,P3y,P3z)
		A3 = P0z - 3*P1z + 3*P2z - P3z
		B3 = 3*P1z - 6*P2z + 3*P3z
		C3 = 3*P2z - 3*P3z
		D3 = P3z		
		tf = SolveCubic(A3, B3, C3, D3+d-P0z)

		A1 = P0y - 3*P1y + 3*P2y - P3y
		B1 = 3*P1y - 6*P2y + 3*P3y
		C1 = 3*P2y - 3*P3y
		D1 = P3y
		
		width = Spline(A1,B1,C1,D1,tf)
		if (width > BWL):
			BWL = width
		
	return BWL*2
