import numpy as np
import math

def cross(a, b):
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]

def raytraceSphere():
    sphereList = []
    rayList = []
    rayDirection = []
    for i in range(0,3):
        element = float(input("Enter sphere center vector(3 numbers): "))
        sphereList.append(element)
    for i in range(0,3):
        element = float(input("Enter ray vector(3 numbers): "))
        rayList.append(element)
    for i in range(0,3):
        element = float(input("Enter ray direction vector(3 numbers): "))
        rayDirection.append(element)
        
    #Initialize vectors     
    sphereCenterVector = np.array(sphereList)
    rayVector = np.array(rayList)
    rayDirectionVector = np.array(rayDirection)
    #print("Sphere Center Vector: ")
    #print(sphereCenterVector)
    #print("Ray Vector: ")
    #print(rayVector)
    #print("Ray Direction Vector: ")
    #print(rayDirectionVector)
    radius = int(input("Enter sphere radius: "))

    #Angle between ray and vector v1 (from sphere center to camera)
    v1Vector = sphereCenterVector - rayVector
    #print("V1: ")
    #print(v1Vector)
    norm = np.linalg.norm(v1Vector)
    v1Normalized = v1Vector/norm
    #print(v1Normalized)

    dotProduct = v1Normalized.dot(rayDirectionVector)
    theta = math.acos(dotProduct)
    thetaDegrees = math.degrees(theta)
    print("The angle between ray and vector v1 is :")
    print(thetaDegrees)
    print("\n")

    #Distance d - from center of sphere to interpolated ray
    D1 = v1Vector.dot(rayDirectionVector)
    #print(D1)
    d = D1 * math.tan(theta)
    print("Distance d is :")
    print(d)
    print("\n")
    
    #Distance a
    a = math.sqrt(radius**2 - d**2)
    print("Distance a is :")
    print(a)
    print("\n")
    
    #Distance b
    b = D1 - a
    print("Distance b is :")
    print(b)
    print("\n")
    
    #Final ray to intersection point
    I = (rayDirectionVector * b) + rayVector
    print("Intersection is: ")
    print(I)
    print("\n")
    
    #Normal vector at I
    N = I - sphereCenterVector
    norm2 = np.linalg.norm(N)
    N = N / norm2
    print("Normal vector at I is: ")
    print(N)
    print("\n")
    
def rayPlane():
    pointList = []
    normalList= []
    rayList = []
    rayDirection = []
    for i in range(0,3):
        element = float(input("Enter point vector(3 numbers): "))
        pointList.append(element)
    for i in range(0,3):
        element = float(input("Enter normal vector(3 numbers): "))
        normalList.append(element)
    for i in range(0,3):
        element = float(input("Enter ray vector(3 numbers): "))
        rayList.append(element)
    for i in range(0,3):
        element = float(input("Enter ray direction vector(3 numbers): "))
        rayDirection.append(element)
        
    #Initialize vectors 
    pointVector = np.array(pointList)
    normalVector = np.array(normalList)
    rayVector = np.array(rayList)
    rayDirectionVector = np.array(rayDirection)
    
    #Distance from ray origin to plane (a)
    w = pointVector - rayVector
    a = w.dot(normalVector)
    print("Distance a is: ")
    print(a)
    print("\n")
    
    #Distance b
    b = rayDirectionVector.dot(normalVector)
    print("Distance b is: ")
    print(b)
    print("\n")
    
    #Ratio k
    k = a/b
    print("Ratio k is: ")
    print(k)
    print("\n")
    
    #Intersection 
    I = (rayDirectionVector * k) + rayVector
    print("I is: ")
    print(I)
    print("\n")

def bezierSpline():
    p0List= []
    p1List = []
    p2List = []
    for i in range(0,3):
        element = float(input("Enter p0 vector(3 numbers): "))
        p0List.append(element)
    for i in range(0,3):
        element = float(input("Enter p1 vector(3 numbers): "))
        p1List.append(element)
    for i in range(0,3):
        element = float(input("Enter p2 vector(3 numbers): "))
        p2List.append(element)

    #Initialize vectors 
    p0Vector = np.array(p0List)
    p1Vector = np.array(p1List)
    p2Vector = np.array(p2List)

    R1 = ((p2Vector - p0Vector) * (1.0/6.0)) + p1Vector
    L2 = (p1Vector - R1) + p1Vector
    L1 = ((p0Vector - p1Vector) * (1.0/3.0)) + L2
    R2 = ((p2Vector - p1Vector) * (1.0/3.0)) + R1

    print("R1: ")
    print(R1)
    print("\n")

    print("L2: ")
    print(L2)
    print("\n")

    print("L1: ")
    print(L1)
    print("\n")

    print("R2: ")
    print(R2)
    print("\n")

def barycentricCrossProduct():
    p0List = []
    p1List = []
    for i in range(0,3):
        element = float(input("Enter p0 vector(3 numbers): "))
        p0List.append(element)
    for i in range(0,3):
        element = float(input("Enter p1 vector(3 numbers): "))
        p1List.append(element)
        
    p0Vector = np.array(p0List)
    p1Vector = np.array(p1List)
    cp = np.cross(p0Vector,p1Vector)
    print(cp)

    magnitude = np.linalg.norm(cp)
    print("Magnitude of cross product: ")
    print(magnitude)
    print("\n")

def bilinearInterpolation():
    pLeftList = []
    pRightList = []
    pList = []
    pLeftInfoList = []
    pRightInfoList = []
    for i in range(0,3):
        element = float(input("Enter pLeft vector(3 numbers): "))
        pLeftList.append(element)
    for i in range(0,3):
        element = float(input("Enter pRight vector(3 numbers): "))
        pRightList.append(element)
    for i in range(0,3):
        element = float(input("Enter p vector(3 numbers): "))
        pList.append(element)
    for i in range(0,3):
        element = float(input("Enter pLeft info vector(normal/color): "))
        pLeftInfoList.append(element)
    for i in range(0,3):
        element = float(input("Enter pRight vector(normal/color): "))
        pRightInfoList.append(element)
    
    pLeftVector = np.array(pLeftList)
    pRightVector = np.array(pRightList)
    pVector = np.array(pList)
    pLeftInfoVector = np.array(pLeftInfoList)
    pRightInfoVector = np.array(pRightInfoList)
    
    d1 = pVector - pLeftVector
    d1Mag = np.linalg.norm(d1)
    d2 = pVector - pRightVector
    d2Mag = np.linalg.norm(d2)
    d = d1Mag + d2Mag
    
    pInfoVector = ( (d2Mag/d) * pLeftInfoVector) + ( (d1Mag/d) * pRightInfoVector)
    print("P info vector is: ")
    print(pInfoVector)
    print("\n")

def clipping():
    aList = []
    bList = []
    normList = []
    pointList = []
    
    for i in range(0,3):
        element = float(input("Enter a vector(3 numbers): "))
        aList.append(element)
    for i in range(0,3):
        element = float(input("Enter b vector(3 numbers): "))
        bList.append(element)
    for i in range(0,3):
        element = float(input("Enter normal vector(3 numbers): "))
        normList.append(element)
    for i in range(0,3):
        element = float(input("Enter point vector(3 numbers): "))
        pointList.append(element)
    aVector = np.array(aList)
    bVector = np.array(bList)
    normVector = np.array(normList)
    pointVector = np.array(pointList)

    d1 = normVector.dot(aVector - pointVector)
    d2 = normVector.dot(bVector - pointVector)
    t = d1/(d1-d2)
    I = aVector + t*(bVector - aVector)
    print("Intersection point is: ")
    print(I)
    print("\n")

def cameraTrans():
    # input: camPos, tarPos, upVector (I think upvector is always 0,0,1) true ask tommy he knows
    upVector = np.array([0,0,1])
    t = []
    c = []
    for i in range(0,3):
        element = float(input("Enter Camera Position vector(3 numbers): "))
        c.append(element)
    for i in range(0,3):
        element = float(input("Enter Target Position vector(3 numbers): "))
        t.append(element)
    # finding and normalizing w
    cl = np.array(t) - np.array(c)
    w = (np.array(c) - np.array(t)) / np.linalg.norm(np.array(c) - np.array(t))
    u = np.cross(upVector, w)
    v = np.cross(w, u)
    print("w: ", w, "\n")
    print("u: ", u, "\n")
    print("v: ", v, "\n")
    print("c: ", cl, "\n")
    arr = np.array(( (u[0], v[0], w[0], cl[0]),
                     (u[1], v[1], w[1], cl[1]),
                     (u[2], v[2], w[2], cl[2]),
                     (0, 0, 0, 1)))
    print(arr)
    #just print each line of the matrix one by one not each element so like u[0], v[0] w[0] changeList[0]

def quadratic(p1,p2,p3,t):
    p12 = np.array(p2) - np.array(p1)
    p23 = np.array(p3) - np.array(p2)
    c1 = (p12 * t) + np.array(p1)
    c2 = (p23 * t) + np.array(p2)
    c3 = c2 - c1
    return (c3 * t) + c1

def cubic(p1,p2,p3,p4,t):
    p12 = np.array(p2) - np.array(p1)
    p23 = np.array(p3) - np.array(p2)
    p34 = np.array(p4) - np.array(p3)
    c1 = (p12 * t) + np.array(p1)
    c2 = (p23 * t) + np.array(p2)
    c3 = (p34 * t) + np.array(p3)
    return quadratic(c1,c2,c3,t)

def quartic(p1,p2,p3,p4,p5,t):
    p12 = np.array(p2) - np.array(p1)
    p23 = np.array(p3) - np.array(p2)
    p34 = np.array(p4) - np.array(p3)
    p45 = np.array(p5) - np.array(p4)
    c1 = (p12 * t) + np.array(p1)
    c2 = (p23 * t) + np.array(p2)
    c3 = (p34 * t) + np.array(p3)
    c4 = (p45 * t) + np.array(p4)
    return cubic(c1,c2,c3,c4,t)
    
def bezierCurve():
    p1 = []
    p2 = []
    p3 = []
    for i in range(0,3):
        element = float(input("Enter p1 (3 numbers): "))
        p1.append(element)
    for i in range(0,3):
        element = float(input("Enter p2 (3 numbers): "))
        p2.append(element)
    for i in range(0,3):
        element = float(input("Enter p3 (3 numbers): "))
        p3.append(element)
    t = float(input("time t: "))
    numPoints = (input("Number of points: "))
    
    if numPoints == "3":
        res = quadratic(p1,p2,p3,t)
        print("result:", res)
    elif numPoints == "4":
        p4 = []
        for i in range(0,3):
            element = float(input("Enter p4 (3 numbers): "))
            p4.append(element)
        res = cubic(p1,p2,p3,p4,t)
        print("result:", res)
    elif numPoints == "5":
        p4 = []
        p5 = []
        for i in range(0,3):
            element = float(input("Enter p4 (3 numbers): "))
            p4.append(element)
        for i in range(0,3):
            element = float(input("Enter p5 (3 numbers): "))
            p5.append(element)
        res = quartic(p1,p2,p3,p4,p5,t)
        print("result: ", res)

def torusDistance():
    torusC=[]
    torusN=[]
    rayO = []
    rayD = []
     
     
    for i in range(0,3):
        element = float(input("Enter Torus Centre (3 numbers): "))
        torusC.append(element)
    for i in range(0,3):
        element = float(input("Enter Torus Normal (3 numbers): "))
        torusN.append(element)
    for i in range(0,3):
        element = float(input("Enter ray origin (3 numbers): "))
        rayO.append(element)
    for i in range(0,3):
        element = float(input("Enter ray direction (3 numbers): "))
        rayD.append(element)
    torusR1 = float(input("Enter torusR1 distance between centre of torus and centre of tube direction: "))
    torusR2 = float(input("Enter torusR2 tube radius: "))

    kPre = (np.array(rayO)-np.array(torusC))
    k    = kPre.dot(np.array(torusN)) 
    pPre = k*(np.array(torusN))
    p    = np.array(rayO) - pPre
    print ("p: ", p)
    bigM = (p-np.array(torusC))/np.linalg.norm(p-np.array(torusC))
    m = bigM*(torusR1)+np.array(torusC)
    print ("m: ", m)
    D = np.linalg.norm(m-np.array(rayO))-torusR2
    print("D: ", D)

def rayBox():
    rayOriginList = []
    rayDirectionList = []
    minList = []
    maxList = []

    for i in range(0,2):
        element = float(input("Enter ray origin (2 numbers): "))
        rayOriginList.append(element)
    for i in range(0,2):
        element = float(input("Enter ray direction (2 numbers): "))
        rayDirectionList.append(element)
    for i in range(0,2):
        element = float(input("Enter min (2 numbers): "))
        minList.append(element)
    for i in range(0,2):
        element = float(input("Enter max (2 numbers): "))
        maxList.append(element)

    #Initialize vectors
    rayDir = np.array(rayDirectionList)
    rayOri = np.array(rayOriginList)
    minimum = np.array(minList)
    maximum = np.array(maxList)

    #minY
    R = rayDir * 999
    a = abs(minimum[1] - rayOri[1])
    k = a / abs((R[1] - rayOri[1]))
    minY = ((R - rayOri) * k) + rayOri

    #maxY
    a = abs(maximum[1] - rayOri[1])
    k = a / abs((R[1] - rayOri[1]))
    maxY = ((R - rayOri) * k) + rayOri

    #minX
    a = abs(minimum[0] - rayOri[0])
    k = a / abs((R[0] - rayOri[0]))
    minX = ((R - rayOri) * k) + rayOri

    #maxX
    a = abs(maximum[0] - rayOri[0])
    k = a / abs((R[0] - rayOri[0]))
    maxX = ((R - rayOri) * k) + rayOri

    print("Min Y: ")
    print(minY)
    print("Max Y: ")
    print(maxY)
    print("Min X: ")
    print(minX)
    print("Max X: ")
    print(maxX)



#Main Menu
while True:
    print ("1. Raytracing a sphere")
    print ("2. Ray plane intersection")
    print ("3. Bezier Spline")
    print ("4. Barycentric Cross Product")
    print ("5. Bilinear Interpolation")
    print ("6. Clipping")
    print ("7. Camera Transformation")
    print ("8. Bezier Curve")
    print ("9. Distance to Torus")
    print ("10. ray box intersection")
    ans = input("Which computation would you like to do today?\n")
    if ans == "1":
        raytraceSphere()
    elif ans == "2":
        rayPlane()
    elif ans == "3":
        bezierSpline()
    elif ans == "4":
        barycentricCrossProduct()
    elif ans == "5":
        bilinearInterpolation()
    elif ans == "6":
        clipping()
    elif ans == "7":
        cameraTrans()
    elif ans == "8":
        bezierCurve()
    elif ans == "9":
        torusDistance()
    elif ans == "10":
        rayBox()
