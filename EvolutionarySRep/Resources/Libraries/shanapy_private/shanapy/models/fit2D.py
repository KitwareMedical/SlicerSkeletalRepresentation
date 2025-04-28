import vtk
import math
import numpy as np
import matplotlib
import scipy
from scipy.interpolate import splprep, splev
import matplotlib.pyplot as plt
from shapely.geometry import Polygon, Point

# fits a circle through 3 points
# finds slope of lines perpendicular to the 2 lines found previously
# draws perpendicular line from the bisection of the first two lines. See where those perpendicular lines intersect
# intersection point is the center of the circle being fit
# Radius is the distance between that center point and one of the original 3 points
# fits a circle through 3 points: Center (middle pt), pLeft (left pt), pRight (right point)
def get_circle(center, pLeft, pRight):
    diff1 = (pLeft[0]-center[0], pLeft[1]-center[1])
    diff2 = (pRight[0]-center[0], pRight[1]-center[1])
    # finds slope of the lines from the middle point to the left point (then get negative reciprocal) to get perpendicular line slope
    slope1 = -1*diff1[0] / diff1[1]
    # finds slope of the lines from the middle point to the right point (then get negative reciprocal) to get perpendicular line slope
    slope2 = -1*diff2[0] / diff2[1]
    # finds the middle of the line from center to pLeft
    mid1 = ((center[0]+pLeft[0])/2, (center[1]+pLeft[1])/2)
    # finds the middle of the line from center to pRight
    mid2 = ((center[0]+pRight[0])/2, (center[1]+pRight[1])/2)
    
    # need to find where the perpendicular bisectors of the two lines intersect
    # int1 is the y-intercept of the first perpendicular bisector
    int1 = mid1[1] - (slope1*mid1[0])
    # int2 is the y-intercept of the second perpendicular bisector
    int2 = mid2[1] - (slope2*mid2[0])
    
    # solving for x-val of the intersection point (this is the x-val of the center of the circle)
    finX = (int2-int1)/(slope1-slope2)
    # solving for y-val of the intersection point (this is the y-val of the center of the circle)
    finY = slope1*finX+int1
    mid = (finX, finY)
    # calculating circle radius
    dist = math.sqrt((mid[0]-center[0])**2 + (mid[1]-center[1])**2)
    return mid, dist


# creating 2x2 second moment matrix
# eigenvalues of the matrix are the lengths of the best fitting ellipse radii
#eigenvectors are the directions of the best fitting ellipse radii
def check_ellipse(xVals, yVals):
    N = len(xVals)
    xmean, ymean = xVals.mean(), yVals.mean()
    xVals = xVals - xmean
    yVals = yVals - ymean
    U, S, V = np.linalg.svd(np.stack((xVals, yVals)))

    transform = np.sqrt(2/N) * U.dot(np.diag(S))   # transformation matrix
    w = [math.hypot(transform[0,0], transform[1,0]), math.hypot(transform[0,1], transform[1,1])]
    v = [[transform[0,0]/w[0], transform[0,1]/w[1]], [transform[1,0]/w[0], transform[1,1]/w[1]]]
    return np.asarray(w), np.asarray(v)

# "solve" finds the nearest point on the ellipse given another point
# takes in the lengths of the two ellipse axes and a point
# Assumes that the ellipse is centered at the  origin and that the semi_major axis is along the x-axis
# ***so if you want to find the nearest point on an ellipse to another given point, note the method assumptions
# may need to translate and rotate point before using this method
def solve(semi_major, semi_minor, p):  
    px = abs(p[0])
    py = abs(p[1])

    tx = 0.707
    ty = 0.707

    a = semi_major
    b = semi_minor

    for x in range(0, 3):
        x = a * tx
        y = b * ty

        ex = (a*a - b*b) * tx**3 / a
        ey = (b*b - a*a) * ty**3 / b

        rx = x - ex
        ry = y - ey

        qx = px - ex
        qy = py - ey

        r = math.hypot(ry, rx)
        q = math.hypot(qy, qx)

        tx = min(1, max(0, (qx * r / q + ex) / a))
        ty = min(1, max(0, (qy * r / q + ey) / b))
        t = math.hypot(ty, tx)
        tx /= t 
        ty /= t 

    return [math.copysign(a * tx, p[0]), math.copysign(b * ty, p[1])]

def gen2Dboundary(inPts = None):
    angle = 0
    id = 0
    CenterX = 10
    CenterY = 30
    r1 = 40
    r2 = 25
    rotation = 0

    points = []
    while (angle < 2 * math.pi):
        point = [r1 * math.cos(angle), r2 * math.sin(angle)]
        if angle > 0 and angle <= math.pi/4:
            point[1] = point[1] - math.sin(angle)*5
        if angle > math.pi/4 and angle <= math.pi/2:
            point[1] = point[1] - (math.cos(angle)*5)
    #         point[0] = point[0] - math.cos(angle)*5
        if angle > math.pi/2 and angle <= 3*math.pi/4:
            point[1] = point[1] - (math.cos(angle)*5)
        if angle > 3*math.pi/4 and angle <= (math.pi):
            point[1] = point[1] + (math.sin(angle)*5)
        x1 = math.cos(rotation) * (point[0]) - math.sin(rotation) * (point[1])
        y1 = math.sin(rotation) * (point[0]) + math.cos(rotation) * (point[1])
        point[0] = x1 + CenterX
        point[1] = y1 + CenterY
        
        
        points.append(point)
        angle = angle + (math.pi / 20)
        id += 1
    points.append(points[0])

    npPoints = np.array(points)
    if inPts is not None:
        inPts.append(inPts[0])
        npPoints = np.array(inPts)

    tck, u = splprep(npPoints.T, u=None, s=0.0, per=1)

    u_new = np.linspace(u.min(), u.max(), 150)
    x_new, y_new = splev(u_new, tck, der=0)

    # plt.plot(x_new, y_new, 'b')
    # plt.axis('scaled')
    # plt.show()

    # x and y are the x-values and y-values of the discrete points sampled along the test case
    # x[i] and y[i] make up the i-th point
    x= x_new[0:-1]
    y= y_new[0:-1]
    orgX = x
    orgY = y
    return orgX, orgY

def curvatureFlow(orgX, orgY):
    x = orgX.copy()
    y = orgY.copy()

    # the way curvature is calculated at a certain point is fitting a circle through that point and a point to the left and right
    # sampling variable determines the interval between the current point and the points to its left and right
    sampling = max(math.floor(len(x) / 40),1)

    # mean curvature flow stops when the shape is near that of an ellipse
    # prevDist is calculating the distance between the points and that of its best fitted ellipse
    prevDist = 0

    # radLen gets length of the major axis of best fitting ellipse, direc is getting the direction of the major axis
    radLen, direc = check_ellipse(x,y)
    coeff = 20

    # diffeoX and Y are keeping track of all the x and y values over every iteration of the curvature flow
    diffeoX = []
    diffeoY = []

    # adding the initial x and y values to diffeoX and diffeoY
    for i in range(0, len(x)):
        diffeoX.append([x[i]])
        diffeoY.append([y[i]])

    # iterating the curvature flow a max 100 times
    for j in range(0,400):
        newPs = []
        for i in range(0, len(x)):
            point0 = (x[i], y[i])
            point1 = (x[i-sampling], y[i-sampling])
            point2 = (x[(i+sampling)%len(x)], y[(i+sampling)%len(x)])
            
            # fitting circle through points
            mid, rad = get_circle(point0, point1, point2)
            
            # getting curvature
            curv = 1/rad
            vect = (mid[0]-point0[0], mid[1]-point0[1])
            mag = math.sqrt(vect[0]**2 + vect[1]**2)
            distance = (curv/mag)*0.1
            disp = (vect[0]*distance, vect[1]*distance)
            
            # moving the current point by curvature value towards the center of the fitted circle
            newP = [disp[0]+point0[0], disp[1]+point0[1]]
            newPs.append(newP)
            
            # keeping track of new x and y values in diffeoX and diffeoY
            x[i] = disp[0]+point0[0]
            y[i] = disp[1]+point0[1]
            diffeoX[i].append(x[i])
            diffeoY[i].append(y[i])
        
        # after every iteration, fit an ellipse to the set of points
        radLen, direc = check_ellipse(x,y)
        
        # calculating the sum of the distances between points on the deformed ellipse and the best fitted ellipse
        # looping through all the points on the deformed ellipse and finding distance from that point to the closest point on the best fit ellipse
        # elDist is the sum of the distances
        elDist = 0
        for i in range(0, len(x)):
            # translating and rotating the points 
            # the "solve" method to find nearest point on ellipse assumes ellipse is at origin and major-axis is along x-axis
            vector_1 = [1, 0]
            vector_2 = [direc[0,0], direc[1,0]]
            unit_vector_1 = vector_1 / np.linalg.norm(vector_1)
            unit_vector_2 = vector_2 / np.linalg.norm(vector_2)
            dot_product = np.dot(unit_vector_1, unit_vector_2)
            angle = np.arccos(dot_product)
            if vector_2[1] > 0:
                angle = -1*angle
        
            avgX = np.mean(x)
            avgY = np.mean(y)
            transX = x[i] - avgX
            transY = y[i] - avgY
            qx = math.cos(angle) * (transX) - math.sin(angle) * (transY)
            qy = math.sin(angle) * (transX) + math.cos(angle) * (transY)
            
            # getting closest point on ellipse
            ptOnEll = solve(radLen[0], radLen[1], [qx, qy])
            
            # summing up distances
            elDist = elDist + math.sqrt((ptOnEll[0]-qx)**2 + (ptOnEll[1]-qy)**2)
        # breaking the curvature flow processes 
        # breaks when the sum of the differences between the points and ellipse is greater than prvious iteration value
        if prevDist != 0:
            # divides by radius length of ellipse to acount for the curvature flow making the shape smaller
            # if we don't divide then the sum of distances would just keep getting smaller
    #         print(prevDist)
            if elDist/radLen[0] > prevDist:
                break
        prevDist = elDist/radLen[0]

    # everything below is to find the closest point on the best fitting ellipse for each of the points on deformed ellipse
    rotPtX = []
    rotPtY = []
    vector_1 = [1, 0]
    vector_2 = [direc[0,0], direc[1,0]]
    unit_vector_1 = vector_1 / np.linalg.norm(vector_1)
    unit_vector_2 = vector_2 / np.linalg.norm(vector_2)
    dot_product = np.dot(unit_vector_1, unit_vector_2)
    angle = np.arccos(dot_product)
    if vector_2[1] > 0:
        angle = -1*angle
    finRotAng = -1*angle

    avgX = np.mean(x)
    avgY = np.mean(y)

    for i in range(0, len(x)):
        transX = x[i] - avgX
        transY = y[i] - avgY
        qx = math.cos(angle) * (transX) - math.sin(angle) * (transY)
        qy = math.sin(angle) * (transX) + math.cos(angle) * (transY)

        ptOnEll = solve(radLen[0], radLen[1], [qx, qy])
        angle2 = angle * -1

        qx2 = math.cos(angle2) * (ptOnEll[0]) - math.sin(angle2) * (ptOnEll[1])
        qy2 = math.sin(angle2) * (ptOnEll[0]) + math.cos(angle2) * (ptOnEll[1])
        qx2 = qx2 + avgX
        qy2 = qy2 + avgY
        
        # rotPtX and rotPtY are the closest points on the best fit ellipse to the points on deformed ellipse
        rotPtX.append(qx2)
        rotPtY.append(qy2)
    
    newPs.append(newPs[0])
    npPoints = np.array(newPs)

    x, y = npPoints.T
    tck, u = splprep(npPoints.T, u=None, s=0.0, per=1)

    u_new = np.linspace(u.min(), u.max(), 100)
    x_new5, y_new5 = splev(u_new, tck, der=0)

    # plt.plot(x_new5, y_new5, 'b')

    # newXs = x_new[0:-1]
    # newYs = y_new[0:-1]

    rotAng = math.atan(direc[1][0]/direc[0][0])
    meanX = np.mean(x)
    meanY = np.mean(y)
    aMaj = radLen[0]
    aMin = radLen[1]


    fitEl = []
    rotRad = 0
    # plotting the best fitting ellipse
    while (rotRad < 2 * math.pi):
        elP = [aMaj * math.cos(rotRad), aMin * math.sin(rotRad)]
        s = math.sin(rotAng)
        c = math.cos(rotAng)
        
        xnew = elP[0] * c - elP[1] * s;
        ynew = elP[0] * s + elP[1] * c;
        
        elP = [xnew+meanX, ynew+meanY]
        
        fitEl.append(elP)
        rotRad = rotRad + (math.pi / 20)
        

    fitEl.append(fitEl[0])

    npfitEl = np.array(fitEl)

    tck, u = splprep(npfitEl.T, u=None, s=0.0, per=1)

    u_new = np.linspace(u.min(), u.max(), 150)
    xnewE, ynewE = splev(u_new, tck, der=0)

    # plt.plot(xnewE, ynewE, 'r')
    # plt.axis('scaled')
    # plt.show()

    return diffeoX, diffeoY, radLen, finRotAng, [meanX, meanY], rotPtX, rotPtY

def ell_srep(diffeoX, diffeoY, radLen, finRotAng, means):
    centX = means[0]
    centY = means[1]
    majA = radLen[0]
    minA = radLen[1] 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
    angle = 0
    id = 0
    radCurv = minA**2/majA
    points = []

    while (angle < 2 * math.pi + (math.pi / 20)):
        point = [majA * math.cos(angle), minA * math.sin(angle)]
        points.append(point)
        angle = angle + (math.pi / 20)
        id += 1
    points[40][1] = points[0][1]
    npPoints = np.array(points)
    x, y = npPoints.T

    tck, u = splprep(npPoints.T, u=None, s=0.0, per=1)
    u_new = np.linspace(u.min(), u.max(), 100)
    x_new, y_new = splev(u_new, tck, der=0)

    # fitting the s-reps
    # the skeleton start the radial curvature in from the "vertex" of the ellipse and ends a radial curv val before the other vertex
    # skeletonXs are the x-values at the ends of the skeleton
    # skeletonYs are the y-values at the ends of the skeleton
    skeletonXs = [-1*majA + radCurv, majA - radCurv]
    # print(centX, centY)
    skeletonYs = [0, 0]

    x1 = math.cos(finRotAng) * (skeletonXs[0]) - math.sin(finRotAng) * (skeletonYs[0])
    y1 = math.sin(finRotAng) * (skeletonXs[0]) + math.cos(finRotAng) * (skeletonYs[0])
    skeletonXs[0] = x1
    skeletonYs[0] = y1

    x2 = math.cos(finRotAng) * (skeletonXs[1]) - math.sin(finRotAng) * (skeletonYs[1])
    y2 = math.sin(finRotAng) * (skeletonXs[1]) + math.cos(finRotAng) * (skeletonYs[1])
    skeletonXs[1] = x2
    skeletonYs[1] = y2

    skeletonXs = [skeletonXs[0]+centX, skeletonXs[1]+centX]
    skeletonYs = [skeletonYs[0]+centY, skeletonYs[1]+centY]
    skelLen = math.sqrt((skeletonXs[1]-skeletonXs[0])**2 + (skeletonYs[1]-skeletonYs[0])**2) /2

    # number of skeletal points
    numSpoke = 20
    theta = -1*math.pi
    sPoints = []
    preTpt = []
    # print(math.pi / numSpoke)
    # sampling points on the skeleton
    while (theta < (math.pi / numSpoke) - 0.000000001):
        sPoint = [skelLen * (math.cos(theta)+1)+ -1*majA + radCurv, 0]
        preTpt.append([sPoint[0], sPoint[1]])
        sPX = math.cos(finRotAng) * (sPoint[0]) - math.sin(finRotAng) * (sPoint[1])
        sPY = math.sin(finRotAng) * (sPoint[0]) + math.cos(finRotAng) * (sPoint[1])
        sPoint[0] = sPX
        sPoint[1] = sPY
        sPoint = [sPoint[0]+centX, sPoint[1]+centY]
        # print(theta)
        #sPoints are the skeletal points
        sPoints.append(sPoint)
        theta = theta + (math.pi / numSpoke)


        
    sPoints = np.array(sPoints)
    sx, sy = sPoints.T
    interSkelX = np.linspace(sx[0], sx[-1], num=100)
    interSkelY = np.linspace(sy[0], sy[-1], num=100)
    # print(sx)
    b_pointsU = []
    b_pointsD = []

    # looping through the skeletal points
    for i in range(0,len(sx)):
        # finding the nearest point on ellipse from the skeletal point
        # solve returns the boundary point of the up-spoke at a skeletal point
        # the boundary point of the down spoke will have the same x-val and -1*y-val
        b_pt = solve(majA, minA,preTpt[i])
        if i == 0:
            b_pt = [-1*majA, 0]
        if i == len(sx)-1:
            b_pt = [majA, 0]
        oldP = b_pt.copy()
        
        # rotating and translating the points back to fit the ellipse
        downBx = math.cos(finRotAng) * (b_pt[0]) - math.sin(finRotAng) * (-1*b_pt[1])
        downBy = math.sin(finRotAng) * (b_pt[0]) + math.cos(finRotAng) * (-1*b_pt[1])
        b_ptX = math.cos(finRotAng) * (b_pt[0]) - math.sin(finRotAng) * (b_pt[1])
        b_ptY = math.sin(finRotAng) * (b_pt[0]) + math.cos(finRotAng) * (b_pt[1])
        b_pt[0] = b_ptX
        b_pt[1] = b_ptY
        b_pt[0] = b_pt[0]+centX
        b_pt[1] = b_pt[1]+centY
        downBx = downBx+centX
        downBy = downBy+centY
        
        # saving all the up boundary points and the down points in b_pointsU and b_pointsD
        b_pointsU.append([b_pt[0], b_pt[1]])
        if i != 0 and i != len(sx)-1:
            b_pointsD.append([downBx, downBy])
        x_vals = [sx[i], b_pt[0]]
        x_vals2 = [sx[i], downBx]
        y_vals = [sy[i], b_pt[1]]
        y_vals2 = [sy[i], downBy]
        # plt.plot(x_vals, y_vals, 'k')
        if i != 0 and i != len(sx)-1:
            pass
            # plt.plot(x_vals2, y_vals2, 'k')
    # print(tSy)
    # plt.plot(centX, centY, 'ro')

    # plt.plot(skeletonXs, skeletonYs,'r')
    # plt.axis('scaled')
    # plt.show()
    return b_pointsU, b_pointsD, sPoints, interSkelX, interSkelY

def inv_curv(b_pointsU, b_pointsD, sPoints, diffeoX, diffeoY, rotPtX, rotPtY, interSkelX, interSkelY):
    invPtsU = b_pointsU.copy()
    invPtsD = b_pointsD.copy()
    skelPt = sPoints.copy()

    # looping through all the iterations of the curvature flow (we are doing the inverse of the flow)
    for j in range(0,len(diffeoX[0])):
        movingX = []
        movingY = []
        source_pts = vtk.vtkPoints()
        target_pts = vtk.vtkPoints()
        for i in range(0,len(rotPtX)):
            pt = []
            if j == 0:  
                pt = [rotPtX[i], rotPtY[i], 0]
            else:
                pt = [diffeoX[i][-1*j], diffeoY[i][-1*j], 0]
            
            # source pts for TPS starts off as best fit ellipse
            # next iteration it will be the last diffeomorphism
            # keep going until it is the second diffeomorphism
            source_pts.InsertNextPoint(pt)

            pt = [diffeoX[i][-1-j], diffeoY[i][-1-j], 0]
            movingX.append(diffeoX[i][-1-j])
            movingY.append(diffeoY[i][-1-j])
            
            # target points is the diffeomorphism before the diffeomorphism of the source points
            # starts as last and goes to first
            target_pts.InsertNextPoint(pt)
        source_pts.Modified()
        target_pts.Modified()

        ### Interpolate deformation with thin-plate-spline
        tps = vtk.vtkThinPlateSplineTransform()
        tps.SetSourceLandmarks(source_pts)
        tps.SetTargetLandmarks(target_pts)
        tps.SetBasisToR()
        tps.Modified()
        sXs = []
        sYs = []
        
        
        # finding interpolated skel points
        for i in range(0, len(interSkelX)):
            intSkelPt = tps.TransformPoint([interSkelX[i], interSkelY[i], 0])
            interSkelX[i] = intSkelPt[0]
            interSkelY[i] = intSkelPt[1]
        
        for i in range(0, len(sPoints)):
            sp = [skelPt[i][0], skelPt[i][1], 0]
            up = [invPtsU[i][0], invPtsU[i][1], 0]
            if i != 0 and i != len(sPoints)-1:
                down = [invPtsD[i-1][0], invPtsD[i-1][1], 0]
            
            # finding new skeletal points
            sk_pt = tps.TransformPoint(sp)
            skelPt[i][0] = sk_pt[0]
            skelPt[i][1] = sk_pt[1]
            sXs.append(sk_pt[0])
            sYs.append(sk_pt[1])
            
            # finding new boundary points for the up spokes
            newUp = tps.TransformPoint(up)
            invPtsU[i][0] = newUp[0]
            invPtsU[i][1] = newUp[1]
            
            # new boundary points of bottom spokes
            if i != 0 and i != len(sPoints)-1:
                newDown = tps.TransformPoint(down)
                invPtsD[i-1][0] = newDown[0]
                invPtsD[i-1][1] = newDown[1]
            
            # plotting the final spokes
            if j == len(diffeoX[0])-1:
                xUp = [sk_pt[0], newUp[0]]
                yUp = [sk_pt[1], newUp[1]]
                # plt.plot(xUp, yUp,'k')
                if i != 0 and i != len(sPoints)-1:
                    xDown = [sk_pt[0], newDown[0]]
                    yDown = [sk_pt[1], newDown[1]]
                    # plt.plot(xDown, yDown,'k')
            
    movingX.append(movingX[0])
    movingY.append(movingY[0])
    # plt.plot(interSkelX, interSkelY,'r')
    # plt.plot(movingX, movingY, 'b')
    # plt.show()
    return sXs, sYs, invPtsU, invPtsD

def refine(sXs, sYs, invPtsU, invPtsD, movingX, movingY):
    radii_arr = []
    dirVecs = []
    movingX.append(movingX[0])
    movingY.append(movingY[0])
    for i in range(0,len(sXs)):
        radmag = math.sqrt((sXs[i]-invPtsU[i][0])**2 + (sYs[i]-invPtsU[i][1])**2)
        vect = [(invPtsU[i][0] - sXs[i])/radmag, (invPtsU[i][1] - sYs[i])/radmag]
        # adding unit vector directions of the up spokes
        dirVecs.append(vect)
        # adding lengths of the up spokes
        radii_arr.append(radmag)
        
        if i != 0 and i != len(sXs)-1:
            radmag = math.sqrt((sXs[i]-invPtsD[i-1][0])**2 + (sYs[i]-invPtsD[i-1][1])**2)
            vect = [(invPtsD[i-1][0] - sXs[i])/radmag, (invPtsD[i-1][1] - sYs[i])/radmag]
            # adding unit vector directions of the down spokes
            dirVecs.append(vect)
            # adding lengths of the down spokes
            radii_arr.append(radmag)
    plotx = np.array(movingX)
    ploty = np.array(movingY)
    tups = np.array((plotx,ploty)).T
    shapePts = list(map(lambda x, y:(x,y), movingX, movingY))
    polyPlot = Polygon(shapePts)
    def obj_func(radii, grad=None):
        total_loss = 0
        for i, radius in enumerate(radii):
            ind = math.ceil((i)/2)
            base_pt   = [sXs[ind], sYs[ind]]
            direc = [dirVecs[i][0] * radius, dirVecs[i][1] * radius]
            bdry_pt   = [base_pt[0]+direc[0], base_pt[1]+direc[1]]
            point = Point(bdry_pt[0],bdry_pt[1])
            dist = polyPlot.exterior.distance(point)
            total_loss += dist ** 2
        return total_loss

    import nlopt
    opt = nlopt.opt(nlopt.LN_NEWUOA, len(radii_arr))
    opt.set_min_objective(obj_func)
    opt.set_maxeval(450)
    opLen = opt.optimize(radii_arr)

    boundPts = []
    for i in range(0, len(sXs)):
        upx = sXs[i] + dirVecs[2*i-1][0]*opLen[2*i-1]
        upy = sYs[i] + dirVecs[2*i-1][1]*opLen[2*i-1]
        xUp = [sXs[i], upx]
        yUp = [sYs[i], upy]
        if i == 0:
            upx = sXs[i] + dirVecs[0][0]*opLen[0]
            upy = sYs[i] + dirVecs[0][1]*opLen[0]
            xUp = [sXs[i], upx]
            yUp = [sYs[i], upy]
        
        boundPts.append([upx, upy])
    
        if i != 0 and i != len(sXs)-1:
            downx = sXs[i] + dirVecs[2*i][0]*opLen[2*i]
            downy = sYs[i] + dirVecs[2*i][1]*opLen[2*i]
            boundPts.append([downx, downy])
            xDown = [sXs[i], downx]
            yDown = [sYs[i], downy]
            # plt.plot(xDown, yDown,'k')
        
        
    #     plt.plot(xUp, yUp,'k')
    # plt.plot(movingX, movingY, 'b')
    # plt.plot(sXs, sYs,'r')
    # plt.axis('scaled')
    # plt.show()
    
    angle = -1*math.pi/2
    min_loss = 1
    min_ang = 0
    min_len = 0
    def obj_func(angles, grad=None):
        total_loss = 0
        for i, angle in enumerate(angles):
            ind = math.ceil((i)/2)
            base_pt   = [sXs[ind], sYs[ind]]
            bdry_pt   = [dirVecs[i][0]*opLen[i], dirVecs[i][1]*opLen[i]]

            s = math.sin(angle)
            c = math.cos(angle)

            xnew = bdry_pt[0] * c - bdry_pt[1] * s
            ynew = bdry_pt[0] * s + bdry_pt[1] * c

            bdry_pt[0] = xnew + sXs[ind]
            bdry_pt[1] = ynew + sYs[ind]
            preFin = bdry_pt.copy()

            vec = [bdry_pt[0]-base_pt[0],bdry_pt[1]-base_pt[1]]
            unitV = vec / np.linalg.norm(vec)

            point = Point(bdry_pt[0],bdry_pt[1])
            while polyPlot.exterior.distance(point) > 0.00000001:

                point = Point(bdry_pt[0],bdry_pt[1])
                dispV = unitV * polyPlot.exterior.distance(point)
                if polyPlot.distance(point) == 0:
                    bdry_pt[0] = bdry_pt[0] + dispV[0]
                    bdry_pt[1] = bdry_pt[1] + dispV[1]

                else:
                    bdry_pt[0] = bdry_pt[0] - dispV[0]
                    bdry_pt[1] = bdry_pt[1] - dispV[1]

            sLength = math.hypot(bdry_pt[0]-sXs[ind], bdry_pt[1]-sYs[ind])

            shortDist = -1
            ptInd = -1
            for j in range(0, len(movingX)):
                ptDist = math.hypot(bdry_pt[0]-movingX[j], bdry_pt[1]-movingY[j])
                if ptDist < shortDist or shortDist == -1:
                    shortDist = ptDist
                    ptInd = j

            ldist = math.hypot(bdry_pt[0]-movingX[ptInd-1], bdry_pt[1]-movingY[ptInd-1])
            if ptInd == 0:
                ldist = math.hypot(bdry_pt[0]-movingX[ptInd-2], bdry_pt[1]-movingY[ptInd-2])
            rdist = math.hypot(bdry_pt[0]-movingX[(ptInd+1)%len(movingX)], bdry_pt[1]-movingY[(ptInd+1)%len(movingX)])
            if ptInd == len(movingX)-1:
                rdist = math.hypot(bdry_pt[0]-movingX[1], bdry_pt[1]-movingY[1])
            left = False
            if ldist < rdist:
                left = True

            vector_1 = [bdry_pt[0]-base_pt[0], bdry_pt[1]-base_pt[1]]
            if left:
                vector_2 = [movingX[ptInd] - movingX[ptInd-1], movingY[ptInd] - movingY[ptInd-1]]
                if ptInd == 0:
                    vector_2 = [movingX[ptInd] - movingX[ptInd-2], movingY[ptInd] - movingY[ptInd-2]]
            else:
                vector_2 = [movingX[(ptInd+1)%len(movingX)] - movingX[ptInd], movingY[(ptInd+1)%len(movingX)] - movingY[ptInd]]
                if ptInd == len(movingX)-1:
                    vector_2 = [movingX[1] - movingX[ptInd], movingY[1] - movingY[ptInd]]
            unit_vector_1 = vector_1 / np.linalg.norm(vector_1)
            unit_vector_2 = vector_2 / np.linalg.norm(vector_2)
            dot_product = np.dot(unit_vector_1, unit_vector_2)
            ang = np.arccos(dot_product)
            total_loss = total_loss + (1 - abs(math.sin(ang)))
            if angle < -1*math.pi/16 or angle > math.pi/16:
                total_loss = total_loss + 100
        return total_loss
    
    import nlopt
    angles = np.zeros((len(sXs)-1)*2)
    opt = nlopt.opt(nlopt.LN_NEWUOA, len(angles))
    opt.set_min_objective(obj_func)
    opt.set_maxeval(600)
    finAng = opt.optimize(angles)

    finBdryPts = []
    for i, angle in enumerate(finAng):
        ind = math.ceil((i)/2)
        base_pt   = [sXs[ind], sYs[ind]]
        bdry_pt   = [dirVecs[i][0]*opLen[i], dirVecs[i][1]*opLen[i]]

        s = math.sin(angle)
        c = math.cos(angle)

        xnew = bdry_pt[0] * c - bdry_pt[1] * s
        ynew = bdry_pt[0] * s + bdry_pt[1] * c

        bdry_pt[0] = xnew + sXs[ind]
        bdry_pt[1] = ynew + sYs[ind]

        vec = [bdry_pt[0]-base_pt[0],bdry_pt[1]-base_pt[1]]
        unitV = vec / np.linalg.norm(vec)

        point = Point(bdry_pt[0],bdry_pt[1])
        while polyPlot.exterior.distance(point) > 0.000000001:
            point = Point(bdry_pt[0],bdry_pt[1])
            dispV = unitV * polyPlot.exterior.distance(point)

            if polyPlot.distance(point) == 0:
                bdry_pt[0] = bdry_pt[0] + dispV[0]
                bdry_pt[1] = bdry_pt[1] + dispV[1]
            else:
                bdry_pt[0] = bdry_pt[0] - dispV[0]
                bdry_pt[1] = bdry_pt[1] - dispV[1]
        spokeXs = [base_pt[0], bdry_pt[0]]
        spokeYs = [base_pt[1], bdry_pt[1]]
        finBdryPts.append(bdry_pt)
        # plt.plot(spokeXs, spokeYs,'k')
        
        shortDist = -1
        ptInd = -1
        for j in range(0, len(movingX)):
            ptDist = math.hypot(bdry_pt[0]-movingX[j], bdry_pt[1]-movingY[j])
            if ptDist < shortDist or shortDist == -1:
                shortDist = ptDist
                ptInd = j

        ldist = math.hypot(bdry_pt[0]-movingX[ptInd-1], bdry_pt[1]-movingY[ptInd-1])
        if ptInd == 0:
            ldist = math.hypot(bdry_pt[0]-movingX[ptInd-2], bdry_pt[1]-movingY[ptInd-2])
        rdist = math.hypot(bdry_pt[0]-movingX[(ptInd+1)%len(movingX)], bdry_pt[1]-movingY[(ptInd+1)%len(movingX)])
        if ptInd == len(movingX)-1:
            rdist = math.hypot(bdry_pt[0]-movingX[1], bdry_pt[1]-movingY[1])
        left = False
        if ldist < rdist:
            left = True

        vector_1 = [bdry_pt[0]-base_pt[0], bdry_pt[1]-base_pt[1]]
        if left:
            vector_2 = [movingX[ptInd] - movingX[ptInd-1], movingY[ptInd] - movingY[ptInd-1]]
            if ptInd == 0:
                vector_2 = [movingX[ptInd] - movingX[ptInd-2], movingY[ptInd] - movingY[ptInd-2]]
        else:
            vector_2 = [movingX[(ptInd+1)%len(movingX)] - movingX[ptInd], movingY[(ptInd+1)%len(movingX)] - movingY[ptInd]]
            if ptInd == len(movingX)-1:
                vector_2 = [movingX[1] - movingX[ptInd], movingY[1] - movingY[ptInd]]

        unit_vector_1 = vector_1 / np.linalg.norm(vector_1)
        unit_vector_2 = vector_2 / np.linalg.norm(vector_2)
        dot_product = np.dot(unit_vector_1, unit_vector_2)
        ang = np.arccos(dot_product)
        # print(ang)
    # plt.plot(movingX, movingY, 'b')
    # plotX = sXs[0:-1]
    # plotX.append(sXs[-1]+3.5)
    # plotX[0] = plotX[0]-3.5
    # plotY = sYs
    # plt.plot(sXs, sYs,'r')

    base_pt1   = [sXs[0], sYs[0]]
    bdry_pt1   = [sXs[0]+dirVecs[0][0]*opLen[0], sYs[0]+dirVecs[0][1]*opLen[0]]

    base_pt2   = [sXs[-1], sYs[-1]]
    bdry_pt2   = [sXs[-1]+dirVecs[-1][0]*opLen[-1], sYs[-1]+dirVecs[-1][1]*opLen[-1]]

    spokeXs = [base_pt1[0], bdry_pt1[0]]
    spokeYs = [base_pt1[1], bdry_pt1[1]]
    # finBdryPts.append(bdry_pt1)
    spokeXs = [base_pt2[0], bdry_pt2[0]]
    spokeYs = [base_pt2[1], bdry_pt2[1]]
    # plt.axis('scaled')
    # plt.show()
    return finBdryPts



def run_sim(points = None):
    if points is not None:
        x, y = gen2Dboundary(points)
    else:
        x, y = gen2Dboundary()
    diffeoX, diffeoY, radLen, finRotAng, means, rotPtX, rotPtY = curvatureFlow(x,y)
    b_pointsU, b_pointsD, sPoints, interSkelX, interSkelY = ell_srep(diffeoX, diffeoY, radLen, finRotAng, means)
    sXs, sYs, invPtsU, invPtsD = inv_curv(b_pointsU, b_pointsD, sPoints, diffeoX, diffeoY, rotPtX, rotPtY, interSkelX, interSkelY)
    b_pts = refine(sXs, sYs, invPtsU, invPtsD, [i[0] for i in diffeoX], [i[0] for i in diffeoY])
    return sXs, sYs, b_pts
# run_sim()
