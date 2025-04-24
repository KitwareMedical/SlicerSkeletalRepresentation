import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import math
import numpy as np
from scipy.spatial import ConvexHull, convex_hull_plot_2d
from scipy.interpolate import splprep, splev
from shanapy.models import fit2D as fit2d
import vtk

def curvedSrep(sharedPts):

    # These constants are to create random data for the sake of this example
    # N_POINTS = 20
    # TARGET_X_SLOPE = 1
    # TARGET_y_SLOPE = 5
    # TARGET_OFFSET  = 1
    # EXTENTS = 5
    # NOISE = 5

    # # Create random data to test code if you dont have data.
    # # In your solution, you would provide your own xs, ys, and zs data.
    # xs = [np.random.uniform(2*EXTENTS)-EXTENTS for i in range(N_POINTS)]
    # ys = [np.random.uniform(2*EXTENTS)-EXTENTS for i in range(N_POINTS)]
    # zs = []
    # for i in range(N_POINTS):
    #     zs.append(xs[i]*TARGET_X_SLOPE + \
    #             ys[i]*TARGET_y_SLOPE + \
    #             TARGET_OFFSET + np.random.normal(scale=NOISE))


    N_POINTS = sharedPts.GetNumberOfPoints()
    xs = []
    ys = []
    zs = []
    for i in range(0, sharedPts.GetNumberOfPoints()):
        sbPt = sharedPts.GetPoint(i)
        xs.append(sbPt[0])
        ys.append(sbPt[1])
        zs.append(sbPt[2])

    # plot raw data
    plt.figure()
    ax = plt.subplot(111, projection='3d')
    # ax.scatter(xs, ys, zs, color='r')


    # fitting a plane through the group of points
    tmp_A = []
    tmp_b = []
    for i in range(len(xs)):
        tmp_A.append([xs[i], ys[i], 1])
        tmp_b.append(zs[i])
    b = np.matrix(tmp_b).T
    A = np.matrix(tmp_A)

    fit = (A.T * A).I * A.T * b
    errors = b - A * fit
    residual = np.linalg.norm(errors)

    
    # equation of plane fit[0]*x + fit[1]*y + fit[2] = z
    print("solution: %f x + %f y + %f = z" % (fit[0], fit[1], fit[2]))
    # print("errors: \n", errors)
    print("residual:", residual)
    print([fit[0], fit[1], fit[2]])

    # z intercept of the plane
    zInt = fit[2,0]

    # getting normal vector to the plane
    vec1 = np.array([0, 1, fit[2,0] + fit[1,0]]) - np.array([0, 0, fit[2,0]])
    vec2 = np.array([1, 0, fit[2,0] + fit[0,0]]) - np.array([0, 0, fit[2,0]])
    cross = np.cross(vec1, vec2)
    norm = np.array([fit[0,0], fit[1,0], -1])
    print(cross)
    # ax.plot([0, fit[0,0]], [0,fit[1,0]], [fit[2,0], fit[2,0]+-1], 'b')
    # ax.plot([0, 1], [0,0], [fit[2,0], fit[2,0]+fit[0,0]], 'b')
    # ax.plot([0, cross[0]], [0, cross[1]], [fit[2,0], cross[2]], 'b')

    c_hat = cross / np.linalg.norm(norm)
    print(np.dot(c_hat, vec2))

    A = np.array([[c_hat[2]+(c_hat[1]**2)*(1-c_hat[2]), c_hat[1]*-1*c_hat[0]*(1-c_hat[2]), -1*c_hat[0]*math.hypot(c_hat[0], c_hat[1])], 
                [c_hat[1]*-1*c_hat[0]*(1-c_hat[2]), c_hat[2]+(c_hat[0]**2)*(1-c_hat[2]), -1*c_hat[1]*math.hypot(c_hat[0], c_hat[1])],
                [c_hat[0]*math.hypot(c_hat[0], c_hat[1]), c_hat[1]*math.hypot(c_hat[0], c_hat[1]), c_hat[2]]])
    print("z-axis")
    print(np.dot(A, c_hat))

    # ax.scatter(0, 0, 0, color='r')

    # projecting all the initial points on the best fitting plane
    planePts = []
    for i in range(0, len(xs)):
        disp = np.array([xs[i], ys[i], zs[i]]) - np.array([0, 0, fit[2,0]])
        dotp = np.dot(c_hat, disp)
        newPt = np.array([xs[i], ys[i], zs[i]]) - (c_hat*dotp)
        mat2 = np.array([[newPt[0]], 
                [newPt[1]],
                [newPt[2] - fit[2,0]]])
        planePts.append([newPt[0], newPt[1], newPt[2] - fit[2,0]])


    # rotating the plane into the XY plane
    trip = np.array([xs[0], ys[0], zs[0]])
    mat1 = np.array([[1, 1, 1], 
                [1, 1, 1],
                [1, 1, 1]])

    xzVec = np.array([0,1])
    xz_norm = np.array([c_hat[0], c_hat[2]])

    unit_vector_1 = xzVec / np.linalg.norm(xzVec)
    unit_vector_2 = xz_norm / np.linalg.norm(xz_norm)
    dot_product = np.dot(unit_vector_1, unit_vector_2)
    angleXZ = np.arccos(dot_product)

    s = math.sin(angleXZ);
    c = math.cos(angleXZ);
    xCoef = xz_norm[0] * c - xz_norm[1] * s;
    zCoef = xz_norm[0] * s + xz_norm[1] * c;

    if abs(xCoef-0) > 10**-12:
        angleXZ = 2*math.pi -angleXZ
        s = math.sin(angleXZ);
        c = math.cos(angleXZ);
        xCoef = xz_norm[0] * c - xz_norm[1] * s;
        zCoef = xz_norm[0] * s + xz_norm[1] * c;

    print("test")
    print([xz_norm])
    print([xCoef, zCoef])
    print(angleXZ)

    yzVec = np.array([0,1])
    yz_norm = np.array([c_hat[1], zCoef])

    unit_vector_1 = yzVec / np.linalg.norm(yzVec)
    unit_vector_2 = yz_norm / np.linalg.norm(yz_norm)
    dot_product = np.dot(unit_vector_1, unit_vector_2)
    angleYZ = np.arccos(dot_product)

    newYZnorm = [c_hat[1], zCoef]

    s = math.sin(angleYZ);
    c = math.cos(angleYZ);
    yCoef = newYZnorm[0] * c - newYZnorm[1] * s;
    zCoef = newYZnorm[0] * s + newYZnorm[1] * c;

    if abs(yCoef-0) > 10**-12:
        angleYZ = 2*math.pi -angleYZ
        s = math.sin(angleYZ);
        c = math.cos(angleYZ);
        yCoef = newYZnorm[0] * c - newYZnorm[1] * s;
        zCoef = newYZnorm[0] * s + newYZnorm[1] * c;


    print("test2")
    print(newYZnorm)
    print([yCoef, zCoef])
    print(angleXZ)


    # rotating the points on the best fitting plane into the XY plane
    xyPts = []
    for i in range(0, len(planePts)):
        s = math.sin(angleXZ);
        c = math.cos(angleXZ);

        xnew = planePts[i][0] * c - planePts[i][2] * s;
        znew = planePts[i][0] * s + planePts[i][2] * c;

        newYZpt = [planePts[i][1], znew]

        s = math.sin(angleYZ);
        c = math.cos(angleYZ);

        ynew = newYZpt[0] * c - newYZpt[1] * s;
        znew = newYZpt[0] * s + newYZpt[1] * c;

        xyPts.append([xnew, ynew])

    # getting the convex hull of the group of points
    xyPts = np.array(xyPts)
    hull = ConvexHull(xyPts)
    verts = hull.vertices
    hullPts = []
    for i in range(0, len(verts)):
        # ax.plot([xyPts[verts[i], 0], xyPts[verts[(i+1)%len(verts)], 0]], [xyPts[verts[i], 1], xyPts[verts[(i+1)%len(verts)], 1]], [0, 0], 'b')
        hullPts.append([xyPts[verts[i], 0], xyPts[verts[i], 1]])

    # fitting 2D flat s-rep to the convex hull of the points on the XY plane
    sXs, sYs, b_pts = fit2d.run_sim(hullPts)

    # hullPts.append(hullPts)
    npPoints = np.array(hullPts)

    tck, u = splprep(npPoints.T, u=None, s=0.0, per=1)

    u_new = np.linspace(u.min(), u.max(), 150)
    x_new, y_new = splev(u_new, tck, der=0)

    # ax.scatter(x_new, y_new, color='b')
    print(b_pts[0])

    # sampling points along the spokes of the 2D s-rep
    sampBpts = []
    numSamp = 4
    for i in range(0, len(sXs)):
        if i == 0:
            length = math.hypot(b_pts[0][0] - sXs[0], b_pts[0][1] - sYs[0])
            unitV = [(b_pts[0][0] - sXs[0])/length, (b_pts[0][1] - sYs[0])/length]
            spokePts = []
            for j in range(1, numSamp+1):
                spokePts.append([sXs[0] + unitV[0]*length*j/numSamp, sYs[0] + unitV[1]*length*j/numSamp])
            sampBpts.append(spokePts)
        elif i == len(sXs) - 1:
            length = math.hypot(b_pts[-1][0] - sXs[-1], b_pts[-1][1] - sYs[-1])
            unitV = [(b_pts[-1][0] - sXs[-1])/length, (b_pts[-1][1] - sYs[-1])/length]
            spokePts = []
            for j in range(1, numSamp+1):
                spokePts.append([sXs[-1] + unitV[0]*length*j/numSamp, sYs[-1] + unitV[1]*length*j/numSamp])
            sampBpts.append(spokePts)
        else:
            length = math.hypot(b_pts[2*i-1][0] - sXs[i], b_pts[2*i-1][1] - sYs[i])
            unitV = [(b_pts[2*i-1][0] - sXs[i])/length, (b_pts[2*i-1][1] - sYs[i])/length]
            spokePts = []
            for j in range(1, numSamp+1):
                spokePts.append([sXs[i] + unitV[0]*length*j/numSamp, sYs[i] + unitV[1]*length*j/numSamp])
            sampBpts.append(spokePts)

            length = math.hypot(b_pts[2*i][0] - sXs[i], b_pts[2*i-1][1] - sYs[i])
            unitV = [(b_pts[2*i][0] - sXs[i])/length, (b_pts[2*i][1] - sYs[i])/length]
            spokePts = []
            for j in range(1, numSamp+1):
                spokePts.append([sXs[i] + unitV[0]*length*j/numSamp, sYs[i] + unitV[1]*length*j/numSamp])
            sampBpts.append(spokePts)

    # ax.scatter(sXs, sYs, color='r')
    # for i in range(0, len(sXs)):
    #     if i == 0:
    #         ax.plot([sXs[0], b_pts[0][0]], [sYs[0], b_pts[0][1]], [0, 0], 'k')
    #         for j in range(0, numSamp):
    #             ax.scatter(sampBpts[i][j][0], sampBpts[i][j][1], color='k')
    #     elif i == len(sXs) -1:
    #         ax.plot([sXs[-1], b_pts[-1][0]], [sYs[-1], b_pts[-1][1]], [0, 0], 'k')
    #         for j in range(0, numSamp):
    #             ax.scatter(sampBpts[-1][j][0], sampBpts[-1][j][1], color='k')
    #     else:
    #         ax.plot([sXs[i], b_pts[2*i - 1][0]], [sYs[i], b_pts[2*i-1][1]], [0, 0], 'k')
    #         ax.plot([sXs[i], b_pts[2*i][0]], [sYs[i], b_pts[2*i][1]], [0, 0], 'k')
    #         for j in range(0, numSamp):
    #             ax.scatter(sampBpts[i*2-1][j][0], sampBpts[i*2-1][j][1], color='k')
    #             ax.scatter(sampBpts[i*2][j][0], sampBpts[i*2][j][1], color='k')

    # rotating the points back into the best fitting plane
    rXs = []
    rYs = []
    rZs = []
    rSamPts = []
    for i in range(0, len(sXs)):
        inAngYZ = angleYZ*-1
        s = math.sin(inAngYZ)
        c = math.cos(inAngYZ)
        inYZ = [sYs[i], 0]
        ynew = inYZ[0] * c - inYZ[1] * s
        znew = inYZ[0] * s + inYZ[1] * c

        inXZ = [sXs[i], znew]
        inAngXZ = angleXZ*-1
        s = math.sin(inAngXZ)
        c = math.cos(inAngXZ)
        xnew = inXZ[0] * c - inXZ[1] * s
        znew = inXZ[0] * s + inXZ[1] * c

        rXs.append(xnew)
        rYs.append(ynew)
        rZs.append(znew+zInt)
        # ax.scatter(xnew, ynew, znew+zInt, color='k')

        if i == 0:
            rSpoke = []
            for j in range(0, numSamp):
                s = math.sin(inAngYZ)
                c = math.cos(inAngYZ)
                inYZ = [sampBpts[i][j][1], 0]
                ynew = inYZ[0] * c - inYZ[1] * s
                znew = inYZ[0] * s + inYZ[1] * c

                inXZ = [sampBpts[i][j][0], znew]
                s = math.sin(inAngXZ)
                c = math.cos(inAngXZ)
                xnew = inXZ[0] * c - inXZ[1] * s
                znew = inXZ[0] * s + inXZ[1] * c
                rSpoke.append([xnew, ynew, znew+zInt])
                # ax.scatter(xnew, ynew, znew+zInt, color='k')
            rSamPts.append(rSpoke)
            # ax.plot([rXs[i], rSpoke[-1][0]], [rYs[i], rSpoke[-1][1]], [rZs[i], rSpoke[-1][2]], 'r')
        elif i == len(sXs)-1:
            rSpoke = []
            for j in range(0, numSamp):
                s = math.sin(inAngYZ)
                c = math.cos(inAngYZ)
                inYZ = [sampBpts[-1][j][1], 0]
                ynew = inYZ[0] * c - inYZ[1] * s
                znew = inYZ[0] * s + inYZ[1] * c

                inXZ = [sampBpts[-1][j][0], znew]
                s = math.sin(inAngXZ)
                c = math.cos(inAngXZ)
                xnew = inXZ[0] * c - inXZ[1] * s
                znew = inXZ[0] * s + inXZ[1] * c
                rSpoke.append([xnew, ynew, znew+zInt])
                # ax.scatter(xnew, ynew, znew+zInt, color='k')
            rSamPts.append(rSpoke)
            # ax.plot([rXs[i], rSpoke[-1][0]], [rYs[i], rSpoke[-1][1]], [rZs[i], rSpoke[-1][2]], 'r')
        else:
            rSpoke = []
            for j in range(0, numSamp):
                s = math.sin(inAngYZ)
                c = math.cos(inAngYZ)
                inYZ = [sampBpts[2*i - 1][j][1], 0]
                ynew = inYZ[0] * c - inYZ[1] * s
                znew = inYZ[0] * s + inYZ[1] * c

                inXZ = [sampBpts[2*i - 1][j][0], znew]
                s = math.sin(inAngXZ)
                c = math.cos(inAngXZ)
                xnew = inXZ[0] * c - inXZ[1] * s
                znew = inXZ[0] * s + inXZ[1] * c
                rSpoke.append([xnew, ynew, znew+zInt])
                # ax.scatter(xnew, ynew, znew+zInt, color='k')
            rSamPts.append(rSpoke)
            # ax.plot([rXs[i], rSpoke[-1][0]], [rYs[i], rSpoke[-1][1]], [rZs[i], rSpoke[-1][2]], 'r')

            rSpoke = []
            for j in range(0, numSamp):
                s = math.sin(inAngYZ)
                c = math.cos(inAngYZ)
                inYZ = [sampBpts[2*i][j][1], 0]
                ynew = inYZ[0] * c - inYZ[1] * s
                znew = inYZ[0] * s + inYZ[1] * c

                inXZ = [sampBpts[2*i][j][0], znew]
                s = math.sin(inAngXZ)
                c = math.cos(inAngXZ)
                xnew = inXZ[0] * c - inXZ[1] * s
                znew = inXZ[0] * s + inXZ[1] * c
                rSpoke.append([xnew, ynew, znew+zInt])
                # ax.scatter(xnew, ynew, znew+zInt, color='k')
            rSamPts.append(rSpoke)
            # ax.plot([rXs[i], rSpoke[-1][0]], [rYs[i], rSpoke[-1][1]], [rZs[i], rSpoke[-1][2]], 'r')
    # ax.plot(rXs, rYs, rZs, 'r')

    source_pts = vtk.vtkPoints()
    target_pts = vtk.vtkPoints()

    for i in range(0,N_POINTS):
        source_pts.InsertNextPoint(planePts[i])
        target_pts.InsertNextPoint(xs[i], ys[i], zs[i])
    source_pts.Modified()
    target_pts.Modified()

    ### Interpolate deformation with thin-plate-spline
    tps = vtk.vtkThinPlateSplineTransform()
    tps.SetSourceLandmarks(source_pts)
    tps.SetTargetLandmarks(target_pts)
    tps.SetBasisToR()
    tps.Modified()

    iXs = []
    iYs = []
    iZs = []
    iSamPts = []
    srepPts = vtk.vtkPoints()

    # interpolating the points from the best fitting plane into their original space
    for i in range(0, len(sXs)):
        intPt = tps.TransformPoint([rXs[i], rYs[i], rZs[i]])
        iXs.append(intPt[0])
        iYs.append(intPt[1])
        iZs.append(intPt[2])

        if i == 0:
            iSpoke = []
            srepPts.InsertNextPoint([iXs[i], iYs[i], iZs[i]])
            for j in range(0, numSamp):
                intPt = tps.TransformPoint(rSamPts[i][j])
                iSpoke.append([intPt[0], intPt[1], intPt[2]])
                srepPts.InsertNextPoint([intPt[0], intPt[1], intPt[2]])
                # ax.scatter(intPt[0], intPt[1], intPt[2], color='k')
                if j == 0:
                    pass
                    # ax.plot([iXs[i], iSpoke[0][0]], [iYs[i], iSpoke[0][1]], [iZs[i], iSpoke[0][2]], 'r')
                else:
                    pass
                    # ax.plot([iSpoke[j-1][0], iSpoke[j][0]], [iSpoke[j-1][1], iSpoke[j][1]], [iSpoke[j-1][2], iSpoke[j][2]], 'r')
            iSamPts.append(iSpoke)
        elif i == len(sXs)-1:
            iSpoke = []
            srepPts.InsertNextPoint([iXs[i], iYs[i], iZs[i]])
            for j in range(0, numSamp):
                intPt = tps.TransformPoint(rSamPts[-1][j])
                iSpoke.append([intPt[0], intPt[1], intPt[2]])
                srepPts.InsertNextPoint([intPt[0], intPt[1], intPt[2]])
                # ax.scatter(intPt[0], intPt[1], intPt[2], color='k')
                if j == 0:
                    pass
                    # ax.plot([iXs[i], iSpoke[0][0]], [iYs[i], iSpoke[0][1]], [iZs[i], iSpoke[0][2]], 'r')
                else:
                    pass
                    # ax.plot([iSpoke[j-1][0], iSpoke[j][0]], [iSpoke[j-1][1], iSpoke[j][1]], [iSpoke[j-1][2], iSpoke[j][2]], 'r')
            iSamPts.append(iSpoke)
        else:
            iSpoke = []
            srepPts.InsertNextPoint([iXs[i], iYs[i], iZs[i]])
            for j in range(0, numSamp):
                intPt = tps.TransformPoint(rSamPts[2*i-1][j])
                iSpoke.append([intPt[0], intPt[1], intPt[2]])
                srepPts.InsertNextPoint([intPt[0], intPt[1], intPt[2]])
                # ax.scatter(intPt[0], intPt[1], intPt[2], color='k')
                if j == 0:
                    pass
                    # ax.plot([iXs[i], iSpoke[0][0]], [iYs[i], iSpoke[0][1]], [iZs[i], iSpoke[0][2]], 'r')
                else:
                    pass
                    # ax.plot([iSpoke[j-1][0], iSpoke[j][0]], [iSpoke[j-1][1], iSpoke[j][1]], [iSpoke[j-1][2], iSpoke[j][2]], 'r')
            iSamPts.append(iSpoke)

            iSpoke = []
            srepPts.InsertNextPoint([iXs[i], iYs[i], iZs[i]])
            for j in range(0, numSamp):
                intPt = tps.TransformPoint(rSamPts[2*i][j])
                iSpoke.append([intPt[0], intPt[1], intPt[2]])
                srepPts.InsertNextPoint([intPt[0], intPt[1], intPt[2]])
                # ax.scatter(intPt[0], intPt[1], intPt[2], color='k')
                if j == 0:
                    pass
                    # ax.plot([iXs[i], iSpoke[0][0]], [iYs[i], iSpoke[0][1]], [iZs[i], iSpoke[0][2]], 'r')
                else:
                    pass
                    # ax.plot([iSpoke[j-1][0], iSpoke[j][0]], [iSpoke[j-1][1], iSpoke[j][1]], [iSpoke[j-1][2], iSpoke[j][2]], 'r')
            iSamPts.append(iSpoke)
    # ax.plot(iXs, iYs, iZs, 'r')
    srepSet = vtk.vtkPolyData()
    srepSet.SetPoints(srepPts)
    srepSet.Modified()
    # writer2 = vtk.vtkPolyDataWriter()
    # writer2.SetInputData(srepSet)
    # writer2.SetFileName('pre2dsrep.vtk')
    # writer2.Write()
    





    #     boundPts.append([upx, upy])

    #     if i != 0 and i != len(sXs)-1:
    #         downx = sXs[i] + dirVecs[2*i][0]*opLen[2*i]
    #         downy = sYs[i] + dirVecs[2*i][1]*opLen[2*i]
    #         boundPts.append([downx, downy])
    #         xDown = [sXs[i], downx]
    #         yDown = [sYs[i], downy]
            # plt.plot(xDown, yDown,'k')



    # # vec1 = np.array([1, 1, 1]) 
    # # vec2 = np.array([-1, 1, -1])
    # # cross = np.cross(vec1, vec2)
    # # print(cross)
    # # ax.plot([0, 1], [0,1], [0, 1], 'b')
    # # ax.plot([0, -1], [0,1], [0, -1], 'b')
    # # ax.plot([0, cross[0]], [0, cross[1]], [0, cross[2]], 'b')

    # # ax.scatter(cross[0], cross[1], cross[2]+fit[2], color='b')

    # ax.plot([0, 10], [0,0], [0, 0], 'r')

    # plot plane
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    X,Y = np.meshgrid(np.arange(xlim[0], xlim[1]),
                    np.arange(ylim[0], ylim[1]))
    Z = np.zeros(X.shape)
    for r in range(X.shape[0]):
        for c in range(X.shape[1]):
            Z[r,c] = fit[0] * X[r,c] + fit[1] * Y[r,c] + fit[2]
    ax.plot_wireframe(X,Y,Z, color='b')

    # xlim = ax.get_xlim()
    # ylim = ax.get_ylim()
    # X,Y = np.meshgrid(np.arange(xlim[0], xlim[1]),
    #                   np.arange(ylim[0], ylim[1]))
    # Z = np.zeros(X.shape)
    # for r in range(X.shape[0]):
    #     for c in range(X.shape[1]):
    #         Z[r,c] = (0) * X[r,c] + (0) * Y[r,c]
    # ax.plot_wireframe(X,Y,Z, color='k')

    # ax.set_xlabel('x')
    # ax.set_ylabel('y')
    # ax.set_zlabel('z')
    # plt.ylim(-30, 30)
    # plt.xlim(-30, 30)
    # ax.set_zlim(-30,30)
    # plt.show()
    # return iXs, iYs, iZs, iSamPts, rXs, rYs, rZs, rSamPts
    return iXs, iYs, iZs, iSamPts
