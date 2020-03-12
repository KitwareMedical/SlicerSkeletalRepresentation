/*==============================================================================

  Program: 3D Slicer

  Portions (c) Copyright Brigham and Women's Hospital (BWH) All Rights Reserved.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

==============================================================================*/

// SkeletalRepresentationRefiner Logic includes
#include "vtkSlicerSkeletalRepresentationRefinerLogic.h"
#include <stdlib.h>
// MRML includes
#include <vtkMRMLScene.h>
#include <vtkMRMLModelNode.h>
#include <vtkMRMLModelDisplayNode.h>
#include <vtkMRMLDisplayNode.h>
#include <vtkMRMLMarkupsDisplayNode.h>
#include <vtkMRMLMarkupsFiducialNode.h>
#include <vtkMRMLMarkupsNode.h>
#include "vtkSlicerMarkupsLogic.h"
// VTK includes
#include <vtkIntArray.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>

#include <vtkPoints.h>
#include <vtkLine.h>
#include <vtkQuad.h>
#include <vtkTriangle.h>
#include <vtkPolyData.h>
#include <vtkFloatArray.h>
#include <vtkTubeFilter.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkExtractSurface.h>
#include <vtkImageData.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkPCANormalEstimation.h>
#include <vtkKochanekSpline.h>
#include <vtkParametricSpline.h>
#include <vtkCellLocator.h>
#include <vtksys/SystemTools.hxx>

#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkSignedDistance.h>
#include <vtkCleanPolyData.h>
#include <vtkPointSource.h>
#include <vtkCurvatures.h>
#include <vtkPointLocator.h>
#include <vtkAppendPolyData.h>
#include <vtkLookupTable.h>
#include <vtkParametricFunctionSource.h>
#include <vtkImplicitPolyDataDistance.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLDataParser.h>
#include <vtkSurfaceReconstructionFilter.h>
#include <vtkColorTransferFunction.h>
#include <vtkProgrammableSource.h>
#include <vtkContourFilter.h>
#include <vtkDistancePolyDataFilter.h>
#include <vtkReverseSense.h>
#include "vtkSlicerSkeletalRepresentationInterpolater.h"
#include "vtkSrep.h"
#include "vtkSpoke.h"
#include "newuoa.h"
#include "vtkPolyData2ImageData.h"
#include "vtkApproximateSignedDistanceMap.h"
#include "vtkGradientDistanceFilter.h"
#include <vtkBoundingBox.h>
#include <vtkMRMLProceduralColorNode.h>
// STD includes
#include <cassert>
const double voxelSpacing = 0.005;
const std::string newFilePrefix = "refined_";
//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkSlicerSkeletalRepresentationRefinerLogic);

//----------------------------------------------------------------------------
vtkSlicerSkeletalRepresentationRefinerLogic::vtkSlicerSkeletalRepresentationRefinerLogic()
{
}

//----------------------------------------------------------------------------
vtkSlicerSkeletalRepresentationRefinerLogic::~vtkSlicerSkeletalRepresentationRefinerLogic()
{
}

//----------------------------------------------------------------------------
void vtkSlicerSkeletalRepresentationRefinerLogic::PrintSelf(ostream& os, vtkIndent indent)
{
    this->Superclass::PrintSelf(os, indent);
}

void vtkSlicerSkeletalRepresentationRefinerLogic::SetImageFileName(const std::string &imageFilePath)
{
    mTargetMeshFilePath = imageFilePath;
    // visualize the input surface mesh
    vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
    reader->SetFileName(imageFilePath.c_str());
    reader->Update();

    vtkSmartPointer<vtkPolyData> surface = reader->GetOutput();
    Visualize(surface, "Input surface mesh", 0, 0, 0);
}

void vtkSlicerSkeletalRepresentationRefinerLogic::SetSrepFileName(const std::string &srepFilePath)
{
    mSrepFilePath = srepFilePath;
    int nRows = 0, nCols = 0;
    double crestShift = 0.0;
    std::string up, down, crest;
    ParseHeader(srepFilePath, &nRows, &nCols, &crestShift, &up, &down, &crest);
    if(nRows == 0 || nCols == 0)
    {
        std::cerr << "The s-rep model is empty." << std::endl;
        return;
    }
    std::vector<double> up_radii, down_radii, up_dirs, down_dirs, up_skeletalPoints, down_skeletalPoints;
    Parse(up, mCoeffArray, up_radii, up_dirs, up_skeletalPoints);

    vtkSrep *srep = new vtkSrep(nRows, nCols, up_radii, up_dirs, up_skeletalPoints);
    if(srep->IsEmpty())
    {
        std::cerr << "The s-rep model is empty." << std::endl;
        delete srep;
        srep = nullptr;
        return;
    }
    vtkSmartPointer<vtkPolyData> upSrepPoly = vtkSmartPointer<vtkPolyData>::New();
    ConvertSpokes2PolyData(srep->GetAllSpokes(), upSrepPoly);
    Visualize(upSrepPoly, "up spokes", 0, 1, 1);

    Parse(down, mCoeffArray, down_radii, down_dirs, down_skeletalPoints);

    vtkSrep *downSrep = new vtkSrep(nRows, nCols, down_radii, down_dirs, down_skeletalPoints);
    if(downSrep->IsEmpty())
    {
        std::cerr << "The s-rep model is empty." << std::endl;
        delete downSrep;
        downSrep = nullptr;
        return;
    }
    vtkSmartPointer<vtkPolyData> downSrepPoly = vtkSmartPointer<vtkPolyData>::New();
    ConvertSpokes2PolyData(downSrep->GetAllSpokes(), downSrepPoly);
    Visualize(downSrepPoly, "down spokes", 1, 0, 1);

    std::vector<vtkSpoke*> crestSpokes, reorderedCrest;
    ParseCrest(crest, crestSpokes);

    vtkSmartPointer<vtkPolyData> crestSrepPoly = vtkSmartPointer<vtkPolyData>::New();
    ConvertSpokes2PolyData(crestSpokes, crestSrepPoly);
    Visualize(crestSrepPoly, "crest spokes", 1, 0, 0);

    vtkSmartPointer<vtkPolyData> foldPoly = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPoints> foldCurvePts = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> foldCurveCell = vtkSmartPointer<vtkCellArray>::New();

    for (int i = 0; i < crestSpokes.size() - 1; ++i) {
        double basePoint[3], bdryPoint[3];
        crestSpokes[i]->GetSkeletalPoint(basePoint);
        crestSpokes[i+1]->GetSkeletalPoint(bdryPoint);
        int id0 = foldCurvePts->InsertNextPoint(basePoint);
        int id1 = foldCurvePts->InsertNextPoint(bdryPoint);
        vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
        line->GetPointIds()->SetId(0, id0);
        line->GetPointIds()->SetId(1, id1);
        foldCurveCell->InsertNextCell(line);

        if(i == crestSpokes.size() - 2) {
            vtkSmartPointer<vtkLine> line_end = vtkSmartPointer<vtkLine>::New();
            line_end->GetPointIds()->SetId(0, id1);
            line_end->GetPointIds()->SetId(1, 0);
            foldCurveCell->InsertNextCell(line_end);
        }
    }
    foldPoly->SetPoints(foldCurvePts);
    foldPoly->SetPolys(foldCurveCell);
    foldPoly->Modified();
    Visualize(foldPoly, "Fold curve", 1, 1, 0);
    // show fold curve
//    ReorderCrestSpokes(nRows, nCols, crestSpokes, reorderedCrest);

}

void vtkSlicerSkeletalRepresentationRefinerLogic::SetOutputPath(const string &outputPath)
{
    mOutputPath = outputPath;
}

void vtkSlicerSkeletalRepresentationRefinerLogic::Refine(double stepSize, double endCriterion, int maxIter, int interpolationLevel)
{
    mFirstCost = true;
    // 1. parse file
    const std::string headerFileName = mSrepFilePath;
    int nRows = 0, nCols = 0;
    double crestShift = 0.0;
    std::string up, down, crest;
    ParseHeader(headerFileName, &nRows, &nCols, &crestShift, &up, &down, &crest);

    if(nRows == 0 || nCols == 0)
    {
        std::cerr << "The s-rep model is empty." << std::endl;
        return;
    }

    mNumCols = nCols;
    mNumRows = nRows;

    // Prepare signed distance image
    AntiAliasSignedDistanceMap(mTargetMeshFilePath);

    // Compute transformation matrix from srep to image coordinate system, namely, unit cube cs.
    TransformSrep(headerFileName);

    // make tuples of interpolation positions (u,v)
    mInterpolationLevel = interpolationLevel;
    mInterpolatePositions.clear();
    double tol = 1e-6;
    int shares = static_cast<int>(pow(2, interpolationLevel));
    double interval = double(1.0 / shares);
    for(int i = 0; i <= shares; ++i)
    {
        for(int j = 0; j <= shares; ++j)
        {
            double u = i * interval;
            double v = j * interval;
            // no interpolation at corners
            if((abs(u) < tol && abs(v) < tol) || (abs(u) < tol && abs(v-1) < tol)
                    || (abs(u-1) < tol && abs(v) < tol) || (abs(u-1) < tol && abs(v-1) < tol))
                continue;
            std::pair<double, double> uv = make_pair(u, v);
            mInterpolatePositions.push_back(uv);
        }
    }
    if(interpolationLevel == 0)
    {
        mInterpolatePositions.push_back(std::pair<double, double>(0, 0));
    }
    // Hide other nodes.
    HideNodesByClass("vtkMRMLModelNode");

    // Refine up spokes
    mUpSpokes = RefinePartOfSpokes(up, stepSize, endCriterion, maxIter);

    // Refine down spokes
    mDownSpokes = RefinePartOfSpokes(down, stepSize, endCriterion, maxIter);

    RefineCrestSpokes(crest, stepSize, endCriterion, maxIter);

    // Update header file
    std::string newHeaderFileName;
    UpdateHeader(headerFileName, mOutputPath, &newHeaderFileName);

    ShowImpliedBoundary(interpolationLevel, newHeaderFileName);
}
void vtkSlicerSkeletalRepresentationRefinerLogic::InterpolateSrep(int interpolationLevel,
                                                                  std::string& srepFileName)
{
    // Hide other nodes.
    HideNodesByClass("vtkMRMLModelNode");
    std::vector<vtkSpoke*> temp;

    // 1. Parse the model into a parameter array that needs to be optimized
    int nRows = 0, nCols = 0;
    std::string up, down, crest;
    double crestShift = 0.0;
    ParseHeader(srepFileName, &nRows, &nCols, &crestShift, &up, &down, &crest);
    if(nRows == 0 || nCols == 0)
    {
        std::cerr << "The s-rep model is empty." << std::endl;
        return;
    }
    InterpolateSrep(interpolationLevel, nRows, nCols, down, crest, temp);
    std::vector<vtkSpoke*> interpolatedSpokes, upSpokes, downSpokes, crestSpokes;
    vtkSmartPointer<vtkPolyData> wireFrame = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPoints> pts, foldCurvePts;
    vtkSmartPointer<vtkCellArray> quads, foldCurveCell;
    pts = vtkSmartPointer<vtkPoints>::New();
    foldCurvePts = vtkSmartPointer<vtkPoints>::New();
    quads = vtkSmartPointer<vtkCellArray>::New();
    foldCurveCell = vtkSmartPointer<vtkCellArray>::New();

    vtkSmartPointer<vtkPolyData> foldCurve = vtkSmartPointer<vtkPolyData>::New();

    // connect implied boundary for up spokes
    vtkSmartPointer<vtkPolyData> northPolyData = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPolyData> southPolyData = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPolyData> crestPolyData = vtkSmartPointer<vtkPolyData>::New();
    std::vector<vtkSpoke*> borderUpSpoke, borderDownSpoke, repeatedInterp;
    ConnectImpliedBoundaryPts(interpolationLevel, nRows, nCols, up,
                              borderUpSpoke, interpolatedSpokes,repeatedInterp, upSpokes,northPolyData);

}

void vtkSlicerSkeletalRepresentationRefinerLogic::InterpolateSrep(int interpolationLevel, int nRows, int nCols,
                                                                  std::string& up, std::string& crest,
                                                                  std::vector<vtkSpoke*> &interpolatedSpokes)
{
    std::vector<vtkSpoke *> interpolatedCrestSpokes;
    InterpolateSrep(interpolationLevel, nRows, nCols, up, crest, interpolatedSpokes, interpolatedCrestSpokes);
}

void vtkSlicerSkeletalRepresentationRefinerLogic::InterpolateSrep(int interpolationLevel, int nRows, int nCols,
                                                                  string &srepFileName,
                                                                  string &crestFileName,
                                                                  std::vector<vtkSpoke *> &interpolatedSpokes,
                                                                  std::vector<vtkSpoke *> &interpolatedCrestSpokes)
{
    std::vector<double> coeffArrayUp, radiiUp, dirsUp, skeletalPointsUp;
    Parse(srepFileName, coeffArrayUp, radiiUp, dirsUp, skeletalPointsUp);
    vtkSrep *srep = new vtkSrep(nRows, nCols, radiiUp, dirsUp, skeletalPointsUp);
    if(srep->IsEmpty())
    {
        std::cerr << "The s-rep model is empty." << std::endl;
        delete srep;
        srep = nullptr;
        return;
    }
    // 1.1 interpolate and visualize for verification
    // collect neighboring spokes around corners
    vtkSlicerSkeletalRepresentationInterpolater interpolater;

    int shares = static_cast<int>(pow(2, interpolationLevel));
    double interval = static_cast<double>(1.0/ shares);
    std::vector<double> steps;
    std::vector<vtkSpoke *> repeatedInterpInterior;
    for(int i = 0; i <= shares; ++i)
    {
        steps.push_back(i * interval);
    }

    int nCrestPoints = nRows*2 + (nCols-2)*2;
    int numCols = 1+static_cast<int>(floor(nRows/2)); // steps from crest point to the skeletal point

    int rowOffset = 0;//static_cast<int>(floor(nRows/2));
    for(int r = 0; r < nCrestPoints; ++r)
    {
        for(int c = 0; c < numCols-1; ++c)
        {
            vtkSpoke *cornerSpokes[4];

            double  dXdu11[3], dXdv11[3],
                    dXdu12[3], dXdv12[3],
                    dXdu21[3], dXdv21[3],
                    dXdu22[3], dXdv22[3];

            for(size_t i = 0; i < steps.size(); ++i)
            {
                for(size_t j = 0; j < steps.size(); ++j)
                {
                    cornerSpokes[0] = srep->GetSpoke(r,c);
                    cornerSpokes[1] = srep->GetSpoke(r+1, c);
                    cornerSpokes[2] = srep->GetSpoke(r+1, c+1);
                    cornerSpokes[3] = srep->GetSpoke(r, c+ 1);

                    ComputeDerivative(skeletalPointsUp, r, c, nRows, nCols, dXdu11, dXdv11);
                    ComputeDerivative(skeletalPointsUp, r+1, c, nRows, nCols, dXdu21, dXdv21);
                    ComputeDerivative(skeletalPointsUp, r, c+1, nRows, nCols, dXdu12, dXdv12);
                    ComputeDerivative(skeletalPointsUp, r+1, c+1, nRows, nCols, dXdu22, dXdv22);
                    if(r == nCrestPoints - 1) {
                        ComputeDerivative(skeletalPointsUp, r, c, nRows, nCols, dXdu11, dXdv11);
                        ComputeDerivative(skeletalPointsUp, 0, c, nRows, nCols, dXdu21, dXdv21);
                        ComputeDerivative(skeletalPointsUp, r, c+1, nRows, nCols, dXdu12, dXdv12);
                        ComputeDerivative(skeletalPointsUp, 0, c+1, nRows, nCols, dXdu22, dXdv22);
                    }
                    interpolater.SetCornerDxdu(dXdu11,
                                               dXdu21,
                                               dXdu22,
                                               dXdu12);
                    interpolater.SetCornerDxdv(dXdv11,
                                               dXdv21,
                                               dXdv22,
                                               dXdv12);

                    vtkSpoke* in1 = new vtkSpoke;
                    interpolater.Interpolate(double(steps[i]), double(steps[j]), cornerSpokes, in1);
                    repeatedInterpInterior.push_back(in1);
                    if((c != numCols - 2 && j == steps.size() -1) || // edges from inside to outside
                            (c == 0 && r >= nCrestPoints/2 && r < nCrestPoints-rowOffset && j == 0)) {//bottom half at the innermost level
                        // avoid redundent interpolation on edges

                    }
                    else {
                        interpolatedSpokes.push_back(in1);
                    }

                }
            }
        }
    }

    vtkSmartPointer<vtkPolyData> upSpokes_polyData = vtkSmartPointer<vtkPolyData>::New();
    ConvertSpokes2PolyData(interpolatedSpokes, upSpokes_polyData);
    Visualize(upSpokes_polyData, "Interpolated", 1, 1, 1);

    vtkSmartPointer<vtkPolyData> primarySpokes = vtkSmartPointer<vtkPolyData>::New();
    ConvertSpokes2PolyData(srep->GetAllSpokes(), primarySpokes);
    Visualize(primarySpokes, "Primary", 1, 0, 0);

    std::vector<vtkSpoke*> crestSpokes, topCrest;
    ParseCrest(crestFileName, crestSpokes);

    std::vector<vtkSpoke *> tempSpokes;
    InterpolateCrest(crestSpokes, repeatedInterpInterior, interpolationLevel, nRows, nCols, interpolatedCrestSpokes, tempSpokes);

    vtkSmartPointer<vtkPolyData> crestSpokes_poly = vtkSmartPointer<vtkPolyData>::New();
    ConvertSpokes2PolyData(interpolatedCrestSpokes, crestSpokes_poly);
    Visualize(crestSpokes_poly, "Crest", 0, 0, 1);

    vtkSmartPointer<vtkPolyData> crestSpokes_primary = vtkSmartPointer<vtkPolyData>::New();
    ConvertSpokes2PolyData(crestSpokes, crestSpokes_primary);
    Visualize(crestSpokes_primary, "Crest Primary", 0, 1, 1);

    // delete pointers
    delete srep;
}

void vtkSlicerSkeletalRepresentationRefinerLogic::SetWeights(double wtImageMatch, double wtNormal, double wtSrad)
{
    mWtImageMatch = wtImageMatch;
    mWtNormalMatch = wtNormal;
    mWtSrad = wtSrad;
}

double vtkSlicerSkeletalRepresentationRefinerLogic::operator ()(double *coeff)
{
    double cost = 0.0;
    cost = EvaluateObjectiveFunction(coeff);
    return cost;
}
double vtkSlicerSkeletalRepresentationRefinerLogic::EvaluateObjectiveFunction(double *coeff)
{
    // TODO: SHOULD add progress bar here.
    //std::cout << "Current iteration: " << iterNum++ << std::endl;

    if(mSrep == nullptr)
    {
        std::cerr << "The srep pointer in the refinement is nullptr." << std::endl;
        return -100000.0;
    }

    // this temporary srep is constructed to compute the cost function value
    // The original srep should not be changed by each iteration
    vtkSrep *tempSrep = new vtkSrep();
    tempSrep->DeepCopy(*mSrep);
    tempSrep->Refine(coeff);
    double imageDist = 0.0, normal = 0.0, srad = 0.0;
    int paramDim = static_cast<int>(mCoeffArray.size());
    int spokeNum = paramDim / 4;

    int nCrestPoints = mNumRows*2 + (mNumCols-2)*2;
    int num_cols = 1+static_cast<int>(floor(mNumRows/2)); // steps from crest point to the skeletal point

    // 1. Compute image match from all spokes and those spokes affected by them
    for(int i = 0; i < spokeNum; ++i)
    {
        int r = i / num_cols;
        int c = i % num_cols;
        vtkSpoke *thisSpoke = tempSrep->GetSpoke(r, c);

        // compute distance for this spoke
        imageDist += ComputeDistance(thisSpoke, &normal);

        for(auto it = mInterpolatePositions.begin(); it != mInterpolatePositions.end(); ++it)
        {
            double u = (*it).first;
            double v = (*it).second;

            // For each spoke at the corner of the srep,
            // its neighbors are all spokes in one quad
            if(r == 0 && c == 0)
            {
                // left-top corner
                imageDist += TotalDistOfLeftTopSpoke(tempSrep, u, v, r, c, &normal);
            }
            else if(r == 0 && c == num_cols - 1)
            {
                // right-top corner
                imageDist += TotalDistOfRightTopSpoke(tempSrep, u, v, r,c, &normal);
            }
            else if(r == nCrestPoints - 1 && c == 0)
            {
                // left-bot corner
                imageDist += TotalDistOfLeftBotSpoke(tempSrep, u, v, r,c, &normal);
            }
            else if(r == nCrestPoints - 1 && c == num_cols - 1)
            {
                // right-bot corner
                imageDist += TotalDistOfRightBotSpoke(tempSrep, u, v, r, c, &normal);
            }
            // For each spoke on the edge of the srep,
            // its neighbors are all spokes in two quads
            else if(r == 0)
            {
                // top edge in middle
                imageDist += TotalDistOfRightTopSpoke(tempSrep, u, v, r, c, &normal);
                imageDist += TotalDistOfLeftTopSpoke(tempSrep, u, v, r, c, &normal);
            }
            else if(r == nCrestPoints - 1)
            {
                // bot edge in middle
                imageDist += TotalDistOfRightBotSpoke(tempSrep, u, v, r, c, &normal);
                imageDist += TotalDistOfLeftBotSpoke(tempSrep, u, v, r, c, &normal);
            }
            else if(c == 0)
            {
                // left edge in middle
                imageDist += TotalDistOfLeftBotSpoke(tempSrep, u, v, r, c, &normal);
                imageDist += TotalDistOfLeftTopSpoke(tempSrep, u, v, r, c, &normal);
            }
            else if(c == num_cols - 1)
            {
                // right edge in middle
                imageDist += TotalDistOfRightBotSpoke(tempSrep, u, v, r, c, &normal);
                imageDist += TotalDistOfRightTopSpoke(tempSrep, u, v, r, c, &normal);
            }
            // for each spoke in the middle of the srep,
            // obtain image distance and normal from all interpolated spoke in 4 quads around it
            else
            {
                imageDist += TotalDistOfRightBotSpoke(tempSrep, u, v, r, c, &normal);
                imageDist += TotalDistOfRightTopSpoke(tempSrep, u, v, r, c, &normal);

                imageDist += TotalDistOfLeftBotSpoke(tempSrep, u, v, r, c, &normal);
                imageDist += TotalDistOfLeftTopSpoke(tempSrep, u, v, r, c, &normal);
            }
        }

    }

    // 2. compute srad penalty
    srad = ComputeRSradPenalty(tempSrep);

    if(mFirstCost)
    {
        // this log helps to adjust the weights of three terms
        std::cout << "ImageMatch:" << imageDist << ", normal:" << normal << ", srad:" << srad << std::endl;
        mFirstCost = false;
    }

    delete tempSrep;
    return mWtImageMatch * imageDist + mWtNormalMatch * normal + mWtSrad * srad;
}

void vtkSlicerSkeletalRepresentationRefinerLogic::AntiAliasSignedDistanceMap(const std::string &meshFileName)
{
    // 1. convert poly data to image data
    vtkPolyData2ImageData polyDataConverter;
    vtkSmartPointer<vtkImageData> img = vtkSmartPointer<vtkImageData>::New();

    // this conversion already put the image into the unit-cube
    polyDataConverter.Convert(meshFileName, img);

    vtkSmartPointer<vtkImageData> antiAliasedImage = vtkSmartPointer<vtkImageData>::New();

    vtkApproximateSignedDistanceMap ssdGenerator;
    ssdGenerator.Convert(img, mAntiAliasedImage);

    // 4. compute normals of the image everywhere
    vtkGradientDistanceFilter gradDistFilter;
    gradDistFilter.Filter(mAntiAliasedImage, mGradDistImage);
}

void vtkSlicerSkeletalRepresentationRefinerLogic::TransformSrep(const std::string &headerFile)
{

    int nRows = 0, nCols = 0;
    std::string up, down, crest;
    double crestShift = 0.0;
    ParseHeader(headerFile, &nRows, &nCols, &crestShift, &up, &down, &crest);

    if(nRows == 0 || nCols == 0)
    {
        std::cerr << "The s-rep model is empty." << std::endl;
        return;
    }

    std::vector<double> radiiUp, dirsUp, skeletalPointsUp, coeffUp;
    Parse(up, coeffUp, radiiUp, dirsUp, skeletalPointsUp);

    vtkSrep *srep = new vtkSrep(nRows, nCols, radiiUp, dirsUp, skeletalPointsUp);
    if(srep->IsEmpty())
    {
        std::cerr << "The s-rep model is empty." << std::endl;
        delete srep;
        srep = nullptr;
        return;
    }

    std::vector<double> radiiDown, dirsDown, skeletalPointsDown, coeffDown;
    Parse(down, coeffDown, radiiDown, dirsDown, skeletalPointsDown);
    srep->AddSpokes(radiiDown, dirsDown, skeletalPointsDown);

    std::vector<double> radiiCrest, dirsCrest, skeletalPointsCrest, coeffCrest;
    Parse(crest, coeffCrest, radiiCrest, dirsCrest, skeletalPointsCrest);
    srep->AddSpokes(radiiCrest, dirsCrest, skeletalPointsCrest);

    TransformSrep2ImageCS(srep, mTransformationMat);

    delete srep;
}

void vtkSlicerSkeletalRepresentationRefinerLogic::ShowImpliedBoundary(int interpolationLevel,
                                                                      const string &srepFileName)
{
    // Hide other nodes.
    HideNodesByClass("vtkMRMLModelNode");

    vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
    reader->SetFileName(mTargetMeshFilePath.c_str());
    reader->Update();
    vtkSmartPointer<vtkPolyData> inputMesh = reader->GetOutput();
    // 1. Parse the model into a parameter array that needs to be optimized
    int nRows = 0, nCols = 0;
    std::string up, down, crest;
    double crestShift = 0.0;
    ParseHeader(srepFileName, &nRows, &nCols, &crestShift, &up, &down, &crest);

    std::vector<vtkSpoke*> interpolatedSpokes, upSpokes, downSpokes, crestSpokes;
    vtkSmartPointer<vtkPolyData> wireFrame = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPoints> pts, foldCurvePts;
    vtkSmartPointer<vtkCellArray> quads, foldCurveCell;
    pts = vtkSmartPointer<vtkPoints>::New();
    foldCurvePts = vtkSmartPointer<vtkPoints>::New();
    quads = vtkSmartPointer<vtkCellArray>::New();
    foldCurveCell = vtkSmartPointer<vtkCellArray>::New();

    vtkSmartPointer<vtkPolyData> foldCurve = vtkSmartPointer<vtkPolyData>::New();

    // connect implied boundary for up spokes
    vtkSmartPointer<vtkPolyData> northPolyData = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPolyData> southPolyData = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPolyData> crestPolyData = vtkSmartPointer<vtkPolyData>::New();
    std::vector<vtkSpoke*> borderUpSpoke, borderDownSpoke, interpolatedCrestSpokes, repeatedInterpUp, repeatedInterpDown;
    ConnectImpliedBoundaryPts(interpolationLevel, nRows, nCols, up,
                              borderUpSpoke, interpolatedSpokes, repeatedInterpUp,upSpokes,northPolyData);

    //cout << "north polydata num points:" << northPolyData->GetNumberOfPoints() << endl;
    ConnectImpliedBoundaryPts(interpolationLevel, nRows, nCols, down,
                              borderDownSpoke,
                              interpolatedSpokes, repeatedInterpDown, downSpokes, southPolyData);

    std::vector<vtkSpoke *> tempSpokes;
    ParseCrest(crest, crestSpokes);
    InterpolateCrest(crestSpokes, repeatedInterpUp, interpolationLevel, nRows, nCols, interpolatedCrestSpokes, tempSpokes);

    vtkSmartPointer<vtkAppendPolyData> appendFilter =
      vtkSmartPointer<vtkAppendPolyData>::New();
    appendFilter->AddInputData(northPolyData);
    appendFilter->AddInputData(southPolyData);

    ConnectImpliedCrest(interpolationLevel,
                        borderUpSpoke,
                        borderDownSpoke,
                        interpolatedCrestSpokes,
                        crestPolyData);
    appendFilter->AddInputData(crestPolyData);
    vtkSmartPointer<vtkPolyData> impliedBoundary = vtkSmartPointer<vtkPolyData>::New();
    appendFilter->Update();
    vtkSmartPointer<vtkCleanPolyData> cleanFilter =
        vtkSmartPointer<vtkCleanPolyData>::New();
    cleanFilter->SetInputConnection(appendFilter->GetOutputPort());
    cleanFilter->Update();
    vtkSmartPointer<vtkPolyData> retImpliedBoundary = appendFilter->GetOutput();

    vtkSmartPointer<vtkPolyData> retiledMesh = vtkSmartPointer<vtkPolyData>::New();
    RetileMesh(inputMesh, retImpliedBoundary, retiledMesh);
    vtkSmartPointer<vtkSmoothPolyDataFilter> smoothFilter =
            vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
    smoothFilter->SetInputData(retiledMesh);
    smoothFilter->SetNumberOfIterations(10);
    smoothFilter->SetRelaxationFactor(0.05); // 0.01
    smoothFilter->FeatureEdgeSmoothingOff();
    smoothFilter->BoundarySmoothingOn();
    smoothFilter->Update();
    impliedBoundary->DeepCopy(smoothFilter->GetOutput());
    Visualize(impliedBoundary, "Implied boundary", 1, 1, 0);

}

void vtkSlicerSkeletalRepresentationRefinerLogic::CLIRefine(const std::string &srepFileName,
                                                            const std::string &imgFileName,
                                                            const std::string &outputPath,
                                                            double stepSize, double endCriterion,
                                                            int maxIter, double wtImg,
                                                            double wtNormal, double wtSrad,
                                                            int interpolationLevel)
{
    SetSrepFileName(srepFileName);
    SetImageFileName(imgFileName);
    SetOutputPath(outputPath);
    SetWeights(wtImg, wtNormal, wtSrad);
    Refine(stepSize, endCriterion, maxIter, interpolationLevel);
}

void vtkSlicerSkeletalRepresentationRefinerLogic::ComputeDerivative(std::vector<double> skeletalPoints,
                                                                    int intr, int intc, int nRows,
                                                                    int intCols, double *dXdu, double *dXdv)
{
    int nCrestPoints = nRows*2 + (intCols-2)*2;
    int numCols = 1+static_cast<int>(floor(nRows/2)); // steps from crest point to the skeletal point

    // 0-based index of elements if arranged in array
    size_t nCols = static_cast<size_t>(numCols);
    size_t r = static_cast<size_t>(intr);
    size_t c= static_cast<size_t>(intc);
    size_t id = static_cast<size_t>(r * nCols + c);
    double head[3], tail[3];
    double factor = 0.5;
    if(r == 0)
    {
        // first row
        // forward difference, next row/col - current row/col
        head[0] = skeletalPoints[(id+nCols)*3];
        head[1] = skeletalPoints[(id+nCols)*3+1];
        head[2] = skeletalPoints[(id+nCols)*3+2];

        tail[0] = skeletalPoints[(id)*3];
        tail[1] = skeletalPoints[(id)*3+1];
        tail[2] = skeletalPoints[(id)*3+2];
        factor = 1.0;
    }
    else if(r == static_cast<size_t>(nCrestPoints - 1))
    {
        // last row
        // backward difference
        tail[0] = skeletalPoints[(id-nCols)*3];
        tail[1] = skeletalPoints[(id-nCols)*3+1];
        tail[2] = skeletalPoints[(id-nCols)*3+2];

        head[0] = skeletalPoints[(id)*3];
        head[1] = skeletalPoints[(id)*3+1];
        head[2] = skeletalPoints[(id)*3+2];

        head[0] = skeletalPoints[(c)*3];
        head[1] = skeletalPoints[(c)*3+1];
        head[2] = skeletalPoints[(c)*3+2];
        factor = 1.0;
    }
    else
    {
        // otherwise, center difference
        head[0] = skeletalPoints[(id+nCols)*3];
        head[1] = skeletalPoints[(id+nCols)*3+1];
        head[2] = skeletalPoints[(id+nCols)*3+2];

        tail[0] = skeletalPoints[(id-nCols)*3];
        tail[1] = skeletalPoints[(id-nCols)*3+1];
        tail[2] = skeletalPoints[(id-nCols)*3+2];
        factor = 0.5;
    }
    ComputeDiff(head, tail, factor, dXdu);

    if(c == 0)
    {
        // first col
        head[0] = skeletalPoints[(id+1)*3];
        head[1] = skeletalPoints[(id+1)*3+1];
        head[2] = skeletalPoints[(id+1)*3+2];

        tail[0] = skeletalPoints[(id)*3];
        tail[1] = skeletalPoints[(id)*3+1];
        tail[2] = skeletalPoints[(id)*3+2];
        factor = 1.0;
    }
    else if(c == nCols - 1)
    {
        // last col
        // backward difference
        tail[0] = skeletalPoints[(id-1)*3];
        tail[1] = skeletalPoints[(id-1)*3+1];
        tail[2] = skeletalPoints[(id-1)*3+2];

        head[0] = skeletalPoints[(id)*3];
        head[1] = skeletalPoints[(id)*3+1];
        head[2] = skeletalPoints[(id)*3+2];
        factor = 1.0;
    }
    else
    {
        // otherwise, center difference
        head[0] = skeletalPoints[(id+1)*3];
        head[1] = skeletalPoints[(id+1)*3+1];
        head[2] = skeletalPoints[(id+1)*3+2];

        tail[0] = skeletalPoints[(id-1)*3];
        tail[1] = skeletalPoints[(id-1)*3+1];
        tail[2] = skeletalPoints[(id-1)*3+2];
        factor = 0.5;
    }
    ComputeDiff(head, tail, factor, dXdv);

}

void vtkSlicerSkeletalRepresentationRefinerLogic::ConvertSpokes2PolyData(std::vector<vtkSpoke *>input, vtkPolyData *output)
{
    vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> arrows = vtkSmartPointer<vtkCellArray>::New();

    size_t numSpokes = input.size();
    for(size_t i = 0; i < numSpokes; ++i)
    {
        vtkSpoke* currSpoke = input[i];
        double basePt[3], bdryPt[3], dir[3];
        currSpoke->GetSkeletalPoint(basePt);
        currSpoke->GetBoundaryPoint(bdryPt);
        currSpoke->GetDirection(dir);
        vtkIdType id0 = (pts->InsertNextPoint(basePt[0], basePt[1], basePt[2]));
        vtkIdType id1 = pts->InsertNextPoint(bdryPt[0], bdryPt[1], bdryPt[2]);

        vtkSmartPointer<vtkLine> currLine = vtkSmartPointer<vtkLine>::New();
        currLine->GetPointIds()->SetId(0, id0);
        currLine->GetPointIds()->SetId(1, id1);
        arrows->InsertNextCell(currLine);
    }
    output->SetPoints(pts);
    output->SetLines(arrows);
    output->Modified();
}

void vtkSlicerSkeletalRepresentationRefinerLogic::SaveSpokes2Vtp(std::vector<vtkSpoke *> input, const string &path)
{
    vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkDoubleArray> spokeDirection = vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> spokeLengths = vtkSmartPointer<vtkDoubleArray>::New();

    size_t numSpokes = input.size();
    spokeLengths->SetNumberOfComponents(1);
    spokeLengths->SetName("spokeLength");

    spokeDirection->SetNumberOfComponents(3);
    spokeDirection->SetName("spokeDirection");

    for(size_t i = 0; i < numSpokes; ++i)
    {
        vtkSpoke* currSpoke = input[i];
        double basePt[3], dir[3], radius;
        currSpoke->GetSkeletalPoint(basePt);
        radius = currSpoke->GetRadius();
        currSpoke->GetDirection(dir);
        pts->InsertNextPoint(basePt[0], basePt[1], basePt[2]);
        spokeDirection->InsertNextTuple(dir);
        spokeLengths->InsertNextTuple1(radius);
    }
    vtkSmartPointer<vtkPolyData> output = vtkSmartPointer<vtkPolyData>::New();
    output->SetPoints(pts);

    output->GetPointData()->AddArray(spokeDirection);
    output->GetPointData()->SetActiveVectors("spokeDirection");
    output->GetPointData()->AddArray(spokeLengths);
    output->GetPointData()->SetActiveScalars("spokeLength");

    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(path.c_str());
    writer->SetInputData(output);
    writer->Update();
}

void vtkSlicerSkeletalRepresentationRefinerLogic::TransSpokes2PolyData(std::vector<vtkSpoke *>input, vtkPolyData *output)
{
    vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> arrows = vtkSmartPointer<vtkCellArray>::New();

    for(size_t i = 0; i < input.size(); ++i)
    {
        vtkSpoke* currSpoke = input[i];
        double basePt[3], bdryPt[3];
        currSpoke->GetSkeletalPoint(basePt);
        currSpoke->GetBoundaryPoint(bdryPt);

        basePt[0] = basePt[0] * mTransformationMat[0][0] + mTransformationMat[3][0];
        basePt[1] = basePt[1] * mTransformationMat[1][1] + mTransformationMat[3][1];
        basePt[2] = basePt[2] * mTransformationMat[2][2] + mTransformationMat[3][2];

        bdryPt[0] = bdryPt[0] * mTransformationMat[0][0] + mTransformationMat[3][0];
        bdryPt[1] = bdryPt[1] * mTransformationMat[1][1] + mTransformationMat[3][1];
        bdryPt[2] = bdryPt[2] * mTransformationMat[2][2] + mTransformationMat[3][2];

        vtkIdType id0 = pts->InsertNextPoint(basePt[0], basePt[1], basePt[2]);
        vtkIdType id1 = pts->InsertNextPoint(bdryPt[0], bdryPt[1], bdryPt[2]);

        vtkSmartPointer<vtkLine> currLine = vtkSmartPointer<vtkLine>::New();
        currLine->GetPointIds()->SetId(0, id0);
        currLine->GetPointIds()->SetId(1, id1);
        arrows->InsertNextCell(currLine);
    }
    output->SetPoints(pts);
    output->SetLines(arrows);
}

void vtkSlicerSkeletalRepresentationRefinerLogic::VisualizePoints(vtkPoints *input)
{
    vtkMRMLScene *scene = this->GetMRMLScene();
    if(!scene)
    {
        vtkErrorMacro(" Invalid scene");
        return;
    }

    vtkSmartPointer<vtkMRMLMarkupsFiducialNode> fidNode = vtkSmartPointer<vtkMRMLMarkupsFiducialNode>::New();

    vtkSmartPointer<vtkMRMLMarkupsDisplayNode> fidDisplayNode = vtkSmartPointer<vtkMRMLMarkupsDisplayNode>::New();
    scene->AddNode(fidDisplayNode);
    fidDisplayNode->SetGlyphScale(0.01);
    fidDisplayNode->SetSelectedColor(1.0, 1.0, 0.0);
    fidDisplayNode->SetTextScale(0.0);
    fidDisplayNode->SetName("surface points");
    scene->AddNode(fidNode);
    fidNode->SetAndObserveDisplayNodeID(fidDisplayNode->GetID());


    fidNode->SetLocked(true);
    for (int i = 0; i < input->GetNumberOfPoints(); ++i) {
        double pt[3];
        input->GetPoint(i, pt);
        fidNode->AddFiducial(pt[0], pt[1], pt[2]);
    }

}

//---------------------------------------------------------------------------
void vtkSlicerSkeletalRepresentationRefinerLogic::SetMRMLSceneInternal(vtkMRMLScene * newScene)
{
  vtkNew<vtkIntArray> events;
  events->InsertNextValue(vtkMRMLScene::NodeAddedEvent);
  events->InsertNextValue(vtkMRMLScene::NodeRemovedEvent);
  events->InsertNextValue(vtkMRMLScene::EndBatchProcessEvent);
  this->SetAndObserveMRMLSceneEventsInternal(newScene, events.GetPointer());
}

//-----------------------------------------------------------------------------
void vtkSlicerSkeletalRepresentationRefinerLogic::RegisterNodes()
{
  assert(this->GetMRMLScene() != nullptr);
}

//---------------------------------------------------------------------------
void vtkSlicerSkeletalRepresentationRefinerLogic::UpdateFromMRMLScene()
{
  assert(this->GetMRMLScene() != nullptr);
}

//---------------------------------------------------------------------------
void vtkSlicerSkeletalRepresentationRefinerLogic
::OnMRMLSceneNodeAdded(vtkMRMLNode* vtkNotUsed(node))
{
}

//---------------------------------------------------------------------------
void vtkSlicerSkeletalRepresentationRefinerLogic
::OnMRMLSceneNodeRemoved(vtkMRMLNode* vtkNotUsed(node))
{
}
void vtkSlicerSkeletalRepresentationRefinerLogic::Parse(const std::string &modelFileName, std::vector<double> &coeffArray,
                                                        std::vector<double> &radii, std::vector<double> &dirs, std::vector<double> &skeletalPoints)
{

    vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
    reader->SetFileName(modelFileName.c_str());
    reader->Update();

    vtkSmartPointer<vtkPolyData> spokesPolyData = reader->GetOutput();
    Visualize(spokesPolyData, "Skeletal sheet", 0, 0, 0);
    vtkSmartPointer<vtkPointData> spokesPointData = spokesPolyData->GetPointData();
    int numOfArrays = spokesPointData->GetNumberOfArrays();
    vtkIdType numOfSpokes = spokesPolyData->GetNumberOfPoints();

    if(numOfSpokes == 0 || numOfArrays == 0)
    {
        return;
    }

    // including Ux, Uy, Uz, r
    vtkSmartPointer<vtkDoubleArray> spokeRadii = vtkDoubleArray::SafeDownCast(spokesPointData->GetArray("spokeLength"));
    vtkSmartPointer<vtkDoubleArray> spokeDirs = vtkDoubleArray::SafeDownCast(spokesPointData->GetArray("spokeDirection"));

    for(int i = 0; i < numOfSpokes; ++i)
    {
        int idxDir = i * 3; // Ux, Uy, Uz

        // coefficients (dirs + radii) for newuoa
        // the coefficient for radii is the exponential value, initially 0
        coeffArray.push_back(spokeDirs->GetValue(idxDir+0));
        coeffArray.push_back(spokeDirs->GetValue(idxDir+1));
        coeffArray.push_back(spokeDirs->GetValue(idxDir+2));
        //coeffArray.push_back(spokeRadii->GetValue(i));
        coeffArray.push_back(0);

        // data for spokes
        radii.push_back(spokeRadii->GetValue(i));

        dirs.push_back(spokeDirs->GetValue(idxDir+0));
        dirs.push_back(spokeDirs->GetValue(idxDir+1));
        dirs.push_back(spokeDirs->GetValue(idxDir+2));

        double tempSkeletalPoint[3];
        int basePointId = i;
        spokesPolyData->GetPoint(basePointId, tempSkeletalPoint);
        skeletalPoints.push_back(tempSkeletalPoint[0]);
        skeletalPoints.push_back(tempSkeletalPoint[1]);
        skeletalPoints.push_back(tempSkeletalPoint[2]);
    }
}

void vtkSlicerSkeletalRepresentationRefinerLogic::ParseHeader(const std::string &headerFileName, int *nRows, int *nCols,
                                                               double *shift, std::string* upFileName,
                                                              std::string* downFileName, std::string* crestFileName)
{
    vtkSmartPointer<vtkXMLDataParser> parser = vtkSmartPointer<vtkXMLDataParser>::New();

    parser->SetFileName(headerFileName.c_str());
    parser->SetIgnoreCharacterData(0);

    if( parser->Parse() == 1)
    {
        vtkXMLDataElement *root = parser->GetRootElement();
        int numElements = root->GetNumberOfNestedElements();;
        for(int i = 0; i < numElements; ++i)
        {
            int r, c;
            double crestShift;
            char *pEnd;
            vtkXMLDataElement *e = root->GetNestedElement(i);
            std::string estimatePath;
            estimatePath = vtksys::SystemTools::GetFilenamePath(headerFileName) + "/";
            std::vector<std::string> components;
            components.push_back(estimatePath);

            char* eName = e->GetName();
            if(strcmp(eName, "nRows") == 0)
            {
                r = static_cast<int>(strtol(e->GetCharacterData(), &pEnd, 10));
                *nRows = r;
            }
            else if(strcmp(eName, "nCols") == 0)
            {
                c = static_cast<int>(strtol(e->GetCharacterData(), &pEnd, 10));
                *nCols = c;
            }
            else if(strcmp(eName, "upSpoke") == 0)
            {
                *upFileName = e->GetCharacterData();
                // some file paths are relative path, others are absolute path
                if(!vtksys::SystemTools::FileIsFullPath(*upFileName))
                {
                    components.push_back(*upFileName);
                    *upFileName = vtksys::SystemTools::JoinPath(components);

                }
                // change to relative path
//                *upFileName = estimatePath+ "up.vtp";
                //cout <<"header file " << headerFileName <<  "use the final up file name:" << *upFileName << endl;
            }
            else if(strcmp(eName, "downSpoke")==0)
            {
                *downFileName = e->GetCharacterData();
                if(!vtksys::SystemTools::FileIsFullPath(*downFileName))
                {
                    components.push_back(*downFileName);
                    *downFileName = vtksys::SystemTools::JoinPath(components);
                }
//                // change to relative path
//                *downFileName = estimatePath + "down.vtp";
            }
            else if(strcmp(eName, "crestSpoke") == 0)
            {
                *crestFileName = e->GetCharacterData();
                if(!vtksys::SystemTools::FileIsFullPath(*crestFileName))
                {
                    components.push_back(*crestFileName);
                    *crestFileName = vtksys::SystemTools::JoinPath(components);
                }
//                *crestFileName =estimatePath + "crest.vtp";
            }
            else if(strcmp(eName, "crestShift") == 0)
            {
                crestShift = atof(e->GetCharacterData());
                *shift = crestShift;
            }
        }
    }

}

void vtkSlicerSkeletalRepresentationRefinerLogic::UpdateHeader(const string &headerFileName,
                                                               const string &outputFilePath,
                                                               std::string *newHeaderFileName)
{
    vtkSmartPointer<vtkXMLDataParser> parser = vtkSmartPointer<vtkXMLDataParser>::New();

    parser->SetFileName(headerFileName.c_str());
    parser->SetIgnoreCharacterData(0);
    if( parser->Parse() == 1)
    {
        vtkXMLDataElement *root = parser->GetRootElement();
        int numElements = root->GetNumberOfNestedElements();;
        std::string newUpFileName, newDownFileName, newCrestFileName;

        std::string estimatePath;
        estimatePath = vtksys::SystemTools::GetFilenamePath(headerFileName) + "/";

        int nRows = 0, nCols = 0;
        for(int i = 0; i < numElements; ++i)
        {
            int r, c;
            char *pEnd;
            vtkXMLDataElement *e = root->GetNestedElement(i);
            std::string estimatePath;
            estimatePath = vtksys::SystemTools::GetFilenamePath(headerFileName) + "/";
            std::vector<std::string> components;
            components.push_back(estimatePath);
            char* eName = e->GetName();
            if(strcmp(eName, "nRows") == 0)
            {
                r = static_cast<int>(strtol(e->GetCharacterData(), &pEnd, 10));
                nRows = r;
            }
            else if(strcmp(eName, "nCols") == 0)
            {
                c = static_cast<int>(strtol(e->GetCharacterData(), &pEnd, 10));
                nCols = c;
            }
            else if(strcmp(eName, "upSpoke") == 0)
            {
                std::string oldFile(e->GetCharacterData());
                // some file paths are relative path, others are absolute path
                oldFile = vtksys::SystemTools::GetFilenameName(oldFile);
                newUpFileName = /*outputFilePath + */newFilePrefix + oldFile;
            }
            else if(strcmp(eName, "downSpoke")==0)
            {
                std::string oldFile(e->GetCharacterData());
                // some file paths are relative path, others are absolute path
                oldFile = vtksys::SystemTools::GetFilenameName(oldFile);
                newDownFileName = /*outputFilePath + */newFilePrefix + oldFile;
            }
            else if(strcmp(eName, "crestSpoke")==0)
            {
                std::string oldFile(e->GetCharacterData());
                // some file paths are relative path, others are absolute path
                oldFile = vtksys::SystemTools::GetFilenameName(oldFile);
                newCrestFileName = /*outputFilePath + */newFilePrefix + oldFile;
            }

        }
        std::stringstream output;

        output<<"<s-rep>"<<std::endl;
        output<<"  <nRows>"<<nRows<<"</nRows>"<<std::endl;
        output<<"  <nCols>"<<nCols<<"</nCols>"<<std::endl;
        output<<"  <meshType>Quad</meshType>"<< std::endl;
        output<<"  <color>"<<std::endl;
        output<<"    <red>0</red>"<<std::endl;
        output<<"    <green>0.5</green>"<<std::endl;
        output<<"    <blue>0</blue>"<<std::endl;
        output<<"  </color>"<<std::endl;
        output<<"  <isMean>False</isMean>"<<std::endl;
        output<<"  <meanStatPath/>"<<std::endl;
        output<<"  <upSpoke>"<< newUpFileName<<"</upSpoke>"<<std::endl;
        output<<"  <downSpoke>"<< newDownFileName << "</downSpoke>"<<std::endl;
        output<<"  <crestSpoke>"<< newCrestFileName << "</crestSpoke>"<<std::endl;
        output<<"</s-rep>"<<std::endl;

        std::string oldHeader = vtksys::SystemTools::GetFilenameName(headerFileName);
        oldHeader = outputFilePath + "//" + newFilePrefix + oldHeader;
        std::string header_file(oldHeader);
        std::ofstream out_file;
        out_file.open(header_file);
        out_file << output.rdbuf();
        out_file.close();
        *newHeaderFileName = header_file;
    }
}

void vtkSlicerSkeletalRepresentationRefinerLogic::ComputeDiff(double *head, double *tail, double factor, double *output)
{
    output[0] = factor * (head[0] - tail[0]);
    output[1] = factor * (head[1] - tail[1]);
    output[2] = factor * (head[2] - tail[2]);
}

double vtkSlicerSkeletalRepresentationRefinerLogic::ComputeDistance(vtkSpoke *theSpoke, double *normalMatch)
{
    // 1. Transform the boundary point to image cs. by applying [x, y, z, 1] * mTransformationMat
    double pt[3];
    theSpoke->GetBoundaryPoint(pt);

    pt[0] = pt[0] * mTransformationMat[0][0] + mTransformationMat[3][0];
    pt[1] = pt[1] * mTransformationMat[1][1] + mTransformationMat[3][1];
    pt[2] = pt[2] * mTransformationMat[2][2] + mTransformationMat[3][2];

    pt[0] /= voxelSpacing;
    pt[1] /= voxelSpacing;
    pt[2] /= voxelSpacing;

    int x = static_cast<int>(pt[0]+0.5);
    int y = static_cast<int>(pt[1]+0.5);
    int z = static_cast<int>(pt[2]+0.5);

    int maxX = static_cast<int>(1 / voxelSpacing - 1);
    int maxY = static_cast<int>(1 / voxelSpacing - 1);
    int maxZ = static_cast<int>(1 / voxelSpacing - 1);

    if(x > maxX) x = maxX;
    if(y > maxY) y = maxY;
    if(z > maxZ) z = maxZ;

    if(x < 0) x = 0;
    if(y < 0) y = 0;
    if(z < 0) z = 0;

    if(mAntiAliasedImage == nullptr)
    {
        std::cerr << "The image in this RefinerLogic instance is empty." << std::endl;
        return -10000.0;
    }
    RealImage::IndexType pixelIndex = {{x,y,z}};
    float dist = mAntiAliasedImage->GetPixel(pixelIndex);

    if(mGradDistImage == nullptr)
    {
        return -10000.0;
    }

    VectorImage::IndexType indexGrad;
    indexGrad[0] = x;
    indexGrad[1] = y;
    indexGrad[2] = z;

    VectorImage::PixelType grad = mGradDistImage->GetPixel(indexGrad);
    double normalVector[3];
    normalVector[0] = static_cast<double>(grad[0]);
    normalVector[1] = static_cast<double>(grad[1]);
    normalVector[2] = static_cast<double>(grad[2]);
    // normalize the normal vector
    vtkMath::Normalize(normalVector);

    double spokeDir[3];
    theSpoke->GetDirection(spokeDir);
    double dotProduct = vtkMath::Dot(normalVector, spokeDir);
    double distSqr = static_cast<double>(dist * dist);

    // The normal match (between [0,1]) is scaled by the distance so that the overall term is comparable
    *normalMatch = *normalMatch + distSqr * (1 - dotProduct);
    // return square of distance
    return distSqr;
}

void vtkSlicerSkeletalRepresentationRefinerLogic::Visualize(vtkPolyData *model, const std::string &modelName,
                                                            double r, double g, double b, bool isVisible)
{
    vtkMRMLScene *scene = this->GetMRMLScene();
    if(!scene)
    {
        vtkErrorMacro(" Invalid scene");
        return;
    }

    // model node
    vtkSmartPointer<vtkMRMLModelNode> modelNode;
    modelNode = vtkSmartPointer<vtkMRMLModelNode>::New();
    modelNode->SetScene(scene);
    modelNode->SetName(modelName.c_str());
    modelNode->SetAndObservePolyData(model);

    // display node
    vtkSmartPointer<vtkMRMLModelDisplayNode> displayModelNode;

    displayModelNode = vtkSmartPointer<vtkMRMLModelDisplayNode>::New();
    if(displayModelNode == nullptr)
    {
        vtkErrorMacro("displayModelNode is NULL");
        return;
    }
    displayModelNode->SetColor(r,g,b);
    displayModelNode->SetScene(scene);
    displayModelNode->SetLineWidth(2.0);
    displayModelNode->SetBackfaceCulling(0);
    displayModelNode->SetRepresentation(1);

    if(isVisible)
    {
        // make the 1st mesh after flow visible
        displayModelNode->SetVisibility(1);
    }
    else
    {
        displayModelNode->SetVisibility(0);
    }

    scene->AddNode(displayModelNode);
    modelNode->AddAndObserveDisplayNodeID(displayModelNode->GetID());

    scene->AddNode(modelNode);

}

void vtkSlicerSkeletalRepresentationRefinerLogic::HideNodesByClass(const string &className)
{
    vtkSmartPointer<vtkCollection> modelNodes = this->GetMRMLScene()->GetNodesByClass(className.c_str());
    modelNodes->InitTraversal();
    for(int i = 0; i < modelNodes->GetNumberOfItems(); i++)
    {
        vtkSmartPointer<vtkMRMLModelNode> thisModelNode = vtkMRMLModelNode::SafeDownCast(modelNodes->GetNextItemAsObject());
        vtkSmartPointer<vtkMRMLModelDisplayNode> displayNode;
        displayNode = thisModelNode->GetModelDisplayNode();
        if(displayNode == nullptr)
        {
            continue;
        }

        displayNode->SetVisibility(0);

    }

}

void vtkSlicerSkeletalRepresentationRefinerLogic::TransformSrep2ImageCS(vtkSrep *input, double mat4x4[][4])
{
    if(input->IsEmpty())
    {
        mat4x4 = nullptr;
        return;
    }
    // 1. Find the bounding box of boundary
    std::vector<vtkSpoke *> spokes = input->GetAllSpokes();
    vtkSmartPointer<vtkPoints> boundaryPts =
            vtkSmartPointer<vtkPoints>::New();
    for(size_t i = 0; i < spokes.size(); ++i)
    {
        double pt[3];
        spokes[i]->GetBoundaryPoint(pt);
        boundaryPts->InsertNextPoint(pt);
    }

    double bounds[6];
    boundaryPts->GetBounds(bounds);
    double xrange = bounds[1] - bounds[0];
    double yrange = bounds[3] - bounds[2];
    double zrange = bounds[5] - bounds[4];

    // the new bounding box keep the ratios between x, y, z
    double xrangeTrans, yrangeTrans, zrangeTrans;
    if(xrange >= yrange && xrange >= zrange)
    {
        xrangeTrans = 1.0;
        yrangeTrans = yrange / xrange;
        zrangeTrans = zrange / xrange;

    }
    else if (yrange >= xrange && yrange >= zrange)
    {
        xrangeTrans = xrange / yrange;
        yrangeTrans = 1.0;
        zrangeTrans = zrange / yrange;
    }
    else if (zrange >= xrange && zrange >= yrange)
    {
        xrangeTrans = xrange / zrange;
        yrangeTrans = yrange / zrange;
        zrangeTrans = 1.0;
    }
    else {
        xrangeTrans = 1.0;
        yrangeTrans = 1.0;
        zrangeTrans = 1.0;
    }

    // the origin of new bounding box, which is centered at (0.5, 0.5,0.5)
    double xoriginTrans, yoriginTrans, zoriginTrans;
    xoriginTrans = 0.5 - xrangeTrans / 2;
    yoriginTrans = 0.5 - yrangeTrans / 2;
    zoriginTrans = 0.5 - zrangeTrans / 2;

    // scale factors to unit cube
    mat4x4[0][0] = xrangeTrans / xrange;
    mat4x4[1][1] = yrangeTrans / yrange;
    mat4x4[2][2] = zrangeTrans / zrange;

    // tranlate amount
    mat4x4[3][0] = xoriginTrans - xrangeTrans * bounds[0] / xrange;
    mat4x4[3][1] = yoriginTrans - yrangeTrans * bounds[2] / yrange;
    mat4x4[3][2] = zoriginTrans - zrangeTrans * bounds[4] / zrange;

    // others are 0
    mat4x4[0][1] = 0; mat4x4[0][2] = 0; mat4x4[0][3] = 0;
    mat4x4[1][0] = 0; mat4x4[1][2] = 0; mat4x4[1][3] = 0;
    mat4x4[2][0] = 0; mat4x4[2][1] = 0; mat4x4[2][3] = 0;
    mat4x4[3][3] = 1; // the bottom-right corner has to be 1 to multiply with another transform matrix
}

void vtkSlicerSkeletalRepresentationRefinerLogic::ConnectImpliedBoundaryPts(int interpolationLevel,int nRows, int nCols,
                                                               const string &srepFileName,
                                                                            std::vector<vtkSpoke *> &borderSpokes,
                                                                            std::vector<vtkSpoke *> &interpolatedSpokes,
                                                                            std::vector<vtkSpoke *> &repeatedInterps,
                                                                            std::vector<vtkSpoke*>& primary,
                                                                            vtkPolyData* impliedPolyData
                                                               )
{
    std::vector<double> coeffArray, radii, dirs, skeletalPoints;
    Parse(srepFileName, coeffArray, radii, dirs, skeletalPoints);

    if(nRows == 0 || nCols == 0)
    {
        std::cerr << "The s-rep model is empty." << std::endl;
        return;
    }

    vtkSrep *srep = new vtkSrep(nRows, nCols, radii, dirs, skeletalPoints);
    if(srep->IsEmpty())
    {
        std::cerr << "The s-rep model is empty." << std::endl;
        delete srep;
        srep = nullptr;
        return;
    }
    std::vector<vtkSpoke *> pSpokes = srep->GetAllSpokes();
    for (size_t i = 0; i < pSpokes.size(); ++i) {
        vtkSpoke *s = new vtkSpoke(*pSpokes[i]);
        primary.push_back(s);
    }
    vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> quads = vtkSmartPointer<vtkCellArray>::New();

    // 1.1 interpolate and visualize for verification
    // collect neighboring spokes around corners
    vtkSlicerSkeletalRepresentationInterpolater interpolater;

    int shares = static_cast<int>(pow(2, interpolationLevel));
    double interval = static_cast<double>(1.0/ shares);
    std::vector<double> steps;

    for(int i = 0; i <= shares; ++i)
    {
        steps.push_back(i * interval);
    }
    vtkSmartPointer<vtkPolyData> surfaceMesh = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPoints> surfacePts = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> surfaceCells = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkAppendPolyData> appendFilter =
      vtkSmartPointer<vtkAppendPolyData>::New();

    int nCrestPoints = nRows*2 + (nCols-2)*2;
    int numCols = 1+static_cast<int>(floor(nRows/2)); // steps from crest point to the skeletal point
    int numPerRow = (numCols - 1) * (steps.size()-1) + 1;
    for(int r = 0; r < nCrestPoints; ++r)
    {
        // interpolate between top and bottom in a grid
        for(size_t i = 0; i < steps.size()-1; ++i)
        {
            // interpolate in the next grid
            for(int c = 0; c < numCols-1; ++c)
            {
                vtkSpoke *cornerSpokes[4];

                double  dXdu11[3], dXdv11[3],
                        dXdu12[3], dXdv12[3],
                        dXdu21[3], dXdv21[3],
                        dXdu22[3], dXdv22[3];

                // interpolate between left and right in a grid
                for(size_t j = 0; j < steps.size(); ++j)
                {
                    cornerSpokes[0] = srep->GetSpoke(r,c);
                    cornerSpokes[1] = srep->GetSpoke(r+1, c);
                    cornerSpokes[2] = srep->GetSpoke(r+1, c+1);
                    cornerSpokes[3] = srep->GetSpoke(r, c+ 1);
                    if(r == nCrestPoints - 1) {
                        ComputeDerivative(skeletalPoints, r, c, nRows, nCols, dXdu11, dXdv11);
                        ComputeDerivative(skeletalPoints, 0, c, nRows, nCols, dXdu21, dXdv21);
                        ComputeDerivative(skeletalPoints, r, c+1, nRows, nCols, dXdu12, dXdv12);
                        ComputeDerivative(skeletalPoints, 0, c+1, nRows, nCols, dXdu22, dXdv22);
                    }
                    else {
                        ComputeDerivative(skeletalPoints, r, c, nRows, nCols, dXdu11, dXdv11);
                        ComputeDerivative(skeletalPoints, r+1, c, nRows, nCols, dXdu21, dXdv21);
                        ComputeDerivative(skeletalPoints, r, c+1, nRows, nCols, dXdu12, dXdv12);
                        ComputeDerivative(skeletalPoints, r+1, c+1, nRows, nCols, dXdu22, dXdv22);
                    }
                    interpolater.SetCornerDxdu(dXdu11,
                                               dXdu21,
                                               dXdu22,
                                               dXdu12);
                    interpolater.SetCornerDxdv(dXdv11,
                                               dXdv21,
                                               dXdv22,
                                               dXdv12);

                    vtkSpoke* in1 = new vtkSpoke;
                    interpolater.Interpolate(double(steps[i]), double(steps[j]), cornerSpokes, in1);
                    repeatedInterps.push_back(in1);
                    if((c != numCols - 2 && j == steps.size() -1)) {
                        // avoid redundent
                        continue;
                    }

                    interpolatedSpokes.push_back(in1);
                    double pt[3];
                    in1->GetBoundaryPoint(pt);
                    vtkIdType currId = surfacePts->InsertNextPoint(pt);
                    int outerPtId, topRightId, rightId;
                    outerPtId = currId + 1;
                    if(r == nCrestPoints - 1 && i == steps.size() - 2 && j < steps.size() - 1)
                    {
                        // connect the gap from the last radial line to the 0th line
                        rightId = c * (steps.size() - 1) + j;
                        topRightId = rightId + 1;
                        vtkSmartPointer<vtkTriangle> topT = vtkSmartPointer<vtkTriangle>::New();

                        topT->GetPointIds()->SetId(0, currId);
                        topT->GetPointIds()->SetId(1, topRightId);
                        topT->GetPointIds()->SetId(2, outerPtId);
                        surfaceCells->InsertNextCell(topT);

                        // bot right triangle
                        vtkSmartPointer<vtkTriangle> botT = vtkSmartPointer<vtkTriangle>::New();
                        botT->GetPointIds()->SetId(0, currId);
                        botT->GetPointIds()->SetId(1, rightId); // right point
                        botT->GetPointIds()->SetId(2, topRightId);
                        surfaceCells->InsertNextCell(botT);
                    }
                    else {
                        if(c == numCols - 2 && j == steps.size() - 1)
                        {
                            // do nothing for the near top line
                            borderSpokes.push_back(in1);
                        }
                        else if (j == steps.size() - 1) {

                        }
                        else {
                            topRightId = outerPtId + numPerRow;
                            rightId = currId + numPerRow;
                            vtkSmartPointer<vtkTriangle> topT = vtkSmartPointer<vtkTriangle>::New();

                            topT->GetPointIds()->SetId(0, currId);
                            topT->GetPointIds()->SetId(1, topRightId);
                            topT->GetPointIds()->SetId(2, outerPtId);
                            surfaceCells->InsertNextCell(topT);

                            // bot right triangle
                            vtkSmartPointer<vtkTriangle> botT = vtkSmartPointer<vtkTriangle>::New();
                            botT->GetPointIds()->SetId(0, currId);
                            botT->GetPointIds()->SetId(1, rightId); // right point
                            botT->GetPointIds()->SetId(2, topRightId);
                            surfaceCells->InsertNextCell(botT);

                        }
                    }

                }

            }

        }
    }
    surfaceMesh->SetPoints(surfacePts);
    surfaceMesh->SetPolys(surfaceCells);
    surfaceMesh->Modified();
    impliedPolyData->DeepCopy(surfaceMesh);

}

void vtkSlicerSkeletalRepresentationRefinerLogic::ConnectImpliedCrest(int interpolationLevel,
                                                                   std::vector<vtkSpoke *> upInterpSpokes,
                                                                   std::vector<vtkSpoke *> downInterpSpokes,
                                                                   std::vector<vtkSpoke *> crestInterpSpokes,
                                                                   vtkPolyData* crestPoly)
{
    int shares = static_cast<int>(pow(2, interpolationLevel));
    double interval = static_cast<double>((1.0/ shares));

    vtkSmartPointer<vtkPoints> crestPoints = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkAppendPolyData> appendFilter =
      vtkSmartPointer<vtkAppendPolyData>::New();
    double ptCrest[3], ptUp[3], ptDown[3], ptInterp[3], ptSkeletal[3], du[9];

    vtkSmartPointer<vtkPolyData> interpS = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkCellArray> interpCell = vtkSmartPointer<vtkCellArray>::New();

    vtkSmartPointer<vtkPolyData> interpDownS = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkCellArray> interpDownCell = vtkSmartPointer<vtkCellArray>::New();

    vtkSmartPointer<vtkPolyData> interpCrestS = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkCellArray> interpCrestCell = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPoints> interpSPts = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkPoints> interpSBrdy = vtkSmartPointer<vtkPoints>::New();

    vtkSmartPointer<vtkPolyData>  triangleCrestPoly = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkCellArray> triangleCrestCell = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPoints>    triangleCrestPts = vtkSmartPointer<vtkPoints>::New();

    for(size_t i = 0; i < crestInterpSpokes.size(); ++i)
    {
        crestInterpSpokes[i]->GetBoundaryPoint(ptCrest);
        crestPoints->InsertNextPoint(ptCrest);
        vtkSmartPointer<vtkPoints> radialCurve = vtkSmartPointer<vtkPoints>::New();
        upInterpSpokes[i]->GetBoundaryPoint(ptUp);
        downInterpSpokes[i]->GetBoundaryPoint(ptDown);
        radialCurve->InsertNextPoint(ptUp);
        radialCurve->InsertNextPoint(ptCrest);
        radialCurve->InsertNextPoint(ptDown);
        vtkSmartPointer<vtkParametricSpline> splineRadial =
                vtkSmartPointer<vtkParametricSpline>::New();
        splineRadial->SetPoints(radialCurve);
        vtkSmartPointer<vtkParametricFunctionSource> functionSourceRadial =
                vtkSmartPointer<vtkParametricFunctionSource>::New();
        functionSourceRadial->SetParametricFunction(splineRadial);
        functionSourceRadial->Update();
        // share the base point among all other interpolated spokes
        crestInterpSpokes[i]->GetSkeletalPoint(ptSkeletal);
        // interpolate along the spline
        for(int j = 1; j < shares; ++j)
        {
            double uInterp = j * interval;
            double u[3] = {uInterp, uInterp, uInterp};
            splineRadial->Evaluate(u, ptInterp, du);
            vtkIdType id0 = interpSPts->InsertNextPoint(ptInterp);
            vtkIdType id1 = interpSPts->InsertNextPoint(ptSkeletal);
            interpSBrdy->InsertNextPoint(ptInterp);
            vtkSmartPointer<vtkLine> interpSLine = vtkSmartPointer<vtkLine>::New();
            interpSLine->GetPointIds()->SetId(0, id0);
            interpSLine->GetPointIds()->SetId(1, id1);
            interpCell->InsertNextCell(interpSLine);

        }

        // for triangular mesh of crest region. Duplicate interpolation at j = 0 and j = shares
        for(int j = 0; j <= shares; ++j) {
            double uInterp = j * interval;
            double u[3] = {uInterp, uInterp, uInterp};
            splineRadial->Evaluate(u, ptInterp, du);
            vtkIdType id_curr = triangleCrestPts->InsertNextPoint(ptInterp);
            int id_left = 0;
            if(i > 0 && j >0) {
                // there are # (shares) interpolated points across crest
                id_left = id_curr - (shares+1);
            }
            else if (j > 0){
                // the first cross crest line should connect with the last one
                id_left = j + (shares+1) * (crestInterpSpokes.size() - 1);
            }
            else continue;
            if(j > 0) {
                // right bot triangle
                vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
                triangle->GetPointIds()->SetId(0, id_left);
                triangle->GetPointIds()->SetId(1, id_curr-1);// righttop point
                triangle->GetPointIds()->SetId(2, id_curr); // right bot point
                triangleCrestCell->InsertNextCell(triangle);

                 //left top triangle
                vtkSmartPointer<vtkTriangle> triangle2 = vtkSmartPointer<vtkTriangle>::New();
                triangle2->GetPointIds()->SetId(0, id_left);
                triangle2->GetPointIds()->SetId(1, id_left-1);
                triangle2->GetPointIds()->SetId(2, id_curr-1);
                triangleCrestCell->InsertNextCell(triangle2);
            }
        }

        appendFilter->AddInputData(functionSourceRadial->GetOutput());
    }
    interpS->SetPoints(interpSPts);
    interpS->SetPolys(interpCell);
    Visualize(interpS, "Interpolated", 0, 0, 1, false);

    triangleCrestPoly->SetPoints(triangleCrestPts);
    triangleCrestPoly->SetPolys(triangleCrestCell);

    double bounds[6];
    interpS->GetBounds(bounds);

    // connect points along crest. There should be #share-1 curves in total
    for(int i = 1; i < shares; ++i)
    {
        vtkSmartPointer<vtkParametricSpline> splineAlongCrest =
                vtkSmartPointer<vtkParametricSpline>::New();
        vtkSmartPointer<vtkPoints> crestSplinePts = vtkSmartPointer<vtkPoints>::New();

        for(int j = 0; j < interpSBrdy->GetNumberOfPoints(); ++j)
        {
            if((j+1) % (shares-1) != i)
            {
                continue;
            }
            double pt[3];
            interpSBrdy->GetPoint(j, pt);
            crestSplinePts->InsertNextPoint(pt);
        }
        if(crestSplinePts->GetNumberOfPoints() > 0)
        {
            splineAlongCrest->SetPoints(crestSplinePts);
            vtkSmartPointer<vtkParametricFunctionSource> functionSourceAlong =
                    vtkSmartPointer<vtkParametricFunctionSource>::New();
            functionSourceAlong->SetParametricFunction(splineAlongCrest);
            functionSourceAlong->Update();
            appendFilter->AddInputData(functionSourceAlong->GetOutput());
        }

    }
    vtkSmartPointer<vtkParametricSpline> splineCrest =
            vtkSmartPointer<vtkParametricSpline>::New();
    splineCrest->SetPoints(crestPoints);
    vtkSmartPointer<vtkParametricFunctionSource> functionSourceCrest =
            vtkSmartPointer<vtkParametricFunctionSource>::New();
    functionSourceCrest->SetParametricFunction(splineCrest);
    functionSourceCrest->Update();
    appendFilter->AddInputData(functionSourceCrest->GetOutput());
    appendFilter->Update();

    // Remove any duplicate points.
    vtkSmartPointer<vtkCleanPolyData> cleanFilter =
            vtkSmartPointer<vtkCleanPolyData>::New();
    cleanFilter->SetInputConnection(appendFilter->GetOutputPort());
    cleanFilter->Update();
    vtkSmartPointer<vtkPolyData> crestConnectPoly = cleanFilter->GetOutput();
    crestPoly->DeepCopy(triangleCrestPoly);

}
void vtkSlicerSkeletalRepresentationRefinerLogic::ConnectImpliedCrest(int interpolationLevel,
                                                                      int nRows, int nCols,
                                                                      const std::string &crest,
                                                                      std::vector<vtkSpoke*> &upSpokes,
                                                                      std::vector<vtkSpoke*> &downSpokes,
                                                                      vtkPolyData* crestPoly)
{
    std::vector<vtkSpoke*> crestSpokes, topCrest;
    ParseCrest(crest, crestSpokes);

    std::vector<vtkSpoke *> upInterpSpokes, downInterpSpokes, crestInterpSpokes, tempInterp, reorderedCrest;
    InterpolateCrest(crestSpokes, upSpokes, interpolationLevel, nRows, nCols, crestInterpSpokes, upInterpSpokes);

    InterpolateCrest(crestSpokes, downSpokes, interpolationLevel, nRows, nCols, tempInterp, downInterpSpokes);

    vtkSmartPointer<vtkPolyData> upInterpSpokesPoly = vtkSmartPointer<vtkPolyData>::New();
    ConvertSpokes2PolyData(upInterpSpokes, upInterpSpokesPoly);
    Visualize(upInterpSpokesPoly, "up spokes", 0, 0,0, true);

    vtkSmartPointer<vtkPolyData> downInterpSpokesPoly = vtkSmartPointer<vtkPolyData>::New();
    ConvertSpokes2PolyData(downInterpSpokes, downInterpSpokesPoly);
    Visualize(downInterpSpokesPoly, "down spokes", 0, 0,0, true);

    vtkSmartPointer<vtkPolyData> crestInterpSpokesPoly = vtkSmartPointer<vtkPolyData>::New();
    ConvertSpokes2PolyData(crestInterpSpokes, crestInterpSpokesPoly);
    Visualize(crestInterpSpokesPoly, "crest spokes", 0, 0,0, true);

    // shares between up spoke to resp. down spoke
    int shares = static_cast<int>(pow(2, interpolationLevel));
    double interval = static_cast<double>((1.0/ shares));

    vtkSmartPointer<vtkPoints> crestPoints = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkAppendPolyData> appendFilter =
      vtkSmartPointer<vtkAppendPolyData>::New();
    double ptCrest[3], ptUp[3], ptDown[3], ptInterp[3], ptSkeletal[3], du[9];

    vtkSmartPointer<vtkPolyData> interpS = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkCellArray> interpCell = vtkSmartPointer<vtkCellArray>::New();

    vtkSmartPointer<vtkPolyData> interpDownS = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkCellArray> interpDownCell = vtkSmartPointer<vtkCellArray>::New();

    vtkSmartPointer<vtkPolyData> interpCrestS = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkCellArray> interpCrestCell = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPoints> interpSPts = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkPoints> interpSBrdy = vtkSmartPointer<vtkPoints>::New();

    vtkSmartPointer<vtkPolyData>  triangleCrestPoly = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkCellArray> triangleCrestCell = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPoints>    triangleCrestPts = vtkSmartPointer<vtkPoints>::New();

    for(size_t i = 0; i < crestInterpSpokes.size(); ++i)
    {
        crestInterpSpokes[i]->GetBoundaryPoint(ptCrest);
        crestPoints->InsertNextPoint(ptCrest);
        vtkSmartPointer<vtkPoints> radialCurve = vtkSmartPointer<vtkPoints>::New();
        upInterpSpokes[i]->GetBoundaryPoint(ptUp);
        downInterpSpokes[i]->GetBoundaryPoint(ptDown);
        radialCurve->InsertNextPoint(ptUp);
        radialCurve->InsertNextPoint(ptCrest);
        radialCurve->InsertNextPoint(ptDown);
        vtkSmartPointer<vtkParametricSpline> splineRadial =
                vtkSmartPointer<vtkParametricSpline>::New();
        splineRadial->SetPoints(radialCurve);
        vtkSmartPointer<vtkParametricFunctionSource> functionSourceRadial =
                vtkSmartPointer<vtkParametricFunctionSource>::New();
        functionSourceRadial->SetParametricFunction(splineRadial);
        functionSourceRadial->Update();
        // share the base point among all other interpolated spokes
        crestInterpSpokes[i]->GetSkeletalPoint(ptSkeletal);
        // interpolate along the spline
        for(int j = 1; j < shares; ++j)
        {
            double uInterp = j * interval;
            double u[3] = {uInterp, uInterp, uInterp};
            splineRadial->Evaluate(u, ptInterp, du);
            vtkIdType id0 = interpSPts->InsertNextPoint(ptInterp);
            vtkIdType id1 = interpSPts->InsertNextPoint(ptSkeletal);
            interpSBrdy->InsertNextPoint(ptInterp);
            vtkSmartPointer<vtkLine> interpSLine = vtkSmartPointer<vtkLine>::New();
            interpSLine->GetPointIds()->SetId(0, id0);
            interpSLine->GetPointIds()->SetId(1, id1);
            interpCell->InsertNextCell(interpSLine);

        }

        // for triangular mesh of crest region. Duplicate interpolation at j = 0 and j = shares
        for(int j = 0; j <= shares; ++j) {
            double uInterp = j * interval;
            double u[3] = {uInterp, uInterp, uInterp};
            splineRadial->Evaluate(u, ptInterp, du);
            vtkIdType id_curr = triangleCrestPts->InsertNextPoint(ptInterp);
            int id_left = 0;
            if(i > 0 && j >0) {
                // there are # (shares) interpolated points across crest
                id_left = id_curr - (shares+1) - 1;
            }
            else if (j > 0){
                // the first cross crest line should connect with the last one
                id_left = j-1 + (shares+1) * (crestInterpSpokes.size() - 2);
            }
            else continue;


            if(j > 0) {
                // right bot triangle
                vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
                triangle->GetPointIds()->SetId(0, id_left);
                triangle->GetPointIds()->SetId(1, id_curr-1);// righttop point
                triangle->GetPointIds()->SetId(2, id_curr); // right bot point
                triangleCrestCell->InsertNextCell(triangle);

                // left top triangle
//                vtkSmartPointer<vtkTriangle> triangle2 = vtkSmartPointer<vtkTriangle>::New();
//                triangle2->GetPointIds()->SetId(0, id_left);
//                triangle2->GetPointIds()->SetId(1, id_left-1);
//                triangle2->GetPointIds()->SetId(2, id_curr-1);
//                triangleCrestCell->InsertNextCell(triangle2);
            }
        }

        appendFilter->AddInputData(functionSourceRadial->GetOutput());
    }
    interpS->SetPoints(interpSPts);
    interpS->SetPolys(interpCell);
    Visualize(interpS, "Interpolated", 0, 0, 1, false);

    triangleCrestPoly->SetPoints(triangleCrestPts);
    triangleCrestPoly->SetPolys(triangleCrestCell);

    double bounds[6];
    interpS->GetBounds(bounds);

    // connect points along crest. There should be #share-1 curves in total
    for(int i = 1; i < shares; ++i)
    {
        vtkSmartPointer<vtkParametricSpline> splineAlongCrest =
                vtkSmartPointer<vtkParametricSpline>::New();
        vtkSmartPointer<vtkPoints> crestSplinePts = vtkSmartPointer<vtkPoints>::New();

        for(int j = 0; j < interpSBrdy->GetNumberOfPoints(); ++j)
        {
            if((j+1) % (shares-1) != i)
            {
                continue;
            }
            double pt[3];
            interpSBrdy->GetPoint(j, pt);
            crestSplinePts->InsertNextPoint(pt);
        }
        if(crestSplinePts->GetNumberOfPoints() > 0)
        {
            splineAlongCrest->SetPoints(crestSplinePts);
            vtkSmartPointer<vtkParametricFunctionSource> functionSourceAlong =
                    vtkSmartPointer<vtkParametricFunctionSource>::New();
            functionSourceAlong->SetParametricFunction(splineAlongCrest);
            functionSourceAlong->Update();
            appendFilter->AddInputData(functionSourceAlong->GetOutput());
        }

    }
    vtkSmartPointer<vtkParametricSpline> splineCrest =
            vtkSmartPointer<vtkParametricSpline>::New();
    splineCrest->SetPoints(crestPoints);
    vtkSmartPointer<vtkParametricFunctionSource> functionSourceCrest =
            vtkSmartPointer<vtkParametricFunctionSource>::New();
    functionSourceCrest->SetParametricFunction(splineCrest);
    functionSourceCrest->Update();
    appendFilter->AddInputData(functionSourceCrest->GetOutput());
    appendFilter->Update();

    // Remove any duplicate points.
    vtkSmartPointer<vtkCleanPolyData> cleanFilter =
            vtkSmartPointer<vtkCleanPolyData>::New();
    cleanFilter->SetInputConnection(appendFilter->GetOutputPort());
    cleanFilter->Update();
    vtkSmartPointer<vtkPolyData> crestConnectPoly = cleanFilter->GetOutput();
    //Visualize(crestConnectPoly, "Implied crest", 0, 1, 1);
    crestPoly->DeepCopy(triangleCrestPoly);
}
void vtkSlicerSkeletalRepresentationRefinerLogic::ConnectFoldCurve(const std::vector<vtkSpoke *> &edgeSpokes,
                                                                   vtkPoints *foldCurvePts, vtkCellArray *foldCurveCell)
{
    if(edgeSpokes.empty())
    {
        return;
    }
    vtkIdType id1 = 0;
    for (size_t i = 0; i < edgeSpokes.size()-1; i+=2) {
        double pt0[3], pt1[3];
        edgeSpokes[i]->GetSkeletalPoint(pt0);
        edgeSpokes[i+1]->GetSkeletalPoint(pt1);
        vtkIdType id0 = foldCurvePts->InsertNextPoint(pt0);
        id1 = foldCurvePts->InsertNextPoint(pt1);

        vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
        line->GetPointIds()->SetId(0, id0);
        line->GetPointIds()->SetId(1, id1);
        foldCurveCell->InsertNextCell(line);
    }
    // connect first and last point to close this curve
    vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
    line->GetPointIds()->SetId(0, id1);
    line->GetPointIds()->SetId(1, 0);
    foldCurveCell->InsertNextCell(line);
}

std::vector<vtkSpoke*>& vtkSlicerSkeletalRepresentationRefinerLogic::RefinePartOfSpokes(const string &srepFileName, double stepSize, double endCriterion, int maxIter)
{
    mCoeffArray.clear();
    std::vector<double> radii, dirs, skeletalPoints;
    Parse(srepFileName, mCoeffArray, radii, dirs, skeletalPoints);

    vtkSrep *srep = new vtkSrep(mNumRows, mNumCols, radii, dirs, skeletalPoints);
    if(srep->IsEmpty())
    {
        std::cerr << "The s-rep model is empty." << std::endl;
        delete srep;
        srep = nullptr;
        std::vector<vtkSpoke*> emptyRet;
        return emptyRet;
    }

    // total number of parameters that need to optimize
    size_t paramDim = mCoeffArray.size();
    double coeff[paramDim];
    for(size_t i = 0; i < paramDim; ++i)
    {
        coeff[i] = mCoeffArray[i];
    }

    mSrep = srep;
    auto orig_spokes = srep->GetAllSpokes();
    vtkSmartPointer<vtkPolyData> origSrep = vtkSmartPointer<vtkPolyData>::New();
    ConvertSpokes2PolyData(orig_spokes, origSrep);

    Visualize(origSrep, "Before refinement", 1, 0, 0);

    mFirstCost = true;
    // 2. Invoke newuoa to optimize
    min_newuoa(static_cast<int>(paramDim), coeff, *this, stepSize, endCriterion, maxIter);

    // Re-evaluate the cost
    mFirstCost = true;
    EvaluateObjectiveFunction(coeff);

    // 3. Visualize the refined srep
    srep->Refine(coeff);
    vtkSmartPointer<vtkPolyData> refinedSrep = vtkSmartPointer<vtkPolyData>::New();
    ConvertSpokes2PolyData(srep->GetAllSpokes(), refinedSrep);
    Visualize(refinedSrep, "Refined interior spokes", 0, 1, 1);
    // write to vtp file
    std::string outputFile(mOutputPath);
    std::string fileName = vtksys::SystemTools::GetFilenameName(srepFileName);
    outputFile = outputFile + "//" + newFilePrefix + fileName;
    SaveSpokes2Vtp(orig_spokes, outputFile);
    return srep->GetAllSpokes();
}

void vtkSlicerSkeletalRepresentationRefinerLogic::RefineCrestSpokes(const string &crest, double stepSize,
                                                                    double vtkNotUsed(endCriterion), int maxIter)
{
    // Show original crest spokes
    std::vector<vtkSpoke*> crestSpokes, topCrest;
    ParseCrest(crest, crestSpokes);

    vtkSmartPointer<vtkPolyData> crestSrep = vtkSmartPointer<vtkPolyData>::New();
    ConvertSpokes2PolyData(crestSpokes, crestSrep);

    Visualize(crestSrep, "Crest before refinement", 1, 1, 0);

    vtkSmartPointer<vtkPolyDataReader> meshReader = vtkSmartPointer<vtkPolyDataReader>::New();
    meshReader->SetFileName(mTargetMeshFilePath.c_str());
    meshReader->Update();
    vtkSmartPointer<vtkPolyData> mesh = meshReader->GetOutput();
    vtkSmartPointer<vtkImplicitPolyDataDistance> implicitPolyDataDistance =
      vtkSmartPointer<vtkImplicitPolyDataDistance>::New();
    implicitPolyDataDistance->SetInput(mesh);
    for (size_t i = 0; i < crestSpokes.size(); ++i) {
        OptimizeCrestSpokeLength(implicitPolyDataDistance, crestSpokes[i], stepSize, maxIter);
    }
    // set crest radii to the reciprocal of crest curvature
    vtkSmartPointer<vtkCurvatures> curvaturesFilter =
        vtkSmartPointer<vtkCurvatures>::New();
    curvaturesFilter->SetInputData(mesh);
    curvaturesFilter->SetCurvatureTypeToMaximum();
    curvaturesFilter->Update();

    vtkSmartPointer<vtkDoubleArray> MC =
        vtkDoubleArray::SafeDownCast(curvaturesFilter->GetOutput()->GetPointData()->GetArray("Maximum_Curvature"));

    if(MC == nullptr) {
        std::cerr << "error in getting max curvature" << std::endl;
        return;
    }

    curvaturesFilter->SetCurvatureTypeToMinimum();
    curvaturesFilter->Update();

    vtkSmartPointer<vtkDoubleArray> MinC =
        vtkDoubleArray::SafeDownCast(curvaturesFilter->GetOutput()->GetPointData()->GetArray("Minimum_Curvature"));
    if(MinC == nullptr)
    {
        std::cout << "error in getting min curvature" << std::endl;
        return;
    }
    // find the nearest point id on the mesh
    vtkSmartPointer<vtkPointLocator> locator = vtkPointLocator::New();
    locator->SetDataSet(mesh);
    locator->BuildLocator();
    for (size_t i = 0; i < crestSpokes.size(); ++i) {
        double bdryPt[3];
        crestSpokes[i]->GetBoundaryPoint(bdryPt);
        vtkIdType idNearest = locator->FindClosestPoint(bdryPt);
        double curr_max = MC->GetValue(idNearest);
        double curr_min = MinC->GetValue(idNearest);
        double rCrest = 1 / max(abs(curr_max), abs(curr_min));
        double rDiff = crestSpokes[i]->GetRadius() - rCrest;
        if(rDiff <= 0) continue;
        // move skeletal point of this crest outward by rDiff
        double skeletalPt[3], u[3];
        crestSpokes[i]->GetDirection(u);
        crestSpokes[i]->GetSkeletalPoint(skeletalPt);
        crestSpokes[i]->SetSkeletalPoint(skeletalPt[0] + u[0] * rDiff, skeletalPt[1] + u[1] * rDiff, skeletalPt[2] + u[2] * rDiff);
        crestSpokes[i]->SetRadius(rCrest);
    }
    // 3. Visualize the refined srep
    vtkSmartPointer<vtkPolyData> refinedSrep = vtkSmartPointer<vtkPolyData>::New();
    ConvertSpokes2PolyData(crestSpokes, refinedSrep);
    Visualize(refinedSrep, "Refined crest", 0, 1, 1);

    // write to vtp file
    std::string outputFile(mOutputPath);
    std::string fileName = vtksys::SystemTools::GetFilenameName(crest);
    outputFile = outputFile + "//" + newFilePrefix + fileName;
    SaveSpokes2Vtp(crestSpokes, outputFile);
    if(mSrep != nullptr)
    {
        delete mSrep;
        mSrep = nullptr;
    }
}

double vtkSlicerSkeletalRepresentationRefinerLogic::TotalDistOfLeftTopSpoke(vtkSrep *tempSrep,
                                                                            double u, double v,
                                                                            int r, int c,
                                                                            double *normalMatch)
{
    vtkSlicerSkeletalRepresentationInterpolater interpolater;
    vtkSpoke *cornerSpokes[4];
    double imageDist = 0.0;
    cornerSpokes[0] = tempSrep->GetSpoke(r, c);
    cornerSpokes[1] = tempSrep->GetSpoke(r+1, c);
    cornerSpokes[2] = tempSrep->GetSpoke(r+1, c+1);
    cornerSpokes[3] = tempSrep->GetSpoke(r, c+1);
    double  dXdu11[3], dXdv11[3],
            dXdu12[3], dXdv12[3],
            dXdu21[3], dXdv21[3],
            dXdu22[3], dXdv22[3];
    std::vector<double> skeletalPts = tempSrep->GetAllSkeletalPoints();
    int nRows = tempSrep->GetNumRows();
    int nCols = tempSrep->GetNumCols();
    ComputeDerivative(skeletalPts, r, c, nRows, nCols, dXdu11, dXdv11);
    ComputeDerivative(skeletalPts, r+1, c, nRows, nCols, dXdu21, dXdv21);
    ComputeDerivative(skeletalPts, r+1, c+1, nRows, nCols, dXdu22, dXdv22);
    ComputeDerivative(skeletalPts, r, c+1, nRows, nCols, dXdu12, dXdv12);

    interpolater.SetCornerDxdu(dXdu11,
                               dXdu21,
                               dXdu22,
                               dXdu12);
    interpolater.SetCornerDxdv(dXdv11,
                               dXdv21,
                               dXdv22,
                               dXdv12);
    vtkSpoke interpolatedSpoke;
    interpolater.Interpolate(u, v, cornerSpokes, &interpolatedSpoke);

    // compute the ssd for this interpolated spoke
    imageDist += ComputeDistance(&interpolatedSpoke, normalMatch);

    return imageDist;
}

double vtkSlicerSkeletalRepresentationRefinerLogic::TotalDistOfRightTopSpoke(vtkSrep *tempSrep,
                                                                             double u, double v,
                                                                             int r, int c,
                                                                             double *normalMatch)
{
    vtkSlicerSkeletalRepresentationInterpolater interpolater;
    vtkSpoke *cornerSpokes[4];
    double imageDist = 0.0;
    cornerSpokes[0] = tempSrep->GetSpoke(r, c-1);
    cornerSpokes[1] = tempSrep->GetSpoke(r+1, c-1);
    cornerSpokes[2] = tempSrep->GetSpoke(r+1, c);
    cornerSpokes[3] = tempSrep->GetSpoke(r, c);
    double  dXdu11[3], dXdv11[3],
            dXdu12[3], dXdv12[3],
            dXdu21[3], dXdv21[3],
            dXdu22[3], dXdv22[3];
    std::vector<double> skeletalPts = tempSrep->GetAllSkeletalPoints();
    int nRows = tempSrep->GetNumRows();
    int nCols = tempSrep->GetNumCols();
    ComputeDerivative(skeletalPts, r, c-1, nRows, nCols, dXdu11, dXdv11);
    ComputeDerivative(skeletalPts, r+1, c-1, nRows, nCols, dXdu21, dXdv21);
    ComputeDerivative(skeletalPts, r+1, c, nRows, nCols, dXdu22, dXdv22);
    ComputeDerivative(skeletalPts, r, c, nRows, nCols, dXdu12, dXdv12);

    interpolater.SetCornerDxdu(dXdu11,
                               dXdu21,
                               dXdu22,
                               dXdu12);
    interpolater.SetCornerDxdv(dXdv11,
                               dXdv21,
                               dXdv22,
                               dXdv12);
    vtkSpoke interpolatedSpoke;
    interpolater.Interpolate(u, v, cornerSpokes, &interpolatedSpoke);

    // compute the ssd & normal match for this interpolated spoke
    imageDist += ComputeDistance(&interpolatedSpoke, normalMatch);
    return imageDist;
}

double vtkSlicerSkeletalRepresentationRefinerLogic::TotalDistOfLeftBotSpoke(vtkSrep *tempSrep,
                                                                            double u, double v,
                                                                            int r, int c,
                                                                            double *normalMatch)
{
    vtkSlicerSkeletalRepresentationInterpolater interpolater;
    vtkSpoke *cornerSpokes[4];
    double imageDist = 0.0;
    cornerSpokes[0] = tempSrep->GetSpoke(r-1, c);
    cornerSpokes[1] = tempSrep->GetSpoke(r, c);
    cornerSpokes[2] = tempSrep->GetSpoke(r, c+1);
    cornerSpokes[3] = tempSrep->GetSpoke(r-1, c+1);
    double  dXdu11[3], dXdv11[3],
            dXdu12[3], dXdv12[3],
            dXdu21[3], dXdv21[3],
            dXdu22[3], dXdv22[3];
    std::vector<double> skeletalPts = tempSrep->GetAllSkeletalPoints();
    int nRows = tempSrep->GetNumRows();
    int nCols = tempSrep->GetNumCols();
    ComputeDerivative(skeletalPts, r-1, c, nRows, nCols, dXdu11, dXdv11);
    ComputeDerivative(skeletalPts, r, c, nRows, nCols, dXdu21, dXdv21);
    ComputeDerivative(skeletalPts, r, c+1, nRows, nCols, dXdu22, dXdv22);
    ComputeDerivative(skeletalPts, r-1, c+1, nRows, nCols, dXdu12, dXdv12);

    interpolater.SetCornerDxdu(dXdu11,
                               dXdu21,
                               dXdu22,
                               dXdu12);
    interpolater.SetCornerDxdv(dXdv11,
                               dXdv21,
                               dXdv22,
                               dXdv12);
    vtkSpoke interpolatedSpoke;
    interpolater.Interpolate(u, v, cornerSpokes, &interpolatedSpoke);

    // compute the ssd for this interpolated spoke
    imageDist += ComputeDistance(&interpolatedSpoke, normalMatch);
    return imageDist;
}

double vtkSlicerSkeletalRepresentationRefinerLogic::TotalDistOfRightBotSpoke(vtkSrep *tempSrep,
                                                                             double u, double v,
                                                                             int r, int c,
                                                                             double *normalMatch)
{
    vtkSlicerSkeletalRepresentationInterpolater interpolater;
    vtkSpoke *cornerSpokes[4];
    double imageDist = 0.0;
    cornerSpokes[0] = tempSrep->GetSpoke(r-1, c-1);
    cornerSpokes[1] = tempSrep->GetSpoke(r, c-1);
    cornerSpokes[2] = tempSrep->GetSpoke(r, c);
    cornerSpokes[3] = tempSrep->GetSpoke(r-1, c);
    double  dXdu11[3], dXdv11[3],
            dXdu12[3], dXdv12[3],
            dXdu21[3], dXdv21[3],
            dXdu22[3], dXdv22[3];
    std::vector<double> skeletalPts = tempSrep->GetAllSkeletalPoints();
    int nRows = tempSrep->GetNumRows();
    int nCols = tempSrep->GetNumCols();
    ComputeDerivative(skeletalPts, r-1, c-1, nRows, nCols, dXdu11, dXdv11);
    ComputeDerivative(skeletalPts, r, c-1, nRows, nCols, dXdu21, dXdv21);
    ComputeDerivative(skeletalPts, r, c, nRows, nCols, dXdu22, dXdv22);
    ComputeDerivative(skeletalPts, r-1, c, nRows, nCols, dXdu12, dXdv12);

    interpolater.SetCornerDxdu(dXdu11,
                               dXdu21,
                               dXdu22,
                               dXdu12);
    interpolater.SetCornerDxdv(dXdv11,
                               dXdv21,
                               dXdv22,
                               dXdv12);
    vtkSpoke interpolatedSpoke;
    interpolater.Interpolate(u, v, cornerSpokes, &interpolatedSpoke);

    // compute the ssd for this interpolated spoke
    imageDist += ComputeDistance(&interpolatedSpoke, normalMatch);
    return imageDist;
}

double vtkSlicerSkeletalRepresentationRefinerLogic::ComputeRSradPenalty(vtkSrep *input)
{
    double penalty = 0.0;
    // Interpolate
    if(input->IsEmpty())
    {
        std::cerr << "The s-rep model is empty in computing rSrad." << std::endl;
        return 0.0;
    }
    // 1.1 interpolate and visualize for verification
    // collect neighboring spokes around corners
    vtkSlicerSkeletalRepresentationInterpolater interpolater;
    int nRows = input->GetNumRows();
    int nCols = input->GetNumCols();
    int nCrestPoints = nRows*2 + (nCols-2)*2;
    int numCols = 1+static_cast<int>(floor(nRows/2)); // steps from crest point to the skeletal point


    std::vector<vtkSpoke*> interpolatedSpokes;
    std::vector<vtkSpoke*> primarySpokes;
    for(int r = 0; r < nCrestPoints; ++r)
    {
        for(int c = 0; c < numCols; ++c)
        {
            vtkSpoke* thisSpoke = input->GetSpoke(r, c);

            std::vector<vtkSpoke*> neighborU, neighborV;
            if(r == 0 && c == 0)
            {
                // top left corner
                FindTopLeftNeigbors(r, c, input, neighborU, neighborV);
                thisSpoke->SetNeighborU(neighborU, true);
                thisSpoke->SetNeighborV(neighborV, true);
            }
            else if(r == 0 && c == numCols - 1)
            {
                // top right corner
                FindTopRightNeigbors(r, c, input, neighborU, neighborV);
                thisSpoke->SetNeighborU(neighborU, true);
                thisSpoke->SetNeighborV(neighborV, false);
            }
            else if(r == 0)
            {
                // top edge
                FindTopRightNeigbors(r, c, input, neighborU, neighborV);
                FindTopLeftNeigbors(r, c, input, neighborU, neighborV);
                neighborU.pop_back();
                thisSpoke->SetNeighborU(neighborU, true);
                thisSpoke->SetNeighborV(neighborV, false);
            }
            else if(r == nCrestPoints - 1 && c == 0)
            {
                // left bot corner
                FindBotLeftNeigbors(r, c, input, neighborU, neighborV);
                thisSpoke->SetNeighborU(neighborU, false);
                thisSpoke->SetNeighborV(neighborV, true);
            }
            else if(r == nCrestPoints - 1 && c == numCols -1)
            {
                // right bot corner
                FindBotRightNeigbors(r, c, input, neighborU, neighborV);
                thisSpoke->SetNeighborU(neighborU, false);
                thisSpoke->SetNeighborV(neighborV, false);
            }
            else if(r == nCrestPoints - 1)
            {
                // bot edge
                FindBotRightNeigbors(r, c, input, neighborU, neighborV);
                FindBotLeftNeigbors(r, c, input, neighborU, neighborV);
                neighborU.pop_back();
                thisSpoke->SetNeighborU(neighborU, false);
                thisSpoke->SetNeighborV(neighborV, false);
            }
            else if(c == 0)
            {
                // left edge
                FindBotLeftNeigbors(r, c, input, neighborU, neighborV);
                FindTopLeftNeigbors(r, c, input, neighborU, neighborV);
                neighborV.pop_back();
                thisSpoke->SetNeighborU(neighborU, false);
                thisSpoke->SetNeighborV(neighborV, true);
            }
            else if(c == numCols - 1)
            {
                // right edge
                FindBotRightNeigbors(r, c, input, neighborU, neighborV);
                FindTopRightNeigbors(r, c, input, neighborU, neighborV);
                neighborV.pop_back();
                thisSpoke->SetNeighborU(neighborU, false);
                thisSpoke->SetNeighborV(neighborV, false);
            }
            else {
                // interior
                FindBotRightNeigbors(r, c, input, neighborU, neighborV);
                FindTopLeftNeigbors(r, c, input, neighborU, neighborV);
                thisSpoke->SetNeighborU(neighborU, false);
                thisSpoke->SetNeighborV(neighborV, false);
            }

            primarySpokes.push_back(thisSpoke);
        }

    }

    //compute the penalty
    double step = mInterpolatePositions[0].second;
    for(size_t i = 0; i < primarySpokes.size(); ++i)
    {
        double thisPenalty = primarySpokes[i]->GetRSradPenalty(step);
        penalty += thisPenalty;
    }
    return penalty;
}

void vtkSlicerSkeletalRepresentationRefinerLogic::FindTopLeftNeigbors(int r, int c,
                                                               vtkSrep* input,
                                                               std::vector<vtkSpoke *> &neighborU,
                                                               std::vector<vtkSpoke *> &neighborV)
{
    vtkSlicerSkeletalRepresentationInterpolater interpolater;
    int nRows = input->GetNumRows();
    int nCols = input->GetNumCols();
    vtkSpoke *cornerSpokes[4];

    double  dXdu11[3], dXdv11[3],
            dXdu12[3], dXdv12[3],
            dXdu21[3], dXdv21[3],
            dXdu22[3], dXdv22[3];
    cornerSpokes[0] = input->GetSpoke(r, c);
    cornerSpokes[1] = input->GetSpoke(r+1, c);
    cornerSpokes[2] = input->GetSpoke(r+1, c+1);
    cornerSpokes[3] = input->GetSpoke(r, c+ 1);
    std::vector<double> skeletalPts = input->GetAllSkeletalPoints();
    ComputeDerivative(skeletalPts, r, c, nRows, nCols, dXdu11, dXdv11);
    ComputeDerivative(skeletalPts, r+1, c, nRows, nCols, dXdu21, dXdv21);
    ComputeDerivative(skeletalPts, r, c+1, nRows, nCols, dXdu12, dXdv12);
    ComputeDerivative(skeletalPts, r+1, c+1, nRows, nCols, dXdu22, dXdv22);

    interpolater.SetCornerDxdu(dXdu11,
                               dXdu21,
                               dXdu22,
                               dXdu12);
    interpolater.SetCornerDxdv(dXdv11,
                               dXdv21,
                               dXdv22,
                               dXdv12);

    vtkSpoke* in1 = new vtkSpoke;
    vtkSpoke* in2 = new vtkSpoke;
    double stepV = mInterpolatePositions[0].second;
    double stepU = stepV;
    interpolater.Interpolate(stepU, 0, cornerSpokes, in1);
    neighborU.push_back(in1);
    interpolater.Interpolate(0, stepV, cornerSpokes, in2);
    neighborV.push_back(in2);
}

void vtkSlicerSkeletalRepresentationRefinerLogic::FindTopRightNeigbors(int r, int c, vtkSrep *input, std::vector<vtkSpoke *> &neighborU, std::vector<vtkSpoke *> &neighborV)
{
    vtkSlicerSkeletalRepresentationInterpolater interpolater;
    int nRows = input->GetNumRows();
    int nCols = input->GetNumCols();
    vtkSpoke *cornerSpokes[4];
    std::vector<double> skeletalPts = input->GetAllSkeletalPoints();
    double  dXdu11[3], dXdv11[3],
            dXdu12[3], dXdv12[3],
            dXdu21[3], dXdv21[3],
            dXdu22[3], dXdv22[3];
    cornerSpokes[0] = input->GetSpoke(r, c- 1);
    cornerSpokes[1] = input->GetSpoke(r+1, c-1);
    cornerSpokes[2] = input->GetSpoke(r+1, c);
    cornerSpokes[3] = input->GetSpoke(r, c);
    ComputeDerivative(skeletalPts, r, c-1, nRows, nCols, dXdu11, dXdv11);
    ComputeDerivative(skeletalPts, r+1, c-1, nRows, nCols, dXdu21, dXdv21);
    ComputeDerivative(skeletalPts, r+1, c, nRows, nCols, dXdu12, dXdv12);
    ComputeDerivative(skeletalPts, r, c, nRows, nCols, dXdu22, dXdv22);

    interpolater.SetCornerDxdu(dXdu11,
                               dXdu21,
                               dXdu22,
                               dXdu12);
    interpolater.SetCornerDxdv(dXdv11,
                               dXdv21,
                               dXdv22,
                               dXdv12);

    vtkSpoke* in1 = new vtkSpoke;
    vtkSpoke* in2 = new vtkSpoke;
    double stepV = mInterpolatePositions[0].second;
    double stepU = stepV;
    interpolater.Interpolate(stepU, 1, cornerSpokes, in1);
    neighborU.push_back(in1);
    interpolater.Interpolate(0, 1-stepV, cornerSpokes, in2);
    neighborV.push_back(in2);
}

void vtkSlicerSkeletalRepresentationRefinerLogic::FindBotLeftNeigbors(int r, int c,
                                                                      vtkSrep *input,
                                                                      std::vector<vtkSpoke *> &neighborU,
                                                                      std::vector<vtkSpoke *> &neighborV)
{
    vtkSlicerSkeletalRepresentationInterpolater interpolater;
    int nRows = input->GetNumRows();
    int nCols = input->GetNumCols();
    vtkSpoke *cornerSpokes[4];
    std::vector<double> skeletalPts = input->GetAllSkeletalPoints();
    double  dXdu11[3], dXdv11[3],
            dXdu12[3], dXdv12[3],
            dXdu21[3], dXdv21[3],
            dXdu22[3], dXdv22[3];
    cornerSpokes[0] = input->GetSpoke(r-1, c);
    cornerSpokes[1] = input->GetSpoke(r, c);
    cornerSpokes[2] = input->GetSpoke(r, c+1);
    cornerSpokes[3] = input->GetSpoke(r-1, c+1);
    ComputeDerivative(skeletalPts, r-1, c, nRows, nCols, dXdu11, dXdv11);
    ComputeDerivative(skeletalPts, r, c, nRows, nCols, dXdu21, dXdv21);
    ComputeDerivative(skeletalPts, r, c+1, nRows, nCols, dXdu12, dXdv12);
    ComputeDerivative(skeletalPts, r-1, c+1, nRows, nCols, dXdu22, dXdv22);

    interpolater.SetCornerDxdu(dXdu11,
                               dXdu21,
                               dXdu22,
                               dXdu12);
    interpolater.SetCornerDxdv(dXdv11,
                               dXdv21,
                               dXdv22,
                               dXdv12);

    vtkSpoke* in1 = new vtkSpoke;
    vtkSpoke* in2 = new vtkSpoke;
    double stepV = mInterpolatePositions[0].second;
    double stepU = stepV;
    interpolater.Interpolate(1-stepU, 0, cornerSpokes, in1);
    neighborU.push_back(in1);
    interpolater.Interpolate(0, stepV, cornerSpokes, in2);
    neighborV.push_back(in2);
}

void vtkSlicerSkeletalRepresentationRefinerLogic::FindBotRightNeigbors(int r, int c,
                                                                       vtkSrep *input,
                                                                       std::vector<vtkSpoke *> &neighborU,
                                                                       std::vector<vtkSpoke *> &neighborV)
{
    vtkSlicerSkeletalRepresentationInterpolater interpolater;
    int nRows = input->GetNumRows();
    int nCols = input->GetNumCols();
    vtkSpoke *cornerSpokes[4];
    std::vector<double> skeletalPts = input->GetAllSkeletalPoints();
    double  dXdu11[3], dXdv11[3],
            dXdu12[3], dXdv12[3],
            dXdu21[3], dXdv21[3],
            dXdu22[3], dXdv22[3];
    cornerSpokes[0] = input->GetSpoke(r-1, c-1);
    cornerSpokes[1] = input->GetSpoke(r, c-1);
    cornerSpokes[2] = input->GetSpoke(r, c);
    cornerSpokes[3] = input->GetSpoke(r-1, c);
    ComputeDerivative(skeletalPts, r-1, c-1, nRows, nCols, dXdu11, dXdv11);
    ComputeDerivative(skeletalPts, r, c-1, nRows, nCols, dXdu21, dXdv21);
    ComputeDerivative(skeletalPts, r, c, nRows, nCols, dXdu12, dXdv12);
    ComputeDerivative(skeletalPts, r-1, c, nRows, nCols, dXdu22, dXdv22);

    interpolater.SetCornerDxdu(dXdu11,
                               dXdu21,
                               dXdu22,
                               dXdu12);
    interpolater.SetCornerDxdv(dXdv11,
                               dXdv21,
                               dXdv22,
                               dXdv12);

    vtkSpoke* in1 = new vtkSpoke;
    vtkSpoke* in2 = new vtkSpoke;
    double stepV = mInterpolatePositions[0].second;
    double stepU = stepV;
    interpolater.Interpolate(1-stepU, 1, cornerSpokes, in1);
    neighborU.push_back(in1);
    interpolater.Interpolate(1, 1-stepV, cornerSpokes, in2);
    neighborV.push_back(in2);
}
void vtkSlicerSkeletalRepresentationRefinerLogic::ParseCrest(const string &crestFileName,
                                                             std::vector<vtkSpoke*> &crestSpokes)
{
    vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
    reader->SetFileName(crestFileName.c_str());
    reader->Update();


    vtkSmartPointer<vtkPolyData> spokesPolyData = reader->GetOutput();
//    Visualize(spokesPolyData, "Fold curve", 1, 1, 0);
    vtkSmartPointer<vtkPointData> spokesPointData = spokesPolyData->GetPointData();
    int numOfArrays = spokesPointData->GetNumberOfArrays();
    vtkIdType numOfSpokes = spokesPolyData->GetNumberOfPoints();

    if(numOfSpokes == 0 || numOfArrays == 0)
    {
        return;
    }

    // including Ux, Uy, Uz, r
    vtkSmartPointer<vtkDoubleArray> spokeRadii = vtkDoubleArray::SafeDownCast(spokesPointData->GetArray("spokeLength"));
    vtkSmartPointer<vtkDoubleArray> spokeDirs = vtkDoubleArray::SafeDownCast(spokesPointData->GetArray("spokeDirection"));

    for(int i = 0; i < numOfSpokes; ++i)
    {
        int idxDir = i * 3; // Ux, Uy, Uz

        vtkSpoke* crestSpoke = new vtkSpoke;
        crestSpoke->SetRadius(spokeRadii->GetValue(i));
        double u[3];
        u[0] = spokeDirs->GetValue(idxDir+0); u[1] = spokeDirs->GetValue(idxDir+1); u[2] = spokeDirs->GetValue(idxDir+2);
        crestSpoke->SetDirection(u);

        double tempSkeletalPoint[3];
        spokesPolyData->GetPoint(i, tempSkeletalPoint);
        crestSpoke->SetSkeletalPoint(tempSkeletalPoint[0], tempSkeletalPoint[1], tempSkeletalPoint[2]);
        crestSpokes.push_back(crestSpoke);
    }
}

// interpolate the crest along clock-wise direction
void vtkSlicerSkeletalRepresentationRefinerLogic::InterpolateCrest(std::vector<vtkSpoke *> &crestSpoke,
                                                                   std::vector<vtkSpoke *> &interiorSpokes,
                                                                   int interpolationLevel, int nRows,
                                                                   int nCols, std::vector<vtkSpoke*> &crest,
                                                                   std::vector<vtkSpoke*> &interior)
{
       std::vector<double> skeletalPts;
       for(size_t i = 0; i < interiorSpokes.size(); ++i)
       {
           double pt[3];
           interiorSpokes[i]->GetSkeletalPoint(pt);
           skeletalPts.push_back(pt[0]);
           skeletalPts.push_back(pt[1]);
           skeletalPts.push_back(pt[2]);
       }
       vtkSlicerSkeletalRepresentationInterpolater interpolater;

       int shares = static_cast<int>(pow(2, interpolationLevel));
       double interval = static_cast<double>((1.0/ shares));
       std::vector<double> steps;

       for(int i = 0; i <= shares; ++i)
       {
           steps.push_back(i * interval);
       }

       vtkSpoke *cornerSpokes[4];

       double  dXdu11[3], dXdv11[3],
               dXdu12[3], dXdv12[3],
               dXdu21[3], dXdv21[3],
               dXdu22[3], dXdv22[3];
       int nCrestPoints = nRows*2 + (nCols-2)*2;
       int numCols = 1+static_cast<int>(floor(nRows/2)); // steps from crest point to the skeletal point

       for(int i = 0; i < nCrestPoints; ++i) {
           int id_next_row = i + 1;
           if(i == nCrestPoints - 1) id_next_row = 0;
           ComputeDerivative(skeletalPts, i, numCols-1, nRows, nCols, dXdu21, dXdv21);
           ComputeDerivative(skeletalPts, id_next_row, numCols-1, nRows, nCols, dXdu22, dXdv22);
           dXdu11[0] = dXdu21[0];
           dXdu11[1] = dXdu21[1];
           dXdu11[2] = dXdu21[2];
           dXdv11[0] = dXdv21[0];
           dXdv11[1] = dXdv21[1];
           dXdv11[2] = dXdv21[2];

           dXdu12[0] = dXdu22[0];
           dXdu12[1] = dXdu22[1];
           dXdu12[2] = dXdu22[2];
           dXdv12[0] = dXdv22[0];
           dXdv12[1] = dXdv22[1];
           dXdv12[2] = dXdv22[2];
           interpolater.SetCornerDxdu(dXdu11,
                                      dXdu21,
                                      dXdu22,
                                      dXdu12);
           interpolater.SetCornerDxdv(dXdv11,
                                      dXdv21,
                                      dXdv22,
                                      dXdv12);
           size_t sti = static_cast<size_t>(i);
           cornerSpokes[0] = crestSpoke[sti];
           cornerSpokes[1] = interiorSpokes[(i+1) * numCols - 1];
           cornerSpokes[2] = NULL;
           if(i == nCrestPoints - 1) {
               cornerSpokes[2] = interiorSpokes[numCols - 1];
               cornerSpokes[3] = crestSpoke[0];
           }
           else {
               cornerSpokes[2] = interiorSpokes[(i+2) * numCols - 1];
               cornerSpokes[3] = crestSpoke[sti+1];
           }

           for(size_t sj = 0; sj < steps.size()-1; ++sj)
           {
               vtkSpoke* in1 = new vtkSpoke;
               interpolater.Interpolate(0.0, double(steps[sj]), cornerSpokes, in1);
               crest.push_back(in1);

               vtkSpoke* in2 = new vtkSpoke;
               interpolater.Interpolate(1.0, double(steps[sj]), cornerSpokes, in2);
               interior.push_back(in2);
           }
       }
}

void vtkSlicerSkeletalRepresentationRefinerLogic::ReorderCrestSpokes(int nRows, int nCols,
                                                                     std::vector<vtkSpoke *> &input,
                                                                     std::vector<vtkSpoke *> &output)
{
    vtkSmartPointer<vtkPoints> foldCurvePts = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> foldCurveCell = vtkSmartPointer<vtkCellArray>::New();
    // top row
    for(int i = 0; i< nCols; ++i) {
        output.push_back(input[i]);
    }
    // right col
    for(int i = nCols + 1; i < (nRows-2)*2+nCols; i+=2) {
        output.push_back(input[i]);
    }
    // right bot
    //output.push_back(input[input.size()-1]);
    // bot row
    for(int i = input.size()-1; i >= (nRows-2)*2+nCols; --i) {
        output.push_back(input[i]);
    }

    // left col
    for(int i = (nRows-2)*2+nCols - 2; i >= nCols; i-=2) {
        output.push_back(input[i]);
    }
    //close circle
    output.push_back(input[0]);


    for(int i = 0; i < input.size()-1; ++i) {
        double pt0[3], pt1[3];
        input[i]->GetSkeletalPoint(pt0);
        input[i+1]->GetSkeletalPoint(pt1);
        vtkIdType id0 = foldCurvePts->InsertNextPoint(pt0);
        vtkIdType id1 = foldCurvePts->InsertNextPoint(pt1);

        vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
        line->GetPointIds()->SetId(0, id0);
        line->GetPointIds()->SetId(1, id1);
        foldCurveCell->InsertNextCell(line);
    }
    vtkSmartPointer<vtkPolyData> foldPoly = vtkSmartPointer<vtkPolyData>::New();
    foldPoly->SetPoints(foldCurvePts);
    foldPoly->SetPolys(foldCurveCell);
    Visualize(foldPoly, "fold curve", 1, 1, 0);
}

void vtkSlicerSkeletalRepresentationRefinerLogic::OptimizeCrestSpokeLength(vtkImplicitPolyDataDistance *distanceFunction,
                                                                           vtkSpoke *targetSpoke, double stepSize, int maxIter)
{
    // 1. Transform the boundary point to image cs. by applying [x, y, z, 1] * mTransformationMat
    double pt[3];
    targetSpoke->GetBoundaryPoint(pt);
    double epsilon = 1e-5;

    double dist = distanceFunction->FunctionValue(pt);
    double oldR = targetSpoke->GetRadius();
    double newR = oldR;
    int iter = 0;
    double oldDist = dist;
    // 2. iteratively update
    while(abs(dist) > epsilon)
    {
        double newBdry[3];
        if(dist > 0)
        {
            // if the spoke is too long, shorten it
            newR -= stepSize;
            targetSpoke->SetRadius(newR);
            targetSpoke->GetBoundaryPoint(newBdry);
            dist = distanceFunction->FunctionValue(newBdry);
            if(oldDist * dist < 0) {
                // if the spoke change from outside to inside, decay the learning rate
                stepSize /= 10;
                oldDist = dist;
            }
        }
        else {
            // enlongate the spoke
            newR += stepSize;
            targetSpoke->SetRadius(newR);
            targetSpoke->GetBoundaryPoint(newBdry);
            dist = distanceFunction->FunctionValue(newBdry);
            // if the spoke change from outside to inside, decay the learning rate
            if(oldDist * dist < 0) {
                stepSize /= 10;
                oldDist = dist;
            }
        }

        iter++;
        if(iter > maxIter)
        {
            break;
        }
    }
}

void vtkSlicerSkeletalRepresentationRefinerLogic::Transform2ImageCS(double *ptInput, int *ptOutput)
{
    ptInput[0] = ptInput[0] * mTransformationMat[0][0] + mTransformationMat[3][0];
    ptInput[1] = ptInput[1] * mTransformationMat[1][1] + mTransformationMat[3][1];
    ptInput[2] = ptInput[2] * mTransformationMat[2][2] + mTransformationMat[3][2];

    ptInput[0] /= voxelSpacing;
    ptInput[1] /= voxelSpacing;
    ptInput[2] /= voxelSpacing;

    int x = static_cast<int>(ptInput[0]+0.5);
    int y = static_cast<int>(ptInput[1]+0.5);
    int z = static_cast<int>(ptInput[2]+0.5);

    int maxX = static_cast<int>(1 / voxelSpacing - 1);
    int maxY = static_cast<int>(1 / voxelSpacing - 1);
    int maxZ = static_cast<int>(1 / voxelSpacing - 1);

    if(x > maxX) x = maxX;
    if(y > maxY) y = maxY;
    if(z > maxZ) z = maxZ;

    if(x < 0) x = 0;
    if(y < 0) y = 0;
    if(z < 0) z = 0;

    ptOutput[0] = x;
    ptOutput[1] = y;
    ptOutput[2] = z;
}

void vtkSlicerSkeletalRepresentationRefinerLogic::RetileMesh(vtkPolyData *targetMesh, vtkPolyData *impliedMesh,
                                                             vtkPolyData* retiledMesh)
{
    vtkSmartPointer<vtkPoints> newMeshPts = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkPolyData> newImpliedMesh = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkCellLocator> cellLocator = vtkSmartPointer<vtkCellLocator>::New();
    cellLocator->SetDataSet(impliedMesh);
    cellLocator->BuildLocator();
    for(int i = 0; i < targetMesh->GetNumberOfPoints(); ++i) {
        double pt[3], closestPt[3];
        targetMesh->GetPoint(i, pt);
        vtkIdType cellId;
        int subId;
        double d; // unsigned distances
        cellLocator->FindClosestPoint(pt, closestPt, cellId, subId, d);
        newMeshPts->InsertNextPoint(closestPt);
    }
    newImpliedMesh->SetPolys(targetMesh->GetPolys());
    newImpliedMesh->SetPoints(newMeshPts);
    newImpliedMesh->Modified();
    //cout << "before clean point num:" << newImpliedMesh->GetNumberOfPoints() << endl;
    vtkSmartPointer<vtkCleanPolyData> cleanFilter = vtkSmartPointer<vtkCleanPolyData>::New();
    cleanFilter->SetInputData(newImpliedMesh);
    cleanFilter->Update();
    retiledMesh->DeepCopy(newImpliedMesh);
    //cout << "after clean point num:" << retiledMesh->GetNumberOfPoints() << endl;

}
