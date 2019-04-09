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
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkImageData.h>
#include <vtksys/SystemTools.hxx>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLDataParser.h>
#include "vtkSlicerSkeletalRepresentationInterpolater.h"
#include "vtkSrep.h"
#include "vtkSpoke.h"
#include "newuoa.h"
#include "vtkPolyData2ImageData.h"
#include "vtkApproximateSignedDistanceMap.h"
#include "vtkGradientDistanceFilter.h"

// STD includes
#include <cassert>
const double voxelSpacing = 0.005;
const std::string newFilePrefix = "/refined_";
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
}

void vtkSlicerSkeletalRepresentationRefinerLogic::SetSrepFileName(const std::string &srepFilePath)
{
    mSrepFilePath = srepFilePath;
}

void vtkSlicerSkeletalRepresentationRefinerLogic::SetOutputPath(const string &outputPath)
{
    mOutputPath = outputPath;
}

int iterNum = 0;
void vtkSlicerSkeletalRepresentationRefinerLogic::Refine(double stepSize, double endCriterion, int maxIter, int interpolationLevel)
{
    mFirstCost = true;
    iterNum = 0;
    // 1. parse file
    const std::string headerFileName = mSrepFilePath;
    int nRows = 0, nCols = 0;
    std::string up, down, crest;
    ParseHeader(headerFileName, &nRows, &nCols, &up, &down, &crest);
    
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
    mInterpolatePositions.clear();
    double tol = 1e-6;
    int shares = pow(2, interpolationLevel);
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
    
    // Refine up spokes
    RefinePartOfSpokes(up, stepSize, endCriterion, maxIter);
    
    // Refine down spokes
    RefinePartOfSpokes(down, stepSize, endCriterion, maxIter);
    
    // Show crest spokes
    
    // Update header file
    UpdateHeader(headerFileName, mOutputPath);
    
}

void vtkSlicerSkeletalRepresentationRefinerLogic::InterpolateSrep(int interpolationLevel, std::string& srepFileName)
{
    // Hide other nodes.
    HideNodesByClass("vtkMRMLModelNode");
    
    // 1. Parse the model into a parameter array that needs to be optimized
    int nRows = 0, nCols = 0;
    std::string up, down, crest;
    ParseHeader(srepFileName, &nRows, &nCols, &up, &down, &crest);
    
    std::vector<double> coeffArrayUp, radiiUp, dirsUp, skeletalPointsUp;
    Parse(up, coeffArrayUp, radiiUp, dirsUp, skeletalPointsUp);
        
    if(nRows == 0 || nCols == 0)
    {
        std::cerr << "The s-rep model is empty." << std::endl;
        return;
    }
    
    vtkSrep *srep = new vtkSrep(nRows, nCols, radiiUp, dirsUp, skeletalPointsUp);
    if(srep->IsEmpty())
    {
        std::cerr << "The s-rep model is empty." << std::endl;
        delete srep;
        srep = NULL;
        return;
    }
    // 1.1 interpolate and visualize for verification
    // collect neighboring spokes around corners
    vtkSlicerSkeletalRepresentationInterpolater interpolater;
    
    int shares = pow(2, interpolationLevel);
    double interval = (double)(1.0/ shares);
    std::vector<double> steps;
    steps.push_back(0.0);
    for(int i = 1; i <= shares; ++i)
    {
        steps.push_back(i * interval);
    }
    steps.push_back(1.0);
    
    std::vector<vtkSpoke*> interpolatedSpokes;
    for(int r = 0; r < nRows-1; ++r)
    {
        for(int c = 0; c < nCols-1; ++c)
        {
            vtkSpoke *cornerSpokes[4];
            
            double  dXdu11[3], dXdv11[3], 
                    dXdu12[3], dXdv12[3], 
                    dXdu21[3], dXdv21[3], 
                    dXdu22[3], dXdv22[3];
            
            for(int i = 0; i < steps.size(); ++i)
            {
                for(int j = 0; j < steps.size(); ++j)
                {
                    cornerSpokes[0] = srep->GetSpoke(r,c);
                    cornerSpokes[1] = srep->GetSpoke(r+1, c);
                    cornerSpokes[2] = srep->GetSpoke(r+1, c+1);
                    cornerSpokes[3] = srep->GetSpoke(r, c+ 1);
                    
                    ComputeDerivative(skeletalPointsUp, r, c, nRows, nCols, dXdu11, dXdv11);
                    ComputeDerivative(skeletalPointsUp, r+1, c, nRows, nCols, dXdu21, dXdv21);
                    ComputeDerivative(skeletalPointsUp, r, c+1, nRows, nCols, dXdu12, dXdv12);
                    ComputeDerivative(skeletalPointsUp, r+1, c+1, nRows, nCols, dXdu22, dXdv22);
                    
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
                    interpolatedSpokes.push_back(in1);
                       
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
    
    if(mSrep == NULL)
    {
        std::cerr << "The srep pointer in the refinement is NULL." << std::endl;
        return -100000.0;
    }
    
    // this temporary srep is constructed to compute the cost function value
    // The original srep should not be changed by each iteration
    vtkSrep *tempSrep = new vtkSrep();
    tempSrep->DeepCopy(*mSrep);
    tempSrep->Refine(coeff);
    double imageDist = 0.0, normal = 0.0, srad = 0.0;
    int paramDim = mCoeffArray.size();
    int spokeNum = paramDim / 4;
    // 1. Compute image match from all spokes and those spokes affected by them
    for(int i = 0; i < spokeNum; ++i)
    {
        int r = i / mNumCols;
        int c = i % mNumCols;
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
            else if(r == 0 && c == mNumCols - 1)
            {
                // right-top corner
                imageDist += TotalDistOfRightTopSpoke(tempSrep, u, v, r,c, &normal);
            }
            else if(r == mNumRows - 1 && c == 0)
            {
                // left-bot corner
                imageDist += TotalDistOfLeftBotSpoke(tempSrep, u, v, r,c, &normal);
            }
            else if(r == mNumRows - 1 && c == mNumCols - 1)
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
            else if(r == mNumRows - 1)
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
            else if(c == mNumCols - 1)
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
    
    // 2. convert image data to signed distance map
//    vtkImage2SignedDistanceMap imageConverter;
//    vtkSmartPointer<vtkImageData> signedDistanceImage = vtkSmartPointer<vtkImageData>::New();
//    imageConverter.Convert(img, signedDistanceImage);
    
//    // 3. anti-aliasing
//    vtkAntiAlias aa;
//    vtkSmartPointer<vtkImageData> antiAliasedImage = vtkSmartPointer<vtkImageData>::New();
//    aa.Filter(signedDistanceImage, antiAliasedImage);
    
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
    ParseHeader(headerFile, &nRows, &nCols, &up, &down, &crest);
    
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
        srep = NULL;
        return;
    }
    
    std::vector<double> radiiDown, dirsDown, skeletalPointsDown, coeffDown;
    Parse(down, coeffDown, radiiDown, dirsDown, skeletalPointsDown);
    srep->AddSpokes(radiiDown, dirsDown, skeletalPointsDown);
    
    std::vector<double> radiiCrest, dirsCrest, skeletalPointsCrest, coeffCrest;
    Parse(crest, coeffCrest, radiiCrest, dirsCrest, skeletalPointsCrest);
    srep->AddSpokes(radiiCrest, dirsCrest, skeletalPointsCrest);
    
    TransformSrep2ImageCS(srep, mTransformationMat);
    
    vtkSmartPointer<vtkPolyData> primarySpokes = vtkSmartPointer<vtkPolyData>::New();
    ConvertSpokes2PolyData(srep->GetAllSpokes(), primarySpokes);
    Visualize(primarySpokes, "Primary", 1, 0, 0);
    
    vtkSmartPointer<vtkPolyData> transUpPrimary = vtkSmartPointer<vtkPolyData>::New();
    TransSpokes2PolyData(srep->GetAllSpokes(), transUpPrimary);
    Visualize(transUpPrimary, "Trans_Primary", 1, 0, 0);
    delete srep;
}

void vtkSlicerSkeletalRepresentationRefinerLogic::ComputeDerivative(std::vector<double> skeletalPoints, int r, int c, int nRows, int nCols, double *dXdu, double *dXdv)
{
    // 0-based index of elements if arranged in array
    int id = r * nCols + c;
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
    else if(r == nRows - 1)
    {
        // last row
        // backward difference
        tail[0] = skeletalPoints[(id-nCols)*3];
        tail[1] = skeletalPoints[(id-nCols)*3+1];
        tail[2] = skeletalPoints[(id-nCols)*3+2];
        
        head[0] = skeletalPoints[(id)*3];
        head[1] = skeletalPoints[(id)*3+1];
        head[2] = skeletalPoints[(id)*3+2];
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
  
    for(int i = 0; i < input.size(); ++i)
    {
        vtkSpoke* currSpoke = input[i];
        double basePt[3], bdryPt[3];
        currSpoke->GetSkeletalPoint(basePt);
        currSpoke->GetBoundaryPoint(bdryPt);
        int id0 = (pts->InsertNextPoint(basePt[0], basePt[1], basePt[2]));
        int id1 = pts->InsertNextPoint(bdryPt[0], bdryPt[1], bdryPt[2]);
        
        vtkSmartPointer<vtkLine> currLine = vtkSmartPointer<vtkLine>::New();
        currLine->GetPointIds()->SetId(0, id0);
        currLine->GetPointIds()->SetId(1, id1);
        arrows->InsertNextCell(currLine);
    }
    output->SetPoints(pts);
    output->SetLines(arrows);
}

void vtkSlicerSkeletalRepresentationRefinerLogic::TransSpokes2PolyData(std::vector<vtkSpoke *>input, vtkPolyData *output)
{
    vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> arrows = vtkSmartPointer<vtkCellArray>::New();
  
    for(int i = 0; i < input.size(); ++i)
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
        
        int id0 = pts->InsertNextPoint(basePt[0], basePt[1], basePt[2]);
        int id1 = pts->InsertNextPoint(bdryPt[0], bdryPt[1], bdryPt[2]);
        
        vtkSmartPointer<vtkLine> currLine = vtkSmartPointer<vtkLine>::New();
        currLine->GetPointIds()->SetId(0, id0);
        currLine->GetPointIds()->SetId(1, id1);
        arrows->InsertNextCell(currLine);
    }
    output->SetPoints(pts);
    output->SetLines(arrows);
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
  assert(this->GetMRMLScene() != 0);
}

//---------------------------------------------------------------------------
void vtkSlicerSkeletalRepresentationRefinerLogic::UpdateFromMRMLScene()
{
  assert(this->GetMRMLScene() != 0);
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
    vtkSmartPointer<vtkPointData> spokesPointData = spokesPolyData->GetPointData();
    int numOfArrays = spokesPointData->GetNumberOfArrays();
    int numOfSpokes = spokesPolyData->GetNumberOfPoints();

    if(numOfSpokes == 0 || numOfArrays == 0)
    {
        return;
    }

    // including Ux, Uy, Uz, r
    vtkSmartPointer<vtkDoubleArray> spokeRadii = vtkDoubleArray::SafeDownCast(spokesPointData->GetArray("spokeLength"));
    vtkSmartPointer<vtkDoubleArray> spokeDirs = vtkDoubleArray::SafeDownCast(spokesPointData->GetArray("spokeDirection"));
    
    for(int i = 0; i < numOfSpokes; ++i)
    {
        int idx = i * 4; // Ux, Uy, Uz,r
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
        spokesPolyData->GetPoint(i, tempSkeletalPoint);
        skeletalPoints.push_back(tempSkeletalPoint[0]);
        skeletalPoints.push_back(tempSkeletalPoint[1]);
        skeletalPoints.push_back(tempSkeletalPoint[2]);
    }
}

void vtkSlicerSkeletalRepresentationRefinerLogic::ParseHeader(const std::string &headerFileName, int *nRows, int *nCols,
                                                               std::string* upFileName, std::string* downFileName, std::string* crestFileName)
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
            char *pEnd; 
            vtkXMLDataElement *e = root->GetNestedElement(i);
            std::string estimatePath;
            estimatePath = vtksys::SystemTools::GetFilenamePath(headerFileName) + "/";
            std::vector<std::string> components;
            components.push_back(estimatePath);
            
            char* eName = e->GetName();
            if(strcmp(eName, "nRows") == 0)
            {
                r = strtol(e->GetCharacterData(), &pEnd, 10);
                *nRows = r;
            }
            else if(strcmp(eName, "nCols") == 0)
            {
                c = strtol(e->GetCharacterData(), &pEnd, 10);
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
            }
            else if(strcmp(eName, "downSpoke")==0)
            {
                *downFileName = e->GetCharacterData();
                if(!vtksys::SystemTools::FileIsFullPath(*downFileName))
                {
                    components.push_back(*downFileName);
                    *downFileName = vtksys::SystemTools::JoinPath(components);
                }
            }
            else if(strcmp(eName, "crestSpoke") == 0)
            {
                *crestFileName = e->GetCharacterData();
                if(!vtksys::SystemTools::FileIsFullPath(*crestFileName))
                {
                    components.push_back(*crestFileName);
                    *crestFileName = vtksys::SystemTools::JoinPath(components);
                }
            }
            
        }
    }
    
}

void vtkSlicerSkeletalRepresentationRefinerLogic::UpdateHeader(const string &headerFileName, const string &outputFilePath)
{
    vtkSmartPointer<vtkXMLDataParser> parser = vtkSmartPointer<vtkXMLDataParser>::New();
    
    parser->SetFileName(headerFileName.c_str());
    parser->SetIgnoreCharacterData(0);
    if( parser->Parse() == 1)
    {
        vtkXMLDataElement *root = parser->GetRootElement(); 
        int numElements = root->GetNumberOfNestedElements();;
        std::string newUpFileName, newDownFileName, newCrestFileName;
        
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
                r = strtol(e->GetCharacterData(), &pEnd, 10);
                nRows = r;
            }
            else if(strcmp(eName, "nCols") == 0)
            {
                c = strtol(e->GetCharacterData(), &pEnd, 10);
                nCols = c;
            }
            else if(strcmp(eName, "upSpoke") == 0)
            {
                std::string oldFile(e->GetCharacterData());
                // some file paths are relative path, others are absolute path
                oldFile = vtksys::SystemTools::GetFilenameName(oldFile);
                newUpFileName = outputFilePath + newFilePrefix + oldFile;
            }
            else if(strcmp(eName, "downSpoke")==0)
            {
                std::string oldFile(e->GetCharacterData());
                // some file paths are relative path, others are absolute path
                oldFile = vtksys::SystemTools::GetFilenameName(oldFile);
                newDownFileName = outputFilePath + newFilePrefix + oldFile;
            }
            else if(strcmp(eName, "crestSpoke")==0)
            {
                std::string oldFile(e->GetCharacterData());
                // some file paths are relative path, others are absolute path
                oldFile = vtksys::SystemTools::GetFilenameName(oldFile);
                newCrestFileName = outputFilePath + oldFile;
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
        oldHeader = outputFilePath + newFilePrefix + oldHeader;
        std::string header_file(oldHeader);
        std::ofstream out_file;
        out_file.open(header_file);
        out_file << output.rdbuf();
        out_file.close();
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
    
    int maxX = 1 / voxelSpacing - 1;
    int maxY = 1 / voxelSpacing - 1;
    int maxZ = 1 / voxelSpacing - 1;
    
    if(x > maxX) x = maxX;
    if(y > maxY) y = maxY;
    if(z > maxZ) z = maxZ;
    
    if(x < 0) x = 0;
    if(y < 0) y = 0;
    if(z < 0) z = 0;
    
    if(mAntiAliasedImage == NULL)
    {
        std::cerr << "The image in this RefinerLogic instance is empty." << std::endl;
        return -10000.0;
    }
    RealImage::IndexType pixelIndex = {{x,y,z}};
    float dist = mAntiAliasedImage->GetPixel(pixelIndex);
    
    if(mGradDistImage == NULL)
    {
        return -10000.0;
    }
    
    VectorImage::IndexType indexGrad;
    indexGrad[0] = x;
    indexGrad[1] = y;
    indexGrad[2] = z;
    
    VectorImage::PixelType grad = mGradDistImage->GetPixel(indexGrad);
    double normalVector[3];
    normalVector[0] = grad[0]; normalVector[1] = grad[1]; normalVector[2] = grad[2];
    // normalize the normal vector
    vtkMath::Normalize(normalVector);
    
    double spokeDir[3]; 
    theSpoke->GetDirection(spokeDir);
    double dotProduct = vtkMath::Dot(normalVector, spokeDir);
    double distSqr = dist * dist;
    
    // The normal match (between [0,1]) is scaled by the distance so that the overall term is comparable
    *normalMatch = *normalMatch + distSqr * (1 - dotProduct);
    // return square of distance
    return distSqr;
}

void vtkSlicerSkeletalRepresentationRefinerLogic::Visualize(vtkPolyData *model, const std::string &modelName, double r, double g, double b)
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
    if(displayModelNode == NULL)
    {
        vtkErrorMacro("displayModelNode is NULL");
        return;
    }
    displayModelNode->SetColor(r,g,b);
    displayModelNode->SetScene(scene);
    displayModelNode->SetLineWidth(2.0);
    displayModelNode->SetBackfaceCulling(0);
    displayModelNode->SetRepresentation(1);
    
    if(true)
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
        if(displayNode == NULL)
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
        mat4x4 = NULL;
        return;
    }
    // 1. Find the bounding box of boundary
    std::vector<vtkSpoke *> spokes = input->GetAllSpokes();
    vtkSmartPointer<vtkPoints> boundaryPts = 
            vtkSmartPointer<vtkPoints>::New();
    for(int i = 0; i < spokes.size(); ++i)
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

void vtkSlicerSkeletalRepresentationRefinerLogic::RefinePartOfSpokes(const string &srepFileName, double stepSize, double endCriterion, int maxIter)
{
    mCoeffArray.clear();
    std::vector<double> radii, dirs, skeletalPoints;
    Parse(srepFileName, mCoeffArray, radii, dirs, skeletalPoints);
    
    vtkSrep *srep = new vtkSrep(mNumRows, mNumCols, radii, dirs, skeletalPoints);
    if(srep->IsEmpty())
    {
        std::cerr << "The s-rep model is empty." << std::endl;
        delete srep;
        srep = NULL;
        return;
    }
    
    // total number of parameters that need to optimize
    int paramDim = mCoeffArray.size();
    double coeff[paramDim];
    for(int i = 0; i < paramDim; ++i)
    {
        coeff[i] = mCoeffArray[i];
    }
    
    mSrep = srep;
    vtkSmartPointer<vtkPolyData> origSrep = vtkSmartPointer<vtkPolyData>::New();
    ConvertSpokes2PolyData(srep->GetAllSpokes(), origSrep);
    // Hide other nodes.
    HideNodesByClass("vtkMRMLModelNode");
    Visualize(origSrep, "Before refinement", 1, 0, 0);
 
    // 2. Invoke newuoa to optimize
    min_newuoa(paramDim, coeff, *this, stepSize, endCriterion, maxIter);
    
    // Re-evaluate the cost
    mFirstCost = true;
    EvaluateObjectiveFunction(coeff);
    
    // 3. Visualize the refined srep
    srep->Refine(coeff);
    vtkSmartPointer<vtkPolyData> refinedSrep = vtkSmartPointer<vtkPolyData>::New();
    ConvertSpokes2PolyData(srep->GetAllSpokes(), refinedSrep);
    // Hide other nodes.
    //HideNodesByClass("vtkMRMLModelNode");
    Visualize(refinedSrep, "Refined", 0, 1, 0);
    
    // write to vtp file
    std::string outputFile(mOutputPath);
    std::string fileName = vtksys::SystemTools::GetFilenameName(srepFileName);
    outputFile = outputFile + newFilePrefix + fileName;
    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(outputFile.c_str());
    writer->SetInputData(refinedSrep);
    writer->Update();
    if(mSrep != NULL)
    {
        delete mSrep;
        mSrep = NULL;
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
    
    std::vector<vtkSpoke*> interpolatedSpokes;
    std::vector<vtkSpoke*> primarySpokes;
    for(int r = 0; r < nRows; ++r)
    {
        for(int c = 0; c < nCols; ++c)
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
            else if(r == 0 && c == nCols - 1)
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
            else if(r == nRows - 1 && c == 0)
            {
                // left bot corner
                FindBotLeftNeigbors(r, c, input, neighborU, neighborV);
                thisSpoke->SetNeighborU(neighborU, false);
                thisSpoke->SetNeighborV(neighborV, true);
            }
            else if(r == nRows - 1 && c == nCols -1)
            {
                // right bot corner
                FindBotRightNeigbors(r, c, input, neighborU, neighborV);
                thisSpoke->SetNeighborU(neighborU, false);
                thisSpoke->SetNeighborV(neighborV, false);
            }
            else if(r == nRows - 1)
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
            else if(c == nCols - 1)
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
    for(int i = 0; i < primarySpokes.size(); ++i)
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
