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
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLDataParser.h>
#include "vtkSlicerSkeletalRepresentationInterpolater.h"
#include "vtkSrep.h"
#include "vtkSpoke.h"
#include "newuoa.h"
#include "vtkPolyData2ImageData.h"
#include "vtkImage2SignedDistanceMap.h"
#include "vtkAntiAlias.h"
#include "vtkApproximateSignedDistanceMap.h"

// STD includes
#include <cassert>

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
    mImageFilePath = imageFilePath;
}

void vtkSlicerSkeletalRepresentationRefinerLogic::SetSrepFileName(const std::string &srepFilePath)
{
    mSrepFilePath = srepFilePath;
}

void vtkSlicerSkeletalRepresentationRefinerLogic::Refine(double stepSize, double endCriterion, int maxIter)
{
    // 1. parse file
    const std::string srepFileName = "/playpen/ra_job/SlicerSkeletalRepresentation/SkeletalRepresentationVisualizer/Testing/test_data/hippo/up.vtp";
    const std::string headerFileName = "/playpen/ra_job/SlicerSkeletalRepresentation/SkeletalRepresentationVisualizer/Testing/test_data/hippo/header.xml";
    
    mCoeffArray.clear();
    std::vector<double> radii, dirs, skeletalPoints;
    Parse(srepFileName, mCoeffArray, radii, dirs, skeletalPoints);
    
    int nRows = 0, nCols = 0;
    ParseHeader(headerFileName, &nRows, &nCols);
    
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
    
    mNumRows = nRows; mNumCols = nCols;
    mSrep = srep;
    
    // make tuples of interpolation positions (u,v)
    mInterpolatePositions.clear();
    int interpolationLevel = 3;
    double interval = double(1.0 / (1+interpolationLevel));
    for(int i = 0; i < interpolationLevel + 1; ++i)
    {
        for(int j = 0; j < interpolationLevel + 1; ++j)
        {
            // no interpolation at corners
            if((i==0 && j == 0) || (i==0 && j==1) || (i==1 && j==0) || (i==1 && j == 1))
                continue;
            std::pair<double, double> uv = make_pair(i * interval, j * interval);
            mInterpolatePositions.insert(uv);
        }
    }
    
    // 2. Invoke newuoa to optimize
    min_newuoa(paramDim, coeff, *this, stepSize, endCriterion, maxIter);
    
    // 3. Visualize the refined srep
    srep->Refine(coeff);
    vtkSmartPointer<vtkPolyData> refinedSrep = vtkSmartPointer<vtkPolyData>::New();
    convertSpokes2PolyData(srep->GetAllSpokes(), refinedSrep);
    // Hide other nodes.
    HideNodesByClass("vtkMRMLModelNode");
    Visualize(refinedSrep, "Refined", 0, 1, 0);
}

void vtkSlicerSkeletalRepresentationRefinerLogic::InterpolateSrep(int interpolationLevel, std::string& srepFileName)
{
    // Hide other nodes.
    HideNodesByClass("vtkMRMLModelNode");
    srepFileName = "/playpen/ra_job/SlicerSkeletalRepresentation/SkeletalRepresentationVisualizer/Testing/test_data/hippo/mod_up.vtp";
    const std::string downFileName = "/playpen/ra_job/SlicerSkeletalRepresentation/SkeletalRepresentationVisualizer/Testing/test_data/hippo/down.vtp";
    const std::string crestFileName = "/playpen/ra_job/SlicerSkeletalRepresentation/SkeletalRepresentationVisualizer/Testing/test_data/hippo/crest.vtp";
    const std::string headerFileName = "/playpen/ra_job/SlicerSkeletalRepresentation/SkeletalRepresentationVisualizer/Testing/test_data/hippo/header.xml";
    // 1. Parse the model into a parameter array that needs to be optimized
    std::vector<double> coeffArrayUp, radiiUp, dirsUp, skeletalPointsUp;
    Parse(srepFileName, coeffArrayUp, radiiUp, dirsUp, skeletalPointsUp);
    
    int nRows = 0, nCols = 0;
    ParseHeader(headerFileName, &nRows, &nCols);
    
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
    //interpolationLevel = 3;
    double interval = 1.0/ (double)(interpolationLevel + 1);
    std::vector<double> steps;
    steps.push_back(0.0);
    for(int i = 1; i <= interpolationLevel; ++i)
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
                    
                    computeDerivative(skeletalPointsUp, r, c, nRows, nCols, dXdu11, dXdv11);
                    computeDerivative(skeletalPointsUp, r+1, c, nRows, nCols, dXdu21, dXdv21);
                    computeDerivative(skeletalPointsUp, r, c+1, nRows, nCols, dXdu12, dXdv12);
                    computeDerivative(skeletalPointsUp, r+1, c+1, nRows, nCols, dXdu22, dXdv22);
                    
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
    convertSpokes2PolyData(interpolatedSpokes, upSpokes_polyData);
    Visualize(upSpokes_polyData, "Interpolated", 1, 1, 1);
    
    vtkSmartPointer<vtkPolyData> primarySpokes = vtkSmartPointer<vtkPolyData>::New();
    convertSpokes2PolyData(srep->GetAllSpokes(), primarySpokes);
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
    if(mInterpolatePositions.empty())
    {
        std::cerr << "The interpolation pairs in the refinement are empty." << std::endl;
        return -100000.0;
    }
    if(mSrep == NULL)
    {
        std::cerr << "The srep pointer in the refinement is NULL." << std::endl;
        return -100000.0;
    }
    double imageDist = 0.0, normal = 0.0, srad = 0.0;
    int paramDim = mCoeffArray.size();
    int spokeNum = paramDim / 4;
    vtkSlicerSkeletalRepresentationInterpolater interpolater;
    
    // 1. Compute image match from all spokes and those spokes affected by them
    double totalDist = 0.0; // total distance square from implied boundary to real boundary of image
    for(int i = 0; i < spokeNum; ++i)
    {
        int r = i / mNumCols;
        int c = i % mNumCols;
        vtkSpoke *thisSpoke = mSrep->GetSpoke(r, c);
        
        // compute distance for this spoke
        
        vtkSpoke *cornerSpokes[4];
        for(auto it = mInterpolatePositions.begin(); it != mInterpolatePositions.end(); ++it)
        {
            double u = (*it).first;
            double v = (*it).second;
            
            // For each spoke at the corner of the srep,
            // its neighbors are all spokes in one quad
            if(r == 0 && c == 0)
            {
                cornerSpokes[0] = thisSpoke;
                cornerSpokes[0] = mSrep->GetSpoke(r+1, c);
                cornerSpokes[0] = mSrep->GetSpoke(r+1, c+1);
                cornerSpokes[0] = mSrep->GetSpoke(r, c+1);
                vtkSpoke *interpolatedSpoke;
                interpolater.Interpolate(u, v, cornerSpokes, interpolatedSpoke);
                
                // compute the ssd for this interpolated spoke
            }
            else if(r == 0 && c == mNumCols - 1)
            {
                
            }
            else if(r == mNumRows - 1 && c == 0)
            {
                
            }
            else if(r == mNumRows - 1 && c == mNumCols - 1)
            {
                
            }
            // For each spoke on the edge of the srep,
            // its neighbors are all spokes in two quads
            else if(r == 0)
            {
                
            }
            else if(r == mNumRows - 1)
            {
                
            }
            else if(c == 0)
            {
                
            }
            else if(c == mNumCols - 1)
            {
                
            }
            // for each spoke in the middle of the srep,
            // obtain image distance and normal from all interpolated spoke in 4 quads around it
            else
            {
                
            }
        }
        
    }
    
    if(mFirstCost)
    {
        // this log help to adjust the weights of three terms
        std::cout << "ImageMatch:" << imageDist << ", normal:" << normal << ", srad:" << srad << std::endl;
        mFirstCost = false;
    }
    return mWtImageMatch * imageDist + mWtNormalMatch * normal + mWtSrad * srad;
}

void vtkSlicerSkeletalRepresentationRefinerLogic::AntiAliasSignedDistanceMap(const std::string &meshFileName)
{
    // 1. convert poly data to image data
    vtkPolyData2ImageData polyDataConverter;
    vtkSmartPointer<vtkImageData> img = vtkSmartPointer<vtkImageData>::New();
            
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
    ssdGenerator.Convert(img, antiAliasedImage);
    // 4. output / save the result
    mAntiAliasedImage->DeepCopy(antiAliasedImage);
}

void vtkSlicerSkeletalRepresentationRefinerLogic::computeDerivative(std::vector<double> skeletalPoints, int r, int c, int nRows, int nCols, double *dXdu, double *dXdv)
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
        head[0] = skeletalPoints[(id-nCols)*3];
        head[1] = skeletalPoints[(id-nCols)*3+1];
        head[2] = skeletalPoints[(id-nCols)*3+2];
        
        tail[0] = skeletalPoints[(id)*3];
        tail[1] = skeletalPoints[(id)*3+1];
        tail[2] = skeletalPoints[(id)*3+2];
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
    computeDiff(head, tail, factor, dXdu);
    
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
        head[0] = skeletalPoints[(id-1)*3];
        head[1] = skeletalPoints[(id-1)*3+1];
        head[2] = skeletalPoints[(id-1)*3+2];
        
        tail[0] = skeletalPoints[(id)*3];
        tail[1] = skeletalPoints[(id)*3+1];
        tail[2] = skeletalPoints[(id)*3+2];
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
    computeDiff(head, tail, factor, dXdv);
    
}

void vtkSlicerSkeletalRepresentationRefinerLogic::convertSpokes2PolyData(std::vector<vtkSpoke *>input, vtkPolyData *output)
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

void vtkSlicerSkeletalRepresentationRefinerLogic::ParseHeader(const std::string &headerFileName, int *nRows, int *nCols)
{
    vtkSmartPointer<vtkXMLDataParser> parser = vtkSmartPointer<vtkXMLDataParser>::New();
    
    parser->SetFileName(headerFileName.c_str());
    parser->SetIgnoreCharacterData(0);
    
    if( parser->Parse() == 1)
    {
        vtkXMLDataElement *root = parser->GetRootElement(); 
        int numElements = root->GetNumberOfNestedElements();;
        if(numElements > 1)
        {
            vtkXMLDataElement *eRow = root->GetNestedElement(0);
            char *pEnd;
            int r, c;
            r = strtol(eRow->GetCharacterData(), &pEnd, 10);
            *nRows = r;
            
            vtkXMLDataElement *eCol = root->GetNestedElement(1);
            c = strtol(eCol->GetCharacterData(), &pEnd, 10);
            *nCols = c;
            
        }
    }
}

void vtkSlicerSkeletalRepresentationRefinerLogic::computeDiff(double *head, double *tail, double factor, double *output)
{
    output[0] = factor * (head[0] - tail[0]);
    output[1] = factor * (head[1] - tail[1]);
    output[2] = factor * (head[2] - tail[2]);
}

void vtkSlicerSkeletalRepresentationRefinerLogic::computeDistance(vtkSpoke *theSpoke)
{
    
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
