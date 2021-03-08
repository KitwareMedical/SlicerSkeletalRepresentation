#include "itkImageFileWriter.h"

#include "itkSmoothingRecursiveGaussianImageFilter.h"

#include "itkPluginUtilities.h"

#include "itkThinPlateSplineKernelTransform.h"
#include "itkResampleImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkBSplineTransform.h"
#include "itkLandmarkBasedTransformInitializer.h"
#include "itkVector.h"
#include "itkLandmarkDisplacementFieldSource.h"
#include "itkWarpImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkMesh.h"
#include "itkTriangleMeshToBinaryImageFilter.h"
#include "itkDisplacementFieldTransform.h"

#include "vtkAppendPolyData.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkXMLDataElement.h"
#include "vtkXMLDataParser.h"
#include "vtkXMLPolyDataReader.h"
#include "vtkXMLPolyDataWriter.h"

#include "vtksys/SystemTools.hxx"

#include <cstdlib>
#include <fstream>
#include <iostream>

#include "InitSRepFromSPHARMCLP.h"

// Use an anonymous namespace to keep class types and function names
// from colliding when module is used as shared object module.  Every
// thing should be in an anonymous namespace except for the module
// entry point, e.g. main()
//
namespace
{
  const unsigned int Dimension = 3;
  using ImageType = itk::Image<float, Dimension>;
  using FixedPointSetType = itk::PointSet<float, Dimension>;
  using MeshType = itk::Mesh<float, Dimension>;
  using DisplacementFieldTransformType = itk::DisplacementFieldTransform<float, 3>;

  ImageType::Pointer meshToImage(MeshType::Pointer mesh)
  {
    auto bounds = mesh->GetBoundingBox()->GetBounds();

    ImageType::SpacingType spacing;
    spacing[0] = (bounds[1] - bounds[0]) / 90;
    spacing[1] = (bounds[3] - bounds[2]) / 90;
    spacing[2] = (bounds[5] - bounds[4]) / 90;

    ImageType::PointType origin;
    origin[0] = bounds[0] - 5*spacing[0];
    origin[1] = bounds[2] - 5*spacing[1];
    origin[2] = bounds[4] - 5*spacing[2];
    
    ImageType::SizeType size;
    size[0] = 100;
    size[1] = 100;
    size[2] = 100;

    using MeshToImageType = itk::TriangleMeshToBinaryImageFilter<MeshType,ImageType>;
    auto meshToImage = MeshToImageType::New();
    meshToImage->SetInput(mesh);
    meshToImage->SetOrigin(origin);
    meshToImage->SetSpacing(spacing);
    meshToImage->SetSize(size);
    meshToImage->Update();

    auto image = meshToImage->GetOutput();
    return image;
  }

  MeshType::Pointer polyDataToMesh(vtkSmartPointer<vtkPolyData> pd)
  {
    auto mesh = MeshType::New();

    int numberOfPoints = pd->GetNumberOfPoints();
    mesh->GetPoints()->Reserve( numberOfPoints );

    for (int p = 0; p < numberOfPoints; p++)
    {
      double* point = pd->GetPoint(p);
      mesh->SetPoint( p, MeshType::PointType( point ) );
    }

    int numberOfTriangles = pd->GetNumberOfCells();
    mesh->GetCells()->Reserve( numberOfTriangles );

    typedef MeshType::CellType CellType;
    typedef itk::TriangleCell<CellType> TriangleCellType;

    for (int c = 0; c < numberOfTriangles; c++)
    {
      unsigned long pointIds[3];
      pointIds[0] = pd->GetCell(c)->GetPointIds()->GetId(0);
      pointIds[1] = pd->GetCell(c)->GetPointIds()->GetId(1);
      pointIds[2] = pd->GetCell(c)->GetPointIds()->GetId(2);
      MeshType::CellAutoPointer itkCell;
      TriangleCellType *tcell = new TriangleCellType;
      TriangleCellType::PointIdentifier itkPts[3];
      for (int ii = 0; ii < 3; ++ii)
      {
        itkPts[ii] = static_cast<TriangleCellType::PointIdentifier>(pointIds[ii]);
      }
      tcell->SetPointIds(itkPts);
      itkCell.TakeOwnership(tcell);
      mesh->SetCell(c,itkCell);
    }  

    return mesh;
  }

  vtkSmartPointer<vtkPolyData> warpSrepPart(vtkSmartPointer<vtkPolyData> pd, DisplacementFieldTransformType::Pointer transform)
  {
      auto num_arrays = pd->GetPointData()->GetNumberOfArrays();
      auto spoke_len_array = pd->GetPointData()->GetArray("spokeLength");
      auto spoke_dir_array = pd->GetPointData()->GetArray("spokeDirection");

      auto warped_pd = vtkSmartPointer<vtkPolyData>::New();
      warped_pd->DeepCopy(pd);

      for (int i = 0; i < pd->GetNumberOfPoints(); i++)
      {
          auto pt = pd->GetPoint(i);
          auto len = spoke_len_array->GetComponent(i, 0);
          auto dir = spoke_dir_array->GetTuple(i);

          std::cout << i << std::endl;
          std::cout << "pt: " << pt[0] << "," << pt[1] << "," << pt[2] << std::endl;
          std::cout << "len: " << len << std::endl;
          std::cout << "dir: " << dir[0] << "," << dir[1] << "," << dir[2] << std::endl;

          auto bdry_pt = new double[3];
          bdry_pt[0] = pt[0] + len * dir[0];
          bdry_pt[1] = pt[1] + len * dir[1];
          bdry_pt[2] = pt[2] + len * dir[2];

          std::cout << "bdry_pt: " << bdry_pt[0] << "," << bdry_pt[1] << "," << bdry_pt[2] << std::endl;

          // Warp both points
          auto warped_pt = transform->TransformPoint(pt);
          auto warped_bdry_pt = transform->TransformPoint(bdry_pt);

          auto warped_dir = new double[3];
          warped_dir[0] = warped_bdry_pt[0] - warped_pt[0];
          warped_dir[1] = warped_bdry_pt[1] - warped_pt[1];
          warped_dir[2] = warped_bdry_pt[2] - warped_pt[2];

          auto warped_len = sqrt(warped_dir[0] * warped_dir[0] + warped_dir[1] * warped_dir[1] + warped_dir[2] * warped_dir[2]);

          warped_dir[0] /= warped_len;
          warped_dir[1] /= warped_len;
          warped_dir[2] /= warped_len;

          warped_pd->GetPoints()->SetPoint(i, warped_pt.GetDataPointer());
          warped_pd->GetPointData()->GetArray("spokeLength")->SetComponent(i, 0, warped_len);
          warped_pd->GetPointData()->GetArray("spokeDirection")->SetTuple(i, warped_dir);
      }

      return warped_pd;
  }

  int DoIt(int argc, char * argv[])
  {
    PARSE_ARGS;

    // Read in template and target meshes

    auto template_reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
    template_reader->SetFileName(ellipsoid.c_str());
    template_reader->Update();
    auto template_mesh = template_reader->GetOutput();

    auto target_reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
    target_reader->SetFileName(target.c_str());
    target_reader->Update();
    auto target_mesh = target_reader->GetOutput();

    const unsigned int Dimension = 3;

    // Combine meshes and set up bounding box image

    auto append = vtkSmartPointer<vtkAppendPolyData>::New();
    append->AddInputData(template_mesh);
    append->AddInputData(target_mesh);
    append->Update();

    auto merged_mesh = polyDataToMesh(append->GetOutput());
    auto bounding_image = meshToImage(merged_mesh);

    // Set up interpolation 

    typedef   itk::Vector< float, Dimension >             FieldVectorType;
    typedef   itk::Image< FieldVectorType, Dimension >   DisplacementFieldType;

    typedef itk::LandmarkDisplacementFieldSource<DisplacementFieldType> DisplacementSourceType;
    DisplacementSourceType::Pointer deformer = DisplacementSourceType::New();
    deformer->SetOutputSpacing(bounding_image->GetSpacing());
    deformer->SetOutputOrigin(bounding_image->GetOrigin());
    deformer->SetOutputRegion(bounding_image->GetLargestPossibleRegion());
    deformer->SetOutputDirection(bounding_image->GetDirection());

    std::cout << deformer->GetKernelTransform()->GetNameOfClass() << std::endl;

    DisplacementSourceType::KernelTransformType::Pointer kernel = deformer->GetKernelTransform();
    kernel->SetStiffness(0.1);
    std::cout << "Stiffness: " << kernel->GetStiffness() << std::endl;

    typedef DisplacementSourceType::LandmarkContainer ContainerType;
    typedef DisplacementSourceType::LandmarkPointType PointType;

    ContainerType::Pointer sourcePointsContainer = ContainerType::New();
    ContainerType::Pointer targetPointsContainer = ContainerType::New();

    for (size_t i = 0; i < template_mesh->GetNumberOfPoints(); i++)
    {
      auto source_point = template_mesh->GetPoint(i);
      PointType p;
      p[0] = -1 * source_point[0];
      p[1] = -1 * source_point[1];
      p[2] = source_point[2];
      sourcePointsContainer->push_back(p);
    }

    for (size_t i = 0; i < target_mesh->GetNumberOfPoints(); i++)
    {
      auto target_point = target_mesh->GetPoint(i);
      PointType p;
      p[0] = -1 * target_point[0];// + (rand() / double(RAND_MAX) ) - 0.5;
      p[1] = -1 * target_point[1];// + (rand() / double(RAND_MAX) ) - 0.5;
      p[2] = target_point[2];// + (rand() / double(RAND_MAX) ) - 0.5;
      targetPointsContainer->push_back(p);
    }

    double minx, miny, minz;
    double maxx, maxy, maxz;

    minx = -1 * template_mesh->GetPoint(0)[0];
    miny = -1 * template_mesh->GetPoint(0)[1];
    minz = template_mesh->GetPoint(0)[2];

    maxx = -1 * template_mesh->GetPoint(0)[0];
    maxy = -1 * template_mesh->GetPoint(0)[1];
    maxz = template_mesh->GetPoint(0)[2];

    for (int i = 1; i < template_mesh->GetNumberOfPoints(); i++)
    {
      if (-1 * template_mesh->GetPoint(i)[0] < minx)
        minx = -1 * template_mesh->GetPoint(i)[0];

      if (-1 * template_mesh->GetPoint(i)[1] < miny)
        miny = -1 * template_mesh->GetPoint(i)[1];

      if (template_mesh->GetPoint(i)[2] < minz)
        minz = template_mesh->GetPoint(i)[2];

      if (-1 * template_mesh->GetPoint(i)[0] > maxx)
        maxx = -1 * template_mesh->GetPoint(i)[0];

      if (-1 * template_mesh->GetPoint(i)[1] > maxy)
        maxy = -1 * template_mesh->GetPoint(i)[1];

      if (template_mesh->GetPoint(i)[2] > maxz)
        maxz = template_mesh->GetPoint(i)[2];
    }

    for (int i = 1; i < template_mesh->GetNumberOfPoints(); i++)
    {
      if (-1 * target_mesh->GetPoint(i)[0] < minx)
        minx = -1 * target_mesh->GetPoint(i)[0];

      if (-1 * target_mesh->GetPoint(i)[1] < miny)
        miny = -1 * target_mesh->GetPoint(i)[1];

      if (target_mesh->GetPoint(i)[2] < minz)
        minz = target_mesh->GetPoint(i)[2];

      if (-1 * target_mesh->GetPoint(i)[0] > maxx)
        maxx = -1 * target_mesh->GetPoint(i)[0];

      if (-1 * target_mesh->GetPoint(i)[1] > maxy)
        maxy = -1 * target_mesh->GetPoint(i)[1];

      if (target_mesh->GetPoint(i)[2] > maxz)
        maxz = target_mesh->GetPoint(i)[2];
    }

    typename ImageType::IndexType minIndex, maxIndex;
    PointType minPoint, maxPoint;

    minPoint[0] = minx - 5 * bounding_image->GetSpacing()[0];
    minPoint[1] = miny - 5 * bounding_image->GetSpacing()[1];
    minPoint[2] = minz - 5 * bounding_image->GetSpacing()[2];

    maxPoint[0] = maxx + 5 * bounding_image->GetSpacing()[0];
    maxPoint[1] = maxy + 5 * bounding_image->GetSpacing()[1];
    maxPoint[2] = maxz + 5 * bounding_image->GetSpacing()[2];

    bounding_image->TransformPhysicalPointToIndex(minPoint, minIndex);
    bounding_image->TransformPhysicalPointToIndex(maxPoint, maxIndex);

    typename ImageType::RegionType region;
    region.SetIndex(minIndex);
    region.SetUpperIndex(maxIndex);

    std::cout << minIndex << std::endl;
    std::cout << maxIndex << std::endl;

    deformer->SetSourceLandmarks(sourcePointsContainer.GetPointer());
    deformer->SetTargetLandmarks(targetPointsContainer.GetPointer());

    deformer->SetOutputRegion(region);
    deformer->Update();
    //  deformer->UpdateLargestPossibleRegion();

    DisplacementFieldType::Pointer field = deformer->GetOutput();

    auto transform = DisplacementFieldTransformType::New();
    transform->SetDisplacementField(field);

    // Warp template mesh as sanity check
    auto warped_mesh = vtkSmartPointer<vtkPolyData>::New();
    warped_mesh->DeepCopy(template_mesh);

    for (int i = 0; i < warped_mesh->GetNumberOfPoints(); i++)
    {
        auto p = warped_mesh->GetPoint(i);
        auto p_t = transform->TransformPoint(p);
        warped_mesh->GetPoints()->SetPoint(i, p_t.GetDataPointer());
    }

    auto outpath = out_dir + "/test.vtp";

    std::cout << outpath << std::endl;

    auto pdw = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    pdw->SetFileName(outpath.c_str());
    pdw->SetInputData(warped_mesh);
    pdw->Update();

    std::cout << "Wrote" << std::endl;
    
    // Read in and warp template s-rep 
    vtkSmartPointer<vtkXMLDataParser> parser = vtkSmartPointer<vtkXMLDataParser>::New();

    parser->SetFileName(srep.c_str());
    parser->SetIgnoreCharacterData(0);

    std::string upFileName, downFileName, crestFileName;
    std::string out_upFileName, out_downFileName, out_crestFileName;
    int r, c;

    if( parser->Parse() == 1)
    {
        vtkXMLDataElement *root = parser->GetRootElement();
        int numElements = root->GetNumberOfNestedElements();;
        for(int i = 0; i < numElements; ++i)
        {
            double crestShift;
            char *pEnd;
            vtkXMLDataElement *e = root->GetNestedElement(i);
            std::string estimatePath;
            estimatePath = vtksys::SystemTools::GetFilenamePath(srep) + "/";
            std::vector<std::string> components;
            components.push_back(estimatePath);
            auto out_estimatePath = vtksys::SystemTools::GetFilenamePath(outpath) + "/";
            std::cout << "out_estimatePath: " << out_estimatePath << std::endl;
            std::vector<std::string> out_components;
            out_components.push_back(out_estimatePath);

            char* eName = e->GetName();
            if (strcmp(eName, "nRows") == 0)
            {
                r = static_cast<int>(strtol(e->GetCharacterData(), &pEnd, 10));
            }
            else if (strcmp(eName, "nCols") == 0)
            {
                c = static_cast<int>(strtol(e->GetCharacterData(), &pEnd, 10));
            }
            if(strcmp(eName, "upSpoke") == 0)
            {
                upFileName = e->GetCharacterData();
                if(!vtksys::SystemTools::FileIsFullPath(upFileName))
                {
                    components.push_back(upFileName);
                    out_components.push_back(upFileName);

                    upFileName = vtksys::SystemTools::JoinPath(components);
                }
                else
                {
                    out_components.push_back(vtksys::SystemTools::GetFilenameName(upFileName));
                }
                out_upFileName = vtksys::SystemTools::JoinPath(out_components);
            }
            else if(strcmp(eName, "downSpoke")==0)
            {
                downFileName = e->GetCharacterData();
                if(!vtksys::SystemTools::FileIsFullPath(downFileName))
                {
                    components.push_back(downFileName);
                    out_components.push_back(downFileName);

                    downFileName = vtksys::SystemTools::JoinPath(components);
                }
                else
                {
                    out_components.push_back(vtksys::SystemTools::GetFilenameName(downFileName));
                }
                out_downFileName = vtksys::SystemTools::JoinPath(out_components);
            }
            else if(strcmp(eName, "crestSpoke") == 0)
            {
                crestFileName = e->GetCharacterData();
                if(!vtksys::SystemTools::FileIsFullPath(crestFileName))
                {
                    components.push_back(crestFileName);
                    out_components.push_back(crestFileName);

                    crestFileName = vtksys::SystemTools::JoinPath(components);
                }
                else
                {
                    out_components.push_back(vtksys::SystemTools::GetFilenameName(crestFileName));
                }
                out_crestFileName = vtksys::SystemTools::JoinPath(out_components);
            }
        }
    }

    std::cout << upFileName << std::endl;
    std::cout << downFileName << std::endl;
    std::cout << crestFileName << std::endl;
    std::cout << out_upFileName << std::endl;

    // Process spokes
    auto spokesReader = vtkSmartPointer<vtkXMLPolyDataReader>::New();

    // Up
    spokesReader->SetFileName(upFileName.c_str());
    spokesReader->Update();
    auto up_pd = spokesReader->GetOutput();
    auto warped_up = warpSrepPart(up_pd, transform);

    pdw->SetFileName(out_upFileName.c_str());
    pdw->SetInputData(warped_up);
    pdw->Update();

    // Down
    spokesReader->SetFileName(downFileName.c_str());
    spokesReader->Update();
    auto down_pd = spokesReader->GetOutput();
    auto warped_down = warpSrepPart(down_pd, transform);

    pdw->SetFileName(out_downFileName.c_str());
    pdw->SetInputData(warped_down);
    pdw->Update();

    // Crest
    spokesReader->SetFileName(crestFileName.c_str());
    spokesReader->Update();
    auto crest_pd = spokesReader->GetOutput();
    auto warped_crest = warpSrepPart(crest_pd, transform);

    pdw->SetFileName(out_crestFileName.c_str());
    pdw->SetInputData(warped_crest);
    pdw->Update();

    // Write warped srep header
    auto warped_header_file = out_dir + "/header.xml";

    ofstream warped_header;
    warped_header.open(warped_header_file);
    warped_header << "<s-rep>" << std::endl;
    warped_header << "  <nRows>" << r << "</nRows>" << std::endl;
    warped_header << "  <nCols>" << c << "</nCols>" << std::endl;
    warped_header << "  <meshType>" << "Quad" << "</meshType>" << std::endl;
    warped_header << "  <color>" << std::endl;
    warped_header << "    <red>" << "0" << "</red>" << std::endl;
    warped_header << "    <green>" << "0.5" << "</green>" << std::endl;
    warped_header << "    <blue>" << "0" << "</blue>" << std::endl;
    warped_header << "  </color>" << std::endl;
    warped_header << "  <isMean>" << "False" << "</isMean>" << std::endl;
    warped_header << "  <meanStatPath/>" << std::endl;
    warped_header << "  <upSpoke>" << out_upFileName << "</upSpoke>" << std::endl;
    warped_header << "  <downSpoke>" << out_downFileName << "</downSpoke>" << std::endl;
    warped_header << "  <crestSpoke>" << out_crestFileName << "</crestSpoke>" << std::endl;
    warped_header << "</s-rep>" << std::endl;
    warped_header.close();

    return EXIT_SUCCESS;
  }

} // end of anonymous namespace

int main(int argc, char * argv[])
{
  PARSE_ARGS;

  itk::ImageIOBase::IOPixelType     pixelType;
  itk::ImageIOBase::IOComponentType componentType;

  try
  {
    DoIt(argc, argv);
  }

  catch (itk::ExceptionObject & excep)
  {
    std::cerr << argv[0] << ": exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
