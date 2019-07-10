#include "vtkImage2SignedDistanceMap.h"

#include "itkVTKImageToImageFilter.h"

#include <itkImageRegionConstIterator.h>
#include <itkConstNeighborhoodIterator.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkConstantPadImageFilter.h>
#include <itkSignedDanielssonDistanceMapImageFilter.h>
#include "itkImageToVTKImageFilter.h"
#include <itkBinaryBallStructuringElement.h>
#include <itkBinaryDilateImageFilter.h>
#include <itkBinaryErodeImageFilter.h>
#include <itkErodeObjectMorphologyImageFilter.h>
#include <itkDilateObjectMorphologyImageFilter.h>
#include <itkRegionOfInterestImageFilter.h>
#include <stack>
#include <math.h>
#include <iomanip>
#include "vtkImageMagnitude.h"
vtkImage2SignedDistanceMap::vtkImage2SignedDistanceMap()
{
    
}


void vtkImage2SignedDistanceMap::Convert(vtkImageData *input, vtkSmartPointer<vtkImageData> output)
{
    std::cout << "Started to convert image to signed distance map..." << std::endl;
    using FilterType = itk::VTKImageToImageFilter< ImageType >;
   
    vtkSmartPointer< vtkImageMagnitude > magnitude =
       vtkSmartPointer< vtkImageMagnitude >::New();
    magnitude->SetInputData(input );
    magnitude->Update();
    FilterType::Pointer filter = FilterType::New();
    filter->SetInput(magnitude->GetOutput());
    
    try
    {
        filter->Update();
        
        ImagePointer image = filter->GetOutput();
        typedef itk::NeighborhoodIterator<ImageType> IteratorType;
        IteratorType::RadiusType radius;
        radius.Fill(1);
        using namespace std;
        cout << "Finding all boundary voxels ... " << endl;
        stack<ImageType::IndexType> boundary;                // A stack to keep track of all changing boundary voxels.
        IteratorType it(radius, image, image->GetRequestedRegion());
        bool inBounds = false;
        for( it.GoToBegin(); !it.IsAtEnd(); ++it ) {
            Pixel p = it.GetCenterPixel();
            for( unsigned int ni = 0; ni != it.Size(); ++ni ) {
                // If voxel has any non-matching neighbor, then it's an edge voxel. Add it to the
                // processing stack.
                if( p != it.GetPixel(ni, inBounds) && inBounds ) {
                    boundary.push(it.GetIndex());
                    break;
                }
            }
        }
        cout << boundary.size() << " voxels" << endl;
        cout << "Eroding pimples and dimples ... " << endl;
        unsigned int pts = 0;
        unsigned int changed = 0;
        while( !boundary.empty() ) {
            // Note central pixel type.
            it.SetLocation( boundary.top() );
            boundary.pop();
            pts++;
            const Pixel p = it.GetCenterPixel();
            // If central pixel has less than x neighbors of the same type, then reverse its
            // classification (foreground becomes background and vice-versa).
            // Strictly speaking this should be a user-controllable parameter, but 9 is a good
            // number as it represents an entire plane of pixels. Anything higher than that
            // might result in a line sitting on top of a plane to be eroded off -- we do not
            // want that to happen.
            const unsigned int threshold = 9;
            unsigned int count = 0;
            for( int xi = -1; xi <= 1; ++xi ) {
                for( int yi = -1; yi <= 1; ++yi ) {
                    for( int zi = -1; zi <= 1; ++zi ) {
                        const IteratorType::OffsetType np = {{xi,yi,zi}};
                        if( it.GetPixel( np ) == p ) {
                            count++;
                            if( count > threshold ) {
                                goto doneWithCheckingThisVoxel;
                            }
                            bool hasANeighbor = false;
                            // Check if this pixel has a neighbor.
                            for( int xni = max(xi-1,-1); xni <= min(xi+1,+1); ++xni ) {
                                for( int yni = max(yi-1,-1); yni <= min(yi+1,+1); ++yni ) {
                                    for( int zni = max(zi-1,-1); zni <= min(zi+1,+1); ++zni ) {
                                        if(xni != 0 && yni != 0 && zni != 0 ) {
                                            const IteratorType::OffsetType np = {{xni,yni,zni}};
                                            if( p == it.GetPixel(np, inBounds) && inBounds ) {
                                                hasANeighbor = true;
                                                goto doneWithCheckingIfThisVoxelHasANeighbor;
                                            }
                                        }
                                    }
                                }
                            }
doneWithCheckingIfThisVoxelHasANeighbor:
                            if( !hasANeighbor ) {
                                // This pixel has no other connected neighbor. The only
                                // connection goes through the central pixel. Therefore, under
                                // no circumstance, erode the central pixel.
                                count = threshold+1;
                                goto doneWithCheckingThisVoxel;
                            }
                        }
                    }
                }
            }
doneWithCheckingThisVoxel:
            if( count <= threshold ) {
                it.SetCenterPixel( (p) ? 0 : 1 );
                changed++;
                // Push all neighbors with same pixel type on the boundary stack.
                for( unsigned int ni = 0; ni != it.Size(); ++ni ) {
                    if( p == it.GetPixel(ni, inBounds) && inBounds ) {
                        boundary.push(it.GetIndex(ni));
                    }
                }
            }
        }
        image = CloseImage( image, image->GetSpacing()[0] - 0.0001 );
        // Do morphological operations primarily in the xy plane. Along the z axis, only get rid of
        // edges that are sticking out.
        // This is not a fundamentally good idea as they can erode away single voxel wide regions.
        //image = openImage( image, image->GetSpacing()[0] - 0.0001 );
        
        //                writeImage<ImageType>(image, clearedImage.c_str() );
        
        itk::SignedDanielssonDistanceMapImageFilter<ImageType,ImageType>::Pointer ddmFilter =
                itk::SignedDanielssonDistanceMapImageFilter<ImageType,ImageType>::New();
        ddmFilter->SetInput(image);
        ddmFilter->Update();
        typedef itk::ImageFileWriter< ImageType > WriterType;
        WriterType::Pointer writer = WriterType::New();
        writer->SetFileName("/playpen/workspace/newuoa/testsdm.mhd");
        writer->SetInput(ddmFilter->GetOutput());
        writer->Update();
        // Convert itk image back to vtk image data
        using BackFilterType = itk::ImageToVTKImageFilter< ImageType >;
        BackFilterType::Pointer backFilter = BackFilterType::New();
        backFilter->SetInput(ddmFilter->GetOutput());
        backFilter->Update();
        
        output->DeepCopy(backFilter->GetOutput());
        std::cout << "Finished generating distance map." << std::endl;
        
    }
    catch( itk::ExceptionObject & error )
    {
        std::cerr << "Error: " << error << std::endl;
        return;
    }
    
}

vtkImage2SignedDistanceMap::ImagePointer vtkImage2SignedDistanceMap::CloseImage(vtkImage2SignedDistanceMap::ImagePointer image, double radius)
{
    // The ball structuring element:
    typedef itk::BinaryBallStructuringElement<Pixel, 3> Kernel;
    // The dilation filter
    //typedef DilateObjectMorphologyImageFilter<ImageType, ImageType, Kernel> DilateFilter;
    typedef itk::BinaryDilateImageFilter<ImageType, ImageType, Kernel> DilateFilter;
    // The erosion filter
    //typedef ErodeObjectMorphologyImageFilter<ImageType, ImageType, Kernel> ErodeFilter;
    typedef itk::BinaryErodeImageFilter<ImageType, ImageType, Kernel> ErodeFilter;

    // Create the structuring element:
    ImageType::SpacingType spacing = image->GetSpacing();
    Kernel ball;
    Kernel::SizeType ballSize;
    ballSize[0] = (unsigned int)ceil(radius/spacing[0]);
    ballSize[1] = (unsigned int)ceil(radius/spacing[1]);
    ballSize[2] = (unsigned int)ceil(radius/spacing[2]);
    cout << "Creating structuring element ["
             << ballSize[0] << ", " << ballSize[1] << ", " << ballSize[2]
             << "] ... " << endl;
    ball.SetRadius(ballSize);
    ball.CreateStructuringElement();
    cout << "done" << endl;
//        writeImage<ImageType>(image, "original.mhd");

    // Before dilating, pad the image, so that structures near the edges do not get messed up.
    typedef itk::ConstantPadImageFilter<ImageType, ImageType> PadFilter;

    PadFilter::Pointer padFilter = PadFilter::New();
    padFilter->SetInput(image);
    padFilter->SetConstant(0);
    const unsigned long pad = std::max(std::max(ballSize[0],ballSize[1]),ballSize[2]) + 1;
    const unsigned long padding[3] = { pad, pad, pad };
    padFilter->SetPadLowerBound(padding);
    padFilter->SetPadUpperBound(padding);
    padFilter->Update();

    // Now do the close operation
    cout << "Dilating ... " << endl;
    DilateFilter::Pointer closeDilate        = DilateFilter::New();
    //closeDilate->SetObjectValue(1);
    closeDilate->SetForegroundValue(1);
    closeDilate->SetKernel(ball);
    closeDilate->SetInput(padFilter->GetOutput());
    closeDilate->Update();
    cout << "done" << endl;
//        writeImage<ImageType>(closeDilate->GetOutput(), "dilated.mhd");

    cout << "Eroding ... " << endl;
    ErodeFilter::Pointer closeErode        = ErodeFilter::New();
    //closeErode->SetObjectValue(1);
    closeErode->SetForegroundValue(1);
    closeErode->SetKernel(ball);
    closeErode->SetInput(closeDilate->GetOutput());
    closeErode->Update();
    cout << "done" << endl;
//        writeImage<ImageType>(closeErode->GetOutput(), "closed.mhd");

    // Remove the extra padding that we had added.
    ImageType::RegionType desiredRegion;
    // ConstantPadImageFilter is crazy -- it creates images with negative indices.
    ImageType::IndexType idxLower = {{ 0, 0, 0 }};
    desiredRegion.SetIndex( idxLower );
    desiredRegion.SetSize( image->GetLargestPossibleRegion().GetSize() );
    // Setup the region of interest filter
    typedef itk::RegionOfInterestImageFilter<ImageType, ImageType> ROIFilter;
    ROIFilter::Pointer roi        = ROIFilter::New();
    roi->SetInput( closeErode->GetOutput() );
    roi->SetRegionOfInterest( desiredRegion );
    roi->Update();

    return roi->GetOutput();
}
