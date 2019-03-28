#include "vtkAntiAlias.h"
#include "itkVTKImageToImageFilter.h"
#include "itkImageToVTKImageFilter.h"
#include <itkImage.h>
#include "itkImageFileReader.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageFileWriter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkCastImageFilter.h"
#include <itkNumericTraits.h>
#include <time.h>
#include "vtkImageMagnitude.h"

#include <vector>

vtkAntiAlias::vtkAntiAlias()
{
    
}

void vtkAntiAlias::Filter(vtkImageData *input, vtkSmartPointer<vtkImageData> output)
{
    std::cout << "Started anti-aliasing ..." << std::endl;
    using InputImageType = itk::Image<unsigned char, 3 >;
    
    using ImageType = itk::Image< float, 3 >;
    using FilterType = itk::VTKImageToImageFilter< InputImageType >;
    
    vtkSmartPointer< vtkImageMagnitude > magnitude =
       vtkSmartPointer< vtkImageMagnitude >::New();
    magnitude->SetInputData(input );
    magnitude->Update();
    
    FilterType::Pointer filter = FilterType::New();
    filter->SetInput(magnitude->GetOutput());
    try
    {
        filter->Update();
        
        using CastFilterType = itk::CastImageFilter< InputImageType, ImageType >;
        CastFilterType::Pointer castFilter = CastFilterType::New();
        castFilter->SetInput( filter->GetOutput() );
        
        ImageType::Pointer image = castFilter->GetOutput();
        
        int xmax = image->GetLargestPossibleRegion().GetSize()[0];
        int ymax = image->GetLargestPossibleRegion().GetSize()[1];
        int zmax = image->GetLargestPossibleRegion().GetSize()[2];
        
        ImageType::IndexType start;
        start[0] = 0;
        start[1] = 0;
        start[2] = 0;
        
        ImageType::SizeType size;
        size[0] = xmax;
        size[1] = ymax;
        size[2] = zmax;
        
        std::vector<int> xs;
        std::vector<int> ys;
        std::vector<int> zs;
        
        // Need to extract raw image values from ITK data structure, there is probably a better way to do this
        // ineighbors is kind of confusing at first: ineighbors[x,y,z] holds the linear index into xs, ys, zs arrays
        // of x,y,z: this is so you can easily find out which neighbors of x,y,z are in the computation band by checking
        // to see if ineighbors[x0,y0,z0] == -1
        
        float *phi = new float[xmax*ymax*zmax];
        int *ineighbors = new int[xmax*ymax*zmax];
        
        float val = 0.0;
        ImageType::IndexType pixelIndex;
        
        double dt = 0.0625;
        double bandwidth = 3;
        
        // We don't want to waste time by optimizing voxels which don't matter for the surface (0 level set).
        
        std::cout << "Finding band of computation...";
        
        itk::ImageRegionIterator<ImageType> imageIterator(image,image->GetLargestPossibleRegion());
        
        while (!imageIterator.IsAtEnd())
        {
            val = imageIterator.Get();
            pixelIndex = imageIterator.GetIndex();
            
            phi[ pixelIndex[0] + pixelIndex[1]*xmax + pixelIndex[2]*xmax*ymax ] = val;
            
            
            if (val > -bandwidth && val < bandwidth)
            {
                xs.push_back(pixelIndex[0]);
                ys.push_back(pixelIndex[1]);
                zs.push_back(pixelIndex[2]);
                ineighbors[pixelIndex[0] + pixelIndex[1]*xmax + pixelIndex[2]*xmax*ymax] = xs.size() - 1;
            }
            else
            {
                ineighbors[pixelIndex[0] + pixelIndex[1]*xmax + pixelIndex[2]*xmax*ymax] = -1;
            }
            
            ++imageIterator;
        }
        
        int num = xs.size();
        
        std::cout << "done. Number of points in band: " << num << " of " << xmax * ymax * zmax << std::endl;
        
        float *phi0 = new float[num];
        
        std::cout << "Building temp images...";
        
        for (int i = 0; i < num; i++)
        {
            phi0[i] = phi[xs[i] + ys[i]*xmax + zs[i]*xmax*ymax];
        }
        
        std::cout << "done." << std::endl;
        
        float phix;
        float phiy;
        float phiz;
        float phixx;
        float phiyy;
        float phizz;
        float phixy;
        float phiyz;
        float phixz;
        float phix2;
        float phiy2;
        float phiz2;
        float denom, disc;
        
        float phinew;
        float *H = new float[num];
        float *K = new float[num];
        float *K1 = new float[num];
        float *K2 = new float[num];
        float *dH = new float[num];
        float *dK1 = new float[num];
        float *dK2 = new float[num];
        float *update = new float[num];
        float *ndenom = new float[num];
        
        float Hsum, K1sum, K2sum;
        
        
        int x, y, z;
        int count;
        
        std::cout << "Iterating: ";
        
        for (int iter = 0; iter < 1000; iter++)
        {
            if (iter % 100 == 0)
                std::cout << iter << "..." << std::flush;
            
#pragma omp parallel for
            for (int i = 0; i < num; i++)
            {
                x = xs[i];
                y = ys[i];
                z = zs[i];
                
                phix = 0.5 * (phi[(x+1) + (y)*xmax + (z)*xmax*ymax] - phi[(x-1) + (y)*xmax + (z)*xmax*ymax]);
                phiy = 0.5 * (phi[(x) + (y+1)*xmax + (z)*xmax*ymax] - phi[(x) + (y-1)*xmax + (z)*xmax*ymax]);
                phiz = 0.5 * (phi[(x) + (y)*xmax + (z+1)*xmax*ymax] - phi[(x) + (y)*xmax + (z-1)*xmax*ymax]);
                
                phixx = phi[(x-1) + (y)*xmax + (z)*xmax*ymax] - 2*phi[(x) + (y)*xmax + (z)*xmax*ymax] + phi[(x+1) + (y)*xmax + (z)*xmax*ymax];
                phiyy = phi[(x) + (y-1)*xmax + (z)*xmax*ymax] - 2*phi[(x) + (y)*xmax + (z)*xmax*ymax] + phi[(x) + (y+1)*xmax + (z)*xmax*ymax];
                phizz = phi[(x) + (y)*xmax + (z-1)*xmax*ymax] - 2*phi[(x) + (y)*xmax + (z)*xmax*ymax] + phi[(x) + (y)*xmax + (z+1)*xmax*ymax];
                
                phixy = 0.25 * (phi[(x+1) + (y+1)*xmax + (z)*xmax*ymax] + phi[(x-1) + (y-1)*xmax + (z)*xmax*ymax]
                        - phi[(x+1) + (y-1)*xmax + (z)*xmax*ymax] - phi[(x-1) + (y+1)*xmax + (z)*xmax*ymax]);
                phiyz = 0.25 * (phi[(x) + (y+1)*xmax + (z+1)*xmax*ymax] + phi[(x) + (y-1)*xmax + (z-1)*xmax*ymax]
                        - phi[(x) + (y+1)*xmax + (z-1)*xmax*ymax] - phi[(x) + (y-1)*xmax + (z+1)*xmax*ymax]);
                phixz = 0.25 * (phi[(x+1) + (y)*xmax + (z+1)*xmax*ymax] + phi[(x-1) + (y)*xmax + (z-1)*xmax*ymax]
                        - phi[(x+1) + (y)*xmax + (z-1)*xmax*ymax] - phi[(x-1) + (y)*xmax + (z+1)*xmax*ymax]);
                
                phix2 = phix * phix;
                phiy2 = phiy * phiy;
                phiz2 = phiz * phiz;
                
                denom = phix2 + phiy2 + phiz2;
                ndenom[i] = sqrt(denom);
                
                if (denom < 0.0002)
                {
                    H[i] = 0;
                    K[i] = 0;
                }
                else
                {
                    H[i] = ( (phiyy+phizz)*phix2 + (phixx+phizz)*phiy2 + (phixx+phiyy)*phiz2 -
                             2*phix*phiy*phixy - 2*phix*phiz*phixz - 2*phiy*phiz*phiyz ) / ( ndenom[i] * ndenom[i] * ndenom[i] );
                    
                    K[i] = ( phix2*(phiyy*phizz - phiyz*phiyz) + phiy2*(phixx*phizz - phixz*phixz) + phiz2*(phixx*phiyy - phixy*phixy) +
                             2*( phix*phiy*(phixz*phiyz - phixy*phizz) + phiy*phiz*(phixy*phixz - phiyz*phixx)
                                 + phix*phiz*(phixy*phiyz - phixz*phiyy) ) ) / (denom * denom);
                }
                
                disc = (H[i]*H[i] - 2*K[i]);
                
                if (disc > 0)
                {
                    K1[i] = 0.5 * (H[i] + sqrt(disc));
                    K2[i] = 0.5 * (H[i] - sqrt(disc));
                }
                else
                {
                    K1[i] = 0;
                    K2[i] = 0;
                }
            }
            
            // Up to this is parallelizable, need to start new block after so all H are computed
#pragma omp parallel for
            for (int i = 0; i < num; i++)
            {
                x = xs[i];
                y = ys[i];
                z = zs[i];
                
                Hsum = 0.0;
                K1sum = 0.0;
                K2sum = 0.0;
                count = 0;
                
                // Check all 6 neighbors for current voxel -- this could probably be cleaner
                
                if (ineighbors[(x-1) + (y)*xmax + (z)*xmax*ymax] >= 0)
                {
                    Hsum += H[ineighbors[(x-1) + (y)*xmax + (z)*xmax*ymax]];
                    K1sum = K1[ineighbors[(x-1) + (y)*xmax + (z)*xmax*ymax]];
                    K2sum = K2[ineighbors[(x-1) + (y)*xmax + (z)*xmax*ymax]];
                    count++;
                }
                if (ineighbors[(x+1) + (y)*xmax + (z)*xmax*ymax] >= 0)
                {
                    Hsum += H[ineighbors[(x+1) + (y)*xmax + (z)*xmax*ymax]];
                    K1sum += K1[ineighbors[(x+1) + (y)*xmax + (z)*xmax*ymax]];
                    K2sum += K2[ineighbors[(x+1) + (y)*xmax + (z)*xmax*ymax]];
                    count++;
                }
                if (ineighbors[(x) + (y-1)*xmax + (z)*xmax*ymax] >= 0)
                {
                    Hsum += H[ineighbors[(x) + (y-1)*xmax + (z)*xmax*ymax]];
                    K1sum += K1[ineighbors[(x) + (y-1)*xmax + (z)*xmax*ymax]];
                    K2sum += K2[ineighbors[(x) + (y-1)*xmax + (z)*xmax*ymax]];
                    count++;
                }
                if (ineighbors[(x) + (y+1)*xmax + (z)*xmax*ymax] >= 0)
                {
                    Hsum += H[ineighbors[(x) + (y+1)*xmax + (z)*xmax*ymax]];
                    K1sum += K1[ineighbors[(x) + (y+1)*xmax + (z)*xmax*ymax]];
                    K2sum += K2[ineighbors[(x) + (y+1)*xmax + (z)*xmax*ymax]];
                    count++;
                }
                if (ineighbors[(x) + (y)*xmax + (z-1)*xmax*ymax] >= 0)
                {
                    Hsum += H[ineighbors[(x) + (y)*xmax + (z-1)*xmax*ymax]];
                    K1sum += K1[ineighbors[(x) + (y)*xmax + (z-1)*xmax*ymax]];
                    K2sum += K2[ineighbors[(x) + (y)*xmax + (z-1)*xmax*ymax]];
                    count++;
                }
                if (ineighbors[(x) + (y)*xmax + (z+1)*xmax*ymax] >= 0)
                {
                    Hsum += H[ineighbors[(x) + (y)*xmax + (z+1)*xmax*ymax]];
                    K1sum += K1[ineighbors[(x) + (y)*xmax + (z+1)*xmax*ymax]];
                    K2sum += K2[ineighbors[(x) + (y)*xmax + (z+1)*xmax*ymax]];
                    count++;
                }
                
                if (count > 0)
                {
                    dH[i] = H[i] - (Hsum / count);
                    dK1[i] = K1[i] - (K1sum / count);
                    dK2[i] = K2[i] - (K2sum / count);
                }
                else
                {
                    dH[i] = 0;
                    dK1[i] = 0;
                    dK2[i] = 0;
                }
                
                if (K[i] < 0)
                    update[i] = dH[i];
                else // K >= 0
                {
                    if (H[i] > 0) // K >= 0 & H > 0
                        update[i] = (dK1[i] > 0 ? dK1[i] : 0); // update = max(dK1,0)
                    else // K >= 0 & H < 0
                        update[i] = (dK2[i] < 0 ? dK2[i] : 0); // update = min(dK2,0)
                }
                
                if (update[i] > 100)
                    update[i] = 100;
                if (update[i] < -100)
                    update[i] = -100;
                
                update[i] = update[i] * ndenom[i];
                
                phinew = phi[x + y*xmax + z*xmax*ymax] + dt * update[i];
                
                if (phinew > phi0[i] + 0.15)
                    phinew = phi0[i] + 0.15;
                if (phinew < phi0[i] - 0.15)
                    phinew = phi0[i] - 0.15;
                
                phi[x + y*xmax + z*xmax*ymax] = phinew;
                
            }
        }
        
        std::cout << "Iterations done, writing output." << std::endl;
        
        pixelIndex[0] = 0;
        pixelIndex[1] = 0;
        pixelIndex[2] = 0;
        
        imageIterator.SetIndex(pixelIndex);
        
        while(!imageIterator.IsAtEnd())
        {
            pixelIndex = imageIterator.GetIndex();
            
            imageIterator.Set(phi[pixelIndex[0] + pixelIndex[1]*xmax + pixelIndex[2]*xmax*ymax]);
            
            ++imageIterator;
        }
        
        typedef itk::ImageFileWriter< ImageType > WriterType;
        WriterType::Pointer writer = WriterType::New();
        writer->SetFileName("/playpen/workspace/newuoa/test.mhd");
        writer->SetInput(image);
        writer->Update();
        std::cout << "Finished anti-aliasing." << std::endl;
        // Convert itk image back to vtk image data
        using BackFilterType = itk::ImageToVTKImageFilter< ImageType >;
        BackFilterType::Pointer backFilter = BackFilterType::New();
        backFilter->SetInput(image);
        backFilter->Update();
        
        output->DeepCopy(backFilter->GetOutput());
        
    }
    catch( itk::ExceptionObject & error )
    {
        std::cerr << "Error: " << error << std::endl;
        return;
    }
}
